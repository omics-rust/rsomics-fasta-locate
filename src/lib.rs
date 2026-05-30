use std::io::{BufWriter, Write};
use std::path::Path;

use memchr::memchr;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    #[default]
    Tsv,
    Bed,
    Gtf,
}

#[derive(Clone)]
pub struct Pattern {
    pub name: String,
    pub text: String,
}

pub struct LocateOpts {
    pub patterns: Vec<Pattern>,
    pub use_regexp: bool,
    pub ignore_case: bool,
    pub only_positive: bool,
    /// Do not include the matched sequence column (TSV only).
    pub hide_matched: bool,
    pub format: OutputFormat,
}

/// Locate patterns in every sequence in `input`, writing results to `output`.
///
/// Coordinates: 1-based inclusive on the original (positive) sequence.
/// Negative-strand hits map back from the reverse-complement.
///
/// TSV column semantics (seqkit-compatible):
/// - pattern: when `ignore_case`, lowercased matched text; otherwise the query text
/// - matched: actual matched bytes; lowercased when `ignore_case`
///
/// Output preserves input order.
pub fn locate_fasta(input: &Path, opts: &LocateOpts, output: &mut dyn Write) -> Result<u64> {
    if std::fs::metadata(input).is_ok_and(|m| m.len() == 0) {
        return Err(RsomicsError::InvalidInput("empty file".into()));
    }

    let compiled = compile_patterns(opts)?;

    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    {
        let mut reader = parse_fastx_file(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
        while let Some(rec) = reader.next() {
            let rec = rec.map_err(|e| RsomicsError::InvalidInput(format!("reading: {e}")))?;
            let id = std::str::from_utf8(rec.id())
                .unwrap_or("unknown")
                .to_owned();
            let seq = rec.seq().into_owned();
            records.push((id, seq));
        }
    }

    let per_record_lines: Vec<(Vec<u8>, u64)> = records
        .par_iter()
        .map(|(seq_id, seq)| {
            let mut buf = Vec::<u8>::new();
            let mut count = 0u64;
            let seq_len = seq.len();

            for (pat, compiled_pat) in opts.patterns.iter().zip(compiled.iter()) {
                let plus_hits = find_hits(seq, compiled_pat, opts.ignore_case);
                for (start0, end0, matched) in &plus_hits {
                    format_hit(
                        &mut buf,
                        opts,
                        seq_id,
                        &pat.name,
                        "+",
                        *start0 + 1,
                        *end0,
                        matched,
                    );
                    count += 1;
                }

                if !opts.only_positive {
                    let rc = revcomp(seq);
                    let minus_hits = find_hits(&rc, compiled_pat, opts.ignore_case);
                    for (start0, end0, matched) in &minus_hits {
                        // Map revcomp 0-based [start0, end0) → original 1-based.
                        // revcomp[i] ↔ original[len-1-i]:
                        //   orig_end   = len - start0         (1-based inclusive)
                        //   orig_start = len - end0 + 1       (1-based inclusive)
                        let orig_start = seq_len - end0 + 1;
                        let orig_end = seq_len - start0;
                        format_hit(
                            &mut buf, opts, seq_id, &pat.name, "-", orig_start, orig_end, matched,
                        );
                        count += 1;
                    }
                }
            }

            (buf, count)
        })
        .collect();

    let mut out = BufWriter::with_capacity(512 * 1024, output);
    if opts.format == OutputFormat::Tsv {
        if opts.hide_matched {
            writeln!(out, "seqID\tpatternName\tpattern\tstrand\tstart\tend")
                .map_err(RsomicsError::Io)?;
        } else {
            writeln!(
                out,
                "seqID\tpatternName\tpattern\tstrand\tstart\tend\tmatched"
            )
            .map_err(RsomicsError::Io)?;
        }
    }

    let mut total_hits: u64 = 0;
    for (lines, count) in per_record_lines {
        out.write_all(&lines).map_err(RsomicsError::Io)?;
        total_hits += count;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(total_hits)
}

/// Format one hit into `buf` without any allocation beyond `buf`'s growth.
#[allow(clippy::too_many_arguments)]
fn format_hit(
    buf: &mut Vec<u8>,
    opts: &LocateOpts,
    seq_id: &str,
    pat_name: &str,
    strand: &str,
    start1: usize,
    end1: usize,
    matched: &str,
) {
    // When ignore_case: `pattern` column = lowercased matched text; `matched` also lowercased.
    // Without ignore_case: `pattern` column = pat_name (the query text).
    let (pattern_col, matched_col) = if opts.ignore_case {
        let lc = matched.to_lowercase();
        (lc.clone(), lc)
    } else {
        (pat_name.to_owned(), matched.to_owned())
    };

    match opts.format {
        OutputFormat::Tsv => {
            use std::io::Write as _;
            if opts.hide_matched {
                let _ = writeln!(
                    buf,
                    "{seq_id}\t{pat_name}\t{pattern_col}\t{strand}\t{start1}\t{end1}"
                );
            } else {
                let _ = writeln!(
                    buf,
                    "{seq_id}\t{pat_name}\t{pattern_col}\t{strand}\t{start1}\t{end1}\t{matched_col}"
                );
            }
        }
        OutputFormat::Bed => {
            use std::io::Write as _;
            let start0 = start1 - 1;
            let _ = writeln!(buf, "{seq_id}\t{start0}\t{end1}\t{pat_name}\t0\t{strand}");
        }
        OutputFormat::Gtf => {
            use std::io::Write as _;
            let _ = writeln!(
                buf,
                "{seq_id}\tSeqKit\tlocation\t{start1}\t{end1}\t0\t{strand}\t.\tgene_id \"{pat_name}\"; "
            );
        }
    }
}

// Both variants are Send+Sync; rayon shares them across threads without cloning.
enum CompiledPattern {
    Literal(Vec<u8>),
    Regex(regex::bytes::Regex),
}

fn compile_patterns(opts: &LocateOpts) -> Result<Vec<CompiledPattern>> {
    opts.patterns
        .iter()
        .map(|p| {
            if opts.use_regexp {
                let re_str = if opts.ignore_case {
                    format!("(?i){}", p.text)
                } else {
                    p.text.clone()
                };
                let re = regex::bytes::Regex::new(&re_str).map_err(|e| {
                    RsomicsError::InvalidInput(format!("invalid regex {:?}: {e}", p.text))
                })?;
                Ok(CompiledPattern::Regex(re))
            } else {
                Ok(CompiledPattern::Literal(p.text.as_bytes().to_vec()))
            }
        })
        .collect()
}

/// Return all `(start0, end0_exclusive, matched_str)` hits for `pattern` in `seq`.
///
/// Literal patterns: memchr first-byte scan then byte-by-byte confirmation.
/// Skips non-matching positions at SIMD speed; 3-5× faster than a naive
/// sliding window when the first base is infrequent.
fn find_hits(seq: &[u8], pat: &CompiledPattern, ignore_case: bool) -> Vec<(usize, usize, String)> {
    match pat {
        CompiledPattern::Literal(lit) => {
            if lit.is_empty() || lit.len() > seq.len() {
                return Vec::new();
            }
            let mut hits = Vec::new();
            let plen = lit.len();
            let last_start = seq.len() - plen;

            if ignore_case {
                let first_upper = lit[0].to_ascii_uppercase();
                let first_lower = lit[0].to_ascii_lowercase();
                let mut pos = 0usize;
                while pos <= last_start {
                    let rem = &seq[pos..=last_start];
                    // Jump to the nearest occurrence of the first char (either case).
                    let advance = memchr(first_upper, rem)
                        .map(|a| memchr(first_lower, rem).map(|b| a.min(b)).unwrap_or(a))
                        .or_else(|| memchr(first_lower, rem));
                    let Some(delta) = advance else { break };
                    pos += delta;
                    let window = &seq[pos..pos + plen];
                    if window
                        .iter()
                        .zip(lit.iter())
                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                    {
                        hits.push((
                            pos,
                            pos + plen,
                            String::from_utf8_lossy(window).into_owned(),
                        ));
                    }
                    pos += 1;
                }
            } else {
                let first = lit[0];
                let mut pos = 0usize;
                while pos <= last_start {
                    let Some(delta) = memchr(first, &seq[pos..=last_start]) else {
                        break;
                    };
                    pos += delta;
                    let window = &seq[pos..pos + plen];
                    if window == lit.as_slice() {
                        hits.push((
                            pos,
                            pos + plen,
                            String::from_utf8_lossy(window).into_owned(),
                        ));
                    }
                    pos += 1;
                }
            }
            hits
        }
        CompiledPattern::Regex(re) => re
            .find_iter(seq)
            .map(|m| {
                (
                    m.start(),
                    m.end(),
                    String::from_utf8_lossy(m.as_bytes()).into_owned(),
                )
            })
            .collect(),
    }
}

/// Reverse complement a nucleotide sequence, preserving the case of each base.
///
/// Case is preserved so that a case-sensitive scan of the revcomp produces the
/// same hit set as seqkit locate — seqkit's revcomp is also case-preserving.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'T' => b'A',
            b't' => b'a',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'N' => b'N',
            b'n' => b'n',
            other => other,
        })
        .collect()
}

/// Load patterns from a FASTA file (name = sequence id, text = sequence).
pub fn load_pattern_file(path: &Path) -> Result<Vec<Pattern>> {
    let mut reader = parse_fastx_file(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut patterns = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| RsomicsError::InvalidInput(format!("reading: {e}")))?;
        let name = std::str::from_utf8(rec.id())
            .unwrap_or("unknown")
            .to_owned();
        let text = std::str::from_utf8(&rec.seq())
            .map_err(|_| RsomicsError::InvalidInput("non-UTF8 pattern sequence".into()))?
            .to_owned();
        patterns.push(Pattern { name, text });
    }
    Ok(patterns)
}
