use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fasta_locate::{LocateOpts, OutputFormat, Pattern, load_pattern_file, locate_fasta};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-fasta-locate",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input FASTA file.
    pub input: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    /// Pattern/motif to locate (may be repeated).
    #[arg(short = 'p', long = "pattern")]
    patterns: Vec<String>,

    /// Pattern/motif file in FASTA format (name = pattern name, seq = pattern).
    #[arg(short = 'f', long = "pattern-file")]
    pattern_file: Option<PathBuf>,

    /// Treat patterns as regular expressions.
    #[arg(short = 'r', long)]
    use_regexp: bool,

    /// Case-insensitive matching.
    #[arg(short = 'i', long)]
    ignore_case: bool,

    /// Search only the positive strand.
    #[arg(short = 'P', long)]
    only_positive_strand: bool,

    /// Do not output the matched sequence column.
    #[arg(short = 'M', long)]
    hide_matched: bool,

    /// Output in BED6 format.
    #[arg(long)]
    bed: bool,

    /// Output in GTF format.
    #[arg(long)]
    gtf: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        if self.patterns.is_empty() && self.pattern_file.is_none() {
            return Err(RsomicsError::InvalidInput(
                "at least one of -p/--pattern or -f/--pattern-file is required".into(),
            ));
        }

        let mut patterns: Vec<Pattern> = self
            .patterns
            .iter()
            .map(|p| Pattern {
                name: p.clone(),
                text: p.clone(),
            })
            .collect();

        if let Some(ref pf) = self.pattern_file {
            patterns.extend(load_pattern_file(pf)?);
        }

        let format = if self.gtf {
            OutputFormat::Gtf
        } else if self.bed {
            OutputFormat::Bed
        } else {
            OutputFormat::Tsv
        };

        let opts = LocateOpts {
            patterns,
            use_regexp: self.use_regexp,
            ignore_case: self.ignore_case,
            only_positive: self.only_positive_strand,
            hide_matched: self.hide_matched,
            format,
        };

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };

        let hits = locate_fasta(&self.input, &opts, &mut out)?;

        if !self.common.quiet {
            eprintln!("{hits} hits");
        }

        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Locate subsequences/motifs in FASTA files — seqkit locate port.",
    origin: Some(Origin {
        upstream: "seqkit locate",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1002/imt2.191"),
    }),
    usage_lines: &["<input.fasta> -p <pattern> [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('p'),
                long: "pattern",
                aliases: &[],
                value: Some("<pattern>"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Pattern/motif (may be repeated).",
                why_default: None,
            },
            FlagSpec {
                short: Some('f'),
                long: "pattern-file",
                aliases: &[],
                value: Some("<file.fa>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "FASTA file of named patterns.",
                why_default: None,
            },
            FlagSpec {
                short: Some('r'),
                long: "use-regexp",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Treat patterns as regular expressions.",
                why_default: None,
            },
            FlagSpec {
                short: Some('i'),
                long: "ignore-case",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Case-insensitive matching.",
                why_default: None,
            },
            FlagSpec {
                short: Some('P'),
                long: "only-positive-strand",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Search only the positive strand.",
                why_default: None,
            },
            FlagSpec {
                short: Some('M'),
                long: "hide-matched",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Omit the matched-sequence column.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Find ATG start codons",
            command: "rsomics-fasta-locate genome.fa -p ATG",
        },
        Example {
            description: "Locate primer on positive strand, BED output",
            command: "rsomics-fasta-locate genome.fa -p ATGCATGC -P --bed",
        },
        Example {
            description: "Regex search for ORF starts",
            command: "rsomics-fasta-locate genome.fa -r -p 'ATG(.{3})+'",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
