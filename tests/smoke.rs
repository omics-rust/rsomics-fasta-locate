use std::path::Path;
use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fasta-locate"))
}

fn fixture() -> &'static Path {
    Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/small.fa"
    ))
}

#[test]
fn finds_atgc_both_strands() {
    let out = bin().args(["-p", "ATGC"]).arg(fixture()).output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    // Header + at least one hit
    assert!(s.contains("seqID\tpatternName"));
    let hits: Vec<&str> = s.lines().filter(|l| !l.starts_with("seqID")).collect();
    assert!(!hits.is_empty(), "expected at least one hit");
    // Both strands present
    assert!(
        hits.iter()
            .any(|h| h.contains('\t') && h.split('\t').nth(3) == Some("+"))
    );
    assert!(
        hits.iter()
            .any(|h| h.contains('\t') && h.split('\t').nth(3) == Some("-"))
    );
}

#[test]
fn positive_strand_only_flag() {
    let out = bin()
        .args(["-p", "ATGC", "-P"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    let hits: Vec<&str> = s.lines().skip(1).collect();
    for hit in &hits {
        let strand = hit.split('\t').nth(3).unwrap_or("");
        assert_eq!(strand, "+", "expected only + strand hits");
    }
}

#[test]
fn bed_output_format() {
    let out = bin()
        .args(["-p", "ATGC", "--bed"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    // BED: no header; each line has 6 tab-separated fields
    for line in s.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 6, "BED line must have 6 fields: {line}");
        // BED start is 0-based: parse as usize
        let _start: usize = fields[1].parse().expect("BED start must be integer");
        let _end: usize = fields[2].parse().expect("BED end must be integer");
    }
}

#[test]
fn no_pattern_exits_nonzero() {
    let out = bin().arg(fixture()).output().unwrap();
    assert!(!out.status.success(), "missing -p should fail");
}

#[test]
fn hide_matched_column() {
    let out = bin()
        .args(["-p", "ATGC", "-M"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    // Header should have 6 columns (no matched)
    let header = s.lines().next().unwrap();
    assert_eq!(header.split('\t').count(), 6);
}

#[test]
fn regex_mode() {
    let out = bin()
        .args(["-r", "-p", "ATG."])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    let hits: Vec<&str> = s.lines().skip(1).collect();
    assert!(!hits.is_empty());
}

#[test]
fn no_hits_returns_empty_body() {
    let out = bin().args(["-p", "ZZZZZ"]).arg(fixture()).output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    // Only header line
    assert_eq!(s.lines().count(), 1);
}

#[test]
fn ignore_case_finds_lowercase_seq() {
    // chr3 has "atgcATGC"; without -i only uppercase ATGC matches from pos 5
    let out_sensitive = bin()
        .args(["-p", "ATGC", "-P"])
        .arg(fixture())
        .output()
        .unwrap();
    let out_insensitive = bin()
        .args(["-p", "ATGC", "-i", "-P"])
        .arg(fixture())
        .output()
        .unwrap();
    let n_sensitive = String::from_utf8_lossy(&out_sensitive.stdout)
        .lines()
        .filter(|l| l.contains("chr3"))
        .count();
    let n_insensitive = String::from_utf8_lossy(&out_insensitive.stdout)
        .lines()
        .filter(|l| l.contains("chr3"))
        .count();
    assert!(
        n_insensitive > n_sensitive,
        "ignore-case should find more chr3 hits"
    );
}
