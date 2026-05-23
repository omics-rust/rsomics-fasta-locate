use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fasta-locate"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small.fa")
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn seqkit_version() -> Option<String> {
    let out = Command::new("seqkit").arg("version").output().ok()?;
    let s = String::from_utf8_lossy(&out.stdout);
    // "seqkit v2.13.0" → extract the version number
    s.split_whitespace()
        .find(|w| w.starts_with('v'))
        .map(|v| v.trim_start_matches('v').to_owned())
}

/// Run our binary with the given args and return stdout bytes.
fn our_stdout(args: &[&str]) -> Vec<u8> {
    let out = bin()
        .args(args)
        .arg(fixture())
        .output()
        .expect("rsomics-fasta-locate failed to run");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

/// Run seqkit locate with the given args and return stdout bytes.
fn seqkit_stdout(args: &[&str]) -> Vec<u8> {
    let out = Command::new("seqkit")
        .args(["locate"])
        .args(args)
        .arg(fixture())
        .output()
        .expect("seqkit failed to run");
    assert!(out.status.success());
    out.stdout
}

// Byte-identical compat for default TSV output.
// seqkit v2.9 and v2.13 have identical locate output for plain literal patterns.
#[test]
fn compat_default_tsv() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ver = seqkit_version().unwrap_or_default();
    eprintln!("seqkit version: {ver}");

    let ours = our_stdout(&["-p", "ATGC"]);
    let theirs = seqkit_stdout(&["-p", "ATGC"]);
    assert_eq!(
        ours,
        theirs,
        "TSV mismatch.\nours:\n{}\ntheirs:\n{}",
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

#[test]
fn compat_positive_strand_only() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = our_stdout(&["-p", "ATGC", "-P"]);
    let theirs = seqkit_stdout(&["-p", "ATGC", "-P"]);
    assert_eq!(ours, theirs);
}

#[test]
fn compat_ignore_case() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = our_stdout(&["-p", "ATGC", "-i"]);
    let theirs = seqkit_stdout(&["-p", "ATGC", "-i"]);
    assert_eq!(
        ours,
        theirs,
        "ignore-case mismatch.\nours:\n{}\ntheirs:\n{}",
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

#[test]
fn compat_bed_output() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = our_stdout(&["-p", "ATGC", "--bed"]);
    let theirs = seqkit_stdout(&["-p", "ATGC", "--bed"]);
    assert_eq!(ours, theirs);
}

#[test]
fn compat_regex() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = our_stdout(&["-r", "-p", "ATG[CGT]"]);
    let theirs = seqkit_stdout(&["-r", "-p", "ATG[CGT]"]);
    assert_eq!(
        ours,
        theirs,
        "regex mismatch.\nours:\n{}\ntheirs:\n{}",
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

// Golden-file tests (no seqkit needed — CI-safe).
#[test]
fn golden_default() {
    let golden =
        std::fs::read(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small_atgc.tsv"))
            .expect("golden file missing");
    let ours = our_stdout(&["-p", "ATGC"]);
    assert_eq!(
        ours,
        golden,
        "golden mismatch.\nours:\n{}",
        String::from_utf8_lossy(&ours)
    );
}

#[test]
fn golden_positive_strand() {
    let golden = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small_atgc_plus.tsv"),
    )
    .expect("golden file missing");
    let ours = our_stdout(&["-p", "ATGC", "-P"]);
    assert_eq!(ours, golden);
}

#[test]
fn golden_ignore_case() {
    let golden = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small_atgc_icase.tsv"),
    )
    .expect("golden file missing");
    let ours = our_stdout(&["-p", "ATGC", "-i"]);
    assert_eq!(
        ours,
        golden,
        "ignore-case golden mismatch.\nours:\n{}",
        String::from_utf8_lossy(&ours)
    );
}

#[test]
fn golden_bed() {
    let golden =
        std::fs::read(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small_atgc.bed"))
            .expect("golden file missing");
    let ours = our_stdout(&["-p", "ATGC", "--bed"]);
    assert_eq!(ours, golden);
}
