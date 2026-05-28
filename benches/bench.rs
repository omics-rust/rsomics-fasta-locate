use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_fasta_locate(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fasta-locate");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fa = manifest.join("tests/golden/small.fa");
    c.bench_function("rsomics-fasta-locate golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([fa.to_str().unwrap(), "ATG"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fasta_locate);
criterion_main!(benches);
