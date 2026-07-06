//! Criterion benchmark suite.
//!
//! At M0 this benches the vle-thermo smoke path only, so `cargo bench` and
//! `cargo test --all-targets` (which compiles benches) work from the start.
//! The real suite — block-Thomas, one Naphtali–Sandholm iteration, a full
//! debutanizer solve — lands from M6 on (PLAN §9), regression-guarded like
//! vle.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use stages_thermo::thermo::smoke_bubble_temperature;

fn bench_smoke_bubble_temperature(c: &mut Criterion) {
    c.bench_function("smoke_bubble_temperature", |b| {
        b.iter(|| black_box(smoke_bubble_temperature().unwrap()));
    });
}

criterion_group!(benches, bench_smoke_bubble_temperature);
criterion_main!(benches);
