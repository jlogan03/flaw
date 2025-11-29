use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use flaw::sos;

pub fn criterion_benchmark(c: &mut Criterion) {
    let cutoff_ratio = 0.05;

    let mut butter4_f64 = sos::butter4::<f64>(cutoff_ratio).unwrap();
    c.bench_function("sos butter4 f64", |b| b.iter(|| {black_box(butter4_f64.update(1.0))} ));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
