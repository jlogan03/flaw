use criterion::{Criterion, criterion_group, criterion_main};
use flaw::{butter2, polynomial_fractional_delay, sos};
use std::hint::black_box;

pub fn criterion_benchmark(c: &mut Criterion) {
    let cutoff_ratio = 0.05;

    let mut butter4_f64 = sos::butter4::<f64>(cutoff_ratio).unwrap();
    c.bench_function("sos butter4 f64", |b| {
        let input = black_box(1.0);
        b.iter(|| black_box(butter4_f64.update(input)))
    });

    let mut butter2_f32 = butter2(cutoff_ratio).unwrap();
    c.bench_function("statespace butter2 f32", |b| {
        let input = black_box(1.0);
        b.iter(|| black_box(butter2_f32.update(input)))
    });

    let mut butter2_bank = [butter2_f32; 20];
    c.bench_function("statespace butter2 f32 bank 20x", |b| {
        let input = black_box([1.0; 20]);
        b.iter(|| {
            black_box(
                butter2_bank
                    .iter_mut()
                    .zip(input.iter())
                    .for_each(|(f, v)| {
                        f.update(*v);
                    }),
            )
        })
    });

    let mut frac_delay3_f32 = polynomial_fractional_delay::<3, f32>(0.5);
    c.bench_function("fractional delay order3 f32", |b| {
        b.iter(|| {
            let input = black_box(1.0);
            black_box(frac_delay3_f32.update(input))
        })
    });

    let mut frac_delay3_bank = [frac_delay3_f32; 20];
    c.bench_function("fractional delay order3 f32 bank 20x", |b| {
        let input = black_box([1.0; 20]);
        b.iter(|| {
            black_box(
                frac_delay3_bank
                    .iter_mut()
                    .zip(input.iter())
                    .for_each(|(f, v)| {
                        f.update(*v);
                    }),
            )
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
