use brickworks_rs::c_wrapper::dist::Dist as DistC;
use brickworks_rs::native::dist::Dist;
use criterion::{Criterion, criterion_group, criterion_main};

const N_CHANNELS: usize = 2;
const N_SAMPLES: usize = 512;
const SAMPLE_RATE: f32 = 48_000.0;

fn benchmarks(c: &mut Criterion) {
    let mut sine = [0.0; N_SAMPLES];

    let freq = 1000.0;
    (0..N_SAMPLES).for_each(|sample| {
        let t = sample as f32 / SAMPLE_RATE;
        sine[sample] = (2.0 * std::f32::consts::PI * freq * t).sin();
    });

    let input_x: [&[f32]; N_CHANNELS] = [&sine, &sine.clone()];
    let y_ch0: [f32; N_SAMPLES] = [0.0; N_SAMPLES];
    let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut y_ch0.clone(), &mut y_ch0.clone()];
    let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut y_ch0.clone(), &mut y_ch0.clone()];

    let mut group = c.benchmark_group("Dist Test");

    group.bench_function("Rust Dist", |b| {
        b.iter(|| {
            let mut dist = Dist::<N_CHANNELS>::new();
            dist.set_sample_rate(SAMPLE_RATE);
            dist.set_distortion(0.5);
            dist.set_volume(0.8);
            dist.set_tone(0.6);

            dist.reset(None, None);
            dist.process(&input_x, &mut rust_y, N_SAMPLES);
        })
    });

    group.bench_function("C Wrapper Dist", |b| {
        b.iter(|| {
            let mut dist = DistC::<N_CHANNELS>::new();
            dist.set_sample_rate(SAMPLE_RATE);
            dist.set_distortion(0.5);
            dist.set_volume(0.8);
            dist.set_tone(0.6);

            dist.reset(None, None);
            dist.process(&input_x, &mut c_y, N_SAMPLES);
        })
    });
    group.finish();
}

criterion_group!(benches, benchmarks);
criterion_main!(benches);
