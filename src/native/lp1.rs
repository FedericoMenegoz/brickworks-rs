use core::f32;

use crate::native::{
    math::{minf, rcpf, tanf},
    one_pole::{OnePoleCoeffs, OnePoleState},
};

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};

pub struct LP1<const N_CHANNELS: usize> {
    coeffs: LP1Coeffs<N_CHANNELS>,
    states: [LP1State; N_CHANNELS],
}

impl<const N_CHANNELS: usize> LP1<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let coeffs = LP1Coeffs::new();
        Self {
            coeffs,
            states: [LP1State::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.coeffs.reset_state(&mut self.states[channel], x0);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.coeffs.reset_state(&mut self.states[channel], x0);
            }),
        }
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        self.coeffs.reset_state_multi(&mut self.states, x0, y0);
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.coeffs.process_multi(&mut self.states, x, y, n_samples);
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        self.coeffs.set_cutoff(value);
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.coeffs.set_prewarp_at_cutoff(value);
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        self.coeffs.set_prewarp_freq(value);
    }
}

impl<const N_CHANNELS: usize> Default for LP1<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct LP1Coeffs<const N_CHANNELS: usize> {
    // Sub-components
    pub(crate) smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    pub(crate) smooth_cutoff_state: OnePoleState,
    pub(crate) smooth_prewarp_freq_state: OnePoleState,

    // Coefficients
    pub(crate) t_k: f32,

    pub(crate) t: f32,
    pub(crate) x_x: f32,
    pub(crate) x_x_z1: f32,
    pub(crate) y_x: f32,

    // Parameters
    pub(crate) cutoff: f32,
    pub(crate) prewarp_k: f32,
    pub(crate) prewarp_freq: f32,
}

impl<const N_CHANNELS: usize> LP1Coeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut smooth_coeffs = OnePoleCoeffs::<N_CHANNELS>::new();
        smooth_coeffs.set_tau(0.005);
        smooth_coeffs.set_sticky_thresh(1e-3);
        let cutoff = 1e3;
        let prewarp_k = 1.0;
        let prewarp_freq = 1e3;

        Self {
            smooth_coeffs,
            smooth_cutoff_state: OnePoleState::new(),
            smooth_prewarp_freq_state: OnePoleState::new(),
            t_k: 0.0,
            t: 0.0,
            x_x: 0.0,
            x_x_z1: 0.0,
            y_x: 0.0,
            cutoff,
            prewarp_k,
            prewarp_freq,
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_positive(sample_rate);
            debug_assert_is_finite(sample_rate);
        }

        self.smooth_coeffs.set_sample_rate(sample_rate);
        self.smooth_coeffs.reset_coeffs();

        self.t_k = f32::consts::PI / sample_rate;
    }

    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.smooth_coeffs
            .reset_state(&mut self.smooth_cutoff_state, self.cutoff);
        self.smooth_coeffs.reset_state(
            &mut self.smooth_prewarp_freq_state,
            self.prewarp_freq + self.prewarp_k * (self.cutoff - self.prewarp_freq),
        );
        self.do_update_coeffs(true);
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut LP1State, x0: f32) -> f32 {
        debug_assert!(x0.is_finite(), "value must be finite, got {}", x0);

        let y = x0;
        state.y_z1 = x0;
        state.x_z1 = 0.0;

        debug_assert!(y.is_finite(), "value must be finite, got {}", y);

        y
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [LP1State; N_CHANNELS],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.reset_state(&mut states[channel], x0[channel]);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut states[channel], x0[channel]);
            }),
        }
    }

    // Not implemented yet: C version only contained assertions
    // need to revisit which assertions from the C version make sense to keep in Rust
    // #[inline(always)]
    // pub fn update_coeffs_ctrl(&mut self) {
    //     todo!()
    // }

    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.do_update_coeffs(false);
    }

    #[inline(always)]
    pub fn process1(&mut self, state: &mut LP1State, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let x_b = self.x_x * (x - state.y_z1) - self.x_x_z1 * state.x_z1;
        let y = x - self.y_x * x_b;
        state.y_z1 = y;
        state.x_z1 = x_b;

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn process(&mut self, state: &mut LP1State, x: &[f32], y: &mut [f32], n_samples: usize) {
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample]);
        });
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [LP1State],
        x: &[&[f32]],
        y: &mut [&mut [f32]],
        n_samples: usize,
    ) {
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            (0..N_CHANNELS).for_each(|channel| {
                y[channel][sample] = self.process1(&mut states[channel], x[channel][sample])
            });
        });
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite());
            debug_assert_range(1e-6..=1e12, value);
        }

        self.cutoff = value;
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.prewarp_k = if value { 1.0 } else { 0.0 };
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(
                value.is_finite(),
                "value must be in range [1e-6, 1e12], got inf"
            );
            debug_assert_range(1e-6..=1e12, value);
        }

        self.prewarp_freq = value;
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    #[inline(always)]
    pub fn coeffs_is_valid(&mut self) -> bool {
        todo!()
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    #[inline(always)]
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }

    // Private
    #[inline(always)]
    fn do_update_coeffs(&mut self, force: bool) {
        let prewarp_freq = self.prewarp_freq + self.prewarp_k * (self.cutoff - self.prewarp_freq);
        let mut prewarp_freq_cur = self.smooth_prewarp_freq_state.get_y_z1();
        let mut cutoff_cur = self.smooth_cutoff_state.get_y_z1();
        let prewarp_freq_changed = force || prewarp_freq != prewarp_freq_cur;
        let cutoff_changed = force || self.cutoff != cutoff_cur;

        if prewarp_freq_changed || cutoff_changed {
            if prewarp_freq_changed {
                prewarp_freq_cur = self
                    .smooth_coeffs
                    .process1_sticky_rel(&mut self.smooth_prewarp_freq_state, prewarp_freq);
                self.t = tanf(minf(
                    self.t_k * prewarp_freq_cur,
                    1.567_654_7, /*34141306*/
                )); // max = 0.499 * fs
            }
            if cutoff_changed {
                cutoff_cur = self
                    .smooth_coeffs
                    .process1_sticky_rel(&mut self.smooth_cutoff_state, self.cutoff);
                self.y_x = rcpf(cutoff_cur);
            }
            let k = cutoff_cur * rcpf(cutoff_cur * self.t + prewarp_freq_cur);
            self.x_x = k * prewarp_freq_cur;
            self.x_x_z1 = k * self.t;
        }
    }
}

impl<const N_CHANNELS: usize> Default for LP1Coeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
#[derive(Default, Debug, Clone, Copy)]
pub struct LP1State {
    y_z1: f32,
    x_z1: f32,
}

impl LP1State {
    pub fn get_y_z1(&self) -> f32 {
        self.y_z1
    }

    pub fn get_x_z1(&self) -> f32 {
        self.x_z1
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_lp1_coeffs as LP1CoeffsWrapper, bw_lp1_state, lp1::LP1 as LP1Wrapper},
        native::one_pole::tests::assert_one_pole_coeffs,
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type LP1T = LP1<N_CHANNELS>;
    type LP1WrapperT = LP1Wrapper<N_CHANNELS>;

    #[test]
    pub fn new() {
        let rust_lp1 = LP1T::new();
        let c_lp1 = LP1WrapperT::new();

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[test]
    pub fn set_sample_rate() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[should_panic(expected = "value must be non negative, got -0.1")]
    #[test]
    pub fn set_sample_rate_invalid() {
        let mut rust_lp1 = LP1T::new();

        rust_lp1.set_sample_rate(-0.1);
    }

    #[test]
    pub fn reset_none() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        rust_lp1.reset(0.0, None);
        c_lp1.reset(0.0, None);

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[test]
    pub fn reset_some() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_lp1.reset(1.0, Some(&mut rust_y0));
        c_lp1.reset(1.0, Some(&mut c_y0));

        assert_lp1(&rust_lp1, &c_lp1);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    pub fn reset_multi() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.5, 0.5];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_lp1.reset_multi(&x0, Some(&mut rust_y0));
        c_lp1.reset_multi(&x0, Some(&mut c_y0));

        assert_lp1(&rust_lp1, &c_lp1);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    pub fn process1() {
        let mut rust_coeffs = LP1Coeffs::<N_CHANNELS>::new();
        let mut c_coeffs = LP1CoeffsWrapper::new();
        let mut rust_state = LP1State::default();
        let mut c_state = bw_lp1_state::default();

        rust_coeffs.set_sample_rate(SAMPLE_RATE);
        c_coeffs.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.5;

        rust_coeffs.reset_coeffs();
        c_coeffs.reset_coeffs();

        rust_coeffs.reset_state(&mut rust_state, x0);
        c_coeffs.reset_state(&mut c_state, x0);

        let rust_y = rust_coeffs.process1(&mut rust_state, x0);
        let c_y = c_coeffs.process1(&mut c_state, x0);
        println!("{rust_y} and {c_y}");
        assert_lp1_coeffs(&rust_coeffs, &c_coeffs);
        assert_eq!(rust_y, c_y);
        assert_eq!(rust_state.x_z1, c_state.X_z1);
        assert_eq!(rust_state.y_z1, c_state.y_z1);
    }

    #[test]
    pub fn process() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 756.21;
        let prewarp_freq = 780.00;

        rust_lp1.set_cutoff(cutoff);
        rust_lp1.set_prewarp_at_cutoff(true);
        rust_lp1.set_prewarp_freq(prewarp_freq);

        c_lp1.set_cutoff(cutoff);
        c_lp1.set_prewarp_at_cutoff(true);
        c_lp1.set_prewarp_freq(prewarp_freq);

        rust_lp1.reset(0.0, None);
        c_lp1.reset(0.0, None);

        let y_ch: Box<dyn Fn() -> [f32; 8]> = Box::new(|| std::array::from_fn(|_| 0.0));

        let mut rust_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];
        let mut c_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];

        rust_lp1.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_lp1.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_lp1(&rust_lp1, &c_lp1);
        assert_eq!(rust_y, c_y);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_lp1.reset(0.0, Some(&mut rust_y0));
        c_lp1.reset(0.0, Some(&mut c_y0));

        assert_lp1(&rust_lp1, &c_lp1);
        assert_eq!(rust_y0, c_y0);

        rust_lp1.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_lp1.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);
        assert_lp1(&rust_lp1, &c_lp1);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    pub fn set_cutoff() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 756.21;

        rust_lp1.set_cutoff(cutoff);
        c_lp1.set_cutoff(cutoff);

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 0")]
    #[test]
    pub fn set_cutoff_invalid() {
        let mut rust_lp1 = LP1T::new();

        let cutoff = 0.0;

        rust_lp1.set_cutoff(cutoff);
    }

    #[test]
    pub fn set_prewarp_at_cutoff() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        rust_lp1.set_prewarp_at_cutoff(true);
        c_lp1.set_prewarp_at_cutoff(true);

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[test]
    pub fn set_prewarp_freq() {
        let mut rust_lp1 = LP1T::new();
        let mut c_lp1 = LP1WrapperT::new();

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        c_lp1.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 1000.0;

        rust_lp1.set_prewarp_at_cutoff(true);
        c_lp1.set_prewarp_at_cutoff(true);

        rust_lp1.set_prewarp_freq(prewarp_freq);
        c_lp1.set_prewarp_freq(prewarp_freq);

        assert_lp1(&rust_lp1, &c_lp1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got inf")]
    #[test]
    pub fn set_prewarp_freq_invalid() {
        let mut rust_lp1 = LP1T::new();
        let prewarp_freq = f32::INFINITY;

        rust_lp1.set_sample_rate(SAMPLE_RATE);
        rust_lp1.set_prewarp_at_cutoff(true);
        rust_lp1.set_prewarp_freq(prewarp_freq);
    }

    fn assert_lp1<const N_CHANNELS: usize>(
        rust_lp1: &LP1<N_CHANNELS>,
        c_lp1: &LP1Wrapper<N_CHANNELS>,
    ) {
        assert_lp1_coeffs(&rust_lp1.coeffs, &c_lp1.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(
                rust_lp1.states[channel].get_x_z1(),
                c_lp1.states[channel].X_z1
            );
            assert_eq!(
                rust_lp1.states[channel].get_y_z1(),
                c_lp1.states[channel].y_z1
            );
        });
    }

    pub(crate) fn assert_lp1_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &LP1Coeffs<N_CHANNELS>,
        c_coeffs: &LP1CoeffsWrapper,
    ) {
        assert_eq!(rust_coeffs.t_k, c_coeffs.t_k);
        assert_eq!(rust_coeffs.t, c_coeffs.t);
        assert_eq!(rust_coeffs.x_x, c_coeffs.X_x);
        assert_eq!(rust_coeffs.x_x_z1, c_coeffs.X_X_z1);
        assert_eq!(rust_coeffs.y_x, c_coeffs.y_X);
        assert_eq!(rust_coeffs.cutoff, c_coeffs.cutoff);
        assert_eq!(rust_coeffs.prewarp_k, c_coeffs.prewarp_k);
        assert_eq!(rust_coeffs.prewarp_freq, c_coeffs.prewarp_freq);

        assert_one_pole_coeffs(&rust_coeffs.smooth_coeffs, &c_coeffs.smooth_coeffs);
        assert_eq!(
            rust_coeffs.smooth_cutoff_state.get_y_z1(),
            c_coeffs.smooth_cutoff_state.y_z1
        );
        assert_eq!(
            rust_coeffs.smooth_prewarp_freq_state.get_y_z1(),
            c_coeffs.smooth_prewarp_freq_state.y_z1
        );
    }
}
