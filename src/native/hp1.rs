use crate::native::lp1::{LP1Coeffs, LP1State};

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
pub struct HP1<const N_CHANNELS: usize> {
    coeffs: HP1Coeffs<N_CHANNELS>,
    states: [HP1State; N_CHANNELS],
}

impl<const N_CHANNELS: usize> HP1<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let coeffs = HP1Coeffs::<N_CHANNELS>::new();
        Self {
            coeffs,
            states: [HP1State::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32]>) {
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
    pub fn process(&mut self, x: &[&[f32]], y: &mut [&mut [f32]], n_samples: usize) {
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

impl<const N_CHANNELS: usize> Default for HP1<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
pub struct HP1Coeffs<const N_CHANNELS: usize> {
    //Sub-component
    lp1_coeffs: LP1Coeffs<N_CHANNELS>,
}

impl<const N_CHANNELS: usize> HP1Coeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let lp1_coeffs = LP1Coeffs::<N_CHANNELS>::new();
        Self { lp1_coeffs }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_positive(sample_rate);
            debug_assert_is_finite(sample_rate);
        }

        self.lp1_coeffs.set_sample_rate(sample_rate);
    }

    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.lp1_coeffs.reset_coeffs();
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut HP1State, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        let lp = self.lp1_coeffs.reset_state(&mut state.lp1_state, x0);
        let y = x0 - lp;

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [HP1State; N_CHANNELS],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut states[channel], x0[channel]);
                y[channel] = 0.0;
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

    // }

    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.lp1_coeffs.update_coeffs_audio();
    }

    #[inline(always)]
    pub fn process1(&mut self, state: &mut HP1State, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let lp = self.lp1_coeffs.process1(&mut state.lp1_state, x);
        let y = x - lp;

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn process(&mut self, state: &mut HP1State, x: &[f32], y: &mut [f32], n_samples: usize) {
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample]);
        });
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [HP1State],
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

        self.lp1_coeffs.set_cutoff(value);
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.lp1_coeffs.set_prewarp_at_cutoff(value);
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

        self.lp1_coeffs.set_prewarp_freq(value);
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
}

impl<const N_CHANNELS: usize> Default for HP1Coeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Default, Clone, Copy, Debug, PartialEq)]
pub struct HP1State {
    // Sub-components
    lp1_state: LP1State,
}

#[cfg(test)]
pub(crate) mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_hp1_coeffs, bw_hp1_state, hp1::HP1 as HP1Wrapper},
        native::lp1::tests::assert_lp1_coeffs,
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type HP1T = HP1<N_CHANNELS>;
    type HP1WrapperT = HP1Wrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_hp1 = HP1T::new();
        let c_hp1 = HP1WrapperT::new();

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[should_panic(expected = "value must be non negative, got -0.1")]
    #[test]
    fn set_sample_rate_invalid() {
        let mut rust_hp1 = HP1T::new();

        rust_hp1.set_sample_rate(-0.1);
    }

    #[test]
    fn reset_none() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        rust_hp1.reset(0.0, None);
        c_hp1.reset(0.0, None);

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[test]
    fn reset_some() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_hp1.reset(1.0, Some(&mut rust_y0));
        c_hp1.reset(1.0, Some(&mut c_y0));

        assert_hp1(&rust_hp1, &c_hp1);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn reset_multi() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.5, 0.5];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_hp1.reset_multi(&x0, Some(&mut rust_y0));
        c_hp1.reset_multi(&x0, Some(&mut c_y0));

        assert_hp1(&rust_hp1, &c_hp1);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 756.21;
        let prewarp_freq = 780.00;

        rust_hp1.set_cutoff(cutoff);
        rust_hp1.set_prewarp_at_cutoff(true);
        rust_hp1.set_prewarp_freq(prewarp_freq);

        c_hp1.set_cutoff(cutoff);
        c_hp1.set_prewarp_at_cutoff(true);
        c_hp1.set_prewarp_freq(prewarp_freq);

        rust_hp1.reset(0.0, None);
        c_hp1.reset(0.0, None);

        let y_ch: Box<dyn Fn() -> [f32; 8]> = Box::new(|| std::array::from_fn(|_| 0.0));

        let mut rust_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];
        let mut c_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];

        rust_hp1.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_hp1.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_hp1(&rust_hp1, &c_hp1);
        assert_eq!(rust_y, c_y);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_hp1.reset(0.0, Some(&mut rust_y0));
        c_hp1.reset(0.0, Some(&mut c_y0));

        assert_hp1(&rust_hp1, &c_hp1);
        assert_eq!(rust_y0, c_y0);

        rust_hp1.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_hp1.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);
        assert_hp1(&rust_hp1, &c_hp1);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    fn set_cutoff() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 756.21;

        rust_hp1.set_cutoff(cutoff);
        c_hp1.set_cutoff(cutoff);

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 0")]
    #[test]
    fn set_cutoff_invalid() {
        let mut rust_hp1 = HP1T::new();

        let cutoff = 0.0;

        rust_hp1.set_cutoff(cutoff);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        rust_hp1.set_prewarp_at_cutoff(true);
        c_hp1.set_prewarp_at_cutoff(true);

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut rust_hp1 = HP1T::new();
        let mut c_hp1 = HP1WrapperT::new();

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        c_hp1.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 1000.0;

        rust_hp1.set_prewarp_at_cutoff(true);
        c_hp1.set_prewarp_at_cutoff(true);

        rust_hp1.set_prewarp_freq(prewarp_freq);
        c_hp1.set_prewarp_freq(prewarp_freq);

        assert_hp1(&rust_hp1, &c_hp1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got inf")]
    #[test]
    fn set_prewarp_freq_invalid() {
        let mut rust_hp1 = HP1T::new();
        let prewarp_freq = f32::INFINITY;

        rust_hp1.set_sample_rate(SAMPLE_RATE);
        rust_hp1.set_prewarp_at_cutoff(true);
        rust_hp1.set_prewarp_freq(prewarp_freq);
    }

    fn assert_hp1<const N_CHANNELS: usize>(
        rust_hp1: &HP1<N_CHANNELS>,
        c_hp1: &HP1Wrapper<N_CHANNELS>,
    ) {
        assert_hp1_coeffs(&rust_hp1.coeffs, &c_hp1.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_hp1_state(&rust_hp1.states[channel], &c_hp1.states[channel]);
        });
    }

    pub(crate) fn assert_hp1_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &HP1Coeffs<N_CHANNELS>,
        c_coeffs: &bw_hp1_coeffs,
    ) {
        assert_lp1_coeffs(&rust_coeffs.lp1_coeffs, &c_coeffs.lp1_coeffs);
    }

    pub(crate) fn assert_hp1_state(rust_state: &HP1State, c_state: &bw_hp1_state) {
        assert_eq!(rust_state.lp1_state.get_x_z1(), c_state.lp1_state.X_z1);
        assert_eq!(rust_state.lp1_state.get_y_z1(), c_state.lp1_state.y_z1);
    }
}
