use crate::native::{
    one_pole::{OnePoleCoeffs, OnePoleState},
    satur::SaturState,
};

pub struct LP1<const N_CHANNELS: usize> {
    coeffs: LP1Coeffs<N_CHANNELS>,
    states: [LP1State; N_CHANNELS],
}

impl<const N_CHANNELS: usize> LP1<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for LP1<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct LP1Coeffs<const N_CHANNELS: usize> {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    smooth_cutoff_state: OnePoleState,
    smooth_prewarp_freq_state: OnePoleState,

    // Coefficients
    t_k: f32,

    t: f32,
    x_x: f32,
    x_x_z1: f32,
    y_x: f32,

    // Parameters
    cutoff: f32,
    prewarp_k: f32,
    prewarp_freq: f32,
}

impl<const N_CHANNELS: usize> LP1Coeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut LP1State, x0: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [SaturState; N_CHANNELS],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn process1(&mut self, state: &mut LP1State, x: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn process(&mut self, state: &mut LP1State, x: &[f32], y: &mut [f32], n_samples: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        state: &mut [LP1State],
        x: &[&[f32]],
        y: &mut [&mut [f32]],
        n_samples: f32,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, vaulue: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, vaulue: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, vaulue: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn coeffs_is_valid(&mut self) -> bool {
        todo!()
    }

    #[inline(always)]
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for LP1Coeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

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
mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_lp1_coeffs as LP1CoeffsWrapper, lp1::LP1 as LP1Wrapper},
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

    #[should_panic(expected = "value must positive, got -0.1")]
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

    fn assert_lp1_coeffs<const N_CHANNELS: usize>(
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
