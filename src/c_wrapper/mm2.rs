use super::*;

pub struct MM2<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_mm2_coeffs,
    pub(crate) states: [bw_mm2_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> MM2<N_CHANNELS> {
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
    pub fn set_q(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_coeff_x(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_coeff_lp(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_coeff_bp(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_coeff_hp(&mut self, value: f32) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::math::INVERSE_2_PI;
    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [&[1.0, 1.0], &[0.0, 0.0]];
    const N_SAMPLES: usize = 2;

    type MM2T = MM2<N_CHANNELS>;

    #[test]
    fn new() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let tau_default = 0.005;
        let cutoff;
        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(mm2.coeffs.gain_x_coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(mm2.coeffs.gain_lp_coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(mm2.coeffs.gain_bp_coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(mm2.coeffs.gain_hp_coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(mm2.coeffs.gain_x_coeffs.gain, 1.);
        assert_eq!(mm2.coeffs.gain_lp_coeffs.gain, 0.);
        assert_eq!(mm2.coeffs.gain_bp_coeffs.gain, 0.);
        assert_eq!(mm2.coeffs.gain_hp_coeffs.gain, 0.);
    }

    #[test]
    fn set_sample_rate() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        assert_eq!(
            mm2.coeffs.svf_coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
        assert_eq!(
            mm2.coeffs.gain_x_coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
        assert_eq!(
            mm2.coeffs.gain_lp_coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
        assert_eq!(
            mm2.coeffs.gain_bp_coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
        assert_eq!(
            mm2.coeffs.gain_hp_coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
    }

    #[test]
    fn reset_none() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
        let x0 = 0.0;
        mm2.reset(x0, None);
    }

    #[test]
    fn reset_some() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn reset_multi() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn process() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_cutoff() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_q() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_coeff_x() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_coeff_lp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_coeff_bp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }

    #[test]
    fn set_coeff_hp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
    }
}
