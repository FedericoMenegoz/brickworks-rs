use crate::native::one_pole::{OnePoleCoeffs, OnePoleState, StickyMode};

pub struct Gain<const N_CHANNELS: usize> {
    coeffs: GainCoeffs<N_CHANNELS>,
}

impl<const N_CHANNELS: usize> Gain<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!();
    }
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!();
    }
    #[inline(always)]
    pub fn reset(&mut self) {
        todo!();
    }
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        todo!();
    }
    #[inline(always)]
    pub fn set_gain_lin(&mut self, value: f32) {
        todo!();
    }
    #[inline(always)]
    pub fn set_gain_db(&mut self, value: f32) {
        todo!();
    }
    #[inline(always)]
    pub fn set_smooth_tau(&mut self, value: f32) {
        todo!();
    }
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        todo!();
    }
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: StickyMode) {
        todo!();
    }
    #[inline(always)]
    pub fn get_gain_lin(&self) -> f32 {
        todo!();
    }
    #[inline(always)]
    pub fn get_gain_cur(&self) -> f32 {
        todo!();
    }
}

impl<const N_CHANNELS: usize> Default for Gain<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct GainCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    smooth_state: OnePoleState,

    // Parameters
    gain: f32,
}

impl<const N_CHANNELS: usize> GainCoeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_reset_coeffs(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_update_coeffs_ctrl(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_update_coeffs_audio(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_update_coeffs_audio_sticky_abs(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_update_coeffs_audio_sticky_rel(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_process1(&mut self, x: f32) -> f32 {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_process(&mut self, x: &[f32], y: &mut [f32], n_samples: usize) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_process_multi(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut f32; N_CHANNELS],
        n_samples: usize,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_gain_lin(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_gain_db(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_smooth_tau(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_sticky_thresh(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_set_sticky_mode(&mut self, value: StickyMode) {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_get_gain_lin(&self) -> f32 {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_get_gain_cur(&self) -> f32 {
        todo!()
    }

    #[inline(always)]
    pub fn bw_gain_coeffs_is_valid() -> bool {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for GainCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::gain::Gain as GainWrapper, native::one_pole::tests::assert_one_pole_coeffs,
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [&[1.0, 1.0], &[0.0, 0.0]];
    const N_SAMPLES: usize = 2;

    type GainT = Gain<N_CHANNELS>;
    type GainWrapperT = GainWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_gain = GainT::new();
        let c_gain = GainWrapperT::new();

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut rust_gain = GainT::new();
        let mut c_gain = GainWrapperT::new();

        rust_gain.set_sample_rate(SAMPLE_RATE);
        c_gain.set_sample_rate(SAMPLE_RATE);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be a positive number, got inf")]
    fn set_sample_rate_invalid() {
        let mut rust_gain = GainT::new();
        rust_gain.set_sample_rate(f32::INFINITY);
    }

    #[test]
    fn reset() {
        const TAU: f32 = 0.00005;
        let sticky_tresh = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();
        let x0 = [0.5; N_CHANNELS];

        c_gain.set_smooth_tau(TAU);
        c_gain.set_sticky_thresh(sticky_tresh);
        c_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        rust_gain.set_sticky_thresh(sticky_tresh);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        c_gain.reset();
        rust_gain.reset();

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn process() {
        const TAU: f32 = 0.00005;
        let sticky_tresh = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_smooth_tau(TAU);
        c_gain.set_sticky_thresh(sticky_tresh);
        c_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        rust_gain.set_sticky_thresh(sticky_tresh);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        c_gain.reset();
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);
        rust_gain.reset();
        rust_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_gain(&rust_gain, &c_gain);

        (0..N_CHANNELS).for_each(|channel| {
            (0..N_SAMPLES).for_each(|sample| {
                assert_eq!(rust_y[channel][sample], c_y[channel][sample]);
            });
        });
    }

    #[test]
    fn set_gain_lin_valid() {
        const GAIN_LIN: f32 = 0.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_LIN);
        c_gain.set_gain_lin(GAIN_LIN);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be ???, got inf")]
    fn set_gain_lin_invalid() {
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_gain_lin(f32::INFINITY);
    }

    #[test]
    fn set_gain_db_valid() {
        const GAIN_DB: f32 = 12.0;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_DB);
        c_gain.set_gain_lin(GAIN_DB);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be ???, got 770.7")]
    fn set_gain_db_invalid() {
        const GAIN_DB: f32 = 770.7;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_gain_db(GAIN_DB);
    }

    #[test]
    fn set_smooth_tau_valid() {
        const TAU: f32 = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        c_gain.set_smooth_tau(TAU);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_smooth_tau_invalid() {
        const TAU: f32 = -1.0;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_gain_db(TAU);
    }

    #[test]
    fn set_sticky_thresh_valid() {
        const STICKY_THRESH: f32 = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_sticky_thresh(STICKY_THRESH);
        c_gain.set_sticky_thresh(STICKY_THRESH);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got -1")]
    fn set_sticky_thresh_invalid() {
        const STICKY_THRESH: f32 = -1.0;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_sticky_thresh(STICKY_THRESH);
    }

    #[test]
    fn set_sticky_mode() {
        const STICKY_MODE: StickyMode = StickyMode::Rel;
        const STICKY_THRESH: f32 = 0.01;

        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_sticky_thresh(STICKY_THRESH);
        c_gain.set_sticky_thresh(STICKY_THRESH);

        rust_gain.set_sticky_mode(STICKY_MODE);
        c_gain.set_sticky_mode(STICKY_MODE as u32);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn get_gain_lin() {
        const GAIN_DB: f32 = 3.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_DB);
        c_gain.set_gain_lin(GAIN_DB);

        rust_gain.reset();
        c_gain.reset();

        rust_gain.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        rust_gain.reset();
        c_gain.reset();

        assert_eq!(&rust_gain.get_gain_lin(), &c_gain.get_gain_lin());
        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn get_gain_cur() {
        const GAIN_LIN: f32 = 0.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_LIN);
        c_gain.set_gain_lin(GAIN_LIN);

        rust_gain.reset();
        c_gain.reset();

        rust_gain.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        rust_gain.reset();
        c_gain.reset();

        assert_eq!(&rust_gain.get_gain_cur(), &c_gain.get_gain_cur());
        assert_gain(&rust_gain, &c_gain);
    }

    fn assert_gain<const N_CHANNELS: usize>(
        rust_gain: &Gain<N_CHANNELS>,
        c_gain: &GainWrapper<N_CHANNELS>,
    ) {
        assert_one_pole_coeffs::<N_CHANNELS>(
            &rust_gain.coeffs.smooth_coeffs,
            &c_gain.coeffs.smooth_coeffs,
        );
        assert_eq!(
            rust_gain.coeffs.smooth_state.get_y_z1(),
            c_gain.coeffs.smooth_state.y_z1
        );
        assert_eq!(rust_gain.coeffs.gain, c_gain.coeffs.gain);
    }
}
