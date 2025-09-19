use crate::native::one_pole::{OnePoleCoeffs, OnePoleState};

pub struct Satur<const N_CHANNELS: usize> {
    pub(crate) coeffs: SaturCoeffs<N_CHANNELS>,
    pub(crate) states: [SaturState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Satur<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32]>) {
        todo!()
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }

    #[inline(always)]
    pub fn process(&mut self, x: &[&[f32]; N_CHANNELS], y: &mut [&mut [f32]], n_samples: usize) {
        todo!()
    }

    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for Satur<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct SaturCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    smooth_bias_state: OnePoleState,
    smooth_gain_state: OnePoleState,

    // Coefficients
    bias_dc: f32,
    inv_gain: f32,

    // Parameters
    bias: f32,
    gain: f32,
    gain_compensation: bool,
}

impl<const N_CHANNELS: usize> SaturCoeffs<N_CHANNELS> {
    #[inline(always)]
    fn tanhf() -> f32 {
        todo!()
    }
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }
    #[inline(always)]
    pub fn set_sample_rate() {
        todo!()
    }
    #[inline(always)]
    pub fn do_update_coeffs() {
        todo!()
    }
    #[inline(always)]
    pub fn reset_coeffs() {
        todo!()
    }
    #[inline(always)]
    pub fn reset_state() -> f32 {
        todo!()
    }
    #[inline(always)]
    pub fn reset_state_multi() {
        todo!()
    }
    #[inline(always)]
    pub fn update_coeffs_ctrl() {
        todo!()
    }
    #[inline(always)]
    pub fn update_coeffs_audio() {
        todo!()
    }
    #[inline(always)]
    pub fn process1() -> f32 {
        todo!()
    }
    #[inline(always)]
    pub fn process1_comp() -> f32 {
        todo!()
    }
    #[inline(always)]
    pub fn process() {
        todo!()
    }
    #[inline(always)]
    pub fn process_multi() {
        todo!()
    }
    #[inline(always)]
    pub fn set_bias() {
        todo!()
    }
    #[inline(always)]
    pub fn set_gain() {
        todo!()
    }
    #[inline(always)]
    pub fn set_gain_compensation() {
        todo!()
    }
    #[inline(always)]
    pub fn coeffs_is_valid() -> bool {
        todo!()
    }
    #[inline(always)]
    pub fn state_is_valid() -> bool {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for SaturCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct SaturState {
    x_z1: f32,
    f_z1: f32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        c_wrapper::{bw_satur_coeffs, satur::Satur as SaturWrapper},
        native::one_pole::tests::assert_one_pole_coeffs,
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type SaturT = Satur<N_CHANNELS>;
    type SaturWrapperT = SaturWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn set_sample_rate() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn reset() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn reset_multi() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn process() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn set_bias() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn set_gain() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    #[test]
    fn set_gain_compensation() {
        let rust_satur = SaturT::new();
        let c_satur = SaturWrapperT::new();
    }

    fn assert_satur_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &SaturCoeffs<N_CHANNELS>,
        c_coeffs: &bw_satur_coeffs,
    ) {
        // Coefficients
        assert_eq!(rust_coeffs.bias_dc, c_coeffs.bias_dc);
        assert_eq!(rust_coeffs.inv_gain, c_coeffs.inv_gain);

        // Parameters
        assert_eq!(rust_coeffs.bias, c_coeffs.bias);
        assert_eq!(rust_coeffs.gain, c_coeffs.gain);
        assert_eq!(
            rust_coeffs.gain_compensation,
            c_coeffs.gain_compensation != 0
        );
    }

    fn assert_satur<const N_CHANNELS: usize>(
        rust_satur: &Satur<N_CHANNELS>,
        c_satur: &SaturWrapper<N_CHANNELS>,
    ) {
        // Sub-components
        assert_one_pole_coeffs(
            &rust_satur.coeffs.smooth_coeffs,
            &c_satur.coeffs.smooth_coeffs,
        );
        assert_eq!(
            &rust_satur.coeffs.smooth_bias_state.get_y_z1(),
            &c_satur.coeffs.smooth_bias_state.y_z1
        );
        assert_eq!(
            &rust_satur.coeffs.smooth_gain_state.get_y_z1(),
            &c_satur.coeffs.smooth_gain_state.y_z1
        );
    }
}
