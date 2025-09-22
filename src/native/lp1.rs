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
