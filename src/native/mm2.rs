pub struct MM2<const N_CHANNELS: usize> {
    pub(crate) coeffs: MM2Coeffs<N_CHANNELS>,
    pub(crate) states: [MM2State; N_CHANNELS],
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

impl<const N_CHANNELS: usize> Default for MM2<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct MM2Coeffs<const N_CHANNELS: usize> {}

pub struct MM2State {}
