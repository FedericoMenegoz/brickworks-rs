use std::ptr::null_mut;

use super::*;

pub struct MM2<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_mm2_coeffs,
    pub(crate) states: [bw_mm2_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> MM2<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_mm2_coeffs::default();
        unsafe {
            bw_mm2_init(&mut coeffs);
        }
        Self {
            coeffs,
            states: [bw_mm2_state::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_mm2_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_mm2_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(val) => (0..N_CHANNELS).for_each(|channel| {
                    val[channel] =
                        bw_mm2_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_mm2_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
            }
        }
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_mm2_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_mm2_reset_coeffs(&mut self.coeffs);
            bw_mm2_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_mm2_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_mm2_process_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_cutoff(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_Q(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        unsafe {
            bw_mm2_set_prewarp_at_cutoff(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_prewarp_freq(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_coeff_x(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_coeff_x(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_coeff_lp(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_coeff_lp(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_coeff_bp(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_coeff_bp(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_coeff_hp(&mut self, value: f32) {
        unsafe {
            bw_mm2_set_coeff_hp(&mut self.coeffs, value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::math::INVERSE_2_PI;
    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    const N_SAMPLES: usize = 8;

    type MM2T = MM2<N_CHANNELS>;

    #[test]
    fn new() {
        let mm2 = MM2T::new();

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
        assert_eq!(
            mm2.coeffs.state,
            bw_mm2_coeffs_state_bw_mm2_coeffs_state_init
        );
    }

    #[test]
    fn set_sample_rate() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
        let fs_2pi = INVERSE_2_PI * SAMPLE_RATE;

        assert_eq!(mm2.coeffs.svf_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(mm2.coeffs.gain_x_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(mm2.coeffs.gain_lp_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(mm2.coeffs.gain_bp_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(mm2.coeffs.gain_hp_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(
            mm2.coeffs.state,
            bw_mm2_coeffs_state_bw_mm2_coeffs_state_set_sample_rate
        );
    }

    #[test]
    fn reset_none() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
        let x0 = 0.0;
        mm2.reset(x0, None);
        assert_eq!(
            mm2.coeffs.state,
            bw_mm2_coeffs_state_bw_mm2_coeffs_state_reset_coeffs
        );
        assert_eq!(mm2.coeffs.reset_id, mm2.states[0].coeffs_reset_id);
    }

    #[test]
    fn reset_some() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;
        let mut y0 = [0.0, 0.0];

        mm2.reset(x0, Some(&mut y0));
        assert_eq!(
            mm2.coeffs.state,
            bw_mm2_coeffs_state_bw_mm2_coeffs_state_reset_coeffs
        );
        assert_eq!(mm2.coeffs.reset_id, mm2.states[0].coeffs_reset_id);
    }

    #[test]
    fn reset_multi() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];
        let mut y0 = [0.0, 0.0];

        mm2.reset_multi(&x0, Some(&mut y0));
        assert_eq!(
            mm2.coeffs.state,
            bw_mm2_coeffs_state_bw_mm2_coeffs_state_reset_coeffs
        );
        assert_eq!(mm2.coeffs.reset_id, mm2.states[1].coeffs_reset_id);
    }

    #[test]
    fn process() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
        let x0 = [0.0, 0.0];
        let mut y0 = [0.0, 0.0];

        mm2.reset_multi(&x0, Some(&mut y0));

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        mm2.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_reset_coeffs);
        assert_eq!(mm2.coeffs.reset_id, mm2.states[0].coeffs_reset_id);
    }

    #[test]
    fn set_cutoff() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);
        let cutoff = 41.20;
        mm2.set_cutoff(cutoff);

        assert_eq!(mm2.coeffs.svf_coeffs.cutoff, cutoff);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_q() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let q = 1.4;
        mm2.set_q(q);

        assert_eq!(mm2.coeffs.svf_coeffs.Q, q);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        mm2.set_prewarp_at_cutoff(true);
        assert_eq!(mm2.coeffs.svf_coeffs.prewarp_k, 1.0);

        mm2.set_prewarp_at_cutoff(false);
        assert_eq!(mm2.coeffs.svf_coeffs.prewarp_k, 0.0);

        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let freq = 185.0;
        mm2.set_prewarp_freq(freq);
        assert_eq!(mm2.coeffs.svf_coeffs.prewarp_freq, freq);

        mm2.set_prewarp_at_cutoff(false);
        assert_eq!(mm2.coeffs.svf_coeffs.prewarp_k, 0.0);

        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_coeff_x() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let gain = 0.9;
        mm2.set_coeff_x(gain);

        assert_eq!(mm2.coeffs.gain_x_coeffs.gain, gain);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_coeff_lp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let gain = 0.9;
        mm2.set_coeff_lp(gain);

        assert_eq!(mm2.coeffs.gain_lp_coeffs.gain, gain);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_coeff_bp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let gain = 0.9;
        mm2.set_coeff_bp(gain);

        assert_eq!(mm2.coeffs.gain_bp_coeffs.gain, gain);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_coeff_hp() {
        let mut mm2 = MM2T::new();
        mm2.set_sample_rate(SAMPLE_RATE);

        let gain = 0.9;
        mm2.set_coeff_hp(gain);

        assert_eq!(mm2.coeffs.gain_hp_coeffs.gain, gain);
        assert!(mm2.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }
}
