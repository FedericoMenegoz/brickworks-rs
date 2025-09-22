use std::ptr::null_mut;

use super::*;

pub struct LP1<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_lp1_coeffs,
    pub(crate) states: [bw_lp1_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> LP1<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_lp1_coeffs::default();
        unsafe {
            bw_lp1_init(&mut coeffs);
        }

        Self {
            coeffs,
            states: [bw_lp1_state::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_lp1_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_lp1_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] =
                        bw_lp1_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_lp1_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
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
        let states_ptrs: [*mut bw_lp1_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_lp1_reset_coeffs(&mut self.coeffs);
            bw_lp1_reset_state_multi(
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
        let states_ptrs: [*mut bw_lp1_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_lp1_process_multi(
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
            bw_lp1_set_cutoff(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        unsafe {
            bw_lp1_set_prewarp_at_cutoff(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        unsafe {
            bw_lp1_set_prewarp_freq(&mut self.coeffs, value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::math::INVERSE_2_PI;
    use std::f32::consts::PI;

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    const N_SAMPLES: usize = 8;

    type LP1T = LP1<N_CHANNELS>;

    #[test]
    pub fn new() {
        let lp1 = LP1T::new();

        let tau_default = 0.005;
        let sticky_tresh_default = 1e-3;
        let cutoff;
        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(lp1.coeffs.cutoff, 1e3);
        assert_eq!(lp1.coeffs.prewarp_k, 1.0);
        assert_eq!(lp1.coeffs.prewarp_freq, 1e3);
        assert_eq!(lp1.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(lp1.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(lp1.coeffs.smooth_coeffs.sticky_thresh, sticky_tresh_default);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_init
        );
    }

    #[test]
    pub fn set_sample_rate() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        assert_eq!(lp1.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
        assert_eq!(lp1.coeffs.t_k, PI / SAMPLE_RATE);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_set_sample_rate
        );
    }

    #[test]
    pub fn reset_none() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;
        lp1.reset(x0, None);

        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        )
    }

    #[test]
    pub fn reset_some() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = 4.0;
        let mut y0 = [0.0, 0.0];
        lp1.reset(x0, Some(&mut y0));

        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(y0[channel], x0);
            assert_eq!(lp1.states[channel].X_z1, 0.0);
            assert_eq!(lp1.states[channel].y_z1, x0);
        });
    }

    #[test]
    pub fn reset_multi() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.1, 0.2];
        let mut y0 = [0.0, 0.0];
        lp1.reset_multi(&x0, Some(&mut y0));

        assert_eq!(y0, x0);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(lp1.states[channel].X_z1, 0.0);
            assert_eq!(lp1.states[channel].y_z1, x0[channel]);
        });
    }

    #[test]
    pub fn process() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];

        lp1.reset_multi(&x0, None);

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        lp1.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(lp1.coeffs.state >= bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs);
        assert_eq!(lp1.coeffs.reset_id, lp1.states[0].coeffs_reset_id);
    }

    #[test]
    pub fn set_cutoff() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 493.883;

        lp1.set_cutoff(cutoff);

        assert_eq!(lp1.coeffs.cutoff, cutoff);
    }

    #[test]
    pub fn set_prewarp_at_cutoff() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        lp1.set_prewarp_at_cutoff(true);
        assert_eq!(lp1.coeffs.prewarp_k, 1.0);

        lp1.set_prewarp_at_cutoff(false);
        assert_eq!(lp1.coeffs.prewarp_k, 0.0);
    }

    #[test]
    pub fn set_prewarp_freq() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 185.0;

        lp1.set_prewarp_freq(prewarp_freq);

        assert_eq!(lp1.coeffs.prewarp_freq, prewarp_freq);
    }
}
