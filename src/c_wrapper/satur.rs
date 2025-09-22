use super::*;
use std::ptr::null_mut;

pub struct Satur<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_satur_coeffs,
    pub(crate) states: [bw_satur_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Satur<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_satur_coeffs::default();
        unsafe {
            bw_satur_init(&mut coeffs);
        }
        Self {
            coeffs,
            states: [bw_satur_state::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_satur_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32]>) {
        unsafe {
            bw_satur_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] =
                        bw_satur_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_satur_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
            }
        }
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_satur_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_satur_reset_coeffs(&mut self.coeffs);
            bw_satur_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
    }

    #[inline(always)]
    pub fn process(&mut self, x: &[&[f32]; N_CHANNELS], y: &mut [&mut [f32]], n_samples: usize) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_satur_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_satur_process_multi(
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
    pub fn set_bias(&mut self, value: f32) {
        unsafe {
            bw_satur_set_bias(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        unsafe {
            bw_satur_set_gain(&mut self.coeffs, value);
        }
    }

    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        unsafe {
            bw_satur_set_gain_compensation(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }
}

impl bw_satur_coeffs {
    #[inline(always)]
    pub(crate) fn tanhf(x: f32) -> f32 {
        unsafe { bw_satur_tanhf(x) }
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

    type SaturT = Satur<N_CHANNELS>;

    #[test]
    fn new() {
        let satur = SaturT::new();

        assert_eq!(satur.coeffs.bias, 0.0);
        assert_eq!(satur.coeffs.gain, 1.0);
        assert_eq!(satur.coeffs.gain_compensation, 0);
        assert_eq!(
            satur.coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_init
        );
    }

    #[test]
    fn set_sample_rate() {
        let mut satur = SaturT::new();
        satur.set_sample_rate(SAMPLE_RATE);

        assert_eq!(
            satur.coeffs.smooth_coeffs.fs_2pi,
            INVERSE_2_PI * SAMPLE_RATE
        );
        assert_eq!(
            satur.coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_set_sample_rate
        );
    }

    #[test]
    fn reset_none() {
        let mut satur = SaturT::new();
        satur.set_sample_rate(SAMPLE_RATE);
        let x0 = 0.0;
        satur.reset(x0, None);
        assert_eq!(
            satur.coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_reset_coeffs
        );
        assert_eq!(satur.coeffs.reset_id, satur.states[0].coeffs_reset_id);
    }

    #[test]
    fn reset_some() {
        let mut satur = SaturT::new();
        satur.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;
        let mut y0 = [0.0, 0.0];

        satur.reset(x0, Some(&mut y0));
        assert_eq!(
            satur.coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_reset_coeffs
        );
        assert_eq!(satur.coeffs.reset_id, satur.states[0].coeffs_reset_id);
    }

    #[test]
    fn reset_multi() {
        let mut satur = SaturT::new();
        satur.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];
        let mut y0 = [0.0, 0.0];

        satur.reset_multi(&x0, Some(&mut y0));
        assert_eq!(
            satur.coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_reset_coeffs
        );
        assert_eq!(satur.coeffs.reset_id, satur.states[1].coeffs_reset_id);
    }

    #[test]
    fn process() {
        let mut satur = SaturT::new();
        satur.set_sample_rate(SAMPLE_RATE);
        let x0 = [0.0, 0.0];

        satur.reset_multi(&x0, None);

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        satur.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(satur.coeffs.state >= bw_satur_coeffs_state_bw_satur_coeffs_state_reset_coeffs);
        assert_eq!(satur.coeffs.reset_id, satur.states[0].coeffs_reset_id);
    }

    #[test]
    fn set_bias() {
        let mut satur = SaturT::new();
        let bias = 0.2;
        satur.set_bias(bias);

        assert_eq!(satur.coeffs.bias, bias);
    }

    #[test]
    fn set_gain() {
        let mut satur = SaturT::new();
        let gain = 5.0;
        satur.set_gain(gain);

        assert_eq!(satur.coeffs.gain, gain);
    }

    #[test]
    fn set_gain_compensation() {
        let mut satur = SaturT::new();

        satur.set_gain_compensation(true);

        assert!(satur.coeffs.gain_compensation != 0);
    }
}
