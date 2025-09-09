use crate::c_wrapper::{
    bw_clip_set_bias, bw_clip_set_gain, bw_clip_set_gain_compensation, bw_clip_set_sample_rate,
    bw_one_pole_coeffs, bw_one_pole_init, bw_one_pole_set_sticky_thresh, bw_one_pole_set_tau,
    utils::make_array,
};

use super::{bw_clip_coeffs, bw_clip_state};

pub(crate) struct Clip<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_clip_coeffs,
    pub(crate) state: [bw_clip_state; N_CHANNELS],
    pub(crate) states_p: [bw_clip_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Clip<N_CHANNELS> {
    pub fn new() -> Self {
        let mut clip = Clip {
            coeffs: bw_clip_coeffs {
                smooth_coeffs: bw_one_pole_coeffs {
                    ..Default::default()
                },
                smooth_bias_state: Default::default(),
                smooth_gain_state: Default::default(),
                bias_dc: Default::default(),
                inv_gain: Default::default(),
                bias: 0.,
                gain: 1.,
                gain_compensation: 0,
            },
            state: make_array::<bw_clip_state, N_CHANNELS>(),
            states_p: make_array::<bw_clip_state, N_CHANNELS>(),
        };

        unsafe {
            bw_one_pole_init(&mut clip.coeffs.smooth_coeffs);
            bw_one_pole_set_tau(&mut clip.coeffs.smooth_coeffs, 0.005);
            bw_one_pole_set_sticky_thresh(&mut clip.coeffs.smooth_coeffs, 0.001);
        }

        clip
    }

    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_clip_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    pub fn reset(&mut self, x0: &[f32], mut y0: Option<&mut [f32]>) {}
    pub fn process(
        &mut self,
        x: &[Vec<f32>],
        y: Option<&mut [Option<&mut [f32]>]>,
        n_samples: usize,
    ) {
    }

    pub fn set_bias(&mut self, value: f32) {
        unsafe {
            bw_clip_set_bias(&mut self.coeffs, value);
        }
    }

    pub fn set_gain(&mut self, value: f32) {
        unsafe {
            bw_clip_set_gain(&mut self.coeffs, value);
        }
    }

    pub fn set_gain_compensation(&mut self, value: bool) {
        unsafe {
            bw_clip_set_gain_compensation(&mut self.coeffs, value as i8);
        }
    }
}

impl Default for bw_clip_state {
    fn default() -> Self {
        Self {
            x_z1: 0.0,
            F_z1: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::c_wrapper::bw_rcpf;
    use std::f32::consts::PI;

    const N_CHANNELS: usize = 2;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    #[test]
    fn new() {
        let mut clip = Clip::<N_CHANNELS>::new();
        let mut cutoff: f32;

        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(0.005);
        }

        assert_eq!(clip.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(clip.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(clip.coeffs.smooth_coeffs.sticky_thresh, 0.001);
        assert_eq!(clip.coeffs.bias, 0.);
        assert_eq!(clip.coeffs.gain, 1.);
        assert_eq!(clip.coeffs.gain_compensation, 0);
    }

    #[test]
    fn set_sample_rate() {
        const SAMPLE_RATE: f32 = 48_000.0;
        let mut clip = Clip::<N_CHANNELS>::new();

        clip.set_sample_rate(SAMPLE_RATE);
        assert_eq!(clip.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE)
    }

    #[test]
    fn set_bias_default() {
        let clip = Clip::<N_CHANNELS>::new();
        // Default value: 0.f.
        let BIAS = 0.0;

        assert_eq!(clip.coeffs.bias, BIAS);
    }

    #[test]
    fn set_bias_in_range() {
        let mut clip = Clip::<N_CHANNELS>::new();
        let BIAS = 200_000.0;
        clip.set_bias(BIAS);

        assert_eq!(clip.coeffs.bias, BIAS);
    }

    #[test]
    fn set_gain_default() {
        let mut clip = Clip::<N_CHANNELS>::new();
        // Default value: 1.f.
        let GAIN = 1.;

        assert_eq!(clip.coeffs.gain, GAIN);
    }

    #[test]
    fn set_gain_in_range() {
        let mut clip = Clip::<N_CHANNELS>::new();
        let GAIN = 200_000.0;
        clip.set_gain(GAIN);

        assert_eq!(clip.coeffs.gain, GAIN);
    }

    #[test]
    fn set_gain_compensation() {
        let mut clip = Clip::<N_CHANNELS>::new();
        let GAIN_COMPENSATION = true;

        clip.set_gain_compensation(GAIN_COMPENSATION);

        assert_eq!(clip.coeffs.gain_compensation != 0, GAIN_COMPENSATION);
    }
}
