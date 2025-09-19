use super::*;
use crate::c_wrapper::utils::make_array;

#[derive(Debug)]
pub struct Clip<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_clip_coeffs,
    pub(crate) states: [bw_clip_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Clip<N_CHANNELS> {
    pub fn new() -> Self {
        let mut clip = Clip {
            coeffs: bw_clip_coeffs::default(),
            states: make_array::<bw_clip_state, N_CHANNELS>(),
        };

        unsafe {
            bw_clip_init(&mut clip.coeffs);
        }
        println!("{:?}", clip);
        clip
    }

    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_clip_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    pub fn reset(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_clip_reset_coeffs(&mut self.coeffs);
            if let Some(out) = y0 {
                (0..N_CHANNELS).for_each(|channel| {
                    out[channel] = bw_clip_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0[channel],
                    );
                });
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    bw_clip_reset_state(&mut self.coeffs, &mut self.states[channel], x0[channel]);
                });
            }
        }
    }
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let mut state_ptrs: [*mut bw_clip_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut bw_clip_state);
        unsafe {
            bw_clip_process_multi(
                &mut self.coeffs,
                state_ptrs.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
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

impl bw_clip_coeffs {
    #[cfg(test)]
    pub(crate) fn process1(&mut self, state: &mut bw_clip_state, x: f32) -> f32 {
        unsafe { bw_clip_process1(self, state, x) }
    }

    #[cfg(test)]
    pub(crate) fn process1_comp(&mut self, state: &mut bw_clip_state, x: f32) -> f32 {
        unsafe { bw_clip_process1_comp(self, state, x) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::c_wrapper::{
        bw_clipf, bw_one_pole_get_y_z1, bw_one_pole_process1_sticky_abs,
        bw_one_pole_process1_sticky_rel, bw_rcpf,
    };
    use std::f32::consts::PI;

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;
    const GAIN: f32 = 2.0;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [&[1.0, 1.0], &[0.0, 0.0]];
    const N_SAMPLES: usize = 2;

    type ClipT = Clip<N_CHANNELS>;

    #[test]
    fn new() {
        let clip = ClipT::new();
        let cutoff: f32;
        let tau_default = 0.005;

        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(clip.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(clip.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(clip.coeffs.smooth_coeffs.sticky_thresh, 0.001);
        assert_eq!(clip.coeffs.bias, 0.);
        assert_eq!(clip.coeffs.gain, 1.);
        assert_eq!(clip.coeffs.gain_compensation, 0);
        assert_eq!(clip.coeffs.state, bw_clip_coeffs_state_bw_clip_coeffs_state_init);
    }
    
    #[test]
    fn set_sample_rate() {
        let mut clip = ClipT::new();
        
        clip.set_sample_rate(SAMPLE_RATE);
        assert_eq!(clip.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
        assert_eq!(clip.coeffs.state, bw_clip_coeffs_state_bw_clip_coeffs_state_set_sample_rate);
    }

    #[test]
    fn set_bias_default() {
        let clip = ClipT::new();
        // Default value: 0.f.
        let BIAS = 0.0;

        assert_eq!(clip.coeffs.bias, BIAS);
    }

    #[test]
    fn set_bias_in_range() {
        let mut clip = ClipT::new();
        let BIAS = 200_000.0;
        clip.set_bias(BIAS);

        assert_eq!(clip.coeffs.bias, BIAS);
    }

    #[test]
    fn set_gain_default() {
        let clip = ClipT::new();
        // Default value: 1.f.
        let gain = 1.;
        assert_eq!(clip.coeffs.gain, gain);
    }

    #[test]
    fn set_gain_in_range() {
        let mut clip = ClipT::new();
        clip.set_gain(GAIN);

        assert_eq!(clip.coeffs.gain, GAIN);
    }

    #[test]
    fn set_gain_compensation() {
        let mut clip = ClipT::new();
        let GAIN_COMPENSATION = true;

        clip.set_gain_compensation(GAIN_COMPENSATION);

        assert_eq!(clip.coeffs.gain_compensation != 0, GAIN_COMPENSATION);
    }

    #[test]
    fn reset() {
        let mut clip = ClipT::new();
        let x0: [f32; N_CHANNELS] = [6.0, 2.0];
        let mut out: [f32; N_CHANNELS] = [3.0, 4.0];

        clip.set_sample_rate(SAMPLE_RATE);
        clip.set_gain_compensation(true);
        clip.reset(&x0, Some(&mut out));

        let inv_gain;
        let bias_dc;
        let y;
        unsafe {
            inv_gain = bw_rcpf(bw_one_pole_process1_sticky_abs(
                &clip.coeffs.smooth_coeffs,
                &mut clip.coeffs.smooth_gain_state,
                clip.coeffs.gain,
            ));
            bias_dc = bw_one_pole_process1_sticky_rel(
                &clip.coeffs.smooth_coeffs,
                &mut clip.coeffs.smooth_bias_state,
                clip.coeffs.bias,
            );
            let x = bw_one_pole_get_y_z1(&clip.coeffs.smooth_gain_state) * x0[0]
                + bw_one_pole_get_y_z1(&clip.coeffs.smooth_bias_state);
            let yb = bw_clipf(x, -1., 1.);
            y = if clip.coeffs.gain_compensation != 0 {
                clip.coeffs.inv_gain
            } else {
                1.
            } * (yb - clip.coeffs.bias_dc)
        }

        assert_eq!(out[0], y);
        assert_eq!(clip.coeffs.inv_gain, inv_gain);
        assert_eq!(clip.coeffs.bias_dc, bias_dc);
    }

    #[test]
    fn process() {
        let mut clip = ClipT::new();
        clip.set_sample_rate(SAMPLE_RATE);
        let x0 = [0.0, 0.0];

        let mut out_0: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut out_1: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut y: [&mut [f32]; 2] = [&mut out_0, &mut out_1];

        clip.set_gain_compensation(true);
        clip.reset(&x0, None);
        clip.process(&PULSE_INPUT, &mut y, N_SAMPLES);

        assert!(clip.coeffs.state >= bw_clip_coeffs_state_bw_clip_coeffs_state_reset_coeffs);
    }
}
