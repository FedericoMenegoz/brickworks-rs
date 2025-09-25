//! **Antialiased hard clipper** with parametric bias and gain (compensation) and output bias removal.
//!
//! In other words this implements (approximately)
//! ```text
//! y(n) = clip(gain * x(n) + bias, -1, 1) - clip(bias, -1, 1)
//! ```
//!
//! with antialiasing and optionally dividing the output by gain.
//!
//! As a side effect, antialiasing causes attenuation at higher frequencies
//! (about 3 dB at 0.5 × Nyquist frequency and rapidly increasing at higher
//! frequencies).
//!
//! The antialiasing technique used here is described in
//! > J. D. Parker, V. Zavalishin, and E. Le Bivic, "Reducing the Aliasing of Nonlinear Waveshaping Using Continuous-Time Convolution", Proc. 19th Intl. Conf. Digital Audio Effects (DAFx-16), pp. 137-144, Brno, Czech Republic, September 2016.
//! # Example
//! ```
//! use brickworks_rs::c_wrapper::clip::Clip;
//!
//! // Create a stereo (2-channel) distortion processor
//! const N_CHANNELS: usize = 2;
//! let mut clip = Clip::<N_CHANNELS>::new();
//!
//! // Configure processor
//! clip.set_sample_rate(48_000.0);
//! clip.set_bias(11.0);
//! clip.set_gain(0.5);
//! clip.set_gain_compensation(true);
//!
//! // Prepare input and output buffers
//! let input_left: Vec<f32> = vec![0.0; 512];
//! let input_right: Vec<f32> = vec![0.0; 512];
//! let mut output_left: Vec<f32> = vec![0.0; 512];
//! let mut output_right: Vec<f32> = vec![0.0; 512];
//!
//! // Build channel slices
//! let inputs = [&input_left[..], &input_right[..]];
//! let mut outputs = [&mut output_left[..], &mut output_right[..]];
//!
//! // Reset before processing
//! clip.reset(0.0, None);
//! // Process 512 samples
//! clip.process(&inputs, &mut outputs, 512);
//! ```
//!
//! # Notes
//! This module provides Rust bindings to the original C implementation.
//! For a fully native Rust implementation with the same interface,
//! see [crate::native::clip].
//! Original C library by [Orastron](https://www.orastron.com/algorithms/bw_clip).
use std::ptr::null_mut;

use super::*;
use crate::c_wrapper::utils::make_array;
/// Antialiased hard clipper with parametric bias and gain (compensation) and output bias removal.
///
/// Manages both the clipper coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
///
/// # Usage
///
/// ```rust
/// use brickworks_rs::c_wrapper::clip::Clip;
///
/// // Stereo filter
/// const N_CHANNELS: usize = 2;
/// let mut clip = Clip::<N_CHANNELS>::new();
///
/// clip.set_sample_rate(44_100.0);
/// clip.set_bias(0.3);
/// clip.set_gain(0.8);
/// clip.reset(0.0, None);
/// // dist.process(...)
/// ```
#[derive(Debug)]
pub struct Clip<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_clip_coeffs,
    pub(crate) states: [bw_clip_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Clip<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
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
    /// Sets the sample rate (Hz).
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_clip_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value x0.
    ///
    /// The corresponding initial output values are written into the y0 array, if it is Some.
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_clip_reset_coeffs(&mut self.coeffs);
            if let Some(out) = y0 {
                (0..N_CHANNELS).for_each(|channel| {
                    out[channel] =
                        bw_clip_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                });
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    bw_clip_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                });
            }
        }
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value in the x0 array.
    ///
    /// The corresponding initial output values are written into the y0 array, if is Some.
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_clip_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_clip_reset_coeffs(&mut self.coeffs);
            bw_clip_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and fills the
    /// first `n_samples` of the `N_CHANNELS` output buffers `y`, while using and updating
    /// both the common coeffs and each of the `N_CHANNELS` states (control and audio rate).
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
    /// Sets the input bias `value`.
    ///
    /// Valid range: [-1e12, 1e12].
    ///
    /// Default value: 0.0.
    pub fn set_bias(&mut self, value: f32) {
        unsafe {
            bw_clip_set_bias(&mut self.coeffs, value);
        }
    }
    /// Sets the gain `value`.
    ///
    /// Valid range: [1e-12, 1e12].
    ///
    /// Default value: 1.0.
    pub fn set_gain(&mut self, value: f32) {
        unsafe {
            bw_clip_set_gain(&mut self.coeffs, value);
        }
    }
    /// Sets whether the output should be divided by gain (true) or not (false).
    ///
    /// Default value: false (off).
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
        assert_eq!(
            clip.coeffs.state,
            bw_clip_coeffs_state_bw_clip_coeffs_state_init
        );
    }

    #[test]
    fn set_sample_rate() {
        let mut clip = ClipT::new();

        clip.set_sample_rate(SAMPLE_RATE);
        assert_eq!(clip.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
        assert_eq!(
            clip.coeffs.state,
            bw_clip_coeffs_state_bw_clip_coeffs_state_set_sample_rate
        );
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
        clip.reset_multi(&x0, Some(&mut out));

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
        clip.reset_multi(&x0, None);
        clip.process(&PULSE_INPUT, &mut y, N_SAMPLES);

        assert!(clip.coeffs.state >= bw_clip_coeffs_state_bw_clip_coeffs_state_reset_coeffs);
    }
}
