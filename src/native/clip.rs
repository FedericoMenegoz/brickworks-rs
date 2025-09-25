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
//! ```text
//! J. D. Parker, V. Zavalishin, and E. Le Bivic, "Reducing the Aliasing of
//! Nonlinear Waveshaping Using Continuous-Time Convolution", Proc. 19th Intl.
//! Conf. Digital Audio Effects (DAFx-16), pp. 137-144, Brno, Czech Republic,
//! September 2016.
//! ```
//!
//! # Example
//! ```
//! use brickworks_rs::native::clip::Clip;
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
//! ## Notes
//! This module provides a native Rust implementation, but the same interface is
//! also available via bindings to the original C library at [crate::c_wrapper::clip].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_clip).
#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
use crate::native::math::{clipf, rcpf};
use crate::native::one_pole::{OnePoleCoeffs, OnePoleState};
/// Antialiased hard clipper with parametric bias and gain (compensation) and output bias removal.
///
/// Manages both the clipper coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`ClipCoeffs`]
/// - [`ClipState`]
/// # Usage
///
/// ```rust
/// use brickworks_rs::native::clip::Clip;
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
pub struct Clip<const N_CHANNELS: usize> {
    coeffs: ClipCoeffs<N_CHANNELS>,
    states: [ClipState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Clip<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        Clip {
            coeffs: ClipCoeffs::new(),
            states: [ClipState::default(); N_CHANNELS],
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, value: f32) {
        self.coeffs.set_sample_rate(value);
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value x0.
    ///
    /// The corresponding initial output values are written into the y0 array, if it is Some.
    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(y0_values) => {
                (0..N_CHANNELS).for_each(|channel| {
                    y0_values[channel] = self.coeffs.reset_state(&mut self.states[channel], x0);
                });
            }
            None => {
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(&mut self.states[channel], x0);
                });
            }
        }
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value in the x0 array.
    ///
    /// The corresponding initial output values are written into the y0 array, if is Some.
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        self.coeffs.reset_state_multi(&mut self.states, x0, y0);
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and fills the
    /// first `n_samples` of the `N_CHANNELS` output buffers `y`, while using and updating
    /// both the common coeffs and each of the `N_CHANNELS` states (control and audio rate).
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.coeffs.process_multi(&mut self.states, x, y, n_samples);
    }
    /// Sets the input bias `value`.
    ///
    /// Valid range: [-1e12, 1e12].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        self.coeffs.set_bias(value);
    }
    /// Sets the gain `value`.
    ///
    /// Valid range: [1e-12, 1e12].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        self.coeffs.set_gain(value);
    }
    /// Sets whether the output should be divided by gain (true) or not (false).
    ///
    /// Default value: false (off).
    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        self.coeffs.set_gain_compensation(value);
    }
}

impl<const N_CHANNELS: usize> Default for Clip<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
pub struct ClipCoeffs<const N_CHANNELS: usize> {
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

impl<const N_CHANNELS: usize> ClipCoeffs<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = ClipCoeffs {
            smooth_coeffs: OnePoleCoeffs::default(),
            smooth_bias_state: OnePoleState::new(),
            smooth_gain_state: OnePoleState::new(),
            bias_dc: 0.0,
            inv_gain: 0.0,
            bias: 0.0,
            gain: 1.0,
            gain_compensation: false,
        };
        coeffs.smooth_coeffs.set_tau(0.005);
        coeffs.smooth_coeffs.set_sticky_thresh(1e-3);

        coeffs
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_positive(value);
            debug_assert_is_finite(value);
        }
        self.smooth_coeffs.set_sample_rate(value);
        self.smooth_coeffs.reset_coeffs();
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.smooth_coeffs
            .reset_state(&mut self.smooth_bias_state, self.bias);
        self.smooth_coeffs
            .reset_state(&mut self.smooth_gain_state, self.gain);
        self.do_update_coeffs(true);
    }
    /// Resets the given state to its initial values using the initial
    /// input value `x0`.
    ///
    /// Returns the corresponding initial output value.
    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut ClipState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        let x = self.smooth_gain_state.get_y_z1() * x0 + self.smooth_bias_state.get_y_z1();
        let a = x.abs();
        let f = if a > 1.0 { a - 0.5 } else { 0.5 * a * a };
        let yb = clipf(x, -1.0, 1.0);
        let y = (if self.gain_compensation {
            self.inv_gain
        } else {
            1.0
        }) * (yb - self.bias_dc);
        state.x_z1 = x;
        state.f_z1 = f;

        debug_assert!(y.is_finite());

        y
    }
    /// Resets each of the `N_CHANNELS` states to its initial values using
    /// the corresponding initial input value in the `x0` array.
    ///
    /// The corresponding initial output values are written into the `y0` array,
    /// if not None.
    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [ClipState; N_CHANNELS],
        x0: &[f32],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        if let Some(y0_values) = y0 {
            (0..N_CHANNELS).for_each(|channel| {
                y0_values[channel] = self.reset_state(&mut states[channel], x0[channel]);
            });
        } else {
            (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut states[channel], x0[channel]);
            });
        }
    }

    // Not implemented yet: C version only contained assertions
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Triggers control-rate update of coefficients.
    // #[inline(always)]
    // pub fn update_coeffs_ctrl(&mut self) {
    //     todo!()
    // }

    /// Triggers audio-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.do_update_coeffs(false);
    }
    /// Processes one input sample `x`, while using and updating state.
    ///
    /// Assumes that gain compensation is disabled.
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1(&mut self, state: &mut ClipState, mut x: f32) -> f32 {
        debug_assert!(x.is_finite());

        x = self.smooth_gain_state.get_y_z1() * x + self.smooth_bias_state.get_y_z1();
        let a = x.abs();
        let f = if a > 1.0 { a - 0.5 } else { 0.5 * a * a };
        let d = x - state.x_z1;
        let yb = if d * d < 1e-6 {
            clipf(0.5 * (x + state.x_z1), -1.0, 1.0)
        } else {
            (f - state.f_z1) * rcpf(d)
        };
        let y = yb - self.bias_dc;
        state.x_z1 = x;
        state.f_z1 = f;

        debug_assert!(y.is_finite());
        y
    }
    /// Processes one input sample `x`, while using and updating state.
    ///
    /// Assumes that gain compensation is enabled.
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_comp(&mut self, state: &mut ClipState, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let y = self.inv_gain * self.process1(state, x);

        debug_assert!(y.is_finite());
        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and fills the first
    /// `n_samples` of the output buffer `y`, while using and updating both coeffs
    /// and state (control and audio rate).
    #[inline(always)]
    pub fn process(&mut self, state: &mut ClipState, x: &[f32], y: &mut [f32], n_samples: usize) {
        if self.gain_compensation {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                y[sample] = self.process1_comp(state, x[sample]);
            });
        } else {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                y[sample] = self.process1(state, x[sample]);
            });
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and fills the
    /// first `n_samples` of the `N_CHANNELS` output buffers `y`, while using and updating
    /// both the common coeffs and each of the `N_CHANNELS` states (control and audio rate).
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [ClipState; N_CHANNELS],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        if self.gain_compensation {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel][sample] =
                        self.process1_comp(&mut states[channel], x[channel][sample]);
                });
            });
        } else {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel][sample] = self.process1(&mut states[channel], x[channel][sample]);
                });
            });
        }
    }
    /// Sets the input bias `value`.
    ///
    /// Valid range: [-1e12, 1e12].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite());
            debug_assert_range(-1e12..=1e12, value);
        }

        self.bias = value;
    }
    /// Sets the gain `value`.
    ///
    /// Valid range: [1e-12, 1e12].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite());
            debug_assert_range(1e-12..=1e12, value);
        }

        self.gain = value;
    }
    /// Sets whether the output should be divided by gain (true) or not (false).
    ///
    /// Default value: false (off).
    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        self.gain_compensation = value;
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    /// Tries to determine whether coeffs is valid and returns true if it seems to
    /// be the case and false if it is certainly not. False positives are possible,
    /// false negatives are not.
    /// # Notes
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn coeffs_is_valid(&mut self) -> bool {
        todo!()
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    /// Tries to determine whether state is valid and returns true if it seems to
    /// be the case and false if it is certainly not. False positives are possible,
    /// false negatives are not.
    /// # Notes
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }

    // Private
    fn do_update_coeffs(&mut self, force: bool) {
        let mut bias_cur = self.smooth_bias_state.get_y_z1();
        if force || self.bias != bias_cur {
            bias_cur = self
                .smooth_coeffs
                .process1_sticky_abs(&mut self.smooth_bias_state, self.bias);
            self.bias_dc = clipf(bias_cur, -1.0, 1.0);
        }
        let mut gain_cur = self.smooth_gain_state.get_y_z1();
        if force || self.gain != gain_cur {
            gain_cur = self
                .smooth_coeffs
                .process1_sticky_rel(&mut self.smooth_gain_state, self.gain);
            self.inv_gain = rcpf(gain_cur);
        }
    }
}

impl<const N_CHANNELS: usize> Default for ClipCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Internal state and related.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ClipState {
    x_z1: f32,
    f_z1: f32,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::Clip;
    use crate::{
        c_wrapper::{bw_clip_coeffs, bw_clip_state, clip::Clip as ClipWrapper},
        native::{
            clip::{ClipCoeffs, ClipState},
            one_pole::tests::assert_one_pole_coeffs,
        },
    };
    use std::f32;
    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    type Clip2 = Clip<N_CHANNELS>;
    type ClipWrapper2 = ClipWrapper<N_CHANNELS>;

    #[test]
    fn new_clip() {
        let rust_clip = Clip2::new();
        let c_clip = ClipWrapper2::new();

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut rust_clip = Clip2::new();
        let mut c_clip = ClipWrapper2::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_sample_rate_must_be_finite() {
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(f32::INFINITY);
    }

    #[test]
    #[should_panic(expected = "value must be in range [-1e12, 1e12], got -1e13")]
    fn set_bias_invalid() {
        let mut rust_clip = Clip2::new();
        rust_clip.set_bias(-1e13);
    }

    #[test]
    #[should_panic(expected = "value must be in range [1e-12, 1e12], got 0e0")]
    fn set_gain_invalid() {
        let mut rust_clip = Clip2::new();
        rust_clip.set_gain(0.0);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_sample_rate_must_be_positive() {
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(-1.);
    }

    #[test]
    fn reset_none() {
        let bias = 5.0;
        let gain = 10.0;
        let x0 = [10.0, 11.0];
        let mut c_clip = ClipWrapper2::new();
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset_multi(&x0, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset_multi(&x0, None);

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn reset_some() {
        let bias = 0.5;
        let gain = 0.01;
        let x0 = [0.1, 0.02];
        let mut rust_y0 = [0.2, 0.7];
        let mut c_y0 = [0.9, 0.7];
        let mut c_clip = ClipWrapper2::new();
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset_multi(&x0, Some(&mut c_y0));

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset_multi(&x0, Some(&mut rust_y0));

        assert_clip(&rust_clip, &c_clip);
        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y0[channel], c_y0[channel]));
    }

    #[test]
    fn process1() {
        let bias = 0.5;
        let gain = 0.01;

        let x = [0.1, 0.02];
        let mut rust_y = [0.2, 0.7];
        let mut c_y = [0.9, 0.7];

        let mut c_clip = ClipWrapper2::new();
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset_multi(&x, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset_multi(&x, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_clip
                .coeffs
                .process1(&mut c_clip.states[channel], x[channel]);
            rust_y[channel] = rust_clip
                .coeffs
                .process1(&mut rust_clip.states[channel], x[channel]);
        });

        assert_clip(&rust_clip, &c_clip);
        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y[channel], c_y[channel]));
    }

    #[test]
    fn process1_comp() {
        let bias = 0.5;
        let gain = 0.01;

        let x = [0.1, 0.02];
        let mut rust_y = [0.2, 0.7];
        let mut c_y = [0.9, 0.7];

        let mut c_clip = ClipWrapper2::new();
        let mut rust_clip = Clip2::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(true);
        c_clip.reset_multi(&x, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(true);
        rust_clip.reset_multi(&x, None);

        (0..N_CHANNELS).for_each(|channel| {
            rust_y[channel] = rust_clip
                .coeffs
                .process1_comp(&mut rust_clip.states[channel], x[channel]);
            c_y[channel] = c_clip
                .coeffs
                .process1_comp(&mut c_clip.states[channel], x[channel]);
            assert_clip(&rust_clip, &c_clip);
        });

        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y[channel], c_y[channel]));
    }

    #[test]
    fn process_with_comp() {
        let n_samples = 2;
        let bias = 0.5;
        let gain = 0.01;
        let x_ch0 = [0.1, 0.02];
        let x_ch1 = [0.1, 0.02];
        let x: [&[f32]; N_CHANNELS] = [&x_ch0, &x_ch1];

        let mut rust_y_ch0 = [0.2, 0.7];
        let mut rust_y_ch1 = [0.2, 0.7];
        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut rust_y_ch0, &mut rust_y_ch1];

        let mut c_y_ch0 = [0.2, 0.7];
        let mut c_y_ch1 = [0.2, 0.7];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut c_y_ch0, &mut c_y_ch1];

        let mut rust_clip = Clip2::new();
        rust_clip.set_sample_rate(SAMPLE_RATE);
        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(true);
        rust_clip.reset_multi(&[0.0; N_CHANNELS], None);

        let mut c_clip = ClipWrapper2::new();
        c_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(true);
        c_clip.reset_multi(&[0.0; N_CHANNELS], None);

        rust_clip.process(&x, &mut rust_y, n_samples);
        c_clip.process(&x, &mut c_y, n_samples);

        (0..N_CHANNELS).for_each(|channel| {
            (0..n_samples).for_each(|sample| {
                assert_eq!(rust_y[channel][sample], c_y[channel][sample]);
            });
        });
        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn process_without_comp() {
        let n_samples = 2;
        let bias = 0.5;
        let gain = 0.01;
        let x_ch0 = [0.1, 0.02];
        let x_ch1 = [0.1, 0.02];
        let x: [&[f32]; N_CHANNELS] = [&x_ch0, &x_ch1];

        let mut rust_y_ch0 = [0.2, 0.7];
        let mut rust_y_ch1 = [0.2, 0.7];
        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut rust_y_ch0, &mut rust_y_ch1];

        let mut c_y_ch0 = [0.2, 0.7];
        let mut c_y_ch1 = [0.2, 0.7];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut c_y_ch0, &mut c_y_ch1];

        let mut rust_clip = Clip2::new();
        rust_clip.set_sample_rate(SAMPLE_RATE);
        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(false); // it is default, but to be clear
        rust_clip.reset_multi(&[0.0; N_CHANNELS], None);

        let mut c_clip = ClipWrapper2::new();
        c_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(false); // it is default, but to be clear
        c_clip.reset_multi(&[0.0; N_CHANNELS], None);

        rust_clip.process(&x, &mut rust_y, n_samples);
        c_clip.process(&x, &mut c_y, n_samples);

        (0..N_CHANNELS).for_each(|channel| {
            (0..n_samples).for_each(|sample| {
                assert_eq!(rust_y[channel][sample], c_y[channel][sample]);
            });
        });
        assert_clip(&rust_clip, &c_clip);
    }

    fn assert_clip<const N_CHANNELS: usize>(
        rust_clip: &Clip<N_CHANNELS>,
        c_clip: &ClipWrapper<N_CHANNELS>,
    ) {
        assert_clip_coeffs(&rust_clip.coeffs, &c_clip.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_clip_state(&rust_clip.states[channel], &c_clip.states[channel]);
        });
    }

    pub(crate) fn assert_clip_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &ClipCoeffs<N_CHANNELS>,
        c_coeffs: &bw_clip_coeffs,
    ) {
        let pre_message = "clip.coeff.";
        let post_message = "does not match";
        assert_eq!(
            rust_coeffs.bias_dc, c_coeffs.bias_dc,
            "{}bias_dc {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.inv_gain, c_coeffs.inv_gain,
            "{}inv_gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.bias, c_coeffs.bias,
            "{}bias {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain, c_coeffs.gain,
            "{}gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain_compensation,
            c_coeffs.gain_compensation != 0,
            "{}gain_compensation {}",
            pre_message,
            post_message
        );
        assert_one_pole_coeffs(&rust_coeffs.smooth_coeffs, &c_coeffs.smooth_coeffs);
        assert_eq!(
            rust_coeffs.smooth_bias_state.get_y_z1(),
            c_coeffs.smooth_bias_state.y_z1
        );
        assert_eq!(
            rust_coeffs.smooth_gain_state.get_y_z1(),
            c_coeffs.smooth_gain_state.y_z1
        );
    }

    pub(crate) fn assert_clip_state(rust_state: &ClipState, c_state: &bw_clip_state) {
        assert_eq!(rust_state.x_z1, c_state.x_z1);
        assert_eq!(rust_state.f_z1, c_state.F_z1);
    }
}
