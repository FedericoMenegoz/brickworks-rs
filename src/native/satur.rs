//! Antialiased tanh-based saturation with parametric bias and gain (compensation) and
//! output bias removal.
//!
//! In other words this implements (approximately):
//!
//! ```text
//! y(n) = tanh(gain * x(n) + bias) - tanh(bias)
//! ```
//!
//! with antialiasing and optionally dividing the output by gain.
//!
//! As a side effect, antialiasing causes attenuation at higher frequencies (about
//! 3 dB at 0.5 × Nyquist frequency and rapidly increasing at higher frequencies).
//!
//! # Examples
//! ```rust
//! use brickworks_rs::native::satur::*;
//!
//! const N_CHANNELS: usize = 2;
//! const SAMPLE_RATE: f32 = 44_100.0;
//! const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
//!     &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!     &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//! ];
//! const N_SAMPLES: usize = 8;
//!
//! fn main() {
//!     let mut satur = Satur::new();
//!     satur.set_sample_rate(SAMPLE_RATE);
//!
//!     let x0 = [0.0, 0.0];
//!
//!     satur.reset_multi(&x0, None);
//!
//!     let mut y: [&mut [f32]; 2] = [
//!         &mut [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
//!         &mut [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
//!     ];
//!     satur.process(&PULSE_INPUT, &mut y, N_SAMPLES);
//! }
//! ```
//!
//! # Notes
//!
//! The antialiasing technique used here is described in:
//!
//! J. D. Parker, V. Zavalishin, and E. Le Bivic, "Reducing the Aliasing of Nonlinear
//! Waveshaping Using Continuous-Time Convolution", Proc. 19th Intl. Conf. Digital
//! Audio Effects (DAFx-16), pp. 137-144, Brno, Czech Republic, September 2016.
//!
//! This module provides a native Rust implementation, but the same interface is
//! also available via bindings to the original C library at [`crate::c_wrapper::satur`].
//!
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_satur).
use crate::native::{
    math::{clipf, rcpf},
    one_pole::{OnePoleCoeffs, OnePoleState},
};

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
/// Antialiased tanh-based saturation with parametric bias and gain (compensation) and
/// output bias removal.
///
/// 
/// This struct manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`SaturCoeffs`] 
/// - [`SaturState`] 
///
/// # Usage
/// ```rust 
/// use brickworks_rs::native::satur::Satur;
/// const N_CHANNELS: usize = 2;
/// let mut satur = Satur::<N_CHANNELS>::new();
/// satur.set_sample_rate(48_000.0);
/// satur.set_bias(1.0);
/// satur.set_gain(0.8);
/// satur.set_gain_compensation(true);
/// // process audio with satur.process(...)
/// ```
pub struct Satur<const N_CHANNELS: usize> {
    pub(crate) coeffs: SaturCoeffs<N_CHANNELS>,
    pub(crate) states: [SaturState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Satur<N_CHANNELS> {
    /// Creates a new instance with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        Satur {
            coeffs: SaturCoeffs::new(),
            states: [SaturState::default(); N_CHANNELS],
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value x0.
    ///
    /// The corresponding initial output values are written into the y0 array, if it is Some.
    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(y) => {
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel] = self.coeffs.reset_state(&mut self.states[channel], x0);
                });
            }
            None => {
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(&mut self.states[channel], x0);
                });
            }
        }
    }
    /// Resets the satur's coeffs and each of the `N_CHANNELS` states to its initial values
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
    /// both the common coeffs and each of the N_CHANNELS states (control and audio rate).
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.coeffs.process_multi(&mut self.states, x, y, n_samples);
    }
    /// Sets the input bias value.
    ///
    /// Valid range: [-1e12, 1e12].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        self.coeffs.set_bias(value);
    }
    /// Sets the gain value.
    ///
    /// Valid range: [1e-12, 1e12].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        self.coeffs.set_gain(value);
    }
    /// Sets whether the output should be divided by gain (`true`) or not (`false`).
    ///
    /// Default value: `false` (off).
    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        self.coeffs.set_gain_compensation(value);
    }
}

impl<const N_CHANNELS: usize> Default for Satur<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
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
    /// Creates a new instance with all field initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut smooth_coeffs = OnePoleCoeffs::<N_CHANNELS>::new();
        let smooth_bias_state = OnePoleState::new();
        let smooth_gain_state = OnePoleState::new();

        smooth_coeffs.set_tau(0.005);
        smooth_coeffs.set_sticky_thresh(1e-3);

        Self {
            smooth_coeffs,
            smooth_bias_state,
            smooth_gain_state,
            bias_dc: 0.0,
            inv_gain: 0.0,
            bias: 0.0,
            gain: 1.0,
            gain_compensation: false,
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_positive(sample_rate);
            debug_assert_is_finite(sample_rate);
        }

        self.smooth_coeffs.set_sample_rate(sample_rate);
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
    /// Resets the given state to its initial values using the initial input value x0.
    ///
    /// Returns the corresponding initial output value.
    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut SaturState, x0: f32) -> f32 {
        #[cfg(debug_assertions)]
        debug_assert_is_finite(x0);

        let x = self.smooth_gain_state.get_y_z1() * x0 + self.smooth_bias_state.get_y_z1();
        let ax = x.abs();
        // let f = if ax >= 2.115287308554551  { ax - 0.6847736211329452 } else {ax * ax * ((0.00304518315009429 * ax - 0.09167437770414569) * ax + 0.5)};
        let f = if ax >= 2.1152873 {
            ax - 0.684_773_6
        } else {
            ax * ax * ((0.003_045_183_1 * ax - 0.091_674_38) * ax + 0.5)
        };
        let yb = Self::tanhf(x);
        let y = (if self.gain_compensation {
            self.inv_gain
        } else {
            1.0
        }) * (yb - self.bias_dc);
        state.x_z1 = x;
        state.f_z1 = f;

        #[cfg(debug_assertions)]
        debug_assert_is_finite(y);

        y
    }
    /// Resets each of the `N_CHANNELS` states to its initial values using
    /// the corresponding initial input value in the `x0` array.
    ///
    /// The corresponding initial output values are written into the `y0` array, if not None.
    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [SaturState],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.reset_state(&mut states[channel], x0[channel]);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut states[channel], x0[channel]);
            }),
        }
    }

    // Not implemented yet: C version only contained assertions
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Triggers control-rate update of coefficients.
    // #[inline(always)]
    // pub fn update_coeffs_ctrl(&mut self, value: f32) {
    //     todo!()
    // }
    /// Triggers audio-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.do_update_coeffs(false);
    }
    /// Processes a single input sample `x`, updating the provided `state`.
    /// Assumes that gain compensation is disabled.
    ///
    /// # Returns
    /// The corresponding output sample.
    /// The actual gain compensation parameter value is ignored.
    #[inline(always)]
    pub fn process1(&mut self, state: &mut SaturState, mut x: f32) -> f32 {
        #[cfg(debug_assertions)]
        debug_assert_is_finite(x);
        x = self.smooth_gain_state.get_y_z1() * x + self.smooth_bias_state.get_y_z1();
        let ax = x.abs();
        // let f = if ax >= 2.115287308554551  { ax - 0.6847736211329452 } else { ax * ax * ((0.00304518315009429 * ax - 0.09167437770414569) * ax + 0.5) };
        let f = if ax >= 2.115_287_3 {
            ax - 0.684_773_6
        } else {
            ax * ax * ((0.003_045_183_1 * ax - 0.091_674_38) * ax + 0.5)
        };
        let d = x - state.x_z1;
        let yb = if d * d < 1e-6 {
            Self::tanhf(0.5 * (x + state.x_z1))
        } else {
            (f - state.f_z1) * rcpf(d)
        };
        let y = yb - self.bias_dc;
        state.x_z1 = x;
        state.f_z1 = f;

        #[cfg(debug_assertions)]
        debug_assert_is_finite(y);

        y
    }
    /// Processes a single input sample `x`, updating the provided `state`.
    /// Assumes that gain compensation is enabled.
    ///
    /// # Returns
    /// The corresponding output sample.
    /// The actual gain compensation parameter value is ignored.
    #[inline(always)]
    pub fn process1_comp(&mut self, state: &mut SaturState, x: f32) -> f32 {
        #[cfg(debug_assertions)]
        debug_assert_is_finite(x);

        let y = self.inv_gain * self.process1(state, x);

        #[cfg(debug_assertions)]
        debug_assert_is_finite(y);

        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and fills the first
    /// `n_samples` of the output buffer `y`, while using and updating both coeffs
    ///  and state (control and audio rate).
    #[inline(always)]
    pub fn process(&mut self, state: &mut SaturState, x: &[f32], y: &mut [f32], n_samples: usize) {
        if self.gain_compensation {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                y[sample] = self.process1_comp(state, x[sample])
            });
        } else {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                y[sample] = self.process1(state, x[sample])
            });
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x`
    /// and fills the first `n_samples` of the `N_CHANNELS` output buffers `y`,
    /// while using and updating both the common coeffs and each of the `N_CHANNELS`
    /// states (control and audio rate).
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [SaturState],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        // Not implemented yet: C version only contained assertions
        // self.update_coeffs_ctrl()

        if self.gain_compensation {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel][sample] =
                        self.process1_comp(&mut states[channel], x[channel][sample])
                });
            });
        } else {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel][sample] = self.process1(&mut states[channel], x[channel][sample])
                });
            });
        }
    }
    /// Sets the input bias value.
    ///
    /// Valid range: [-1e12, 1e12].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(-1e12..=1e12, value);
        }

        self.bias = value;
    }
    /// Sets the gain value.
    ///
    /// Valid range: [1e-12, 1e12].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(-1e12..=1e12, value);
        }

        self.gain = value;
    }
    /// Sets whether the output should be divided by gain (`true`) or not (`false`).
    ///
    /// Default value: `false` (off).
    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        self.gain_compensation = value;
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    /// Tries to determine whether coeffs is valid and returns true if it seems to
    /// be the case and `false` if it is certainly not. False positives are possible,
    /// false negatives are not.
    ///
    /// # Note
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn coeffs_is_valid() -> bool {
        todo!()
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    /// Tries to determine whether state is valid and returns `true` if it seems to
    /// be the case and `false` if it is certainly not. False positives are possible,
    /// false negatives are not.
    ///
    /// # Note
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn state_is_valid() -> bool {
        todo!()
    }

    // Private
    #[inline(always)]
    fn do_update_coeffs(&mut self, force: bool) {
        let mut bias_cur = self.smooth_bias_state.get_y_z1();
        if force || self.bias != bias_cur {
            bias_cur = self
                .smooth_coeffs
                .process1_sticky_abs(&mut self.smooth_bias_state, self.bias);
            self.bias_dc = Self::tanhf(bias_cur)
        }
        let mut gain_cur = self.smooth_gain_state.get_y_z1();
        if force || self.gain != gain_cur {
            gain_cur = self
                .smooth_coeffs
                .process1_sticky_rel(&mut self.smooth_gain_state, self.gain);
            self.inv_gain = rcpf(gain_cur);
        }
    }

    #[inline(always)]
    fn tanhf(x: f32) -> f32 {
        // let xm = clipf(x, -2.115287308554551, 2.115287308554551);
        let xm = clipf(x, -2.115_287_3, 2.115_287_3);
        let axm = xm.abs();
        // xm * axm * (0.01218073260037716 * axm - 0.2750231331124371) + xm;c
        xm * axm * (0.012_180_733 * axm - 0.275_023_13) + xm
    }
}

impl<const N_CHANNELS: usize> Default for SaturCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Internal state and related.
#[derive(Default, Clone, Copy, Debug, PartialEq)]
pub struct SaturState {
    x_z1: f32,
    f_z1: f32,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::{
        c_wrapper::{
            bw_satur_coeffs as SaturCoeffsWrapper, bw_satur_state, satur::Satur as SaturWrapper,
        },
        native::one_pole::tests::{assert_one_pole_coeffs, assert_one_pole_state},
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

        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);

        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn reset_none() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);
        rust_satur.reset(0.0, None);
        c_satur.reset(0.0, None);

        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn reset_some() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_satur.reset(1.0, Some(&mut rust_y0));
        c_satur.reset(1.0, Some(&mut c_y0));

        assert_satur(&rust_satur, &c_satur);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn reset_multi() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.5, 0.5];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_satur.reset_multi(&x0, Some(&mut rust_y0));
        c_satur.reset_multi(&x0, Some(&mut c_y0));

        assert_satur(&rust_satur, &c_satur);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();
        let bias = 1.3;
        let gain = 0.98;

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);

        rust_satur.set_bias(bias);
        rust_satur.set_gain(gain);

        c_satur.set_bias(bias);
        c_satur.set_gain(gain);

        rust_satur.reset(0.0, None);
        c_satur.reset(0.0, None);

        let y_ch: Box<dyn Fn() -> [f32; 8]> = Box::new(|| std::array::from_fn(|_| 0.0));

        let mut rust_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];
        let mut c_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];

        rust_satur.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_satur.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_satur(&rust_satur, &c_satur);
        assert_eq!(rust_y, c_y);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_satur.reset(0.0, Some(&mut rust_y0));
        c_satur.reset(0.0, Some(&mut c_y0));

        assert_satur(&rust_satur, &c_satur);
        assert_eq!(rust_y0, c_y0);

        rust_satur.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_satur.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);
        assert_satur(&rust_satur, &c_satur);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    fn set_bias() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();
        let bias = 1.3;

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);
        rust_satur.set_bias(bias);
        c_satur.set_bias(bias);

        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn set_gain() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();
        let gain = 1.1;

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);
        rust_satur.set_gain(gain);
        c_satur.set_gain(gain);

        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn set_gain_compensation() {
        let mut rust_satur = SaturT::new();
        let mut c_satur = SaturWrapperT::new();

        rust_satur.set_sample_rate(SAMPLE_RATE);
        c_satur.set_sample_rate(SAMPLE_RATE);
        rust_satur.set_gain_compensation(true);
        c_satur.set_gain_compensation(true);
        assert_satur(&rust_satur, &c_satur);

        rust_satur.set_gain_compensation(false);
        c_satur.set_gain_compensation(false);
        assert_satur(&rust_satur, &c_satur);
    }

    #[test]
    fn tanhf() {
        let x = [1.2, 0.4, -0.45, 0.08746328];
        x.iter().for_each(|x| {
            assert_eq!(
                SaturCoeffs::<N_CHANNELS>::tanhf(*x),
                SaturCoeffsWrapper::tanhf(*x)
            )
        });
    }

    fn assert_satur<const N_CHANNELS: usize>(
        rust_satur: &Satur<N_CHANNELS>,
        c_satur: &SaturWrapper<N_CHANNELS>,
    ) {
        assert_satur_coeffs(&rust_satur.coeffs, &c_satur.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_satur_state(&rust_satur.states[channel], &c_satur.states[channel]);
        });
    }

    pub(crate) fn assert_satur_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &SaturCoeffs<N_CHANNELS>,
        c_coeffs: &SaturCoeffsWrapper,
    ) {
        // Sub-components
        assert_one_pole_coeffs(&rust_coeffs.smooth_coeffs, &c_coeffs.smooth_coeffs);
        assert_one_pole_state(&rust_coeffs.smooth_bias_state, &c_coeffs.smooth_bias_state);
        assert_one_pole_state(&rust_coeffs.smooth_gain_state, &c_coeffs.smooth_gain_state);
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

    pub(crate) fn assert_satur_state(rust_state: &SaturState, c_state: &bw_satur_state) {
        assert_eq!(rust_state.f_z1, c_state.F_z1);
        assert_eq!(rust_state.x_z1, c_state.x_z1);
    }
}
