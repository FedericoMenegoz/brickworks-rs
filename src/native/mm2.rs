//! **Second-order multimode filter**.
//!
//! It implements an approximation of the Laplace-domain transfer function
//! ```text
//! H(s) = coeff_x + (coeff_hp s^2 + 2 pi fc s coeff_bp + (2 pi fc)^2 coeff_lp) / (s^2 + 2 pi fc / Q s + (2 pi fc)^2)
//! ```
//!
//! where `fc` is the cutoff frequency and `Q` is the quality factor.
//!
//! ## Example
//!
//! ```rust
//! use brickworks_rs::native::mm2::MM2;
//!
//! const N_CHANNELS: usize = 2;
//! const N_SAMPLES: usize = 8;
//!
//! fn main() {
//!     // Create a stereo MM2 filter
//!     let mut mm2 = MM2::<N_CHANNELS>::new();
//!
//!     // Configure sample rate
//!     mm2.set_sample_rate(44_100.0);
//!
//!     // Set filter parameters
//!     mm2.set_cutoff(2_000.0);   // cutoff frequency in Hz
//!     mm2.set_q(0.707);          // Q factor
//!
//!     // Mix coefficients
//!     mm2.set_coeff_lp(0.9);
//!     mm2.set_coeff_x(0.3);
//!
//!     // Reset states with initial input = silence
//!     mm2.reset(0.0, None);
//!
//!     // Example input: stereo pulse at first sample
//!     let input: [&[f32]; N_CHANNELS] = [
//!         &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!         &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!     ];
//!
//!     // Output buffers
//!     let mut output_l = [0.0; N_SAMPLES];
//!     let mut output_r = [0.0; N_SAMPLES];
//!     let mut outputs: [&mut [f32]; N_CHANNELS] = [&mut output_l, &mut output_r];
//!
//!     // Process audio block
//!     mm2.process(&input, &mut outputs, N_SAMPLES);
//!
//!     println!("Left channel:  {:?}", outputs[0]);
//!     println!("Right channel: {:?}", outputs[1]);
//! }
//! ```
//! ## Notes
//! This module provides a native Rust implementation, but the same interface is
//! also available via bindings to the original C library at [crate::c_wrapper::mm2].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_mm2).
use crate::native::{
    gain::GainCoeffs,
    svf::{SVFCoeffs, SVFState},
};

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_positive, debug_assert_range};
/// Multi–Mode 2–pole filter (`MM2`) with per–channel state and coefficients.
///
/// Manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`MM2Coeffs`]
/// - [`MM2State`]
///
/// # Usage
///
/// ```rust
/// use brickworks_rs::native::mm2::MM2;
///
/// // Stereo filter
/// let mut mm2 = MM2::<2>::new();
///
/// mm2.set_sample_rate(44_100.0);
/// mm2.set_cutoff(2_000.0);
/// mm2.set_q(0.707);
/// mm2.set_coeff_x(0.8);
/// mm2.set_coeff_bp(0.2);
/// mm2.reset(0.0, None);
/// // mm2.process(...)
/// ```
pub struct MM2<const N_CHANNELS: usize> {
    pub(crate) coeffs: MM2Coeffs<N_CHANNELS>,
    pub(crate) states: [MM2State; N_CHANNELS],
}

impl<const N_CHANNELS: usize> MM2<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let coeffs = MM2Coeffs::<N_CHANNELS>::new();
        Self {
            coeffs,
            states: [MM2State::default(); N_CHANNELS],
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
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.coeffs.reset_state(&mut self.states[channel], x0);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.coeffs.reset_state(&mut self.states[channel], x0);
            }),
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
    /// Sets the cutoff frequency value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        self.coeffs.set_cutoff(value);
    }
    /// Sets the quality factor to the given value.
    ///
    /// Valid range: [1e-6, 1e6].
    ///
    /// Default value: 0.5.
    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        self.coeffs.set_q(value);
    }
    /// Sets whether bilinear transform prewarping frequency should match the
    /// cutoff frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.coeffs.set_prewarp_at_cutoff(value);
    }
    /// Sets the prewarping frequency value (Hz).
    ///
    /// Only used when the `prewarp_at_cutoff` parameter is off and however
    /// internally limited to avoid instability.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        self.coeffs.set_prewarp_freq(value);
    }
    /// Sets the input mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_coeff_x(&mut self, value: f32) {
        self.coeffs.set_coeff_x(value);
    }
    /// Sets the lowpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_lp(&mut self, value: f32) {
        self.coeffs.set_coeff_lp(value);
    }
    /// Sets the bandpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_bp(&mut self, value: f32) {
        self.coeffs.set_coeff_bp(value);
    }
    /// Sets the highpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_hp(&mut self, value: f32) {
        self.coeffs.set_coeff_hp(value);
    }
}

impl<const N_CHANNELS: usize> Default for MM2<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
pub struct MM2Coeffs<const N_CHANNELS: usize> {
    // Sub-components
    svf_coeffs: SVFCoeffs<N_CHANNELS>,
    gain_x_coeffs: GainCoeffs<N_CHANNELS>,
    gain_lp_coeffs: GainCoeffs<N_CHANNELS>,
    gain_bp_coeffs: GainCoeffs<N_CHANNELS>,
    gain_hp_coeffs: GainCoeffs<N_CHANNELS>,
}

impl<const N_CHANNELS: usize> MM2Coeffs<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let svf_coeffs = SVFCoeffs::<N_CHANNELS>::new();
        let mut gain_x_coeffs = GainCoeffs::<N_CHANNELS>::new();
        let mut gain_lp_coeffs = GainCoeffs::<N_CHANNELS>::new();
        let mut gain_bp_coeffs = GainCoeffs::<N_CHANNELS>::new();
        let mut gain_hp_coeffs = GainCoeffs::<N_CHANNELS>::new();

        let default_tau = 0.005;
        let x_gain = 1.0;
        let passes_gain = 0.0;

        gain_x_coeffs.set_smooth_tau(default_tau);
        gain_lp_coeffs.set_smooth_tau(default_tau);
        gain_bp_coeffs.set_smooth_tau(default_tau);
        gain_hp_coeffs.set_smooth_tau(default_tau);
        gain_x_coeffs.set_gain_lin(x_gain);
        gain_lp_coeffs.set_gain_lin(passes_gain);
        gain_bp_coeffs.set_gain_lin(passes_gain);
        gain_hp_coeffs.set_gain_lin(passes_gain);

        Self {
            svf_coeffs,
            gain_x_coeffs,
            gain_lp_coeffs,
            gain_bp_coeffs,
            gain_hp_coeffs,
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(sample_rate.is_finite());
            debug_assert_positive(sample_rate);
        }

        self.svf_coeffs.set_sample_rate(sample_rate);
        self.gain_x_coeffs.set_sample_rate(sample_rate);
        self.gain_lp_coeffs.set_sample_rate(sample_rate);
        self.gain_bp_coeffs.set_sample_rate(sample_rate);
        self.gain_hp_coeffs.set_sample_rate(sample_rate);
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.svf_coeffs.reset_coeffs();
        self.gain_x_coeffs.reset_coeffs();
        self.gain_lp_coeffs.reset_coeffs();
        self.gain_bp_coeffs.reset_coeffs();
        self.gain_hp_coeffs.reset_coeffs();
    }
    /// Resets the given state to its initial values using the initial
    /// input value `x0`.
    ///
    /// Returns the corresponding initial output value.
    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut MM2State, x0: f32) -> f32 {
        let (mut lp, mut bp, mut hp) = (0.0, 0.0, 0.0);
        self.svf_coeffs
            .reset_state(&mut state.svf_state, x0, &mut lp, &mut bp, &mut hp);

        let y = self.gain_x_coeffs.get_gain_cur() * x0
            + self.gain_lp_coeffs.get_gain_cur() * lp
            + self.gain_bp_coeffs.get_gain_cur() * bp
            + self.gain_hp_coeffs.get_gain_cur() * hp;

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
        state: &mut [MM2State],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.reset_state(&mut state[channel], x0[channel]);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut state[channel], x0[channel]);
            }),
        }
    }
    /// Triggers control-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        // Not implemented yet as only asserting:
        // see native::svf::SVFCoeffs.update_coeff_ctrl()
        // self.svf_coeffs.update_coeffs_ctrl();

        self.gain_x_coeffs.update_coeffs_ctrl();
        self.gain_lp_coeffs.update_coeffs_ctrl();
        self.gain_bp_coeffs.update_coeffs_ctrl();
        self.gain_hp_coeffs.update_coeffs_ctrl();
    }
    /// Triggers audio-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.svf_coeffs.update_coeffs_audio();
        self.gain_x_coeffs.update_coeffs_audio();
        self.gain_lp_coeffs.update_coeffs_audio();
        self.gain_bp_coeffs.update_coeffs_audio();
        self.gain_hp_coeffs.update_coeffs_audio();
    }
    /// Processes one input sample `x`, while using and updating state.
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1(&mut self, state: &mut MM2State, x: f32) -> f32 {
        debug_assert!(x.is_finite(), "value must be finite, got {}", x);
        let (mut lp, mut bp, mut hp) = (0.0, 0.0, 0.0);

        self.svf_coeffs
            .process1(&mut state.svf_state, x, &mut lp, &mut bp, &mut hp);

        let vx = self.gain_x_coeffs.process1(x);
        let vlp = self.gain_lp_coeffs.process1(lp);
        let vbp = self.gain_bp_coeffs.process1(bp);
        let vhp = self.gain_hp_coeffs.process1(hp);
        let y = vx + vlp + vbp + vhp;

        debug_assert!(y.is_finite());

        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and fills the first
    /// `n_samples` of the output buffer `y`, while using and updating both coeffs
    /// and state (control and audio rate).
    #[inline(always)]
    pub fn process(&mut self, state: &mut MM2State, x: &[f32], y: &mut [f32], n_samples: usize) {
        self.update_coeffs_ctrl();

        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample]);
        });
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and fills the
    /// first `n_samples` of the `N_CHANNELS` output buffers `y`, while using and updating
    /// both the common coeffs and each of the `N_CHANNELS` states (control and audio rate).
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [MM2State],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.update_coeffs_ctrl();

        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            (0..N_CHANNELS).for_each(|channel| {
                y[channel][sample] = self.process1(&mut states[channel], x[channel][sample]);
            });
        });
    }
    /// Sets the cutoff frequency value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite(), "value must be finite, got {}", value);
            debug_assert_range(1e-6..=1e12, value);
        }

        self.svf_coeffs.set_cutoff(value);
    }
    /// Sets the quality factor to the given value.
    ///
    /// Valid range: [1e-6, 1e6].
    ///
    /// Default value: 0.5.
    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite(), "value must be finite, got {}", value);
            debug_assert_range(1e-6..=1e6, value);
        }

        self.svf_coeffs.set_q(value);
    }
    /// Sets whether bilinear transform prewarping frequency should match the
    /// cutoff frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.svf_coeffs.set_prewarp_at_cutoff(value);
    }
    /// Sets the prewarping frequency value (Hz).
    ///
    /// Only used when the `prewarp_at_cutoff` parameter is off and however
    /// internally limited to avoid instability.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite(), "value must be finite, got {}", value);
            debug_assert_range(1e-6..=1e12, value);
        }

        self.svf_coeffs.set_prewarp_freq(value);
    }
    /// Sets the input mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_coeff_x(&mut self, value: f32) {
        debug_assert!(value.is_finite(), "value must be finite, got {}", value);

        self.gain_x_coeffs.set_gain_lin(value);
    }
    /// Sets the lowpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_lp(&mut self, value: f32) {
        debug_assert!(value.is_finite(), "value must be finite, got {}", value);

        self.gain_lp_coeffs.set_gain_lin(value);
    }
    /// Sets the bandpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_bp(&mut self, value: f32) {
        debug_assert!(value.is_finite(), "value must be finite, got {}", value);

        self.gain_bp_coeffs.set_gain_lin(value);
    }
    /// Sets the highpass mode coefficient value.
    ///
    /// value must be finite.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_coeff_hp(&mut self, value: f32) {
        debug_assert!(value.is_finite(), "value must be finite, got {}", value);

        self.gain_hp_coeffs.set_gain_lin(value);
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
        todo!();
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
        todo!();
    }
}

impl<const N_CHANNELS: usize> Default for MM2Coeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

/// Internal state and related.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct MM2State {
    // Sub-components
    svf_state: SVFState,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::{
        c_wrapper::{bw_mm2_coeffs, bw_mm2_state, mm2::MM2 as MM2Wrapper},
        native::{
            gain::tests::assert_gain_coeffs,
            svf::tests::{assert_svf_coeffs, assert_svf_states},
        },
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type MM2T = MM2<N_CHANNELS>;
    type MM2WrapperT = MM2Wrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_mm2 = MM2T::new();
        let c_mm2 = MM2WrapperT::new();

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn reset_none() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;

        rust_mm2.reset(x0, None);
        c_mm2.reset(x0, None);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn reset_some() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_mm2.reset(x0, Some(&mut rust_y0));
        c_mm2.reset(x0, Some(&mut c_y0));

        assert_mm2(&rust_mm2, &c_mm2);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn reset_multi() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.045, 0.045];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_mm2.reset_multi(&x0, Some(&mut rust_y0));
        c_mm2.reset_multi(&x0, Some(&mut c_y0));

        assert_mm2(&rust_mm2, &c_mm2);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];

        rust_mm2.set_cutoff(2000.0);
        rust_mm2.set_coeff_bp(0.7);
        rust_mm2.set_coeff_x(0.5);
        rust_mm2.set_coeff_lp(0.9);

        c_mm2.set_cutoff(2000.0);
        c_mm2.set_coeff_bp(0.7);
        c_mm2.set_coeff_x(0.5);
        c_mm2.set_coeff_lp(0.9);

        rust_mm2.reset(0.0, None);
        c_mm2.reset(0.0, None);

        rust_mm2.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_mm2.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_mm2(&rust_mm2, &c_mm2);
        (0..N_CHANNELS).for_each(|channel| {
            (0..N_SAMPLES).for_each(|sample| {
                assert_eq!(
                    rust_y[channel][sample], c_y[channel][sample],
                    "Mismatch at ch {channel}, sample {sample}"
                );
            });
        });
    }

    #[test]
    fn set_cutoff() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let cutoff = 261.63;

        rust_mm2.set_cutoff(cutoff);
        c_mm2.set_cutoff(cutoff);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_q() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let q = 1.4;

        rust_mm2.set_q(q);
        c_mm2.set_q(q);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_prewarp() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let cutoff = 2000.0;
        rust_mm2.set_cutoff(cutoff);
        c_mm2.set_cutoff(cutoff);

        let prewarp_freq = 2200.0;
        rust_mm2.set_prewarp_freq(prewarp_freq);
        c_mm2.set_prewarp_freq(prewarp_freq);

        rust_mm2.set_prewarp_at_cutoff(false);
        c_mm2.set_prewarp_at_cutoff(false);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_coeff_x() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let coeff_x = 0.5;
        rust_mm2.set_coeff_x(coeff_x);
        c_mm2.set_coeff_x(coeff_x);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_coeff_lp() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let coeff_lp = 0.9;
        rust_mm2.set_coeff_lp(coeff_lp);
        c_mm2.set_coeff_lp(coeff_lp);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_coeff_bp() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let coeff_bp = 0.7;
        rust_mm2.set_coeff_bp(coeff_bp);
        c_mm2.set_coeff_bp(coeff_bp);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    #[test]
    fn set_coeff_hp() {
        let mut rust_mm2 = MM2T::new();
        let mut c_mm2 = MM2WrapperT::new();

        rust_mm2.set_sample_rate(SAMPLE_RATE);
        c_mm2.set_sample_rate(SAMPLE_RATE);

        let coeff_hp = 0.6;
        rust_mm2.set_coeff_hp(coeff_hp);
        c_mm2.set_coeff_hp(coeff_hp);

        assert_mm2(&rust_mm2, &c_mm2);
    }

    fn assert_mm2<const N_CHANNELS: usize>(
        rust_mm2: &MM2<N_CHANNELS>,
        c_mm2: &MM2Wrapper<N_CHANNELS>,
    ) {
        assert_mm2_coeffs(&rust_mm2.coeffs, &c_mm2.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_mm2_state(&rust_mm2.states[channel], &c_mm2.states[channel]);
        });
    }

    pub(crate) fn assert_mm2_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &MM2Coeffs<N_CHANNELS>,
        c_coeffs: &bw_mm2_coeffs,
    ) {
        assert_svf_coeffs(&rust_coeffs.svf_coeffs, &c_coeffs.svf_coeffs);
        assert_gain_coeffs(&rust_coeffs.gain_x_coeffs, &c_coeffs.gain_x_coeffs);
        assert_gain_coeffs(&rust_coeffs.gain_lp_coeffs, &c_coeffs.gain_lp_coeffs);
        assert_gain_coeffs(&rust_coeffs.gain_bp_coeffs, &c_coeffs.gain_bp_coeffs);
        assert_gain_coeffs(&rust_coeffs.gain_hp_coeffs, &c_coeffs.gain_hp_coeffs);
    }

    pub(crate) fn assert_mm2_state(rust_state: &MM2State, c_state: &bw_mm2_state) {
        assert_svf_states(&rust_state.svf_state, &c_state.svf_state);
    }
}
