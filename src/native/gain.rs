//! **Smoothed gain** module with optional sticky gain-reach threshold.
//!
//! # Example
//! ```rust
//! use brickworks_rs::native::gain::Gain;
//! use brickworks_rs::native::one_pole::StickyMode;
//! const N_CHANNELS: usize = 2;
//! const N_SAMPLES: usize = 4;
//! const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
//!         &[1.0, 0.0, 0.0, 0.0],
//!         &[1.0, 0.0, 0.0, 0.0]
//!     ];
//! let mut y: [&mut[f32]; N_CHANNELS] = [
//!         &mut [0.0, 0.0, 0.0, 0.0],
//!         &mut [0.0, 0.0, 0.0, 0.0],
//!     ];
//! let mut gain = Gain::<N_CHANNELS>::new();
//! gain.set_sample_rate(44_100.0);
//! gain.set_gain_lin(0.9);
//! gain.set_smooth_tau(0.2);
//! gain.set_sticky_thresh(1.2);
//! gain.set_sticky_mode(StickyMode::Abs);
//!
//! gain.reset();
//! gain.process(&PULSE_INPUT, &mut y, N_SAMPLES);
//!
//! println!("Output: {:?}", y);
//! ```
//! //!
//! # Notes
//! This module provides a native Rust implementation of the filter, but the same interface is
//! also available via bindings to the original C library at [crate::c_wrapper::gain].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_gain).
#[cfg(debug_assertions)]
use super::common::{debug_assert_positive, debug_assert_range};
use crate::native::{
    math::db2linf,
    one_pole::{OnePoleCoeffs, OnePoleState, StickyMode},
};
/// **Smoothed gain** module with optional sticky gain-reach threshold.
/// This struct manages the gain coefficients for a given number of channels
/// (`N_CHANNELS`).  
///
/// It wraps:
/// - [`GainCoeffs`]
///
/// # Usage
/// ```rust
/// use brickworks_rs::native::gain::Gain;
/// const N_CHANNELS: usize = 2;
/// let mut gain = Gain::<N_CHANNELS>::new();
/// gain.set_sample_rate(48_000.0);
/// gain.set_gain_db(40.0);
/// gain.reset();
/// // process audio with gain.process(...)
/// ```
pub struct Gain<const N_CHANNELS: usize> {
    coeffs: GainCoeffs<N_CHANNELS>,
}

impl<const N_CHANNELS: usize> Gain<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            coeffs: GainCoeffs::new(),
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset(&mut self) {
        self.coeffs.reset_coeffs();
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and
    /// writes the results to the first `n_samples` of the `N_CHANNELS` output
    /// buffers `y`, while updating the shared coefficients (control and audio rate).
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.coeffs.process_multi(x, y, n_samples);
    }
    /// Sets the gain parameter to the given value (linear gain) in coeffs.
    ///
    /// `value` must be finite.
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain_lin(&mut self, value: f32) {
        self.coeffs.set_gain_lin(value);
    }
    /// Sets the gain parameter to the given value (dB) in coeffs.
    ///
    /// `value` must be less than or equal to 770.630.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_gain_db(&mut self, value: f32) {
        self.coeffs.set_gain_db(value);
    }
    /// Sets the smoothing time constant value (s) in coeffs.
    ///
    /// `value` must be non-negative.
    ///
    /// Default value: 0.05.
    #[inline(always)]
    pub fn set_smooth_tau(&mut self, value: f32) {
        self.coeffs.set_smooth_tau(value);
    }
    /// Sets the gain-reach threshold specified by value in coeffs.
    ///
    /// When the difference between the current and the target gain would
    /// fall under such threshold according to the current distance metric
    /// (see [Gain::set_sticky_mode()]), the current gain is forcefully
    /// set to be equal to the target gain value.
    ///
    /// Valid range: [0.0, 1e18].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        self.coeffs.set_sticky_thresh(value);
    }
    /// Sets the current distance metric for sticky behavior to value in coeffs.
    ///
    /// *   [StickyMode::Abs]: absolute gain difference (|current - target|, with
    ///     current and target linear);
    /// *   [StickyMode::Rel]: relative gain difference with respect to target gain
    ///     (|current - target| / |target|, with current and target linear).
    ///
    /// Default value: [StickyMode::Abs].
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: StickyMode) {
        self.coeffs.set_sticky_mode(value);
    }
    /// Returns the current gain parameter value (linear gain) in coeffs.
    #[inline(always)]
    pub fn get_gain_lin(&self) -> f32 {
        self.coeffs.get_gain_lin()
    }
    /// Returns the actual current gain coefficient (linear gain) in coeffs.
    ///
    /// coeffs must be at least in the "reset" state.
    #[inline(always)]
    pub fn get_gain_cur(&self) -> f32 {
        self.coeffs.get_gain_cur()
    }
}

impl<const N_CHANNELS: usize> Default for Gain<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

/// Coefficients and related.
pub struct GainCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    smooth_state: OnePoleState,
    // Parameters
    gain: f32,
}

impl<const N_CHANNELS: usize> GainCoeffs<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut smooth_coeffs = OnePoleCoeffs::new();
        smooth_coeffs.set_tau(0.05);
        Self {
            smooth_coeffs,
            smooth_state: OnePoleState::new(),
            gain: 1.0,
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(
                sample_rate.is_finite(),
                "value must be finite, got {}",
                sample_rate
            );
            debug_assert_positive(sample_rate);
        }
        self.smooth_coeffs.set_sample_rate(sample_rate);
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.smooth_coeffs.reset_coeffs();
        self.smooth_coeffs
            .reset_state(&mut self.smooth_state, self.gain);
    }
    /// Triggers control-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        self.smooth_coeffs.update_coeffs_ctrl();
    }
    /// Triggers audio-rate update of coefficients.
    ///
    /// Assumes that the gain-reach threshold is 0.0.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        // OnePole::update_coeffs_audio() is not implemented yet:
        // C version only contained assertions need to revisit
        // which assertions from the C version make sense to keep in Rust
        // self.smooth_coeffs.update_coeffs_audio();

        self.smooth_coeffs
            .process1(&mut self.smooth_state, self.gain);
    }
    /// Triggers audio-rate update of coefficients.
    ///
    /// Assumes that the gain-reach threshold is not 0.0 and the distance
    /// metric for sticky behavior is set to StickyMode::Abs.
    #[inline(always)]
    pub fn update_coeffs_audio_sticky_abs(&mut self) {
        // OnePole::update_coeffs_audio() is not implemented yet:
        // C version only contained assertions need to revisit
        // which assertions from the C version make sense to keep in Rust
        // self.smooth_coeffs.update_coeffs_audio();

        self.smooth_coeffs
            .process1_sticky_abs(&mut self.smooth_state, self.gain);
    }
    /// Triggers audio-rate update of coefficients.
    ///
    /// Assumes that the gain-reach threshold is not 0.0 and the distance
    /// metric for sticky behavior is set to StickyMode::Rel.
    #[inline(always)]
    pub fn update_coeffs_audio_sticky_rel(&mut self) {
        // OnePole::update_coeffs_audio() is not implemented yet:
        // C version only contained assertions need to revisit
        // which assertions from the C version make sense to keep in Rust
        // self.smooth_coeffs.update_coeffs_audio();

        self.smooth_coeffs
            .process1_sticky_rel(&mut self.smooth_state, self.gain);
    }
    /// Processes one input sample `x` using coeffs and returns the corresponding
    /// output sample.
    #[inline(always)]
    pub fn process1(&mut self, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let y = self.smooth_state.get_y_z1() * x;

        debug_assert!(y.is_finite());

        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and fills the first
    /// `n_samples` of the output buffer `y`, while using and updating coeffs
    /// (control and audio rate).
    #[inline(always)]
    pub fn process(&mut self, x: &[f32], y: &mut [f32], n_samples: usize) {
        self.update_coeffs_ctrl();
        if self.smooth_coeffs.get_sticky_thresh() == 0.0 {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                y[sample] = self.process1(x[sample]);
            });
        } else {
            match self.smooth_coeffs.get_sticky_mode() {
                StickyMode::Abs => {
                    (0..n_samples).for_each(|sample| {
                        self.update_coeffs_audio_sticky_abs();
                        y[sample] = self.process1(x[sample]);
                    });
                }
                StickyMode::Rel => {
                    (0..n_samples).for_each(|sample| {
                        self.update_coeffs_audio_sticky_rel();
                        y[sample] = self.process1(x[sample]);
                    });
                }
            }
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and
    /// writes the results to the first `n_samples` of the `N_CHANNELS` output
    /// buffers `y`, while updating the shared coefficients (control and audio rate).
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.update_coeffs_ctrl();
        if self.smooth_coeffs.get_sticky_thresh() == 0.0 {
            (0..n_samples).for_each(|sample| {
                self.update_coeffs_audio();
                (0..N_CHANNELS).for_each(|channel| {
                    y[channel][sample] = self.process1(x[channel][sample]);
                });
            });
        } else {
            match self.smooth_coeffs.get_sticky_mode() {
                StickyMode::Abs => {
                    (0..n_samples).for_each(|sample| {
                        self.update_coeffs_audio_sticky_abs();
                        (0..N_CHANNELS).for_each(|channel| {
                            y[channel][sample] = self.process1(x[channel][sample]);
                        });
                    });
                }
                StickyMode::Rel => {
                    (0..n_samples).for_each(|sample| {
                        self.update_coeffs_audio_sticky_rel();
                        (0..N_CHANNELS).for_each(|channel| {
                            y[channel][sample] = self.process1(x[channel][sample]);
                        });
                    });
                }
            }
        }
    }
    /// Sets the gain parameter to the given value (linear gain) in coeffs.
    ///
    /// `value` must be finite.
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain_lin(&mut self, value: f32) {
        debug_assert!(value.is_finite(), "value must be finite, got {}", value);

        self.gain = value;
    }
    /// Sets the gain parameter to the given value (dB) in coeffs.
    ///
    /// `value` must be less than or equal to 770.630.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_gain_db(&mut self, value: f32) {
        debug_assert!(
            value <= 770.630,
            "value must be less or equal to 770.630, got {}",
            value
        );
        debug_assert!(!value.is_nan());

        self.gain = db2linf(value);
    }
    /// Sets the smoothing time constant value (s) in coeffs.
    ///
    /// `value` must be non-negative.
    ///
    /// Default value: 0.05.
    #[inline(always)]
    pub fn set_smooth_tau(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(!value.is_nan());
            debug_assert_positive(value);
        }

        self.smooth_coeffs.set_tau(value);
    }
    /// Sets the gain-reach threshold specified by value in coeffs.
    ///
    /// When the difference between the current and the target gain would
    /// fall under such threshold according to the current distance metric
    /// (see [Gain::set_sticky_mode()]), the current gain is forcefully
    /// set to be equal to the target gain value.
    ///
    /// Valid range: [0.0, 1e18].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(!value.is_nan());
            debug_assert_range(0.0..=1e18, value);
        }

        self.smooth_coeffs.set_sticky_thresh(value);
    }
    /// Sets the current distance metric for sticky behavior to value in coeffs.
    ///
    /// *   [StickyMode::Abs]: absolute gain difference (|current - target|, with
    ///     current and target linear);
    /// *   [StickyMode::Rel]: relative gain difference with respect to target gain
    ///     (|current - target| / |target|, with current and target linear).
    ///
    /// Default value: [StickyMode::Abs].
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: StickyMode) {
        self.smooth_coeffs.set_sticky_mode(value);
    }
    /// Returns the current gain parameter value (linear gain) in coeffs.
    #[inline(always)]
    pub fn get_gain_lin(&self) -> f32 {
        self.gain
    }
    /// Returns the actual current gain coefficient (linear gain) in coeffs.
    ///
    /// coeffs must be at least in the "reset" state.
    #[inline(always)]
    pub fn get_gain_cur(&self) -> f32 {
        self.smooth_state.get_y_z1()
    }
    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    /// Tries to determine whether coeffs is valid and returns true if it seems to
    /// be the case and false if it is certainly not. False positives are possible,
    /// false negatives are not.
    /// # Notes
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn coeffs_is_valid() -> bool {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for GainCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_gain_coeffs, gain::Gain as GainWrapper},
        native::one_pole::tests::{assert_one_pole_coeffs, assert_one_pole_state},
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [&[1.0, 1.0], &[0.0, 0.0]];
    const N_SAMPLES: usize = 2;

    type GainT = Gain<N_CHANNELS>;
    type GainWrapperT = GainWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_gain = GainT::new();
        let c_gain = GainWrapperT::new();

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut rust_gain = GainT::new();
        let mut c_gain = GainWrapperT::new();

        rust_gain.set_sample_rate(SAMPLE_RATE);
        c_gain.set_sample_rate(SAMPLE_RATE);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_sample_rate_invalid() {
        let mut rust_gain = GainT::new();
        rust_gain.set_sample_rate(f32::INFINITY);
    }

    #[test]
    fn reset() {
        const TAU: f32 = 0.00005;
        let sticky_tresh = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_smooth_tau(TAU);
        c_gain.set_sticky_thresh(sticky_tresh);
        c_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        rust_gain.set_sticky_thresh(sticky_tresh);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        c_gain.reset();
        rust_gain.reset();

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn process() {
        const TAU: f32 = 0.00005;
        let sticky_tresh = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_smooth_tau(TAU);
        c_gain.set_sticky_thresh(sticky_tresh);
        c_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        rust_gain.set_sticky_thresh(sticky_tresh);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        c_gain.reset();
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);
        rust_gain.reset();
        rust_gain.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);

        assert_gain(&rust_gain, &c_gain);

        (0..N_CHANNELS).for_each(|channel| {
            (0..N_SAMPLES).for_each(|sample| {
                assert_eq!(rust_y[channel][sample], c_y[channel][sample]);
            });
        });
    }

    #[test]
    fn set_gain_lin_valid() {
        const GAIN_LIN: f32 = 0.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_LIN);
        c_gain.set_gain_lin(GAIN_LIN);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_gain_lin_invalid() {
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_gain_lin(f32::INFINITY);
    }

    #[test]
    fn set_gain_db_valid() {
        const GAIN_DB: f32 = 12.0;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_DB);
        c_gain.set_gain_lin(GAIN_DB);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be less or equal to 770.630, got 770.7")]
    fn set_gain_db_invalid() {
        const GAIN_DB: f32 = 770.7;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_gain_db(GAIN_DB);
    }

    #[test]
    fn set_smooth_tau_valid() {
        const TAU: f32 = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_smooth_tau(TAU);
        c_gain.set_smooth_tau(TAU);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_smooth_tau_invalid() {
        const TAU: f32 = -1.0;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_smooth_tau(TAU);
    }

    #[test]
    fn set_sticky_thresh_valid() {
        const STICKY_THRESH: f32 = 0.1;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_sticky_thresh(STICKY_THRESH);
        c_gain.set_sticky_thresh(STICKY_THRESH);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got -1")]
    fn set_sticky_thresh_invalid() {
        const STICKY_THRESH: f32 = -1.0;
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_sticky_thresh(STICKY_THRESH);
    }

    #[test]
    fn set_sticky_mode() {
        const STICKY_THRESH: f32 = 0.01;

        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_sticky_thresh(STICKY_THRESH);
        c_gain.set_sticky_thresh(STICKY_THRESH);

        rust_gain.set_sticky_mode(crate::native::one_pole::StickyMode::Abs);
        c_gain.set_sticky_mode(crate::c_wrapper::one_pole::StickyMode::Abs);

        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn get_gain_lin() {
        const GAIN_DB: f32 = 3.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_DB);
        c_gain.set_gain_lin(GAIN_DB);

        rust_gain.reset();
        c_gain.reset();

        rust_gain.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        rust_gain.reset();
        c_gain.reset();

        assert_eq!(&rust_gain.get_gain_lin(), &c_gain.get_gain_lin());
        assert_gain(&rust_gain, &c_gain);
    }

    #[test]
    fn get_gain_cur() {
        const GAIN_LIN: f32 = 0.5;
        let mut c_gain = GainWrapperT::new();
        let mut rust_gain = GainT::new();

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];

        c_gain.set_sample_rate(SAMPLE_RATE);
        rust_gain.set_sample_rate(SAMPLE_RATE);

        rust_gain.set_gain_lin(GAIN_LIN);
        c_gain.set_gain_lin(GAIN_LIN);

        rust_gain.reset();
        c_gain.reset();

        rust_gain.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_gain.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        rust_gain.reset();
        c_gain.reset();

        assert_eq!(&rust_gain.get_gain_cur(), &c_gain.get_gain_cur());
        assert_gain(&rust_gain, &c_gain);
    }

    fn assert_gain<const N_CHANNELS: usize>(
        rust_gain: &Gain<N_CHANNELS>,
        c_gain: &GainWrapper<N_CHANNELS>,
    ) {
        assert_gain_coeffs(&rust_gain.coeffs, &c_gain.coeffs);
    }

    pub(crate) fn assert_gain_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &GainCoeffs<N_CHANNELS>,
        c_coeffs: &bw_gain_coeffs,
    ) {
        assert_one_pole_coeffs(&rust_coeffs.smooth_coeffs, &c_coeffs.smooth_coeffs);
        assert_one_pole_state(&rust_coeffs.smooth_state, &c_coeffs.smooth_state);
        assert_eq!(rust_coeffs.gain, c_coeffs.gain);
    }
}
