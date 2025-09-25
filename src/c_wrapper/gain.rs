//! **Smoothed gain** module with optional sticky gain-reach threshold.
//!
//! # Example
//! ```rust
//! use brickworks_rs::c_wrapper::gain::Gain;
//! use brickworks_rs::c_wrapper::one_pole::StickyMode;
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
//!
//! # Notes
//! This module provides Rust bindings to the original C implementation.
//! For a fully native Rust implementation with the same interface,
//! see [crate::native::gain].
//! Original C library by [Orastron](https://www.orastron.com/algorithms/bw_gain).
use super::*;
use crate::c_wrapper::one_pole::StickyMode;
/// **Smoothed gain** module with optional sticky gain-reach threshold.
/// This struct manages the gain coefficients for a given number of channels
/// (`N_CHANNELS`).  
///
/// It wraps:
/// - [`GainCoeffs`]
///
/// # Usage
/// ```rust
/// use brickworks_rs::c_wrapper::gain::Gain;
/// const N_CHANNELS: usize = 2;
/// let mut gain = Gain::<N_CHANNELS>::new();
/// gain.set_sample_rate(48_000.0);
/// gain.set_gain_db(40.0);
/// gain.reset();
/// // process audio with gain.process(...)
/// ```
pub struct Gain<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_gain_coeffs,
}

impl<const N_CHANNELS: usize> Gain<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_gain_coeffs::default();
        unsafe {
            bw_gain_init(&mut coeffs);
        }
        Self { coeffs }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_gain_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset(&mut self) {
        unsafe {
            bw_gain_reset_coeffs(&mut self.coeffs);
        }
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
        let x_p: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_p: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        unsafe {
            bw_gain_process_multi(
                &mut self.coeffs,
                x_p.as_ptr(),
                y_p.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }
    /// Sets the gain parameter to the given value (linear gain) in coeffs.
    ///
    /// `value` must be finite.
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_gain_lin(&mut self, value: f32) {
        unsafe {
            bw_gain_set_gain_lin(&mut self.coeffs, value);
        }
    }
    /// Sets the gain parameter to the given value (dB) in coeffs.
    ///
    /// `value` must be less than or equal to 770.630.
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_gain_db(&mut self, value: f32) {
        unsafe {
            bw_gain_set_gain_dB(&mut self.coeffs, value);
        }
    }
    /// Sets the smoothing time constant value (s) in coeffs.
    ///
    /// `value` must be non-negative.
    ///
    /// Default value: 0.05.
    #[inline(always)]
    pub fn set_smooth_tau(&mut self, value: f32) {
        unsafe {
            bw_gain_set_smooth_tau(&mut self.coeffs, value);
        }
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
        unsafe {
            bw_gain_set_sticky_thresh(&mut self.coeffs, value);
        }
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
        unsafe {
            bw_gain_set_sticky_mode(&mut self.coeffs, value.to_bw());
        }
    }
    /// Returns the current gain parameter value (linear gain) in coeffs.
    #[inline(always)]
    pub fn get_gain_lin(&self) -> f32 {
        unsafe { bw_gain_get_gain_lin(&self.coeffs) }
    }
    /// Returns the actual current gain coefficient (linear gain) in coeffs.
    ///
    /// coeffs must be at least in the "reset" state.
    #[inline(always)]
    pub fn get_gain_cur(&self) -> f32 {
        unsafe { bw_gain_get_gain_cur(&self.coeffs) }
    }
}

impl<const N_CHANNELS: usize> Default for Gain<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        c_wrapper::{
            bw_dB2linf, bw_gain_coeffs, bw_gain_init, bw_gain_process1, bw_gain_reset_coeffs,
            bw_gain_set_sample_rate, bw_gain_set_sticky_thresh, bw_gain_update_coeffs_audio,
            bw_gain_update_coeffs_ctrl, bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs,
            bw_rcpf, gain::Gain, one_pole::StickyMode,
        },
        native::math::INVERSE_2_PI,
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    type GainT = Gain<N_CHANNELS>;

    #[test]
    fn new() {
        let gain = GainT::new();

        let cutoff: f32;
        let gain_default = 1.0;
        let tau_default = 0.05;

        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(gain.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(gain.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(gain.coeffs.gain, gain_default);
    }

    #[test]
    fn set_sample_rate() {
        let mut gain = GainT::new();

        gain.set_sample_rate(SAMPLE_RATE);
        assert_eq!(gain.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE)
    }

    #[test]
    fn reset() {
        let gain_par = 1.1;
        let sticky_tresh = 0.1;
        let mut gain = GainT::new();

        gain.set_gain_lin(gain_par);
        gain.set_sticky_thresh(sticky_tresh);
        gain.set_sample_rate(SAMPLE_RATE);
        gain.reset();

        assert!(gain.coeffs.gain == gain_par);
        assert_eq!(gain.coeffs.smooth_coeffs.sticky_thresh, sticky_tresh);
    }

    #[test]
    fn process() {
        const N_SAMPLES: usize = 2;
        // Wrapper
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);
        gain.set_sticky_thresh(0.0);

        let sample_0: [f32; N_CHANNELS] = [1.0, 1.0];
        let sample_1: [f32; N_CHANNELS] = [0.0, 0.0];
        let x: [&[f32]; 2] = [&sample_0, &sample_1];

        let mut out_0: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut out_1: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut y: [&mut [f32]; 2] = [&mut out_0, &mut out_1];

        // C
        let mut coeffs = bw_gain_coeffs::default();

        let mut out_0_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut out_1_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut y_c: [&mut [f32]; 2] = [&mut out_0_c, &mut out_1_c];

        // Process wrapper
        gain.reset();
        gain.process(&x, &mut y, N_SAMPLES);

        // Process C
        unsafe {
            bw_gain_init(&mut coeffs);
            bw_gain_set_sample_rate(&mut coeffs, SAMPLE_RATE);
            bw_gain_set_sticky_thresh(&mut coeffs, 0.0);
            bw_gain_reset_coeffs(&mut coeffs);
            bw_gain_update_coeffs_ctrl(&mut coeffs);
            (0..N_SAMPLES).for_each(|sample| {
                bw_gain_update_coeffs_audio(&mut coeffs);
                (0..N_CHANNELS).for_each(|channel| {
                    y_c[channel][sample] = bw_gain_process1(&coeffs, x[channel][sample])
                });
            });
        }

        (0..N_SAMPLES).for_each(|sample| {
            (0..N_CHANNELS).for_each(|channel| {
                assert_eq!(
                    y[channel][sample], y_c[channel][sample],
                    "sample {sample} and channel {channel} does not match"
                );
            });
        });
    }

    #[test]
    fn set_gain_lin() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);
        let gain_lin = 0.9;
        gain.set_gain_lin(gain_lin);

        assert_eq!(gain.coeffs.gain, gain_lin);
    }

    #[test]
    fn set_gain_db() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);
        let gain_db = 6.0;
        let gain_val;
        unsafe {
            gain_val = bw_dB2linf(gain_db);
        }

        gain.set_gain_db(gain_db);

        assert_eq!(gain.coeffs.gain, gain_val);
    }

    #[test]
    fn set_smooth_tau() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);

        let tau = 0.01;
        let cutoff;
        gain.set_smooth_tau(tau);

        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau);
        }

        assert_eq!(gain.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(gain.coeffs.smooth_coeffs.cutoff_down, cutoff);
    }

    #[test]
    fn set_sticky_thresh() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);

        let sticky_tresh = 0.12;
        gain.set_sticky_thresh(sticky_tresh);

        assert_eq!(gain.coeffs.smooth_coeffs.sticky_thresh, sticky_tresh);
    }

    #[test]
    fn set_sticky_mode() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);

        gain.set_sticky_mode(StickyMode::Abs);

        assert_eq!(
            gain.coeffs.smooth_coeffs.sticky_mode,
            bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs
        );
    }

    #[test]
    fn get_gain_lin() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);
        let gain_db = 6.0;
        let gain_val;
        unsafe {
            gain_val = bw_dB2linf(gain_db);
        }

        gain.set_gain_db(gain_db);

        assert_eq!(gain.get_gain_lin(), gain_val);
    }

    #[test]
    fn get_gain_cur() {
        let mut gain = GainT::new();
        gain.set_sample_rate(SAMPLE_RATE);
        let gain_lin = 0.9;
        let n_samples = 2;

        let x_sample1 = [1.0, 1.0];
        let x_sample2 = [0.0, 0.0];
        let x: [&[f32]; 2] = [&x_sample1, &x_sample2];

        let mut y_sample1 = [0.0, 0.0];
        let mut y_sample2 = [0.0, 0.0];
        let mut y: [&mut [f32]; 2] = [&mut y_sample1, &mut y_sample2];

        gain.set_gain_lin(gain_lin);
        gain.reset();
        gain.process(&x, &mut y, n_samples);
        assert_eq!(gain.get_gain_cur(), gain.coeffs.smooth_state.y_z1);
    }
}
