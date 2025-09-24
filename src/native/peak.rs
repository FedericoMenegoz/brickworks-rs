//! Second-order peak filter with unitary gain at DC and asymptotically as
//! frequency increases.
//!
//! The quality factor of the underlying bandpass filter can be either directly
//! controlled via the Q parameter or indirectly through the bandwidth parameter,
//! which designates the distance in octaves between midpoint gain frequencies,
//! i.e., frequencies with gain = peak gain / 2 in dB terms.
//! The use_bandiwdth parameter allows you to choose which parameterization to use.
//!
//! # Example
//! ```rust
//! use brickworks_rs::native::peak::*;
//!
//! const N_CHANNELS: usize = 2;
//! const N_SAMPLES: usize = 8;
//! const SAMPLE_RATE: f32 = 44_100.0;
//! const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
//!     &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!     &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//! ];
//! fn main() {
//!     let mut peak = Peak::new();
//!     peak.set_sample_rate(SAMPLE_RATE);
//!     peak.set_q(1.4);
//!
//!     let x0 = [0.0, 0.0];
//!
//!     peak.reset_multi(&x0, None);
//!
//!     let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
//!     peak.process(&PULSE_INPUT, &mut y, N_SAMPLES);
//! }
//! ```
//!
//! # Notes
//! This module provides a native Rust implementation, but the same interface is
//! also available via bindings to the original C library at [crate::c_wrapper::peak].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_peak).
//!
use bitflags::bitflags;

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
use crate::native::{
    math::{db2linf, pow2f, rcpf, sqrtf},
    mm2::{MM2Coeffs, MM2State},
};
/// Second-order peak filter with unitary gain at DC and asymptotically 
/// as frequency increases.
/// 
/// This struct manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`PeakCoeffs`] 
/// - [`PeakState`] 
/// 
/// # Usage
/// ```rust 
/// use brickworks_rs::native::peak::Peak;
/// const N_CHANNELS: usize = 2;
/// let mut peak = Peak::<N_CHANNELS>::new();
/// peak.set_sample_rate(48_000.0);
/// peak.set_cutoff(1_000.0);
/// peak.set_q(0.707);
/// peak.set_prewarp_at_cutoff(true);
/// peak.reset(None, None);
/// // process audio with svf.process(...)
/// ```
pub struct Peak<const N_CHANNELS: usize> {
    coeffs: PeakCoeffs<N_CHANNELS>,
    states: [PeakState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Peak<N_CHANNELS> {
    /// Creates a new instance with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            coeffs: PeakCoeffs::new(),
            states: [PeakState::default(); N_CHANNELS],
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
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self
                    .coeffs
                    .reset_state(&mut self.states[channel], x0.unwrap_or(0.0));
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.coeffs
                    .reset_state(&mut self.states[channel], x0.unwrap_or(0.0));
            }),
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
    /// Sets the input cutoff frequency value (Hz).
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
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.coeffs.set_prewarp_at_cutoff(value);
    }
    /// Sets the prewarping frequency value (Hz).
    ///
    /// Only used when the prewarp_at_cutoff parameter is off and however internally
    /// limited to avoid instability.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        self.coeffs.set_prewarp_freq(value);
    }
    /// Sets the peak gain parameter to the given value (linear gain) in coeffs.
    ///
    /// Valid range: [1e-30, 1e30].
    ///
    /// If actually using the bandwidth parameter to control q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_peak_gain_lin(&mut self, value: f32) {
        self.coeffs.set_peak_gain_lin(value);
    }
    /// Sets the peak gain parameter to the given value (dB) in coeffs.
    ///
    /// Valid range: [-600.0, 600.0].
    ///
    /// If actually using the bandwidth parameter to control q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_peak_gain_db(&mut self, value: f32) {
        self.coeffs.set_peak_gain_db(value);
    }
    /// Sets the bandwidth value (octaves) in coeffs.
    ///
    /// Valid range: [1e-6, 90.0].
    ///
    /// If actually using the bandwidth parameter to control q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 2.543_106_6.
    #[inline(always)]
    pub fn set_bandwidth(&mut self, value: f32) {
        self.coeffs.set_bandwidth(value);
    }
    /// Sets whether the quality factor should be controlled via the bandwidth
    /// parameter (true) or via the Q parameter (false).
    ///
    /// Default value: true (use bandwidth parameter).
    #[inline(always)]
    pub fn set_use_bandwidth(&mut self, value: bool) {
        self.coeffs.set_use_bandwidth(value);
    }
}

impl<const N_CHANNELS: usize> Default for Peak<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
pub struct PeakCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    mm2_coeffs: MM2Coeffs<N_CHANNELS>,

    // Coefficients
    bw_k: f32,

    // Parameters
    q: f32,
    peak_gain: f32,
    bandwidth: f32,
    use_bandwidth: bool,
    param_changed: ParamChanged,
}

impl<const N_CHANNELS: usize> PeakCoeffs<N_CHANNELS> {
    /// Creates a new instance with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            mm2_coeffs: MM2Coeffs::<N_CHANNELS>::new(),
            bw_k: 0.0,
            q: 0.5,
            peak_gain: 1.0,
            bandwidth: 2.543_106_6, /*06327224*/
            use_bandwidth: true,
            param_changed: ParamChanged::all(),
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

        self.mm2_coeffs.set_sample_rate(sample_rate);
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.param_changed = ParamChanged::all();
        self.update_mm2_params();
        self.mm2_coeffs.reset_coeffs();
    }
    /// Resets the given state to its initial values using the initial input value `x0`.
    ///
    /// # Returns
    /// The corresponding initial output value.
    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut PeakState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        let y = self.mm2_coeffs.reset_state(&mut state.mm2_state, x0);

        debug_assert!(y.is_finite());

        y
    }
    /// Resets each of the `N_CHANNELS` states to its initial values using the
    /// corresponding input values in `x0`.
    ///
    /// The output values are written into `y0` if it is not `None`.
    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [PeakState],
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
    /// Triggers control-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        self.update_mm2_params();
        self.mm2_coeffs.update_coeffs_ctrl();
    }
    /// Triggers audio-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.mm2_coeffs.update_coeffs_audio();
    }
    /// Processes one input sample `x`, while using and updating the given state.
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1(&mut self, state: &mut PeakState, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let y = self.mm2_coeffs.process1(&mut state.mm2_state, x);

        debug_assert!(y.is_finite());

        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and fills the first
    /// `n_samples` of the output buffer `y`, while using and updating both coeffs
    /// and state (control and audio rate).
    #[inline(always)]
    pub fn process(&mut self, state: &mut PeakState, x: &[f32], y: &mut [f32], n_samples: usize) {
        self.update_coeffs_ctrl();
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample]);
        });
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and fills the
    /// first `n_samples` of the `N_CHANNELS` output buffers `y`, while using and updating
    /// both the common coeffs and each of the N_CHANNELS states (control and audio rate).
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [PeakState],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.update_coeffs_ctrl();
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            (0..N_CHANNELS).for_each(|channel| {
                y[channel][sample] = self.process1(&mut states[channel], x[channel][sample])
            });
        });
    }
    /// Sets the input cutoff frequency value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e12, value);
        }

        self.mm2_coeffs.set_cutoff(value);
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
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e6, value);
        }

        if self.q != value {
            self.q = value;
            self.param_changed |= ParamChanged::Q
        }
    }
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.mm2_coeffs.set_prewarp_at_cutoff(value);
    }
    /// Sets the prewarping frequency value (Hz).
    ///
    /// Only used when the prewarp_at_cutoff parameter is off and however internally
    /// limited to avoid instability.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e12, value);
        }
        self.mm2_coeffs.set_prewarp_freq(value);
    }
    /// Sets the peak gain parameter to the given value (linear gain) in coeffs.
    ///
    /// Valid range: [1e-30, 1e30].
    ///
    /// If actually using the bandwidth parameter to control Q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_peak_gain_lin(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-30..=1e30, value);
        }

        if self.peak_gain != value {
            self.peak_gain = value;
            self.param_changed |= ParamChanged::PEAK_GAIN;
        }
    }
    /// Sets the peak gain parameter to the given value (dB) in coeffs.
    ///
    /// Valid range: [-600.0, 600.0].
    ///
    /// If actually using the bandwidth parameter to control q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_peak_gain_db(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(-600.0..=600.0, value);
        }
        self.set_peak_gain_lin(db2linf(value));
    }
    /// Sets the bandwidth value (octaves) in coeffs.
    ///
    /// Valid range: [1e-6, 90.0].
    ///
    /// If actually using the bandwidth parameter to control q, by the time
    /// reset_*(), update_coeffs_*(), or process*() is called,
    /// sqrtf(pow2f(bandwidth) * peak_gain) * rcpf(pow2f(bandwidth) - 1.0)
    /// must be in [1e-6, 1e6].
    ///
    /// Default value: 2.543_106_6.
    #[inline(always)]
    pub fn set_bandwidth(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=90.0, value);
        }

        if self.bandwidth != value {
            self.bandwidth = value;
            self.param_changed |= ParamChanged::BANDWIDTH;
        }
    }
    /// Sets whether the quality factor should be controlled via the bandwidth
    /// parameter (true) or via the Q parameter (false).
    ///
    /// Default value: true (use bandwidth parameter).
    #[inline(always)]
    pub fn set_use_bandwidth(&mut self, value: bool) {
        if self.use_bandwidth != value {
            self.use_bandwidth = value;
            self.param_changed |= ParamChanged::Q | ParamChanged::BANDWIDTH;
        }
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
    pub fn coeffs_is_valid(&mut self) -> bool {
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
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }

    // Private
    #[inline(always)]
    fn update_mm2_params(&mut self) {
        if !self.param_changed.is_empty() {
            if self.use_bandwidth {
                if self
                    .param_changed
                    .intersects(ParamChanged::PEAK_GAIN | ParamChanged::BANDWIDTH)
                {
                    if self.param_changed.contains(ParamChanged::BANDWIDTH) {
                        self.bw_k = pow2f(self.bandwidth);
                    }
                    let q = sqrtf(self.bw_k * self.peak_gain) * rcpf(self.bw_k - 1.0);
                    self.mm2_coeffs.set_q(q);
                    self.mm2_coeffs
                        .set_coeff_bp((self.peak_gain - 1.0) * rcpf(q));
                }
            } else if self
                .param_changed
                .intersects(ParamChanged::PEAK_GAIN | ParamChanged::Q)
            {
                if self.param_changed.contains(ParamChanged::Q) {
                    self.mm2_coeffs.set_q(self.q);
                }
                self.mm2_coeffs
                    .set_coeff_bp((self.peak_gain - 1.0) * rcpf(self.q));
            }

            self.param_changed = ParamChanged::empty();
        }
    }
}

impl<const N_CHANNELS: usize> Default for PeakCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    struct ParamChanged: u32 {
        const Q = 1;
        const PEAK_GAIN = 1<<1;
        const BANDWIDTH = 1<<2;

        // added this as peak_initialization test fail caused by param_changed
        const _ = !0;
    }
}

/// Internal state and related.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct PeakState {
    //Sub-components
    mm2_state: MM2State,
}

#[cfg(test)]
pub(crate) mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_peak_coeffs, bw_peak_state, peak::Peak as PeakWrapper},
        native::mm2::tests::{assert_mm2_coeffs, assert_mm2_state},
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type PeakT = Peak<N_CHANNELS>;
    type PeakWrapperT = PeakWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_peak = PeakT::new();
        let c_peak = PeakWrapperT::new();

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn reset_none() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();

        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.reset(None, None);
        c_peak.reset(None, None);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn reset_some() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();

        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.3;
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_peak.reset(Some(x0), Some(&mut rust_y0));
        c_peak.reset(Some(x0), Some(&mut c_y0));

        assert_peak(&rust_peak, &c_peak);
        assert_eq!(rust_y0, c_y0);
    }

    #[test]
    fn reset_multi() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.045, 0.045];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_peak.reset_multi(&x0, Some(&mut rust_y0));
        c_peak.reset_multi(&x0, Some(&mut c_y0));

        assert_peak(&rust_peak, &c_peak);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let bw = 39.0;
        let cutoff = 456.1;
        let prewarp_freq = 460.1;
        let gain_db = 56.0;

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];

        rust_peak.set_bandwidth(bw);
        rust_peak.set_cutoff(cutoff);
        rust_peak.set_peak_gain_db(gain_db);
        rust_peak.set_prewarp_at_cutoff(true);
        rust_peak.set_prewarp_freq(prewarp_freq);
        rust_peak.set_use_bandwidth(true);

        c_peak.set_bandwidth(bw);
        c_peak.set_cutoff(cutoff);
        c_peak.set_peak_gain_db(gain_db);
        c_peak.set_prewarp_at_cutoff(true);
        c_peak.set_prewarp_freq(prewarp_freq);
        c_peak.set_use_bandwidth(true);

        rust_peak.reset(None, None);
        c_peak.reset(None, None);
        rust_peak.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_peak.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_peak(&rust_peak, &c_peak);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    fn set_cutoff() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let cutoff = 261.63;

        rust_peak.set_cutoff(cutoff);
        c_peak.set_cutoff(cutoff);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got -1e0")]
    #[test]
    fn set_cutoff_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_cutoff(-1.0);
    }

    #[test]
    fn set_q() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let q = 12.0;

        rust_peak.set_q(q);
        c_peak.set_q(q);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e6], got 1e7")]
    #[test]
    fn set_q_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_q(1e7);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_prewarp_at_cutoff(true);
        c_peak.set_prewarp_at_cutoff(true);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 340.12;
        rust_peak.set_prewarp_freq(prewarp_freq);
        c_peak.set_prewarp_freq(prewarp_freq);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 1e13")]
    #[test]
    fn set_prewarp_freq_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_prewarp_freq(1e13);
    }

    #[test]
    fn set_peak_gain_lin() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let gain_lin = 1.76;

        rust_peak.set_peak_gain_lin(gain_lin);
        c_peak.set_peak_gain_lin(gain_lin);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-30, 1e30], got 0e0")]
    #[test]
    fn set_peak_gain_lin_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_peak_gain_lin(0.0);
    }

    #[test]
    fn set_peak_gain_db() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let gain_db = 450.0;

        rust_peak.set_peak_gain_lin(gain_db);
        c_peak.set_peak_gain_lin(gain_db);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [-6e2, 6e2], got 6.001e2")]
    #[test]
    fn set_peak_gain_db_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_peak_gain_db(600.1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 9e1], got 9.01e1")]
    #[test]
    fn set_bandwidth_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_bandwidth(90.1);
    }

    #[test]
    fn set_bandwidth() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        let bw = 34.12;

        rust_peak.set_peak_gain_lin(bw);
        c_peak.set_peak_gain_lin(bw);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_use_bandwidth() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_use_bandwidth(true);
        c_peak.set_use_bandwidth(true);

        assert_peak(&rust_peak, &c_peak);
    }

    fn assert_peak<const N_CHANNELS: usize>(
        rust_peak: &Peak<N_CHANNELS>,
        c_peak: &PeakWrapper<N_CHANNELS>,
    ) {
        assert_peak_coeffs(&rust_peak.coeffs, &c_peak.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_peak_state(&rust_peak.states[channel], &c_peak.states[channel]);
        });
    }

    pub(crate) fn assert_peak_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &PeakCoeffs<N_CHANNELS>,
        c_coeffs: &bw_peak_coeffs,
    ) {
        // Sub-components
        assert_mm2_coeffs(&rust_coeffs.mm2_coeffs, &c_coeffs.mm2_coeffs);

        assert_eq!(rust_coeffs.bw_k, c_coeffs.bw_k);
        assert_eq!(rust_coeffs.q, c_coeffs.Q);
        assert_eq!(rust_coeffs.peak_gain, c_coeffs.peak_gain);
        assert_eq!(rust_coeffs.bandwidth, c_coeffs.bandwidth);
        assert_eq!(rust_coeffs.use_bandwidth, c_coeffs.use_bandwidth != 0,);
        assert_eq!(
            rust_coeffs.param_changed.bits(),
            c_coeffs.param_changed as u32
        );
    }

    pub(crate) fn assert_peak_state(rust_state: &PeakState, c_state: &bw_peak_state) {
        assert_mm2_state(&rust_state.mm2_state, &c_state.mm2_state);
    }
}
