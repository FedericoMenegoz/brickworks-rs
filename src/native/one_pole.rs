//! **One-pole (6 dB/oct) lowpass filter** with unitary DC gain, separate attack and decay time constants, and sticky target-reach threshold.
//!
//! This is better suited to implement smoothing than [crate::native::lp1].
//!
//! # Example
//! ```rust
//! use brickworks_rs::native::one_pole::*;
//!
//! const CUTOFF: f32 = 100_000.0;
//! const STICKY_THRESH: f32 = 0.9;
//! const N_CHANNELS: usize = 2;
//! const N_SAMPLES: usize = 1;
//! const SAMPLE_RATE: f32 = 48_000.0;
//!
//! // Create a new OnePole filter instance for N_CHANNELS
//! let mut one_pole = OnePole::<N_CHANNELS>::new();
//!
//! // Input signal: one sample per channel
//! let x:[&[f32]; N_CHANNELS] = [&[1.0], &[0.0]];
//!
//! // Output buffer, same shape as input
//! let mut y_ch1 = [0.0];
//! let mut y_ch2 = [0.0];
//! let mut y: [Option<&mut[f32]>; N_CHANNELS] = [Some(&mut y_ch1), Some(&mut y_ch2)];
//!
//! // Configure the filter
//! one_pole.set_sample_rate(SAMPLE_RATE);
//! one_pole.set_cutoff(CUTOFF);
//! one_pole.set_sticky_mode(StickyMode::Rel);
//! one_pole.set_sticky_thresh(STICKY_THRESH);
//!
//! // Initialize the filter state for each channel
//! one_pole.reset(None, None);
//!
//! // Process one sample per channel
//! one_pole.process(&x, Some(&mut y), N_SAMPLES);
//!
//! // Output the filtered result
//! println!("Filtered output: {:?}", y);
//! ```
//!
//! # Notes
//! This module provides a native Rust implementation of the filter, but the same interface is
//! also available via bindings to the original C library at [crate::c_wrapper::one_pole].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_one_pole).
use super::math::rcpf;
#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
use crate::native::math::{INVERSE_2_PI, NANO};
use bitflags::bitflags;
/// One-pole (6 dB/oct) lowpass filter with unitary DC gain, separate attack and decay time constants, and sticky target-reach threshold.
/// This struct manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`OnePoleCoeffs`]
/// - [`OnePoleState`]
///
/// # Usage
/// ```rust
/// use brickworks_rs::native::one_pole::{OnePole, StickyMode};
/// const N_CHANNELS: usize = 2;
/// let mut op = OnePole::<N_CHANNELS>::new();
/// op.set_sample_rate(48_000.0);
/// op.set_cutoff(1_000.0);
/// op.set_sticky_mode(StickyMode::Rel);
/// op.set_sticky_thresh(0.1);
/// // process audio with op.process(...)
/// ```
#[derive(Debug)]
pub struct OnePole<const N_CHANNELS: usize> {
    coeffs: OnePoleCoeffs<N_CHANNELS>,
    states: [OnePoleState; N_CHANNELS],
}
/// Distance metrics for sticky behavior.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum StickyMode {
    /// `StickyMode::Abs`: absolute difference ( `|out - in|` );
    Abs,
    /// `StickyMode::Rel`: relative difference with respect to input (`|out - in| / |in|`).
    Rel,
}
impl<const N_CHANNELS: usize> OnePole<N_CHANNELS> {
    /// Creates a new instance with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        OnePole {
            coeffs: OnePoleCoeffs::new(),
            states: [OnePoleState { y_z1: 0.0 }; N_CHANNELS],
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }
    /// Resets the states and coeffs for all channels to the initial input `x0`,
    /// or to 0 if `x0` is not provided.
    /// If `y0` is provided, the resulting initial outputs are stored in it.
    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(out) => {
                (0..N_CHANNELS).for_each(|channel| {
                    out[channel] = self.states[channel].reset(x0.unwrap_or(0.));
                });
            }
            None => {
                (0..N_CHANNELS).for_each(|channel| {
                    self.states[channel].reset(x0.unwrap_or(0.));
                });
            }
        }
    }
    /// Resets the state and coefficients for all channels using the provided initial
    /// input values.
    ///
    /// Both the coefficients and all channel states are reset.
    /// If `y0` is `Some`, the initial outputs are written into it.
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        self.coeffs.reset_state_multi(&mut self.states, x0, y0);
    }
    /// Processes a block of input samples using a one-pole per channel.
    ///
    /// This method applies a one-pole low-pass filter to `n_samples` of input per channel.
    /// It supports multiple filtering modes based on whether the filter is *asymmetric*
    /// (`cutoff_up != cutoff_down`), and whether it uses *sticky*
    /// behavior (which prevents updates when the change is below a threshold).
    ///
    /// ### Parameters
    /// - `x`: A slice of input sample vectors, one per channel.
    /// - `y`: An optional mutable slice of output buffers, one per channel.
    /// - `n_samples`: Number of samples to process per channel.
    ///
    /// ### Filter Behavior
    /// The internal logic chooses the processing mode based on the current coefficients:
    ///
    /// - `process1`: Symmetric filtering with no sticky behavior.
    /// - `process1_sticky_abs`: Symmetric filtering with sticky behavior using an absolute difference metric.
    /// - `process1_sticky_rel`: Symmetric filtering with sticky behavior using a relative difference metric.
    /// - `process1_asym`: Asymmetric filtering with no sticky behavior.
    /// - `process1_asym_sticky_abs`: Asymmetric filtering with sticky behavior using an absolute difference metric.
    /// - `process1_asym_sticky_rel`: Asymmetric filtering with sticky behavior using a relative difference metric.
    ///
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        n_samples: usize,
    ) {
        self.coeffs.process_multi(&mut self.states, x, y, n_samples);
    }
    /// Sets both the upgoing (attack) and downgoing (decay) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling both `set_cutoff_up()` and `set_cutoff_down()` with same
    /// value or calling `set_tau()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        self.coeffs.set_cutoff(value);
    }
    /// Sets the upgoing (attack) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_tau_up()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_cutoff_up(&mut self, value: f32) {
        self.coeffs.set_cutoff_up(value);
    }
    /// Sets the downgoing (decay) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_tau_down()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_cutoff_down(&mut self, value: f32) {
        self.coeffs.set_cutoff_down(value);
    }
    /// Sets both the upgoing (attack) and downgoing (decay) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling both `set_tau_up()` and `set_tau_down()` with same
    /// value or calling set_cutoff() with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_tau(&mut self, value: f32) {
        self.coeffs.set_tau(value);
    }
    /// Sets the upgoing (attack) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_cutoff_up()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_tau_up(&mut self, value: f32) {
        self.coeffs.set_tau_up(value);
    }
    /// Sets both the downgoing (decay) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_cutoff_down()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    #[inline(always)]
    pub fn set_tau_down(&mut self, value: f32) {
        self.coeffs.set_tau_down(value);
    }
    /// Sets the target-reach threshold specified by value in `OnePole<N_CHANNELS>`.
    ///
    /// When the difference between the output and the input would fall under such threshold according
    /// to the current distance metric (see [StickyMode]), the output is forcefully set to be equal to the input value.
    ///
    ///
    /// ### Parameters
    /// - `value`: sticky threshhold, default is `0.0` and its valid range is [0.0, 1e18].
    /// - `value`: sticky threshhold, default is `0.0` and its valid range is [0.0, 1e18].
    ///
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        self.coeffs.set_sticky_thresh(value);
    }
    /// Sets the current distance metric for sticky behavior to value in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: sticky mode, default is [StickyMode::Abs]
    /// - `value`: sticky mode, default is [StickyMode::Abs]
    ///
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: StickyMode) {
        self.coeffs.set_sticky_mode(value);
    }
    /// Returns the current target-reach threshold in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_sticky_thresh(&self) -> f32 {
        self.coeffs.get_sticky_thresh()
    }
    /// Returns the current distance metric for sticky behavior in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_sticky_mode(&self) -> StickyMode {
        self.coeffs.get_sticky_mode()
    }
    /// Returns the last output sample as stored in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_y_z1(&self, channel: usize) -> f32 {
        self.states[channel].y_z1
    }
}

impl<const N_CHANNELS: usize> Default for OnePole<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
#[derive(Clone, Debug, Copy)]
pub struct OnePoleCoeffs<const N_CHANNELS: usize> {
    fs_2pi: f32,
    m_a1u: f32,
    m_a1d: f32,
    st2: f32,
    cutoff_up: f32,
    cutoff_down: f32,
    sticky_thresh: f32,
    sticky_mode: StickyMode,
    param_changed: ParamChanged,
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    struct ParamChanged: u32 {
        const CUTOFF_UP = 1;
        const CUTOFF_DOWN = 1<<1;
        const STICKY_TRESH = 1<<2;

        // added this as one_pole_initialization test fail caused by param_changed
        const _ = !0;
    }
}

/// Internal state and related.
#[derive(Debug, Clone, Copy)]
pub struct OnePoleState {
    y_z1: f32,
}
impl OnePoleState {
    #[inline(always)]
    pub(crate) fn new() -> Self {
        Self { y_z1: 0.0 }
    }

    #[inline(always)]
    pub(crate) fn reset(&mut self, x0: f32) -> f32 {
        self.y_z1 = x0;
        x0
    }

    #[inline(always)]
    pub(crate) fn get_y_z1(&self) -> f32 {
        self.y_z1
    }
}

impl Default for OnePoleState {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N_CHANNELS: usize> OnePoleCoeffs<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            fs_2pi: Default::default(),
            m_a1u: Default::default(),
            m_a1d: Default::default(),
            st2: Default::default(),
            cutoff_up: f32::INFINITY,
            cutoff_down: f32::INFINITY,
            sticky_thresh: 0.0,
            sticky_mode: StickyMode::Abs,
            param_changed: ParamChanged::all(),
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
        self.fs_2pi = INVERSE_2_PI * sample_rate;
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.param_changed = ParamChanged::all();
        self.do_update_coeffs_ctrl();
    }
    /// Resets the given state to its initial values using the initial input value `x0`.
    ///
    /// Returns the corresponding initial output value.
    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut OnePoleState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());
        state.reset(x0)
    }
    /// Resets each of the `N_CHANNELS` states to its initial values using the
    /// corresponding input values in `x0`.
    ///
    /// The output values are written into `y0` if it is not `None`.
    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [OnePoleState],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        // No need to check states are in different addresses cause it is
        // enforced by design
        if let Some(out) = y0 {
            (0..N_CHANNELS).for_each(|channel| {
                out[channel] = states[channel].reset(x0[channel]);
            });
        } else {
            (0..N_CHANNELS).for_each(|channel| {
                states[channel].reset(x0[channel]);
            });
        }
    }
    /// Triggers control-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        self.do_update_coeffs_ctrl();
    }

    // Not implemented yet: C version only contained assertions
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Triggers audio-rate update of coefficients.
    // #[inline(always)]
    // fn update_coeffs_audio(&self) {
    //     todo!()
    // }

    /// Processes a single input sample `x`, updating the provided `state`.
    /// Assumes that the upgoing and downgoing cutoff/tau are equal, and the
    /// target-reach threshold is `0.0`.
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let y = x + self.m_a1u * (state.get_y_z1() - x);
        state.y_z1 = y;
        y
    }
    /// Processes a single input sample `x` with sticky absolute threshold behavior,
    /// updating the provided `state`.
    /// Assumes upgoing and downgoing cutoff/tau are equal, the target-reach
    /// threshold is not `0.0`, and sticky mode is absolute (`StickyMode::Abs`).
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_sticky_abs(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let mut y = x + self.m_a1u * (state.get_y_z1() - x);
        let d = y - x;
        if d * d <= self.st2 {
            y = x
        }
        state.y_z1 = y;

        y
    }
    /// Processes a single input sample `x` with sticky relative threshold behavior,
    /// updating the provided `state`.
    /// Assumes upgoing and downgoing cutoff/tau are equal, the target-reach
    /// threshold is not `0.0`, and sticky mode is relative (`StickyMode::Rel`).
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_sticky_rel(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let mut y = x + self.m_a1u * (state.get_y_z1() - x);
        let d = y - x;
        if d * d <= self.st2 * x * x {
            y = x
        }
        state.y_z1 = y;

        y
    }
    /// Processes a single input sample `x`, updating the provided `state`.
    /// Assumes upgoing and downgoing cutoff/tau may differ and the target-reach threshold is `0.0`.
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_asym(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let y_z1 = state.get_y_z1();
        let ma1 = if x >= y_z1 { self.m_a1u } else { self.m_a1d };

        let y = x + ma1 * (y_z1 - x);
        state.y_z1 = y;

        y
    }
    /// Processes a single input sample `x` with sticky absolute threshold,
    /// updating the provided `state`.
    /// Assumes upgoing and downgoing cutoff/tau may differ, the target-reach
    /// threshold is not `0.0`, and sticky mode is absolute (`StickyMode::Abs`).
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_asym_sticky_abs(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let y_z1 = state.get_y_z1();
        let ma1 = if x >= y_z1 { self.m_a1u } else { self.m_a1d };

        let mut y = x + ma1 * (y_z1 - x);
        let d = y - x;
        if d * d <= self.st2 {
            y = x;
        }
        state.y_z1 = y;

        y
    }
    /// Processes a single input sample `x` with sticky relative threshold,
    /// updating the provided `state`.
    /// Assumes upgoing and downgoing cutoff/tau may differ, the target-reach
    /// threshold is not `0.0`, and sticky mode is relative (`StickyMode::Rel`).
    ///
    /// Returns the corresponding output sample.
    #[inline(always)]
    pub fn process1_asym_sticky_rel(&mut self, state: &mut OnePoleState, x: f32) -> f32 {
        let y_z1 = state.get_y_z1();
        let ma1 = if x >= y_z1 { self.m_a1u } else { self.m_a1d };

        let mut y = x + ma1 * (y_z1 - x);
        let d = y - x;
        if d * d <= self.st2 * x * x {
            y = x;
        }
        state.y_z1 = y;

        y
    }
    /// Processes the first `n_samples` of the input buffer `x` and writes the results
    /// to the first `n_samples` of the output buffer `y`, while updating both
    /// coefficients and state (control and audio rate).
    ///
    #[inline(always)]
    pub fn process(
        &mut self,
        state: &mut OnePoleState,
        x: &[f32],
        y: Option<&mut [f32]>,
        n_samples: usize,
    ) {
        self.update_coeffs_ctrl();
        match y {
            Some(y_values) => {
                if self.is_asym() {
                    if self.is_sticky() {
                        match self.get_sticky_mode() {
                            StickyMode::Abs => {
                                (0..n_samples).for_each(|sample| {
                                    y_values[sample] =
                                        self.process1_asym_sticky_abs(state, x[sample]);
                                });
                            }
                            StickyMode::Rel => {
                                (0..n_samples).for_each(|sample| {
                                    y_values[sample] =
                                        self.process1_asym_sticky_rel(state, x[sample]);
                                });
                            }
                        }
                    } else {
                        (0..n_samples).for_each(|sample| {
                            y_values[sample] = self.process1_asym(state, x[sample]);
                        });
                    }
                } else if self.is_sticky() {
                    match self.get_sticky_mode() {
                        StickyMode::Abs => {
                            (0..n_samples).for_each(|sample| {
                                y_values[sample] = self.process1_sticky_abs(state, x[sample]);
                            });
                        }
                        StickyMode::Rel => {
                            (0..n_samples).for_each(|sample| {
                                y_values[sample] = self.process1_sticky_rel(state, x[sample]);
                            });
                        }
                    }
                } else {
                    (0..n_samples).for_each(|sample| {
                        y_values[sample] = self.process1(state, x[sample]);
                    });
                }
            }
            None => {
                if self.is_asym() {
                    if self.is_sticky() {
                        match self.get_sticky_mode() {
                            StickyMode::Abs => {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_asym_sticky_abs(state, x[sample]);
                                });
                            }
                            StickyMode::Rel => {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_asym_sticky_rel(state, x[sample]);
                                });
                            }
                        }
                    } else {
                        (0..n_samples).for_each(|sample| {
                            self.process1_asym(state, x[sample]);
                        });
                    }
                } else if self.is_sticky() {
                    match self.get_sticky_mode() {
                        StickyMode::Abs => {
                            (0..n_samples).for_each(|sample| {
                                self.process1_sticky_abs(state, x[sample]);
                            });
                        }
                        StickyMode::Rel => {
                            (0..n_samples).for_each(|sample| {
                                self.process1_sticky_rel(state, x[sample]);
                            });
                        }
                    }
                } else {
                    (0..n_samples).for_each(|sample| {
                        self.process1(state, x[sample]);
                    });
                }
            }
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and writes the
    /// results to the first `n_samples` of the `N_CHANNELS` output buffers `y`,
    /// while updating the shared coefficients and each channel's state (control and audio rate).
    ///
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [OnePoleState; N_CHANNELS],
        x: &[&[f32]; N_CHANNELS],
        y: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        n_samples: usize,
    ) {
        // As for reset_multi no need to check states are in
        // different addresses cause it is enforced by design
        // Still need to investigate if the other assertions
        // are sanity check only needed in c

        self.update_coeffs_ctrl();
        if let Some(y_values) = y {
            if self.is_asym() {
                if self.is_sticky() {
                    match self.sticky_mode {
                        StickyMode::Abs => {
                            (0..N_CHANNELS).for_each(|channel| {
                                if let Some(y_value) = y_values[channel].as_deref_mut() {
                                    (0..n_samples).for_each(|sample| {
                                        y_value[sample] = self.process1_asym_sticky_abs(
                                            &mut states[channel],
                                            x[channel][sample],
                                        );
                                    });
                                } else {
                                    (0..n_samples).for_each(|sample| {
                                        self.process1_asym_sticky_abs(
                                            &mut states[channel],
                                            x[channel][sample],
                                        );
                                    });
                                }
                            });
                        }
                        StickyMode::Rel => {
                            (0..N_CHANNELS).for_each(|channel| {
                                if let Some(y_value) = y_values[channel].as_deref_mut() {
                                    (0..n_samples).for_each(|sample| {
                                        y_value[sample] = self.process1_asym_sticky_rel(
                                            &mut states[channel],
                                            x[channel][sample],
                                        );
                                    });
                                } else {
                                    (0..n_samples).for_each(|sample| {
                                        self.process1_asym_sticky_rel(
                                            &mut states[channel],
                                            x[channel][sample],
                                        );
                                    });
                                }
                            });
                        }
                    }
                } else {
                    (0..N_CHANNELS).for_each(|channel| {
                        if let Some(y_value) = y_values[channel].as_deref_mut() {
                            (0..n_samples).for_each(|sample| {
                                y_value[sample] =
                                    self.process1_asym(&mut states[channel], x[channel][sample]);
                            });
                        } else {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym(&mut states[channel], x[channel][sample]);
                            });
                        }
                    });
                }
            } else if self.is_sticky() {
                match self.sticky_mode {
                    StickyMode::Abs => {
                        (0..N_CHANNELS).for_each(|channel| {
                            if let Some(y_value) = y_values[channel].as_deref_mut() {
                                (0..n_samples).for_each(|sample| {
                                    y_value[sample] = self.process1_sticky_abs(
                                        &mut states[channel],
                                        x[channel][sample],
                                    );
                                });
                            } else {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_sticky_abs(
                                        &mut states[channel],
                                        x[channel][sample],
                                    );
                                });
                            }
                        });
                    }
                    StickyMode::Rel => {
                        (0..N_CHANNELS).for_each(|channel| {
                            if let Some(y_value) = y_values[channel].as_deref_mut() {
                                (0..n_samples).for_each(|sample| {
                                    y_value[sample] = self.process1_sticky_rel(
                                        &mut states[channel],
                                        x[channel][sample],
                                    );
                                });
                            } else {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_sticky_rel(
                                        &mut states[channel],
                                        x[channel][sample],
                                    );
                                });
                            }
                        });
                    }
                }
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    if let Some(y_value) = y_values[channel].as_deref_mut() {
                        (0..n_samples).for_each(|sample| {
                            y_value[sample] =
                                self.process1(&mut states[channel], x[channel][sample]);
                        });
                    } else {
                        (0..n_samples).for_each(|sample| {
                            self.process1(&mut states[channel], x[channel][sample]);
                        });
                    }
                });
            }
        } else if self.is_asym() {
            if self.is_sticky() {
                match self.sticky_mode {
                    StickyMode::Abs => {
                        (0..N_CHANNELS).for_each(|channel| {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym_sticky_abs(
                                    &mut states[channel],
                                    x[channel][sample],
                                );
                            });
                        });
                    }
                    StickyMode::Rel => {
                        (0..N_CHANNELS).for_each(|channel| {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym_sticky_rel(
                                    &mut states[channel],
                                    x[channel][sample],
                                );
                            });
                        });
                    }
                }
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    (0..n_samples).for_each(|sample| {
                        self.process1_asym(&mut states[channel], x[channel][sample]);
                    });
                });
            }
        } else if self.is_sticky() {
            match self.sticky_mode {
                StickyMode::Abs => {
                    (0..N_CHANNELS).for_each(|channel| {
                        (0..n_samples).for_each(|sample| {
                            self.process1_sticky_abs(&mut states[channel], x[channel][sample]);
                        });
                    });
                }
                StickyMode::Rel => {
                    (0..N_CHANNELS).for_each(|channel| {
                        (0..n_samples).for_each(|sample| {
                            self.process1_sticky_rel(&mut states[channel], x[channel][sample]);
                        });
                    });
                }
            }
        } else {
            (0..N_CHANNELS).for_each(|channel| {
                (0..n_samples).for_each(|sample| {
                    self.process1(&mut states[channel], x[channel][sample]);
                });
            });
        }
    }
    /// Sets both the upgoing (attack) and downgoing (decay) cutoff frequency to `value` (Hz) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: cutoff frequency in Hz. Must be positive.
    ///
    /// # Notes
    /// Equivalent to calling `set_cutoff_up(value)` and `set_cutoff_down(value)`.
    /// Can also be derived from `set_tau(value)` as `1 / (2 * pi * value)` (net of numerical errors).
    /// Default: `INFINITY`.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        self.set_cutoff_up(value);
        self.set_cutoff_down(value);
    }
    /// Sets the upgoing (attack) cutoff frequency to `value` (Hz) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: cutoff frequency in Hz. Must be positive.
    ///
    /// # Notes
    /// Equivalent to calling `set_tau_up(value)` as `1 / (2 * pi * value)` (net of numerical errors).
    /// Default: `INFINITY`.
    #[inline(always)]
    pub fn set_cutoff_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        if self.cutoff_up != value {
            self.cutoff_up = value;
            self.param_changed |= ParamChanged::CUTOFF_UP;
        }
    }
    /// Sets the downgoing (decay) cutoff frequency to `value` (Hz) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: cutoff frequency in Hz. Must be positive.
    ///
    /// # Notes
    /// Equivalent to calling `set_tau_down(value)` as `1 / (2 * pi * value)` (net of numerical errors).
    /// Default: `INFINITY`.
    #[inline(always)]
    pub fn set_cutoff_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        if self.cutoff_down != value {
            self.cutoff_down = value;
            self.param_changed |= ParamChanged::CUTOFF_DOWN;
        }
    }
    /// Sets both the upgoing (attack) and downgoing (decay) time constants to `value` (seconds) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: time constant in seconds. Must be non-negative.
    ///
    /// # Notes
    /// Equivalent to calling `set_tau_up(value)` and `set_tau_down(value)`, or `set_cutoff(1 / (2 * pi * value))` (net of numerical errors).
    /// Default: `0.0`.
    #[inline(always)]
    pub fn set_tau(&mut self, value: f32) {
        self.set_tau_up(value);
        self.set_tau_down(value);
    }
    /// Sets the upgoing (attack) time constant to `value` (seconds) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: time constant in seconds. Must be non-negative.
    ///
    /// # Notes
    /// Equivalent to `set_cutoff_up(1 / (2 * pi * value))` (net of numerical errors).
    /// Default: `0.0`.
    #[inline(always)]
    pub fn set_tau_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        let cutoff = if value < NANO {
            f32::INFINITY
        } else {
            INVERSE_2_PI * rcpf(value)
        };
        self.set_cutoff_up(cutoff);
    }
    /// Sets the downgoing (decay) time constant to `value` (seconds) in the coefficients.
    ///
    /// # Parameters
    /// - `value`: time constant in seconds. Must be non-negative.
    ///
    /// # Notes
    /// Equivalent to `set_cutoff_down(1 / (2 * pi * value))` (net of numerical errors).
    /// Default: `0.0`.
    #[inline(always)]
    pub fn set_tau_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        let cutoff = if value < NANO {
            f32::INFINITY
        } else {
            INVERSE_2_PI * rcpf(value)
        };
        self.set_cutoff_down(cutoff);
    }
    /// Sets the target-reach threshold in the coefficients.
    ///
    /// # Parameters
    /// - `value`: threshold value in the range `[0.0, 1e18]`.
    ///
    /// # Notes
    /// When the difference between the output and input falls below this threshold
    /// according to the current sticky mode, the output is forced to equal the input.
    /// Default: `0.0`.
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_range(0.0..=1.0e18, value);
        if self.sticky_thresh != value {
            self.sticky_thresh = value;
            self.param_changed |= ParamChanged::STICKY_TRESH;
        }
    }
    /// Sets the current distance metric for sticky behavior in the coefficients.
    ///
    /// # Parameters
    /// - `value`: the sticky mode to use.
    ///
    /// # Notes
    /// Default: `StickyMode::Abs`.
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: StickyMode) {
        self.sticky_mode = value;
    }
    /// Returns the current target-reach threshold from the coefficients.
    #[inline(always)]
    pub fn get_sticky_thresh(&self) -> f32 {
        self.sticky_thresh
    }
    /// Returns the current sticky mode from the coefficients.
    #[inline(always)]
    pub fn get_sticky_mode(&self) -> StickyMode {
        self.sticky_mode
    }
    /// Returns the last output sample stored in the state.
    ///
    /// # Parameters
    /// - `state`: the `OnePoleState` to query.
    ///
    /// Returns the last output sample as `f32`.
    #[inline(always)]
    pub fn get_y_z1(&self, state: OnePoleState) -> f32 {
        state.get_y_z1()
    }
    /// Checks whether the coefficients are valid.
    ///
    /// # Notes
    /// <div class="warning">Not implemented yet!</div>
    #[inline(always)]
    pub fn coeffs_is_valid(&self) {
        todo!()
    }
    /// Checks whether the given state is valid.
    ///
    /// # Parameters
    /// - `state`: the `OnePoleState` to check.
    #[inline(always)]
    pub fn state_is_valid(&self, state: OnePoleState) -> bool {
        state.get_y_z1().is_finite()
    }

    // Private
    #[inline(always)]
    fn do_update_coeffs_ctrl(&mut self) {
        if !self.param_changed.is_empty() {
            if self.param_changed.contains(ParamChanged::CUTOFF_UP) {
                // tau < 1 ns is instantaneous for any practical purpose
                self.m_a1u = if self.cutoff_up > 1.591_549_4e8 {
                    0.0
                } else {
                    self.fs_2pi * rcpf(self.fs_2pi + self.cutoff_up)
                }
            }
            if self.param_changed.contains(ParamChanged::CUTOFF_DOWN) {
                // as before
                self.m_a1d = if self.cutoff_down > 1.591_549_4e8 {
                    0.0
                } else {
                    self.fs_2pi * rcpf(self.fs_2pi + self.cutoff_down)
                }
            }
            if self.param_changed.contains(ParamChanged::STICKY_TRESH) {
                self.st2 = self.sticky_thresh * self.sticky_thresh;
            }
            self.param_changed = ParamChanged::empty()
        }
    }
    // Helper functions for clarity sake
    #[inline(always)]
    fn is_asym(&self) -> bool {
        self.m_a1d != self.m_a1u
    }
    #[inline(always)]
    fn is_sticky(&self) -> bool {
        self.st2 != 0.
    }
}

impl<const N_CHANNELS: usize> Default for OnePoleCoeffs<N_CHANNELS> {
    #[inline(always)]
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::f32;
    #[cfg(debug_assertions)]
    use std::panic;

    use super::*;
    use crate::c_wrapper::{
        one_pole::OnePole as OnePoleWrapper, one_pole::StickyMode as StikyModeWrapper, *,
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    type OnePoleT = OnePole<N_CHANNELS>;
    type OnePoleCoeffsT = OnePoleCoeffs<N_CHANNELS>;
    type OnePoleWrapperT = OnePoleWrapper<N_CHANNELS>;

    #[test]
    fn one_pole_initialization() {
        let c_one_pole = OnePoleWrapperT::new();
        let rust_one_pole = OnePoleT::new();

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn set_sample_rate() {
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        assert_eq!(rust_one_pole.coeffs.fs_2pi, c_one_pole.coeffs.fs_2pi)
    }

    #[test]
    fn set_cutoff() {
        const CUTOFF: f32 = 1000.0;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_cutoff(CUTOFF);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_cutoff_negative() {
        const CUTOFF: f32 = -1.0;
        let mut rust_one_pole = OnePoleT::new();

        rust_one_pole.set_cutoff(CUTOFF);
    }

    #[test]
    fn set_cutoff_up() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_up(CUTOFF);
        rust_one_pole.set_cutoff_up(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_up, CUTOFF);
        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn set_cutoff_down() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_down(CUTOFF);
        rust_one_pole.set_cutoff_down(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_down, CUTOFF);
        assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
    }

    #[test]
    fn set_tau() {
        const CUTOFF: f32 = 150.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn set_tau_up() {
        const CUTOFF: f32 = 1500.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau_up(TAU);
        rust_one_pole.set_tau_up(TAU);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn set_tau_down() {
        const CUTOFF: f32 = 10_000.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau_down(TAU);
        rust_one_pole.set_tau_down(TAU);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        const TAU: f32 = 0.1e-9;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_tau_negative() {
        let mut rust_one_pole = OnePoleT::new();

        rust_one_pole.set_tau(-1.);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.set_sticky_thresh(sticky_tresh);
        rust_one_pole.set_sticky_thresh(sticky_tresh);

        assert!(rust_one_pole.coeffs.param_changed.bits() & BW_ONE_POLE_PARAM_STICKY_THRESH != 0);
        assert_eq!(rust_one_pole.get_sticky_thresh(), sticky_tresh);
        assert_eq!(
            rust_one_pole.get_sticky_thresh(),
            c_one_pole.get_sticky_thresh()
        );
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got -1e0")]
    fn set_sticky_tresh_negative() {
        let mut rust_one_pole = OnePoleT::new();

        rust_one_pole.set_sticky_thresh(-1.);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got 1.1e18")]
    fn set_sticky_tresh_too_high() {
        let mut rust_one_pole = OnePoleT::new();

        rust_one_pole.set_sticky_thresh(1.1e18);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let c_mode = StikyModeWrapper::Abs;
        let rust_mode = StickyMode::Abs;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.set_sticky_mode(c_mode);
        rust_one_pole.set_sticky_mode(rust_mode);

        assert_eq!(
            rust_one_pole.get_sticky_mode() as u32,
            c_one_pole.get_sticky_mode() as u32
        );
        assert_eq!(rust_one_pole.get_sticky_mode(), rust_mode);
    }

    #[test]
    fn set_sticky_mode_rel() {
        let c_mode = StikyModeWrapper::Rel;
        let rust_mode = StickyMode::Rel;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.set_sticky_mode(c_mode);
        rust_one_pole.set_sticky_mode(rust_mode);

        assert_eq!(
            rust_one_pole.get_sticky_mode() as u32,
            c_one_pole.get_sticky_mode() as u32
        );
        assert_eq!(rust_one_pole.get_sticky_mode(), rust_mode);
    }

    #[test]
    fn reset_none() {
        const CUTOFF: f32 = 1200.0;
        let sticky_tresh = 0.1;
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();
        let x0 = [0.5; N_CHANNELS];

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_tresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_tresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(&x0, None);
        rust_one_pole.reset_multi(&x0, None);

        assert_one_pole(&rust_one_pole, &c_one_pole);
    }

    #[test]
    fn reset_with_input_and_output() {
        const CUTOFF: f32 = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N_CHANNELS];
        let mut c_y0_output = [0.0; N_CHANNELS];
        let mut rust_y0_output = [0.0; N_CHANNELS];
        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_thresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_thresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(&x0_input, Some(&mut c_y0_output));
        rust_one_pole.reset_multi(&x0_input, Some(&mut rust_y0_output));

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.states[i].y_z1, x0_input[i]);
            assert_eq!(rust_one_pole.states[i].y_z1, c_one_pole.states[i].y_z1);
            assert_eq!(rust_y0_output[i], x0_input[i]);
        }
    }
    #[test]
    fn process1() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF: f32 = 100.;
        const STICKY_TRESH: f32 = 0.;
        let x0 = [45.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_abs() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF: f32 = 100_000.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [1.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(StickyMode::Abs);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(StikyModeWrapper::Abs);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1_sticky_abs(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1_sticky_abs(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_rel() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF: f32 = 100_000.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [1.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(StickyMode::Rel);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(StikyModeWrapper::Rel);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1_sticky_rel(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1_sticky_rel(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_asym() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF_UP: f32 = 100.;
        const CUTOFF_DOWN: f32 = 300.;
        const STICKY_TRESH: f32 = 0.;
        let x0 = [0.5, 0.32];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff_up(CUTOFF_UP);
        rust_one_pole.set_cutoff_down(CUTOFF_DOWN);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1_asym(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1_asym(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_asym_sticy_abs() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF_UP: f32 = 100.;
        const CUTOFF_DOWN: f32 = 300.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [0.5, 0.32];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff_up(CUTOFF_UP);
        rust_one_pole.set_cutoff_down(CUTOFF_DOWN);
        rust_one_pole.set_sticky_mode(StickyMode::Abs);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(StikyModeWrapper::Abs);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1_asym_sticky_abs(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1_asym_sticky_abs(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_asym_sticky_rel() {
        let mut rust_one_pole = OnePoleT::new();
        let mut c_one_pole = OnePoleWrapperT::new();
        const CUTOFF_UP: f32 = 100.;
        const CUTOFF_DOWN: f32 = 300.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [0.5, 0.32];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff_up(CUTOFF_UP);
        rust_one_pole.set_cutoff_down(CUTOFF_DOWN);
        rust_one_pole.set_sticky_mode(StickyMode::Rel);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset_multi(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(StikyModeWrapper::Rel);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole
                .coeffs
                .process1_asym_sticky_rel(&mut c_one_pole.states[channel], x0[channel]);
            rust_y[channel] = rust_one_pole
                .coeffs
                .process1_asym_sticky_rel(&mut rust_one_pole.states[channel], x0[channel]);

            assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process_with_y() {
        const N_CHANNELS: usize = 2;
        const N_SAMPLES: usize = 4;
        const CUTOFF: f32 = 1000.0;
        const SAMPLE_RATE: f32 = 48000.0;

        let x0: [&[f32]; N_CHANNELS] = [&[1.0, 2.0, 3.0, 4.0], &[0.5, 1.5, 2.5, 3.5]];

        let mut c_buf0 = [0.0; N_SAMPLES];
        let mut c_buf1 = [0.0; N_SAMPLES];
        let mut rust_buf0 = [0.0; N_SAMPLES];
        let mut rust_buf1 = [0.0; N_SAMPLES];

        let mut c_output: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut c_buf0), Some(&mut c_buf1)];
        let mut rust_output: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut rust_buf0), Some(&mut rust_buf1)];

        let mut c_one_pole = OnePoleWrapperT::new();
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.reset(&[0.0; N_CHANNELS], None);

        c_one_pole.process(&x0, Some(&mut c_output), N_SAMPLES);

        let mut rust_one_pole = OnePoleT::new();
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.reset_multi(&[0.0; N_CHANNELS], None);

        rust_one_pole.process(&x0, Some(&mut rust_output), N_SAMPLES);

        for ch in 0..N_CHANNELS {
            for sample in 0..N_SAMPLES {
                assert_eq!(
                    rust_output[ch].as_ref().unwrap()[sample],
                    c_output[ch].as_ref().unwrap()[sample]
                );
                println!(
                    "C output: {}\nRust output: {}",
                    c_output[ch].as_ref().unwrap()[sample],
                    rust_output[ch].as_ref().unwrap()[sample]
                );
            }
        }
    }

    #[test]
    fn get_y_z1() {
        const N_CHANNELS: usize = 2;
        const SAMPLE_RATE: f32 = 44100.0;
        const CUTOFF: f32 = 1000.0;
        const N_SAMPLES: usize = 4;

        let mut c_one_pole = OnePoleWrapperT::new();
        let mut rust_one_pole = OnePoleT::new();

        let x0: [f32; N_CHANNELS] = [0.0, 0.0];

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_mode(StikyModeWrapper::Rel);
        c_one_pole.reset(&x0, None);

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(StickyMode::Rel);
        c_one_pole.reset(&x0, None);

        let input: [&[f32]; N_CHANNELS] = [&[1.0, 2.0, 3.0, 4.0], &[0.5, 1.5, 2.5, 3.5]];

        let mut c_buf0 = [0.0, 0.1, 0.2, 0.3];
        let mut c_buf1 = [1.0, 1.1, 1.2, 1.3];
        let mut rust_buf0 = [0.0, 0.1, 0.2, 0.3];
        let mut rust_buf1 = [1.0, 1.1, 1.2, 1.3];

        let mut c_output: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut c_buf0), Some(&mut c_buf1)];
        let mut rust_output: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut rust_buf0), Some(&mut rust_buf1)];

        c_one_pole.process(&input, Some(&mut c_output), N_SAMPLES);

        rust_one_pole.process(&input, Some(&mut rust_output), N_SAMPLES);

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.get_y_z1(i), c_one_pole.get_y_z1(i));
            assert_eq!(
                rust_one_pole.get_y_z1(i),
                rust_output[i].as_ref().unwrap()[3]
            );
        }
    }

    // Do we need c sanity checks??
    // By design I can not insert a non valid value
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got NaN")]
    fn set_cutoff_to_nan() {
        let result = panic::catch_unwind(|| {
            let mut one_pole_coeffs = OnePoleCoeffsT::default();
            one_pole_coeffs.set_cutoff_down(f32::NAN);
        });
        assert!(result.is_err());

        let result = panic::catch_unwind(|| {
            let mut one_pole_coeffs = OnePoleCoeffsT::default();
            one_pole_coeffs.set_cutoff_up(f32::NAN);
        });
        assert!(result.is_err());

        let mut one_pole_coeffs = OnePoleCoeffsT::default();
        one_pole_coeffs.set_cutoff(f32::NAN);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got inf")]
    fn set_sticky_thresh_to_infinite() {
        let mut one_pole_coeffs = OnePoleCoeffsT::default();

        one_pole_coeffs.set_sticky_thresh(f32::INFINITY);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_sample_rate_must_be_finite() {
        let mut one_pole_coeffs = OnePoleCoeffsT::default();

        one_pole_coeffs.set_sample_rate(f32::INFINITY);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_sample_rate_must_be_positive() {
        let mut rust_one_pole = OnePoleT::new();

        rust_one_pole.set_sample_rate(-1.);
    }

    pub(crate) fn assert_one_pole_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &OnePoleCoeffs<N_CHANNELS>,
        c_coeffs: &bw_one_pole_coeffs,
    ) {
        let pre_message = "one_pole.coeff.";
        let post_message = "does not match";
        assert_eq!(
            rust_coeffs.fs_2pi, c_coeffs.fs_2pi,
            "{}fs_2pi {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.m_a1u, c_coeffs.mA1u,
            "{}mA1u {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.m_a1d, c_coeffs.mA1d,
            "{}mA1d {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.st2, c_coeffs.st2,
            "{}st2 {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.cutoff_up, c_coeffs.cutoff_up,
            "{}cutoff_up {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.cutoff_down, c_coeffs.cutoff_down,
            "{}cutoff_down {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.sticky_thresh, c_coeffs.sticky_thresh,
            "{}sticky_thresh {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.sticky_mode as u32, c_coeffs.sticky_mode,
            "{}sticky_mode {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.param_changed.bits(),
            c_coeffs.param_changed as u32,
            // to check only the used flags, see ParamChanged struct
            // c_coeffs.param_changed as u32 & (0b111),
            "{}param_changed {}",
            pre_message,
            post_message
        );
    }

    pub(crate) fn assert_one_pole_state(rust_states: &OnePoleState, c_states: &bw_one_pole_state) {
        let pre_message = "one_pole.states";
        let post_message = "does not match";
        assert_eq!(
            rust_states.y_z1, c_states.y_z1,
            "{pre_message}.y_z1 {post_message}"
        )
    }

    pub(crate) fn assert_one_pole<const N_CHANNELS: usize>(
        rust_one_pole: &OnePole<N_CHANNELS>,
        c_one_pole: &OnePoleWrapper<N_CHANNELS>,
    ) {
        assert_one_pole_coeffs(&rust_one_pole.coeffs, &c_one_pole.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_one_pole_state(&rust_one_pole.states[channel], &c_one_pole.states[channel]);
        });
    }
}
