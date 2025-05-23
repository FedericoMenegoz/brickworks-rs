use bitflags::bitflags;

use super::math::rcpf;
use crate::native::common::{INVERSE_2_PI, NANO};
#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};

/// One-pole filter: Rust native port of the original implementation by [Orastron](https://www.orastron.com/algorithms/bw_one_pole).
///
/// One-pole (6 dB/oct) lowpass filter with unitary DC gain, separate attack and decay time constants, and sticky target-reach threshold.
///
/// This is better suited to implement smoothing than bw_lp1.
///
/// # Example
/// ```
/// use brickworks_rs::native::one_pole::*;
///
/// const CUTOFF: f32 = 100_000.0;
/// const STICKY_THRESH: f32 = 0.9;
/// const SAMPLE_RATE: f32 = 48_100.0;
/// const N_CHANNELS: usize = 2;
/// const N_SAMPLES: usize = 1;
/// fn main() {
///     // Create a new OnePole filter instance for N_CHANNELS
///     let mut one_pole = OnePole::<N_CHANNELS>::new();
///
///     // Input signal: one sample per channel
///     let x = vec![vec![1.0], vec![33.0]];
///
///     // Output buffer, same shape as input
///     let mut y = vec![vec![0.0], vec![0.0]];
///
///     // Configure the filter
///     one_pole.set_sample_rate(SAMPLE_RATE);
///     one_pole.set_cutoff(CUTOFF);
///     one_pole.set_sticky_mode(OnePoleStickyMode::Rel);
///     one_pole.set_sticky_thresh(STICKY_THRESH);
///
///     // Initialize the filter state for each channel
///     one_pole.reset(&[0.0, 0.0], None);
///     // !!! This is ugly yeah...need to investigate
///     let mut y_wrapped: Vec<Option<&mut [f32]>> =
///         y.iter_mut().map(|ch| Some(&mut **ch)).collect();
///
///     // Process one sample per channel
///     one_pole.process(&x, Some(&mut y_wrapped), N_SAMPLES);
///
///     // Output the filtered result
///     println!("Filtered output: {:?}", y);
/// }
///
/// ```
/// 
///# Notes
/// This module provides a native Rust implementation of the filter, but the same interface is 
/// also available via bindings to the original C library at [crate::c_wrapper::one_pole].
/// 
#[derive(Debug)]
pub struct OnePole<const N_CHANNELS: usize> {
    coeffs: OnePoleCoeffs,
    states: Vec<OnePoleState>,
    _states_p: Vec<OnePoleState>, // BW_RESTRICT to check what is for
}
/// Distance metrics for sticky behavior.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OnePoleStickyMode {
    /// `OnePoleStickyMode::Abs`: absolute difference ( `|out - in|` );
    Abs,
    /// `OnePoleStickyMode::Rel`: relative difference with respect to input (`|out - in| / |in|`).
    Rel,
}

impl<const N_CHANNELS: usize> OnePole<N_CHANNELS> {
    /// Creates a new `OnePole` filter with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        OnePole {
            coeffs: Default::default(),
            states: vec![OnePoleState { y_z1: 0.0 }; N_CHANNELS],
            _states_p: vec![OnePoleState { y_z1: 0.0 }; N_CHANNELS],
        }
    }
    /// Sets the filter's sample rate, required for accurate time constant conversion.
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }
    /// Resets the filter state for all channels using the initial input `x0`.
    /// If `y0` is provided, the resulting initial outputs are stored in it.
    #[inline(always)]
    pub fn reset(&mut self, x0: &[f32], y0: Option<&mut [f32]>) {
        self.coeffs.reset_coeffs();
        self.reset_state_multi(x0, y0);
    }
    /// Processes a block of input samples using a one-pole filter per channel.
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
        x: &[Vec<f32>],
        y: Option<&mut [Option<&mut [f32]>]>,
        n_samples: usize,
    ) {
        self.process_multi(x, y, n_samples);
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
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
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
    /// to the current distance metric (see [OnePoleStickyMode]), the output is forcefully set to be equal to the input value.
    ///
    ///
    /// ### Parameters
    /// - `value`: sticky threshhold, default is `f32::0.0` and its valid range is [0.f, 1e18f].
    ///
    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        self.coeffs.set_sticky_thresh(value);
    }
    /// Sets the current distance metric for sticky behavior to value in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: sticky mode, default is [OnePoleStickyMode::Abs]
    ///
    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: OnePoleStickyMode) {
        self.coeffs.set_sticky_mode(value);
    }
    /// Returns the current target-reach threshold in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_sticky_thresh(&self) -> f32 {
        self.coeffs.get_sticky_thresh()
    }
    /// Returns the current distance metric for sticky behavior in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_sticky_mode(&self) -> OnePoleStickyMode {
        self.coeffs.get_sticky_mode()
    }
    /// Returns the last output sample as stored in `OnePole<N_CHANNELS>`.
    #[inline(always)]
    pub fn get_y_z1(&self, channel: usize) -> f32 {
        self.states[channel].y_z1
    }

    // private methods
    #[inline(always)]
    fn reset_state_multi(&mut self, x0: &[f32], y0: Option<&mut [f32]>) {
        // No need to check states are in different addresses cause it is
        // enforced by design
        if let Some(out) = y0 {
            (0..N_CHANNELS).for_each(|channel| {
                out[channel] = self.states[channel].reset(x0[channel]);
            });
        } else {
            (0..N_CHANNELS).for_each(|channel| {
                self.states[channel].reset(x0[channel]);
            });
        }
    }

    #[inline(always)]
    fn process1(&mut self, x: f32, channel: usize) -> f32 {
        let y = x + self.coeffs.m_a1u * (self.get_y_z1(channel) - x);
        self.states[channel].y_z1 = y;
        y
    }

    #[inline(always)]
    fn process1_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        let mut y = x + self.coeffs.m_a1u * (self.get_y_z1(channel) - x);
        let d = y - x;
        if d * d <= self.coeffs.st2 {
            y = x
        }
        self.states[channel].y_z1 = y;

        y
    }

    #[inline(always)]
    fn process1_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        let mut y = x + self.coeffs.m_a1u * (self.get_y_z1(channel) - x);
        let d = y - x;
        if d * d <= self.coeffs.st2 * x * x {
            y = x
        }
        self.states[channel].y_z1 = y;

        y
    }

    #[inline(always)]
    fn process1_asym(&mut self, x: f32, channel: usize) -> f32 {
        let y_z1 = self.get_y_z1(channel);
        let ma1 = if x >= y_z1 {
            self.coeffs.m_a1u
        } else {
            self.coeffs.m_a1d
        };

        let y = x + ma1 * (y_z1 - x);
        self.states[channel].y_z1 = y;

        y
    }

    #[inline(always)]
    fn process1_asym_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        let y_z1 = self.get_y_z1(channel);
        let ma1 = if x >= y_z1 {
            self.coeffs.m_a1u
        } else {
            self.coeffs.m_a1d
        };

        let mut y = x + ma1 * (y_z1 - x);
        let d = y - x;
        if d * d <= self.coeffs.st2 {
            y = x;
        }
        self.states[channel].y_z1 = y;

        y
    }

    #[inline(always)]
    fn process1_asym_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        let y_z1 = self.get_y_z1(channel);
        let ma1 = if x >= y_z1 {
            self.coeffs.m_a1u
        } else {
            self.coeffs.m_a1d
        };

        let mut y = x + ma1 * (y_z1 - x);
        let d = y - x;
        if d * d <= self.coeffs.st2 * x * x {
            y = x;
        }
        self.states[channel].y_z1 = y;

        y
    }

    #[inline(always)]
    fn process_multi(
        &mut self,
        x: &[Vec<f32>],
        y: Option<&mut [Option<&mut [f32]>]>,
        n_samples: usize,
    ) {
        // As for reset_multi no need to check states are in
        // different addresses cause it is enforced by design
        // Still need to investigate if the other assertions
        // are sanity check only needed in c

        self.coeffs.update_coeffs_ctrl();
        if let Some(y_values) = y {
            if self.is_asym() {
                if self.is_sticky() {
                    match self.coeffs.sticky_mode {
                        OnePoleStickyMode::Abs => {
                            (0..N_CHANNELS).for_each(|channel| {
                                if let Some(y_value) = y_values[channel].as_deref_mut() {
                                    (0..n_samples).for_each(|sample| {
                                        y_value[sample] = self
                                            .process1_asym_sticky_abs(x[channel][sample], channel);
                                    });
                                } else {
                                    (0..n_samples).for_each(|sample| {
                                        self.process1_asym_sticky_abs(x[channel][sample], channel);
                                    });
                                }
                            });
                        }
                        OnePoleStickyMode::Rel => {
                            (0..N_CHANNELS).for_each(|channel| {
                                if let Some(y_value) = y_values[channel].as_deref_mut() {
                                    (0..n_samples).for_each(|sample| {
                                        y_value[sample] = self
                                            .process1_asym_sticky_rel(x[channel][sample], channel);
                                    });
                                } else {
                                    (0..n_samples).for_each(|sample| {
                                        self.process1_asym_sticky_rel(x[channel][sample], channel);
                                    });
                                }
                            });
                        }
                    }
                } else {
                    (0..N_CHANNELS).for_each(|channel| {
                        if let Some(y_value) = y_values[channel].as_deref_mut() {
                            (0..n_samples).for_each(|sample| {
                                y_value[sample] = self.process1_asym(x[channel][sample], channel);
                            });
                        } else {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym(x[channel][sample], channel);
                            });
                        }
                    });
                }
            } else if self.is_sticky() {
                match self.coeffs.sticky_mode {
                    OnePoleStickyMode::Abs => {
                        (0..N_CHANNELS).for_each(|channel| {
                            if let Some(y_value) = y_values[channel].as_deref_mut() {
                                (0..n_samples).for_each(|sample| {
                                    y_value[sample] =
                                        self.process1_sticky_abs(x[channel][sample], channel);
                                });
                            } else {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_sticky_abs(x[channel][sample], channel);
                                });
                            }
                        });
                    }
                    OnePoleStickyMode::Rel => {
                        (0..N_CHANNELS).for_each(|channel| {
                            if let Some(y_value) = y_values[channel].as_deref_mut() {
                                (0..n_samples).for_each(|sample| {
                                    y_value[sample] =
                                        self.process1_sticky_rel(x[channel][sample], channel);
                                });
                            } else {
                                (0..n_samples).for_each(|sample| {
                                    self.process1_sticky_rel(x[channel][sample], channel);
                                });
                            }
                        });
                    }
                }
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    if let Some(y_value) = y_values[channel].as_deref_mut() {
                        (0..n_samples).for_each(|sample| {
                            y_value[sample] = self.process1(x[channel][sample], channel);
                        });
                    } else {
                        (0..n_samples).for_each(|sample| {
                            self.process1(x[channel][sample], channel);
                        });
                    }
                });
            }
        } else if self.is_asym() {
            if self.is_sticky() {
                match self.coeffs.sticky_mode {
                    OnePoleStickyMode::Abs => {
                        (0..N_CHANNELS).for_each(|channel| {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym_sticky_abs(x[channel][sample], channel);
                            });
                        });
                    }
                    OnePoleStickyMode::Rel => {
                        (0..N_CHANNELS).for_each(|channel| {
                            (0..n_samples).for_each(|sample| {
                                self.process1_asym_sticky_rel(x[channel][sample], channel);
                            });
                        });
                    }
                }
            } else {
                (0..N_CHANNELS).for_each(|channel| {
                    (0..n_samples).for_each(|sample| {
                        self.process1_asym(x[channel][sample], channel);
                    });
                });
            }
        } else if self.is_sticky() {
            match self.coeffs.sticky_mode {
                OnePoleStickyMode::Abs => {
                    (0..N_CHANNELS).for_each(|channel| {
                        (0..n_samples).for_each(|sample| {
                            self.process1_sticky_abs(x[channel][sample], channel);
                        });
                    });
                }
                OnePoleStickyMode::Rel => {
                    (0..N_CHANNELS).for_each(|channel| {
                        (0..n_samples).for_each(|sample| {
                            self.process1_sticky_rel(x[channel][sample], channel);
                        });
                    });
                }
            }
        } else {
            (0..N_CHANNELS).for_each(|channel| {
                (0..n_samples).for_each(|sample| {
                    self.process1(x[channel][sample], channel);
                });
            });
        }
    }
    // Not implementing the process() defined with BW_CXX_NO_ARRAY

    // Helper functions for clarity sake
    #[inline(always)]
    fn is_asym(&self) -> bool {
        self.coeffs.m_a1d != self.coeffs.m_a1u
    }
    #[inline(always)]
    fn is_sticky(&self) -> bool {
        self.coeffs.st2 != 0.
    }
}

impl<const N_CHANNELS: usize> Default for OnePole<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug, Copy)]
struct OnePoleCoeffs {
    fs_2pi: f32,
    m_a1u: f32,
    m_a1d: f32,
    st2: f32,
    cutoff_up: f32,
    cutoff_down: f32,
    sticky_thresh: f32,
    sticky_mode: OnePoleStickyMode,
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

#[derive(Debug, Clone, Copy)]
struct OnePoleState {
    y_z1: f32,
}

impl OnePoleState {
    #[inline(always)]
    fn reset(&mut self, x0: f32) -> f32 {
        self.y_z1 = x0;
        x0
    }
}

impl OnePoleCoeffs {
    #[inline(always)]
    fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_positive(sample_rate);
            debug_assert_is_finite(sample_rate);
        }
        self.fs_2pi = INVERSE_2_PI * sample_rate;
    }

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

    #[inline(always)]
    fn reset_coeffs(&mut self) {
        self.param_changed = ParamChanged::all();
        self.do_update_coeffs_ctrl();
    }

    #[inline(always)]
    fn update_coeffs_ctrl(&mut self) {
        self.do_update_coeffs_ctrl();
    }

    // This is only asserting
    // #[inline(always)]
    // fn update_coeffs_audio(&self) {
    //     todo!()
    // }

    #[inline(always)]
    fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        self.set_cutoff_up(value);
        self.set_cutoff_down(value);
    }

    #[inline(always)]
    fn set_cutoff_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        if self.cutoff_up != value {
            self.cutoff_up = value;
            self.param_changed |= ParamChanged::CUTOFF_UP;
        }
    }

    #[inline(always)]
    fn set_cutoff_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        if self.cutoff_down != value {
            self.cutoff_down = value;
            self.param_changed |= ParamChanged::CUTOFF_DOWN;
        }
    }

    #[inline(always)]
    fn set_tau(&mut self, value: f32) {
        self.set_tau_up(value);
        self.set_tau_down(value);
    }

    #[inline(always)]
    fn set_tau_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        let cutoff = if value < NANO {
            f32::INFINITY
        } else {
            INVERSE_2_PI * rcpf(value)
        };
        self.set_cutoff_up(cutoff);
    }

    #[inline(always)]
    fn set_tau_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        let cutoff = if value < NANO {
            f32::INFINITY
        } else {
            INVERSE_2_PI * rcpf(value)
        };
        self.set_cutoff_down(cutoff);
    }

    #[inline(always)]
    fn set_sticky_thresh(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_range(0., 1.0e18, value);
        if self.sticky_thresh != value {
            self.sticky_thresh = value;
            self.param_changed |= ParamChanged::STICKY_TRESH;
        }
    }

    #[inline(always)]
    fn set_sticky_mode(&mut self, value: OnePoleStickyMode) {
        self.sticky_mode = value;
    }

    #[inline(always)]
    fn get_sticky_thresh(&self) -> f32 {
        self.sticky_thresh
    }

    #[inline(always)]
    fn get_sticky_mode(&self) -> OnePoleStickyMode {
        self.sticky_mode
    }
}

impl Default for OnePoleCoeffs {
    #[inline(always)]
    fn default() -> Self {
        Self {
            fs_2pi: Default::default(),
            m_a1u: Default::default(),
            m_a1d: Default::default(),
            st2: Default::default(),
            cutoff_up: f32::INFINITY,
            cutoff_down: f32::INFINITY,
            sticky_thresh: 0.0,
            sticky_mode: OnePoleStickyMode::Abs,
            param_changed: ParamChanged::all(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f32;
    #[cfg(debug_assertions)]
    use std::panic;

    use super::*;
    use crate::c_wrapper::{one_pole::OnePole as OnePoleWrapper, one_pole::OnePoleStickyMode as WrapperStikyMode, *};

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    #[test]
    fn update_coeffs_ctrl_all_changed() {
        let cutoff = 100.;
        let sticky_thresh = 0.01;

        let mut rust_coeffs = OnePoleCoeffs {
            cutoff_up: cutoff,
            cutoff_down: cutoff,
            sticky_thresh: sticky_thresh,
            param_changed: ParamChanged::all(),
            ..Default::default()
        };
        let mut c_coeffs = bw_one_pole_coeffs {
            cutoff_up: Default::default(),
            cutoff_down: Default::default(),
            sticky_thresh: Default::default(),
            param_changed: Default::default(),
            fs_2pi: Default::default(),
            mA1u: Default::default(),
            mA1d: Default::default(),
            st2: Default::default(),
            sticky_mode: Default::default(),
        };

        unsafe {
            bw_one_pole_init(&mut c_coeffs);
            c_coeffs.cutoff_up = cutoff;
            c_coeffs.cutoff_down = cutoff;
            c_coeffs.sticky_thresh = sticky_thresh;
            bw_one_pole_update_coeffs_ctrl(&mut c_coeffs);
        }
        rust_coeffs.update_coeffs_ctrl();

        assert_coeffs_rust_c(rust_coeffs, c_coeffs);
    }

    #[test]
    fn update_coeffs_ctrl_nothing_changed() {
        let mut rust_coeffs = OnePoleCoeffs {
            param_changed: ParamChanged::empty(),
            ..Default::default()
        };
        let mut c_coeffs = bw_one_pole_coeffs {
            cutoff_up: Default::default(),
            cutoff_down: Default::default(),
            sticky_thresh: Default::default(),
            param_changed: 0,
            fs_2pi: Default::default(),
            mA1u: Default::default(),
            mA1d: Default::default(),
            st2: Default::default(),
            sticky_mode: Default::default(),
        };

        unsafe {
            bw_one_pole_init(&mut c_coeffs);
            bw_one_pole_do_update_coeffs_ctrl(&mut c_coeffs);
        }
        rust_coeffs.do_update_coeffs_ctrl();

        assert_coeffs_rust_c(rust_coeffs, c_coeffs);
    }

    #[test]
    fn one_pole_initialization() {
        let c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let rust_one_pole = OnePole::<N_CHANNELS>::new();

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_sample_rate() {
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        assert_eq!(rust_one_pole.coeffs.fs_2pi, c_one_pole.coeffs.fs_2pi)
    }

    #[test]
    fn set_cutoff() {
        const CUTOFF: f32 = 1000.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_cutoff(CUTOFF);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_cutoff_negative() {
        const CUTOFF: f32 = -1.0;
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_cutoff(CUTOFF);
    }

    #[test]
    fn set_cutoff_up() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_up(CUTOFF);
        rust_one_pole.set_cutoff_up(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_up, CUTOFF);
        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_cutoff_down() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_down(CUTOFF);
        rust_one_pole.set_cutoff_down(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_down, CUTOFF);
        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau() {
        const CUTOFF: f32 = 150.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_up() {
        const CUTOFF: f32 = 1500.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau_up(TAU);
        rust_one_pole.set_tau_up(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_down() {
        const CUTOFF: f32 = 10_000.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau_down(TAU);
        rust_one_pole.set_tau_down(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        const TAU: f32 = 0.1e-9;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_tau_negative() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_tau(-1.);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

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
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_sticky_thresh(-1.);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got 1.1e18")]
    fn set_sticky_tresh_too_high() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_sticky_thresh(1.1e18);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let c_mode = WrapperStikyMode::Abs;
        let rust_mode = OnePoleStickyMode::Abs;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

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
        let c_mode = WrapperStikyMode::Rel;
        let rust_mode = OnePoleStickyMode::Rel;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

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
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let x0 = [0.5; N_CHANNELS];

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_tresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_tresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(&x0, None);
        rust_one_pole.reset(&x0, None);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn reset_with_input_and_output() {
        const CUTOFF: f32 = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N_CHANNELS];
        let mut c_y0_output = [0.0; N_CHANNELS];
        let mut rust_y0_output = [0.0; N_CHANNELS];
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_thresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_thresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(&x0_input, Some(&mut c_y0_output));
        rust_one_pole.reset(&x0_input, Some(&mut rust_y0_output));

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.states[i].y_z1, x0_input[i]);
            assert_eq!(rust_one_pole.states[i].y_z1, c_one_pole.states[i].y_z1);
            assert_eq!(rust_y0_output[i], x0_input[i]);
        }
    }
    #[test]
    fn process1() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        const CUTOFF: f32 = 100.;
        const STICKY_TRESH: f32 = 0.;
        let x0 = [45.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_abs() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        const CUTOFF: f32 = 100_000.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [1.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Abs);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(WrapperStikyMode::Abs);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1_sticky_abs(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1_sticky_abs(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_rel() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        const CUTOFF: f32 = 100_000.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [1.0, 33.0];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Rel);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(WrapperStikyMode::Rel);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1_sticky_rel(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1_sticky_rel(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_asym() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
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
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1_asym(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1_asym(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_sticky_asym_sticy_abs() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        const CUTOFF_UP: f32 = 100.;
        const CUTOFF_DOWN: f32 = 300.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [0.5, 0.32];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff_up(CUTOFF_UP);
        rust_one_pole.set_cutoff_down(CUTOFF_DOWN);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Abs);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(WrapperStikyMode::Abs);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1_asym_sticky_abs(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1_asym_sticky_abs(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
            assert_eq!(
                rust_one_pole.get_y_z1(channel),
                c_one_pole.get_y_z1(channel)
            );
            assert_eq!(rust_y[channel], c_y[channel]);
        });
    }

    #[test]
    fn process1_asym_sticky_rel() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        const CUTOFF_UP: f32 = 100.;
        const CUTOFF_DOWN: f32 = 300.;
        const STICKY_TRESH: f32 = 0.9;
        let x0 = [0.5, 0.32];
        let mut c_y = [0.0, 0.0];
        let mut rust_y = [0.0, 0.0];

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff_up(CUTOFF_UP);
        rust_one_pole.set_cutoff_down(CUTOFF_DOWN);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Rel);
        rust_one_pole.set_sticky_thresh(STICKY_TRESH);
        rust_one_pole.reset(&x0, None);

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff_up(CUTOFF_UP);
        c_one_pole.set_cutoff_down(CUTOFF_DOWN);
        c_one_pole.set_sticky_thresh(STICKY_TRESH);
        c_one_pole.set_sticky_mode(WrapperStikyMode::Rel);
        c_one_pole.reset(&x0, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_one_pole.process1_asym_sticky_rel(x0[channel], channel);
            rust_y[channel] = rust_one_pole.process1_asym_sticky_rel(x0[channel], channel);

            assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
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

        let x0 = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];

        let mut c_buf0 = [0.0; N_SAMPLES];
        let mut c_buf1 = [0.0; N_SAMPLES];
        let mut rust_buf0 = [0.0; N_SAMPLES];
        let mut rust_buf1 = [0.0; N_SAMPLES];

        let mut c_output: [&mut [f32]; N_CHANNELS] = [&mut c_buf0, &mut c_buf1];
        let mut rust_output: [&mut [f32]; N_CHANNELS] = [&mut rust_buf0, &mut rust_buf1];

        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.reset(&[0.0; N_CHANNELS], None);

        let mut c_output_wrapped: Vec<Option<&mut [f32]>> =
            c_output.iter_mut().map(|ch| Some(&mut **ch)).collect();
        c_one_pole.process(&x0, Some(&mut c_output_wrapped), N_SAMPLES);

        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.reset(&[0.0; N_CHANNELS], None);

        let mut output_wrapped: Vec<Option<&mut [f32]>> =
            rust_output.iter_mut().map(|ch| Some(&mut **ch)).collect();
        rust_one_pole.process(&x0, Some(&mut output_wrapped), N_SAMPLES);

        for ch in 0..N_CHANNELS {
            for sample in 0..N_SAMPLES {
                assert_eq!(rust_output[ch][sample], c_output[ch][sample]);
                println!(
                    "C output: {}\nRust output: {}",
                    c_output[ch][sample], rust_output[ch][sample]
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

        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_mode(WrapperStikyMode::Rel);

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Rel);

        let input = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];

        let mut c_buf0 = [0.0, 0.1, 0.2, 0.3];
        let mut c_buf1 = [1.0, 1.1, 1.2, 1.3];
        let mut rust_buf0 = [0.0, 0.1, 0.2, 0.3];
        let mut rust_buf1 = [1.0, 1.1, 1.2, 1.3];

        let mut c_output: [&mut [f32]; N_CHANNELS] = [&mut c_buf0, &mut c_buf1];
        let mut rust_output: [&mut [f32]; N_CHANNELS] = [&mut rust_buf0, &mut rust_buf1];

        let mut c_output_wrapped: Vec<Option<&mut [f32]>> =
            c_output.iter_mut().map(|ch| Some(&mut **ch)).collect();
        c_one_pole.process(&input, Some(&mut c_output_wrapped), N_SAMPLES);

        let mut output_wrapped: Vec<Option<&mut [f32]>> =
            rust_output.iter_mut().map(|ch| Some(&mut **ch)).collect();
        rust_one_pole.process(&input, Some(&mut output_wrapped), N_SAMPLES);

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.get_y_z1(i), c_one_pole.get_y_z1(i));
            assert_eq!(rust_one_pole.get_y_z1(i), rust_output[i][3]);
        }
    }

    // Do we need c sanity checks??
    // By design I can not insert a non valid value
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got NaN")]
    fn set_cutoff_to_nan() {
        let result = panic::catch_unwind(|| {
            let mut one_pole_coeffs = OnePoleCoeffs::default();
            one_pole_coeffs.set_cutoff_down(f32::NAN);
        });
        assert!(result.is_err());

        let result = panic::catch_unwind(|| {
            let mut one_pole_coeffs = OnePoleCoeffs::default();
            one_pole_coeffs.set_cutoff_up(f32::NAN);
        });
        assert!(result.is_err());

        let mut one_pole_coeffs = OnePoleCoeffs::default();
        one_pole_coeffs.set_cutoff(f32::NAN);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got inf")]
    fn set_sticky_thresh_to_infinite() {
        let mut one_pole_coeffs = OnePoleCoeffs::default();

        one_pole_coeffs.set_sticky_thresh(f32::INFINITY);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_sample_rate_invalid() {
        let mut one_pole_coeffs = OnePoleCoeffs::default();

        one_pole_coeffs.set_sample_rate(f32::INFINITY);
    }

    fn assert_coeffs_rust_c(rust_coeffs: OnePoleCoeffs, c_coeffs: bw_one_pole_coeffs) {
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
}
