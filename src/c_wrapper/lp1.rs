//! **First-order lowpass filter** (6 dB/oct) with unitary DC gain.
//!
//! This is better suited to filtering actual audio than [crate::c_wrapper::one_pole].
//!
//! # Example
//! ```rust
//! use brickworks_rs::c_wrapper::lp1::LP1;
//!
//! const N_CHANNELS: usize = 2;
//! const N_SAMPLES: usize = 8;
//! let mut filter = LP1::<N_CHANNELS>::new();
//!
//! // Set sample rate
//! filter.set_sample_rate(44_100.0);
//!
//! // Set cutoff and prewarp
//! filter.set_cutoff(1_000.0);
//! filter.set_prewarp_at_cutoff(true);
//!
//! // Reset filter states
//! filter.reset(0.0, None);
//!
//!
//! // Example input: multi-channel array of samples
//! let input: [&[f32]; N_CHANNELS] = [&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!                                    &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
//!
//! let mut output_buffers: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];
//!
//! // Process the samples
//! filter.process(&input, &mut output_buffers, N_SAMPLES);
//!
//! println!("Filtered output: {:?}", output_buffers);
//! ```
//! # Notes
//! This module provides Rust bindings to the original C implementation.
//! For a fully native Rust implementation with the same interface,
//! see [crate::native::lp1].
//! Original C library by [Orastron](https://www.orastron.com/algorithms/bw_lp1).
use super::*;
use std::ptr::null_mut;
/// First-order lowpass filter (6 dB/oct) with unitary DC gain.
///
/// This struct manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`LP1Coeffs`]
/// - [`LP1State`]
///
/// # Usage
/// ```rust
/// use brickworks_rs::c_wrapper::lp1::LP1;
/// const N_CHANNELS: usize = 2;
/// let mut lp1 = LP1::<N_CHANNELS>::new();
/// lp1.set_sample_rate(48_000.0);
/// lp1.set_cutoff(1_000.0);
/// lp1.set_prewarp_at_cutoff(false);
/// lp1.set_prewarp_freq(990.0);
/// // process audio with lp1.process(...)
/// ```
pub struct LP1<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_lp1_coeffs,
    pub(crate) states: [bw_lp1_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> LP1<N_CHANNELS> {
    /// Creates a new instance with default parameters and zeroed state.
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_lp1_coeffs::default();
        unsafe {
            bw_lp1_init(&mut coeffs);
        }

        Self {
            coeffs,
            states: [bw_lp1_state::default(); N_CHANNELS],
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_lp1_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }
    /// Resets the states and coeffs for all channels to the initial input `x0`,
    /// or to 0 if `x0` is not provided.
    /// If `y0` is provided, the resulting initial outputs are stored in it.
    #[inline(always)]
    pub fn reset(&mut self, x0: f32, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_lp1_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] =
                        bw_lp1_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_lp1_reset_state(&mut self.coeffs, &mut self.states[channel], x0);
                }),
            }
        }
    }
    /// Resets the state and coefficients for all channels using the provided initial
    /// input values.
    ///
    /// Both the coefficients and all channel states are reset.
    /// If `y0` is `Some`, the initial outputs are written into it.
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_lp1_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_lp1_reset_coeffs(&mut self.coeffs);
            bw_lp1_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x` and writes the
    /// results to the first `n_samples` of the `N_CHANNELS` output buffers `y`,
    /// while updating the shared coefficients and each channel's state (control and audio rate).
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_lp1_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_lp1_process_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }
    /// Sets the cutoff frequency to the given value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_lp1_set_cutoff(&mut self.coeffs, value);
        }
    }
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        unsafe {
            bw_lp1_set_prewarp_at_cutoff(&mut self.coeffs, if value { 1 } else { 0 });
        }
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
        unsafe {
            bw_lp1_set_prewarp_freq(&mut self.coeffs, value);
        }
    }
}

impl bw_lp1_coeffs {
    #[cfg(test)]
    pub(crate) fn new() -> Self {
        let mut coeffs = Self::default();
        unsafe {
            bw_lp1_init(&mut coeffs);
        }
        coeffs
    }

    #[cfg(test)]
    pub(crate) fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_lp1_set_sample_rate(self, sample_rate);
        }
    }

    #[cfg(test)]
    pub(crate) fn process1(&mut self, state: &mut bw_lp1_state, x: f32) -> f32 {
        unsafe { bw_lp1_process1(self, state, x) }
    }

    #[cfg(test)]
    pub(crate) fn reset_coeffs(&mut self) {
        unsafe {
            bw_lp1_reset_coeffs(self);
        }
    }

    #[cfg(test)]
    pub(crate) fn reset_state(&mut self, state: &mut bw_lp1_state, x0: f32) -> f32 {
        unsafe { bw_lp1_reset_state(self, state, x0) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::math::INVERSE_2_PI;
    use std::f32::consts::PI;

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    const N_SAMPLES: usize = 8;

    type LP1T = LP1<N_CHANNELS>;

    #[test]
    fn new() {
        let lp1 = LP1T::new();

        let tau_default = 0.005;
        let sticky_tresh_default = 1e-3;
        let cutoff;
        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(lp1.coeffs.cutoff, 1e3);
        assert_eq!(lp1.coeffs.prewarp_k, 1.0);
        assert_eq!(lp1.coeffs.prewarp_freq, 1e3);
        assert_eq!(lp1.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(lp1.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(lp1.coeffs.smooth_coeffs.sticky_thresh, sticky_tresh_default);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_init
        );
    }

    #[test]
    fn set_sample_rate() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        assert_eq!(lp1.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
        assert_eq!(lp1.coeffs.t_k, PI / SAMPLE_RATE);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_set_sample_rate
        );
    }

    #[test]
    fn reset_none() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.0;
        lp1.reset(x0, None);

        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        )
    }

    #[test]
    fn reset_some() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = 4.0;
        let mut y0 = [0.0, 0.0];
        lp1.reset(x0, Some(&mut y0));

        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(y0[channel], x0);
            assert_eq!(lp1.states[channel].X_z1, 0.0);
            assert_eq!(lp1.states[channel].y_z1, x0);
        });
    }

    #[test]
    fn reset_multi() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.1, 0.2];
        let mut y0 = [0.0, 0.0];
        lp1.reset_multi(&x0, Some(&mut y0));

        assert_eq!(y0, x0);
        assert_eq!(
            lp1.coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(lp1.states[channel].X_z1, 0.0);
            assert_eq!(lp1.states[channel].y_z1, x0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];

        lp1.reset_multi(&x0, None);

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        lp1.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(lp1.coeffs.state >= bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs);
        assert_eq!(lp1.coeffs.reset_id, lp1.states[0].coeffs_reset_id);
    }

    #[test]
    fn set_cutoff() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let cutoff = 493.883;

        lp1.set_cutoff(cutoff);

        assert_eq!(lp1.coeffs.cutoff, cutoff);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        lp1.set_prewarp_at_cutoff(true);
        assert_eq!(lp1.coeffs.prewarp_k, 1.0);

        lp1.set_prewarp_at_cutoff(false);
        assert_eq!(lp1.coeffs.prewarp_k, 0.0);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut lp1 = LP1T::new();
        lp1.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 185.0;

        lp1.set_prewarp_freq(prewarp_freq);

        assert_eq!(lp1.coeffs.prewarp_freq, prewarp_freq);
    }
}
