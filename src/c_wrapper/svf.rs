//! State variable filter (2nd order, 12 dB/oct) model with separated lowpass,
//! bandpass, and highpass outputs.
//!
//! # Example
//! ```rust
//! use brickworks_rs::c_wrapper::svf::*;
//!
//! const SAMPLE_RATE: f32 = 48_000.0;
//! const N_CHANNELS: usize = 2;
//!
//! fn main() {
//!     let mut svf = SVF::new();
//!     svf.set_sample_rate(SAMPLE_RATE);
//!     let cutoff = 185.0;
//!     let q = 0.707;
//!     let prewarpfreq = 185.0;
//!     let n_samples = 2;
//!
//!     let x_ch0 = [1.0, 0.0];
//!     let x_ch1 = [1.0, 0.0];
//!     let x: [&[f32]; 2] = [&x_ch0, &x_ch1];
//!
//!     let mut y_lp_ch0 = [0.0, 0.0];
//!     let mut y_lp_ch1 = [0.0, 0.0];
//!     let mut y_lp: [Option<&mut [f32]>; N_CHANNELS] =
//!         [Some(&mut y_lp_ch0), Some(&mut y_lp_ch1)];
//!
//!     svf.set_sample_rate(SAMPLE_RATE);
//!     svf.set_cutoff(cutoff);
//!     svf.set_prewarp_at_cutoff(true);
//!     svf.set_q(q);
//!     svf.set_prewarp_freq(prewarpfreq);
//!     svf.reset(0.0, None, None, None);
//!
//!
//!     svf.process(
//!         &x,
//!         Some(&mut y_lp),
//!         None,
//!         None,
//!         n_samples,
//!     );
//! }
//! ```
//! # Notes
//! This module provides Rust bindings to the original C implementation.
//! For a fully native Rust implementation with the same interface,
//! see [crate::native::svf].
//! Original C library by [Orastron](https://www.orastron.com/algorithms/bw_svf).
//!
use super::*;
use crate::c_wrapper::utils::{from_opt_to_raw, make_array};
use std::ptr::null_mut;
#[derive(Debug)]
/// State Variable Filter (SVF) with multiple channels.
///
/// The SVF provides lowpass, bandpass, and highpass outputs for each channel.
/// Filter parameters are shared across channels, while each channel maintains
/// its own state.
/// 
/// # Usage
/// ```rust 
/// use brickworks_rs::c_wrapper::svf::SVF;
/// const N_CHANNELS: usize = 2;
/// let mut svf = SVF::<N_CHANNELS>::new();
/// svf.set_sample_rate(48_000.0);
/// svf.set_cutoff(1_000.0);
/// svf.set_q(0.707);
/// // process audio with svf.process(...)
/// ```
pub struct SVF<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_svf_coeffs,
    pub(crate) states: [bw_svf_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SVF<N_CHANNELS> {
    /// Creates a new `SVF` filter with default parameters and zeroed state.
    pub fn new() -> Self {
        let mut coeffs = bw_svf_coeffs::default();
        unsafe {
            bw_svf_init(&mut coeffs);
        }
        Self {
            coeffs,
            states: make_array::<bw_svf_state, N_CHANNELS>(),
        }
    }
    /// Sets the filter's sample rate (Hz).
    pub fn set_sample_rate(&mut self, value: f32) {
        unsafe {
            bw_svf_set_sample_rate(&mut self.coeffs, value);
        }
    }
    /// Resets the given state to its initial values using the given coeffs and the
    /// initial input value x_0.
    ///
    /// The corresponding initial lowpass, bandpass, and highpass output values are
    /// put into y_lp_0, y_bp_0, and y_hp_0 respectively if they are Some.
    pub fn reset(
        &mut self,
        x0: f32,
        mut y_lp0: Option<&mut [f32; N_CHANNELS]>,
        mut y_bp0: Option<&mut [f32; N_CHANNELS]>,
        mut y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        unsafe { bw_svf_reset_coeffs(&mut self.coeffs) }
        // This is more readable, but it's actually slower than the c version
        // need to refactor as in native::svf::SVFCoeffs::reset_state_multi
        for channel in 0..N_CHANNELS {
            let mut v_lp = 0.0;
            let mut v_bp = 0.0;
            let mut v_hp = 0.0;

            let y_lp = match &mut y_lp0 {
                Some(array) => &mut array[channel],
                None => &mut v_lp,
            };
            let y_bp = match &mut y_bp0 {
                Some(array) => &mut array[channel],
                None => &mut v_bp,
            };
            let y_hp = match &mut y_hp0 {
                Some(array) => &mut array[channel],
                None => &mut v_hp,
            };

            unsafe {
                bw_svf_reset_state(
                    &mut self.coeffs,
                    &mut self.states[channel],
                    x0,
                    y_lp,
                    y_bp,
                    y_hp,
                )
            };
        }
    }
    /// Resets each of the n_channels states to its initial values using the given
    /// coeffs and the corresponding initial input value in the x_0 array.
    ///
    /// The corresponding initial lowpass, bandpass, and highpass output values are
    /// put into the y_lp0, y_bp0, and y_hp0 arrays, respectively, if they are Some.
    pub fn reset_multi(
        &mut self,
        x0: &[f32; N_CHANNELS],
        y_lp0: Option<&mut [f32; N_CHANNELS]>,
        y_bp0: Option<&mut [f32; N_CHANNELS]>,
        y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        unsafe {
            bw_svf_reset_coeffs(&mut self.coeffs);
        }

        let state_ptrs: [*mut bw_svf_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);

        let y_lp_p = match y_lp0 {
            Some(array) => array.as_mut_ptr(),
            None => null_mut(),
        };
        let y_bp_p = match y_bp0 {
            Some(array) => array.as_mut_ptr(),
            None => null_mut(),
        };
        let y_hp_p = match y_hp0 {
            Some(array) => array.as_mut_ptr(),
            None => null_mut(),
        };
        unsafe {
            bw_svf_reset_state_multi(
                &mut self.coeffs,
                state_ptrs.as_ptr(),
                x0.as_ptr(),
                y_lp_p,
                y_bp_p,
                y_hp_p,
                N_CHANNELS,
            );
        }
    }
    /// Processes the first `n_sample` of the `N_CHANNELS` input buffers `x` and fills
    /// the first `n_sample` of the `N_CHANNELS` output buffers `y_lp` (lowpass),
    /// `y_bp` (bandpass) and `y_hp`(highpass) if they are Some, while using and
    /// updating both the common coeffs and each of the `N_CHANNELS` states (control
    /// and audio rate).
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y_lp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_bp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_hp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        n_samples: usize,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_lp_ptrs: [*mut f32; N_CHANNELS] = from_opt_to_raw(y_lp);
        let mut y_bp_ptrs: [*mut f32; N_CHANNELS] = from_opt_to_raw(y_bp);
        let mut y_hp_ptrs: [*mut f32; N_CHANNELS] = from_opt_to_raw(y_hp);

        let mut states_ptrs: [*mut bw_svf_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut bw_svf_state);

        unsafe {
            bw_svf_process_multi(
                &mut self.coeffs,
                states_ptrs.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_lp_ptrs.as_mut_ptr(),
                y_bp_ptrs.as_mut_ptr(),
                y_hp_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }
    /// Sets the cutoff frequency to the given value (Hz) in SVF.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    pub fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_svf_set_cutoff(&mut self.coeffs, value);
        }
    }
    /// Sets the quality factor to the given value in SVF.
    ///
    /// Valid range: [1e-6, 1e6].
    ///
    /// Default value: 0.5.
    pub fn set_q(&mut self, value: f32) {
        unsafe {
            bw_svf_set_Q(&mut self.coeffs, value);
        }
    }
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (non-0) or not (0).
    ///
    /// Default value: non-0 (on).
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        unsafe {
            bw_svf_set_prewarp_at_cutoff(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }
    /// Sets the prewarping frequency value (Hz) in coeffs.
    ///
    /// Only used when the prewarp_at_cutoff parameter is off and however internally
    /// limited to avoid instability.
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    pub fn set_prewarp_freq(&mut self, value: f32) {
        unsafe {
            bw_svf_set_prewarp_freq(&mut self.coeffs, value);
        }
    }
}

impl<const N_CHANNELS: usize> Default for SVF<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

impl bw_svf_coeffs {
    // Wrapping these to test them against the native ones
    #[cfg(test)]
    pub(crate) fn reset_coeffs(&mut self) {
        unsafe {
            bw_svf_reset_coeffs(self);
        }
    }

    #[cfg(test)]
    pub(crate) fn process1(
        &mut self,
        state: &mut bw_svf_state,
        x: f32,
        y_lp: &mut f32,
        y_bp: &mut f32,
        y_hp: &mut f32,
    ) {
        unsafe {
            bw_svf_process1(self, state, x, y_lp, y_bp, y_hp);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::c_wrapper::{
        bw_rcpf, bw_svf_coeffs, bw_svf_init, bw_svf_process1, bw_svf_reset_coeffs,
        bw_svf_reset_state_multi, bw_svf_set_cutoff, bw_svf_set_prewarp_at_cutoff,
        bw_svf_set_prewarp_freq, bw_svf_set_sample_rate, bw_svf_state, bw_svf_update_coeffs_audio,
        svf::SVF,
    };
    use std::{f32::consts::PI, ptr::null_mut};

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    type SVFTest = SVF<N_CHANNELS>;

    #[test]
    fn new() {
        let svf = SVFTest::new();
        let tau_default = 0.005;
        let sticky_tresh_default = 1e-3;
        let cutoff;
        unsafe {
            cutoff = INVERSE_2_PI * bw_rcpf(tau_default);
        }

        assert_eq!(svf.coeffs.cutoff, 1e3);
        assert_eq!(svf.coeffs.Q, 0.5);
        assert_eq!(svf.coeffs.prewarp_freq, 1e3);
        assert_eq!(svf.coeffs.prewarp_k, 1.0);
        assert_eq!(svf.coeffs.smooth_coeffs.cutoff_up, cutoff);
        assert_eq!(svf.coeffs.smooth_coeffs.cutoff_down, cutoff);
        assert_eq!(svf.coeffs.smooth_coeffs.sticky_thresh, sticky_tresh_default);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut svf = SVFTest::new();
        svf.set_sample_rate(SAMPLE_RATE);

        assert_eq!(svf.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
        assert_eq!(svf.coeffs.t_k, PI / SAMPLE_RATE);
        assert_eq!(svf.coeffs.prewarp_freq_max, 0.499 * SAMPLE_RATE);
    }

    #[test]
    fn reset() {
        let mut svf = SVFTest::new();
        let x0 = 0.1;
        svf.set_sample_rate(SAMPLE_RATE);
        svf.reset(x0, None, None, None);

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(svf.states[channel].lp_z1, x0);
            assert_eq!(svf.states[channel].bp_z1, 0.);
            assert_eq!(svf.states[channel].hp_z1, 0.);
            assert_eq!(svf.states[channel].cutoff_z1, svf.coeffs.cutoff);
        });

        svf.set_cutoff(1000.0);
        let mut y_lp0: [f32; 2] = [10.0, 11.0];
        let mut y_bp0: [f32; 2] = [20.0, 22.0];
        let mut y_hp0: [f32; 2] = [30.0, 33.0];

        svf.reset(x0, Some(&mut y_lp0), Some(&mut y_bp0), Some(&mut y_hp0));
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(svf.states[channel].lp_z1, x0);
            assert_eq!(svf.states[channel].bp_z1, 0.);
            assert_eq!(svf.states[channel].hp_z1, 0.);
            assert_eq!(svf.states[channel].cutoff_z1, svf.coeffs.cutoff);
        });
    }

    #[test]
    fn reset_multi() {
        let mut svf = SVFTest::new();
        let x0: [f32; N_CHANNELS] = [0.1, 0.2];

        svf.set_sample_rate(SAMPLE_RATE);
        svf.reset_multi(&x0, None, None, None);

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(svf.states[channel].lp_z1, x0[channel]);
            assert_eq!(svf.states[channel].bp_z1, 0.);
            assert_eq!(svf.states[channel].hp_z1, 0.);
            assert_eq!(svf.states[channel].cutoff_z1, svf.coeffs.cutoff);
        });

        svf.set_cutoff(1000.0);
        let mut y_lp0: [f32; 2] = [10.0, 11.0];
        let mut y_bp0: [f32; 2] = [20.0, 22.0];
        let mut y_hp0: [f32; 2] = [30.0, 33.0];

        svf.reset_multi(&x0, Some(&mut y_lp0), Some(&mut y_bp0), Some(&mut y_hp0));
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(svf.states[channel].lp_z1, x0[channel]);
            assert_eq!(svf.states[channel].bp_z1, 0.);
            assert_eq!(svf.states[channel].hp_z1, 0.);
            assert_eq!(svf.states[channel].cutoff_z1, svf.coeffs.cutoff);
        });
    }

    #[test]
    fn process() {
        let mut svf = SVFTest::new();

        let sample_0: [f32; N_CHANNELS] = [6.0, 2.0];
        let sample_1: [f32; N_CHANNELS] = [6.0, 2.0];
        let x: [&[f32]; 2] = [&sample_0, &sample_1];

        // Wrapper data
        let mut y_lp_0: [f32; N_CHANNELS] = [1.0, 2.0];
        let mut y_lp_1: [f32; N_CHANNELS] = [3.0, 4.0];

        let mut y_bp_0: [f32; N_CHANNELS] = [5.0, 6.0];
        let mut y_bp_1: [f32; N_CHANNELS] = [7.0, 8.0];

        let mut y_hp_0: [f32; N_CHANNELS] = [9.0, 10.0];
        let mut y_hp_1: [f32; N_CHANNELS] = [11.0, 12.0];

        // C data
        let mut c_y_lp_0: [f32; N_CHANNELS] = y_lp_0.clone();
        let mut c_y_lp_1: [f32; N_CHANNELS] = y_lp_1.clone();

        let mut c_y_bp_0: [f32; N_CHANNELS] = y_bp_0.clone();
        let mut c_y_bp_1: [f32; N_CHANNELS] = y_bp_1.clone();

        let mut c_y_hp_0: [f32; N_CHANNELS] = y_hp_0.clone();
        let mut c_y_hp_1: [f32; N_CHANNELS] = y_hp_1.clone();

        // Prepare output for wrapper
        let mut y_lp: [Option<&mut [f32]>; 2] = [Some(&mut y_lp_0), Some(&mut y_lp_1)];
        let mut y_bp: [Option<&mut [f32]>; 2] = [Some(&mut y_bp_0), Some(&mut y_bp_1)];
        let mut y_hp: [Option<&mut [f32]>; 2] = [Some(&mut y_hp_0), Some(&mut y_hp_1)];

        // Prepare output for C
        let mut c_y_lp: [&mut [f32]; 2] = [&mut c_y_lp_0, &mut c_y_lp_1];
        let mut c_y_bp: [&mut [f32]; 2] = [&mut c_y_bp_0, &mut c_y_bp_1];
        let mut c_y_hp: [&mut [f32]; 2] = [&mut c_y_hp_0, &mut c_y_hp_1];

        // Parameters
        let cutoff = 1000.0;
        let prewarp_freq = 800.0;
        let n_samples = 2;
        // Wrapper
        svf.set_sample_rate(SAMPLE_RATE);
        svf.set_cutoff(cutoff);
        svf.set_prewarp_at_cutoff(true);
        svf.set_prewarp_freq(prewarp_freq);
        svf.reset_multi(&[0.0, 0.0], None, None, None);
        svf.process(
            &x,
            Some(&mut y_lp),
            Some(&mut y_bp),
            Some(&mut y_hp),
            n_samples,
        );

        // Process C
        let mut c_coeffs = bw_svf_coeffs::default();
        let mut c_state_ch0 = bw_svf_state::default();
        let mut c_state_ch1 = bw_svf_state::default();
        let c_states: [*mut bw_svf_state; N_CHANNELS] = [&mut c_state_ch0, &mut c_state_ch1];

        unsafe {
            bw_svf_init(&mut c_coeffs);
            bw_svf_set_sample_rate(&mut c_coeffs, SAMPLE_RATE);
            bw_svf_set_cutoff(&mut c_coeffs, cutoff);
            bw_svf_set_prewarp_freq(&mut c_coeffs, prewarp_freq);
            bw_svf_set_prewarp_at_cutoff(&mut c_coeffs, 1);
            bw_svf_reset_coeffs(&mut c_coeffs);
            bw_svf_reset_state_multi(
                &mut c_coeffs,
                c_states.as_ptr(),
                [0.0, 0.0].as_ptr(),
                null_mut(),
                null_mut(),
                null_mut(),
                N_CHANNELS,
            );

            (0..n_samples).for_each(|sample| {
                bw_svf_update_coeffs_audio(&mut c_coeffs);
                (0..N_CHANNELS).for_each(|channel| {
                    bw_svf_process1(
                        &c_coeffs,
                        c_states[channel],
                        x[channel][sample],
                        &mut c_y_lp[channel][sample],
                        &mut c_y_bp[channel][sample],
                        &mut c_y_hp[channel][sample],
                    );
                });
            });
        }

        (0..n_samples).for_each(|sample| {
            (0..N_CHANNELS).for_each(|channel| {
                assert_eq!(
                    y_lp[channel].as_ref().unwrap()[sample],
                    c_y_lp[channel][sample],
                    "low pass sample {sample} and channel {channel} does not match"
                );
                assert_eq!(
                    y_bp[channel].as_ref().unwrap()[sample],
                    c_y_bp[channel][sample],
                    "band pass sample {sample} and channel {channel} does not match"
                );
                assert_eq!(
                    y_hp[channel].as_ref().unwrap()[sample],
                    c_y_hp[channel][sample],
                    "high pass sample {sample} and channel {channel} does not match"
                );
            });
        });
    }

    #[test]
    fn set_cutoff_valid() {
        let mut svf = SVFTest::new();
        let cutoff = 124.3;

        svf.set_cutoff(cutoff);

        assert_eq!(svf.coeffs.cutoff, cutoff);
    }

    #[test]
    fn set_q_valid() {
        let mut svf = SVFTest::new();
        let q = 20.3;

        svf.set_q(q);

        assert_eq!(svf.coeffs.Q, q);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut svf = SVFTest::new();

        svf.set_prewarp_at_cutoff(true);

        assert_eq!(svf.coeffs.prewarp_k, 1.0);

        svf.set_prewarp_at_cutoff(false);

        assert_eq!(svf.coeffs.prewarp_k, 0.0);
    }

    #[test]
    fn set_prewarp_freq_valid() {
        let mut svf = SVFTest::new();
        let prewarp_freq = 100.9;

        svf.set_prewarp_freq(prewarp_freq);

        assert_eq!(svf.coeffs.prewarp_freq, prewarp_freq);
    }
}
