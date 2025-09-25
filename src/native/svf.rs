//! **State variable filter** (2nd order, 12 dB/oct) model with separated lowpass,
//! bandpass, and highpass outputs.
//!
//! # Example
//! ```rust
//! use brickworks_rs::native::svf::*;
//!
//! const SAMPLE_RATE: f32 = 48_000.0;
//! const N_CHANNELS: usize = 2;
//!
//! let mut svf = SVF::new();
//! svf.set_sample_rate(SAMPLE_RATE);
//! let cutoff = 185.0;
//! let q = 0.707;
//! let prewarpfreq = 185.0;
//! let n_samples = 2;
//!
//! let x_ch0 = [1.0, 0.0];
//! let x_ch1 = [1.0, 0.0];
//! let x: [&[f32]; 2] = [&x_ch0, &x_ch1];
//!
//! let mut y_lp_ch0 = [0.0, 0.0];
//! let mut y_lp_ch1 = [0.0, 0.0];
//! let mut y_lp: [Option<&mut [f32]>; N_CHANNELS] =
//!     [Some(&mut y_lp_ch0), Some(&mut y_lp_ch1)];
//!
//! svf.set_sample_rate(SAMPLE_RATE);
//! svf.set_cutoff(cutoff);
//! svf.set_prewarp_at_cutoff(true);
//! svf.set_q(q);
//! svf.set_prewarp_freq(prewarpfreq);
//! svf.reset(0.0, None, None, None);
//!
//!
//! svf.process(
//!     &x,
//!     Some(&mut y_lp),
//!     None,
//!     None,
//!     n_samples,
//! );
//! ```
//! # Notes
//! This module provides a native Rust implementation of the SVF Filter, but
//! the same interface is also available via bindings to the original C library
//! at [crate::c_wrapper::svf].
//! Original implementation by [Orastron](https://www.orastron.com/algorithms/bw_svf).
//!
use crate::native::{
    math::{minf, rcpf, tanf},
    one_pole::{OnePoleCoeffs, OnePoleState},
};
use std::f32::consts::PI;

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_positive, debug_assert_range};
/// State Variable Filter (SVF) with multiple channels.
///
/// The SVF provides lowpass, bandpass, and highpass outputs for each channel.
/// Filter parameters are shared across channels, while each channel maintains
/// its own state.
///
/// This struct manages both the filter coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
/// It wraps:
/// - [`SVFCoeffs`]
/// - [`SVFState`]
///
///
/// # Usage
/// ```rust
/// use brickworks_rs::native::svf::SVF;
/// const N_CHANNELS: usize = 2;
/// let mut svf = SVF::<N_CHANNELS>::new();
/// svf.set_sample_rate(48_000.0);
/// svf.set_cutoff(1_000.0);
/// svf.set_q(0.707);
/// // process audio with svf.process(...)
/// ```
#[derive(Debug)]
pub struct SVF<const N_CHANNELS: usize> {
    coeffs: SVFCoeffs<N_CHANNELS>,
    states: [SVFState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SVF<N_CHANNELS> {
    /// Creates a new `SVF` filter with default parameters and zeroed state.
    pub fn new() -> Self {
        Self {
            coeffs: SVFCoeffs::new(),
            states: [SVFState::new(); N_CHANNELS],
        }
    }
    /// Sets the filter's sample rate (Hz).
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }
    /// Resets the given state to its initial values using the given coeffs and the
    /// initial input value x_0.
    ///
    /// The corresponding initial lowpass, bandpass, and highpass output values are
    /// put into y_lp_0, y_bp_0, and y_hp_0 respectively if they are Some.
    pub fn reset(
        &mut self,
        x0: f32,
        y_lp0: Option<&mut [f32; N_CHANNELS]>,
        y_bp0: Option<&mut [f32; N_CHANNELS]>,
        y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        self.coeffs.reset_coeffs();
        match (y_lp0, y_bp0, y_hp0) {
            (Some(lp0), Some(bp0), Some(hp0)) => {
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut lp0[channel],
                        &mut bp0[channel],
                        &mut hp0[channel],
                    );
                });
            }
            (Some(lp0), Some(bp0), None) => {
                let mut v_hp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut lp0[channel],
                        &mut bp0[channel],
                        &mut v_hp,
                    );
                });
            }
            (Some(lp0), None, Some(hp0)) => {
                let mut v_bp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut lp0[channel],
                        &mut v_bp,
                        &mut hp0[channel],
                    );
                });
            }
            (Some(lp0), None, None) => {
                let (mut v_bp, mut v_hp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut lp0[channel],
                        &mut v_bp,
                        &mut v_hp,
                    );
                });
            }
            (None, Some(bp0), Some(hp0)) => {
                let mut v_lp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut v_lp,
                        &mut bp0[channel],
                        &mut hp0[channel],
                    );
                });
            }
            (None, Some(bp0), None) => {
                let (mut v_lp, mut v_hp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut v_lp,
                        &mut bp0[channel],
                        &mut v_hp,
                    );
                });
            }
            (None, None, Some(hp0)) => {
                let (mut v_lp, mut v_bp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut v_lp,
                        &mut v_bp,
                        &mut hp0[channel],
                    );
                });
            }
            (None, None, None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.coeffs.reset_state(
                        &mut self.states[channel],
                        x0,
                        &mut v_lp,
                        &mut v_bp,
                        &mut v_hp,
                    );
                });
            }
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
        self.coeffs.reset_coeffs();
        self.coeffs
            .reset_state_multi(&mut self.states, x0, y_lp0, y_bp0, y_hp0);
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
        self.coeffs
            .process_multi(&mut self.states, x, y_lp, y_bp, y_hp, n_samples);
    }
    /// Sets the cutoff frequency to the given value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    pub fn set_cutoff(&mut self, value: f32) {
        self.coeffs.set_cutoff(value);
    }
    /// Sets the quality factor to the given value.
    ///
    /// Valid range: [1e-6, 1e6].
    ///
    /// Default value: 0.5.
    pub fn set_q(&mut self, value: f32) {
        self.coeffs.set_q(value);
    }
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (non-0) or not (0).
    ///
    /// Default value: non-0 (on).
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.coeffs.set_prewarp_at_cutoff(value);
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
        self.coeffs.set_prewarp_freq(value);
    }
}

impl<const N_CHANNELS: usize> Default for SVF<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Coefficients and related.
#[derive(Debug)]
pub struct SVFCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs<N_CHANNELS>,
    smooth_cutoff_state: OnePoleState,
    smooth_q_state: OnePoleState,
    smooth_prewarp_freq_state: OnePoleState,

    // Coefficients
    t_k: f32,
    prewarp_freq_max: f32,

    kf: f32,
    kbl: f32,
    k: f32,
    hp_hb: f32,
    hp_x: f32,

    // Parameters
    cutoff: f32,
    q: f32,
    prewarp_k: f32,
    prewarp_freq: f32,
}

impl<const N_CHANNELS: usize> SVFCoeffs<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = Self {
            smooth_coeffs: OnePoleCoeffs::new(),
            smooth_cutoff_state: OnePoleState::new(),
            smooth_q_state: OnePoleState::new(),
            smooth_prewarp_freq_state: OnePoleState::new(),
            t_k: 0.0,
            prewarp_freq_max: 0.0,
            kf: 0.0,
            kbl: 0.0,
            k: 0.0,
            hp_hb: 0.0,
            hp_x: 0.0,
            cutoff: 1e3,
            q: 0.5,
            prewarp_k: 1.0,
            prewarp_freq: 1e3,
        };
        coeffs.smooth_coeffs.set_tau(0.005);
        coeffs.smooth_coeffs.set_sticky_thresh(1e-3);

        coeffs
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(sample_rate.is_finite());
            debug_assert_positive(sample_rate);
        }
        self.smooth_coeffs.set_sample_rate(sample_rate);
        self.smooth_coeffs.reset_coeffs();

        self.t_k = PI / sample_rate;
        self.prewarp_freq_max = 0.499 * sample_rate;
    }
    /// Resets coefficients to assume their target values.
    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.smooth_coeffs
            .reset_state(&mut self.smooth_cutoff_state, self.cutoff);
        self.smooth_coeffs
            .reset_state(&mut self.smooth_q_state, self.q);
        self.smooth_coeffs.reset_state(
            &mut self.smooth_prewarp_freq_state,
            self.prewarp_freq + self.prewarp_k * (self.cutoff - self.prewarp_freq),
        );
        self.do_update_coeffs(true);
    }
    /// Resets the given state to its initial values using the initial input `x0`.
    ///
    /// # Notes
    /// The corresponding initial lowpass, bandpass, and highpass outputs are
    /// stored in `y_lp`, `y_bp`, and `y_hp`.
    #[inline(always)]
    pub fn reset_state(
        &mut self,
        state: &mut SVFState,
        x0: f32,
        y_lp0: &mut f32,
        y_bp0: &mut f32,
        y_hp0: &mut f32,
    ) {
        state.hp_z1 = 0.0;
        state.lp_z1 = x0;
        state.bp_z1 = 0.0;
        state.cutoff_z1 = self.cutoff;
        *y_lp0 = x0;
        *y_bp0 = 0.0;
        *y_hp0 = 0.0;
    }
    /// Resets each of the `N_CHANNELS` states to their initial values using the
    /// provided coefficients and initial inputs.
    ///
    /// # Notes
    /// The corresponding initial lowpass, bandpass, and highpass outputs are
    /// written into `y_lp`, `y_bp`, and `y_hp` arrays if provided.
    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        state: &mut [SVFState; N_CHANNELS],
        x0: &[f32; N_CHANNELS],
        y_lp0: Option<&mut [f32; N_CHANNELS]>,
        y_bp0: Option<&mut [f32; N_CHANNELS]>,
        y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match (y_lp0, y_bp0, y_hp0) {
            (Some(lp0), Some(bp0), Some(hp0)) => {
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut lp0[channel],
                        &mut bp0[channel],
                        &mut hp0[channel],
                    );
                });
            }
            (Some(lp0), Some(bp0), None) => {
                let mut v_hp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut lp0[channel],
                        &mut bp0[channel],
                        &mut v_hp,
                    );
                });
            }
            (Some(lp0), None, Some(hp0)) => {
                let mut v_bp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut lp0[channel],
                        &mut v_bp,
                        &mut hp0[channel],
                    );
                });
            }
            (Some(lp0), None, None) => {
                let (mut v_bp, mut v_hp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut lp0[channel],
                        &mut v_bp,
                        &mut v_hp,
                    );
                });
            }
            (None, Some(bp0), Some(hp0)) => {
                let mut v_lp = 0.0;
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut v_lp,
                        &mut bp0[channel],
                        &mut hp0[channel],
                    );
                });
            }
            (None, Some(bp0), None) => {
                let (mut v_lp, mut v_hp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut v_lp,
                        &mut bp0[channel],
                        &mut v_hp,
                    );
                });
            }
            (None, None, Some(hp0)) => {
                let (mut v_lp, mut v_bp) = (0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut v_lp,
                        &mut v_bp,
                        &mut hp0[channel],
                    );
                });
            }
            (None, None, None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..N_CHANNELS).for_each(|channel| {
                    self.reset_state(
                        &mut state[channel],
                        x0[channel],
                        &mut v_lp,
                        &mut v_bp,
                        &mut v_hp,
                    );
                });
            }
        }
    }

    // Not implemented yet: C version only contained assertions
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Triggers audio-rate update of coefficients in `OnePoleCoeffs`.
    // #[inline(always)]
    // pub fn update_coeffs_ctrl(&mut self) {
    //     todo!()
    // }
    /// Triggers audio-rate update of coefficients.
    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.do_update_coeffs(false);
    }
    /// Processes a single input sample `x`, updating the provided state.
    ///
    /// # Notes
    /// The lowpass, bandpass, and highpass outputs are written into
    /// `y_lp`, `y_bp`, and `y_hp`, respectively.
    #[inline(always)]
    pub fn process1(
        &mut self,
        state: &mut SVFState,
        x: f32,
        y_lp: &mut f32,
        y_bp: &mut f32,
        y_hp: &mut f32,
    ) {
        let kk = self.kf * state.cutoff_z1;
        let lp_xz1 = state.lp_z1 + kk * state.bp_z1;
        let bp_xz1 = state.bp_z1 + kk * state.hp_z1;
        *y_hp = self.hp_x * (x - self.hp_hb * bp_xz1 - lp_xz1);
        *y_bp = bp_xz1 + self.kbl * *y_hp;
        *y_lp = lp_xz1 + self.kbl * *y_bp;
        state.hp_z1 = *y_hp;
        state.lp_z1 = *y_lp;
        state.bp_z1 = *y_bp;
        state.cutoff_z1 = self.smooth_cutoff_state.get_y_z1();

        debug_assert!(y_lp.is_finite());
        debug_assert!(y_lp.is_finite());
        debug_assert!(y_bp.is_finite());
    }
    /// Processes the first `n_samples` of the input buffer `x`.
    ///
    /// # Notes
    /// The coefficients and state are updated during processing.
    #[inline(always)]
    pub fn process(
        &mut self,
        state: &mut SVFState,
        x: &[f32],
        y_lp: Option<&mut [f32]>,
        y_bp: Option<&mut [f32]>,
        y_hp: Option<&mut [f32]>,
        n_samples: usize,
    ) {
        debug_assert!(
            y_lp.as_ref()
                .zip(y_bp.as_ref())
                .is_none_or(|(lp, bp)| lp != bp)
        );
        debug_assert!(
            y_lp.as_ref()
                .zip(y_hp.as_ref())
                .is_none_or(|(lp, hp)| lp != hp)
        );
        debug_assert!(
            y_bp.as_ref()
                .zip(y_hp.as_ref())
                .is_none_or(|(bp, hp)| bp != hp)
        );

        match (y_lp, y_bp, y_hp) {
            (Some(lp), Some(bp), Some(hp)) => {
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(
                        state,
                        x[sample],
                        &mut lp[sample],
                        &mut bp[sample],
                        &mut hp[sample],
                    );
                });
            }
            (Some(lp), Some(bp), None) => {
                let mut v_hp = 0.0;
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(
                        state,
                        x[sample],
                        &mut lp[sample],
                        &mut bp[sample],
                        &mut v_hp,
                    );
                });
            }
            (Some(lp), None, Some(hp)) => {
                let mut v_bp = 0.0;
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(
                        state,
                        x[sample],
                        &mut lp[sample],
                        &mut v_bp,
                        &mut hp[sample],
                    );
                });
            }
            (Some(lp), None, None) => {
                let (mut v_bp, mut v_hp) = (0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(state, x[sample], &mut lp[sample], &mut v_bp, &mut v_hp);
                });
            }
            (None, Some(bp), Some(hp)) => {
                let mut v_lp = 0.0;
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(
                        state,
                        x[sample],
                        &mut v_lp,
                        &mut bp[sample],
                        &mut hp[sample],
                    );
                });
            }
            (None, Some(bp), None) => {
                let (mut v_lp, mut v_hp) = (0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(state, x[sample], &mut v_lp, &mut bp[sample], &mut v_hp);
                });
            }
            (None, None, Some(hp)) => {
                let (mut v_lp, mut v_bp) = (0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(state, x[sample], &mut v_lp, &mut v_bp, &mut hp[sample]);
                });
            }
            (None, None, None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    self.process1(state, x[sample], &mut v_lp, &mut v_bp, &mut v_hp);
                });
            }
        }
    }
    /// Processes the first `n_samples` of the `N_CHANNELS` input buffers `x`.
    ///
    /// # Notes
    /// The coefficients and all channel states are updated during processing.
    /// Output buffers must be provided for all channels.
    #[inline(always)]
    pub fn process_multi(
        &mut self,
        state: &mut [SVFState; N_CHANNELS],
        x: &[&[f32]],
        y_lp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_bp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_hp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        n_samples: usize,
    ) {
        // missing debug
        match (y_lp, y_bp, y_hp) {
            (Some(lp), Some(bp), Some(hp)) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(lp) = &mut lp[channel] {
                            lp[sample] = v_lp
                        }
                        if let Some(bp) = &mut bp[channel] {
                            bp[sample] = v_bp
                        }
                        if let Some(hp) = &mut hp[channel] {
                            hp[sample] = v_hp
                        }
                    })
                });
            }
            (Some(lp), Some(bp), None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(lp) = &mut lp[channel] {
                            lp[sample] = v_lp
                        }
                        if let Some(bp) = &mut bp[channel] {
                            bp[sample] = v_bp
                        }
                    })
                });
            }
            (Some(lp), None, Some(hp)) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(lp) = &mut lp[channel] {
                            lp[sample] = v_lp
                        }
                        if let Some(hp) = &mut hp[channel] {
                            hp[sample] = v_hp
                        }
                    })
                });
            }
            (Some(lp), None, None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(lp) = &mut lp[channel] {
                            lp[sample] = v_lp
                        }
                    })
                });
            }
            (None, Some(bp), Some(hp)) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(bp) = &mut bp[channel] {
                            bp[sample] = v_bp
                        }
                        if let Some(hp) = &mut hp[channel] {
                            hp[sample] = v_hp
                        }
                    })
                });
            }
            (None, Some(bp), None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(bp) = &mut bp[channel] {
                            bp[sample] = v_bp
                        }
                    })
                });
            }
            (None, None, Some(hp)) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                        if let Some(hp) = &mut hp[channel] {
                            hp[sample] = v_hp
                        }
                    })
                });
            }
            (None, None, None) => {
                let (mut v_lp, mut v_bp, mut v_hp) = (0.0, 0.0, 0.0);
                (0..n_samples).for_each(|sample| {
                    self.update_coeffs_audio();
                    (0..N_CHANNELS).for_each(|channel| {
                        self.process1(
                            &mut state[channel],
                            x[channel][sample],
                            &mut v_lp,
                            &mut v_bp,
                            &mut v_hp,
                        );
                    })
                });
            }
        }
    }
    /// Sets the cutoff frequency to the given value (Hz).
    ///
    /// Valid range: [1e-6, 1e12].
    ///
    /// Default value: 1e3.
    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(value.is_finite());
            debug_assert_range(1e-6..=1e12, value);
        }

        self.cutoff = value;
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
            debug_assert!(value.is_finite());
            debug_assert_range(1e-6..=1e6, value);
        }

        self.q = value;
    }
    /// Sets whether bilinear transform prewarping frequency should match the cutoff
    /// frequency (true) or not (false).
    ///
    /// Default value: true (on).
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.prewarp_k = if value { 1.0 } else { 0.0 }
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
            debug_assert!(value.is_finite());
            debug_assert_range(1e-6..=1e12, value);
        }

        self.prewarp_freq = value;
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
    #[inline(always)]
    fn do_update_coeffs(&mut self, force: bool) {
        let prewarp_freq = self.prewarp_freq + self.prewarp_k * (self.cutoff - self.prewarp_freq);
        let mut cutoff_cur = self.smooth_cutoff_state.get_y_z1();
        let mut prewarp_freq_cur = self.smooth_prewarp_freq_state.get_y_z1();
        let mut q_cur = self.smooth_q_state.get_y_z1();
        let cutoff_changed = force || self.cutoff != cutoff_cur;
        let prewarp_freq_changed = force || prewarp_freq != prewarp_freq_cur;
        let q_changed = force || self.q != q_cur;

        if cutoff_changed || prewarp_freq_changed || q_changed {
            if cutoff_changed || prewarp_freq_changed {
                if cutoff_changed {
                    cutoff_cur = self
                        .smooth_coeffs
                        .process1_sticky_rel(&mut self.smooth_cutoff_state, self.cutoff);
                }

                if prewarp_freq_changed {
                    prewarp_freq_cur = self
                        .smooth_coeffs
                        .process1_sticky_rel(&mut self.smooth_prewarp_freq_state, prewarp_freq);
                    let f = minf(prewarp_freq_cur, self.prewarp_freq_max);
                    self.kf = tanf(self.t_k * f) * rcpf(f);
                }

                self.kbl = self.kf * cutoff_cur;
            }

            if q_changed {
                q_cur = self
                    .smooth_coeffs
                    .process1_sticky_abs(&mut self.smooth_q_state, self.q);
                self.k = rcpf(q_cur);
            }

            self.hp_hb = self.k + self.kbl;
            self.hp_x = rcpf(1.0 + self.kbl * self.hp_hb);
        }
    }
}

impl<const N_CHANNELS: usize> Default for SVFCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
/// Internal state and related.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SVFState {
    // State
    hp_z1: f32,
    lp_z1: f32,
    bp_z1: f32,
    cutoff_z1: f32,
}

impl SVFState {
    pub(crate) fn new() -> Self {
        Self {
            hp_z1: 0.0,
            lp_z1: 0.0,
            bp_z1: 0.0,
            cutoff_z1: 0.0,
        }
    }
}

impl Default for SVFState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::SVF;
    use crate::{
        c_wrapper::{bw_svf_coeffs, bw_svf_state, svf::SVF as SVFWrapper},
        native::{
            one_pole::tests::{assert_one_pole_coeffs, assert_one_pole_state},
            svf::{SVFCoeffs, SVFState},
        },
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    type SVFT = SVF<N_CHANNELS>;
    type SVFWrapperT = SVFWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_svf = SVFT::new();
        let c_svf = SVFWrapperT::new();

        assert_svf(&rust_svf, &c_svf);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);

        assert_svf(&rust_svf, &c_svf);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -48000")]
    fn set_sample_rate_must_be_positive() {
        let mut rust_svf = SVFT::new();
        rust_svf.set_sample_rate(-SAMPLE_RATE);
    }

    #[test]
    fn reset_coeffs() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);

        c_svf.coeffs.reset_coeffs();
        rust_svf.coeffs.reset_coeffs();

        assert_svf(&rust_svf, &c_svf);
    }

    #[test]
    fn reset_none() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);

        let x0 = 1.0;
        rust_svf.reset(x0, None, None, None);
        c_svf.reset(x0, None, None, None);

        assert_svf(&rust_svf, &c_svf);
    }

    #[test]
    fn reset_some() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);

        let x0 = 1.0;
        let mut y_lp0 = [1.0, 1.1];
        let mut y_bp0 = [2.0, 2.1];
        let mut y_hp0 = [3.0, 3.1];

        let mut c_y_lp0 = [1.0, 1.1];
        let mut c_y_bp0 = [2.0, 2.1];
        let mut c_y_hp0 = [3.0, 3.1];

        rust_svf.reset(x0, Some(&mut y_lp0), Some(&mut y_bp0), Some(&mut y_hp0));
        c_svf.reset(
            x0,
            Some(&mut c_y_lp0),
            Some(&mut c_y_bp0),
            Some(&mut c_y_hp0),
        );

        assert_svf(&rust_svf, &c_svf);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(y_lp0[channel], c_y_lp0[channel]);
            assert_eq!(y_bp0[channel], c_y_bp0[channel]);
            assert_eq!(y_hp0[channel], c_y_hp0[channel]);
        });
    }

    #[test]
    fn process1() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        let cutoff = 739.99;
        let q = 1.0;
        let prewarpfreq = 740.0;

        let x = 1.0;
        let mut y_lp = [1.0, 1.1];
        let mut y_bp = [2.0, 2.1];
        let mut y_hp = [3.0, 3.1];

        let mut c_y_lp = [1.0, 1.1];
        let mut c_y_bp = [2.0, 2.1];
        let mut c_y_hp = [3.0, 3.1];

        rust_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_cutoff(cutoff);
        rust_svf.set_prewarp_at_cutoff(true);
        rust_svf.set_q(q);
        rust_svf.set_prewarp_freq(prewarpfreq);
        rust_svf.reset(0.0, None, None, None);

        c_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_cutoff(cutoff);
        c_svf.set_prewarp_at_cutoff(true);
        c_svf.set_q(q);
        c_svf.set_prewarp_freq(prewarpfreq);
        c_svf.reset(0.0, None, None, None);

        (0..N_CHANNELS).for_each(|channel| {
            rust_svf.coeffs.process1(
                &mut rust_svf.states[channel],
                x,
                &mut y_lp[channel],
                &mut y_bp[channel],
                &mut y_hp[channel],
            );
            c_svf.coeffs.process1(
                &mut c_svf.states[channel],
                x,
                &mut c_y_lp[channel],
                &mut c_y_bp[channel],
                &mut c_y_hp[channel],
            );
        });

        assert_svf(&rust_svf, &c_svf);

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(y_lp[channel], c_y_lp[channel]);
            assert_eq!(y_bp[channel], c_y_bp[channel]);
            assert_eq!(y_hp[channel], c_y_hp[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);
        let cutoff = 185.0;
        let q = 0.707;
        let prewarpfreq = 185.0;
        let n_samples = 2;

        let x_ch0 = [1.0, 0.0];
        let x_ch1 = [1.0, 0.0];
        let x: [&[f32]; 2] = [&x_ch0, &x_ch1];

        let mut y_lp_ch0 = [0.0, 0.0];
        let mut y_lp_ch1 = [0.0, 0.0];
        let mut y_lp: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut y_lp_ch0), Some(&mut y_lp_ch1)];
        let mut y_bp_ch0 = [0.0, 0.0];
        let mut y_bp_ch1 = [0.0, 0.0];
        let mut y_bp: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut y_bp_ch0), Some(&mut y_bp_ch1)];
        let mut y_hp_ch0 = [0.0, 0.0];
        let mut y_hp_ch1 = [0.0, 0.0];
        let mut y_hp: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut y_hp_ch0), Some(&mut y_hp_ch1)];

        let mut c_y_lp_ch0 = [0.0, 0.0];
        let mut c_y_lp_ch1 = [0.0, 0.0];
        let mut c_y_lp: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut c_y_lp_ch0), Some(&mut c_y_lp_ch1)];
        let mut c_y_bp_ch0 = [0.0, 0.0];
        let mut c_y_bp_ch1 = [0.0, 0.0];
        let mut c_y_bp: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut c_y_bp_ch0), Some(&mut c_y_bp_ch1)];
        let mut c_y_hp_ch0 = [0.0, 0.0];
        let mut c_y_hp_ch1 = [0.0, 0.0];
        let mut c_y_hp: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut c_y_hp_ch0), Some(&mut c_y_hp_ch1)];

        rust_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_cutoff(cutoff);
        rust_svf.set_prewarp_at_cutoff(true);
        rust_svf.set_q(q);
        rust_svf.set_prewarp_freq(prewarpfreq);
        rust_svf.reset(0.0, None, None, None);

        c_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_cutoff(cutoff);
        c_svf.set_prewarp_at_cutoff(true);
        c_svf.set_q(q);
        c_svf.set_prewarp_freq(prewarpfreq);
        c_svf.reset(0.0, None, None, None);

        rust_svf.process(
            &x,
            Some(&mut y_lp),
            Some(&mut y_bp),
            Some(&mut y_hp),
            n_samples,
        );
        c_svf.process(
            &x,
            Some(&mut c_y_lp),
            Some(&mut c_y_bp),
            Some(&mut c_y_hp),
            n_samples,
        );

        assert_svf(&rust_svf, &c_svf);
        (0..N_CHANNELS).for_each(|channel| {
            (0..n_samples).for_each(|sample| {
                assert_eq!(
                    y_lp[channel].as_ref().unwrap()[sample],
                    c_y_lp[channel].as_ref().unwrap()[sample]
                );
                assert_eq!(
                    y_bp[channel].as_ref().unwrap()[sample],
                    c_y_bp[channel].as_ref().unwrap()[sample]
                );
                assert_eq!(
                    y_hp[channel].as_ref().unwrap()[sample],
                    c_y_hp[channel].as_ref().unwrap()[sample]
                );
            })
        });
    }

    #[test]
    fn set_cutoff_valid() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        let cutoff = 1479.98;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_cutoff(cutoff);
        c_svf.set_cutoff(cutoff);

        assert_svf(&rust_svf, &c_svf);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 1e-7")]
    fn set_cutoff_invalid() {
        let mut rust_svf = SVFT::new();
        let cutoff = 1e-7;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_cutoff(cutoff);
    }

    #[test]
    fn set_q() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        let q = 0.707;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_q(q);
        c_svf.set_q(q);

        assert_svf(&rust_svf, &c_svf);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [1e-6, 1e6], got 1e-7")]
    fn set_q_invalid() {
        let mut rust_svf = SVFT::new();
        let q = 1e-7;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_q(q);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);

        rust_svf.set_prewarp_at_cutoff(true);
        c_svf.set_prewarp_at_cutoff(true);

        assert_svf(&rust_svf, &c_svf);
    }

    #[test]
    fn set_prewarp_freq_valid() {
        let mut rust_svf = SVFT::new();
        let mut c_svf = SVFWrapperT::new();
        let prewarp_freq = 92.50;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        c_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_prewarp_freq(prewarp_freq);
        c_svf.set_prewarp_freq(prewarp_freq);

        assert_svf(&rust_svf, &c_svf);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 1e13")]
    fn set_prewarp_freq_invalid() {
        let mut rust_svf = SVFT::new();
        let prewarp_freq = 1e13;

        rust_svf.set_sample_rate(SAMPLE_RATE);
        rust_svf.set_prewarp_freq(prewarp_freq);
    }

    fn assert_svf<const N_CHANNELS: usize>(
        rust_svf: &SVF<N_CHANNELS>,
        c_svf: &SVFWrapper<N_CHANNELS>,
    ) {
        assert_svf_coeffs(&rust_svf.coeffs, &c_svf.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_svf_states(&rust_svf.states[channel], &c_svf.states[channel]);
        });
    }

    pub(crate) fn assert_svf_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &SVFCoeffs<N_CHANNELS>,
        c_coeffs: &bw_svf_coeffs,
    ) {
        let pre_message = "svf.coeff.";
        let post_message = "does not match";
        assert_one_pole_coeffs(&rust_coeffs.smooth_coeffs, &c_coeffs.smooth_coeffs);
        assert_one_pole_state(
            &rust_coeffs.smooth_cutoff_state,
            &c_coeffs.smooth_cutoff_state,
        );
        assert_one_pole_state(&rust_coeffs.smooth_q_state, &c_coeffs.smooth_Q_state);
        assert_one_pole_state(
            &rust_coeffs.smooth_prewarp_freq_state,
            &c_coeffs.smooth_prewarp_freq_state,
        );
        assert_eq!(
            rust_coeffs.t_k, c_coeffs.t_k,
            "{}t_k {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.prewarp_freq_max, c_coeffs.prewarp_freq_max,
            "{}prewarp_freq_max {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.kf, c_coeffs.kf,
            "{}kf {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.kbl, c_coeffs.kbl,
            "{}kbl {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.k, c_coeffs.k,
            "{}k {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.hp_hb, c_coeffs.hp_hb,
            "{}hp_hb {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.hp_x, c_coeffs.hp_x,
            "{}hp_x {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.cutoff, c_coeffs.cutoff,
            "{}cutoff {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.q, c_coeffs.Q,
            "{}Q {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.prewarp_k, c_coeffs.prewarp_k,
            "{}prewarp_k {}",
            pre_message, post_message,
        );
        assert_eq!(
            rust_coeffs.prewarp_freq, c_coeffs.prewarp_freq,
            "{}prewarp_freq {}",
            pre_message, post_message,
        );
    }

    pub fn assert_svf_states(rust_state: &SVFState, c_state: &bw_svf_state) {
        assert_eq!(rust_state.hp_z1, c_state.hp_z1);
        assert_eq!(rust_state.lp_z1, c_state.lp_z1);
        assert_eq!(rust_state.bp_z1, c_state.bp_z1);
        assert_eq!(rust_state.cutoff_z1, c_state.cutoff_z1);
    }
}
