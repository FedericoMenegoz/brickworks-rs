//! **Distortion Effect**
//!
//! Loosely inspired to the "rodent" distortion pedal.
//!
//! # Example
//! ```
//! use brickworks_rs::c_wrapper::dist::Dist;
//!
//! // Create a stereo (2-channel) distortion processor
//! let mut dist = Dist::<2>::new();
//!
//! // Configure processor
//! dist.set_sample_rate(48_000.0);
//! dist.set_distortion(0.8);
//! dist.set_tone(0.5);
//! dist.set_volume(0.9);
//!
//! // Prepare input and output buffers
//! let input_left: Vec<f32> = vec![0.0; 512];
//! let input_right: Vec<f32> = vec![0.0; 512];
//! let mut output_left: Vec<f32> = vec![0.0; 512];
//! let mut output_right: Vec<f32> = vec![0.0; 512];
//!
//! // Build channel slices
//! let inputs = [&input_left[..], &input_right[..]];
//! let mut outputs = [&mut output_left[..], &mut output_right[..]];
//!
//! // Reset state before processing
//! dist.reset(None, None);
//! // Process 512 samples
//! dist.process(&inputs, &mut outputs, 512);
//! ```
//! # Notes
//! This module provides Rust bindings to the original C implementation.
//! For a fully native Rust implementation with the same interface,
//! see [crate::native::dist].
//! Original C library by [Orastron](https://www.orastron.com/algorithms/bw_dist).
use crate::c_wrapper::{
    bw_dist_coeffs, bw_dist_init, bw_dist_process_multi, bw_dist_reset_coeffs, bw_dist_reset_state,
    bw_dist_reset_state_multi, bw_dist_set_distortion, bw_dist_set_sample_rate, bw_dist_set_tone,
    bw_dist_set_volume, bw_dist_state,
};
use std::ptr::null_mut;
/// Distortion effect.
///
/// Manages both the distortion coefficients and the runtime states
/// for a given number of channels (`N_CHANNELS`).  
///
/// # Usage
///
/// ```rust
/// use brickworks_rs::c_wrapper::dist::Dist;
///
/// // Stereo filter
/// const N_CHANNELS: usize = 2;
/// let mut dist = Dist::<N_CHANNELS>::new();
///
/// dist.set_sample_rate(44_100.0);
/// dist.set_distortion(0.3);
/// dist.set_tone(0.8);
/// dist.set_volume(0.5);
/// dist.reset(None, None);
/// // dist.process(...)
/// ```
pub struct Dist<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_dist_coeffs,
    pub(crate) states: [bw_dist_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Dist<N_CHANNELS> {
    /// Creates a new instance with all fields initialized.
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_dist_coeffs::default();

        unsafe {
            bw_dist_init(&mut coeffs);
            Self {
                coeffs,
                states: [bw_dist_state::default(); N_CHANNELS],
            }
        }
    }
    /// Sets the sample rate (Hz).
    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_dist_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }
    /// Resets the coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value x0.
    ///
    /// The corresponding initial output values are written into the y0 array, if it is Some.
    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_dist_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] = bw_dist_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_dist_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
            }
        }
    }
    /// Resets the satur's coeffs and each of the `N_CHANNELS` states to its initial values
    /// using the corresponding initial input value in the x0 array.
    ///
    /// The corresponding initial output values are written into the y0 array, if is Some.
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_dist_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_dist_reset_coeffs(&mut self.coeffs);
            bw_dist_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
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
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_dist_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_dist_process_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }
    /// Sets the distortion (input gain, approximately) to the given `value`.
    ///
    /// Valid range: [0.0 (low distortion), 1.0 (high distortion)].
    ///
    /// Default value: 0.0.
    #[inline(always)]
    pub fn set_distortion(&mut self, value: f32) {
        unsafe {
            bw_dist_set_distortion(&mut self.coeffs, value);
        }
    }
    /// Sets the tone (filter) to the given `value`.
    ///
    /// Valid range: [0.0 (low cutoff), 1.0 (high cutoff)].
    ///
    /// Default value: 0.5.
    #[inline(always)]
    pub fn set_tone(&mut self, value: f32) {
        unsafe {
            bw_dist_set_tone(&mut self.coeffs, value);
        }
    }
    /// Sets the volume (output gain) to the given `value`.
    ///
    /// Valid range: [0.0 (silence), 1.0 (max volume)].
    ///
    /// Default value: 1.0.
    #[inline(always)]
    pub fn set_volume(&mut self, value: f32) {
        unsafe {
            bw_dist_set_volume(&mut self.coeffs, value);
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::{c_wrapper::*, native::math::INVERSE_2_PI};

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    const N_SAMPLES: usize = 8;

    type DistT = Dist<N_CHANNELS>;
    #[test]
    fn new() {
        let dist = DistT::new();
        assert_eq!(
            dist.coeffs.state,
            bw_dist_coeffs_state_bw_dist_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.hp1_coeffs.state,
            bw_hp1_coeffs_state_bw_hp1_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.peak_coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.clip_coeffs.state,
            bw_clip_coeffs_state_bw_clip_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.satur_coeffs.state,
            bw_satur_coeffs_state_bw_satur_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.lp1_coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_init
        );
        assert_eq!(
            dist.coeffs.gain_coeffs.state,
            bw_gain_coeffs_state_bw_gain_coeffs_state_init
        );

        assert_eq!(dist.coeffs.hp1_coeffs.lp1_coeffs.cutoff, 7.0);

        assert_eq!(dist.coeffs.peak_coeffs.mm2_coeffs.svf_coeffs.cutoff, 2e3);
        assert_eq!(dist.coeffs.peak_coeffs.bandwidth, 10.0);

        assert_eq!(dist.coeffs.clip_coeffs.bias, 0.75 / 4.25);
        assert_eq!(dist.coeffs.clip_coeffs.gain, 1.0 / 4.25);
        assert_eq!(dist.coeffs.clip_coeffs.gain_compensation, 1);

        assert_eq!(dist.coeffs.satur_coeffs.gain, 1.0 / 0.7);
        assert_eq!(dist.coeffs.satur_coeffs.gain_compensation, 1);

        assert_eq!(
            dist.coeffs.lp1_coeffs.cutoff,
            475.0 + (20e3 - 475.0) * 0.125
        );
    }

    #[test]
    fn set_sample_rate() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);
        let fs_2pi = INVERSE_2_PI * SAMPLE_RATE;

        assert_eq!(
            dist.coeffs.hp1_coeffs.lp1_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            dist.coeffs
                .peak_coeffs
                .mm2_coeffs
                .svf_coeffs
                .smooth_coeffs
                .fs_2pi,
            fs_2pi
        );
        assert_eq!(dist.coeffs.clip_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(dist.coeffs.satur_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(dist.coeffs.lp1_coeffs.smooth_coeffs.fs_2pi, fs_2pi);
        assert_eq!(dist.coeffs.gain_coeffs.smooth_coeffs.fs_2pi, fs_2pi);

        assert_eq!(
            dist.coeffs.state,
            bw_dist_coeffs_state_bw_dist_coeffs_state_set_sample_rate
        );
    }

    #[test]
    fn reset() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);

        dist.reset(None, None);

        assert_eq!(
            dist.coeffs.state,
            bw_dist_coeffs_state_bw_dist_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.peak_coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.lp1_coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.gain_coeffs.state,
            bw_gain_coeffs_state_bw_gain_coeffs_state_reset_coeffs
        );
    }

    #[test]
    fn reset_multi() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);
        let x0 = [0.31, 0.92];
        let mut y0 = [0.0, 0.0];

        dist.reset_multi(&x0, Some(&mut y0));

        assert_eq!(
            dist.coeffs.state,
            bw_dist_coeffs_state_bw_dist_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.peak_coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.lp1_coeffs.state,
            bw_lp1_coeffs_state_bw_lp1_coeffs_state_reset_coeffs
        );
        assert_eq!(
            dist.coeffs.gain_coeffs.state,
            bw_gain_coeffs_state_bw_gain_coeffs_state_reset_coeffs
        );

        unsafe {
            assert_eq!(
                dist.states[0].hash,
                bw_hash_sdbm("bw_dist_state".as_ptr() as *const i8)
            );
        }

        assert_eq!(dist.states[1].coeffs_reset_id, dist.coeffs.reset_id);
    }

    #[test]
    fn process() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];
        let mut y0 = [0.0, 0.0];

        dist.reset_multi(&x0, Some(&mut y0));

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        dist.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(dist.coeffs.state >= bw_dist_coeffs_state_bw_dist_coeffs_state_reset_coeffs);
        assert_eq!(dist.coeffs.reset_id, dist.states[0].coeffs_reset_id);
    }

    #[test]
    fn set_distortion() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);

        let distortion = 0.9;

        dist.set_distortion(distortion);
        assert_eq!(dist.coeffs.peak_coeffs.peak_gain, unsafe {
            bw_dB2linf(60.0 * distortion)
        });
    }

    #[test]
    fn set_tone() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);

        let tone = 0.5;

        dist.set_tone(tone);

        assert_eq!(
            dist.coeffs.lp1_coeffs.cutoff,
            475. + (20e3 - 475.) * tone * tone * tone
        );
    }

    #[test]
    fn set_volume() {
        let mut dist = DistT::new();
        dist.set_sample_rate(SAMPLE_RATE);

        let volume = 0.8;

        dist.set_volume(volume);

        assert_eq!(dist.coeffs.gain_coeffs.gain, volume * volume * volume);
    }
}
