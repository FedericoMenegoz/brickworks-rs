use crate::native::common::{debug_assert_is_finite, debug_assert_positive};
use crate::native::one_pole::{OnePoleCoeffs, OnePoleState};

pub struct Clip<const N_CHANNELS: usize> {
    coeffs: ClipCoeffs,
    states: [ClipState; N_CHANNELS],
    _states_p: [ClipState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Clip<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        // #[cfg(debug_assertions)]
        // {
        //     debug_assert_positive(sample_rate);
        //     debug_assert_is_finite(sample_rate);
        // }
        // self.coeffs.smooth_coeffs.set_sample_rate(sample_rate);
        todo!()
    }

    #[inline(always)]
    pub fn reset(&mut self, x_0: &[f32; N_CHANNELS], y_0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        n_sample: usize,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn set_bias(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain_compensation(&mut self, value: bool) {
        todo!()
    }

    // Private methods
    #[inline(always)]
    fn do_update_coeffs(&mut self, force: bool) {
        todo!()
    }

    #[inline(always)]
    fn process1(&mut self, x: f32, channel: usize) -> f32 {
        todo!()
    }

    #[inline(always)]
    fn process1_comp(&mut self, x: f32, channel: usize) -> f32 {
        todo!()
    }

    #[inline(always)]
    fn process_multi(&mut self, x: &[&[f32]; N_CHANNELS], y: &mut [&mut [f32]; N_CHANNELS]) {
        todo!()
    }

    /*
    To understand when is needed see line 578 of bw_clip.h:

    static inline void bw_clip_process(
        bw_clip_coeffs * BW_RESTRICT coeffs,
        bw_clip_state * BW_RESTRICT  state,
        const float *                x,
        float *                      y,
        size_t                       n_samples) {...}
    */
    // #[inline(always)]
    // fn process_single(&mut self, x: &[f32], y: &mut [f32]) {
    //     todo!()
    // }
}

impl<const N_CHANNELS: usize> Default for Clip<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
pub struct ClipCoeffs {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs,
    smooth_bias_state: OnePoleState,
    smooth_gain_state: OnePoleState,

    // Coefficients
    bias_dc: f32,
    inv_gain: f32,

    // Parameters
    bias: f32,
    gain: f32,
    gain_compensation: bool,
}

pub struct ClipState {
    x_z1: f32,
    f_z1: f32,
}

#[cfg(test)]
mod tests {
    use std::f32;

    use super::Clip;
    use crate::{
        c_wrapper::{bw_clip_coeffs, bw_clip_do_update_coeffs, clip::Clip as ClipWrapper},
        native::{
            clip::ClipCoeffs,
            one_pole::{OnePoleCoeffs, tests::assert_one_pole_coeffs},
        },
    };
    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    #[test]
    fn clip_initialization() {
        let rust_clip = Clip::<N_CHANNELS>::new();
        let c_clip = ClipWrapper::<N_CHANNELS>::new();

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn set_sample_rate_valid() {
        let mut rust_clip = Clip::<N_CHANNELS>::new();
        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();

        rust_clip.set_sample_rate(SAMPLE_RATE);
        c_clip.set_sample_rate(SAMPLE_RATE);

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    #[should_panic(expected = "value must be finite, got inf")]
    fn set_sample_rate_must_be_finite() {
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        rust_clip.set_sample_rate(f32::INFINITY);
    }

    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_sample_rate_must_be_positive() {
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        rust_clip.set_sample_rate(-1.);
    }

    #[test]
    fn do_update_coeffs() {
        let bias = 1.0;
        let mut rust_clip = Clip::<N_CHANNELS>::new();
        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        rust_clip.set_bias(bias);
        c_clip.set_bias(bias);
        rust_clip.do_update_coeffs(false);
        c_clip.do_update_coeffs(false);

        assert_clip(&rust_clip, &c_clip);

        let gain = 2.0;
        rust_clip.set_gain(gain);
        c_clip.set_gain(gain);
        rust_clip.do_update_coeffs(false);
        c_clip.do_update_coeffs(false);
    }

    #[test]
    fn reset_none() {
        // let mut rust_clip = Clip::<N_CHANNELS>::new()
        let bias = 5.0;
        let gain = 10.0;
        let x0 = [10.0, 11.0];
        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset(&x0, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset(&x0, None);

        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn reset_some() {
        let bias = 0.5;
        let gain = 0.01;
        let x0 = [0.1, 0.02];
        let mut rust_y0 = [0.2, 0.7];
        let mut c_y0 = [0.9, 0.7];
        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset(&x0, Some(&mut c_y0));

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset(&x0, Some(&mut rust_y0));

        assert_clip(&rust_clip, &c_clip);
        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y0[channel], c_y0[channel]));
    }

    #[test]
    fn process1() {
        let bias = 0.5;
        let gain = 0.01;

        let x = [0.1, 0.02];
        let mut rust_y = [0.2, 0.7];
        let mut c_y = [0.9, 0.7];

        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.reset(&x, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.reset(&x, None);

        (0..N_CHANNELS).for_each(|channel| {
            c_y[channel] = c_clip.process1(x[channel], channel);
            rust_y[channel] = rust_clip.process1(x[channel], channel);
        });

        assert_clip(&rust_clip, &c_clip);
        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y[channel], c_y[channel]));
    }

    #[test]
    fn process1_comp() {
        let bias = 0.5;
        let gain = 0.01;

        let x = [0.1, 0.02];
        let mut rust_y = [0.2, 0.7];
        let mut c_y = [0.9, 0.7];

        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        let mut rust_clip = Clip::<N_CHANNELS>::new();

        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(true);
        c_clip.reset(&x, None);

        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(true);
        rust_clip.reset(&x, None);

        (0..N_CHANNELS).for_each(|channel| {
            rust_y[channel] = rust_clip.process1_comp(x[channel], channel);
            c_y[channel] = c_clip.process1_comp(x[channel], channel);
            assert_clip(&rust_clip, &c_clip);
        });

        (0..N_CHANNELS).for_each(|channel| assert_eq!(rust_y[channel], c_y[channel]));
    }

    #[test]
    fn process_with_comp() {
        let n_samples = 2;
        let bias = 0.5;
        let gain = 0.01;
        let x_ch0 = [0.1, 0.02];
        let x_ch1 = [0.1, 0.02];
        let x: [&[f32]; N_CHANNELS] = [&x_ch0, &x_ch1];

        let mut rust_y_ch0 = [0.2, 0.7];
        let mut rust_y_ch1 = [0.2, 0.7];
        let mut rust_y: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut rust_y_ch0), Some(&mut rust_y_ch1)];

        let mut c_y_ch0 = [0.2, 0.7];
        let mut c_y_ch1 = [0.2, 0.7];
        let mut c_y: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut c_y_ch0), Some(&mut c_y_ch1)];

        let mut rust_clip = Clip::<N_CHANNELS>::new();
        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(true);
        rust_clip.reset(&[0.0; N_CHANNELS], None);

        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(true);
        c_clip.reset(&[0.0; N_CHANNELS], None);

        rust_clip.process(&x, Some(&mut c_y), n_samples);
        c_clip.process(&x, Some(&mut rust_y), n_samples);

        (0..N_CHANNELS).for_each(|channel| {
            (0..n_samples).for_each(|sample| {
                assert_eq!(
                    rust_y[channel].as_ref().unwrap()[sample],
                    c_y[channel].as_ref().unwrap()[sample]
                );
            });
        });
        assert_clip(&rust_clip, &c_clip);
    }

    #[test]
    fn process_without_comp() {
        let n_samples = 2;
        let bias = 0.5;
        let gain = 0.01;
        let x_ch0 = [0.1, 0.02];
        let x_ch1 = [0.1, 0.02];
        let x: [&[f32]; N_CHANNELS] = [&x_ch0, &x_ch1];

        let mut rust_y_ch0 = [0.2, 0.7];
        let mut rust_y_ch1 = [0.2, 0.7];
        let mut rust_y: [Option<&mut [f32]>; N_CHANNELS] =
            [Some(&mut rust_y_ch0), Some(&mut rust_y_ch1)];

        let mut c_y_ch0 = [0.2, 0.7];
        let mut c_y_ch1 = [0.2, 0.7];
        let mut c_y: [Option<&mut [f32]>; N_CHANNELS] = [Some(&mut c_y_ch0), Some(&mut c_y_ch1)];

        let mut rust_clip = Clip::<N_CHANNELS>::new();
        rust_clip.set_bias(bias);
        rust_clip.set_gain(gain);
        rust_clip.set_gain_compensation(false); // it is default, but to be clear
        rust_clip.reset(&[0.0; N_CHANNELS], None);

        let mut c_clip = ClipWrapper::<N_CHANNELS>::new();
        c_clip.set_bias(bias);
        c_clip.set_gain(gain);
        c_clip.set_gain_compensation(false); // it is default, but to be clear
        c_clip.reset(&[0.0; N_CHANNELS], None);

        rust_clip.process(&x, Some(&mut c_y), n_samples);
        c_clip.process(&x, Some(&mut rust_y), n_samples);

        (0..N_CHANNELS).for_each(|channel| {
            (0..n_samples).for_each(|sample| {
                assert_eq!(
                    rust_y[channel].as_ref().unwrap()[sample],
                    c_y[channel].as_ref().unwrap()[sample]
                );
            });
        });
        assert_clip(&rust_clip, &c_clip);
    }

    fn assert_clip_coeffs(rust_coeffs: &ClipCoeffs, c_coeffs: &bw_clip_coeffs) {
        let pre_message = "clip.coeff.";
        let post_message = "does not match";
        assert_eq!(
            rust_coeffs.bias_dc, c_coeffs.bias_dc,
            "{}bias_dc {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.inv_gain, c_coeffs.inv_gain,
            "{}inv_gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.bias, c_coeffs.bias,
            "{}bias {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain, c_coeffs.gain,
            "{}gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain_compensation,
            c_coeffs.gain_compensation != 0,
            "{}gain_compensation {}",
            pre_message,
            post_message
        );
        assert_one_pole_coeffs(rust_coeffs.smooth_coeffs, c_coeffs.smooth_coeffs);
        assert_eq!(
            rust_coeffs.smooth_bias_state.y_z1,
            c_coeffs.smooth_bias_state.y_z1
        );
        assert_eq!(
            rust_coeffs.smooth_gain_state.y_z1,
            c_coeffs.smooth_gain_state.y_z1
        );
    }

    fn assert_clip<const N_CHANNELS: usize>(
        rust_clip: &Clip<N_CHANNELS>,
        c_clip: &ClipWrapper<N_CHANNELS>,
    ) {
        assert_clip_coeffs(&rust_clip.coeffs, &c_clip.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_clip.states[channel].x_z1, c_clip.states[channel].x_z1);
            assert_eq!(rust_clip.states[channel].f_z1, c_clip.states[channel].F_z1);
            assert_eq!(
                rust_clip._states_p[channel].x_z1,
                c_clip.states_p[channel].x_z1
            );
            assert_eq!(
                rust_clip._states_p[channel].f_z1,
                c_clip.states_p[channel].F_z1
            );
        });
    }
}
