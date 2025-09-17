use crate::c_wrapper::{bw_gain_coeffs, bw_one_pole_sticky_mode};

pub struct Gain<const N_CHANNELS: usize> {
    coeffs: bw_gain_coeffs,
}

impl<const N_CHANNELS: usize> Gain<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        todo!()
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn reset(&mut self) {
        todo!()
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &[&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain_lin(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_gain_db(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_smooth_tau(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_sticky_thresh(&mut self, value: f32) {
        todo!()
    }

    #[inline(always)]
    pub fn set_sticky_mode(&mut self, value: bw_one_pole_sticky_mode) {
        todo!()
    }

    #[inline(always)]
    pub fn get_gain_lin(&self) -> f32 {
        todo!()
    }

    #[inline(always)]
    pub fn get_gain_cur(&self) -> f32 {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for Gain<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for bw_gain_coeffs {
    fn default() -> Self {
        Self {
            smooth_coeffs: Default::default(),
            smooth_state: Default::default(),
            gain: Default::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        c_wrapper::{
            bw_dB2linf, bw_gain_coeffs, bw_gain_init, bw_gain_process1,
            bw_gain_sticky_mode_bw_gain_sticky_mode_abs, bw_gain_update_coeffs_audio,
            bw_gain_update_coeffs_ctrl, bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs,
            bw_rcpf, gain::Gain,
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
        let tau_default = 0.005;

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
        assert_eq!(gain.coeffs.smooth_coeffs.cutoff_up, f32::INFINITY);
    }

    #[test]
    fn process() {
        const N_SAMPLES: usize = 2;
        // Wrapper
        let mut gain = GainT::new();

        let sample_0: [f32; N_CHANNELS] = [1.0, 1.0];
        let sample_1: [f32; N_CHANNELS] = [0.0, 0.0];
        let x: [&[f32]; 2] = [&sample_0, &sample_1];

        let mut out_0: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut out_1: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut y: [&mut [f32]; 2] = [&mut out_0, &mut out_1];

        // C
        let mut coeffs = bw_gain_coeffs::default();

        let sample_0_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let sample_1_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let x_c: [&[f32]; 2] = [&sample_0_c, &sample_1_c];

        let mut out_0_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut out_1_c: [f32; N_CHANNELS] = [0.0, 0.0];
        let mut y_c: [&mut [f32]; 2] = [&mut out_0_c, &mut out_1_c];

        // Process wrapper
        gain.process(&x, &mut y, N_SAMPLES);

        // Process C
        unsafe {
            bw_gain_init(&mut coeffs);
            bw_gain_update_coeffs_ctrl(&mut coeffs);
            (0..N_SAMPLES).for_each(|sample| {
                bw_gain_update_coeffs_audio(&mut coeffs);
                (0..N_CHANNELS).for_each(|channel| {
                    y_c[channel][sample] = bw_gain_process1(&coeffs, x_c[channel][sample])
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

        let gain_sticky_mode = bw_gain_sticky_mode_bw_gain_sticky_mode_abs;

        gain.set_sticky_mode(gain_sticky_mode);

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
        gain.process(&x, &mut y, n_samples);
        assert_eq!(gain.get_gain_cur(), gain.coeffs.smooth_state.y_z1);
    }
}
