use crate::c_wrapper::{bw_svf_coeffs, bw_svf_state};

pub struct SVF<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_svf_coeffs,
    pub(crate) states: [bw_svf_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SVF<N_CHANNELS> {
    pub fn new() -> Self {
        todo!();
    }

    pub fn set_sample_rate(&mut self, value: f32) {
        todo!()
    }

    pub fn reset(
        &mut self,
        x0: f32,
        y_lp0: Option<&mut [f32; N_CHANNELS]>,
        y_bp0: Option<&mut [f32; N_CHANNELS]>,
        y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        todo!()
    }

    pub fn reset_multi(
        &mut self,
        x0: &[f32; N_CHANNELS],
        y_lp0: Option<&mut [f32; N_CHANNELS]>,
        y_bp0: Option<&mut [f32; N_CHANNELS]>,
        y_hp0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        todo!()
    }

    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y_lp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_bp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
        y_hp: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
    ) {
        todo!()
    }

    pub fn set_cutoff(&mut self, value: f32) {
        todo!()
    }

    pub fn set_q(&mut self, value: f32) {
        todo!()
    }

    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        todo!()
    }

    pub fn set_prewarp_freq(&mut self, value: f32) {
        todo!()
    }
}

impl Default for bw_svf_state {
    fn default() -> Self {
        todo!()
    }
}

impl Default for bw_svf_coeffs {
    fn default() -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::c_wrapper::{
        bw_rcpf, bw_svf_coeffs, bw_svf_init, bw_svf_process1, bw_svf_set_cutoff,
        bw_svf_set_prewarp_at_cutoff, bw_svf_set_prewarp_freq, bw_svf_state,
        bw_svf_update_coeffs_audio, svf::SVF,
    };
    use std::f32::consts::PI;

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;
    const GAIN: f32 = 200_000.0;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    type SVFTest = SVF<N_CHANNELS>;

    #[test]
    fn new() {
        let mut svf = SVFTest::new();
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

        assert_eq!(svf.coeffs.smooth_coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE)
    }

    #[test]
    fn reset() {
        let mut svf = SVFTest::new();
        let x0 = 0.1;
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
        svf.set_cutoff(cutoff);
        svf.set_prewarp_at_cutoff(true);
        svf.set_prewarp_freq(prewarp_freq);
        svf.process(&x, Some(&mut y_lp), Some(&mut y_bp), Some(&mut y_hp));

        // Process C
        let mut c_coeffs = bw_svf_coeffs::default();
        let mut c_states = [bw_svf_state::default(); N_CHANNELS];

        unsafe {
            bw_svf_init(&mut c_coeffs);
            bw_svf_set_cutoff(&mut c_coeffs, cutoff);
            bw_svf_set_prewarp_freq(&mut c_coeffs, prewarp_freq);
            bw_svf_set_prewarp_at_cutoff(&mut c_coeffs, 1);

            (0..n_samples).for_each(|sample| {
                bw_svf_update_coeffs_audio(&mut c_coeffs);
                (0..N_CHANNELS).for_each(|channel| {
                    bw_svf_process1(
                        &c_coeffs,
                        &mut c_states[channel],
                        x[channel][sample],
                        &mut c_y_lp[channel][sample],
                        &mut c_y_bp[channel][sample],
                        &mut c_y_hp[channel][sample],
                    )
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

        svf.set_cutoff(q);

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
