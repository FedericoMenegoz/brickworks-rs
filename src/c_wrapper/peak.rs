use std::ptr::null_mut;

use super::*;

pub struct Peak<const N_CHANNELS: usize> {
    coeffs: bw_peak_coeffs,
    states: [bw_peak_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Peak<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut coeffs = bw_peak_coeffs::default();
        unsafe {
            bw_peak_init(&mut coeffs);
        }

        Self {
            coeffs,
            states: [bw_peak_state::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_peak_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            bw_peak_reset_coeffs(&mut self.coeffs);
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] = bw_peak_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_peak_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
            }
        }
    }
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_peak_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_peak_reset_coeffs(&mut self.coeffs);
            bw_peak_reset_state_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x0.as_ptr(),
                y0_ptrs,
                N_CHANNELS,
            );
        }
    }
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_peak_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_peak_process_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_peak_set_cutoff(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        unsafe {
            bw_peak_set_Q(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        unsafe {
            bw_peak_set_prewarp_at_cutoff(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }
    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        unsafe {
            bw_peak_set_prewarp_freq(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_peak_gain_lin(&mut self, value: f32) {
        unsafe {
            bw_peak_set_peak_gain_lin(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_peak_gain_db(&mut self, value: f32) {
        unsafe {
            bw_peak_set_peak_gain_dB(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_bandwidth(&mut self, value: f32) {
        unsafe {
            bw_peak_set_bandwidth(&mut self.coeffs, value);
        }
    }
    #[inline(always)]
    pub fn set_use_bandwidth(&mut self, value: bool) {
        unsafe {
            bw_peak_set_use_bandwidth(&mut self.coeffs, if value { 1 } else { 0 });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::math::INVERSE_2_PI;

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    const N_SAMPLES: usize = 8;

    type PeakT = Peak<N_CHANNELS>;

    #[test]
    pub fn new() {
        let peak = PeakT::new();

        let q_default = 0.5;
        let peak_gain_default = 1.0;
        let bandwidth_default = 2.543_106_6/*06327224*/;

        assert_eq!(peak.coeffs.Q, q_default);
        assert_eq!(peak.coeffs.peak_gain, peak_gain_default);
        assert_eq!(peak.coeffs.bandwidth, bandwidth_default);
        assert!(peak.coeffs.use_bandwidth != 0);
        assert_eq!(peak.coeffs.param_changed, !0);
        assert_eq!(
            peak.coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_init
        );
    }

    #[test]
    pub fn set_sample_rate() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let fs_2pi = INVERSE_2_PI * SAMPLE_RATE;

        assert_eq!(
            peak.coeffs.mm2_coeffs.svf_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            peak.coeffs.mm2_coeffs.gain_x_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            peak.coeffs.mm2_coeffs.gain_lp_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            peak.coeffs.mm2_coeffs.gain_bp_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            peak.coeffs.mm2_coeffs.gain_hp_coeffs.smooth_coeffs.fs_2pi,
            fs_2pi
        );
        assert_eq!(
            peak.coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_set_sample_rate
        );
    }

    #[test]
    pub fn reset_none() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        peak.reset(None, None);

        assert_eq!(
            peak.coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs
        );
        assert_eq!(peak.coeffs.reset_id, peak.states[0].coeffs_reset_id);
    }

    #[test]
    pub fn reset_some() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let x0 = 4.0;
        let mut y0 = [0.0, 0.0];
        peak.reset(Some(x0), Some(&mut y0));

        assert_eq!(
            peak.coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs
        );

        assert_eq!(peak.coeffs.reset_id, peak.states[0].coeffs_reset_id);
    }

    #[test]
    pub fn reset_multi() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.1, 0.2];
        let mut y0 = [0.0, 0.0];
        peak.reset_multi(&x0, Some(&mut y0));

        assert_eq!(y0, x0);
        assert_eq!(
            peak.coeffs.state,
            bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs
        );
        assert_eq!(peak.coeffs.reset_id, peak.states[1].coeffs_reset_id);
    }

    #[test]
    pub fn process() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.0, 0.0];

        peak.reset_multi(&x0, None);

        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0], &mut [0.0, 0.0]];
        peak.process(&PULSE_INPUT, &mut y, N_SAMPLES);
        assert!(peak.coeffs.state >= bw_peak_coeffs_state_bw_peak_coeffs_state_reset_coeffs);
        assert_eq!(peak.coeffs.reset_id, peak.states[0].coeffs_reset_id);
    }

    #[test]
    pub fn set_cutoff() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let cutoff = 493.883;

        peak.set_cutoff(cutoff);

        assert_eq!(peak.coeffs.mm2_coeffs.svf_coeffs.cutoff, cutoff);
        assert!(peak.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    fn set_q() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let q = 1.4;
        peak.set_q(q);

        assert_eq!(peak.coeffs.Q, q);
        assert!(peak.coeffs.state >= bw_peak_coeffs_state_bw_peak_coeffs_state_init);
    }

    #[test]
    pub fn set_prewarp_at_cutoff() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        peak.set_prewarp_at_cutoff(true);
        assert_eq!(peak.coeffs.mm2_coeffs.svf_coeffs.prewarp_k, 1.0);

        peak.set_prewarp_at_cutoff(false);
        assert_eq!(peak.coeffs.mm2_coeffs.svf_coeffs.prewarp_k, 0.0);

        assert!(peak.coeffs.state >= bw_mm2_coeffs_state_bw_mm2_coeffs_state_init);
    }

    #[test]
    pub fn set_prewarp_freq() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        let freq = 185.0;
        peak.set_prewarp_freq(freq);
        assert_eq!(peak.coeffs.mm2_coeffs.svf_coeffs.prewarp_freq, freq);

        peak.set_prewarp_at_cutoff(false);
        assert_eq!(peak.coeffs.mm2_coeffs.svf_coeffs.prewarp_k, 0.0);

        assert!(peak.coeffs.state >= bw_peak_coeffs_state_bw_peak_coeffs_state_init);
    }

    #[test]
    fn set_peak_gain_lin() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);
        let peak_gain = 1.2;
        peak.set_peak_gain_lin(peak_gain);

        assert_eq!(peak.coeffs.peak_gain, peak_gain);
        assert!(peak.coeffs.param_changed as u32 & BW_PEAK_PARAM_PEAK_GAIN != 0);
    }

    #[test]
    fn set_peak_gain_db() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);
        let peak_gain = 300.0;
        peak.set_peak_gain_db(peak_gain);

        unsafe {
            assert_eq!(peak.coeffs.peak_gain, bw_dB2linf(peak_gain));
        }
        assert!(peak.coeffs.param_changed as u32 & BW_PEAK_PARAM_PEAK_GAIN != 0);
    }

    #[test]
    fn set_bandwidth() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);
        let bandwidth = 40.0;
        peak.set_bandwidth(bandwidth);

        assert_eq!(peak.coeffs.bandwidth, bandwidth);
        assert!(peak.coeffs.param_changed as u32 & BW_PEAK_PARAM_BANDWIDTH != 0)
    }

    #[test]
    fn set_use_bandwidth() {
        let mut peak = PeakT::new();
        peak.set_sample_rate(SAMPLE_RATE);

        peak.set_use_bandwidth(true);

        assert!(peak.coeffs.use_bandwidth != 0);
        assert!(peak.coeffs.param_changed as u32 & (BW_PEAK_PARAM_Q | BW_PEAK_PARAM_BANDWIDTH) != 0)
    }
}
