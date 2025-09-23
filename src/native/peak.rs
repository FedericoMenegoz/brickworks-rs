use bitflags::bitflags;

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};
use crate::native::{
    math::{db2linf, pow2f, rcpf, sqrtf},
    mm2::{MM2Coeffs, MM2State},
};

pub struct Peak<const N_CHANNELS: usize> {
    coeffs: PeakCoeffs<N_CHANNELS>,
    states: [PeakState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Peak<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            coeffs: PeakCoeffs::new(),
            states: [PeakState::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.set_sample_rate(sample_rate);
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self
                    .coeffs
                    .reset_state(&mut self.states[channel], x0.unwrap_or(0.0));
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.coeffs
                    .reset_state(&mut self.states[channel], x0.unwrap_or(0.0));
            }),
        }
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32; N_CHANNELS], y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_coeffs();
        self.coeffs.reset_state_multi(&mut self.states, x0, y0);
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.coeffs.process_multi(&mut self.states, x, y, n_samples);
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        self.coeffs.set_cutoff(value);
    }

    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        self.coeffs.set_q(value);
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.coeffs.set_prewarp_at_cutoff(value);
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        self.coeffs.set_prewarp_freq(value);
    }

    #[inline(always)]
    pub fn set_peak_gain_lin(&mut self, value: f32) {
        self.coeffs.set_peak_gain_lin(value);
    }

    #[inline(always)]
    pub fn set_peak_gain_db(&mut self, value: f32) {
        self.coeffs.set_peak_gain_db(value);
    }

    #[inline(always)]
    pub fn set_bandwidth(&mut self, value: f32) {
        self.coeffs.set_bandwidth(value);
    }

    #[inline(always)]
    pub fn set_use_bandwidth(&mut self, value: bool) {
        self.coeffs.set_use_bandwidth(value);
    }
}

impl<const N_CHANNELS: usize> Default for Peak<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}
pub struct PeakCoeffs<const N_CHANNELS: usize> {
    // Sub-components
    mm2_coeffs: MM2Coeffs<N_CHANNELS>,

    // Coefficients
    bw_k: f32,

    // Parameters
    q: f32,
    peak_gain: f32,
    bandwidth: f32,
    use_bandwidth: bool,
    param_changed: ParamChanged,
}

impl<const N_CHANNELS: usize> PeakCoeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mm2_coeffs = MM2Coeffs::<N_CHANNELS>::new();

        let q = 0.5;
        let peak_gain = 1.0;
        let bandwidth = 2.543_106_6/*06327224*/;
        let use_bandwidth = true;

        let param_changed = ParamChanged::all();

        Self {
            mm2_coeffs,
            bw_k: 0.0,
            q,
            peak_gain,
            bandwidth,
            use_bandwidth,
            param_changed,
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(sample_rate.is_finite());
            debug_assert_positive(sample_rate);
        }

        self.mm2_coeffs.set_sample_rate(sample_rate);
    }

    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.param_changed = ParamChanged::all();
        self.update_mm2_params();
        self.mm2_coeffs.reset_coeffs();
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut PeakState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        let y = self.mm2_coeffs.reset_state(&mut state.mm2_state, x0);

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [PeakState],
        x0: &[f32; N_CHANNELS],
        y0: Option<&mut [f32; N_CHANNELS]>,
    ) {
        match y0 {
            Some(y) => (0..N_CHANNELS).for_each(|channel| {
                y[channel] = self.reset_state(&mut states[channel], x0[channel]);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.reset_state(&mut states[channel], x0[channel]);
            }),
        }
    }

    #[inline(always)]
    pub fn update_coeffs_ctrl(&mut self) {
        self.update_mm2_params();
        self.mm2_coeffs.update_coeffs_ctrl();
    }

    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.mm2_coeffs.update_coeffs_audio();
    }

    #[inline(always)]
    pub fn process1(&mut self, state: &mut PeakState, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let y = self.mm2_coeffs.process1(&mut state.mm2_state, x);

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn process(&mut self, state: &mut PeakState, x: &[f32], y: &mut [f32], n_samples: usize) {
        self.update_coeffs_ctrl();
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample]);
        });
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [PeakState],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_samples: usize,
    ) {
        self.update_coeffs_ctrl();
        (0..n_samples).for_each(|sample| {
            self.update_coeffs_audio();
            (0..N_CHANNELS).for_each(|channel| {
                y[channel][sample] = self.process1(&mut states[channel], x[channel][sample])
            });
        });
    }

    #[inline(always)]
    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e12, value);
        }

        self.mm2_coeffs.set_cutoff(value);
    }

    #[inline(always)]
    pub fn set_q(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e6, value);
        }

        if self.q != value {
            self.q = value;
            self.param_changed |= ParamChanged::Q
        }
    }

    #[inline(always)]
    pub fn set_prewarp_at_cutoff(&mut self, value: bool) {
        self.mm2_coeffs.set_prewarp_at_cutoff(value);
    }

    #[inline(always)]
    pub fn set_prewarp_freq(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=1e12, value);
        }
        self.mm2_coeffs.set_prewarp_freq(value);
    }

    #[inline(always)]
    pub fn set_peak_gain_lin(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-30..=1e30, value);
        }

        if self.peak_gain != value {
            self.peak_gain = value;
            self.param_changed |= ParamChanged::PEAK_GAIN;
        }
    }

    #[inline(always)]
    pub fn set_peak_gain_db(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(-600.0..=600.0, value);
        }
        self.set_peak_gain_lin(db2linf(value));
    }

    #[inline(always)]
    pub fn set_bandwidth(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(1e-6..=90.0, value);
        }

        if self.bandwidth != value {
            self.bandwidth = value;
            self.param_changed |= ParamChanged::BANDWIDTH;
        }
    }

    #[inline(always)]
    pub fn set_use_bandwidth(&mut self, value: bool) {
        if self.use_bandwidth != value {
            self.use_bandwidth = value;
            self.param_changed |= ParamChanged::Q | ParamChanged::BANDWIDTH;
        }
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    #[inline(always)]
    pub fn coeffs_is_valid(&mut self) -> bool {
        todo!()
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    #[inline(always)]
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }

    // Private
    #[inline(always)]
    fn update_mm2_params(&mut self) {
        if !self.param_changed.is_empty() {
            if self.use_bandwidth {
                if self
                    .param_changed
                    .intersects(ParamChanged::PEAK_GAIN | ParamChanged::BANDWIDTH)
                {
                    if self.param_changed.contains(ParamChanged::BANDWIDTH) {
                        self.bw_k = pow2f(self.bandwidth);
                    }
                    let q = sqrtf(self.bw_k * self.peak_gain) * rcpf(self.bw_k - 1.0);
                    self.mm2_coeffs.set_q(q);
                    self.mm2_coeffs
                        .set_coeff_bp((self.peak_gain - 1.0) * rcpf(q));
                }
            } else if self
                .param_changed
                .intersects(ParamChanged::PEAK_GAIN | ParamChanged::Q)
            {
                if self.param_changed.contains(ParamChanged::Q) {
                    self.mm2_coeffs.set_q(self.q);
                }
                self.mm2_coeffs
                    .set_coeff_bp((self.peak_gain - 1.0) * rcpf(self.q));
            }

            self.param_changed = ParamChanged::empty();
        }
    }
}

impl<const N_CHANNELS: usize> Default for PeakCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    struct ParamChanged: u32 {
        const Q = 1;
        const PEAK_GAIN = 1<<1;
        const BANDWIDTH = 1<<2;

        // added this as peak_initialization test fail caused by param_changed
        const _ = !0;
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct PeakState {
    //Sub-components
    mm2_state: MM2State,
}

#[cfg(test)]
mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_peak_coeffs, peak::Peak as PeakWrapper},
        native::mm2::tests::{assert_mm2_coeffs, assert_mm2_state},
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type PeakT = Peak<N_CHANNELS>;
    type PeakWrapperT = PeakWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_peak = PeakT::new();
        let c_peak = PeakWrapperT::new();

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn reset_none() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();

        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.reset(None, None);
        c_peak.reset(None, None);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn reset_some() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();

        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.3;
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_peak.reset(Some(x0), Some(&mut rust_y0));
        c_peak.reset(Some(x0), Some(&mut c_y0));

        assert_peak(&rust_peak, &c_peak);
        assert_eq!(rust_y0, c_y0);
    }

    #[test]
    fn reset_multi() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.045, 0.045];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_peak.reset_multi(&x0, Some(&mut rust_y0));
        c_peak.reset_multi(&x0, Some(&mut c_y0));

        assert_peak(&rust_peak, &c_peak);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let bw = 39.0;
        let cutoff = 456.1;
        let prewarp_freq = 460.1;
        let gain_db = 56.0;

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];

        rust_peak.set_bandwidth(bw);
        rust_peak.set_cutoff(cutoff);
        rust_peak.set_peak_gain_db(gain_db);
        rust_peak.set_prewarp_at_cutoff(true);
        rust_peak.set_prewarp_freq(prewarp_freq);
        rust_peak.set_use_bandwidth(true);

        c_peak.set_bandwidth(bw);
        c_peak.set_cutoff(cutoff);
        c_peak.set_peak_gain_db(gain_db);
        c_peak.set_prewarp_at_cutoff(true);
        c_peak.set_prewarp_freq(prewarp_freq);
        c_peak.set_use_bandwidth(true);

        rust_peak.reset(None, None);
        c_peak.reset(None, None);
        rust_peak.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_peak.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_peak(&rust_peak, &c_peak);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    fn set_cutoff() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let cutoff = 261.63;

        rust_peak.set_cutoff(cutoff);
        c_peak.set_cutoff(cutoff);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got -1e0")]
    #[test]
    fn set_cutoff_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_cutoff(-1.0);
    }

    #[test]
    fn set_q() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let q = 12.0;

        rust_peak.set_q(q);
        c_peak.set_q(q);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e6], got 1e7")]
    #[test]
    fn set_q_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_q(1e7);
    }

    #[test]
    fn set_prewarp_at_cutoff() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_prewarp_at_cutoff(true);
        c_peak.set_prewarp_at_cutoff(true);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_prewarp_freq() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let prewarp_freq = 340.12;
        rust_peak.set_prewarp_freq(prewarp_freq);
        c_peak.set_prewarp_freq(prewarp_freq);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-6, 1e12], got 1e13")]
    #[test]
    fn set_prewarp_freq_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_prewarp_freq(1e13);
    }

    #[test]
    fn set_peak_gain_lin() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let gain_lin = 1.76;

        rust_peak.set_peak_gain_lin(gain_lin);
        c_peak.set_peak_gain_lin(gain_lin);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [1e-30, 1e30], got 0e0")]
    #[test]
    fn set_peak_gain_lin_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_peak_gain_lin(0.0);
    }

    #[test]
    fn set_peak_gain_db() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        let gain_db = 450.0;

        rust_peak.set_peak_gain_lin(gain_db);
        c_peak.set_peak_gain_lin(gain_db);

        assert_peak(&rust_peak, &c_peak);
    }

    #[should_panic(expected = "value must be in range [-6e2, 6e2], got 6.001e2")]
    #[test]
    fn set_peak_gain_db_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        rust_peak.set_peak_gain_db(600.1);
    }

    #[should_panic(expected = "value must be in range [1e-6, 9e1], got 9.01e1")]
    #[test]
    fn set_bandwidth_invalid() {
        let mut rust_peak = PeakT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_bandwidth(90.1);
    }

    #[test]
    fn set_bandwidth() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);
        let bw = 34.12;

        rust_peak.set_peak_gain_lin(bw);
        c_peak.set_peak_gain_lin(bw);

        assert_peak(&rust_peak, &c_peak);
    }

    #[test]
    fn set_use_bandwidth() {
        let mut rust_peak = PeakT::new();
        let mut c_peak = PeakWrapperT::new();
        rust_peak.set_sample_rate(SAMPLE_RATE);
        c_peak.set_sample_rate(SAMPLE_RATE);

        rust_peak.set_use_bandwidth(true);
        c_peak.set_use_bandwidth(true);

        assert_peak(&rust_peak, &c_peak);
    }

    fn assert_peak_coeff<const N_CHANNELS: usize>(
        rust_coeffs: &PeakCoeffs<N_CHANNELS>,
        c_coeffs: &bw_peak_coeffs,
    ) {
        // Sub-components
        assert_mm2_coeffs(&rust_coeffs.mm2_coeffs, &c_coeffs.mm2_coeffs);

        assert_eq!(rust_coeffs.bw_k, c_coeffs.bw_k);
        assert_eq!(rust_coeffs.q, c_coeffs.Q);
        assert_eq!(rust_coeffs.peak_gain, c_coeffs.peak_gain);
        assert_eq!(rust_coeffs.bandwidth, c_coeffs.bandwidth);
        assert_eq!(
            rust_coeffs.use_bandwidth,
            match c_coeffs.use_bandwidth {
                0 => false,
                _ => true,
            }
        );
        assert_eq!(
            rust_coeffs.param_changed.bits(),
            c_coeffs.param_changed as u32
        );
    }

    fn assert_peak<const N_CHANNELS: usize>(
        rust_peak: &Peak<N_CHANNELS>,
        c_peak: &PeakWrapper<N_CHANNELS>,
    ) {
        assert_peak_coeff(&rust_peak.coeffs, &c_peak.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_mm2_state(
                &rust_peak.states[channel].mm2_state,
                &c_peak.states[channel].mm2_state,
            );
        });
    }
}
