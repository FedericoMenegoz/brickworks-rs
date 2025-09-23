use crate::native::{
    clip::{ClipCoeffs, ClipState},
    gain::GainCoeffs,
    hp1::{HP1Coeffs, HP1State},
    lp1::{LP1Coeffs, LP1State},
    peak::{PeakCoeffs, PeakState},
    satur::{SaturCoeffs, SaturState},
};

#[cfg(debug_assertions)]
use super::common::{debug_assert_is_finite, debug_assert_positive, debug_assert_range};

pub struct Dist<const N_CHANNELS: usize> {
    coeffs: DistCoeffs<N_CHANNELS>,
    states: [DistState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> Dist<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            coeffs: DistCoeffs::new(),
            states: [DistState::default(); N_CHANNELS],
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
    pub fn set_distortion(&mut self, value: f32) {
        self.coeffs.set_distortion(value);
    }

    #[inline(always)]
    pub fn set_tone(&mut self, value: f32) {
        self.coeffs.set_tone(value);
    }

    #[inline(always)]
    pub fn set_volume(&mut self, value: f32) {
        self.coeffs.set_volume(value);
    }
}

impl<const N_CHANNELS: usize> Default for Dist<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct DistCoeffs<const N_CHANNELS: usize> {
    hp1_coeffs: HP1Coeffs<N_CHANNELS>,
    peak_coeffs: PeakCoeffs<N_CHANNELS>,
    clip_coeffs: ClipCoeffs<N_CHANNELS>,
    satur_coeffs: SaturCoeffs<N_CHANNELS>,
    lp1_coeffs: LP1Coeffs<N_CHANNELS>,
    gain_coeffs: GainCoeffs<N_CHANNELS>,
}

impl<const N_CHANNELS: usize> DistCoeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new() -> Self {
        let mut hp1_coeffs = HP1Coeffs::new();
        let mut peak_coeffs = PeakCoeffs::new();
        let mut clip_coeffs = ClipCoeffs::new();
        let mut satur_coeffs = SaturCoeffs::new();
        let mut lp1_coeffs = LP1Coeffs::new();
        let gain_coeffs = GainCoeffs::new();
        hp1_coeffs.set_cutoff(7.0);
        peak_coeffs.set_cutoff(2e3);
        peak_coeffs.set_bandwidth(10.0);

        clip_coeffs.set_bias(0.75 / 4.25);
        clip_coeffs.set_gain(1.0 / 4.25);
        clip_coeffs.set_gain_compensation(true);

        satur_coeffs.set_gain(1.0 / 0.7);
        satur_coeffs.set_gain_compensation(true);

        lp1_coeffs.set_cutoff(475.0 + (20e3 - 475.0) * 0.125);

        Self {
            hp1_coeffs,
            peak_coeffs,
            clip_coeffs,
            satur_coeffs,
            lp1_coeffs,
            gain_coeffs,
        }
    }

    #[inline(always)]
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert!(sample_rate.is_finite());
            debug_assert_positive(sample_rate);
        }
        self.hp1_coeffs.set_sample_rate(sample_rate);
        self.peak_coeffs.set_sample_rate(sample_rate);
        self.clip_coeffs.set_sample_rate(sample_rate);
        self.satur_coeffs.set_sample_rate(sample_rate);
        self.lp1_coeffs.set_sample_rate(sample_rate);
        self.gain_coeffs.set_sample_rate(sample_rate);
        self.hp1_coeffs.reset_coeffs();
        self.clip_coeffs.reset_coeffs();
        self.satur_coeffs.reset_coeffs();
    }

    #[inline(always)]
    pub fn reset_coeffs(&mut self) {
        self.peak_coeffs.reset_coeffs();
        self.lp1_coeffs.reset_coeffs();
        self.gain_coeffs.reset_coeffs();
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut DistState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        let mut y = self.hp1_coeffs.reset_state(&mut state.hp1_state, x0);
        y = self.peak_coeffs.reset_state(&mut state.peak_state, y);
        y = self.clip_coeffs.reset_state(&mut state.clip_state, y);
        y = self.satur_coeffs.reset_state(&mut state.satur_state, y);
        y = self.lp1_coeffs.reset_state(&mut state.lp1_state, y);
        y *= self.gain_coeffs.get_gain_cur();

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [DistState; N_CHANNELS],
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
        self.peak_coeffs.update_coeffs_ctrl();

        // Not implemented yet: C version only contained assertions
        // need to revisit which assertions from the C version make sense to keep in Rust
        // self.lp1_coeffs.update_coeffs_ctrl();

        self.gain_coeffs.update_coeffs_ctrl();
    }

    #[inline(always)]
    pub fn update_coeffs_audio(&mut self) {
        self.peak_coeffs.update_coeffs_audio();
        self.lp1_coeffs.update_coeffs_audio();
        self.gain_coeffs.update_coeffs_audio();
    }

    #[inline(always)]
    pub fn process1(&mut self, state: &mut DistState, x: f32) -> f32 {
        debug_assert!(x.is_finite());

        let mut y = self.hp1_coeffs.process1(&mut state.hp1_state, x);
        y = self.peak_coeffs.process1(&mut state.peak_state, y);
        y = self.clip_coeffs.process1_comp(&mut state.clip_state, y);
        y = self.satur_coeffs.process1_comp(&mut state.satur_state, y);
        y = self.lp1_coeffs.process1(&mut state.lp1_state, y);
        y = self.gain_coeffs.process1(y);

        debug_assert!(y.is_finite());

        y
    }

    #[inline(always)]
    pub fn process(&mut self, state: &mut DistState, x: &[f32], y: &mut [f32], n_sample: usize) {
        self.update_coeffs_ctrl();
        (0..n_sample).for_each(|sample| {
            self.update_coeffs_audio();
            y[sample] = self.process1(state, x[sample])
        });
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [DistState; N_CHANNELS],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_sample: usize,
    ) {
        self.update_coeffs_ctrl();
        (0..n_sample).for_each(|sample| {
            self.update_coeffs_audio();
            (0..N_CHANNELS).for_each(|channel| {
                y[channel][sample] = self.process1(&mut states[channel], x[channel][sample])
            });
        });
    }

    #[inline(always)]
    pub fn set_distortion(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(0.0..=1.0, value);
        }

        self.peak_coeffs.set_peak_gain_db(60.0 * value);
    }

    #[inline(always)]
    pub fn set_tone(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(0.0..=1.0, value);
        }

        self.lp1_coeffs
            .set_cutoff(475.0 + (20e3 - 475.0) * value * value * value);
    }

    #[inline(always)]
    pub fn set_volume(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        {
            debug_assert_is_finite(value);
            debug_assert_range(0.0..=1.0, value);
        }

        self.gain_coeffs.set_gain_lin(value * value * value);
    }

    #[inline(always)]
    pub fn coeffs_is_valid(&mut self) -> bool {
        todo!()
    }

    #[inline(always)]
    pub fn state_is_valid(&mut self) -> bool {
        todo!()
    }
}

impl<const N_CHANNELS: usize> Default for DistCoeffs<N_CHANNELS> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Default, PartialEq, Clone, Copy)]
pub struct DistState {
    hp1_state: HP1State,
    peak_state: PeakState,
    clip_state: ClipState,
    satur_state: SaturState,
    lp1_state: LP1State,
}

#[cfg(test)]
pub(crate) mod tests {
    use core::f32;

    use super::*;
    use crate::{
        c_wrapper::{bw_dist_coeffs, bw_dist_state, dist::Dist as DistWrapper},
        native::{
            clip::tests::{assert_clip_coeffs, assert_clip_state},
            gain::tests::assert_gain_coeffs,
            hp1::tests::{assert_hp1_coeffs, assert_hp1_state},
            lp1::tests::{assert_lp1_coeffs, assert_lp1_state},
            peak::tests::{assert_peak_coeffs, assert_peak_state},
            satur::tests::{assert_satur_coeffs, assert_satur_state},
        },
    };

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [&[1.0, 1.0], &[0.0, 0.0]];
    const N_SAMPLES: usize = 2;

    type DistT = Dist<N_CHANNELS>;
    type DistWrapperT = DistWrapper<N_CHANNELS>;

    #[test]
    fn new() {
        let rust_dist = DistT::new();
        let c_dist = DistWrapperT::new();

        assert_dist(&rust_dist, &c_dist);
    }

    #[test]
    fn set_sample_rate() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        assert_dist(&rust_dist, &c_dist);
    }

    #[test]
    fn reset_none() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        rust_dist.reset(None, None);
        c_dist.reset(None, None);

        assert_dist(&rust_dist, &c_dist);
    }

    #[test]
    fn reset_some() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let x0 = 0.3;
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_dist.reset(Some(x0), Some(&mut rust_y0));
        c_dist.reset(Some(x0), Some(&mut c_y0));

        assert_dist(&rust_dist, &c_dist);
        assert_eq!(rust_y0, c_y0);
    }

    #[test]
    fn reset_multi() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let x0 = [0.045, 0.045];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_dist.reset_multi(&x0, Some(&mut rust_y0));
        c_dist.reset_multi(&x0, Some(&mut c_y0));

        assert_dist(&rust_dist, &c_dist);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let distortion = 0.4;
        let volume = 0.9;
        let tone = 1.0;

        let mut rust_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];
        let mut c_y: [&mut [f32]; N_CHANNELS] = [&mut [0.0; N_SAMPLES], &mut [0.0; N_SAMPLES]];

        rust_dist.set_distortion(distortion);
        rust_dist.set_volume(volume);
        rust_dist.set_tone(tone);
        c_dist.set_distortion(distortion);
        c_dist.set_volume(volume);
        c_dist.set_tone(tone);

        rust_dist.reset(None, None);
        c_dist.reset(None, None);
        rust_dist.process(&PULSE_INPUT, &mut rust_y, N_SAMPLES);
        c_dist.process(&PULSE_INPUT, &mut c_y, N_SAMPLES);

        assert_dist(&rust_dist, &c_dist);
        assert_eq!(rust_y, c_y);
    }

    #[test]
    fn set_distortion() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let distortion = 1.0;
        rust_dist.set_distortion(distortion);
        c_dist.set_distortion(distortion);

        assert_dist(&rust_dist, &c_dist);
    }

    #[should_panic(expected = "value must be in range [0e0, 1e0], got -1e0")]
    #[test]
    fn set_distortion_invalid() {
        let mut rust_dist = DistT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);

        rust_dist.set_distortion(-1.0);
    }

    #[test]
    fn set_tone() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let tone = 0.5;

        rust_dist.set_tone(tone);
        c_dist.set_tone(tone);

        assert_dist(&rust_dist, &c_dist);
    }

    #[should_panic(expected = "value must be in range [0e0, 1e0], got 1.1e0")]
    #[test]
    fn set_tone_invalid() {
        let mut rust_dist = DistT::new();
        rust_dist.set_sample_rate(SAMPLE_RATE);

        rust_dist.set_tone(1.1);
    }

    #[test]
    fn set_volume() {
        let mut rust_dist = DistT::new();
        let mut c_dist = DistWrapperT::new();

        rust_dist.set_sample_rate(SAMPLE_RATE);
        c_dist.set_sample_rate(SAMPLE_RATE);

        let volume = 0.1234;

        rust_dist.set_volume(volume);
        c_dist.set_volume(volume);

        assert_dist(&rust_dist, &c_dist);
    }

    #[should_panic(expected = "value must be in range [0e0, 1e0], got 2e0")]
    #[test]
    fn set_volume_invalid() {
        let mut rust_dist = DistT::new();
        rust_dist.set_sample_rate(SAMPLE_RATE);

        rust_dist.set_volume(2.0);
    }

    fn assert_dist<const N_CHANNELS: usize>(
        rust_dist: &Dist<N_CHANNELS>,
        c_dist: &DistWrapper<N_CHANNELS>,
    ) {
        assert_dist_coeffs(&rust_dist.coeffs, &c_dist.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_dist_state(&rust_dist.states[channel], &c_dist.states[channel]);
        });
    }

    fn assert_dist_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &DistCoeffs<N_CHANNELS>,
        c_coeffs: &bw_dist_coeffs,
    ) {
        assert_hp1_coeffs(&rust_coeffs.hp1_coeffs, &c_coeffs.hp1_coeffs);
        assert_peak_coeffs(&rust_coeffs.peak_coeffs, &c_coeffs.peak_coeffs);
        assert_clip_coeffs(&rust_coeffs.clip_coeffs, &c_coeffs.clip_coeffs);
        assert_satur_coeffs(&rust_coeffs.satur_coeffs, &c_coeffs.satur_coeffs);
        assert_lp1_coeffs(&rust_coeffs.lp1_coeffs, &c_coeffs.lp1_coeffs);
        assert_gain_coeffs(&rust_coeffs.gain_coeffs, &c_coeffs.gain_coeffs);
    }

    fn assert_dist_state(rust_state: &DistState, c_state: &bw_dist_state) {
        assert_hp1_state(&rust_state.hp1_state, &c_state.hp1_state);
        assert_peak_state(&rust_state.peak_state, &c_state.peak_state);
        assert_clip_state(&rust_state.clip_state, &c_state.clip_state);
        assert_satur_state(&rust_state.satur_state, &c_state.satur_state);
        assert_lp1_state(&rust_state.lp1_state, &c_state.lp1_state);
    }
}
