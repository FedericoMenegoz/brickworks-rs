use crate::global::{assert_positive, INVERSE_2_PI};
use bitflags::bitflags;

#[allow(dead_code)]
struct OnePole<const N_CHANNELS: usize> {
    coeffs: OnePoleCoeffs,
    states: Vec<OnePoleState>,
    states_p: Vec<OnePoleState>, // BW_RESTRICT to check what is for
}

#[allow(dead_code, unused_mut, unused_variables)]
#[derive(Clone, Debug, Copy)]
struct OnePoleCoeffs {
    pub fs_2pi: f32,
    pub m_a1u: f32,
    pub m_a1d: f32,
    pub st2: f32,
    pub cutoff_up: f32,
    pub cutoff_down: f32,
    pub sticky_thresh: f32,
    pub sticky_mode: OnePoleStickyMode,
    pub param_changed: ParamChanged,
}

#[allow(dead_code, unused_mut, unused_variables)]
#[derive(Debug, PartialEq, Clone, Copy)]
enum OnePoleStickyMode {
    Abs,
    Rel,
}

#[allow(dead_code, unused_mut, unused_variables)]
#[derive(Debug, Clone, Copy)]
struct OnePoleState {
    y_z1: f32,
}

bitflags! {
    #[derive(Clone, Debug, Copy)]
    struct ParamChanged: u32 {
        const CUTOFF_UP = 1;
        const CUTOFF_DOWN = 1<<1;
        const STICKY_TRESH = 1<<2;
    }
}

#[allow(dead_code, unused_mut, unused_variables)]
impl<const N_CHANNELS: usize> OnePole<N_CHANNELS> {
    pub fn new() -> Self {
        OnePole {
            coeffs: Default::default(),
            states: vec![OnePoleState { y_z1: 0.0 }; N_CHANNELS],
            states_p: vec![OnePoleState { y_z1: 0.0 }; N_CHANNELS],
        }
    }

    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        self.coeffs.fs_2pi = INVERSE_2_PI * sample_rate;
    }

    pub fn reset(&mut self, x0: Option<&[f32]>, mut y0: Option<&mut [f32]>) {
        todo!()
    }

    pub fn process(&mut self, x: &[Vec<f32>], y: Option<&mut [&mut [f32]]>, n_samples: usize) {
        todo!()
    }

    pub fn set_cutoff(&mut self, value: f32) {
        assert_positive(value);
        self.set_cutoff_up(value);
        self.set_cutoff_down(value);
    }

    pub fn set_cutoff_up(&mut self, value: f32) {
        assert_positive(value);
        if self.coeffs.cutoff_up != value {
            self.coeffs.cutoff_up = value;
            self.coeffs.param_changed |= ParamChanged::CUTOFF_UP;
        }
    }

    pub fn set_cutoff_down(&mut self, value: f32) {
        assert_positive(value);
        if self.coeffs.cutoff_down != value {
            self.coeffs.cutoff_down = value;
            self.coeffs.param_changed |= ParamChanged::CUTOFF_DOWN;
        }
    }

    pub fn set_tau(&mut self, value: f32) {
        todo!()
    }

    pub fn set_tau_up(&mut self, value: f32) {
        todo!()
    }

    pub fn set_tau_down(&mut self, value: f32) {
        todo!()
    }

    pub fn set_sticky_thresh(&mut self, value: f32) {
        todo!()
    }

    pub fn set_sticky_mode(&mut self, value: OnePoleStickyMode) {
        todo!()
    }

    pub fn get_sticky_thresh(&self) -> f32 {
        todo!()
    }

    pub fn get_sticky_mode(&self) -> OnePoleStickyMode {
        todo!()
    }

    pub fn get_yz1(&self, channel: usize) -> f32 {
        todo!()
    }
}

impl Default for OnePoleCoeffs {
    fn default() -> Self {
        Self {
            fs_2pi: Default::default(),
            m_a1u: Default::default(),
            m_a1d: Default::default(),
            st2: Default::default(),
            cutoff_up: f32::INFINITY,
            cutoff_down: f32::INFINITY,
            sticky_thresh: 0.0,
            sticky_mode: OnePoleStickyMode::Abs,
            param_changed: ParamChanged::all(),
        }
    }
}

#[allow(unused_mut)]
#[allow(unused_assignments)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::c_wrapper::{one_pole_wrapper::*, *};

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;

    #[test]
    fn one_pole_initialization() {
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_sample_rate() {
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        assert_eq!(rust_one_pole.coeffs.fs_2pi, c_one_pole.coeffs.fs_2pi)
    }

    #[test]
    fn set_cutoff() {
        const CUTOFF: f32 = 1000.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_cutoff(CUTOFF);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    #[should_panic(expected = "Value must be non negative, got -1!")]
    fn set_cutoff_negative() {
        const CUTOFF: f32 = -1.0;
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_cutoff(CUTOFF);
    }

    #[test]
    fn set_cutoff_up() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_up(CUTOFF);
        rust_one_pole.set_cutoff_up(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_up, CUTOFF);
        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_cutoff_down() {
        const CUTOFF: f32 = 1200.0;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_cutoff_down(CUTOFF);
        rust_one_pole.set_cutoff_down(CUTOFF);

        assert_eq!(rust_one_pole.coeffs.cutoff_down, CUTOFF);
        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau() {
        const CUTOFF: f32 = 150.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_up() {
        const CUTOFF: f32 = 1500.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;

        c_one_pole.set_tau_up(TAU);
        rust_one_pole.set_tau_up(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_down() {
        const CUTOFF: f32 = 10_000.0;
        const TAU: f32 = INVERSE_2_PI / CUTOFF;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau_down(TAU);
        rust_one_pole.set_tau_down(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        const TAU: f32 = 0.1e-9;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.coeffs.param_changed = 0;
        rust_one_pole.coeffs.param_changed = ParamChanged::empty();

        c_one_pole.set_tau(TAU);
        rust_one_pole.set_tau(TAU);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    #[should_panic(expected = "Value must be non negative, got -1!")]
    fn set_negative_tau() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_tau(-1.);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sticky_thresh(sticky_tresh);
        rust_one_pole.set_sticky_thresh(sticky_tresh);

        assert!(rust_one_pole.coeffs.param_changed.bits() & BW_ONE_POLE_PARAM_STICKY_THRESH != 0);
        assert_eq!(rust_one_pole.get_sticky_thresh(), sticky_tresh);
        assert_eq!(
            rust_one_pole.get_sticky_thresh(),
            c_one_pole.get_sticky_thresh()
        );
    }

    #[test]
    #[should_panic(expected = "Value must be in range [0e0, 1e18], got -1e0!")]
    fn set_sticky_tresh_negative() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_sticky_thresh(-1.);
    }

    #[test]
    #[should_panic(expected = "Value must be in range [0e0, 1e18], got 1.1e18!")]
    fn set_sticky_tresh_too_high() {
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        rust_one_pole.set_sticky_thresh(1.1e18);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let c_mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs;
        let rust_mode = OnePoleStickyMode::Abs;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sticky_mode(c_mode);
        rust_one_pole.set_sticky_mode(rust_mode);

        assert_eq!(
            rust_one_pole.get_sticky_mode() as u32,
            c_one_pole.get_sticky_mode()
        );
        assert_eq!(rust_one_pole.get_sticky_mode(), rust_mode);
    }

    #[test]
    fn set_sticky_mode_rel() {
        let c_mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel;
        let rust_mode = OnePoleStickyMode::Rel;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sticky_mode(c_mode);
        rust_one_pole.set_sticky_mode(rust_mode);

        assert_eq!(
            rust_one_pole.get_sticky_mode() as u32,
            c_one_pole.get_sticky_mode()
        );
        assert_eq!(rust_one_pole.get_sticky_mode(), rust_mode);
    }

    #[test]
    fn reset_none() {
        const CUTOFF: f32 = 1200.0;
        let sticky_tresh = 0.1;
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_tresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_tresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(None, None);
        rust_one_pole.reset(None, None);

        assert_coeffs_rust_c(rust_one_pole.coeffs, c_one_pole.coeffs);
    }

    #[test]
    fn reset_with_input_and_output() {
        const CUTOFF: f32 = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N_CHANNELS];
        let mut c_y0_output = [0.0; N_CHANNELS];
        let mut rust_y0_output = [0.0; N_CHANNELS];
        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_thresh(sticky_thresh);
        c_one_pole.set_sample_rate(SAMPLE_RATE);

        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_thresh(sticky_thresh);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);

        c_one_pole.reset(Some(&x0_input), Some(&mut c_y0_output));
        rust_one_pole.reset(Some(&x0_input), Some(&mut rust_y0_output));

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.states[i].y_z1, x0_input[i]);
            assert_eq!(rust_one_pole.states[i].y_z1, c_one_pole.states[i].y_z1);
            assert_eq!(rust_y0_output[i], x0_input[i]);
        }
    }

    #[test]
    fn process_with_y() {
        const N_CHANNELS: usize = 2;
        const N_SAMPLES: usize = 4;
        const CUTOFF: f32 = 1000.0;

        let input = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut c_output: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0]];
        let mut rust_output: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0]];

        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.reset(None, None);

        c_one_pole.process(&input, Some(&mut c_output), N_SAMPLES);

        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.reset(None, None);

        rust_one_pole.process(&input, Some(&mut rust_output), N_SAMPLES);

        for ch in 0..N_CHANNELS {
            for sample in 0..N_SAMPLES {
                assert_eq!(rust_output[ch][sample], c_output[ch][sample]);
                println!(
                    "C output: {}\nRust output: {}",
                    c_output[ch][sample], rust_output[ch][sample]
                )
            }
        }
    }

    #[test]
    fn get_yz1() {
        const N_CHANNELS: usize = 2;
        const SAMPLE_RATE: f32 = 44100.0;
        const CUTOFF: f32 = 1000.0;
        const N_SAMPLES: usize = 4;

        let mut c_one_pole = OnePoleWrapper::<N_CHANNELS>::new();
        let mut rust_one_pole = OnePole::<N_CHANNELS>::new();

        c_one_pole.set_sample_rate(SAMPLE_RATE);
        c_one_pole.set_cutoff(CUTOFF);
        c_one_pole.set_sticky_mode(bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel);

        rust_one_pole.set_sample_rate(SAMPLE_RATE);
        rust_one_pole.set_cutoff(CUTOFF);
        rust_one_pole.set_sticky_mode(OnePoleStickyMode::Rel);

        let input = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut c_output: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.1, 0.2, 0.3], &mut [1.0, 1.1, 1.2, 1.3]];
        let mut rust_output: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.1, 0.2, 0.3], &mut [1.0, 1.1, 1.2, 1.3]];

        c_one_pole.process(&input, Some(&mut c_output), N_SAMPLES);
        rust_one_pole.process(&input, Some(&mut rust_output), N_SAMPLES);

        for i in 0..N_CHANNELS {
            assert_eq!(rust_one_pole.get_yz1(i), c_one_pole.get_yz1(i));
            assert_eq!(rust_one_pole.get_yz1(i), rust_output[i][3]);
        }
    }

    fn assert_coeffs_rust_c(rust_coeffs: OnePoleCoeffs, c_coeffs: bw_one_pole_coeffs) {
        assert_eq!(rust_coeffs.fs_2pi, c_coeffs.fs_2pi);
        assert_eq!(rust_coeffs.m_a1u, c_coeffs.mA1u);
        assert_eq!(rust_coeffs.m_a1d, c_coeffs.mA1d);
        assert_eq!(rust_coeffs.st2, c_coeffs.st2);
        assert_eq!(rust_coeffs.cutoff_up, c_coeffs.cutoff_up);
        assert_eq!(rust_coeffs.cutoff_down, c_coeffs.cutoff_down);
        assert_eq!(rust_coeffs.sticky_thresh, c_coeffs.sticky_thresh);
        assert_eq!(rust_coeffs.sticky_mode as u32, c_coeffs.sticky_mode);
        assert_eq!(rust_coeffs.param_changed.bits(), c_coeffs.param_changed as u32);
    }
}
