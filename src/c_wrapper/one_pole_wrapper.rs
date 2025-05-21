use std::ptr::null_mut;

#[cfg(debug_assertions)]
use crate::global::{debug_assert_positive, debug_assert_range};

use super::*;

pub struct OnePoleWrapper<const N_CHANNELS: usize> {
    pub coeffs: bw_one_pole_coeffs,
    pub states: [bw_one_pole_state; N_CHANNELS],
    pub states_p: [bw_one_pole_state; N_CHANNELS], // BW_RESTRICT to check what is for
}

impl<const N_CHANNELS: usize> OnePoleWrapper<N_CHANNELS> {
    pub fn new() -> Self {
        let states: [bw_one_pole_state; N_CHANNELS] =
            std::array::from_fn(|_| bw_one_pole_state { y_z1: 0.0 });

        let states_p: [bw_one_pole_state; N_CHANNELS] =
            std::array::from_fn(|_| bw_one_pole_state { y_z1: 0.0 });

        let mut one_pole = OnePoleWrapper {
            coeffs: bw_one_pole_coeffs {
                fs_2pi: Default::default(),
                mA1u: Default::default(),
                mA1d: Default::default(),
                st2: Default::default(),
                cutoff_up: Default::default(),
                cutoff_down: Default::default(),
                sticky_thresh: Default::default(),
                sticky_mode: Default::default(),
                param_changed: Default::default(),
            },
            states,
            states_p,
        };
        unsafe {
            bw_one_pole_init(&mut one_pole.coeffs);
        }
        one_pole
    }

    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_one_pole_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    pub fn reset(&mut self, x0: &[f32], mut y0: Option<&mut [f32]>) {
        unsafe {
            bw_one_pole_reset_coeffs(&mut self.coeffs);
            (0..N_CHANNELS).for_each(|channel| {
                let result = bw_one_pole_reset_state(
                    &mut self.coeffs,
                    &mut self.states[channel],
                    x0[channel],
                );
                if let Some(slice) = y0.as_deref_mut() {
                    slice[channel] = result
                }
            });
        }
    }

    pub fn process(&mut self, x: &[Vec<f32>], y: Option<&mut [&mut [f32]]>, n_samples: usize) {
        unsafe {
            // In case y is None this will be passed
            let null_ptrs: [*mut f32; N_CHANNELS] = [null_mut(); N_CHANNELS];
            // Need to prepare the data into raw pointers for c
            let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
            let y_ptrs: Option<[*mut f32; N_CHANNELS]> =
                y.map(|state| std::array::from_fn(|i| state[i].as_mut_ptr()));
            let mut state_ptrs: [*mut bw_one_pole_state; N_CHANNELS] =
                std::array::from_fn(|i| &mut self.states[i] as *mut _);

            bw_one_pole_process_multi(
                &mut self.coeffs,
                state_ptrs.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.unwrap_or(null_ptrs).as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }

    pub fn set_cutoff(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_cutoff(&mut self.coeffs, value);
        }
    }

    pub fn set_cutoff_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_cutoff_up(&mut self.coeffs, value);
        }
    }

    pub fn set_cutoff_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_cutoff_down(&mut self.coeffs, value);
        }
    }

    pub fn set_tau(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_tau(&mut self.coeffs, value);
        }
    }

    pub fn set_tau_up(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_tau_up(&mut self.coeffs, value);
        }
    }

    pub fn set_tau_down(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_positive(value);
        unsafe {
            bw_one_pole_set_tau_down(&mut self.coeffs, value);
        }
    }

    pub fn set_sticky_thresh(&mut self, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert_range(0., 1.0e18, value);
        unsafe {
            bw_one_pole_set_sticky_thresh(&mut self.coeffs, value);
        }
    }

    pub fn set_sticky_mode(&mut self, value: bw_one_pole_sticky_mode) {
        debug_assert!(
            value == bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs
                || value == bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel,
            "sticky mode can be {bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs} or {bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel}, got {value}"
        );
        unsafe {
            bw_one_pole_set_sticky_mode(&mut self.coeffs, value as u32);
        }
    }

    pub fn get_sticky_thresh(&self) -> f32 {
        unsafe { bw_one_pole_get_sticky_thresh(&self.coeffs) }
    }

    pub fn get_sticky_mode(&self) -> bw_one_pole_sticky_mode {
        unsafe { bw_one_pole_get_sticky_mode(&self.coeffs) }
    }

    pub fn get_y_z1(&self, channel: usize) -> f32 {
        unsafe { bw_one_pole_get_y_z1(&self.states[channel]) }
    }

    pub(crate) fn process1(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1(&mut self.coeffs, &mut self.states[channel], x) }
    }

    pub(crate) fn process1_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_sticky_abs(&mut self.coeffs, &mut self.states[channel], x) }
    }

    pub(crate) fn process1_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_sticky_rel(&mut self.coeffs, &mut self.states[channel], x) }
    }
    pub(crate) fn process1_asym(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_asym(&mut self.coeffs, &mut self.states[channel], x) }
    }
    pub(crate) fn process1_asym_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        unsafe {
            bw_one_pole_process1_asym_sticky_abs(&mut self.coeffs, &mut self.states[channel], x)
        }
    }
    pub(crate) fn process1_asym_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        unsafe {
            bw_one_pole_process1_asym_sticky_rel(&mut self.coeffs, &mut self.states[channel], x)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::global::*;
    use core::f32;
    const N_CHANNELS: usize = 2;
    const BW_RCPF_ERROR: f32 = 0.0013;
    const SAMPLE_RATE: f32 = 48_000.0;

    #[test]
    fn new() {
        let one_pole = OnePoleWrapper::<N_CHANNELS>::new();

        assert_eq!(one_pole.coeffs.cutoff_up, f32::INFINITY);
        assert_eq!(one_pole.coeffs.cutoff_down, f32::INFINITY);
        assert_eq!(one_pole.coeffs.sticky_thresh, 0.);
        assert_eq!(
            one_pole.coeffs.sticky_mode,
            bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs
        );
        assert_eq!(one_pole.coeffs.param_changed, !0); // from brickworks lib coeffs->param_changed = ~0; // useless, just to make compilers happy about uninitialized variables
    }

    #[test]
    fn set_sample_rate() {
        const SAMPLE_RATE: f32 = 48_000.0;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_sample_rate(SAMPLE_RATE);

        assert_eq!(f.coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE)
    }

    #[test]
    fn set_cutoff() {
        const CUTOFF: f32 = 1000.0;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff(CUTOFF);

        assert_eq!(f.coeffs.cutoff_up, CUTOFF);
        assert_eq!(f.coeffs.cutoff_down, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1000")]
    fn set_cutoff_negative() {
        const CUTOFF: f32 = -1000.0;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_cutoff(CUTOFF);
    }

    #[test]
    fn set_cutoff_up() {
        const CUTOFF: f32 = 1200.0;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_up(CUTOFF);

        assert_eq!(f.coeffs.cutoff_up, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_cutoff_down() {
        const CUTOFF: f32 = 1200.0;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_down(CUTOFF);

        assert_eq!(f.coeffs.cutoff_down, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau() {
        let cutoff = 150.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau(tau);

        let relative_error = (f.coeffs.cutoff_down - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR); // see bw_rcpf in bw_math.h
        assert_eq!(f.coeffs.cutoff_down, f.coeffs.cutoff_up);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_tau_up() {
        let cutoff = 1500.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau_up(tau);

        let relative_error = (f.coeffs.cutoff_up - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_tau_down() {
        let cutoff = 10_000.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau_down(tau);

        let relative_error = (f.coeffs.cutoff_down - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        let tau = 0.1e-9;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau(tau);

        assert_eq!(f.coeffs.cutoff_up, f32::INFINITY);
        assert_eq!(f.coeffs.cutoff_down, f32::INFINITY);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN == 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP == 0);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be non negative, got -1")]
    fn set_negative_tau() {
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_tau(-1.);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_sticky_thresh(sticky_tresh);

        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_STICKY_THRESH != 0);
        assert_eq!(f.get_sticky_thresh(), sticky_tresh);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got -1e0")]
    fn set_sticky_tresh_negative() {
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_sticky_thresh(-1.);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "value must be in range [0e0, 1e18], got 1.1e18")]
    fn set_sticky_tresh_too_high() {
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_sticky_thresh(1.1e18);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[test]
    fn set_sticky_mode_rel() {
        let mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "sticky mode can be 0 or 1, got 3")]
    fn set_sticky_mode_not_valid() {
        let mode = 3;
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();

        f.set_sticky_mode(mode);
    }

    #[test]
    fn reset_none() {
        let cutoff = 1200.0;
        let sticky_tresh = 0.1;
        let x0_input = [0.5; N_CHANNELS];
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_tresh);
        f.set_sample_rate(SAMPLE_RATE);
        let inverse = unsafe { f.coeffs.fs_2pi * bw_rcpf(f.coeffs.fs_2pi + cutoff) };

        f.reset(&x0_input, None);

        assert!(f.coeffs.param_changed == 0);
        assert_eq!(f.coeffs.mA1u, inverse);
        assert_eq!(f.coeffs.mA1d, inverse);
        assert_eq!(f.coeffs.st2, sticky_tresh * sticky_tresh);
    }

    #[test]
    fn reset_with_input_and_output() {
        let cutoff = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N_CHANNELS];
        let mut y0_output = [0.0; N_CHANNELS];
        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_thresh);
        f.set_sample_rate(SAMPLE_RATE);

        f.reset(&x0_input, Some(&mut y0_output));

        for i in 0..N_CHANNELS {
            assert_eq!(f.states[i].y_z1, x0_input[i]);
            assert_eq!(y0_output[i], x0_input[i]);
        }
    }

    #[test]
    fn process_with_y() {
        const N_CHANNELS: usize = 2;
        const N_SAMPLES: usize = 4;
        let x0_input = [0.5; N_CHANNELS];

        let input_data = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut output_data: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0]];

        let mut filter = OnePoleWrapper::<N_CHANNELS>::new();
        filter.set_cutoff(1000.0);
        filter.set_sample_rate(44100.0);
        filter.reset(&x0_input, None);

        filter.process(&input_data, Some(&mut output_data), N_SAMPLES);

        for ch in 0..N_CHANNELS {
            for sample in 0..N_SAMPLES {
                println!("Output[{}][{}] = {}", ch, sample, output_data[ch][sample]);
                assert!(output_data[ch][sample].is_finite());
            }
        }
    }

    #[test]
    fn get_y_z1() {
        const N_CHANNELS: usize = 2;
        const SAMPLE_RATE: f32 = 44100.0;
        const CUTOFF: f32 = 1000.0;
        const N_SAMPLES: usize = 4;

        let mut f = OnePoleWrapper::<N_CHANNELS>::new();
        f.set_sample_rate(SAMPLE_RATE);
        f.set_cutoff(CUTOFF);
        f.set_sticky_mode(bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel);
        let input_data = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut output_data: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.1, 0.2, 0.3], &mut [1.0, 1.1, 1.2, 1.3]];
        f.process(&input_data, Some(&mut output_data), N_SAMPLES);

        for i in 0..N_CHANNELS {
            assert_eq!(f.get_y_z1(i), output_data[i][3])
        }
    }
}
