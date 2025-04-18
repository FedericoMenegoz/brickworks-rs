use super::*;

pub(crate) struct OnePoleWrapper<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_one_pole_coeffs,
    pub(crate) states: [bw_one_pole_state; N_CHANNELS],
    pub(crate) states_p: [bw_one_pole_state; N_CHANNELS], // BW_RESTRICT to check what is for
}

impl<const N_CHANNELS: usize> OnePoleWrapper<N_CHANNELS> {
    pub(crate) fn new() -> Self {
        let states: [bw_one_pole_state; N_CHANNELS] = std::array::from_fn(|_| {
            bw_one_pole_state { y_z1: 0.0 }
        });

        let states_p: [bw_one_pole_state; N_CHANNELS] = std::array::from_fn(|_| {
            bw_one_pole_state { y_z1: 0.0 }
        });

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

    pub(crate) fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_one_pole_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }

    pub(crate) fn reset(&mut self, x0: Option<&[f32]>, mut y0: Option<&mut [f32]>) {
        unsafe {
            bw_one_pole_reset_coeffs(&mut self.coeffs);
            (0..N_CHANNELS).for_each(|i| {
                let x0_value = x0.map_or(0.0, |slice| slice[i]);
                let result = bw_one_pole_reset_state(&mut self.coeffs, &mut self.states[i], x0_value);
                if let Some(slice) = y0.as_deref_mut() {
                    slice[i] = result
                }
            });
        }
    }

    pub(crate) fn process(&mut self, x: &[&[f32]], y: &mut [&mut [f32]], n_samples: usize) {
        unsafe {
            let mut x_ptrs: [*const f32; N_CHANNELS] = [std::ptr::null(); N_CHANNELS];
            let mut y_ptrs: [*mut f32; N_CHANNELS] = [std::ptr::null_mut(); N_CHANNELS];
            let mut state_ptrs: [*mut bw_one_pole_state; N_CHANNELS] = [std::ptr::null_mut(); N_CHANNELS];
            (0..N_CHANNELS).for_each(|i| {
                x_ptrs[i] = x[i].as_ptr();
                y_ptrs[i] = y[i].as_mut_ptr();
                state_ptrs[i] = &mut self.states[i]
            });

            bw_one_pole_process_multi(
                &mut self.coeffs,
                state_ptrs.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }

    pub(crate) fn set_cutoff(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}!");
        unsafe {
            bw_one_pole_set_cutoff(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_cutoff_up(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}!");
        unsafe {
            bw_one_pole_set_cutoff_up(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_cutoff_down(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}!");
        unsafe {
            bw_one_pole_set_cutoff_down(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}!");
        unsafe {
            bw_one_pole_set_tau(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau_up(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}!");
        unsafe {
            bw_one_pole_set_tau_up(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau_down(&mut self, value: f32) {
        assert!(value >= 0., "Value must be non negative, got {value}.");
        unsafe {
            bw_one_pole_set_tau_down(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_sticky_thresh(&mut self, value: f32) {
        assert!(
            value >= 0. && value <= 1.0e18,
            "Value must be in range [0.0, 1.0e18], got {value}."
        );
        unsafe {
            bw_one_pole_set_sticky_thresh(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_sticky_mode(&mut self, value: bw_one_pole_sticky_mode) {
        assert!(
            value == bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs
                || value == bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel,
            "Sticky mode can be {bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs} or {bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel}, got {value}"
        );
        unsafe {
            bw_one_pole_set_sticky_mode(&mut self.coeffs, value as u32);
        }
    }

    pub(crate) fn get_sticky_thresh(&self) -> f32 {
        unsafe { bw_one_pole_get_sticky_thresh(&self.coeffs) }
    }

    pub(crate) fn get_sticky_mode(&self) -> bw_one_pole_sticky_mode {
        unsafe { bw_one_pole_get_sticky_mode(&self.coeffs) }
    }

    pub(crate) fn get_yz1(&self, channel: usize) -> f32 {
        unsafe { bw_one_pole_get_y_z1(&self.states[channel]) }
    }
}

#[cfg(test)]
mod tests {
    use core::f32;
    use std::f32::consts::PI;

    use super::*;
    const N: usize = 2;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);
    const BW_RCPF_ERROR: f32 = 0.0013;
    const SAMPLE_RATE: f32 = 48_000.0;

    #[test]
    fn new() {
        let one_pole = OnePoleWrapper::<N>::new();

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
        let mut f = OnePoleWrapper::<N>::new();

        f.set_sample_rate(SAMPLE_RATE);

        assert_eq!(f.coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE);
    }

    #[test]
    fn set_cutoff() {
        const CUT_OFF: f32 = 1000.0;
        let mut f = OnePoleWrapper::<N>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff(CUT_OFF);

        assert_eq!(f.coeffs.cutoff_up, CUT_OFF);
        assert_eq!(f.coeffs.cutoff_down, CUT_OFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    #[should_panic]
    fn set_cutoff_negative() {
        const CUT_OFF: f32 = -1000.0;
        let mut f = OnePoleWrapper::<N>::new();

        f.set_cutoff(CUT_OFF);
    }

    #[test]
    fn set_cutoff_up() {
        const CUT_OFF: f32 = 1200.0;
        let mut f = OnePoleWrapper::<N>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_up(CUT_OFF);

        assert_eq!(f.coeffs.cutoff_up, CUT_OFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_cutoff_down() {
        const CUT_OFF: f32 = 1200.0;
        let mut f = OnePoleWrapper::<N>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_down(CUT_OFF);

        assert_eq!(f.coeffs.cutoff_down, CUT_OFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau() {
        let cutoff = 150.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePoleWrapper::<N>::new();
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
        let mut f = OnePoleWrapper::<N>::new();
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
        let mut f = OnePoleWrapper::<N>::new();
        f.coeffs.param_changed = 0;

        f.set_tau_down(tau);

        let relative_error = (f.coeffs.cutoff_down - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        let tau = 0.1e-9;
        let mut f = OnePoleWrapper::<N>::new();
        f.coeffs.param_changed = 0;

        f.set_tau(tau);

        assert_eq!(f.coeffs.cutoff_up, f32::INFINITY);
        assert_eq!(f.coeffs.cutoff_down, f32::INFINITY);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN == 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP == 0);
    }

    #[test]
    #[should_panic]
    fn set_negative_tau() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_tau(-1.);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut f = OnePoleWrapper::<N>::new();

        f.set_sticky_thresh(sticky_tresh);

        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_STICKY_THRESH != 0);
        assert_eq!(f.get_sticky_thresh(), sticky_tresh);

    }

    #[test]
    #[should_panic]
    fn set_sticky_tresh_negative() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_sticky_thresh(-1.);
    }

    #[test]
    #[should_panic]
    fn set_sticky_tresh_too_high() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_sticky_thresh(1.1e18);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs;
        let mut f = OnePoleWrapper::<N>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[test]
    fn set_sticky_mode_rel() {
        let mode = bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel;
        let mut f = OnePoleWrapper::<N>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[test]
    #[should_panic]
    fn set_sticky_mode_not_valid() {
        let mode = 3;
        let mut f = OnePoleWrapper::<N>::new();

        f.set_sticky_mode(mode);
    }

    #[test]
    fn reset_none() {
        let cutoff = 1200.0;
        let sticky_tresh = 0.1;
        let mut f = OnePoleWrapper::<N>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_tresh);
        f.set_sample_rate(SAMPLE_RATE);
        let inverse= unsafe {
             f.coeffs.fs_2pi * bw_rcpf(f.coeffs.fs_2pi + cutoff)
        };

        f.reset(None, None);

        assert!(f.coeffs.param_changed == 0);
        assert_eq!(f.coeffs.mA1u, inverse);
        assert_eq!(f.coeffs.mA1d, inverse);
        assert_eq!(f.coeffs.st2, sticky_tresh * sticky_tresh);
    }

    #[test]
    fn reset_with_input_and_output() {
        let cutoff = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N];
        let mut y0_output = [0.0; N];
        let mut f = OnePoleWrapper::<N>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_thresh);
        f.set_sample_rate(SAMPLE_RATE);

        f.reset(Some(&x0_input), Some(&mut y0_output));

        for i in 0..N {
            assert_eq!(f.states[i].y_z1, x0_input[i]);
            assert_eq!(y0_output[i], x0_input[i]);
        }
    }

    #[test]
    fn process() {
        const N_SAMPLES: usize = 16;
        let mut f = OnePoleWrapper::<N>::new();
        let x = vec![vec![0.0; N_SAMPLES]; N];
        let mut y = vec![vec![0.0; N_SAMPLES]; N];
        let x_refs: Vec<&[f32]> = x.iter().map(|v| v.as_slice()).collect();
        let mut y_refs: Vec<&mut [f32]> = y.iter_mut().map(|v| v.as_mut_slice()).collect();

        f.process(&x_refs, &mut y_refs, N_SAMPLES);
    }

    #[test]
    fn test_get_yz1() {
        let f = OnePoleWrapper::<N>::new();
        for i in 0..N {
            let _ = f.get_yz1(i);
        }
    }
}
