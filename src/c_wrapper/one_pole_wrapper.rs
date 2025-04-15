use super::*;

pub(crate) struct OnePoleWrapper<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_one_pole_coeffs,
    pub(crate) states: [*mut bw_one_pole_state; N_CHANNELS],
    pub(crate) states_p: [*mut bw_one_pole_state; N_CHANNELS], // BW_RESTRICT to check what is for
}

impl<const N_CHANNELS: usize> OnePoleWrapper<N_CHANNELS> {
    pub(crate) fn new() -> Self {

        let states: [*mut bw_one_pole_state; N_CHANNELS] = std::array::from_fn(|_| {
            let state = Box::new(bw_one_pole_state { y_z1: 0.0 });
            Box::into_raw(state)
        });

        let states_p: [*mut bw_one_pole_state; N_CHANNELS] = std::array::from_fn(|_| {
            let state = Box::new(bw_one_pole_state { y_z1: 0.0 });
            Box::into_raw(state)
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

    pub(crate) fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32]>) {
        unsafe {
            bw_one_pole_reset_coeffs(&mut self.coeffs);
            if let Some(y0_value) = y0 {
                (0..N_CHANNELS).for_each(|i| {
                    y0_value[i] = bw_one_pole_reset_state(
                        &mut self.coeffs,
                    self.states[i],
                        x0.unwrap_or(0.),
                    )
                });
            } else {
                (0..N_CHANNELS).for_each(|i| {
                    bw_one_pole_reset_state(
                        &mut self.coeffs,
                        self.states[i],
                        x0.unwrap_or(0.),
                    );
                });
            }
        }
    }

    pub(crate) fn process(&mut self, x: &[&[f32]], y: &mut [&mut [f32]], n_samples: usize) {
        unsafe {
            
            let mut x_ptrs: [*const f32; N_CHANNELS] = [std::ptr::null(); N_CHANNELS];
            let mut y_ptrs: [*mut f32; N_CHANNELS] = [std::ptr::null_mut(); N_CHANNELS];
            (0..N_CHANNELS).for_each(|i| {
                x_ptrs[i] = x[i].as_ptr();
                y_ptrs[i] = y[i].as_mut_ptr();
            });

            bw_one_pole_process_multi(
                &mut self.coeffs,
                self.states.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }

    pub(crate) fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_cutoff_up(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff_up(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_cutoff_down(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff_down(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau_up(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau_up(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_tau_down(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau_down(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_sticky_thresh(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_sticky_thresh(&mut self.coeffs, value);
        }
    }

    pub(crate) fn set_sticky_mode(&mut self, value: bw_one_pole_sticky_mode) {
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
        unsafe { bw_one_pole_get_y_z1(self.states[channel]) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 2;

    #[test]
    fn new() {
        let _ = OnePoleWrapper::<N>::new();
    }

    #[test]
    fn set_sample_rate() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_sample_rate(48000.0);
    }

    #[test]
    fn set_cutoff() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_cutoff(1000.0);
    }

    #[test]
    fn set_cutoff_up() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_cutoff_up(1200.0);
    }

    #[test]
    fn set_cutoff_down() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_cutoff_down(800.0);
    }

    #[test]
    fn set_tau() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_tau(0.1);
    }

    #[test]
    fn set_tau_up() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_tau_up(0.15);
    }

    #[test]
    fn set_tau_down() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_tau_down(0.05);
    }

    #[test]
    fn set_sticky_thresh() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_sticky_thresh(0.2);
    }

    #[test]
    fn set_sticky_mode() {
        let mut f = OnePoleWrapper::<N>::new();
        f.set_sticky_mode(0);
    }

    #[test]
    fn get_sticky_thresh() {
        let f = OnePoleWrapper::<N>::new();
        let _ = f.get_sticky_thresh();
    }

    #[test]
    fn get_sticky_mode() {
        let f = OnePoleWrapper::<N>::new();
        let _ = f.get_sticky_mode();
    }

    #[test]
    fn reset_none() {
        let mut f = OnePoleWrapper::<N>::new();
        f.reset(Some(0.0), None);
    }

    #[test]
    fn reset_some() {
        let mut f = OnePoleWrapper::<N>::new();
        let mut y = vec![0.0; N];
        f.reset(Some(0.0), Some(&mut y));
    }

    #[test]
    fn process() {
        let mut f = OnePoleWrapper::<N>::new();
        let x = vec![vec![0.0; 16]; N];
        let mut y = vec![vec![0.0; 16]; N];
        let x_refs: Vec<&[f32]> = x.iter().map(|v| v.as_slice()).collect();
        let mut y_refs: Vec<&mut [f32]> = y.iter_mut().map(|v| v.as_mut_slice()).collect();
        
        f.process(&x_refs, &mut y_refs, 16);
    }

    #[test]
    fn test_get_yz1() {
        let f = OnePoleWrapper::<N>::new();
        for i in 0..N {
            let _ = f.get_yz1(i);
        }
    }
}
