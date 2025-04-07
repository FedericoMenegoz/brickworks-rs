struct OnePole<const N_CHANNELS: usize> {
    coeffs: OnePoleCoeffs,
    states: [OnePoleState; N_CHANNELS],
    states_p: [OnePoleState; N_CHANNELS], // BW_RESTRICT to check what is for
}

impl<const N_CHANNELS: usize> OnePole<N_CHANNELS> {
    pub fn new() -> Self {
        todo!()
    }

    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        todo!()
    }

    pub fn reset(&self, x0: Option<f32>, y0: Option<f32>) {
        todo!()
    }

    pub fn process(x: &[&[f32]], y: &[&mut [f32]], n_samples: usize) {
        todo!()
    }

    pub fn set_cutoff(value: f32) {
        todo!()
    }

    pub fn set_cutoff_up(value: f32) {
        todo!()
    }

    pub fn set_cutoff_down(value: f32) {
        todo!()
    }

    pub fn set_tau(value: f32) {
        todo!()
    }

    pub fn set_tau_up(value: f32) {
        todo!()
    }

    pub fn set_tau_down(value: f32) {
        todo!()
    }

    pub fn set_sticky_thresh(value: f32) {
        todo!()
    }

    pub fn set_sticky_mode(value: OnePoleStickyMode) {
        todo!()
    }

    pub fn get_sticky_thresh() -> f32 {
        todo!()
    }

    pub fn get_sticky_mode() -> OnePoleStickyMode {
        todo!()
    }

    pub fn get_yz1(channel: usize) {
        todo!()
    }
}
struct OnePoleCoeffs {
    pub fs_2pi: f32,
    pub m_a1u: f32,
    pub m_a1d: f32,
    pub st2: f32,
    pub cutoff_up: f32,
    pub cutoff_down: f32,
    pub sticky_thresh: f32,
    pub sticky_mode: u32,
    pub param_changed: u32,
}
enum OnePoleStickyMode {
    Abs,
    Rel,
}
struct OnePoleState {
    y_z1: f32,
}

#[allow(unused_mut)]
#[allow(unused_assignments)]
#[cfg(test)]
mod test {
    use std::panic;

    use super::*;
    use crate::c_wrapper::*;

    #[test]
    fn one_pole_initialization() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let rust_one_pole: OnePole<1>;

        // Act
        unsafe { bw_one_pole_init(&mut c_coeffs) }
        rust_one_pole = OnePole::<1>::new();

        // Assert
        assert_coeffs_rust_c(rust_one_pole.coeffs, c_coeffs);
    }

    #[test]
    fn one_pole_set_sample_rate_greater_than_zero() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let sample_rate: f32 = 100.;
        let mut rust_one_pole = OnePole::<1>::new();

        // Act
        unsafe {
            bw_one_pole_set_sample_rate(&mut c_coeffs, sample_rate);
        }
        rust_one_pole.set_sample_rate(sample_rate);

        // Assert
        assert_eq!(rust_one_pole.coeffs.fs_2pi, c_coeffs.fs_2pi);
    }

    #[test]
    fn one_pole_set_sample_rate_zero() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let sample_rate: f32 = 0.;

        // Act
        unsafe {
            let c_result = panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                bw_one_pole_set_sample_rate(&mut c_coeffs, sample_rate);
            }));

            // Initialize rust_one_pole here to avoid the same issue
            let mut rust_one_pole = OnePole::<1>::new();
            let rust_result = panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                rust_one_pole.set_sample_rate(sample_rate)
            }));

            //Assert
            assert!(
                c_result.is_err() && rust_result.is_err(),
                "Expected panic in both c and rust, got c = {:?} and rust = {:?}",
                c_result,
                rust_result
            );
        }
    }

    #[test]
    fn one_pole_do_update_coeffs_ctrl_cutoff_up() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        c_coeffs.param_changed = BW_ONE_POLE_PARAM_CUTOFF_UP as i32;
        let mut rust_one_pole = OnePole::<1>::new();
        rust_one_pole.coeffs.param_changed = BW_ONE_POLE_PARAM_CUTOFF_UP;

        // Act
        unsafe {
            bw_one_pole_do_update_coeffs_ctrl(&mut c_coeffs);
        }
        todo!();
        // rust_one_pole.do_update_coeffs_ctrl();

        // Assert
        // assert_coeffs_rust_c(rust_one_pole.coeffs, c_coeffs);
    }

    #[test]
    fn one_pole_do_update_coeffs_ctrl_cutoff_down() {
        todo!();
        // Arrange
        // Act
        // Assert
    }

    #[test]
    fn one_pole_do_update_coeffs_ctrl_sticky_tresh() {
        todo!();
        // Arrange
        // Act
        // Assert
    }

    #[test]
    fn one_pole_do_update_coeffs_ctrl_all() {
        todo!();
        // Arrange
        // Act
        // Assert
    }

    #[test]
    fn one_pole_reset_coeffs() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let mut rust_one_pole = OnePole::<1>::new();

        // Act
        unsafe {
            bw_one_pole_reset_coeffs(&mut c_coeffs);
        }
        todo!();
        // rust_one_pole.reset_coeffs();

        // // Assert
        // assert_coeffs_rust_c(rust_one_pole.coeffs, c_coeffs);
    }

    #[test]
    fn one_pole_reset_state() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let mut c_state = bw_one_pole_state { y_z1: 0.1 };
        let x_0 = 0.2;
        let mut rust_one_pole = OnePole::<1>::new();
        let c_return;
        // let rust_return;
        // Act
        unsafe {
            c_return = bw_one_pole_reset_state(&mut c_coeffs, &mut c_state, x_0);
        }
        todo!();
        // rust_return = rust_one_pole.reset_coeffs(x_0);

        // // Assert
        // assert_eq!(rust_one_pole.state.y_z1, c_state.y_z1);
        // assert_eq!(rust_return, c_return);
    }

    #[test]
    fn one_pole_reset_state_multi() {
        // Arrange
        let n_channels = 2;
        let mut c_states: Vec<bw_one_pole_state> = (0..n_channels)
            .map(|val| bw_one_pole_state {
                y_z1: val as f32 / 10.,
            })
            .collect();
        let mut c_state_ptrs: Vec<*mut bw_one_pole_state> =
            c_states.iter_mut().map(|s| s as *mut _).collect();
        let c_state_ptrs_ref = c_state_ptrs.as_ptr();
        let mut x_0 = (0..n_channels)
            .map(|val| val as f32 + 0.5)
            .collect::<Vec<f32>>();
        let mut y_0 = (0..n_channels)
            .map(|val| val as f32 + 0.3)
            .collect::<Vec<f32>>();
        let mut rust_one_pole = OnePole::<1>::new();

        // Act
        unsafe {
            bw_one_pole_reset_state_multi(
                &C_COEFFS,
                c_state_ptrs_ref,
                x_0.as_mut_ptr(),
                y_0.as_mut_ptr(),
                n_channels,
            );
        }
        // UNDERSTAND HOW TO DEAL WITH CHANNELS AND STATES ARE THEY PART OF THE ONE POLE?
        // rust_one_pole.reset_coeffs_multi(x_0, y_0);

        // Assert
        todo!()
    }

    // bw_one_pole_update_coeffs_ctrl and bw_one_pole_update_coeffs_audio are
    // only checking if the coeffs are valid
    // need to undestrand about these asserts in the c library
    #[test]
    fn bw_one_pole_update_coeffs_ctrl() {
        todo!()
    }
    #[test]
    fn bw_one_pole_update_coeffs_audio() {
        todo!()
    }

    #[test]
    fn one_pole_process1() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let mut c_state = bw_one_pole_state { y_z1: 0.1 };
        let x = 0.1;
        let mut rust_one_pole: OnePole<1>;
        let c_return;
        // let rust_return;

        // Act
        unsafe {
            c_return = bw_one_pole_process1(&mut c_coeffs, &mut c_state, x);
        }
        todo!();
        // rust_return = rust_one_pole.process1(x);

        // // Assert
        // assert_eq!(rust_return, c_return);
        // assert_eq!(rust_one_pole.state.y_z1, c_state.y_z1);
    }

    #[test]
    fn one_pole_process1_sticky_abs() {
        // Arrange
        let mut c_coeffs = bw_one_pole_coeffs::new();
        let mut c_state = bw_one_pole_state { y_z1: 0.1 };
        let x = 0.1;
        let mut rust_one_pole: OnePole<1>;
        let c_return;
        // let rust_return;

        // Act
        unsafe {
            c_return = bw_one_pole_process1(&mut c_coeffs, &mut c_state, x);
        }
        todo!();

        // rust_return = rust_one_pole.process1(x);

        // // Assert
        // assert_eq!(rust_return, c_return);
        // assert_eq!(rust_one_pole.state.y_z1, c_state.y_z1);
    }

    // Utility functions for testing
    impl bw_one_pole_coeffs {
        fn new() -> Self {
            bw_one_pole_coeffs {
                fs_2pi: 0.0,
                mA1u: 0.0,
                mA1d: 0.0,
                st2: 0.0,
                cutoff_up: 0.0,
                cutoff_down: 0.0,
                sticky_thresh: 0.0,
                sticky_mode: 0,
                param_changed: 0,
            }
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
        assert_eq!(rust_coeffs.sticky_mode, c_coeffs.sticky_mode);
        assert_eq!(rust_coeffs.param_changed, c_coeffs.param_changed as u32);
    }

    // Consts and variable for testing
    const C_COEFFS: bw_one_pole_coeffs = bw_one_pole_coeffs {
        fs_2pi: 0.0,
        mA1u: 0.0,
        mA1d: 0.0,
        st2: 0.0,
        cutoff_up: 0.0,
        cutoff_down: 0.0,
        sticky_thresh: 0.0,
        sticky_mode: 0,
        param_changed: 0,
    };
}
