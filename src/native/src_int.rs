use std::f32::consts::PI;

use crate::native::math::tanf;

pub struct SRCInt<const N_CHANNELS: usize> {
    pub(crate) coeffs: SRCIntCoeffs<N_CHANNELS>,
    pub(crate) states: [SRCIntState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SRCInt<N_CHANNELS> {
    #[inline(always)]
    pub fn new(ratio: i32) -> Self {
        Self {
            coeffs: SRCIntCoeffs::<N_CHANNELS>::new(ratio),
            states: [SRCIntState::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
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
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        self.coeffs.reset_state_multi(&mut self.states, x0, y0);
    }

    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_in_samples: usize,
        n_out_samples: Option<&mut [usize; N_CHANNELS]>,
    ) {
        self.coeffs
            .process_multi(&mut self.states, x, y, n_in_samples, n_out_samples);
    }
}

pub struct SRCIntCoeffs<const N_CHANNELS: usize> {
    ratio: i32,
    b0: f32,
    ma1: f32,
    ma2: f32,
    ma3: f32,
    ma4: f32,
}

impl<const N_CHANNELS: usize> SRCIntCoeffs<N_CHANNELS> {
    #[inline(always)]
    pub fn new(ratio: i32) -> Self {
        debug_assert!(!(-1..=1).contains(&ratio));

        const PI_OVER_TWO: f32 = PI / 2.0;
        const BUTTERWORTH_COEFF: f32 = 2.613_126/*_929752753*/;
        const BUTTERWORTH_COEFF_TIMES_2: f32 = BUTTERWORTH_COEFF * 2.0;
        const OTHER_NUMBER: f32 =  3.414_213_7/*62373095*/;
        const OTHER_NUMBER_TIMES_TWO: f32 = OTHER_NUMBER * 2.0;

        // 4th-degree Butterworth with cutoff at ratio * Nyquist, using bilinear transform w/ prewarping
        let ratio_f = ratio.abs() as f32;
        let t = tanf(PI_OVER_TWO / ratio_f);
        let t2 = t * t;
        let k = 1.0
            / (t * (t * (t * (t + BUTTERWORTH_COEFF) + OTHER_NUMBER) + BUTTERWORTH_COEFF) + 1.0);
        Self {
            ratio,
            b0: k * t2 * t2,
            ma1: k
                * (t * (t2 * (-BUTTERWORTH_COEFF_TIMES_2 - 4.0 * t) + BUTTERWORTH_COEFF_TIMES_2)
                    + 4.0),
            ma2: k * ((OTHER_NUMBER_TIMES_TWO - 6.0 * t2) * t2 - 6.0),
            ma3: k
                * (t * (t2 * (BUTTERWORTH_COEFF_TIMES_2 - 4.0 * t) - BUTTERWORTH_COEFF_TIMES_2)
                    + 4.0),
            ma4: k
                * (t * (t * ((BUTTERWORTH_COEFF - t) * t - OTHER_NUMBER) + BUTTERWORTH_COEFF)
                    - 1.0),
        }
    }

    #[inline(always)]
    pub fn reset_state(&mut self, state: &mut SRCIntState, x0: f32) -> f32 {
        debug_assert!(x0.is_finite());

        if self.ratio < 0 {
            // DF-II
            state.z1 = x0 / (1.0 - self.ma1 - self.ma2 - self.ma3 - self.ma4);
            state.z2 = state.z1;
            state.z3 = state.z2;
            state.z4 = state.z3;
            state.i = 0;
        } else {
            // TDF-II
            let k = 4.0 * self.b0;
            state.z4 = (self.b0 + self.ma4) * x0;
            state.z3 = (k + self.ma3) * x0 + state.z4;
            state.z2 = (6.0 * self.b0 + self.ma2) * x0 + state.z3;
            state.z1 = (k + self.ma1) * x0 + state.z2;
        }

        debug_assert!(x0.is_finite());

        x0
    }

    #[inline(always)]
    pub fn reset_state_multi(
        &mut self,
        states: &mut [SRCIntState; N_CHANNELS],
        x0: &[f32],
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
    pub fn process(
        &mut self,
        state: &mut SRCIntState,
        x: &[f32],
        y: &mut [f32],
        n_in_sample: usize,
    ) -> usize {
        let mut n = 0;
        if self.ratio < 0 {
            (0..n_in_sample).for_each(|sample| {
                // DF-II
                let z0 = x[sample]
                    + self.ma1 * state.z1
                    + self.ma2 * state.z2
                    + self.ma3 * state.z3
                    + self.ma4 * state.z4;
                if state.i == 0 {
                    state.i = -self.ratio;
                    y[n] = self.b0 * (z0 + state.z4 + 4.0 * (state.z1 + state.z3) + 6.0 * state.z2);
                    n += 1;
                }
                state.i -= 1;
                state.z4 = state.z3;
                state.z3 = state.z2;
                state.z2 = state.z1;
                state.z1 = z0;
            });
        } else {
            (0..n_in_sample).for_each(|sample| {
                // TDF-II
                let in_ = self.ratio as f32 * x[sample];
                let v0 = self.b0 * in_;
                let v1 = 4.0 * v0;
                let v2 = 6.0 * v0;
                let mut o = v0 + state.z1;
                state.z1 = v1 + self.ma1 * o + state.z2;
                state.z2 = v2 + self.ma2 * o + state.z3;
                state.z3 = v1 + self.ma3 * o + state.z4;
                state.z4 = v0 + self.ma4 * o;
                y[n] = o;
                n += 1;
                println!("[rust] n = {}", n);
                (1..self.ratio).for_each(|_| {
                    o = state.z1;
                    state.z1 = self.ma1 * o + state.z2;
                    state.z2 = self.ma2 * o + state.z3;
                    state.z3 = self.ma3 * o + state.z4;
                    state.z4 = self.ma4 * o;
                    y[n] = o;
                    n += 1;
                });
            });
        }
        n
    }

    #[inline(always)]
    pub fn process_multi(
        &mut self,
        states: &mut [SRCIntState],
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_in_samples: usize,
        n_out_samples: Option<&mut [usize; N_CHANNELS]>,
    ) {
        match n_out_samples {
            Some(out) => (0..N_CHANNELS).for_each(|channel| {
                out[channel] =
                    self.process(&mut states[channel], x[channel], y[channel], n_in_samples);
            }),
            None => (0..N_CHANNELS).for_each(|channel| {
                self.process(&mut states[channel], x[channel], y[channel], n_in_samples);
            }),
        }
    }

    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Tries to determine whether state is valid and returns `true` if it seems to
    // /// be the case and `false` if it is certainly not. False positives are possible,
    // /// false negatives are not.
    // ///
    // /// # Note
    // /// <div class="warning">Not implemented yet!</div>
    // #[inline(always)]
    // pub fn coeffs_is_valid(&self) {
    //     todo!()
    // }
    // Not implemented yet:
    // need to revisit which assertions from the C version make sense to keep in Rust
    // /// Tries to determine whether state is valid and returns `true` if it seems to
    // /// be the case and `false` if it is certainly not. False positives are possible,
    // /// false negatives are not.
    // ///
    // /// # Note
    // /// <div class="warning">Not implemented yet!</div>
    // #[inline(always)]
    // pub fn state_is_valid(&self, state: &SRCIntState) {
    //     todo!()
    // }
}

#[derive(Default, Clone, Copy, Debug)]
pub struct SRCIntState {
    i: i32,
    z1: f32,
    z2: f32,
    z3: f32,
    z4: f32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::c_wrapper::{
        bw_src_int_coeffs as SRCIntCoeffsWrapper, bw_src_int_state,
        src_int::SRCInt as SRCIntWrapper,
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type SRCIntT = SRCInt<N_CHANNELS>;
    type SRCIntWrapperT = SRCIntWrapper<N_CHANNELS>;

    #[test]
    fn new_with_ratio_positive() {
        let ratio = 2.0;
        let rust_src_int = SRCIntT::new(ratio as i32);
        let c_src_int = SRCIntWrapperT::new(ratio as i32);

        assert_src_int(&rust_src_int, &c_src_int);
    }

    #[test]
    fn new_with_ratio_negative() {
        let ratio = -4.0;
        let rust_src_int = SRCIntT::new(ratio as i32);
        let c_src_int = SRCIntWrapperT::new(ratio as i32);

        assert_src_int(&rust_src_int, &c_src_int);
    }

    #[test]
    fn reset_none() {
        let ratio = -2.0;
        let mut rust_src_int = SRCIntT::new(ratio as i32);
        let mut c_src_int = SRCIntWrapperT::new(ratio as i32);

        rust_src_int.reset(None, None);
        c_src_int.reset(None, None);

        assert_src_int(&rust_src_int, &c_src_int);
    }

    #[test]
    fn reset_some() {
        let ratio = 2.0;
        let mut rust_src_int = SRCIntT::new(ratio as i32);
        let mut c_src_int = SRCIntWrapperT::new(ratio as i32);

        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_src_int.reset(Some(1.0), Some(&mut rust_y0));
        c_src_int.reset(Some(1.0), Some(&mut c_y0));

        assert_src_int(&rust_src_int, &c_src_int);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn reset_multi_down_sampling() {
        let ratio = -3.0;
        let mut rust_src_int = SRCIntT::new(ratio as i32);
        let mut c_src_int = SRCIntWrapperT::new(ratio as i32);

        let x0 = [0.5, 0.5];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_src_int.reset_multi(&x0, Some(&mut rust_y0));
        c_src_int.reset_multi(&x0, Some(&mut c_y0));

        assert_src_int(&rust_src_int, &c_src_int);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn reset_multi_up_sampling() {
        let ratio = 3.0;
        let mut rust_src_int = SRCIntT::new(ratio as i32);
        let mut c_src_int = SRCIntWrapperT::new(ratio as i32);

        let x0 = [0.5, 0.5];
        let mut rust_y0 = [0.0, 0.0];
        let mut c_y0 = [0.0, 0.0];

        rust_src_int.reset_multi(&x0, Some(&mut rust_y0));
        c_src_int.reset_multi(&x0, Some(&mut c_y0));

        assert_src_int(&rust_src_int, &c_src_int);
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(rust_y0[channel], c_y0[channel]);
        });
    }

    #[test]
    fn process_down_sampling() {
        let ratio = -2.0;
        let mut rust_src_int = SRCIntT::new(ratio as i32);
        let mut c_src_int = SRCIntWrapperT::new(ratio as i32);

        rust_src_int.reset(None, None);
        c_src_int.reset(None, None);

        let y_ch: Box<dyn Fn() -> [f32; 8]> = Box::new(|| std::array::from_fn(|_| 0.0));

        let mut rust_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];
        let mut c_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];

        let mut rust_n_out_samples: [usize; N_CHANNELS] = [0, 0];
        let mut c_n_out_samples: [usize; N_CHANNELS] = [0, 0];

        rust_src_int.process(
            &PULSE_INPUT,
            &mut rust_y,
            N_SAMPLES,
            Some(&mut rust_n_out_samples),
        );
        c_src_int.process(
            &PULSE_INPUT,
            &mut c_y,
            N_SAMPLES,
            Some(&mut c_n_out_samples),
        );

        assert_src_int(&rust_src_int, &c_src_int);
        assert_eq!(rust_y, c_y);
        assert_eq!(rust_n_out_samples, c_n_out_samples);
    }

    #[test]
    fn process_up_sampling() {
        const RATIO: f32 = 2.0;
        const N_SAMPLES_NEW: usize = N_SAMPLES * RATIO as usize;
        let mut rust_src_int = SRCIntT::new(RATIO as i32);
        let mut c_src_int = SRCIntWrapperT::new(RATIO as i32);

        rust_src_int.reset(None, None);
        c_src_int.reset(None, None);

        let y_ch: Box<dyn Fn() -> [f32; N_SAMPLES_NEW]> = Box::new(|| std::array::from_fn(|_| 0.0));

        let mut rust_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];
        let mut c_y: [&mut [f32]; 2] = [&mut y_ch(), &mut y_ch()];

        let mut rust_n_out_samples: [usize; N_CHANNELS] = [0, 0];
        let mut c_n_out_samples: [usize; N_CHANNELS] = [0, 0];

        c_src_int.process(
            &PULSE_INPUT,
            &mut c_y,
            N_SAMPLES,
            Some(&mut c_n_out_samples),
        );
        rust_src_int.process(
            &PULSE_INPUT,
            &mut rust_y,
            N_SAMPLES,
            Some(&mut rust_n_out_samples),
        );

        assert_src_int(&rust_src_int, &c_src_int);
        assert_eq!(rust_y, c_y);
        assert_eq!(rust_n_out_samples, c_n_out_samples);
    }

    fn assert_src_int<const N_CHANNELS: usize>(
        rust_src_int: &SRCInt<N_CHANNELS>,
        c_src_int: &SRCIntWrapper<N_CHANNELS>,
    ) {
        assert_src_int_coeffs(&rust_src_int.coeffs, &c_src_int.coeffs);
        (0..N_CHANNELS).for_each(|channel| {
            assert_src_int_state(&rust_src_int.states[channel], &c_src_int.states[channel]);
        });
    }

    fn assert_src_int_coeffs<const N_CHANNELS: usize>(
        rust_coeffs: &SRCIntCoeffs<N_CHANNELS>,
        c_coeffs: &SRCIntCoeffsWrapper,
    ) {
        assert_eq!(rust_coeffs.b0, c_coeffs.b0);
        assert_eq!(rust_coeffs.ma1, c_coeffs.ma1);
        assert_eq!(rust_coeffs.ma2, c_coeffs.ma2);
        assert_eq!(rust_coeffs.ma3, c_coeffs.ma3);
        assert_eq!(rust_coeffs.ma4, c_coeffs.ma4);
        assert_eq!(rust_coeffs.ratio, c_coeffs.ratio);
    }

    fn assert_src_int_state(rust_state: &SRCIntState, c_state: &bw_src_int_state) {
        assert_eq!(rust_state.z1, c_state.z1);
        assert_eq!(rust_state.z2, c_state.z2);
        assert_eq!(rust_state.z3, c_state.z3);
        assert_eq!(rust_state.z4, c_state.z4);
        assert_eq!(rust_state.i, c_state.i);
    }
}
