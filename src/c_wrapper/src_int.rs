use super::*;
use std::ptr::null_mut;

pub struct SRCInt<const N_CHANNELS: usize> {
    pub(crate) coeffs: bw_src_int_coeffs,
    pub(crate) states: [bw_src_int_state; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SRCInt<N_CHANNELS> {
    #[inline(always)]
    pub fn new(ratio: i32) -> Self {
        let mut coeffs = bw_src_int_coeffs::default();
        unsafe {
            bw_src_int_init(&mut coeffs, ratio);
        }
        Self {
            coeffs,
            states: [bw_src_int_state::default(); N_CHANNELS],
        }
    }

    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        unsafe {
            match y0 {
                Some(y) => (0..N_CHANNELS).for_each(|channel| {
                    y[channel] = bw_src_int_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
                None => (0..N_CHANNELS).for_each(|channel| {
                    bw_src_int_reset_state(
                        &mut self.coeffs,
                        &mut self.states[channel],
                        x0.unwrap_or(0.0),
                    );
                }),
            }
        }
    }

    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        let y0_ptrs = match y0 {
            Some(y) => y.as_mut_ptr(),
            None => null_mut(),
        };
        let states_ptrs: [*mut bw_src_int_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {
            bw_src_int_reset_state_multi(
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
        n_in_samples: usize,
        n_out_samples: Option<&mut [usize; N_CHANNELS]>,
    ) {
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let mut y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());
        let states_ptrs: [*mut bw_src_int_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        let out_samples: *mut usize = match n_out_samples {
            Some(out) => out.as_mut_ptr(),
            None => null_mut(),
        };
        unsafe {
            bw_src_int_process_multi(
                &mut self.coeffs,
                states_ptrs.as_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.as_mut_ptr(),
                N_CHANNELS,
                n_in_samples,
                out_samples,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use crate::c_wrapper::{bw_tanf, src_int::SRCInt};

    const N_CHANNELS: usize = 2;
    const SAMPLE_RATE: f32 = 48_000.0;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    #[test]
    fn new() {
        let ratio = 2.0;
        let butterworth_coeff = 2.613_125_930/*752753*/;
        let butterworth_coeff_times_2 = butterworth_coeff * 2.0;
        let other_number =  3.414_213_562/*373095*/;
        let src_int = SRCInt::<N_CHANNELS>::new(ratio as i32);
        let (t, t2, k);
        unsafe {
            t = bw_tanf((PI / 2.0) / ratio);
            t2 = t * t;
            k = 1.0
                / (t * (t * (t * (t + butterworth_coeff) + other_number) + butterworth_coeff)
                    + 1.0);
        }
        assert_eq!(k * t2 * t2, src_int.coeffs.b0);
        assert_eq!(
            k * (t * (t2 * (-butterworth_coeff_times_2 - 4.0 * t) + butterworth_coeff_times_2)
                + 4.0),
            src_int.coeffs.ma1,
        );
        assert_eq!(
            k * ((6.828_427_12474619 - 6.0 * t2) * t2 - 6.0),
            src_int.coeffs.ma2
        );
        assert_eq!(
            k * (t * (t2 * (butterworth_coeff_times_2 - 4.0 * t) - butterworth_coeff_times_2)
                + 4.0),
            src_int.coeffs.ma3
        );
        assert_eq!(
            k * (t * (t * ((butterworth_coeff - t) * t - other_number) + butterworth_coeff) - 1.0),
            src_int.coeffs.ma4
        );
    }

    #[test]
    fn reset() {
        let ratio = -3.0;
        let mut src_int = SRCInt::<N_CHANNELS>::new(ratio as i32);

        let x0 = 0.2;
        let mut y0: [f32; N_CHANNELS] = [0.0, 0.0];

        src_int.reset(Some(x0), Some(&mut y0));

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(
                src_int.states[channel].z1,
                x0 / (1.0
                    - src_int.coeffs.ma1
                    - src_int.coeffs.ma2
                    - src_int.coeffs.ma3
                    - src_int.coeffs.ma4)
            );
            assert_eq!(src_int.states[channel].z2, src_int.states[channel].z1);
            assert_eq!(src_int.states[channel].z3, src_int.states[channel].z1);
            assert_eq!(src_int.states[channel].z4, src_int.states[channel].z1);
            assert_eq!(src_int.states[channel].i, 0);
            assert_eq!(y0[channel], x0);
        });
    }

    #[test]
    fn reset_multi() {
        let ratio = 3.0;
        let mut src_int = SRCInt::<N_CHANNELS>::new(ratio as i32);

        let x0: [f32; N_CHANNELS] = [0.01, 0.02];
        let mut y0: [f32; N_CHANNELS] = [0.0, 0.0];

        src_int.reset_multi(&x0, Some(&mut y0));
        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(y0[channel], x0[channel]);
        });
    }

    #[test]
    fn process() {
        let ratio = 2.0;
        let mut src_int = SRCInt::<N_CHANNELS>::new(ratio as i32);
        let n_in_samples = 2;

        let x: [&[f32]; 2] = [&[0.2, 0.3], &[0.4, 0.5]];
        let mut y: [&mut [f32]; 2] = [&mut [0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0]];
        let mut n_out_samples: [usize; N_CHANNELS] = [0, 0];
        src_int.reset(None, None);
        src_int.process(&x, &mut y, n_in_samples, Some(&mut n_out_samples));

        (0..N_CHANNELS).for_each(|channel| {
            assert_eq!(n_out_samples[channel], n_in_samples * ratio as usize);
        });
    }
}
