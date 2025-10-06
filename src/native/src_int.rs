pub struct SRCInt<const N_CHANNELS: usize> {
    pub(crate) coeffs: SRCIntCoeffs<N_CHANNELS>,
    pub(crate) states: [SRCIntState; N_CHANNELS],
}

impl<const N_CHANNELS: usize> SRCInt<N_CHANNELS> {
    #[inline(always)]
    pub fn new(ratio: i32) -> Self {
        todo!()
    }
    #[inline(always)]
    pub fn reset(&mut self, x0: Option<f32>, y0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }
    #[inline(always)]
    pub fn reset_multi(&mut self, x0: &[f32], y0: Option<&mut [f32; N_CHANNELS]>) {
        todo!()
    }
    #[inline(always)]
    pub fn process(
        &mut self,
        x: &[&[f32]; N_CHANNELS],
        y: &mut [&mut [f32]; N_CHANNELS],
        n_in_samples: usize,
        n_out_samples: Option<&mut [usize; N_CHANNELS]>,
    ) {
        todo!()
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

pub struct SRCIntState {
    i: i32,
    z1: f32,
    z2: f32,
    z3: f32,
    z4: f32,
}

mod tests {
    use super::*;
    use crate::c_wrapper::{
        bw_src_int_coeffs as SRCIntCoeffsWrapper, bw_src_int_state,
        src_int::SRCInt as SRCIntWrapper,
    };

    const N_CHANNELS: usize = 2;
    const N_SAMPLES: usize = 8;
    const SAMPLE_RATE: f32 = 44_100.0;

    const PULSE_INPUT: [&[f32]; N_CHANNELS] = [
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    type SRCIntT = SRCInt<N_CHANNELS>;
    type SRCIntWrapperT = SRCIntWrapper<N_CHANNELS>;
}
