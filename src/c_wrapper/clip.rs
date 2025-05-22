use super::{bw_clip_coeffs, bw_clip_state};


pub struct Clip<const N_CHANNELS: usize> {
    coeffs: bw_clip_coeffs,
    state: [bw_clip_state; N_CHANNELS],
    states_p: [bw_clip_state; N_CHANNELS],
}

impl <const N_CHANNELS: usize> Clip<N_CHANNELS> {
    pub fn new() {}
    pub fn set_sample_rate(sample_rate: f32) {}
    pub fn reset(&mut self, x0: &[f32], mut y0: Option<&mut [f32]>) {}
    pub fn process(
        &mut self,
        x: &[Vec<f32>],
        y: Option<&mut [Option<&mut [f32]>]>,
        n_samples: usize,
    ) {}
    pub fn set_bias(value: f32){}
    pub fn set_gain_compensation(value: bool){}
}

mod tests {
    
}