use crate::native::one_pole::{OnePoleCoeffs, OnePoleState};

pub struct Clip<const N_CHANNELS: usize> {
    coeffs: ClipCoeffs,
    state: [ClipState; N_CHANNELS],
    _state_p: [ClipState; N_CHANNELS],
}

pub struct ClipCoeffs {
    // Sub-components
    smooth_coeffs: OnePoleCoeffs,
    smooth_bias_state: OnePoleState,
    smooth_gain_state: OnePoleState,

    // Coefficients
    bias_dc: f32,
    inv_gain: f32,

    // Parameters
    bias: f32,
    gain: f32,
    gain_compensation: bool,
}

pub struct ClipState {
    x_z1: f32,
    f_z1: f32,
}

#[cfg(test)]
mod tests {
    use crate::{
        c_wrapper::bw_clip_coeffs,
        native::{clip::ClipCoeffs, one_pole::tests::assert_one_pole_coeffs_rust_c},
    };

    fn assert_clip_coeffs_rust_c(rust_coeffs: ClipCoeffs, c_coeffs: bw_clip_coeffs) {
        let pre_message = "clip.coeff.";
        let post_message = "does not match";
        assert_eq!(
            rust_coeffs.bias_dc, c_coeffs.bias_dc,
            "{}bias_dc {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.inv_gain, c_coeffs.inv_gain,
            "{}inv_gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.bias, c_coeffs.bias,
            "{}bias {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain, c_coeffs.gain,
            "{}gain {}",
            pre_message, post_message
        );
        assert_eq!(
            rust_coeffs.gain_compensation,
            c_coeffs.gain_compensation != 0,
            "{}gain_compensation {}",
            pre_message,
            post_message
        );
        assert_one_pole_coeffs_rust_c(rust_coeffs.smooth_coeffs, c_coeffs.smooth_coeffs);
    }
}
