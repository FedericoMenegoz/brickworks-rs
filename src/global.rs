use std::f32::consts::PI;
#[allow(dead_code)]
pub (crate) const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

pub (crate) fn assert_positive(value: f32) {
    assert!(value >= 0.0, "Value must be non negative, got {}!", value);
}

pub (crate) fn assert_range(min: f32, max: f32, value: f32) {
    assert!(
        value >= min && value <= max,
        "Value must be in range [{min:e}, {max:e}], got {value:e}!"
    );
}