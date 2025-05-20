use std::f32::consts::PI;

pub(crate) const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);
pub(crate) const NANO: f32 = 1e-9;

pub(crate) fn assert_positive(value: f32) {
    assert!(value >= 0.0, "Value must be non negative, got {}!", value);
}

pub(crate) fn assert_range(min: f32, max: f32, value: f32) {
    assert!(
        value >= min && value <= max,
        "Value must be in range [{min:e}, {max:e}], got {value:e}!"
    );
}

// Need to contitionally add this whole module when refactored
// the assertion conditionally as well
#[allow(dead_code)]
pub(crate) fn assert_is_finite(value: f32) {
    debug_assert!(value.is_finite(), "value must be finite, got {}", value);
}
