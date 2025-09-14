#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_range};

//  Newton-Raphson reciprocal approsimation
#[inline(always)]
pub fn rcpf(x: f32) -> f32 {
    #[cfg(debug_assertions)]
    {
        debug_assert_is_finite(x);
        if x > 0. {
            debug_assert_range(8.077_936e-28, 1.237_940_1e27, x);
        } else {
            debug_assert_range(-1.237_940_1e27, -8.077_936e-28, x);
        }
    }

    let magic = 0x7ef0e840u32;
    let u = x.to_bits();
    let mut result = f32::from_bits(magic - u);

    result = result + result - x * result * result;
    result = result + result - x * result * result;

    #[cfg(debug_assertions)]
    debug_assert_is_finite(result);
    result
}
/// Returns the minimum of `a` and `b`.
#[inline(always)]
pub fn minf(a: f32, b: f32) -> f32 {
    debug_assert!(!a.is_nan());
    debug_assert!(!b.is_nan());

    let y = a.min(b);

    debug_assert!(!y.is_nan());
    y
}
/// Returns the maximum of `a` and `b`.
#[inline(always)]
pub fn maxf(a: f32, b: f32) -> f32 {
    debug_assert!(!a.is_nan());
    debug_assert!(!b.is_nan());

    let y = a.max(b);

    debug_assert!(!y.is_nan());
    y
}
///    Returns `x` unless it is smaller than `m_small`, in which case it returns `m_small`,
///    or bigger than `m_big`, in which case it returns `m_big`.
///
///    `m_big` must be greater than or equal to `m_small`.
#[inline(always)]
pub fn clipf(x: f32, m_small: f32, m_big: f32) -> f32 {
    debug_assert!(!x.is_nan());
    debug_assert!(!m_small.is_nan());
    debug_assert!(!m_big.is_nan());
    debug_assert!(
        m_big > m_small,
        "m_small must be less than m_big, got {} and {}",
        m_small,
        m_big
    );

    let y = minf(maxf(x, m_small), m_big);

    debug_assert!(!y.is_nan());
    return y;
}
#[cfg(test)]
mod tests {
    use crate::{
        c_wrapper::bw_rcpf,
        native::math::{clipf, rcpf},
    };

    #[test]
    fn rcpf_should_return_same_result_as_bw_rcpf() {
        let values = vec![2.0, 345.0, 8.1e-28, 1.2e27];
        unsafe {
            values.iter().for_each(|value| {
                assert_eq!(rcpf(*value), bw_rcpf(*value));
            });
        }
    }

    #[cfg(debug_assertions)]
    #[should_panic(expected = "value must be in range [8.077936e-28, 1.2379401e27], got 8e-28")]
    #[test]
    fn rcpf_should_panic_with_small_number() {
        let value = 8.0e-28;
        rcpf(value);
    }

    #[cfg(debug_assertions)]
    #[should_panic(expected = "value must be in range [-1.2379401e27, -8.077936e-28], got -1.3e27")]
    #[test]
    fn rcpf_should_panic_with_big_negative_number() {
        let value = -1.3e27;
        rcpf(value);
    }

    #[test]
    fn clipf_valid() {
        assert_eq!(clipf(0.10, 0.20, 1.00), 0.20);
        assert_eq!(clipf(1.10, 0.20, 1.00), 1.00);
        assert_eq!(clipf(0.15, 0.10, 1.00), 0.15);
    }

    #[should_panic(expected = "m_small must be less than m_big, got 0.5 and 0.1")]
    #[test]
    fn clipf_invalid() {
        clipf(0.50, 0.50, 0.10);
    }
}
