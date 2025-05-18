use crate::global::{assert_is_finite, assert_range};

//  Newton-Raphson reciprocal approsimation
#[inline(always)]
pub(crate) fn rcpf(x: f32) -> f32 {
    #[cfg(debug_assertions)]
    {
        assert_is_finite(x);
        if x > 0. {
            assert_range(8.077_936e-28, 1.237_940_1e27, x);
        } else {
            assert_range(-1.237_940_1e27, -8.077_936e-28, x);
        }
    }

    let magic = 0x7ef0e840u32;
    let u = x.to_bits();
    let mut result = f32::from_bits(magic - u);

    result = result + result - x * result * result;
    result = result + result - x * result * result;

    #[cfg(debug_assertions)]
    assert_is_finite(result);
    result
}

#[cfg(test)]
mod tests {
    use crate::{c_wrapper::bw_rcpf, native::math::rcpf};

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
    #[should_panic(expected = "Value must be in range [8.077936e-28, 1.2379401e27], got 8e-28")]
    #[test]
    fn rcpf_should_panic_with_small_number() {
        let value = 8.0e-28;
        rcpf(value);
    }

    #[cfg(debug_assertions)]
    #[should_panic(expected = "Value must be in range [-1.2379401e27, -8.077936e-28], got -1.3e27")]
    #[test]
    fn rcpf_should_panic_with_big_negative_number() {
        let value = -1.3e27;
        rcpf(value);
    }
}
