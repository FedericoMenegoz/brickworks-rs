fn has_inf(x: &[f32]) -> bool {
    x.iter().any(|val| val.is_infinite())
}

fn has_nan(x: &[f32]) -> bool {
    x.iter().any(|val| val.is_nan())
}

fn has_only_finite(x: &[f32]) -> bool {
    x.iter().all(|val| val.is_finite())
}

fn hash_sdbm(s: &str) -> u32 {
    let mut hash: u32 = 0;

    for ch in s.encode_utf16() {
        hash = u32::from(ch)
            .wrapping_add(hash << 6)
            .wrapping_add(hash << 16)
            .wrapping_sub(hash);
    }

    hash
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::c_wrapper::*;

    // IEEE Floating-Point Representation
    const INFTY: f32 = f32::from_bits(0xFF800000);
    const NAN: f32 = f32::from_bits(0x7FC00000);

    // BW_HAS_INF
    #[test]
    fn has_inf_is_false_for_normal_array() {
        let arr = [1.0, 2.0, 3.0];
        assert_eq!(has_inf(&arr), false);
        unsafe {
            assert_eq!(has_inf(&arr), bw_has_inf(arr.as_ptr(), arr.len()) != 0);
        }
    }
    #[test]
    fn has_inf_is_true_if_contains_infinity() {
        let arr = [1.0, INFTY, 3.0];
        assert_eq!(has_inf(&arr), true);
        unsafe {
            assert_eq!(has_inf(&arr), bw_has_inf(arr.as_ptr(), arr.len()) != 0);
        }
    }

    // BW_HAS_NAN
    #[test]
    fn has_nan_is_false_for_normal_array() {
        let arr = [1.0, 2.0, 3.0];
        assert_eq!(has_nan(&arr), false);
        unsafe {
            assert_eq!(has_nan(&arr), bw_has_nan(arr.as_ptr(), arr.len()) != 0);
        }
    }
    #[test]
    fn has_nan_is_true_if_contains_nan() {
        let arr = [1.0, f32::NAN, 3.0];
        assert_eq!(has_nan(&arr), true);
        unsafe {
            assert_eq!(has_nan(&arr), bw_has_nan(arr.as_ptr(), arr.len()) != 0);
        }
    }

    // BW_HAS_ONLY_FINITE
    #[test]
    fn has_only_finite_is_true_for_normal_array() {
        let arr = [1.0, 2.0, 3.0];
        assert_eq!(has_only_finite(&arr), true);
        unsafe {
            assert_eq!(
                has_only_finite(&arr),
                bw_has_only_finite(arr.as_ptr(), arr.len()) != 0
            );
        }
    }
    #[test]
    fn has_only_finite_is_false_if_contains_inf_or_nan() {
        let arr = [1.0, f32::INFINITY, 3.0];
        assert_eq!(has_only_finite(&arr), false);
        unsafe {
            assert_eq!(
                has_only_finite(&arr),
                bw_has_only_finite(arr.as_ptr(), arr.len()) != 0
            );
        }

        let arr_nan = [1.0, f32::NAN, 3.0];
        assert_eq!(has_only_finite(&arr_nan), false);
        unsafe {
            assert_eq!(
                has_only_finite(&arr_nan),
                bw_has_only_finite(arr_nan.as_ptr(), arr_nan.len()) != 0
            );
        }
    }

    // BW_HAS_SDBM
    #[test]
    fn bw_hash_sdbm_is_consistent() {
        let input = "test";
        unsafe {
            assert_eq!(hash_sdbm(&input), bw_hash_sdbm(input.as_ptr() as *const i8));
        }
    }
}
