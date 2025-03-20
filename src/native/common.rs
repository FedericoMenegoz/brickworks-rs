fn is_inf(x: f32) -> bool {
    todo!()
}

fn is_nan(x: f32) -> bool {
    todo!()
}

fn is_finite(x: f32) -> bool {
    todo!()
}

fn has_inf(x: &[f32]) -> bool {
    todo!()
}

fn has_nan(x: &[f32]) -> bool {
    todo!()
}

fn has_only_finite(x: &[f32]) -> bool {
    todo!()
}

fn hash_sdbm(x: &[i8]) -> u32 {
    todo!()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::c_wrapper::*;

    // IEEE Floating-Point Representation
    const INFTY: f32 = f32::from_bits(0xFF800000);
    const NAN: f32 = f32::from_bits(0x7FC00000);

    // BW_IS_INF
    #[test]
    fn is_inf_is_false_for_positive_numbers() {
        assert_eq!(is_inf(2.0), false);
        unsafe {
            assert_eq!(is_inf(2.0), bw_is_inf(2.0) != 0);
        }
    }
    #[test]
    fn is_inf_is_false_for_negative_numbers() {
        assert_eq!(is_inf(-2.0), false);
        unsafe {
            assert_eq!(is_inf(-2.0), bw_is_inf(-2.0) != 0);
        }
    }
    #[test]
    fn is_inf_is_false_for_zero() {
        assert_eq!(is_inf(0.0), false);
        unsafe {
            assert_eq!(is_inf(0.0), bw_is_inf(0.0) != 0);
        }
    }
    #[test]
    fn is_inf_is_true_for_plus_infinity() {
        assert_eq!(is_inf(INFTY), true);
        unsafe {
            assert_eq!(is_inf(INFTY), bw_is_inf(INFTY) != 0);
        }
    }
    #[test]
    fn is_inf_is_true_for_minus_infinity() {
        assert_eq!(is_inf(INFTY), true);
        unsafe {
            assert_eq!(is_inf(INFTY), bw_is_inf(INFTY) != 0);
        }
    }

    // BW_IS_NAN
    #[test]
    fn is_nan_is_true_for_nan() {
        assert_eq!(is_nan(NAN), true);
        unsafe {
            assert_eq!(is_nan(NAN), bw_is_nan(NAN) != 0);
        }
    }
    #[test]
    fn is_nan_is_false_for_number() {
        assert_eq!(is_nan(2.0), false);
        unsafe {
            assert_eq!(is_nan(2.0), bw_is_nan(2.0) != 0);
        }
    }

    // BW_IS_FINITE
    #[test]
    fn is_finite_is_false_for_infinity() {
        assert_eq!(is_finite(INFTY), false);
        unsafe {
            assert_eq!(is_finite(INFTY), bw_is_finite(INFTY) != 0);
        }
    }
    #[test]
    fn is_finite_is_true_for_number() {
        let x = 2.0;
        assert_eq!(is_finite(x), true);
        unsafe {
            assert_eq!(is_finite(x), bw_is_finite(x) != 0);
        }
    }

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
        let input = [3, 2, 1];
        unsafe {
            assert_eq!(hash_sdbm(&input), bw_hash_sdbm(input.as_ptr()));
        }
    }
}
