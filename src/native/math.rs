use core::f32::{self, consts::PI};

pub const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);
pub const NANO: f32 = 1e-9;
pub const PI_OVER_2: f32 = PI / 2.0;

#[cfg(debug_assertions)]
use crate::native::common::{debug_assert_is_finite, debug_assert_range};

//  Newton-Raphson reciprocal approsimation
#[inline(always)]
pub fn rcpf(x: f32) -> f32 {
    #[cfg(debug_assertions)]
    {
        debug_assert_is_finite(x);
        if x > 0. {
            debug_assert_range(8.077_936e-28..=1.237_940_1e27, x);
        } else {
            debug_assert_range(-1.237_940_1e27..=-8.077_936e-28, x);
        }
    }
    let mut result = f32::from_bits(0x7ef0e840u32.wrapping_sub(x.to_bits()));

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
/// Returns `x` unless it is smaller than `m_small`, in which case it returns `m_small`,
/// or bigger than `m_big`, in which case it returns `m_big`.
///
/// `m_big` must be greater than or equal to `m_small`.
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
    y
}
/// Returns an approximation of the tangent of `x`, where `x` is given in
/// radians.
///
/// `x` must be finite and in [-pi/2 + 1e-3f, pi/2 - 1e-3f] + k * pi, where k
/// is any integer number.
///
/// Absolute error < 0.06 or relative error < 0.8%, whatever is worse.
#[inline(always)]
pub fn tanf(mut x: f32) -> f32 {
    debug_assert!(x.is_finite());
    debug_assert!(
        (x - PI * (INVERSE_2_PI * x).floor() <= PI_OVER_2 - 1e-3)
            || (x - PI * (INVERSE_2_PI * x).floor() >= PI_OVER_2 + 1e-3),
        "value must be in range [-pi/2 + 1e-3f, pi/2 - 1e-3f] + k * pi got {}",
        x
    );
    x *= INVERSE_2_PI;
    let y: f32 = sin2pif(x) * rcpf(cos2pif(x));
    debug_assert!(y.is_finite());
    y
}
/// Returns an approximation of the sine of 2 * pi * `x`, where `x` is given
/// in radians.
///
/// `x` must be finite.
///
/// Absolute error < 0.011 or relative error < 1.7%, whatever is worse.
#[inline(always)]
pub fn sin2pif(mut x: f32) -> f32 {
    debug_assert!(x.is_finite());
    x = x - x.floor();
    let xp1 = x + x - 1.0;
    let xp2 = xp1.abs();
    let xp = PI_OVER_2 - PI_OVER_2 * (xp2 + xp2 - 1.0).abs();
    let y = -1.0_f32.copysign(xp1) * (xp + xp * xp * (-0.05738534 - 0.11073982 * xp));
    debug_assert!(y.is_finite());
    y
}
///
/// Returns an approximation of the cosine of 2 * pi * `x`, where `x` is given
/// in radians.
///
/// `x` must be finite.
///
/// Absolute error < 0.011 or relative error < 1.7%, whatever is worse.
#[inline(always)]
pub fn cos2pif(x: f32) -> f32 {
    debug_assert!(x.is_finite());
    let y = sin2pif(x + 0.25);
    debug_assert!(y.is_finite());
    y
}
/// Returns an approximation of 10 raised to the power of `x` / 20 (dB to
/// linear ratio conversion). For `x < -758.5955890732315f` it just returns
/// `0.f`.
///
/// `x` must be less than or equal to `770.630f`.
///
/// Relative error < 0.062%.
#[inline(always)]
pub fn db2linf(x: f32) -> f32 {
    debug_assert!(!x.is_nan());
    debug_assert!(
        x <= 770.630,
        "value must be less or equal to 770.630, got {}",
        x
    );

    let y = pow2f(0.1660964/*047443682*/ * x);

    debug_assert!(y.is_finite());
    y
}

/// Returns an approximation of the square root of `x`.
///
/// `x` must be finite and non-negative.
///
/// Absolute error < 1.09e-19 or relative error < 0.0007%, whatever is worse.
#[inline(always)]
pub fn sqrtf(x: f32) -> f32 {
    debug_assert!(x.is_finite(), "value must be finite, got {}", x);
    debug_assert!(x >= 0.0, "value must be non negative, got {}", x);
    // if x < 1.1754943508222875e-38 { return 0.0 };
    if x < 1.175_494_3e-38 {
        return 0.0;
    };

    let mut v_bits = x.to_bits();

    let i = (v_bits >> 26) & 0x38;
    v_bits += (0x200000e0 << i) & 0xff00_0000;

    let mut v = f32::from_bits(v_bits);

    let r = rcpf(v);

    v_bits = (v.to_bits().wrapping_sub(0x3f82a127) >> 1).wrapping_add(0x3f7d8fc7) & 0x7fff_ffff;
    v = f32::from_bits(v_bits);

    v = v + v * (0.5 - 0.5 * r * v * v);
    v = v + v * (0.5 - 0.5 * r * v * v);

    v_bits = v.to_bits();
    v_bits -= (0x100000f0 << i) & 0xff00_0000;
    v = f32::from_bits(v_bits);

    debug_assert!(v.is_finite());

    v
}
/// Returns an approximation of 2 raised to the power of `x`. For `x < -126.f`
/// it just returns `0.f`.
///
/// `x` must be less than or equal to `127.999f`.
///
/// Relative error < 0.062%.
///
#[inline(always)]
pub fn pow2f(x: f32) -> f32 {
    debug_assert!(!x.is_nan());
    debug_assert!(
        x <= 127.999,
        "value must be less or equal to 127.999, got {}",
        x
    );

    if x < -126.0 {
        return 0.0;
    }

    let xi = x as i32;
    let v_bits = x.to_bits();
    let l = xi - ((v_bits >> 31) & 1) as i32;
    let f = x - l as f32;

    let v_i = (l + 127) << 23;
    let v_f = f32::from_bits(v_i as u32);

    let y =
        //excessive precision v_f + v_f * f * (0.6931471805599453 + f * (0.2274112777602189 + f * 0.07944154167983575));
        v_f + v_f * f * (f32::consts::LN_2 + f * (0.22741129 + f * 0.07944154));

    debug_assert!(y.is_finite());
    y
}
#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use crate::{
        c_wrapper::{
            bw_clipf, bw_cos2pif, bw_dB2linf, bw_pow2f, bw_rcpf, bw_sin2pif, bw_sqrtf, bw_tanf,
        },
        native::math::{clipf, cos2pif, db2linf, pow2f, rcpf, sin2pif, sqrtf, tanf},
    };

    #[test]
    fn rcpf_should_return_same_result_as_bw_rcpf() {
        let values = vec![-0.98, 2.0, 345.0, 8.1e-28, 1.2e27];
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
        unsafe {
            assert_eq!(clipf(0.10, 0.20, 1.00), bw_clipf(0.10, 0.20, 1.00));
            assert_eq!(clipf(1.10, 0.20, 1.00), bw_clipf(1.10, 0.20, 1.00));
            assert_eq!(clipf(0.15, 0.10, 1.00), bw_clipf(0.15, 0.10, 1.00));
        }
    }

    #[should_panic(expected = "m_small must be less than m_big, got 0.5 and 0.1")]
    #[test]
    fn clipf_invalid() {
        clipf(0.50, 0.50, 0.10);
    }

    #[test]
    fn tanf_valid() {
        unsafe {
            assert_eq!(tanf(-PI / 4.0), bw_tanf(-PI / 4.0));
        }
    }

    #[test]
    fn sin2pif_valid() {
        unsafe {
            assert_eq!(sin2pif(-PI / 4.0), bw_sin2pif(-PI / 4.0));
        }
    }

    #[test]
    fn cos2pif_valid() {
        unsafe {
            assert_eq!(cos2pif(-PI / 4.0), bw_cos2pif(-PI / 4.0));
        }
    }

    #[should_panic(
        expected = "value must be in range [-pi/2 + 1e-3f, pi/2 - 1e-3f] + k * pi got -1.5707964"
    )]
    #[test]
    fn tanf_invalid() {
        tanf(-PI / 2.0);
    }

    #[test]
    fn db2linf_valid() {
        unsafe {
            assert_eq!(db2linf(-1000.0), bw_dB2linf(-1000.0));
            assert_eq!(db2linf(10.0), bw_dB2linf(10.0));
            assert_eq!(db2linf(500.0), bw_dB2linf(500.0));
        }
    }

    #[should_panic(expected = "value must be less or equal to 770.630, got 770.6301")]
    #[test]
    fn db2linf_invalid() {
        db2linf(770.6301);
    }

    #[test]
    fn pow2f_valid() {
        let mut mine;
        let mut bw;

        unsafe {
            for i in 0..10 {
                let p = i as f32 + 0.1;
                mine = pow2f(p);
                bw = bw_pow2f(p);
                assert_eq!(mine, bw);
            }
            assert_eq!(pow2f(45.0), bw_pow2f(45.0));
            assert_eq!(pow2f(127.99), bw_pow2f(127.99));
            assert_eq!(pow2f(-1000.0), bw_pow2f(-1000.0));
        }
    }

    #[should_panic(expected = "value must be less or equal to 127.999, got 128")]
    #[test]
    fn pow2f_invalid() {
        pow2f(128.0);
    }

    #[test]
    fn sqrtf_valid() {
        let ns = [23.0, 0.00004, 65.1, 123559.0, 1.1e-38];

        unsafe {
            for val in ns {
                assert_eq!(sqrtf(val), bw_sqrtf(val));
            }
        }
    }
}
