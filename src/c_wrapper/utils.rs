#[inline(always)]
pub(crate) fn make_array<T: Default, const N: usize>() -> [T; N] {
    std::array::from_fn(|_| T::default())
}
