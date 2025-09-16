use std::ptr::null_mut;

#[inline(always)]
pub(crate) fn make_array<T: Default, const N: usize>() -> [T; N] {
    std::array::from_fn(|_| T::default())
}

#[inline(always)]
pub (crate) fn from_opt_to_raw<const N_CHANNELS: usize>(buffer: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>) -> [*mut f32; N_CHANNELS] {
    match buffer {
        Some(buf_ch) => std::array::from_fn(|i| {
            if let Some(ch) = &mut buf_ch[i] {
                ch.as_mut_ptr()
            } else {
                null_mut()
            }
        }),
        None => [null_mut(); N_CHANNELS],
    }
}