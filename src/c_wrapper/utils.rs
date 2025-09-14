#[inline(always)]
pub(crate) fn make_array<T: Default, const N: usize>() -> [T; N] {
    std::array::from_fn(|_| T::default())
}

// Need to prepare the data into raw pointers for c
#[inline(always)]
pub fn prepare_input_and_optional_output_states_ptrs<State, const N_CHANNELS: usize>(
    x: &[&[f32]; N_CHANNELS],
    y: Option<&mut [Option<&mut [f32]>; N_CHANNELS]>,
    states: &mut [State; N_CHANNELS],
) -> (
    [*const f32; N_CHANNELS],
    [*mut f32; N_CHANNELS],
    [*mut State; N_CHANNELS],
) {
    // Input pointers
    let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());

    // Output pointers (null if missing)
    let y_ptrs: [*mut f32; N_CHANNELS] = match y {
        Some(y_channels) => std::array::from_fn(|i| {
            if let Some(buf) = &mut y_channels[i] {
                buf.as_mut_ptr()
            } else {
                std::ptr::null_mut()
            }
        }),
        None => [std::ptr::null_mut(); N_CHANNELS],
    };

    // State pointers
    let state_ptrs: [*mut State; N_CHANNELS] =
        std::array::from_fn(|i| &mut states[i] as *mut State);

    (x_ptrs, y_ptrs, state_ptrs)
}

pub fn prepare_input_and_output_states_ptrs<State, const N_CHANNELS: usize>(
    x: &[&[f32]; N_CHANNELS],
    y: &mut [&mut [f32]; N_CHANNELS],
    states: &mut [State; N_CHANNELS],
) -> (
    [*const f32; N_CHANNELS],
    [*mut f32; N_CHANNELS],
    [*mut State; N_CHANNELS],
) {
    let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());

    let y_ptrs: [*mut f32; N_CHANNELS] = std::array::from_fn(|i| y[i].as_mut_ptr());

    let state_ptrs: [*mut State; N_CHANNELS] =
        std::array::from_fn(|i| &mut states[i] as *mut State);

    (x_ptrs, y_ptrs, state_ptrs)
}
