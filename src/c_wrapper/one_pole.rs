use super::*;
use std::ptr::null_mut;


/// One-pole filter: Rust binding to the C library by [Orastron](https://www.orastron.com/algorithms/bw_one_pole).
///
/// One-pole (6 dB/oct) lowpass filter with unitary DC gain, separate attack and decay time constants, and sticky target-reach threshold.
///
/// This is better suited to implement smoothing than bw_lp1.
///
/// # Example
/// ```
/// use brickworks_rs::c_wrapper::one_pole::*;
///
/// const CUTOFF: f32 = 100_000.0;
/// const STICKY_THRESH: f32 = 0.9;
/// const SAMPLE_RATE: f32 = 48_100.0;
/// const N_CHANNELS: usize = 2;
/// const N_SAMPLES: usize = 1;
/// fn main() {
///     // Create a new OnePole filter instance for N_CHANNELS
///     let mut one_pole = OnePole::<N_CHANNELS>::new();
///
///     // Input signal: one sample per channel
///     let x = vec![vec![1.0], vec![33.0]];
///
///     // Output buffer, same shape as input
///     let mut y = vec![vec![0.0], vec![0.0]];
///
///     // Configure the filter
///     one_pole.set_sample_rate(SAMPLE_RATE);
///     one_pole.set_cutoff(CUTOFF);
///     one_pole.set_sticky_mode(OnePoleStickyMode::Rel);
///     one_pole.set_sticky_thresh(STICKY_THRESH);
///
///     // Initialize the filter state for each channel
///     one_pole.reset(&[0.0, 0.0], None);
///     // !!! This is ugly yeah...need to investigate
///     let mut y_wrapped: Vec<Option<&mut [f32]>> =
///         y.iter_mut().map(|ch| Some(&mut **ch)).collect();
///
///     // Process one sample per channel
///     one_pole.process(&x, Some(&mut y_wrapped), N_SAMPLES);
///
///     // Output the filtered result
///     println!("Filtered output: {:?}", y);
/// }
///
/// ```
/// 
/// # Notes
/// This module provides Rust bindings to the original C implementation.
/// For a fully native Rust implementation with the same interface, see [crate::native::one_pole].
pub struct OnePole<const N_CHANNELS: usize> {
    pub coeffs: bw_one_pole_coeffs,
    pub states: [bw_one_pole_state; N_CHANNELS],
    pub states_p: [bw_one_pole_state; N_CHANNELS], // BW_RESTRICT to check what is for
}

/// Distance metrics for sticky behavior.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OnePoleStickyMode {
    /// `OnePoleStickyMode::Abs`: absolute difference ( `|out - in|` );
    Abs,
    /// `OnePoleStickyMode::Rel`: relative difference with respect to input (`|out - in| / |in|`).
    Rel,
}

impl<const N_CHANNELS: usize> OnePole<N_CHANNELS> {
    /// Creates a new `OnePole` filter with default parameters and zeroed state.
    pub fn new() -> Self {
        let states: [bw_one_pole_state; N_CHANNELS] =
            std::array::from_fn(|_| bw_one_pole_state { y_z1: 0.0 });

        let states_p: [bw_one_pole_state; N_CHANNELS] =
            std::array::from_fn(|_| bw_one_pole_state { y_z1: 0.0 });

        let mut one_pole = OnePole {
            coeffs: bw_one_pole_coeffs {
                fs_2pi: Default::default(),
                mA1u: Default::default(),
                mA1d: Default::default(),
                st2: Default::default(),
                cutoff_up: Default::default(),
                cutoff_down: Default::default(),
                sticky_thresh: Default::default(),
                sticky_mode: Default::default(),
                param_changed: Default::default(),
            },
            states,
            states_p,
        };
        unsafe {
            bw_one_pole_init(&mut one_pole.coeffs);
        }
        one_pole
    }
    /// Sets the filter's sample rate, required for accurate time constant conversion.
    pub fn set_sample_rate(&mut self, sample_rate: f32) {
        unsafe {
            bw_one_pole_set_sample_rate(&mut self.coeffs, sample_rate);
        }
    }
    /// Resets the filter state for all channels using the initial input `x0`.
    /// If `y0` is provided, the resulting initial outputs are stored in it.
    pub fn reset(&mut self, x0: &[f32], mut y0: Option<&mut [f32]>) {
        unsafe {
            bw_one_pole_reset_coeffs(&mut self.coeffs);
            (0..N_CHANNELS).for_each(|channel| {
                let result = bw_one_pole_reset_state(
                    &mut self.coeffs,
                    &mut self.states[channel],
                    x0[channel],
                );
                if let Some(slice) = y0.as_deref_mut() {
                    slice[channel] = result
                }
            });
        }
    }
    /// Processes a block of input samples using a one-pole filter per channel.
    ///
    /// This method applies a one-pole low-pass filter to `n_samples` of input per channel.
    /// It supports multiple filtering modes based on whether the filter is *asymmetric*
    /// (`cutoff_up != cutoff_down`), and whether it uses *sticky*
    /// behavior (which prevents updates when the change is below a threshold).
    ///
    /// ### Parameters
    /// - `x`: A slice of input sample vectors, one per channel.
    /// - `y`: An optional mutable slice of output buffers, one per channel.
    /// - `n_samples`: Number of samples to process per channel.
    ///
    /// ### Filter Behavior
    /// The internal logic chooses the processing mode based on the current coefficients:
    ///
    /// - `process1`: Symmetric filtering with no sticky behavior.
    /// - `process1_sticky_abs`: Symmetric filtering with sticky behavior using an absolute difference metric.
    /// - `process1_sticky_rel`: Symmetric filtering with sticky behavior using a relative difference metric.
    /// - `process1_asym`: Asymmetric filtering with no sticky behavior.
    /// - `process1_asym_sticky_abs`: Asymmetric filtering with sticky behavior using an absolute difference metric.
    /// - `process1_asym_sticky_rel`: Asymmetric filtering with sticky behavior using a relative difference metric.
    ///
    pub fn process(
        &mut self,
        x: &[Vec<f32>],
        y: Option<&mut [Option<&mut [f32]>]>,
        n_samples: usize,
    ) {
        // In case y is None this will be passed
        let null_ptrs: [*mut f32; N_CHANNELS] = [null_mut(); N_CHANNELS];
        // Need to prepare the data into raw pointers for c
        let x_ptrs: [*const f32; N_CHANNELS] = std::array::from_fn(|i| x[i].as_ptr());
        let y_ptrs: Option<[*mut f32; N_CHANNELS]> = y.map(|y_channel| {
            std::array::from_fn(|i| {
                // state[i].unwrap_or(null_mut::<[f32]>()).as_mut_ptr()
                if let Some(y_samples) = y_channel[i].as_mut() {
                    y_samples.as_mut_ptr()
                } else {
                    null_ptrs[0]
                }
            })
        });
        let mut state_ptrs: [*mut bw_one_pole_state; N_CHANNELS] =
            std::array::from_fn(|i| &mut self.states[i] as *mut _);
        unsafe {

            bw_one_pole_process_multi(
                &mut self.coeffs,
                state_ptrs.as_mut_ptr(),
                x_ptrs.as_ptr(),
                y_ptrs.unwrap_or(null_ptrs).as_mut_ptr(),
                N_CHANNELS,
                n_samples,
            );
        }
    }
    /// Sets both the upgoing (attack) and downgoing (decay) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling both `set_cutoff_up()` and `set_cutoff_down()` with same
    /// value or calling `set_tau()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_cutoff(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff(&mut self.coeffs, value);
        }
    }
    /// Sets the upgoing (attack) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_tau_up()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_cutoff_up(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff_up(&mut self.coeffs, value);
        }
    }
    /// Sets the downgoing (decay) cutoff frequency to the given value (Hz) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: cutoff value to be set, default is `f32::INFINITE`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_tau_down()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_cutoff_down(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_cutoff_down(&mut self.coeffs, value);
        }
    }
    /// Sets both the upgoing (attack) and downgoing (decay) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling both `set_tau_up()` and `set_tau_down()` with same
    /// value or calling set_cutoff() with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_tau(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau(&mut self.coeffs, value);
        }
    }
    /// Sets the upgoing (attack) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_cutoff_up()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_tau_up(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau_up(&mut self.coeffs, value);
        }
    }
    /// Sets both the downgoing (decay) time constant to the given value (s) in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: tau value to be set, default is `f32::0.0`, must be non-negative.
    ///
    /// ### Note
    /// This is equivalent to calling `set_cutoff_down()` with value = 1 / (2 * pi * value) (net of numerical errors).
    ///
    pub fn set_tau_down(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_tau_down(&mut self.coeffs, value);
        }
    }
    /// Sets the target-reach threshold specified by value in `OnePole<N_CHANNELS>`.
    ///
    /// When the difference between the output and the input would fall under such threshold according
    /// to the current distance metric (absolute or relative), the output is forcefully set to be equal to the input value.
    ///
    ///
    /// ### Parameters
    /// - `value`: sticky threshhold, default is `f32::0.0` and its valid range is [0.f, 1e18f].
    ///    
    pub fn set_sticky_thresh(&mut self, value: f32) {
        unsafe {
            bw_one_pole_set_sticky_thresh(&mut self.coeffs, value);
        }
    }
    /// Sets the current distance metric for sticky behavior to value in `OnePole<N_CHANNELS>`.
    ///
    /// ### Parameters
    /// - `value`: sticky mode, default is absolute.
    ///
    pub fn set_sticky_mode(&mut self, value: OnePoleStickyMode) {
        unsafe {
            bw_one_pole_set_sticky_mode(&mut self.coeffs, value.to_bw());
        }
    }
    /// Returns the current target-reach threshold in `OnePole<N_CHANNELS>`.
    pub fn get_sticky_thresh(&self) -> f32 {
        unsafe { bw_one_pole_get_sticky_thresh(&self.coeffs) }
    }
    /// Returns the current distance metric for sticky behavior in `OnePole<N_CHANNELS>`.
    pub fn get_sticky_mode(&self) -> OnePoleStickyMode {
        unsafe { OnePoleStickyMode::from_bw(bw_one_pole_get_sticky_mode(&self.coeffs)) }
    }
    /// Returns the last output sample as stored in `OnePole<N_CHANNELS>`.
    pub fn get_y_z1(&self, channel: usize) -> f32 {
        unsafe { bw_one_pole_get_y_z1(&self.states[channel]) }
    }
    // Wrapping these to test them against the native ones
    #[cfg(test)]
    pub(crate) fn process1(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1(&mut self.coeffs, &mut self.states[channel], x) }
    }

    #[cfg(test)]
    pub(crate) fn process1_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_sticky_abs(&mut self.coeffs, &mut self.states[channel], x) }
    }

    #[cfg(test)]
    pub(crate) fn process1_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_sticky_rel(&mut self.coeffs, &mut self.states[channel], x) }
    }

    #[cfg(test)]
    pub(crate) fn process1_asym(&mut self, x: f32, channel: usize) -> f32 {
        unsafe { bw_one_pole_process1_asym(&mut self.coeffs, &mut self.states[channel], x) }
    }

    #[cfg(test)]
    pub(crate) fn process1_asym_sticky_abs(&mut self, x: f32, channel: usize) -> f32 {
        unsafe {
            bw_one_pole_process1_asym_sticky_abs(&mut self.coeffs, &mut self.states[channel], x)
        }
    }

    #[cfg(test)]
    pub(crate) fn process1_asym_sticky_rel(&mut self, x: f32, channel: usize) -> f32 {
        unsafe {
            bw_one_pole_process1_asym_sticky_rel(&mut self.coeffs, &mut self.states[channel], x)
        }
    }
}

impl OnePoleStickyMode {
    fn to_bw(&self) -> bw_one_pole_sticky_mode {
        match self {
            OnePoleStickyMode::Abs => bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs,
            OnePoleStickyMode::Rel => bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel,
        }
    }
    fn from_bw(bw_sticky_mode: bw_one_pole_sticky_mode) -> Self {
        match bw_sticky_mode {
            bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs => OnePoleStickyMode::Abs,
            bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_rel => OnePoleStickyMode::Rel,
            err_val => panic!("non valid bw_one_pole_sticky_mode (0, 1) got {err_val}")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f32;
    use std::f32::consts::PI;
    const N_CHANNELS: usize = 2;
    const BW_RCPF_ERROR: f32 = 0.0013;
    const SAMPLE_RATE: f32 = 48_000.0;
    const INVERSE_2_PI: f32 = 1.0 / (2.0 * PI);

    #[test]
    fn new() {
        let one_pole = OnePole::<N_CHANNELS>::new();

        assert_eq!(one_pole.coeffs.cutoff_up, f32::INFINITY);
        assert_eq!(one_pole.coeffs.cutoff_down, f32::INFINITY);
        assert_eq!(one_pole.coeffs.sticky_thresh, 0.);
        assert_eq!(
            one_pole.coeffs.sticky_mode,
            bw_one_pole_sticky_mode_bw_one_pole_sticky_mode_abs
        );
        assert_eq!(one_pole.coeffs.param_changed, !0); // from brickworks lib coeffs->param_changed = ~0; // useless, just to make compilers happy about uninitialized variables
    }

    #[test]
    fn set_sample_rate() {
        const SAMPLE_RATE: f32 = 48_000.0;
        let mut f = OnePole::<N_CHANNELS>::new();

        f.set_sample_rate(SAMPLE_RATE);

        assert_eq!(f.coeffs.fs_2pi, INVERSE_2_PI * SAMPLE_RATE)
    }

    #[test]
    fn set_cutoff() {
        const CUTOFF: f32 = 1000.0;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff(CUTOFF);

        assert_eq!(f.coeffs.cutoff_up, CUTOFF);
        assert_eq!(f.coeffs.cutoff_down, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_cutoff_up() {
        const CUTOFF: f32 = 1200.0;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_up(CUTOFF);

        assert_eq!(f.coeffs.cutoff_up, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_cutoff_down() {
        const CUTOFF: f32 = 1200.0;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_cutoff_down(CUTOFF);

        assert_eq!(f.coeffs.cutoff_down, CUTOFF);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau() {
        let cutoff = 150.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau(tau);

        let relative_error = (f.coeffs.cutoff_down - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR); // see bw_rcpf in bw_math.h
        assert_eq!(f.coeffs.cutoff_down, f.coeffs.cutoff_up);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_tau_up() {
        let cutoff = 1500.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau_up(tau);

        let relative_error = (f.coeffs.cutoff_up - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP != 0);
    }

    #[test]
    fn set_tau_down() {
        let cutoff = 10_000.0;
        let tau = INVERSE_2_PI / cutoff;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau_down(tau);

        let relative_error = (f.coeffs.cutoff_down - cutoff).abs() * 100.0 / cutoff;
        assert!(relative_error < BW_RCPF_ERROR);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN != 0);
    }

    #[test]
    fn set_tau_small_should_not_change_cutoff() {
        let tau = 0.1e-9;
        let mut f = OnePole::<N_CHANNELS>::new();
        f.coeffs.param_changed = 0;

        f.set_tau(tau);

        assert_eq!(f.coeffs.cutoff_up, f32::INFINITY);
        assert_eq!(f.coeffs.cutoff_down, f32::INFINITY);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_DOWN == 0);
        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_CUTOFF_UP == 0);
    }

    #[test]
    fn set_sticky_thresh() {
        let sticky_tresh = 0.01;
        let mut f = OnePole::<N_CHANNELS>::new();

        f.set_sticky_thresh(sticky_tresh);

        assert!(f.coeffs.param_changed as u32 & BW_ONE_POLE_PARAM_STICKY_THRESH != 0);
        assert_eq!(f.get_sticky_thresh(), sticky_tresh);
    }

    #[test]
    fn set_sticky_mode_abs() {
        let mode = OnePoleStickyMode::Abs;
        let mut f = OnePole::<N_CHANNELS>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[test]
    fn set_sticky_mode_rel() {
        let mode = OnePoleStickyMode::Rel;
        let mut f = OnePole::<N_CHANNELS>::new();

        f.set_sticky_mode(mode);

        assert_eq!(f.get_sticky_mode(), mode);
    }

    #[test]
    fn reset_none() {
        let cutoff = 1200.0;
        let sticky_tresh = 0.1;
        let x0_input = [0.5; N_CHANNELS];
        let mut f = OnePole::<N_CHANNELS>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_tresh);
        f.set_sample_rate(SAMPLE_RATE);
        let inverse = unsafe { f.coeffs.fs_2pi * bw_rcpf(f.coeffs.fs_2pi + cutoff) };

        f.reset(&x0_input, None);

        assert!(f.coeffs.param_changed == 0);
        assert_eq!(f.coeffs.mA1u, inverse);
        assert_eq!(f.coeffs.mA1d, inverse);
        assert_eq!(f.coeffs.st2, sticky_tresh * sticky_tresh);
    }

    #[test]
    fn reset_with_input_and_output() {
        let cutoff = 1000.0;
        let sticky_thresh = 0.2;
        let x0_input = [0.5; N_CHANNELS];
        let mut y0_output = [0.0; N_CHANNELS];
        let mut f = OnePole::<N_CHANNELS>::new();
        f.set_cutoff(cutoff);
        f.set_sticky_thresh(sticky_thresh);
        f.set_sample_rate(SAMPLE_RATE);

        f.reset(&x0_input, Some(&mut y0_output));

        for i in 0..N_CHANNELS {
            assert_eq!(f.states[i].y_z1, x0_input[i]);
            assert_eq!(y0_output[i], x0_input[i]);
        }
    }

    #[test]
    fn process_with_y() {
        const N_CHANNELS: usize = 2;
        const N_SAMPLES: usize = 4;
        let x0_input = [0.5; N_CHANNELS];

        let input_data = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut output_data: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.0, 0.0, 0.0], &mut [0.0, 0.0, 0.0, 0.0]];

        let mut filter = OnePole::<N_CHANNELS>::new();
        filter.set_cutoff(1000.0);
        filter.set_sample_rate(44100.0);
        filter.reset(&x0_input, None);

        filter.process(
            &input_data,
            Some(
                &mut output_data
                    .iter_mut()
                    .map(|ch| Some(&mut **ch))
                    .collect::<Vec<_>>(),
            ),
            N_SAMPLES,
        );

        for ch in 0..N_CHANNELS {
            for sample in 0..N_SAMPLES {
                println!("Output[{}][{}] = {}", ch, sample, output_data[ch][sample]);
                assert!(output_data[ch][sample].is_finite());
            }
        }
    }

    #[test]
    fn get_y_z1() {
        const N_CHANNELS: usize = 2;
        const SAMPLE_RATE: f32 = 44100.0;
        const CUTOFF: f32 = 1000.0;
        const N_SAMPLES: usize = 4;

        let mut f = OnePole::<N_CHANNELS>::new();
        f.set_sample_rate(SAMPLE_RATE);
        f.set_cutoff(CUTOFF);
        f.set_sticky_mode(OnePoleStickyMode::Rel);
        let input_data = [vec![1.0, 2.0, 3.0, 4.0], vec![0.5, 1.5, 2.5, 3.5]];
        let mut output_data: [&mut [f32]; N_CHANNELS] =
            [&mut [0.0, 0.1, 0.2, 0.3], &mut [1.0, 1.1, 1.2, 1.3]];
        // f.process(&input_data, Some(&mut output_data), N_SAMPLES);
        f.process(
            &input_data,
            Some(
                &mut output_data
                    .iter_mut()
                    .map(|ch| Some(&mut **ch))
                    .collect::<Vec<_>>(),
            ),
            N_SAMPLES,
        );

        for i in 0..N_CHANNELS {
            assert_eq!(f.get_y_z1(i), output_data[i][3])
        }
    }
}
