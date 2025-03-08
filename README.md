# Brickworks-rs 

Brickworks-rs is a Rust port of the [Brickworks](https://www.orastron.com/algorithms#collections-bundles) audio DSP library, designed for real-time audio processing. This project provides both a pure Rust implementation and an FFI wrapper around the original C library to compare performance and correctness.

## Installing

## Testing

## Docs

## Roadmap

### Distortion effect [bw_dist](https://www.orastron.com/algorithms/bw_dist)
The submodules below are listed in topological order based on their dependency graph.

#### Common definitions `bw_common`
- [ ] `bw_common` wrapper
- [ ] `bw_common` port skeleton and TDD
- [ ] `bw_common` port implementation
- [ ] `bw_common` testing

#### Mathematical routines `bw_math`
- [ ] `bw_math` wrapper
- [ ] `bw_math` port skeleton and TDD
- [ ] `bw_math` port implementation
- [ ] `bw_math` testing

#### One-pole lowpass filter `bw_one_pole`
- [ ] `bw_one_pole` wrapper
- [ ] `bw_one_pole` port skeleton and TDD
- [ ] `bw_one_pole` port implementation
- [ ] `bw_one_pole` testing

#### Antialiased hard clipper `bw_clip`
- [ ] `bw_clip` wrapper
- [ ] `bw_clip` port skeleton and TDD
- [ ] `bw_clip` port implementation
- [ ] `bw_clip` testing

#### Gain module `bw_gain`
- [ ] `bw_gain` wrapper
- [ ] `bw_gain` port skeleton and TDD
- [ ] `bw_gain` port implementation
- [ ] `bw_gain` testing

#### Antialiased tanh-based saturation `bw_satur`
- [ ] `bw_satur` wrapper
- [ ] `bw_satur` port skeleton and TDD
- [ ] `bw_satur` port implementation
- [ ] `bw_satur` testing

#### First-order lowpass filter `bw_lp1`
- [ ] `bw_lp1` wrapper
- [ ] `bw_lp1` port skeleton and TDD
- [ ] `bw_lp1` port implementation
- [ ] `bw_lp1` testing

#### First-order highpass filter `bw_hp1`
- [ ] `bw_hp1` wrapper
- [ ] `bw_hp1` port skeleton and TDD
- [ ] `bw_hp1` port implementation
- [ ] `bw_hp1` testing

#### State variable filter `bw_svf`
- [ ] `bw_svf` wrapper
- [ ] `bw_svf` port skeleton and TDD
- [ ] `bw_svf` port implementation
- [ ] `bw_svf` testing

#### Second-order multimode filter `bw_mm2`
- [ ] `bw_mm2` wrapper
- [ ] `bw_mm2` port skeleton and TDD
- [ ] `bw_mm2` port implementation
- [ ] `bw_mm2` testing

#### Second-order peak filter `bw_peak`
- [ ] `bw_peak` wrapper
- [ ] `bw_peak` port skeleton and TDD
- [ ] `bw_peak` port implementation
- [ ] `bw_peak` testing

#### Distortion `bw_dist`
- [ ] `bw_dist` wrapper
- [ ] `bw_dist` port skeleton and TDD
- [ ] `bw_dist` port implementation
- [ ] `bw_dist` testing

## License
Brickworks-rs is distributed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) License.
