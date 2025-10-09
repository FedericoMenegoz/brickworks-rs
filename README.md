# Brickworks-rs 

Brickworks-rs is a Rust port of the [Brickworks](https://www.orastron.com/algorithms#collections-bundles) audio DSP library, designed for real-time audio processing. The aim of this project is to provide both a pure Rust implementation and an FFI wrapper around the original C library to compare performance and correctness. 

Check the Project Wiki [here](https://github.com/FedericoMenegoz/brickworks-rs/wiki).

# Building / Benchmarking in release mode (LTO optimizations)
## Linux and macOS
Use environment variables to set Rust compilation flags only for release builds or benchmarks.
```bash
# build --release
RUSTFLAGS="-Clinker-plugin-lto -Clinker=clang -Clink-arg=-fuse-ld=lld" cargo build --release
# bench
RUSTFLAGS="-Clinker-plugin-lto -Clinker=clang -Clink-arg=-fuse-ld=lld" cargo bench
```
## macOS specific error
If you get an error like:
```
error: linking with `clang` failed: exit status: 1
  |
  = note: some arguments are omitted. use `--verbose` to show all linker arguments
  = note: ld64.lld: error: unknown argument '-plugin-opt=O0'
          ld64.lld: error: unknown argument '-plugin-opt=mcpu=penryn'
          clang: error: linker command failed with exit code 1 (use -v to see invocation)
```
This is a known issue with `lld` on macOS. A workaround and discussion can be found here: [rust-lang/rust issue #60059]((https://github.com/rust-lang/rust/issues/60059)).

# Windows
On windows in order to use LTO you need to use `lld-link`:

Use default MSVC linker:
```bash
# Build in release mode
$env:RUSTFLAGS = "-C linker=lld-link" cargo build --release

# Run benchmarks
$env:RUSTFLAGS = "-C linker=lld-link" cargo bench
```

## Example of use
[Here](https://github.com/FedericoMenegoz/Brickworks-rs-plugin) you can find a minimal example of how to build a VST3 and CLAP plugin using [nih_plugin](https://github.com/robbert-vdh/nih-plug.git).

## License
Brickworks-rs is distributed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) License.
