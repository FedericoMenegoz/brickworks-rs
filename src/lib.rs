//! # Brickworks-rs
//!
//! `brickworks-rs` is a Rust port of the original [Brickworks C library](https://www.orastron.com/algorithms), a collection of DSP building blocks for real-time audio processing.
//!
//! ---
//! This crate provides both:
//! - A **native** Rust implementation of the algorithms
//! - A **wrapper** around the original C code using `bindgen`
//!
//! Both implementations expose the same Rust interface.
//!
//! ### See also
//! - [example of use](https://github.com/FedericoMenegoz/Brickworks-rs-plugin)
//! - [project wiki](https://github.com/FedericoMenegoz/brickworks-rs/wiki/Index)
// Silence clippy warnings for the wrapper
#[allow(clippy::all)]
// Silence clippy warnings for the wrapper
pub mod c_wrapper;
pub mod native;
