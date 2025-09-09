// Auto generated code does not follow these...
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
// To silence warnings "`extern` block uses type `u128`, which is not FFI-safe
// see comments here https://stackoverflow.com/questions/63526076/what-to-do-about-warning-extern-block-uses-type-u128-which-is-not-ffi-safe
#![allow(improper_ctypes)]

//! Rust bindings to C library [Orastron](https://www.orastron.com/algorithms)
//!
//!
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
#[allow(dead_code, unused_variables, unused_mut)]
pub mod clip;
pub mod one_pole;

// Helper functions
pub(crate) mod utils;
