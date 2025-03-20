#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

// To silence warnings "`extern` block uses type `u128`, which is not FFI-safe
// see comments here https://stackoverflow.com/questions/63526076/what-to-do-about-warning-extern-block-uses-type-u128-which-is-not-ffi-safe
// #![allow(improper_ctypes)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod test {
    use crate::c_wrapper::*;

    #[test]
    pub fn dummyTestJustChecking() {
        unsafe {
            bw_is_inf(3.0);
        }

        assert!(true)
    }
}
