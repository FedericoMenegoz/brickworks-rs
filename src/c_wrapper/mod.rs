#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

// To silence warnings "`extern` block uses type `u128`, which is not FFI-safe"
// see comments here https://stackoverflow.com/questions/63526076/what-to-do-about-warning-extern-block-uses-type-u128-which-is-not-ffi-safe
// #![allow(improper_ctypes)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub fn prova() {
    unsafe {
        println!("{}",__nan())
    }
}

#[cfg(test)]
mod test{
    use crate::c_wrapper::__nan;


    #[test]
    pub fn dummyTestJustChecking() {
        unsafe {
            println!("{}",__nan())
        }
        assert!(true)
    }
}