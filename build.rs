use std::collections::HashSet;
use std::path::PathBuf;
use bindgen::Builder;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

fn main() {
    let wrapper_path = PathBuf::from("src/c_wrapper/wrapper.h")
        .canonicalize()
        .expect("Cannot canonicalize wrapper path");
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Paths for the static inline wrapper
    let static_filename = "static_inline_wrapper";
    let static_fns_path = out_dir.join(format!("{}.c", static_filename));

    // Generate Rust bindings with bindgen
    let mut builder = Builder::default()
        .header(wrapper_path.to_str().unwrap())
        .generate_comments(false)
        .parse_callbacks(Box::new(IgnoreMacros::new()))
        .derive_default(true)
        .wrap_static_fns(true)
        .wrap_static_fns_path(static_fns_path.to_str().unwrap());

    // Set debug/release macros
    if cfg!(debug_assertions) {
        builder = builder.clang_arg("-DBW_DEBUG_DEEP");
    } else {
        builder = builder.clang_arg("-DBW_NO_DEBUG");
    }

    let bindings = builder.generate().expect("Unable to generate bindings");
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Could not write bindings file");

    // Compile the generated C wrapper using cc crate
    let mut build = cc::Build::new();
    build.file(static_fns_path)
        .include(wrapper_path.parent().unwrap())
        .flag_if_supported("-O3")
        .flag_if_supported("-flto")
        .flag_if_supported("-fPIC");

    // Optional: add debug/release macro defines for C code
    if cfg!(debug_assertions) {
        build.define("BW_DEBUG_DEEP", None);
    } else {
        build.define("BW_NO_DEBUG", None);
    }

    // Compile as a static library
    build.compile(static_filename);

    // Tell Cargo where to find the library (handled automatically by cc)
}

// Ignore certain macros to avoid bindgen parsing errors
const IGNORE_MACROS: [&str; 5] = [
    "FP_INFINITE",
    "FP_NAN",
    "FP_NORMAL",
    "FP_SUBNORMAL",
    "FP_ZERO",
];

#[derive(Debug)]
struct IgnoreMacros(HashSet<String>);

impl ParseCallbacks for IgnoreMacros {
    fn will_parse_macro(&self, name: &str) -> MacroParsingBehavior {
        if self.0.contains(name) {
            MacroParsingBehavior::Ignore
        } else {
            MacroParsingBehavior::Default
        }
    }
}

impl IgnoreMacros {
    fn new() -> Self {
        Self(IGNORE_MACROS.into_iter().map(|s| s.to_owned()).collect())
    }
}
