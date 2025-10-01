use bindgen::Builder;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};
use std::collections::HashSet;
use std::path::PathBuf;

fn main() {
    let wrapper_path = PathBuf::from("src/c_wrapper/wrapper.h");
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Paths for the static inline wrapper
    let static_filename = "static_inline_wrapper";
    let static_fns_path = out_dir.join(format!("{}.c", static_filename));

    // Generate Rust bindings with bindgen
    let mut bind_build = Builder::default()
        .header(wrapper_path.to_str().unwrap())
        .clang_arg("-I./brickworks/include")
        .generate_comments(false)
        .parse_callbacks(Box::new(IgnoreMacros::new()))
        .derive_default(true)
        .wrap_static_fns(true)
        .wrap_static_fns_path(static_fns_path.to_str().unwrap());

    // Compile the generated C wrapper using cc crate
    let mut compiler_build = cc::Build::new();
    compiler_build.file(&static_fns_path).warnings(false);

    compiler_build
        .include(wrapper_path.parent().unwrap())
        // Windows so that it will find the files in wrapper.h
        .include(".");

    // Test and debug mode
    if cfg!(debug_assertions) {
        bind_build = bind_build.clang_arg("-DBW_DEBUG_DEEP");
        compiler_build.define("BW_DEBUG_DEEP", None);
    }
    // Optimization for release and benchmark
    else {
        // Windows msvc toolchain
        if cfg!(target_env = "msvc") {
            compiler_build
                .flag("/O2") // Maximum optimization (speed)
                .flag("/GL"); // Whole program optimization
        } else {
            compiler_build
                .flag_if_supported("-O3")
                .flag_if_supported("-flto")
                .flag_if_supported("-fPIC");
            // .flag_if_supported("-march=native");
        }
        bind_build = bind_build.clang_arg("-DBW_NO_DEBUG");
        compiler_build.define("BW_NO_DEBUG", None);
    }

    let bindings = bind_build.generate().expect("Unable to generate bindings");
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Could not write bindings file");
    // Compile as a static library
    compiler_build.compile(static_filename);
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
