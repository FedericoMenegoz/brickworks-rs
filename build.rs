use bindgen::Builder;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};
use std::collections::HashSet;
use std::path::PathBuf;

fn main() {
    println!("cargo::rerun-if-changed=brickworks/include");
    let wrapper_path = PathBuf::from("src/c_wrapper/wrapper.h");
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // paths for the static inline wrapper
    let static_filename = "static_inline_wrapper";
    let static_fns_path = out_dir.join(format!("{}.c", static_filename));

    // generate Rust bindings with bindgen
    let mut bind_build = Builder::default()
        .header(wrapper_path.to_str().unwrap())
        .clang_arg("-I./brickworks/include")
        .generate_comments(false)
        .parse_callbacks(Box::new(IgnoreMacros::new()))
        .derive_default(true)
        .wrap_static_fns(true)
        .wrap_static_fns_path(static_fns_path.to_str().unwrap());

    // compile the generated C wrapper using cc crate
    let mut compiler_build = cc::Build::new();
    compiler_build.file(&static_fns_path).warnings(false);

    compiler_build
        .include(wrapper_path.parent().unwrap())
        // Windows so that it will find the files in wrapper.h
        .include(".");

    // test and debug mode
    if cfg!(debug_assertions) {
        println!("cargo:warning=debug assertion mode");
        bind_build = bind_build.clang_arg("-DBW_DEBUG_DEEP");
        compiler_build.define("BW_DEBUG_DEEP", None);
    }
    // optimization for release and benchmark
    else {
        println!("cargo:warning=release mode");
        // windows
        if cfg!(target_env = "msvc") {
            println!("cargo:warning=msvc");
            compiler_build.flag("/O2");
        }
        // linux
        else if cfg!(target_os = "linux") {
            println!("cargo:warning=linux");
            compiler_build
                .flag_if_supported("-O")
                .flag_if_supported("-lto=thin")
                .flag_if_supported("-fPIC");
        }
        // macos
        else if cfg!(target_os = "macos") {
            println!("cargo:warning=macos");
            compiler_build
                .flag_if_supported("-O3")
                .flag_if_supported("-flto")
                .flag_if_supported("-march=native")
                .flag_if_supported("-ffast-math")
                .flag_if_supported("-funroll-loops")
                .flag_if_supported("-fomit-frame-pointer")
                .flag_if_supported("-fstrict-aliasing");
        }
        bind_build = bind_build.clang_arg("-DBW_NO_DEBUG");
        compiler_build.define("BW_NO_DEBUG", None);
    }

    let bindings = bind_build.generate().expect("Unable to generate bindings");
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Could not write bindings file");
    // compile as a static library
    compiler_build.compile(static_filename);
}

// ignore certain macros to avoid bindgen parsing errors
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
