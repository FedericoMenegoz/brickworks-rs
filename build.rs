use bindgen::Builder;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};
use std::collections::HashSet;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo::rerun-if-changed=brickworks/include");
    
    let wrapper_path = PathBuf::from("src/c_wrapper/wrapper.h")
        .canonicalize()
        .expect("Cannot canonicalize wrapper path");
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // paths for the static inline wrapper
    // where will be created the static library
    let static_filename = "static_inline_wrapper";
    let static_fns_path = out_dir.join(format!("{}.c", static_filename));

    // object file path (platform-specific extension)
    #[cfg(target_os = "windows")]
    let obj_path = out_dir.join(format!("{}.obj", static_filename));
    #[cfg(target_os = "windows")]
    let lib_path = out_dir.join(format!("{}.lib", static_filename));
    
    // library file path (platform-specific)
    #[cfg(not(target_os = "windows"))]
    let obj_path = out_dir.join(format!("{}.o", static_filename));
    #[cfg(not(target_os = "windows"))]
    let lib_path = out_dir.join(format!("lib{}.a", static_filename));
    
    // generate Rust bindings with bindgen
    let mut bind_build = Builder::default()
        .header(wrapper_path.to_str().unwrap())
        .clang_arg("-I./brickworks/include")
        .generate_comments(false)
        .parse_callbacks(Box::new(IgnoreMacros::new()))
        .derive_default(true)
        .wrap_static_fns(true)
        .wrap_static_fns_path(static_fns_path.to_str().unwrap());

    // prepare clang command for compilation
    let mut clang_cmd = Command::new("clang");
    clang_cmd
        .arg("-c")
        .arg("-o")
        .arg(&obj_path)
        .arg(&static_fns_path)
        .arg("-I./brickworks/include")
        .arg(format!("-I{}", wrapper_path.parent().unwrap().display()))
        .arg("-include")
        .arg(&wrapper_path);

    // debug and release mode configuration
    let is_release = std::env::var("PROFILE")
        .map(|p| p == "release")
        .unwrap_or(false);
    
    if !is_release {
        println!("cargo:warning=DEBUG MODE");
        bind_build = bind_build.clang_arg("-DBW_DEBUG_DEEP");
        clang_cmd.arg("-DBW_DEBUG_DEEP");
        clang_cmd.arg("-O0").arg("-g");
    } else {
        println!("cargo:warning=RELEASE MODE");
        bind_build = bind_build.clang_arg("-DBW_NO_DEBUG");
        clang_cmd.arg("-DBW_NO_DEBUG");

        // platform-specific optimization flags
        if cfg!(target_os = "windows") {
            println!("cargo:warning=Windows");
            clang_cmd
                .arg("-O3")
                .arg("-flto")
                .arg("-target")
                .arg("x86_64-pc-windows-msvc");
        } else if cfg!(target_os = "linux") {
            println!("cargo:warning=Linux");
            clang_cmd
                .arg("-O3")
                .arg("-flto=thin")
                .arg("-fPIC");
        } else if cfg!(target_os = "macos") {
            println!("cargo:warning=macOS");
            clang_cmd
                .arg("-O3")
                .arg("-flto")
                .arg("-march=native");
        }
    }

    // generate bindings
    let bindings = bind_build
        .generate()
        .expect("Unable to generate bindings");
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Could not write bindings file");
    // compile the static wrapper with clang
    let clang_output = clang_cmd.output().expect("Failed to run clang");
    
    if !clang_output.status.success() {
        panic!(
            "Could not compile object file:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&clang_output.stdout),
            String::from_utf8_lossy(&clang_output.stderr)
        );
    }

    // create static library
    #[cfg(target_os = "windows")]
    let lib_output = Command::new("llvm-lib")
        .arg(&obj_path)
        .arg(format!("/OUT:{}", lib_path.display()))
        .output()
        .expect("Failed to run llvm-lib. Make sure LLVM is in your PATH.");

    #[cfg(not(target_os = "windows"))]
    let lib_output = Command::new("ar")
        .arg("rcs")
        .arg(&lib_path)
        .arg(&obj_path)
        .output()
        .expect("Failed to run ar");

    if !lib_output.status.success() {
        panic!(
            "Could not create library file:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&lib_output.stdout),
            String::from_utf8_lossy(&lib_output.stderr)
        );
    }

    // tell cargo where to find the library
    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static={}", static_filename);
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