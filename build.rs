use std::collections::HashSet;
use std::path::PathBuf;
use std::process::Command;

use bindgen::Builder;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

fn main() {
    let wrapper_path = PathBuf::from("src/c_wrapper/wrapper.h")
        .canonicalize()
        .expect("Cannot canonicalize wrapper path")
        .to_str()
        .expect("The wrapper path is not a valid string.")
        .to_owned();
    let output_path = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Path for the static and inline function
    let static_filename = "static_inline_wrapper";
    let static_fns_path = output_path
        .join(static_filename.to_owned() + ".c")
        .to_str()
        .expect("Error joining $OUT_DIR with static fns.")
        .to_owned();

    // This is the path to the object file
    let obj_path = output_path.join(static_filename.to_owned() + ".o");
    // This is the path to the static library file.
    let lib_path = output_path.join("lib".to_owned() + static_filename + ".a");

    // Build the bindings
    let mut builder = Builder::default()
        .header(&wrapper_path)
        .generate_comments(false)
        .parse_callbacks(Box::new(IgnoreMacros::new()))
        .derive_default(true)
        .wrap_static_fns(true)
        .wrap_static_fns_path(&static_fns_path);

    // Prepare clang command for compilation
    let mut clang_cmd = std::process::Command::new("clang");
    clang_cmd
        .arg("-O")
        .arg("-c")
        .arg("-fPIC")
        .arg("-o")
        .arg(&obj_path)
        .arg(&static_fns_path)
        .arg("-include")
        .arg(wrapper_path);

    // To define DEBUG_DEEP when not in release
    if cfg!(debug_assertions) {
        println!("cargo:warning=INFO: DEBUG_DEEP is defined.");
        // Define BW_DEBUG_DEEP for bindgen
        builder = builder.clang_arg("-DBW_DEBUG_DEEP");
        // Define BW_DEBUG_DEEP for clang compilation
        clang_cmd.arg("-D").arg("BW_DEBUG_DEEP");
    }

    let bindings = builder.generate().expect("Unable to generate bindings");

    let clang_output = clang_cmd.output().unwrap();

    if !clang_output.status.success() {
        panic!(
            "Could not compile object file:\n{}",
            String::from_utf8_lossy(&clang_output.stderr)
        );
    }

    // Turn the object file into a static library
    #[cfg(not(target_os = "windows"))]
    let lib_output = Command::new("ar")
        .arg("rcs")
        .arg(lib_path)
        .arg(obj_path)
        .output()
        .unwrap();
    #[cfg(target_os = "windows")]
    let lib_output = Command::new("LIB")
        .arg(obj_path)
        .arg(format!(
            "/OUT:{}",
            output_path
                .join("lib".to_owned() + static_filename + ".lib")
                .display()
        ))
        .output()
        .unwrap();
    if !lib_output.status.success() {
        panic!(
            "Could not emit library file:\n{}",
            String::from_utf8_lossy(&lib_output.stderr)
        );
    }
    // So that cargo can find the brickworks wrapped library for inline static
    println!("cargo:rustc-link-search=native={}", output_path.display());
    // Tell cargo to statically link against the static wrapper
    println!("cargo:rustc-link-lib=static={}", static_filename);

    // Write the rust bindings
    bindings
        .write_to_file(output_path.join("bindings.rs"))
        .expect("Cound not write bindings to the Rust file");
}

// See https://github.com/rust-lang/rust-bindgen/issues/687#issuecomment-450750547
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
