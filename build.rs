use std::path::PathBuf;
use std::process::Command;

use bindgen::Builder;

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
    let bindings = Builder::default()
        .header(&wrapper_path)
        .generate_comments(false)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .wrap_static_fns(true)
        .wrap_static_fns_path(&static_fns_path)
        .generate()
        .expect("Unable to generate bindings");

    // Compile the generated wrappers into an object file
    let clang_output = std::process::Command::new("clang")
        .arg("-O")
        .arg("-c")
        .arg("-o")
        .arg(&obj_path)
        .arg(&static_fns_path)
        .arg("-include")
        .arg(wrapper_path)
        .output()
        .unwrap();

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
            out_dir_path
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
