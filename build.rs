use std::env;
use std::path::PathBuf;
use std::process::Command;

use bindgen::Builder;

fn main() {
    let wrapper_path = "src/c_wrapper/wrapper.hpp";
    let current_dir = env::current_dir().expect("Error unwrapping current_dir()");
    let absolute_wrapper_path = current_dir
        .join(wrapper_path)
        .to_str()
        .expect("Error joining current dir with wrapper path.")
        .to_owned();
    let output_path = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let static_fns_path = output_path
        .join("extern.c")
        .to_str()
        .expect("Error joining $OUT_DIR with static fns.")
        .to_owned();
    // Tell bindgen to generate wrappers for static functions
    let bindings = Builder::default()
        .header(&absolute_wrapper_path)
        .clang_arg("-std=c++11")
        .generate_comments(false)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .wrap_static_fns(true)
        .wrap_static_fns_path("./extern.c")
        // .wrap_static_fns_path(&static_fns_path)
        .generate()
        .expect("Unable to generate bindings");

    // This is the path to the object file
    let obj_path = output_path.join("extern.o");
    // This is the path to the static library file.
    let lib_path = output_path.join("libextern.a");

    // Compile the generated wrappers into an object file
    let clang_output = std::process::Command::new("clang")
        .arg("-O")
        .arg("-c")
        .arg("-o")
        .arg(&obj_path)
        .arg("extern.c")
        // .arg(&static_fns_path)
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
            out_dir_path.join("libextern.lib").display()
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
    // Tell cargo to statically link against the `libextern` static library
    println!("cargo:rustc-link-lib=static=extern");

    // Write the rust bindings
    bindings
        .write_to_file(output_path.join("bindings.rs"))
        .expect("Cound not write bindings to the Rust file");
}
