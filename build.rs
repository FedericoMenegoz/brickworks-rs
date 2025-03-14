use std::env;
use std::path::{Path, PathBuf};

fn main() {
    let src_path = Path::new("src/c_wrapper");

    let bindings = bindgen::Builder::default()
        .header(src_path.join("wrapper.h").as_os_str().to_str().unwrap())
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
