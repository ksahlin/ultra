//! Compile the vendored namfinder into the uLTRA binary.
//!
//! Mirrors upstream's `salib` CMake target source list, plus `main.cpp` (with
//! its `main` renamed out of the way) and `shim.cpp`. See
//! `vendor/namfinder/VENDORING.md`.
//!
//! Requirements: a C++17 compiler and zlib. No CMake, and nothing on PATH at
//! run time -- which is the point, since a missing namfinder is what makes the
//! reference uninstallable on osx-arm64.

use std::path::PathBuf;

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("vendor/namfinder");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Upstream generates these two with CMake's configure_file. They carry one
    // #define each, so we write them rather than depend on CMake.
    std::fs::write(
        out.join("buildconfig.hpp"),
        "#ifndef NAMFINDER_CONFIG_HPP\n#define NAMFINDER_CONFIG_HPP\n#define CMAKE_BUILD_TYPE \"Release\"\n#endif\n",
    )
    .expect("write buildconfig.hpp");
    std::fs::write(
        out.join("version.hpp"),
        "#ifndef NAMFINDER_VERSION_HPP\n#define NAMFINDER_VERSION_HPP\n#include <string>\n#define VERSION_STRING \"0.1.3\"\nstd::string version_string();\n#endif\n",
    )
    .expect("write version.hpp");

    // Exactly upstream's salib list.
    let cpp = [
        "refs", "fastq", "cmdline", "index", "indexparameters", "output", "pc",
        "aln", "nam", "randstrobes", "version", "io",
    ];

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++17")
        .include(root.join("src"))
        .include(root.join("ext"))
        .include(&out)
        .warnings(false);
    for f in cpp {
        build.file(root.join("src").join(format!("{f}.cpp")));
    }
    build.file(root.join("shim.cpp"));
    // namfinder's main.cpp holds run_strobealign; rename its `main` so it does
    // not collide with the Rust entry point.
    build.file(root.join("src/main.cpp")).define("main", "namfinder_cli_main");
    build.compile("namfinder");

    // xxhash is C, not C++, so it needs its own unit.
    cc::Build::new()
        .file(root.join("ext/xxhash.c"))
        .include(root.join("ext"))
        .warnings(false)
        .compile("xxhash");

    // edlib: one .cpp, vendored for the same reason namfinder is -- exact by
    // construction, and it settles the "which location comes first" tie-break
    // that help_functions.edlib_alignment depends on.
    let edlib = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("vendor/edlib");
    cc::Build::new()
        .cpp(true)
        .std("c++11")
        .file(edlib.join("edlib.cpp"))
        .include(&edlib)
        .warnings(false)
        .compile("edlib");
    println!("cargo:rerun-if-changed=vendor/edlib");

    println!("cargo:rustc-link-lib=z");
    println!("cargo:rerun-if-changed=vendor/namfinder");
}
