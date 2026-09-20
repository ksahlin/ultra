// C ABI over namfinder's own entry point.
//
// `run_strobealign` is the function namfinder's `main` calls. Going through it
// rather than reimplementing randstrobe seeding means the argument parsing,
// index construction, NAM finding and output formatting are all namfinder's
// own code -- exact by construction rather than exact by oracle.
//
// build.rs compiles namfinder's main.cpp with `-Dmain=namfinder_cli_main` so
// that its `main` does not collide with the Rust entry point while
// `run_strobealign` stays available.

extern int run_strobealign(int argc, char **argv);

extern "C" int nf_run(int argc, char **argv) {
    return run_strobealign(argc, argv);
}
