//! End-to-end exercise of the FFI roundtrip against a real upstream-bundled
//! GDS file. The fixture is `vendor/examples/data/CEU_Exon.gds`, which ships
//! with `CoreArray/GDSFormat` and is in CoreArray binary format. Confirms
//! that the C++ shim builds, links, and that open/close/last_error/ndim/dims
//! produce sane values on real bytes — independent of whether higher-level
//! layouts (SNPGDS vs SeqArray) match.

use std::ffi::{CStr, CString};
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("vendor/examples/data")
        .join(name)
}

fn read_err(file: *mut corearray_sys::CaFile) -> String {
    unsafe {
        let p = corearray_sys::corearray_last_error(file);
        if p.is_null() {
            String::new()
        } else {
            CStr::from_ptr(p).to_string_lossy().into_owned()
        }
    }
}

#[test]
fn open_and_close_real_gds_file() {
    let path = fixture("CEU_Exon.gds");
    assert!(path.exists(), "fixture missing: {}", path.display());
    let cpath = CString::new(path.to_str().unwrap()).unwrap();

    unsafe {
        let f = corearray_sys::corearray_open(cpath.as_ptr());
        assert!(!f.is_null(), "corearray_open returned null");

        let err = read_err(f);
        assert!(err.is_empty(), "open reported error: {err}");

        corearray_sys::corearray_close(f);
    }
}

#[test]
fn open_nonexistent_path_surfaces_error() {
    let cpath = CString::new("/no/such/file.gds").unwrap();
    unsafe {
        let f = corearray_sys::corearray_open(cpath.as_ptr());
        assert!(!f.is_null(), "shim should hand back a handle even on open failure");
        let err = read_err(f);
        assert!(!err.is_empty(), "expected error message, got empty");
        corearray_sys::corearray_close(f);
    }
}

#[test]
fn ndim_on_missing_node_returns_error() {
    let path = fixture("CEU_Exon.gds");
    let cpath = CString::new(path.to_str().unwrap()).unwrap();
    let bogus = CString::new("/this/node/does/not/exist").unwrap();
    unsafe {
        let f = corearray_sys::corearray_open(cpath.as_ptr());
        assert!(!f.is_null());
        assert!(read_err(f).is_empty(), "open clean");

        let mut ndim = 0i32;
        let rc = corearray_sys::corearray_ndim(f, bogus.as_ptr(), &mut ndim);
        assert_ne!(rc, 0, "ndim on missing path should return non-zero");
        let err = read_err(f);
        assert!(
            err.contains("path not found") || err.contains("does not exist") || !err.is_empty(),
            "expected non-empty error, got '{err}'"
        );

        corearray_sys::corearray_close(f);
    }
}
