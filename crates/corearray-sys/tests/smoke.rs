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

/// Ignored by default. Set `FAVOR_AGDS_PROBE_PATH` to a real production
/// SeqArray (a)GDS file (e.g. one of the 136k CCDG chr*.agds files) and
/// run with `--ignored` to probe the FFI against real data.
///
/// Verifies open + reads `/sample.id` count, `/position` count, and the
/// first three positions. No annotation channels are touched.
#[test]
#[ignore]
fn probe_real_agds_from_env() {
    let path = match std::env::var("FAVOR_AGDS_PROBE_PATH") {
        Ok(p) => p,
        Err(_) => panic!("set FAVOR_AGDS_PROBE_PATH to a real .agds file"),
    };
    let cpath = CString::new(path.as_str()).unwrap();
    let sample_id = CString::new("/sample.id").unwrap();
    let position = CString::new("/position").unwrap();
    let chromosome = CString::new("/chromosome").unwrap();
    let allele = CString::new("/allele").unwrap();
    let geno = CString::new("/genotype/data").unwrap();

    unsafe {
        let f = corearray_sys::corearray_open(cpath.as_ptr());
        assert!(!f.is_null());
        let err = read_err(f);
        assert!(err.is_empty(), "open '{path}' reported: {err}");

        // /sample.id is a 1-d string array.
        let mut nd = 0i32;
        assert_eq!(
            corearray_sys::corearray_ndim(f, sample_id.as_ptr(), &mut nd),
            0,
            "/sample.id ndim: {}",
            read_err(f)
        );
        assert_eq!(nd, 1, "/sample.id should be 1-d");
        let mut sample_dims = [0i64];
        assert_eq!(
            corearray_sys::corearray_dims(f, sample_id.as_ptr(), 1, sample_dims.as_mut_ptr()),
            0,
            "/sample.id dims: {}",
            read_err(f)
        );
        let n_samples = sample_dims[0];
        eprintln!("/sample.id: {n_samples} samples");
        assert!(n_samples > 0, "expected at least one sample");

        // /chromosome and /position are 1-d, length == n_variants.
        let mut pos_dims = [0i64];
        assert_eq!(
            corearray_sys::corearray_dims(f, position.as_ptr(), 1, pos_dims.as_mut_ptr()),
            0,
            "/position dims: {}",
            read_err(f)
        );
        let n_variants = pos_dims[0];
        eprintln!("/position: {n_variants} variants");
        assert!(n_variants > 0);

        // /genotype/data is 3-d: (ploidy, sample, variant).
        assert_eq!(
            corearray_sys::corearray_ndim(f, geno.as_ptr(), &mut nd),
            0,
            "/genotype/data ndim: {}",
            read_err(f)
        );
        eprintln!("/genotype/data ndim: {nd}");
        assert!(nd == 2 || nd == 3, "expected 2- or 3-d genotype array");

        let mut gd = vec![0i64; nd as usize];
        assert_eq!(
            corearray_sys::corearray_dims(f, geno.as_ptr(), nd, gd.as_mut_ptr()),
            0,
            "/genotype/data dims: {}",
            read_err(f)
        );
        eprintln!("/genotype/data dims: {gd:?}");

        // Read first three /position entries.
        let count = 3.min(n_variants as i32);
        let mut pos = vec![0i32; count as usize];
        let s = [0i32];
        let l = [count];
        assert_eq!(
            corearray_sys::corearray_read_int32(
                f,
                position.as_ptr(),
                1,
                s.as_ptr(),
                l.as_ptr(),
                pos.as_mut_ptr(),
            ),
            0,
            "read_int32 /position: {}",
            read_err(f)
        );
        eprintln!("first {count} positions: {pos:?}");

        // SeqArray /genotype/data is 3-d (variant, sample, ploidy) with
        // dims = [n_variants, n_samples, 2]. Read one variant: pin axis 0
        // to 0, take all samples and ploidies. 2-d fallback (sample, variant).
        let (geno_count, gs, gl): (usize, Vec<i32>, Vec<i32>) = if nd == 3 {
            let nsamp = gd[1];
            let ploidy = gd[2];
            ((nsamp * ploidy) as usize, vec![0, 0, 0], vec![1, nsamp as i32, ploidy as i32])
        } else {
            let nsamp = gd[0];
            (nsamp as usize, vec![0, 0], vec![nsamp as i32, 1])
        };
        let mut buf = vec![0i8; geno_count];
        assert_eq!(
            corearray_sys::corearray_read_int8(
                f,
                geno.as_ptr(),
                nd,
                gs.as_ptr(),
                gl.as_ptr(),
                buf.as_mut_ptr(),
            ),
            0,
            "read_int8 /genotype/data: {}",
            read_err(f)
        );
        let nonneg = buf.iter().filter(|&&c| c >= 0).count();
        let neg = buf.iter().filter(|&&c| c < 0).count();
        eprintln!("variant 0 dosage: {nonneg} non-missing alleles, {neg} missing");

        // /chromosome and /allele are 1-d string arrays. Read first one.
        let dummy_cb: corearray_sys::corearray_str_cb = {
            extern "C" fn cb(
                _u: *mut std::ffi::c_void,
                idx: i32,
                ptr: *const std::os::raw::c_char,
                len: usize,
            ) {
                let s = unsafe {
                    std::str::from_utf8(std::slice::from_raw_parts(ptr as *const u8, len))
                        .unwrap_or("?")
                };
                eprintln!("  string[{idx}] = {s:?}");
            }
            cb
        };
        eprintln!("first /chromosome:");
        assert_eq!(
            corearray_sys::corearray_read_string(
                f,
                chromosome.as_ptr(),
                0,
                1,
                dummy_cb,
                std::ptr::null_mut(),
            ),
            0,
            "read /chromosome: {}",
            read_err(f)
        );
        eprintln!("first /allele:");
        assert_eq!(
            corearray_sys::corearray_read_string(
                f,
                allele.as_ptr(),
                0,
                1,
                dummy_cb,
                std::ptr::null_mut(),
            ),
            0,
            "read /allele: {}",
            read_err(f)
        );

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
