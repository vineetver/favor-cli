//! GDS (SeqArray / CoreArray) variant reader.
//!
//! Reads /sample.id, /chromosome, /position, /allele once at open time. Then
//! streams one variant at a time, decoding the bit2-packed /genotype/data
//! slice into a per-variant samples_text buffer that mirrors the VCF reader's
//! contract (FORMAT-prefixed, tab-separated, "GT" column only).
//!
//! /genotype/data is a 3-d int8 array with axes (ploidy, sample, variant).
//! Per-allele codes are 0=ref, 1..k=alt_k, anything < 0 is missing
//! (CoreArray emits the bit-pattern's signed extension; we treat any negative
//! value as missing). Diploid GT is built as "a/b" using these codes; "." for
//! missing alleles.

#![allow(dead_code)] // GdsVariantReader is wired via FormatHandler trait object

use std::ffi::{c_void, CStr, CString};
use std::fmt::Write as _;
use std::os::raw::c_char;
use std::path::Path;

use crate::error::CohortError;

use super::reader::{RawRecord, VariantReader};

/// Compile-time CStr from a byte literal with explicit NUL terminator.
const fn c(bytes: &'static [u8]) -> &'static CStr {
    // SAFETY: caller passes a literal that ends with `\0` and contains no
    // interior NULs.
    unsafe { CStr::from_bytes_with_nul_unchecked(bytes) }
}

const PATH_SAMPLE_ID: &CStr = c(b"/sample.id\0");
const PATH_POSITION: &CStr = c(b"/position\0");
const PATH_CHROMOSOME: &CStr = c(b"/chromosome\0");
const PATH_ALLELE: &CStr = c(b"/allele\0");
const PATH_GENOTYPE_DATA: &CStr = c(b"/genotype/data\0");

pub struct GdsVariantReader {
    file: GdsHandle,
    sample_count: i32,
    variant_count: i32,
    ploidy: i32,
    chromosome: Vec<String>,
    position: Vec<i32>,
    allele: Vec<String>,
    geno_buf: Vec<i8>,
    samples_text: String,
}

impl GdsVariantReader {
    pub fn open(path: &Path) -> Result<Self, CohortError> {
        let file = GdsHandle::open(path)?;

        let sample_count = file.array_dim(PATH_SAMPLE_ID, 0)? as i32;
        let variant_count = file.array_dim(PATH_POSITION, 0)? as i32;

        let geno_dims = file.array_dims(PATH_GENOTYPE_DATA, 3)?;
        let ploidy = geno_dims[0] as i32;
        if geno_dims[1] as i32 != sample_count || geno_dims[2] as i32 != variant_count {
            return Err(CohortError::Input(format!(
                "GDS '{}': /genotype/data dims {:?} disagree with sample_count={} variant_count={}",
                path.display(),
                geno_dims,
                sample_count,
                variant_count
            )));
        }

        let chromosome = file.read_string_array_1d(PATH_CHROMOSOME, 0, variant_count)?;
        let position = file.read_int32_array_1d(PATH_POSITION, 0, variant_count)?;
        let allele = file.read_string_array_1d(PATH_ALLELE, 0, variant_count)?;

        let geno_buf = vec![0i8; (ploidy as usize) * (sample_count as usize)];
        let samples_text =
            String::with_capacity(2 + (sample_count as usize) * (ploidy as usize * 2 + 1));

        Ok(Self {
            file,
            sample_count,
            variant_count,
            ploidy,
            chromosome,
            position,
            allele,
            geno_buf,
            samples_text,
        })
    }

    pub fn sample_count(&self) -> i32 {
        self.sample_count
    }

    pub fn sample_names(&self) -> Result<Vec<String>, CohortError> {
        self.file
            .read_string_array_1d(PATH_SAMPLE_ID, 0, self.sample_count)
    }
}

impl VariantReader for GdsVariantReader {
    fn for_each(
        &mut self,
        f: &mut dyn for<'a> FnMut(RawRecord<'a>) -> Result<(), CohortError>,
    ) -> Result<(), CohortError> {
        for v in 0..self.variant_count {
            let start = [0i32, 0, v];
            let length = [self.ploidy, self.sample_count, 1];
            self.file
                .read_int8_slice(PATH_GENOTYPE_DATA, &start, &length, &mut self.geno_buf)?;

            self.samples_text.clear();
            self.samples_text.push_str("GT");
            for s in 0..(self.sample_count as usize) {
                self.samples_text.push('\t');
                for p in 0..(self.ploidy as usize) {
                    if p > 0 {
                        self.samples_text.push('/');
                    }
                    let code = self.geno_buf[p * (self.sample_count as usize) + s];
                    if code < 0 {
                        self.samples_text.push('.');
                    } else {
                        let _ = write!(self.samples_text, "{}", code);
                    }
                }
            }

            let alleles = self.allele[v as usize].as_str();
            let (ref_allele, alt_alleles) = match alleles.find(',') {
                Some(i) => (&alleles[..i], &alleles[i + 1..]),
                None => (alleles, ""),
            };

            let rec = RawRecord {
                chromosome: self.chromosome[v as usize].as_str(),
                position: self.position[v as usize],
                ref_allele,
                alt_alleles,
                rsid: None,
                qual: None,
                filter: None,
                samples_text: self.samples_text.as_str(),
            };
            f(rec)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Thin RAII wrapper over the corearray-sys FFI. All entry points return
// Result<_, CohortError> with the upstream error string preserved.

struct GdsHandle {
    inner: *mut corearray_sys::CaFile,
    path: String,
}

unsafe impl Send for GdsHandle {}

impl GdsHandle {
    fn open(path: &Path) -> Result<Self, CohortError> {
        let cpath = CString::new(path.as_os_str().to_string_lossy().as_bytes())
            .map_err(|e| CohortError::Input(format!("path contains NUL: {e}")))?;
        let inner = unsafe { corearray_sys::corearray_open(cpath.as_ptr()) };
        if inner.is_null() {
            return Err(CohortError::Input(format!(
                "GDS open: allocation failed for '{}'",
                path.display()
            )));
        }
        let h = GdsHandle {
            inner,
            path: path.display().to_string(),
        };
        // CoreArray surfaces open failures via last_error rather than a null handle.
        let err = h.last_error();
        if !err.is_empty() {
            return Err(CohortError::Input(format!(
                "GDS open '{}': {err}",
                path.display()
            )));
        }
        Ok(h)
    }

    fn last_error(&self) -> String {
        unsafe {
            let p = corearray_sys::corearray_last_error(self.inner);
            if p.is_null() {
                String::new()
            } else {
                std::ffi::CStr::from_ptr(p).to_string_lossy().into_owned()
            }
        }
    }

    fn err<T>(&self, where_: &str) -> Result<T, CohortError> {
        Err(CohortError::Input(format!(
            "GDS {where_} on '{}': {}",
            self.path,
            self.last_error()
        )))
    }

    fn array_dim(&self, path: &std::ffi::CStr, axis: usize) -> Result<i64, CohortError> {
        let mut ndim: i32 = 0;
        let rc = unsafe { corearray_sys::corearray_ndim(self.inner, path.as_ptr(), &mut ndim) };
        if rc != 0 {
            return self.err("ndim");
        }
        if axis as i32 >= ndim {
            return Err(CohortError::Input(format!(
                "GDS '{}': axis {axis} out of range for ndim={ndim} on {}",
                self.path,
                path.to_string_lossy()
            )));
        }
        let mut dims = vec![0i64; ndim as usize];
        let rc = unsafe {
            corearray_sys::corearray_dims(self.inner, path.as_ptr(), ndim, dims.as_mut_ptr())
        };
        if rc != 0 {
            return self.err("dims");
        }
        Ok(dims[axis])
    }

    fn array_dims(
        &self,
        path: &std::ffi::CStr,
        expected_ndim: i32,
    ) -> Result<Vec<i64>, CohortError> {
        let mut ndim: i32 = 0;
        let rc = unsafe { corearray_sys::corearray_ndim(self.inner, path.as_ptr(), &mut ndim) };
        if rc != 0 {
            return self.err("ndim");
        }
        if ndim != expected_ndim {
            return Err(CohortError::Input(format!(
                "GDS '{}': '{}' has ndim={ndim}, expected {expected_ndim}",
                self.path,
                path.to_string_lossy()
            )));
        }
        let mut dims = vec![0i64; ndim as usize];
        let rc = unsafe {
            corearray_sys::corearray_dims(self.inner, path.as_ptr(), ndim, dims.as_mut_ptr())
        };
        if rc != 0 {
            return self.err("dims");
        }
        Ok(dims)
    }

    fn read_int32_array_1d(
        &self,
        path: &std::ffi::CStr,
        start: i32,
        length: i32,
    ) -> Result<Vec<i32>, CohortError> {
        let mut buf = vec![0i32; length as usize];
        let s = [start];
        let l = [length];
        let rc = unsafe {
            corearray_sys::corearray_read_int32(
                self.inner,
                path.as_ptr(),
                1,
                s.as_ptr(),
                l.as_ptr(),
                buf.as_mut_ptr(),
            )
        };
        if rc != 0 {
            return self.err("read_int32");
        }
        Ok(buf)
    }

    fn read_int8_slice(
        &self,
        path: &std::ffi::CStr,
        start: &[i32],
        length: &[i32],
        out: &mut [i8],
    ) -> Result<(), CohortError> {
        let need: usize = length.iter().map(|&n| n as usize).product();
        if out.len() < need {
            return Err(CohortError::Input(format!(
                "GDS read_int8: buffer of {} too small for {need} elements",
                out.len()
            )));
        }
        let rc = unsafe {
            corearray_sys::corearray_read_int8(
                self.inner,
                path.as_ptr(),
                start.len() as i32,
                start.as_ptr(),
                length.as_ptr(),
                out.as_mut_ptr(),
            )
        };
        if rc != 0 {
            return self.err("read_int8");
        }
        Ok(())
    }

    fn read_string_array_1d(
        &self,
        path: &std::ffi::CStr,
        start: i32,
        length: i32,
    ) -> Result<Vec<String>, CohortError> {
        let mut out: Vec<String> = vec![String::new(); length as usize];
        let user_ptr: *mut Vec<String> = &mut out;

        extern "C" fn cb(
            user: *mut c_void,
            idx: i32,
            ptr: *const c_char,
            len: usize,
        ) {
            unsafe {
                let v = &mut *(user as *mut Vec<String>);
                let bytes = std::slice::from_raw_parts(ptr as *const u8, len);
                v[idx as usize] = String::from_utf8_lossy(bytes).into_owned();
            }
        }

        let rc = unsafe {
            corearray_sys::corearray_read_string(
                self.inner,
                path.as_ptr(),
                start,
                length,
                cb,
                user_ptr as *mut c_void,
            )
        };
        if rc != 0 {
            return self.err("read_string");
        }
        Ok(out)
    }
}

impl Drop for GdsHandle {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe { corearray_sys::corearray_close(self.inner) };
            self.inner = std::ptr::null_mut();
        }
    }
}
