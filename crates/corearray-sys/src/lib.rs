//! Raw `extern "C"` bindings to the CoreArray GDS library via the C++ shim
//! in `cpp/shim.cpp`. All functions return 0 on success; non-zero return
//! means the file's last_error string is populated.
//!
//! Higher-level Rust code (`favor::ingest::gds`) builds the streaming
//! `VariantReader` on top of this surface.

#![allow(non_camel_case_types)]

use std::ffi::c_void;
use std::os::raw::{c_char, c_int};

#[repr(C)]
pub struct CaFile {
    _private: [u8; 0],
}

pub type corearray_str_cb = extern "C" fn(
    user: *mut c_void,
    idx: i32,
    ptr: *const c_char,
    len: usize,
);

extern "C" {
    pub fn corearray_open(path: *const c_char) -> *mut CaFile;
    pub fn corearray_close(file: *mut CaFile);

    pub fn corearray_last_error(file: *mut CaFile) -> *const c_char;

    pub fn corearray_ndim(file: *mut CaFile, path: *const c_char, out_ndim: *mut i32) -> c_int;
    pub fn corearray_dims(
        file: *mut CaFile,
        path: *const c_char,
        ndim: i32,
        out_dims: *mut i64,
    ) -> c_int;

    pub fn corearray_read_int32(
        file: *mut CaFile,
        path: *const c_char,
        ndim: i32,
        start: *const i32,
        length: *const i32,
        out_buf: *mut i32,
    ) -> c_int;

    pub fn corearray_read_int8(
        file: *mut CaFile,
        path: *const c_char,
        ndim: i32,
        start: *const i32,
        length: *const i32,
        out_buf: *mut i8,
    ) -> c_int;

    pub fn corearray_read_string(
        file: *mut CaFile,
        path: *const c_char,
        start: i32,
        length: i32,
        cb: corearray_str_cb,
        user: *mut c_void,
    ) -> c_int;
}
