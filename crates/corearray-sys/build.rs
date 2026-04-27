// Compiles the vendored CoreArray GDS C++ library plus a thin C shim
// (`cpp/shim.cpp`) that exposes a flat `extern "C"` surface to Rust.
//
// LZMA support is compiled out via `COREARRAY_NO_LZMA`, which gates lines
// 2164..2784 of dStream.cpp so no liblzma symbols are referenced. SeqArray
// defaults to LZ4_RA / ZIP_RA for genotype storage, so this is the common
// case. LZMA-compressed datasets will fail at read time with a CoreArray
// codec error.

use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vendor = manifest.join("vendor");
    let core = vendor.join("src/CoreArray");
    let geno = vendor.join("src/GenoGDS");
    let lz4 = vendor.join("src/LZ4");
    let zlib = vendor.join("src/ZLIB");
    let include = vendor.join("include");
    let shim = manifest.join("cpp");

    if !vendor.join("src/CoreArray/CoreArray.cpp").exists() {
        panic!(
            "vendor/ submodule not initialized. run: git submodule update --init --recursive"
        );
    }

    let cpp_sources = [
        "CoreArray.cpp",
        "dAllocator.cpp",
        "dAny.cpp",
        "dBase.cpp",
        "dBitGDS.cpp",
        "dEndian.cpp",
        "dFile.cpp",
        "dParallel.cpp",
        "dPlatform.cpp",
        "dRealGDS.cpp",
        "dSerial.cpp",
        "dSparse.cpp",
        "dStrGDS.cpp",
        "dStream.cpp",
        "dStruct.cpp",
        "dVLIntGDS.cpp",
    ];
    let c_core = ["dParallel_Ext.c"];
    let lz4_sources = ["lz4.c", "lz4hc.c", "lz4frame.c", "xxhash.c"];
    let zlib_sources = [
        "adler32.c",
        "compress.c",
        "crc32.c",
        "deflate.c",
        "infback.c",
        "inffast.c",
        "inflate.c",
        "inftrees.c",
        "trees.c",
        "uncompr.c",
        "zutil.c",
    ];

    let mut cxx = cc::Build::new();
    cxx.cpp(true)
        .std("c++11")
        .include(&include)
        .include(&core)
        .include(&geno)
        .define("COREARRAY_NO_LZMA", None)
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-unused-variable")
        .flag_if_supported("-Wno-unused-function")
        .flag_if_supported("-Wno-deprecated-declarations")
        .flag_if_supported("-Wno-misleading-indentation")
        .flag_if_supported("-Wno-class-memaccess")
        .flag_if_supported("-Wno-implicit-fallthrough");
    for f in cpp_sources {
        cxx.file(core.join(f));
    }
    cxx.file(geno.join("GenoGDS.cpp"));
    cxx.file(shim.join("shim.cpp"));
    cxx.compile("corearray_cpp");

    let mut c_lib = cc::Build::new();
    c_lib
        .include(&include)
        .include(&core)
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-implicit-function-declaration");
    for f in c_core {
        c_lib.file(core.join(f));
    }
    for f in lz4_sources {
        c_lib.file(lz4.join(f));
    }
    for f in zlib_sources {
        c_lib.file(zlib.join(f));
    }
    c_lib.compile("corearray_c");

    println!("cargo:rerun-if-changed=cpp/shim.cpp");
    println!("cargo:rerun-if-changed=cpp/shim.h");
    println!("cargo:rerun-if-changed=vendor/src");
    println!("cargo:rerun-if-changed=vendor/include");
}
