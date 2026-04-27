#include "shim.h"
#include "../vendor/src/CoreArray/CoreArray.h"

#include <new>
#include <string>
#include <vector>

using namespace CoreArray;

struct CaFile {
    CdGDSFile gds;
    std::string err;
};

namespace {

CdAllocArray *resolve_array(CaFile *f, const char *path) {
    CdGDSObj *obj = f->gds.Root().Path(UTF8String(path));
    if (!obj) {
        f->err = std::string("path not found: ") + path;
        return nullptr;
    }
    CdAllocArray *arr = dynamic_cast<CdAllocArray *>(obj);
    if (!arr) {
        f->err = std::string("not an allocatable array: ") + path;
        return nullptr;
    }
    return arr;
}

template <typename Fn>
int trap(CaFile *f, Fn &&body) {
    try {
        body();
        return 0;
    } catch (std::exception &e) {
        f->err = e.what();
    } catch (...) {
        f->err = "unknown C++ exception";
    }
    return 1;
}

}  // namespace

extern "C" {

CaFile *corearray_open(const char *path) {
    static bool registered = false;
    if (!registered) {
        RegisterClass();
        registered = true;
    }
    CaFile *f = new (std::nothrow) CaFile;
    if (!f) return nullptr;
    try {
        f->gds.LoadFile(path, true);
        return f;
    } catch (std::exception &e) {
        f->err = e.what();
        // Hand back the handle so callers can read the error message before close.
        return f;
    } catch (...) {
        f->err = "unknown C++ exception during open";
        return f;
    }
}

void corearray_close(CaFile *f) {
    if (!f) return;
    try { f->gds.CloseFile(); } catch (...) {}
    delete f;
}

const char *corearray_last_error(CaFile *f) {
    return f ? f->err.c_str() : "null file handle";
}

int corearray_ndim(CaFile *f, const char *path, int32_t *out_ndim) {
    return trap(f, [&]() {
        CdAllocArray *arr = resolve_array(f, path);
        if (!arr) throw ErrCoreArray("resolve_array failed");
        *out_ndim = arr->DimCnt();
    });
}

int corearray_dims(CaFile *f, const char *path, int32_t ndim, int64_t *out_dims) {
    return trap(f, [&]() {
        CdAllocArray *arr = resolve_array(f, path);
        if (!arr) throw ErrCoreArray("resolve_array failed");
        if (arr->DimCnt() != ndim) throw ErrCoreArray("ndim mismatch");
        for (int i = 0; i < ndim; ++i) out_dims[i] = arr->GetDLen(i);
    });
}

int corearray_read_int32(
    CaFile *f, const char *path, int32_t ndim,
    const int32_t *start, const int32_t *length, int32_t *out_buf) {
    return trap(f, [&]() {
        CdAllocArray *arr = resolve_array(f, path);
        if (!arr) throw ErrCoreArray("resolve_array failed");
        if (arr->DimCnt() != ndim) throw ErrCoreArray("ndim mismatch");
        arr->ReadData(start, length, out_buf, svInt32);
    });
}

int corearray_read_int8(
    CaFile *f, const char *path, int32_t ndim,
    const int32_t *start, const int32_t *length, int8_t *out_buf) {
    return trap(f, [&]() {
        CdAllocArray *arr = resolve_array(f, path);
        if (!arr) throw ErrCoreArray("resolve_array failed");
        if (arr->DimCnt() != ndim) throw ErrCoreArray("ndim mismatch");
        arr->ReadData(start, length, out_buf, svInt8);
    });
}

int corearray_read_string(
    CaFile *f, const char *path, int32_t start, int32_t length,
    corearray_str_cb cb, void *user) {
    return trap(f, [&]() {
        CdAllocArray *arr = resolve_array(f, path);
        if (!arr) throw ErrCoreArray("resolve_array failed");
        if (arr->DimCnt() != 1) throw ErrCoreArray("string array must be 1-d");
        std::vector<UTF8String> tmp(length);
        int32_t s = start, l = length;
        arr->ReadData(&s, &l, tmp.data(), svStrUTF8);
        for (int32_t i = 0; i < length; ++i) {
            const UTF8String &s = tmp[i];
            cb(user, i, s.data(), s.size());
        }
    });
}

}  // extern "C"
