// Flat C surface over CoreArray's C++ API. Each entry point traps any C++
// exception and converts it to a non-zero return code; the message is
// retrievable via corearray_last_error.
//
// All functions return 0 on success.

#ifndef COREARRAY_SHIM_H
#define COREARRAY_SHIM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CaFile CaFile;

CaFile *corearray_open(const char *path);
void corearray_close(CaFile *file);

const char *corearray_last_error(CaFile *file);

// Returns the rank (number of dimensions) of the array at `path`.
int corearray_ndim(CaFile *file, const char *path, int32_t *out_ndim);

// Fills out_dims[0..ndim] with the per-axis sizes.
int corearray_dims(CaFile *file, const char *path, int32_t ndim, int64_t *out_dims);

// Reads a contiguous slab from an N-dim array, casting to int32. Caller passes
// `start` and `length` arrays of `ndim` entries; total elements = product(length).
int corearray_read_int32(
    CaFile *file,
    const char *path,
    int32_t ndim,
    const int32_t *start,
    const int32_t *length,
    int32_t *out_buf);

// Same shape as read_int32 but casts to int8. Used for bit2-packed dosages
// in /genotype/data, where per-allele codes 0/1/2/3 (3=missing) come out as
// signed 8-bit integers; pair adjacent alleles for diploid dosage.
int corearray_read_int8(
    CaFile *file,
    const char *path,
    int32_t ndim,
    const int32_t *start,
    const int32_t *length,
    int8_t *out_buf);

// Reads a 1-d UTF-8 string array slice. For each string in [start, start+length),
// invokes `cb(user, idx, ptr, len)`. Strings are NOT null-terminated; copy if you
// need to keep them past the callback.
typedef void (*corearray_str_cb)(void *user, int32_t idx, const char *ptr, size_t len);
int corearray_read_string(
    CaFile *file,
    const char *path,
    int32_t start,
    int32_t length,
    corearray_str_cb cb,
    void *user);

#ifdef __cplusplus
}
#endif

#endif
