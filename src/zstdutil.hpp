
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <zstd.h>

/** \file
 * Wrappers for Zstandard compression and decompression.
 *
 * TODO: Override xsputn, xsgetn for faster compression?
 * TODO: Move constructors for streams?
 * TODO: is_open(), close() for streams?
 */

namespace vg {

//------------------------------------------------------------------------------

/// Zstandard compression buffer that writes to another stream buffer.
class zstd_compress_buf : public std::streambuf {
public:
    explicit zstd_compress_buf(std::streambuf* inner, int compression_level = ZSTD_defaultCLevel());
    ~zstd_compress_buf();

    zstd_compress_buf(const zstd_compress_buf&) = delete;
    zstd_compress_buf& operator=(const zstd_compress_buf&) = delete;
    zstd_compress_buf(zstd_compress_buf&&) = default;
    zstd_compress_buf& operator=(zstd_compress_buf&&) = default;

protected:
    int_type overflow(int_type ch) override;
    int sync() override;

    std::streambuf* inner;
    ZSTD_CCtx* context;
    std::vector<char> in_buffer;
    std::vector<char> out_buffer;
};

/// Zstandard decompression buffer that reads from another stream buffer.
class zstd_decompress_buf : public std::streambuf {
public:
    explicit zstd_decompress_buf(std::streambuf* inner);
    ~zstd_decompress_buf();

    zstd_decompress_buf(const zstd_decompress_buf&) = delete;
    zstd_decompress_buf& operator=(const zstd_decompress_buf&) = delete;
    zstd_decompress_buf(zstd_decompress_buf&&) = default;
    zstd_decompress_buf& operator=(zstd_decompress_buf&&) = default;

protected:
    int_type underflow() override;

    std::streambuf* inner;
    ZSTD_DCtx* context;
    std::vector<char> in_buffer;
    std::vector<char> out_buffer;
    size_t in_offset;
};

//------------------------------------------------------------------------------

/// Zstandard output file stream.
/// The object cannot be copied or moved.
class zstd_ofstream : public std::ostream {
public:
    explicit zstd_ofstream(const std::string& filename, int compression_level = ZSTD_defaultCLevel()) :
        std::ostream(&buffer),
        inner(filename, std::ios::binary),
        buffer(inner.rdbuf(), compression_level) {}

    zstd_ofstream(const zstd_ofstream&) = delete;
    zstd_ofstream& operator=(const zstd_ofstream&) = delete;
    zstd_ofstream(zstd_ofstream&&) = delete;
    zstd_ofstream& operator=(zstd_ofstream&&) = delete;

protected:
    std::ofstream inner;
    zstd_compress_buf buffer;
};

/// Zstandard input file stream.
/// The object cannot be copied or moved.
class zstd_ifstream : public std::istream {
public:
    explicit zstd_ifstream(const std::string& filename) :
        std::istream(&buffer),
        inner(filename, std::ios::binary),
        buffer(inner.rdbuf()) {}

    zstd_ifstream(const zstd_ifstream&) = delete;
    zstd_ifstream& operator=(const zstd_ifstream&) = delete;
    zstd_ifstream(zstd_ifstream&&) = delete;
    zstd_ifstream& operator=(zstd_ifstream&&) = delete;

protected:
    std::ifstream inner;
    zstd_decompress_buf buffer;
};

//------------------------------------------------------------------------------

} // namespace vg

// TODO: Get rid of these when we have something better.

//
// -*- coding: utf-8-unix; -*-
//  Copyright (c) 2020 Tencent, Inc.
//     All rights reserved.
//
// Date:   2020/11/30 13:45
// File:   util.h
// Desc:
//

namespace zstdutil {

const int DEFAULTCOMPRESSLEVEL = 5;

// if return code not 0 is error
int CompressString(const std::string& src, std::string& dst,
                   int compressionlevel = DEFAULTCOMPRESSLEVEL);

// if return code not 0 is error
int DecompressString(const std::string& src, std::string& dst);

// if return code not 0 is error
int StreamDecompressString(const std::string& src, std::string& dst,
                           int compressionlevel = DEFAULTCOMPRESSLEVEL);

// if return code not 0 is error
int StreamCompressString(const std::string& src, std::string& dst,
                         int compressionlevel = DEFAULTCOMPRESSLEVEL);

}  // namespace zstdutil

//------------------------------------------------------------------------------
