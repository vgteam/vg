#include "zstdutil.hpp"

namespace vg {

//------------------------------------------------------------------------------

zstd_compress_buf::zstd_compress_buf(std::streambuf* inner, int compression_level) :
    inner(inner), context(ZSTD_createCCtx())
{
    this->in_buffer.resize(ZSTD_CStreamInSize());
    this->setp(this->in_buffer.data(), this->in_buffer.data() + this->in_buffer.size());
    this->out_buffer.resize(ZSTD_CStreamOutSize());
    ZSTD_CCtx_setParameter(this->context, ZSTD_c_compressionLevel, compression_level);
}

zstd_compress_buf::~zstd_compress_buf() {
    this->sync();
    ZSTD_freeCCtx(this->context); this->context = nullptr;
}

zstd_compress_buf::int_type zstd_compress_buf::overflow(int_type ch) {
    if (ch != traits_type::eof()) {
        if (this->sync() == -1) {
            return traits_type::eof();
        }
        *this->pptr() = traits_type::to_char_type(ch);
        this->pbump(1);
    }
    return ch;
}

zstd_compress_buf::int_type zstd_compress_buf::sync() {
    if (this->inner == nullptr) {
        throw std::runtime_error("zstd_compress_buf: inner stream buffer is null");
    }

    ZSTD_inBuffer input = { this->pbase(), static_cast<size_t>(this->pptr() - this->pbase()), 0 };
    ZSTD_EndDirective mode = (this->pptr() < this->epptr() ? ZSTD_e_end : ZSTD_e_continue);
    bool finished = false;
    while (!finished) {
        ZSTD_outBuffer output = { this->out_buffer.data(), this->out_buffer.size(), 0 };
        size_t result = ZSTD_compressStream2(this->context, &output, &input, mode);
        if (ZSTD_isError(result)) {
            std::string msg = "zstd_compress_buf: compression failed: " + std::string(ZSTD_getErrorName(result));
            throw std::runtime_error(msg);
        }
        size_t n = this->inner->sputn(this->out_buffer.data(), output.pos);
        if (n != output.pos) {
            throw std::runtime_error("zstd_compress_buf: failed to write compressed data");
        }
        finished = (input.pos >= input.size);
        if (mode == ZSTD_e_end) {
            finished &= (result == 0);
        }
    }

    this->setp(this->in_buffer.data(), this->in_buffer.data() + this->in_buffer.size());
    return 0;
}

//------------------------------------------------------------------------------

zstd_decompress_buf::zstd_decompress_buf(std::streambuf* inner) :
    inner(inner), context(ZSTD_createDCtx())
{
    this->in_buffer.resize(ZSTD_DStreamInSize());
    this->in_offset = this->in_buffer.size();
    this->out_buffer.resize(ZSTD_DStreamOutSize());
}

zstd_decompress_buf::~zstd_decompress_buf() {
    ZSTD_freeDCtx(this->context); this->context = nullptr;
}

zstd_decompress_buf::int_type zstd_decompress_buf::underflow() {
    if (this->gptr() < this->egptr()) {
        return traits_type::to_int_type(*this->gptr());
    }

    // Fill the input buffer if necessary.
    if (this->in_offset >= this->in_buffer.size()) {
        size_t n = this->inner->sgetn(this->in_buffer.data(), this->in_buffer.size());
        this->in_offset = 0;
        this->in_buffer.resize(n);
    }

    // Decompress the data into the input buffer.
    ZSTD_inBuffer input = { this->in_buffer.data(), this->in_buffer.size(), this->in_offset };
    ZSTD_outBuffer output = { this->out_buffer.data(), this->out_buffer.size(), 0 };
    size_t result = ZSTD_decompressStream(this->context, &output, &input);
    if (ZSTD_isError(result)) {
        std::string msg = "zstd_decompress_buf: decompression failed: " + std::string(ZSTD_getErrorName(result));
        throw std::runtime_error(msg);
    }
    this->in_offset = input.pos;

    // Tell the stream to use the output buffer.
    this->setg(this->out_buffer.data(), this->out_buffer.data(), this->out_buffer.data() + output.pos);
    return (output.pos > 0 ? traits_type::to_int_type(*this->gptr()) : traits_type::eof());
}

//------------------------------------------------------------------------------

} // namespace vg

//------------------------------------------------------------------------------

//
// -*- coding: utf-8-unix; -*-
//  Copyright (c) 2020 Tencent, Inc.
//     All rights reserved.
//
// Date:   2020/11/30 13:45
// File:   zstd.cc
// Desc:
//

namespace zstdutil {

int CompressString(const std::string& src, std::string& dst, int compressionlevel) {
  size_t const cBuffSize = ZSTD_compressBound(src.size());
  dst.resize(cBuffSize);
  auto dstp = const_cast<void*>(static_cast<const void*>(dst.c_str()));
  auto srcp = static_cast<const void*>(src.c_str());
  size_t const cSize = ZSTD_compress(dstp, cBuffSize, srcp, src.size(), compressionlevel);
  auto code = ZSTD_isError(cSize);
  if (code) {
    return code;
  }
  dst.resize(cSize);
  return code;
}

int DecompressString(const std::string& src, std::string& dst) {
  size_t const cBuffSize = ZSTD_getFrameContentSize(src.c_str(), src.size());

  if (0 == cBuffSize) {
    return cBuffSize;
  }

  if (ZSTD_CONTENTSIZE_UNKNOWN == cBuffSize) {
    return StreamDecompressString(src, dst);
  }

  if (ZSTD_CONTENTSIZE_ERROR == cBuffSize) {
    return -2;
  }

  dst.resize(cBuffSize);
  auto dstp = const_cast<void*>(static_cast<const void*>(dst.c_str()));
  auto srcp = static_cast<const void*>(src.c_str());
  size_t const cSize = ZSTD_decompress(dstp, cBuffSize, srcp, src.size());
  auto code = ZSTD_isError(cSize);
  if (code) {
    return code;
  }
  dst.resize(cSize);
  return code;
}

int StreamCompressString(const std::string& src, std::string& dst, int compressionlevel) {
  size_t const buffInSize = ZSTD_CStreamInSize();
  std::string buffInTmp;
  buffInTmp.reserve(buffInSize);
  auto buffIn = const_cast<void*>(static_cast<const void*>(buffInTmp.c_str()));

  auto buffOutSize = ZSTD_CStreamOutSize();
  std::string buffOutTmp;
  buffOutTmp.reserve(buffOutSize);
  auto buffOut = const_cast<void*>(static_cast<const void*>(buffOutTmp.c_str()));

  ZSTD_CCtx* const cctx = ZSTD_createCCtx();
  ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionlevel);

  size_t const toRead = buffInSize;
  auto local_pos = 0;
  auto buff_tmp = const_cast<char*>(buffInTmp.c_str());
  for (;;) {
    size_t read = src.copy(buff_tmp, toRead, local_pos);
    local_pos += read;

    int const lastChunk = (read < toRead);
    ZSTD_EndDirective const mode = lastChunk ? ZSTD_e_end : ZSTD_e_continue;

    ZSTD_inBuffer input = {buffIn, read, 0};
    int finished;

    do {
      ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
      size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
      dst.insert(dst.end(), buffOutTmp.begin(), buffOutTmp.begin() + output.pos);
      finished = lastChunk ? (remaining == 0) : (input.pos == input.size);
    } while (!finished);

    if (lastChunk) {
      break;
    }
  }
  
  ZSTD_freeCCtx(cctx);

  return 0;
}

int StreamDecompressString(const std::string& src, std::string& dst, int compressionlevel) {
  size_t const buffInSize = ZSTD_DStreamInSize();
  std::string buffInTmp;
  buffInTmp.reserve(buffInSize);
  auto buffIn = const_cast<void*>(static_cast<const void*>(buffInTmp.c_str()));

  auto buffOutSize = ZSTD_DStreamOutSize();
  std::string buffOutTmp;
  buffOutTmp.reserve(buffOutSize);
  auto buffOut = const_cast<void*>(static_cast<const void*>(buffOutTmp.c_str()));

  ZSTD_DCtx* const dctx = ZSTD_createDCtx();

  size_t const toRead = buffInSize;
  size_t read;
  size_t last_ret = 0;
  size_t local_pos = 0;
  auto buff_tmp = const_cast<char*>(buffInTmp.c_str());

  while ((read = src.copy(buff_tmp, toRead, local_pos))) {
    local_pos += read;
    ZSTD_inBuffer input = {buffIn, read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
      size_t const ret = ZSTD_decompressStream(dctx, &output, &input);
      dst.insert(dst.end(), buffOutTmp.begin(), buffOutTmp.begin() + output.pos);
      last_ret = ret;
    }
  }
  
  ZSTD_freeDCtx(dctx);

  if(last_ret != 0) {
    return -3;
  }

  return 0;
}

}  // namespace zstdutil

//------------------------------------------------------------------------------
