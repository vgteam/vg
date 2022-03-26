//
// -*- coding: utf-8-unix; -*-
//  Copyright (c) 2020 Tencent, Inc.
//     All rights reserved.
//
// Date:   2020/11/30 13:45
// File:   zstd.cc
// Desc:
//

#include "zstdutil.hpp"

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

}  // namespace util
