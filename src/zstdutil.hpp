//
// -*- coding: utf-8-unix; -*-
//  Copyright (c) 2020 Tencent, Inc.
//     All rights reserved.
//
// Date:   2020/11/30 13:45
// File:   util.h
// Desc:
//

#pragma once

#include <string>
#include <zstd.h>

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

}  // namespace util
