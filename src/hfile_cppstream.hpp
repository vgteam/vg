/// \file hfile_cppstream.hpp
/// hFILE* C++ streams wrapper
/// Modeled on https://github.com/samtools/htslib-plugins/blob/master/hfile_mmap.c

// We need to provide a C++ stream plugin for htslib hFILE* files so we can
// connect Protobuf Zero Copy Streams to C++ streams while filtering through
// BGZF file handles.

#ifndef VG_HFILE_CPPSTREAM_HPP_INCLUDED
#define VG_HFILE_CPPSTREAM_HPP_INCLUDED

#include <htslib/hfile.h>

#include <istream>
#include <ostream>

namespace vg {

namespace stream {

/// Wrap a C++ output stream as an hFILE* that can be written by BGZF
hFILE* hfile_wrap(std::ostream& output);

/// Wrap a C++ input stream as an hFILE* that can be read by BGZF
hFILE* hfile_wrap(std::istream& input);

}

}

#endif // VG_HFILE_CPPSTREAM_HPP_INCLUDED
