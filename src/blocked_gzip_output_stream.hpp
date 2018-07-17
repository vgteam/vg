#ifndef VG_BLOCKED_GZIP_OUTPUT_STREAM_HPP_INCLUDED
#define VG_BLOCKED_GZIP_OUTPUT_STREAM_HPP_INCLUDED

#include <google/protobuf/io/zero_copy_stream.h>

#include <htslib/bgzf.h>

namespace vg {

namespace stream {


/// Protobuf-style ZeroCopyOutputStream that writes data in blocked gzip
/// format, and allows interacting with virtual offsets.
class BlockedGzipOutputStream : public ::google::protobuf::io::ZeroCopyOutputStream {

public:
    /// Make a new stream outputting to the given open BGZF file handle.
    /// The stream will own the BGZF file and close it when destructed.
    BlockedGzipOutputStream(BGZF* bgzf_handle);
    
    /// Make a new stream outputting to the given C++ std::ostream, wrapping it
    /// in a BGZF.
    BlockedGzipOutputStream(std::ostream& stream);

    /// Destroy the stream, finishing all writes if necessary.
    virtual ~BlockedGzipOutputStream();
    
    /// Get a buffer to write to. Saves the address of the buffer where data
    /// points, and the size of the buffer where size points. Returns false on
    /// an unrecoverable error, and true if a buffer was gotten. The stream is
    /// responsible for making sure data in the buffer makes it into the
    /// output. The data pointer must be valid until the next write call or
    /// until the stream is destroyed.
    virtual bool Next(void** data, int* size);
    
    /// When called after Next(), mark the last count bytes of the buffer that
    /// Next() produced as not to be written to the output. The user must not
    /// have touched those bytes.
    virtual void BackUp(int count);
    
    /// Get the number of bytes written since the stream was constructed.
    virtual int64_t ByteCount() const;
    
    /// Take the given data at the given address into the stream as written.
    /// Only works if AllowsAliasing() returns true. Returns true on success,
    /// and false on an unrecoverable error.
    virtual bool WriteAliasedRaw(const void * data, int size);
    
    /// Return true if WriteAliasedRaw() is actually available, and false otherwise.
    virtual bool AllowsAliasing() const;
    
private:
    
    /// Actually dump the buffer data to the file, if needed. Sadly, we can't
    /// really be zero-copy because the bgzf library isn't.
    /// Throws on failure.
    void flush();
    
    /// The open BGZF handle being written to
    BGZF* handle;
    
    /// This vector will own the memory we use as our void* buffer.
    std::vector<char> buffer;
    
    /// The number of characters that have been backed up from the end of the buffer
    size_t backed_up;
    
    /// The counter to back ByteCount
    size_t byte_count;
    
};

/// Wrap a C++ input stream as an hFILE* that can be read by BGZF
hFILE* hfile_wrap(std::istream& input);
/// Wrap a C++ output stream as an hFILE* that can be written by BGZF
hFILE* hfile_wrap(std::ostream& output);

}

}


#endif
