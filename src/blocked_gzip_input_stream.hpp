#ifndef VG_BLOCKED_GZIP_INPUT_STREAM_HPP_INCLUDED
#define VG_BLOCKED_GZIP_INPUT_STREAM_HPP_INCLUDED

#include <google/protobuf/io/zero_copy_stream.h>

#include <htslib/bgzf.h>

namespace vg {

namespace stream {


/// Protobuf-style ZeroCopyInputStream that reads data from blocked gzip
/// format, and allows interacting with virtual offsets.
/// Cannot be moved or copied, because the base class can't be moved or copied.
class BlockedGzipInputStream : public ::google::protobuf::io::ZeroCopyInputStream {

public:
    
    // Does not support construction off a raw BGZF because there's no way to
    // force the file cursor to the start of a new block. And because it's a
    // bad API anyway.
    
    /// Make a new stream reading from the given C++ std::istream, wrapping it
    /// in a BGZF. The stream must be at a BGZF block header, since the header
    /// info is peeked.
    BlockedGzipInputStream(std::istream& stream);

    /// Destroy the stream.
    virtual ~BlockedGzipInputStream();
    
    // Explicitly say we can't be copied/moved, to simplify errors.
    BlockedGzipInputStream(const BlockedGzipInputStream& other) = delete;
    BlockedGzipInputStream& operator=(const BlockedGzipInputStream& other) = delete;
    BlockedGzipInputStream(BlockedGzipInputStream&& other) = delete;
    BlockedGzipInputStream& operator=(BlockedGzipInputStream&& other) = delete;
   
    ///////////////////////////////////////////////////////////////////////////
    // ZeroCopyInputStream interface
    ///////////////////////////////////////////////////////////////////////////
   
    /// Get a buffer to read from. Saves the address of the buffer where data
    /// points, and the size of the buffer where size points. Returns false on
    /// an unrecoverable error or EOF, and true if a buffer was gotten. The
    /// data pointer must be valid until the next read call or until the stream
    /// is destroyed.
    virtual bool Next(const void** data, int* size);
    
    /// When called after Next(), mark the last count bytes of the buffer that
    /// Next() produced as not having been read.
    virtual void BackUp(int count);
    
    /// Skip ahead the given number of bytes. Return false if the end of the
    /// stream is reached, or an error occurs. If the end of the stream is hit,
    /// advances to the end of the stream.
    virtual bool Skip(int count);
    
    /// Get the number of bytes read since the stream was constructed.
    virtual int64_t ByteCount() const;
    
    ///////////////////////////////////////////////////////////////////////////
    // BGZF support interface
    ///////////////////////////////////////////////////////////////////////////
    
    /// Return the blocked gzip virtual offset at which the next fresh buffer
    /// returned by Next() will start, or -1 if operating on an untellable
    /// stream like standard input or on a non-blocked file. Note that this
    /// will only get you the position of the next read if anything you are
    /// reading through is fully backed up to the next actually-unread byte.
    /// See Protobuf's CodedInputStream::Trim().
    virtual int64_t Tell() const;
    
    /// Seek to the given virtual offset. Return true if successful, or false
    /// if the backing stream is unseekable, or not blocked. Note that this
    /// will cause problems if something reading from this stream is still
    /// operating on outstanding buffers; Any CodedInputStreams reading from
    /// this stream *must* be destroyed before this function is called.
    virtual bool Seek(int64_t virtual_offset);
    
protected:
    
    /// The open BGZF handle being read from. We use the BGZF's buffer as our
    /// buffer, and its block_offset for our seeks and back-ups.
    BGZF* handle;
    
    /// The counter to back ByteCount
    size_t byte_count;
    
    /// Flag for whether our backing stream is tellable.
    bool know_offset;
    
};

}

}


#endif
