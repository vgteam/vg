#include "blocked_gzip_input_stream.hpp"

#include "hfile_cppstream.hpp"

// We need the hFILE* internals available.
#include <hfile_internal.h>

namespace vg {

namespace stream {

using namespace std;

BlockedGzipInputStream::BlockedGzipInputStream(std::istream& stream) : handle(nullptr), byte_count(0),
    know_offset(false) {
    
    // See where the stream is
    stream.clear();
    auto file_start = stream.tellg();
    bool good = stream.good();
    
    // Wrap the stream in an hFILE*
    hFILE* wrapped = hfile_wrap(stream);
    if (wrapped == nullptr) {
        throw runtime_error("Unable to wrap stream");
    }
    
    // Give ownership of it to a BGZF that reads, which we in turn own.
    handle = bgzf_hopen(wrapped, "r");
    if (handle == nullptr) {
        throw runtime_error("Unable to set up BGZF library on wrapped stream");
    }
    
    if (file_start >= 0 && good && bgzf_compression(handle) == 2) {
        // The stream we are wrapping is seekable, and the data is block-compressed
        
        // We need to make sure BGZF knows where its blocks are starting.
    
        // We just freshly opened the BGZF so it thinks it is at 0.
        
        // Tell the BGZF where its next block is actually starting.
        handle->block_address = file_start;
        
        // Remember the virtual offsets will be valid
        know_offset = true;
    }
}

BlockedGzipInputStream::~BlockedGzipInputStream() {
    // Close the GBZF
    bgzf_close(handle);
}

bool BlockedGzipInputStream::Next(const void** data, int* size) {
    if (handle->block_length != 0 && handle->block_offset != handle->block_length) {
        // We aren't just after a seek, but we also aren't at the end of a block. We backed up.
        
        // The already-read data may have started at an offset. But we don't
        // care, because if we back up by X bytes we always re-read the last X
        // bytes of the block.
    
        // Return the unread part of the BGZF file's buffer
        *data = (void*)((char*)handle->uncompressed_block + handle->block_offset);
        *size = handle->block_length - handle->block_offset;
        
        // Send the offset to the end of the block again
        handle->block_offset = handle->block_length;
        handle->uncompressed_address += *size;
        
#ifdef debug
        cerr << "Re-emit " << *size << " bytes of backed-up data, move to " << handle->block_offset << "/" << handle->block_length << endl;
        cerr << "errcode: " << handle->errcode << endl;
#endif
        
        return true;
    } else {
        // We need new data. Either we did a seek, or we are at the end of the previous block.
        
        if (bgzf_compression(handle) != 2) {
            // We're not BGZF compressed. bgzf_read_block only resets the
            // block_offset (when not seeking) for BGZF files. We have to do it
            // manually.
            handle->block_offset = 0;
        }
        
#ifdef debug
        cerr << "Compression mode: " << bgzf_compression(handle) << endl;
        cerr << "Read next block; offset is " << handle->block_offset << endl;
        if (handle->gz_stream != nullptr) {
            cerr << "\tavail_in: " << handle->gz_stream->avail_in << endl;
            if (handle->gz_stream->avail_in > 0) {
                cerr << "\t\tShould not read from backing file!" << endl;
            }
            cerr << "\tavail_out: " << handle->gz_stream->avail_out << endl;
        }
#endif
        
        // Make the BGZF read the next block
        if (bgzf_read_block(handle) != 0) {
            // We have encountered an error
            
#ifdef debug
            cerr << "Failed to read next block" << endl;
            cerr << "\terrcode: " << handle->errcode << endl;
#endif
            
            return false;
        }
        
#ifdef debug
        cerr << "See a block of length " << handle->block_length << " at " << handle->block_address << endl;
        cerr << "\terrcode: " << handle->errcode << endl;
        if (handle->gz_stream != nullptr) {
            cerr << "\tavail_in: " << handle->gz_stream->avail_in << endl;
            cerr << "\tavail_out: " << handle->gz_stream->avail_out << endl;
        }
#endif
        
        if (handle->block_length == 0) {
            // We have hit EOF
            
#ifdef debug
            cerr << "Next block reports length 0 (EOF)" << endl;
            cerr << "\terrcode: " << handle->errcode << endl;
#endif
            
            return false;
        }
        
        // Otherwise we have data.
        
        if (handle->block_offset > handle->block_length) {
            // We don't have enough data to fulfill the most recent seek. Signal an error.
            
#ifdef debug
            cerr << "Tried to seek to " << handle->block_offset << " but got only " <<  handle->block_length << " bytes" << endl;
            cerr << "\terrcode: " << handle->errcode << endl;
#endif
            
            return false;
        }
        
        // Send out the address and size, accounting for seek offset
        *data = (void*)((char*)handle->uncompressed_block + handle->block_offset);
        *size = handle->block_length - handle->block_offset;
        
        // Record the bytes read
        byte_count += handle->block_length - handle->block_offset;
        
        // Tell the BGZF that the cursor is at the end of the block (because it
        // is; subsequent reads come from there)
        handle->block_offset = handle->block_length;
        handle->uncompressed_address += *size;
        
#ifdef debug
        cerr << "Emit " << *size << " bytes in fresh block" << endl;
        cerr << "\terrcode: " << handle->errcode << endl;
#endif
        
        return true;
    }
}

void BlockedGzipInputStream::BackUp(int count) {
    assert(count <= handle->block_offset);
    handle->block_offset -= count;
    handle->uncompressed_address -= count;
    
#ifdef debug
    cerr << "Back up " << count << " bytes to " << handle->block_offset << "/" << handle->block_length << endl;
#endif
}

bool BlockedGzipInputStream::Skip(int count) {
    // We just implement this in terms of next and back up. There's not really
    // a more efficient way, since we can't do relative seeks.
    
    // We have to support this happening immediately after a seek.
    
    while (count > 0) {
        // Keep nexting until we get the block that is count away from where we are.
        const void* ignored_data;
        int size;
        
        if (!Next(&ignored_data, &size)) {
            // We hit EOF, or had an error.
            return false;
        }
        
        // We accomplished this much skipping.
        count -= size;
    }
    
    if (count < 0) {
        // We went too far. But we know we want to go somewhere in this buffer,
        // or we would have finished the loop before we did.
        BackUp(-count);
        count = 0;
    }
    
    return true;
    
}

int64_t BlockedGzipInputStream::ByteCount() const {
    return byte_count;
}

int64_t BlockedGzipInputStream::Tell() const {
    if (know_offset) {
        // Our virtual offsets are true.
       
        // But since we are happy to leave the BGZF's cursor at the ends of
        // blocks, we have to work out what the real virtual offset should be
        // in that case (byte 0 of the next block)
        if (handle->block_offset == handle->block_length) {
            // We need to know where the next block is
            
            // We don't have bgzf_htell so we fake it.
            // We also manually shift the block address to the right place.
            return htell(handle->fp) << 16;
            
        } else {
            // Since we use the BGZF's internal cursor correctly, we can rely on its tell function.
            return bgzf_tell(handle);
        }
    } else {
        // We don't know where the zero position in the stream was, so we can't
        // trust BGZF's virtual offsets.
        return -1;
    }
}

bool BlockedGzipInputStream::Seek(int64_t virtual_offset) {
    if (!know_offset) {
        // We can't seek
        return false;
    }
    
    // Do the seek.
    // This will set handle->block_length to 0, so we know we need to read the block when we read next.
    if(bgzf_seek(handle, virtual_offset, SEEK_SET) == 0) {
        // The seek succeeded
        return true;
    } else {
        // The seek failed
        // Error info is in handle->errcode and handle->fp->has_errno
        // Log the error
        cerr << "error[vg::BlockedGzipInputStream]: BGZF seek error: errcode: "
            << handle->errcode << " errno: " << handle->fp->has_errno << endl;
        return false;
    }
    
    // We won't find out if there's actually data there until we try to read...
}

}

}

