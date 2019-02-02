#include "blocked_gzip_output_stream.hpp"

#include "hfile_cppstream.hpp"

// We need the hFILE* internals available.
#include <hfile_internal.h>

namespace vg {

namespace stream {

using namespace std;

BlockedGzipOutputStream::BlockedGzipOutputStream(BGZF* bgzf_handle) : handle(bgzf_handle), buffer(), backed_up(0), byte_count(0),
    know_offset(false), end_file(false) {
    
    if (handle->mt) {
        // I don't want to deal with BGZF multithreading, because I'm going to be hacking its internals
        throw runtime_error("Multithreaded BGZF is not supported");
    }
    
    // Force the BGZF to start a new block by flushing the old one, if it exists.
    if (bgzf_flush(handle) != 0) {
        throw runtime_error("Unable to flush BGZF");
    }
   
    // Try seeking the hfile's backend to exactly the position it is at, to get the actual offset.
    // This lets us know if the stream is really seekable/tellable, because htell always works.
    auto cur_pos = (*(handle->fp->backend->seek))(handle->fp, 0, SEEK_CUR);
    if (cur_pos >= 0) {
        // The seek succeeded. We know where we are, and so, we assume, does
        // the hFILE.
        
        // Tell the BGZF where it is (which is at the hFILE's position rather
        // than the backend's, but we know the hFILE position is correct)
        handle->block_address = htell(handle->fp);
        
        // We are backed by a tellable stream
        know_offset = true;
    }
}

BlockedGzipOutputStream::BlockedGzipOutputStream(std::ostream& stream) : handle(nullptr), buffer(), backed_up(0), byte_count(0),
    know_offset(false), end_file(false) {
    
    // Wrap the stream in an hFILE*
    hFILE* wrapped = hfile_wrap(stream);
    if (wrapped == nullptr) {
        throw runtime_error("Unable to wrap stream");
    }
    
    // Give ownership of it to a BGZF that writes, which we in turn own.
    handle = bgzf_hopen(wrapped, "w");
    if (handle == nullptr) {
        throw runtime_error("Unable to set up BGZF library on wrapped stream");
    }
    
    stream.clear();
    auto file_start = stream.tellp();
    if (file_start >= 0 && stream.good()) {
        // The stream we are wrapping is seekable.
        
        // We need to make sure BGZF knows where its blocks are starting.
    
        // No need to flush because we just freshly opened the BGZF
        
        // Tell the BGZF where its next block is actually starting.
        handle->block_address = file_start;
        
        // Remember the virtual offsets will be valid
        know_offset = true;
    }
}

BlockedGzipOutputStream::~BlockedGzipOutputStream() {

#ifdef debug
    cerr << "Destroying BlockedGzipOutputStream" << endl;
#endif

    // Make sure to finish writing before destructing.
    flush();
    
    if (end_file) {
        // Close the file with an EOF block.
#ifdef debug
        cerr << "Close BlockedGzipOutputStream normally with EOF" << endl;
#endif
        bgzf_close(handle);
    } else {
        // Close the BGZF *without* writing an EOF block.
#ifdef debug
        cerr << "Force BlockedGzipOutputStream closed" << endl;
#endif
        force_close();
    }
    
#ifdef debug
    cerr << "BlockedGzipOutputStream destroyed" << endl;
#endif

}

bool BlockedGzipOutputStream::Next(void** data, int* size) {
    try {
        // Dump data if we have it
        flush();
        
        // Allocate some space in the buffer
        buffer.resize(4096);
        
#ifdef debug
        cerr << "Allocate buffer of " << buffer.size() << " bytes " << endl;
#endif
        
        // None of it is backed up
        backed_up = 0;
        
        // Tell the caller where to write
        *data = (void*)&buffer[0];
        *size = buffer.size();
        
        // It worked
        return true;
        
    } catch(exception e) {
        return false;
    }
}

void BlockedGzipOutputStream::BackUp(int count) {
    backed_up += count;
    assert(backed_up <= buffer.size());
    
#ifdef debug
    cerr << "Back up " << count << " bytes to " << (buffer.size() - backed_up) << " still written" << endl;
#endif
}

int64_t BlockedGzipOutputStream::ByteCount() const {
#ifdef debug
    cerr << "Report total bytes written as " << byte_count << endl;
#endif
    return byte_count;
}

bool BlockedGzipOutputStream::WriteAliasedRaw(const void* data, int size) {
    // Not allowed
    return false;
}

bool BlockedGzipOutputStream::AllowsAliasing() const {
    return false;
}


int64_t BlockedGzipOutputStream::Tell() {
    if (know_offset) {
        // Our virtual offsets are true.
        
        // Make sure all data has been sent to BGZF
        flush();
        
        // See where we are now. No de-aliasing is necessary; the BGZF never
        // leaves the cursor past the end of the block when writing, so we
        // always have the cannonical virtual offset.
        return bgzf_tell(handle);
    } else {
        // We don't know where the zero position in the stream was, so we can't
        // trust BGZF's virtual offsets.
        return -1;
    }
}

void BlockedGzipOutputStream::StartFile() {
    // We know since nothing has been written that we are working with a fresh
    // BGZF at what it thinks is virtual offset 0.
    assert(bgzf_tell(handle) == 0);
    know_offset = true;
}

void BlockedGzipOutputStream::EndFile() {
#ifdef debug
    cerr << "Setting BlockedGzipOutputStream end_file flag" << endl;
#endif
    end_file = true;
}

void BlockedGzipOutputStream::flush() {
    // How many bytes are left to write?
    auto outstanding = buffer.size() - backed_up;
    if (outstanding > 0) {
#ifdef debug
        cerr << "Flush " << outstanding << " bytes to BGZF" << endl;
#endif
    
        // Save the buffer
        auto written = bgzf_write(handle, (void*)&buffer[0], outstanding);
        
        if (written != outstanding) {
            // This only happens when there is an error
            throw runtime_error("IO error writing data in BlockedGzipOutputStream");
        }
        
        // Record the actual write
        byte_count += written;
        
        // Make sure we don't try and write the same data twice by scrapping the buffer.
        buffer.resize(0);
        backed_up = 0;
    }
}

void BlockedGzipOutputStream::force_close() {
    // Sneakily close the BGZF file without letting it write an EOF empty block marker.
    
#ifdef debug
    cerr << "Forcing BlockedGzipOutputStream closed without EOF" << endl;
#endif
    
    // Flush the data, which the close function won't do in the path we want it to take
    if (bgzf_flush(handle) != 0) {
        throw runtime_error("Could not flush the BGZF");
    }
    
    // Lie to BGZF and tell it that it did not just write compressed data.
    // This causes close to bypass the EOF block write.
    handle->is_compressed = 0;
    
    // Do the close operation, which does all the other cleanup still.
    bgzf_close(handle);
    handle = nullptr;
}

}

}

