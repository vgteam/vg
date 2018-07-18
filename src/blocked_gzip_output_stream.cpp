#include "blocked_gzip_output_stream.hpp"

// We need to provide a C++ stream plugin for htslib hFILE* files so we can
// connect Protobuf Zero Copy Streams to C++ streams while filtering through
// BGZF file handles.
#include <hfile_internal.h>

#include <errno.h>

#define debug

namespace vg {

namespace stream {

using namespace std;

BlockedGzipOutputStream::BlockedGzipOutputStream(BGZF* bgzf_handle) : handle(bgzf_handle), buffer(), backed_up(0), byte_count(0) {
    // Nothing to do
}

BlockedGzipOutputStream::BlockedGzipOutputStream(std::ostream& stream) : handle(nullptr), buffer(), backed_up(0), byte_count(0) {
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
}

BlockedGzipOutputStream::~BlockedGzipOutputStream() {
    // Make sure to finish writing before destructing.
    flush();
    // Close the GBZF
    bgzf_close(handle);
}

bool BlockedGzipOutputStream::Next(void** data, int* size) {
    try {
        // Dump data if we have it
        flush();
        
        // Allocate some space in the buffer
        buffer.resize(4096);
        
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
}

int64_t BlockedGzipOutputStream::ByteCount() const {
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
    // Make sure all data has been sent to BGZF
    flush();
    
    // See where we are now
    return bgzf_tell(handle);
}

void BlockedGzipOutputStream::flush() {
    // How many bytes are left to write?
    auto outstanding = buffer.size() - backed_up;
    if (outstanding > 0) {
        // Save the buffer
        auto written = bgzf_write(handle, (void*)&buffer[0], outstanding);
        
        if (written != outstanding) {
            // This only happens when there is an error
            throw runtime_error("IO error writing data in BlockedGzipOutputStream");
        }
        
        // Record the actual write
        byte_count += written;
    }
}

// hFILE* C++ streams plugin
// Modeled on https://github.com/samtools/htslib-plugins/blob/master/hfile_mmap.c

/// Define a c-style-inheritance derived struct that holds the hFILE and the
/// stream pointers. Either stream pointer may be null (if we are in the other
/// mode), or both can be non-null and point to the same iostream object.
typedef struct {
    hFILE base;
    istream* input;
    ostream* output;
} hFILE_cppstream;


// Define read, write, seek (which also can tell), flush, and close functions

/// Read data. Return bytes read, or a negative value on error. Set errno on error.
static ssize_t cppstream_read(hFILE *fpv, void *buffer, size_t nbytes) {
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;
    
    if (fp->input == nullptr) {
        // No input stream
        errno = EBADF;
        return -1;
    }
    
    // Read the data and record how much we got
    fp->input->clear();
    fp->input->read((char*) buffer, nbytes);
    ssize_t found = fp->input->gcount();
    
    // TODO: it is not clear how we should signal EOF
    
    if (!fp->input->good() && !fp->input->eof()) {
        // An error happened
        errno = EIO;
        return -1;
    }
    
    // Otherwise the read worked
    return found;
}

/// Write data. Return the number of bytes actually written. Return a negative
/// value and set errno on error.
static ssize_t cppstream_write(hFILE *fpv, const void *buffer, size_t nbytes) {
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;
    
    if (fp->output == nullptr) {
        // No output stream
        errno = EBADF;
        return -1;
    }
    
    // Write the data and record how much we put
    fp->output->clear();
    fp->output->write((char*) buffer, nbytes);
    
    // TODO: it is not clear how we should signal EOF
    
    if (!fp->output->good()) {
        // An error happened
        errno = EIO;
        return -1;
    }
    
    // Otherwise the write worked, and we wrote all the bytes
    return nbytes;
}

/// Seek relative to SEEK_SET (beginning), SEEK_CUR, or SEEK_END. Return the
/// resulting offset from the beginning of the file.
/// Returns a negative value on error.
static off_t cppstream_seek(hFILE *fpv, off_t offset, int whence) {

#ifdef debug
    cerr << "cppstream_seek(" << fpv << ", " << offset << ", " << whence << ")" << endl;
#endif
    
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;
    
    // How are we seeking?
    ios_base::seekdir way;
    switch (whence) {
    case SEEK_SET:
        way = ios_base::beg;
        break;
    case SEEK_CUR:
        way = ios_base::cur;
        break;
    case SEEK_END:
        way = ios_base::end;
        break;
    default:
        errno = EINVAL;
        return -1;
    }
    
    off_t arrived_at = 0;
    
    if (fp->input != nullptr) {
        // Seek the input stream
        fp->input->clear();
        fp->input->seekg(offset, way);
        if (!fp->input->good()) {
            // Seek failed.
            // Assume it is because this is a pipe.
            errno = ESPIPE;
            return -1;
        }
        
        auto reached = fp->input->tellg();
        if (reached == -1) {
            // Definitely a pipe
            errno = ESPIPE;
            return -1;
        }
        
        arrived_at = reached;
    }
    
    if (fp->output != nullptr) {
        // Seek the output stream
        fp->output->clear();
        fp->output->seekp(offset, way);
        if (!fp->output->good()) {
            // Seek failed.
            // Assume it is because this is a pipe.
            errno = ESPIPE;
            return -1;
        }
        
        auto reached = fp->output->tellp();
        if (reached == -1) {
            // Definitely a pipe
            errno = ESPIPE;
            return -1;
        }
        
        if (fp->input != nullptr && reached != arrived_at) {
            // We have two streams and they are out of sync!
            errno = EIO;
            return -1;
        }
        
        arrived_at = reached;
    }
    
#ifdef debug
    cerr << "\t" << arrived_at << endl;
#endif
    
    // We worked!
    return arrived_at;
}

/// Flush the output stream, if we are doing output. Return 0 for success, or a
/// negative number and set errno on error.
static int cppstream_flush(hFILE *fpv) {
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;

    if (fp->output != nullptr) {
        // We have an output stream to flush
        fp->output->clear();
        fp->output->flush();
        
        if (!fp->output->good()) {
            // Flushing did not work
            errno = EIO;
            return -1;
        }
    }
    
    return 0;
}

/// Close the file. Return 0 on success, or a negative number and set errno on
/// failure.
static int cppstream_close(hFILE *fpv) {
    // This is tricky because we don't own the stream. We also can't close generic istreams and ostreams.
    
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;
    
    // Just null out the stream fields. They will be closed when destroyed, and we don't own them. 
    fp->input = nullptr;
    fp->output = nullptr;
    
    return 0;
    
}

/// Define an hFILE backend for cpp streams
static const struct hFILE_backend cppstream_backend = {
    cppstream_read,
    cppstream_write,
    cppstream_seek,
    cppstream_flush,
    cppstream_close
};

hFILE* hfile_wrap(std::istream& input) {
    /// Make the base struct, making sure it knows how big we are
    hFILE_cppstream* fp = (hFILE_cppstream*) hfile_init(sizeof(hFILE_cppstream), "r", 0);
    
    if (fp == nullptr) {
        // Couldn't allocate the file for some reason?
        return nullptr;
    }

    // Do our initialization
    fp->input = &input;
    fp->output = nullptr;
    
    // Set the backend
    fp->base.backend = &cppstream_backend;
    
    // Tell the file that it is starting at the offset that the stream is at
    input.clear();
    auto start_pos = input.tellg();
    if (start_pos < 0 || !input.good()) {
        throw runtime_error("Could not determine initial input position");
    }
    fp->base.offset = start_pos;
    
    // Return the base hFILE*
    return &fp->base;
}

hFILE* hfile_wrap(std::ostream& output) {
    /// Make the base struct, making sure it knows how big we are
    hFILE_cppstream* fp = (hFILE_cppstream*) hfile_init(sizeof(hFILE_cppstream), "w", 0);
    
    if (fp == nullptr) {
        // Couldn't allocate the file for some reason?
        return nullptr;
    }

    // Do our initialization
    fp->input = nullptr;
    fp->output = &output;
    
    // Set the backend
    fp->base.backend = &cppstream_backend;
    
    // Tell the file that it is starting at the offset that the stream is at
    output.clear();
    auto start_pos = output.tellp();
    if (start_pos < 0 || !output.good()) {
        throw runtime_error("Could not determine initial output position");
    }
    fp->base.offset = start_pos;
    
    // Return the base hFILE*
    return &fp->base;
}


}

}

