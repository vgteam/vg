#include "hfile_cppstream.hpp"

// We need the hFILE* internals available.
#include <hfile_internal.h>

#include <errno.h>

#include <iostream>

namespace vg {

namespace stream {

using namespace std;

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
#ifdef debug
    cerr << "cppstream_read(" << fpv << ", " << buffer << ", " << nbytes << ")" << endl;
#endif
    
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
    
    if (!fp->input->good() && !fp->input->eof()) {
        // An error happened
        errno = EIO;
        return -1;
    }
    
#ifdef debug
    cerr << "\tFound " << found << "/" << nbytes << " bytes up to " << fp->input->tellg() << endl;
#endif
    
    // Otherwise the read worked
    return found;
}

/// Write data. Return the number of bytes actually written. Return a negative
/// value and set errno on error.
static ssize_t cppstream_write(hFILE *fpv, const void *buffer, size_t nbytes) {
#ifdef debug
    cerr << "cppstream_write(" << fpv << ", " << buffer << ", " << nbytes << ")" << endl;
#endif
    
    // Cast the hFILE to the derived class
    hFILE_cppstream* fp = (hFILE_cppstream*) fpv;
    
    if (fp->output == nullptr) {
        // No output stream
        errno = EBADF;
        return -1;
    }
    
    // Write the data and record how much we put
    fp->output->clear();
    // Note that the stream always takes all the bytes
    fp->output->write((char*) buffer, nbytes);
    
    if (!fp->output->good()) {
        // An error happened
        errno = EIO;
        return -1;
    }
    
#ifdef debug
    cerr << "Successfully wrote " << nbytes << " bytes" << endl;
#endif
    
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
#ifdef debug
    cerr << "cppstream_flush(" << fpv << ")" << endl;
#endif

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
    
#ifdef debug
    cerr << "cppstream_close(" << fpv << ")" << endl;
#endif
    
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
        // The offset can't be determined, because this isn't a seekable stream.
        // Use a 0 offset.
        start_pos = 0;
        
        // TODO: There's no real way to prevent the hfile from seeking
        // internally in its buffer. 
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
        // The offset can't be determined, because this isn't a seekable stream.
        // Use a 0 offset.
        start_pos = 0;
        
        // TODO: There's no real way to prevent the hfile from seeking
        // internally in its buffer. 
        
    }
    fp->base.offset = start_pos;
    
    
    
    // Return the base hFILE*
    return &fp->base;
}

}

}
