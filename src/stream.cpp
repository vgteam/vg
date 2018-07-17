#include "stream.hpp"

namespace stream {

using namespace std;

BlockedGzipOutputStream::BlockedGzipOutputStream(BGZF* bgzf_handle) : handle(bgzf_handle), buffer(), backed_up(0), byte_count(0) {
    // Nothing to do
}

BlockedGzipOutputStream::~BlockedGzipOutputStream() {
    // Make sure to finish writing before destructing.
    flush();
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


}

