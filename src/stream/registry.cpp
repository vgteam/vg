/**
 * \file registry.cpp
 * Non-template implementations of registry functions.
 */

#include "registry.hpp"
#include "fd_streams.hpp"

#include <unistd.h>

namespace vg {

namespace stream {

using namespace std;

load_function_t wrap_stream_loader(function<void*(istream&)> istream_loader) {
    // Capture the istream-using function by value
    return [istream_loader](const message_sender_function_t& for_each_message) -> void* {
    
        // Make a place to put the final void* result
        void* result = nullptr;
    
        // Open a pipe with an istream and an ostream.
        
        // First we make an array to hold input and output ends of a pipe.
        int pipe_fds[2];
        
        // Then we make the pipe with Unix magic.
        assert(pipe(pipe_fds) == 0);
        
        {
            // Then wrap each side in a C++ stream.
            // We had to write these ourselves.
            FDIstream read_pipe(pipe_fds[0]);
            FDOstream write_pipe(pipe_fds[1]);
            
            // Start a thread to run the real loader function and put its result in our pointer
            
            // In our thread, loop over all the messages with for_each_message
            
            // When we get a message, write it to the loader's stream, blocking.
            
            // After we got all the messages, wait for the loader thread to finish.
        }
        
        // Close up the pipe
        close(pipe_fds[0]);
        close(pipe_fds[1]);
        
        // Return the contents of the shared pointer it should have written to.
        return result;
    };
}


auto Registry::get_tables() -> Tables& {
    // Keep state in a static function local, for on-demand initialization.
    static Tables tables;
    return tables;
}

}

}
