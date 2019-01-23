/**
 * \file registry.cpp
 * Non-template implementations of registry functions.
 */

#include "registry.hpp"
#include "fd_streams.hpp"

#include <unistd.h>
#include <thread>

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
        
        // Start a thread to run the real loader function and put its result in our pointer
        thread run_loader([&]() {
            // Wrap the read FD
            fdistream read_pipe(pipe_fds[0]);
            
            // Run the loader on that stream and save its result
            result = istream_loader(read_pipe);
        });
            
        {
            // Wrap the write FD
            fdostream write_pipe(pipe_fds[1]);
            
            // In our thread, loop over all the messages with for_each_message
            for_each_message([&](const string& message) {
                // When we get a message, write it to the loader's stream, blocking.
                write_pipe << message;
            });
            
        }
        
        // After we got all the messages, and sent all the data, wait for the loader thread to finish.
        run_loader.join();
        
        // Close up the pipe after the streams are gone
        close(pipe_fds[0]);
        close(pipe_fds[1]);
        
        // Return the contents of the shared pointer it should have written to.
        return result;
    };
}

save_function_t wrap_stream_saver(function<void(const void*, ostream&)> ostream_saver) {
    // Capture the ostream-using function by value
    return [ostream_saver](const void* to_save, const message_consumer_function_t& emit_message) {
    
        // Open a pipe with an istream and an ostream.
        
        // First we make an array to hold input and output ends of a pipe.
        int pipe_fds[2];
        
        // Then we make the pipe with Unix magic.
        assert(pipe(pipe_fds) == 0);
        
        // Start a thread to run the real saver function and write to the pipe
        thread run_saver([&]() {
            // Wrap the write FD
            fdostream write_pipe(pipe_fds[1]);
            
            // Run the saver on that stream
            ostream_saver(to_save, write_pipe);
        });
            
        {
            // Wrap the read FD
            fdistream read_pipe(pipe_fds[0]);
            
            // We don't want to emit an empty trailing chunk, but we should
            // emit an empty chunk if there is actually no data in the stream.
            bool emitted_anything = false;
            
            while (read_pipe) {
                // Until the read pipe EOFs or otherwise crashes
                
                // Collect chunks of a certain size
                // TODO: What size is best?
                string chunk(4096, 0);
            
                // Read blockingly.
                // See https://stackoverflow.com/a/1816382
                read_pipe.read(&chunk[0], chunk.size());
                
                // Remember if we got all the characters or not.
                auto bytes_gotten = read_pipe.gcount();
                
                if (bytes_gotten < chunk.size()) {
                    // Shrink to fit
                    chunk.resize(bytes_gotten);
                }
                
                if (bytes_gotten != 0) {
                    // Emit some bytes
                    emit_message(chunk);
                    emitted_anything = true;
                }
            }
            
            if (!emitted_anything) {
                // Nothing came through before EOF.
                // We should still emit *a* message.
                emit_message("");
            }
        }
        
        // After we got all the data and sent all the messages, wait for the
        // writer thread to do whatever it does after writing everything.
        run_saver.join();
        
        // Close up the pipe after the streams are gone
        close(pipe_fds[0]);
        close(pipe_fds[1]);
    };
}


auto Registry::get_tables() -> Tables& {
    // Keep state in a static function local, for on-demand initialization.
    static Tables tables;
    return tables;
}

}

}
