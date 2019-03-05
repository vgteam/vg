/**
 * \file registry.cpp
 * Non-template implementations of registry functions.
 */

#include "registry.hpp"
#include "fdstream.hpp"

#include "register_loader_saver_gcsa.hpp"
#include "register_loader_saver_lcp.hpp"
#include "register_loader_saver_gbwt.hpp"
#include "register_loader_saver_xg.hpp"
#include "register_loader_saver_vg.hpp"
#include "register_loader_saver_minimizer.hpp"
#include "register_loader_saver_snarl_manager.hpp"

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/util/type_resolver.h>
#include <google/protobuf/util/type_resolver_util.h>
#include <google/protobuf/type.pb.h>
#include <google/protobuf/descriptor.pb.h>

#include <unistd.h>
#include <thread>

namespace vg {

namespace stream {

using namespace std;

auto Registry::register_everything() -> bool {

    // Register all the Protobufs
    register_protobuf<Alignment>("GAM");
    register_protobuf<MultipathAlignment>("MGAM");
    register_protobuf<Graph>("VG");
    register_protobuf<Snarl>("SNARL");
    register_protobuf<Pileup>("PILEUP");
    register_protobuf<Locus>("LOCUS");
    register_protobuf<Translation>("TRANS");

    // Register all the stream loader/savers.
    // These all call back to the registry.
    register_loader_saver_gcsa();
    register_loader_saver_lcp();
    register_loader_saver_gbwt();
    register_loader_saver_xg();
    register_loader_saver_vg();
    register_loader_saver_minimizer();
    register_loader_saver_snarl_manager();
    
    return true;
}

// Make sure the register_everything function is statically invoked.
static bool registration_success = Registry::register_everything();

auto Registry::get_tables() -> Tables& {
    // Keep state in a static function local, for on-demand initialization.
    static Tables tables;
    return tables;
}

auto Registry::is_valid_tag(const string& tag) -> bool {
    
    if (tag.size() > MAX_TAG_LENGTH) {
        // Too long to be correct, and definitely longer than we want to hash.
        return false;
    }
    
    // Get our state
    Tables& tables = get_tables();
    
    // Check if the tag is known in any table
    
    if (tables.tag_to_protobuf.count(tag)) {
        return true;
    }
    
    if (tables.tag_to_loader.count(tag)) {
        return true;
    }
    
    // If not, it may be just a raw Protobuf type name, like "vg.Type".
    if (tag.size() < 3 || tag[0] != 'v' || tag[1] != 'g' || tag[2] != '.') {
        // All the auto-assigned tags start with "vg.". So reject this as not an auto-assigned tag.
        return false;
    }
    
    // Sniff if it is known by Protobuf.
    
    ::google::protobuf::Type scratch;
    // The returned status is not really documented.
    // See <https://github.com/protocolbuffers/protobuf/blob/master/src/google/protobuf/stubs/status.h>
    // We need a leading slash or it can't find it.
    auto find_status = get_resolver().ResolveMessageType("/" + tag, &scratch);
    
#ifdef debug
    cerr << "Resolve " << tag << " status " << find_status.ok()  << " " << find_status.error_code() << " " << find_status.error_message() << endl;
#endif
    
    // If we resolved it OK, it can be a tag.
    // Otherwise, it is new and suspicious. Treat it as a message.
    // If it's really some new thing we haven't heard of, we'll hopefully crash.
    return find_status.ok();
}

auto Registry::get_resolver() -> ::google::protobuf::util::TypeResolver& {
    // We will have one shared resolver that we allocate here.
    // The unique_ptr makes sure it gets destructed at the end of the program.
    // And static local initialization is thread safe for free! See <https://stackoverflow.com/a/8102145>
    // Also see <https://stackoverflow.com/a/41651378> for how you even make a resolver.
    static unique_ptr<::google::protobuf::util::TypeResolver> resolver(
        ::google::protobuf::util::NewTypeResolverForDescriptorPool("", ::google::protobuf::DescriptorPool::generated_pool()));
        
    return *resolver;
}

auto wrap_bare_loader(function<void*(istream&)> istream_loader) -> load_function_t {
    // Capture the istream-using function by value
    return [istream_loader](const message_sender_function_t& for_each_message) -> void* {
    
        // Make a place to put the final void* result
        void* result = nullptr;
    
        // Open a pipe with an istream and an ostream.
        
        // First we make an array to hold input and output ends of a pipe.
        int pipe_fds[2];
        
        // Then we make the pipe with Unix magic.
        assert(pipe(pipe_fds) == 0);

#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Starting loader thread" << endl;
#endif
        
        // Start a thread to run the real loader function and put its result in our pointer
        thread run_loader([&]() {
        
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "Loader thread running" << endl;
#endif
        
            {
                // Wrap the read FD
                fdistream read_pipe(pipe_fds[0]);
                
                // Run the loader on that stream and save its result
                result = istream_loader(read_pipe);
            }
            
            // Close it up
            close(pipe_fds[0]);
            
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "Loader thread done" << endl;
#endif
        });
            
        {
            
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "Feeder thread running" << endl;
#endif
        
            {
                // Wrap the write FD
                fdostream write_pipe(pipe_fds[1]);
                
                // In our thread, loop over all the messages with for_each_message
                for_each_message([&](const string& message) {
                
#ifdef debug
                    #pragma omp critical (cerr)
                    cerr << "Feeder thread sending " << message.size() << " bytes" << endl;
#endif
                    
                    // When we get a message, write it to the loader's stream, blocking.
                    write_pipe << message;
                });
            }
            
            // Close the write FD so the reader sees EOF.
            close(pipe_fds[1]);
            
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "Feeder thread done" << endl;
#endif
            
        }
        
        // After we got all the messages, and sent all the data, wait for the loader thread to finish.
        run_loader.join();
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Load threads joined" << endl;
#endif
        
        // Return the contents of the shared pointer it should have written to.
        return result;
    };
}

void with_function_calling_stream(const message_consumer_function_t& emit_message, const function<void(ostream&)>& use_stream) {
    // Open a pipe with an istream and an ostream.
    
    // First we make an array to hold input and output ends of a pipe.
    int pipe_fds[2];
    
    // Then we make the pipe with Unix magic.
    assert(pipe(pipe_fds) == 0);
    
#ifdef debug
    #pragma omp critical (cerr)
    cerr << "Starting saver thread" << endl;
#endif
    
    // Start a thread to do the data schelpping
    thread run_chunker([&]() {
    
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Chunker thread running" << endl;
#endif

        // We don't want to emit an empty trailing chunk, but we should
        // emit an empty chunk if there is actually no data in the stream.
        bool emitted_anything = false;
    
        {
            // Wrap the read FD
            fdistream read_pipe(pipe_fds[0]);
            
            while (read_pipe) {
                // Until the read pipe EOFs or otherwise crashes
                
                // Collect chunks of a certain size
                // TODO: What size is best?
                string chunk(4096 * 1024, 0);
            
#ifdef debug
                #pragma omp critical (cerr)
                cerr << "Chunker thread waiting on pipe" << endl;
#endif
            
                // Read blockingly.
                // See https://stackoverflow.com/a/1816382
                read_pipe.read(&chunk[0], chunk.size());
                
                // Remember if we got all the characters or not.
                auto bytes_gotten = read_pipe.gcount();
                
#ifdef debug
                #pragma omp critical (cerr)
                cerr << "Chunker thread got " << bytes_gotten << " bytes" << endl;
#endif
                
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
            
        }
        
        close(pipe_fds[0]);
        
        if (!emitted_anything) {
            // Nothing came through before EOF.
            // We should still emit *a* message.
            emit_message("");
        }
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Chunker thread done" << endl;
#endif
    

    });
    
#ifdef debug
    #pragma omp critical (cerr)
    cerr << "Saver thread running" << endl;
#endif

    {
        // Wrap the write FD
        fdostream write_pipe(pipe_fds[1]);
        
        // Run the saver on that stream
        use_stream(write_pipe);
    }
    
    // Close the FD so the chunker thread can see EOF and make progress.
    close(pipe_fds[1]);
    
#ifdef debug
    #pragma omp critical (cerr)
    cerr << "Saver thread done" << endl;
#endif
    
    // Now both threads are done.
    run_chunker.join();
    
#ifdef debug
    #pragma omp critical (cerr)
    cerr << "Save threads joined" << endl;
#endif
    
}

auto wrap_bare_saver(function<void(const void*, ostream&)> ostream_saver) -> save_function_t {
    // Capture the ostream-using function by value
    return [ostream_saver](const void* to_save, const message_consumer_function_t& emit_message) {
   
        assert(to_save != nullptr);
        
        // Get ahold of an ostream that calls our emit_message function
        with_function_calling_stream(emit_message, [&ostream_saver, &to_save](ostream& out) {
            // And save to it
            ostream_saver(to_save, out);
        });
        
    };
}



}

}
