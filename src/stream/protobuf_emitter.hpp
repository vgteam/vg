#ifndef VG_STREAM_PROTOBUF_EMITTER_HPP_INCLUDED
#define VG_STREAM_PROTOBUF_EMITTER_HPP_INCLUDED

/**
 * \file protobuf_emitter.hpp
 * Defines an output cursor for writing Protobuf data to files.
 */

#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <functional>
#include <vector>

#include "blocked_gzip_output_stream.hpp"
#include "stream.hpp"

namespace vg {

namespace stream {

using namespace std;

/**
 *
 * Class that wraps an output stream and allows emitting groups of Protobuf
 * objects to it, with internal buffering. Handles finishing the file on its
 * own, and allows tracking of BGZF virtual offsets within a non-seekable
 * stream (as long as the entire stream is controlled by one instance). Cannot
 * be copied, but can be moved.
 *
 * Can call callbacks with the groups emitted and their virtual offsets, for
 * indexing purposes.
 * 
 * Note that the callbacks may be called by the ProtobufEmitter's destructor,
 * so anything they reference needs to outlive the ProtobufEmitter.
 *
 * Not thread-safe. May be more efficient than repeated write/write_buffered
 * calls because a single BGZF stream can be used.
 */
template <typename T>
class ProtobufEmitter {
public:
    /// Constructor
    ProtobufEmitter(std::ostream& out, size_t max_group_size = 1000) :
        group(),
        max_group_size(max_group_size),
        bgzip_out(new BlockedGzipOutputStream(out))
    {
#ifdef debug
        cerr << "Creating ProtobufEmitter" << endl;
#endif
        if (bgzip_out->Tell() == -1) {
            // Say we are starting at the beginnign of the stream, if we don't know where we are.
            bgzip_out->StartFile();
        }
    }
    
    /// Destructor that finishes the file
    ~ProtobufEmitter() {
#ifdef debug
        cerr << "Destroying ProtobufEmitter" << endl;
#endif
        if (bgzip_out.get() != nullptr) {
#ifdef debug
            cerr << "ProtobufEmitter emitting final group" << endl;
#endif
        
            // Before we are destroyed, write stuff out.
            emit_group();
            
#ifdef debug
            cerr << "ProtobufEmitter ending file" << endl;
#endif
            
            // Tell our stream to finish the file (since it hasn't been moved away)
            bgzip_out->EndFile();
        }
        
#ifdef debug
        cerr << "ProtobufEmitter destroyed" << endl;
#endif
    }
    
    // Prohibit copy
    ProtobufEmitter(const ProtobufEmitter& other) = delete;
    ProtobufEmitter& operator=(const ProtobufEmitter& other) = delete;
    // Allow default move
    ProtobufEmitter(ProtobufEmitter&& other) = default;
    ProtobufEmitter& operator=(ProtobufEmitter&& other) = default;
    
    /// Emit the given item.
    /// TODO: Not thread safe.
    void write(T&& item) {
        if (group.size() >= max_group_size) {
            emit_group();
        }
        group.emplace_back(std::move(item));
    }
    
    /// Emit a copy of the given item.
    /// To use when you have something you can't move.
    void write_copy(const T& item) {
        if (group.size() >= max_group_size) {
            emit_group();
        }
        group.push_back(item);
    }
    
    /// Define a type for group emission event listeners
    using listener_t = std::function<void(const vector<T>&, int64_t, int64_t)>;
    
    /// Add an event listener that listens for emitted groups. The listener
    /// will be called with the group buffer, the start virtual offset, and the
    /// past-end virtual offset. Moves the function passed in.
    /// Anything the function uses by reference must outlive this object!
    void on_group(listener_t&& listener) {
        group_handlers.emplace_back(std::move(listener));
    }
    
    /// Actually write out everything in the buffer.
    /// Doesn't actually flush the underlying streams to disk.
    /// Assumes that no more than one group's worht of items are in the buffer.
    void emit_group() {
        if (group.empty()) {
            // Nothing to do
            return;
        }
        
        // We can't write a non-empty buffer if our stream is gone/moved away
        assert(bgzip_out.get() != nullptr);
        
        auto handle = [](bool ok) {
            if (!ok) {
                throw std::runtime_error("stream::ProtobufEmitter::emit_group: I/O error writing protobuf");
            }
        };
    
        // Work out where the group we emit will start
        int64_t virtual_offset = bgzip_out->Tell();
    
        ::google::protobuf::io::CodedOutputStream coded_out(bgzip_out.get());

        // Prefix the group with the number of objects
        coded_out.WriteVarint64(group.size());
        handle(!coded_out.HadError());

        std::string s;
        size_t written = 0;
        for (auto& item : group) {
            handle(item.SerializeToString(&s));
            if (s.size() > MAX_PROTOBUF_SIZE) {
                throw std::runtime_error("stream::ProtobufEmitter::emit_group: message too large error writing protobuf");
            }
            
#ifdef debug
            cerr << "Writing message of " << s.size() << " bytes in group @ " << virtual_offset << endl;
#endif
            
            // And prefix each object with its size
            coded_out.WriteVarint32(s.size());
            handle(!coded_out.HadError());
            coded_out.WriteRaw(s.data(), s.size());
            handle(!coded_out.HadError());
        }
        
        // Work out where we ended
        coded_out.Trim();
        int64_t next_virtual_offset = bgzip_out->Tell();
        
        for (auto& handler : group_handlers) {
            // Report the group to each group handler that is listening
            handler(group, virtual_offset, next_virtual_offset);
        }
        
        // Empty the buffer because everything in it is written
        group.clear();
    }
    
    
private:

    // This is our internal buffer
    vector<T> group;
    // This is how big we let it get before we dump it
    size_t max_group_size;
    // Since Protobuf streams can't be copied or moved, we wrap ours in a uniqueptr_t so we can be moved.
    unique_ptr<BlockedGzipOutputStream> bgzip_out;
    
    // If someone wants to listen in on emitted groups, they can register a handler
    vector<listener_t> group_handlers;

};


/// Produce an std::function that can be invoked with Protobuf objects and save them to the given stream.
/// Easy way to get a dumping callback to feed to something that wants a callback.
/// The passed stream must outlive the resulting function.
template<typename Item>
std::function<void(const Item&)> emit_to(ostream& out) {
    // We are going to be clever and make a lambda capture a shared_ptr to an
    // emitter, so we can have the emitter last as long as the function we
    // return.
    shared_ptr<ProtobufEmitter<Item>> emitter(new ProtobufEmitter<Item>(out));

    return [emitter](const Item& item) {
        // Write out each item.
        // TODO: Set up so we can use the move operation the cursors support
        // Not easy because of https://stackoverflow.com/a/30394755
        emitter->write_copy(item);
    };
}


}

}

#endif
