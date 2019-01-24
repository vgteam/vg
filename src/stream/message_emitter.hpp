#ifndef VG_STREAM_MESSAGE_EMITTER_HPP_INCLUDED
#define VG_STREAM_MESSAGE_EMITTER_HPP_INCLUDED

/**
 * \file message_emitter.hpp
 * Defines an output cursor for writing grouped, type-tagged messages to BGZF files.
 */

#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <functional>
#include <vector>
#include <string>
#include <memory>

#include <google/protobuf/io/coded_stream.h>

#include "blocked_gzip_output_stream.hpp"

namespace vg {

namespace stream {

using namespace std;

/**
 * Class that wraps an output stream and allows emitting groups of binary
 * messages to it, with internal buffering. Handles finishing the file on its
 * own, and allows tracking of BGZF virtual offsets within a non-seekable
 * stream (as long as the entire stream is controlled by one instance). Cannot
 * be copied, but can be moved.
 *
 * Each group consists of:
 * - A 64-bit varint with the number of messages plus 1
 * - A 64-bit varint header/tag length
 * - Header/tag data
 * And then for each message item:
 * - A 64-bit varint message length
 * - Message data
 *
 * This format is designed to be syntactically the same as the old untagged VG
 * Protobuf format, to allow easy sniffing/reading of old files.
 *
 * Can call callbacks with the groups emitted and their virtual offsets, for
 * indexing purposes.
 * 
 * Note that the callbacks may be called by the object's destructor, so
 * anything they reference needs to outlive the object.
 *
 * Not thread-safe. May be more efficient than repeated write/write_buffered
 * calls because a single BGZF stream can be used.
 */
class MessageEmitter {
public:

    /// We refuse to serialize individual messages longer than this size.
    const size_t MAX_MESSAGE_SIZE = 1000000000;

    /// Constructor
    MessageEmitter(std::ostream& out, size_t max_group_size = 1000);
    
    /// Destructor that finishes the file
    ~MessageEmitter();
    
    // Prohibit copy
    MessageEmitter(const MessageEmitter& other) = delete;
    MessageEmitter& operator=(const MessageEmitter& other) = delete;
    // Allow default move
    MessageEmitter(MessageEmitter&& other) = default;
    MessageEmitter& operator=(MessageEmitter&& other) = default;
    
    /// Emit the given message with the given type tag.
    void write(const string& tag, string&& message);
    
    /// Emit a copy of the given message with the given type tag.
    /// To use when you have something you can't move.
    void write_copy(const string& tag, const string& message);
    
    /// Define a type for group emission event listeners.
    /// Arguments are: type tag, group buffer, start virtual offset, and past-end virtual offset.
    using listener_t = std::function<void(const string&, const vector<string>&, int64_t, int64_t)>;
    
    /// Add an event listener that listens for emitted groups. The listener
    /// will be called with the group buffer, the start virtual offset, and the
    /// past-end virtual offset. Moves the function passed in.
    /// Anything the function uses by reference must outlive this object!
    void on_group(listener_t&& listener);
    
    /// Actually write out everything in the buffer.
    /// Doesn't actually flush the underlying streams to disk.
    /// Assumes that no more than one group's worth of messages are in the buffer.
    void emit_group();

private:

    // This is our internal tag string for what is in our buffer.
    string group_tag;
    // This is our internal buffer
    vector<string> group;
    // This is how big we let it get before we dump it
    size_t max_group_size;
    // Since Protobuf streams can't be copied or moved, we wrap ours in a uniqueptr_t so we can be moved.
    unique_ptr<BlockedGzipOutputStream> bgzip_out;
    
    // If someone wants to listen in on emitted groups, they can register a handler
    vector<listener_t> group_handlers;

};

}

}

#endif
