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
#include <list>
#include <mutex>

#include "message_emitter.hpp"
#include "registry.hpp"
#include "../json2pb.h"

namespace vg {

namespace stream {

using namespace std;

/**
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
 * May be more efficient than repeated write/write_buffered calls because a
 * single BGZF stream can be used.
 *
 * Thread-safe to call into. Serialization is done before locking. If a
 * particular order is needed between objects, use the multi-object write
 * functions. Listeners will be called inside the lock, so only one will be in
 * progress at a time.
 */
template <typename T>
class ProtobufEmitter {
public:
    /// Constructor
    ProtobufEmitter(std::ostream& out, size_t max_group_size = 1000);
    
    /// Destructor that finishes the file
    ~ProtobufEmitter();
    
    // Prohibit copy
    ProtobufEmitter(const ProtobufEmitter& other) = delete;
    ProtobufEmitter& operator=(const ProtobufEmitter& other) = delete;
    // Allow default move
    ProtobufEmitter(ProtobufEmitter&& other) = default;
    ProtobufEmitter& operator=(ProtobufEmitter&& other) = default;
    
    /// Emit the given item.
    /// TODO: May not really be any more efficient.
    /// We serialize to string right away in either case.
    void write(T&& item);
    
    /// Emit the given collection of items in order, with no other intervening
    /// items between them.
    void write_many(vector<T>&& ordered_items);
    
    /// Emit a copy of the given item.
    /// To use when you have something you can't move.
    void write_copy(const T& item);
    
    /// Define a type for group emission event listeners.
    /// The arguments are the start virtual offset and the past-end virtual offset.
    using group_listener_t = std::function<void(int64_t, int64_t)>;
    
    /// Add an event listener that listens for emitted groups. The listener
    /// will be called with the start virtual offset, and the
    /// past-end virtual offset. Moves the function passed in.
    /// Anything the function uses by reference must outlive this object!
    void on_group(group_listener_t&& listener);
    
    /// Define a type for message emission event listeners.
    /// This gets called for every message we emit, and then the group listeners get called for the whole group.
    using message_listener_t = std::function<void(const T&)>;
    
    /// Add an event listener that will be called every time a message is emitted.
    void on_message(message_listener_t&& listener);
    
    /// Actually write out everything in the buffer.
    /// Doesn't actually flush the underlying streams to disk.
    /// Assumes that no more than one group's worht of items are in the buffer.
    void emit_group();
    
private:

    /// Mutex to controll access to the backing MessageEmitter.
    /// Also needs to control access to the listener lists.
    mutex out_mutex;

    /// We wrap a MessageEmitter that handles tagged message IO
    MessageEmitter message_emitter;
    
    /// And a single precomputed copy of the tag string to use
    string tag;
    
    /// And all the group handler functions. These need to never move; they are
    /// captured by reference to listeners in our MessageEmitter.
    list<group_listener_t> group_handlers;
    
    /// These we invoke ourselves per message.
    vector<message_listener_t> message_handlers;
    
    /// Make sure the given Protobuf-library bool return value is true, and fail otherwise with a message.
    void handle(bool ok);

};

/// Produce an std::function that can be invoked with Protobuf objects and save them to the given stream.
/// Easy way to get a dumping callback to feed to something that wants a callback.
/// The passed stream must outlive the resulting function.
template<typename Item>
std::function<void(const Item&)> emit_to(ostream& out);

/////////
// Template implementations
/////////

template<typename T>
ProtobufEmitter<T>::ProtobufEmitter(std::ostream& out, size_t max_group_size) : message_emitter(out, max_group_size),
    tag(Registry::get_protobuf_tag<T>()) {
    
    // Nothing to do!
}

template<typename T>
ProtobufEmitter<T>::~ProtobufEmitter() {
#ifdef debug
    cerr << "Destroying ProtobufEmitter" << endl;
#endif
    
    // Emit the final group, so the MessageEmitter is empty when it destructs
    // and doesn't try to call any callbacks.
    // TODO: The whole callback ownership system is weird and should be re-done better somehow.
    emit_group();
    
#ifdef debug
    cerr << "ProtobufEmitter destroyed" << endl;
#endif
}

template<typename T>
auto ProtobufEmitter<T>::write(T&& item) -> void {
    // Grab the item
    T to_encode = std::move(item);
    
    // Encode it to a string
    string encoded;
    handle(to_encode.SerializeToString(&encoded));
    
    // Lock the backing emitter
    lock_guard<mutex> lock(out_mutex);
    
    // Write it with the correct tag.
    message_emitter.write(tag, std::move(encoded));
    
    for (auto& handler : message_handlers) {
        // Fire the handlers in serial
        handler(to_encode);
    }
}

template<typename T>
auto ProtobufEmitter<T>::write_many(vector<T>&& ordered_items) -> void {
    // Grab the items
    vector<T> to_encode = std::move(ordered_items);
    
    // Encode them all to strings
    vector<string> encoded(to_encode.size());
    for (size_t i = 0; i < to_encode.size(); i++) {
        handle(to_encode[i].SerializeToString(&encoded[i]));
    }
    
    // Lock the backing emitter
    lock_guard<mutex> lock(out_mutex);
    
    for (size_t i = 0; i < to_encode.size(); i++) {
        // Write each message with the correct tag.
        message_emitter.write(tag, std::move(encoded[i]));
        
        for (auto& handler : message_handlers) {
            // Fire the handlers in serial
            handler(to_encode[i]);
        }
    }
    
}

template<typename T>
auto ProtobufEmitter<T>::write_copy(const T& item) -> void {
    // Encode it to a string
    string encoded;
    handle(item.SerializeToString(&encoded));
    
#ifdef debug
    cerr << "Write Protobuf: " << pb2json(item) << " to " << encoded.size() << " bytes" << endl;
#endif

    // Lock the backing emitter
    lock_guard<mutex> lock(out_mutex);
    
    // Write it with the correct tag.
    message_emitter.write(tag, std::move(encoded));
    
    for (auto& handler : message_handlers) {
        // Fire the handlers in serial
        handler(item);
    }
}

template<typename T>
auto ProtobufEmitter<T>::on_group(group_listener_t&& listener) -> void {
    // Lock the handler list
    lock_guard<mutex> lock(out_mutex);
    
    // Take custody
    group_handlers.emplace_back(std::move(listener));
    
    // Grab a reference
    auto& owned_listener = group_handlers.back();
    
    // Capture by reference in another listener.
    // TODO: This isn't going to work at all if we want to ever use the same MessageEmitter with multiple ProtobufEmitters...
    message_emitter.on_group([&owned_listener](const string& tag, int64_t start_vo, int64_t past_end_vo) {
        // Call back with the group info.
        owned_listener(start_vo, past_end_vo);
    });
}

template<typename T>
auto ProtobufEmitter<T>::on_message(message_listener_t&& listener) -> void {
    // Lock the handler list
    lock_guard<mutex> lock(out_mutex);

    // Put in the collection to loop through on every message.
    message_handlers.emplace_back(std::move(listener));
}

template<typename T>
auto ProtobufEmitter<T>::emit_group() -> void {
    // Lock the backing emitter
    lock_guard<mutex> lock(out_mutex);

    message_emitter.emit_group();
}

template<typename T>
auto ProtobufEmitter<T>::handle(bool ok) -> void {
    if (!ok) {
        throw std::runtime_error("stream::ProtobufEmitter: could not write Protobuf");
    }
}

template<typename Item>
auto emit_to(ostream& out) -> std::function<void(const Item&)> {
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
