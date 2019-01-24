#ifndef VG_STREAM_PROTOBUF_ITERATOR_HPP_INCLUDED
#define VG_STREAM_PROTOBUF_ITERATOR_HPP_INCLUDED

/**
 * \file protobuf_iterator.hpp
 * Defines a cursor for reading Protobuf messages from files.
 *
 * TODO: This implementation copies messages. Go back to/find a way to do zero-copy stream parsing directly.
 */

#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <functional>
#include <vector>

#include "message_iterator.hpp"
#include "registry.hpp"
#include "../json2pb.h"

namespace vg {

namespace stream {

using namespace std;


/**
 * Refactored stream::for_each function that follows the unidirectional iterator interface.
 * Also supports seeking and telling at the group level in bgzip files.
 * Cannot be copied, but can be moved.
 */
template <typename T>
class ProtobufIterator {
public:
    /// Constructor
    ProtobufIterator(istream& in);
    
    ///////////
    // C++ Iterator Interface
    ///////////
    
    /// Default constructor for an end iterator.
    ProtobufIterator() = default;
    
    /// Get the current item. Caller may move it away.
    /// Only legal to call if we are not an end iterator.
    T& operator*();
    
    /// Get the current item when we are const.
    /// Only legal to call if we are not an end iterator.
    const T& operator*() const;
    
    /// In-place pre-increment to advance the iterator.
    const ProtobufIterator<T>& operator++();
    
    /// Check if two iterators are equal. Since you can only have one on
    /// a stream, this only has two equality classes: iterators that have hit
    /// the end, and iterators that haven't.
    bool operator==(const ProtobufIterator<T>& other) const;
    
    /// Check if two iterators are not equal. Since you can only have one on
    /// a stream, this only has two equality classes: iterators that have hit
    /// the end, and iterators that haven't.
    bool operator!=(const ProtobufIterator<T>& other) const;
    
    /// Returns iterators that act like begin() and end() for a stream containing messages
    static pair<ProtobufIterator<T>, ProtobufIterator<T>> range(istream& in);
    
    ///////////
    // has_next()/take() interface
    ///////////
    
    /// Return true if dereferencing the iterator will produce a valid value, and false otherwise.
    bool has_next() const;
    
    /// Advance the iterator to the next message, or the end if this was the last message.
    /// Basically the same as ++.
    void get_next();
    
    /// Take the current item, which must exist, and advance the iterator to the next one.
    T take();
    
    ///////////
    // File position and seeking
    ///////////
    
    /// Return the virtual offset of the group being currently read (i.e. the
    /// group to which the current message belongs), to seek back to. You can't
    /// seek back to the current message, just to the start of the group.
    /// Returns -1 instead if the underlying file doesn't support seek/tell.
    /// Returns the past-the-end virtual offset of the file if EOF is reached.
    int64_t tell_group() const;
    
    /// Seek to the given virtual offset and start reading the group that is there.
    /// The next value produced will be the first value in that group.
    /// If already at the start of the group at the given virtual offset, does nothing.
    /// Return false if seeking is unsupported or the seek fails.
    bool seek_group(int64_t virtual_offset);
    
private:
    
    /// Wrap a MessageIterator and just do Protobuf parsing on top of that.
    MessageIterator message_it;
    
    /// We always maintain a parsed version of the current message.
    T value;
    
    /// Fill in value, if message_it has a value of an appropriate tag.
    void fill_value();
};

///////////
// Template implementations
///////////

template<typename T>
ProtobufIterator<T>::ProtobufIterator(istream& in) : message_it(in) {
    // Make sure to fill in our value field.
    fill_value();
}

template<typename T>
auto ProtobufIterator<T>::operator*() -> T& {
    return value;
}

template<typename T>
auto ProtobufIterator<T>::operator*() const -> const T& {
    return value;
}

template<typename T>
auto ProtobufIterator<T>::operator++() -> const ProtobufIterator<T>& {
    ++message_it;
    fill_value();
    return *this;
}

template<typename T>
auto ProtobufIterator<T>::operator==(const ProtobufIterator<T>& other) const -> bool {
    return message_it == other.message_it;
}

template<typename T>
auto ProtobufIterator<T>::operator!=(const ProtobufIterator<T>& other) const -> bool {
    return message_it != other.message_it;
}

template<typename T>
auto ProtobufIterator<T>::range(istream& in) -> pair<ProtobufIterator<T>, ProtobufIterator<T>> {
    return make_pair(ProtobufIterator<T>(in), ProtobufIterator<T>());
}

template<typename T>
auto ProtobufIterator<T>::has_next() const -> bool {
    return message_it.has_next();
}

template<typename T>
auto ProtobufIterator<T>::get_next() -> void {
    ++(*this);
}

template<typename T>
auto ProtobufIterator<T>::take() -> T {
    auto temp = std::move(value);
    get_next();
    // Return by value, which gets moved.
    return temp;
}

template<typename T>
auto ProtobufIterator<T>::tell_group() const -> int64_t {
    return message_it.tell_group();
}

template<typename T>
auto ProtobufIterator<T>::seek_group(int64_t virtual_offset) -> bool {
    if (message_it.seek_group(virtual_offset)) {
        fill_value();
        return true;
    }
    return false;
}

template<typename T>
auto ProtobufIterator<T>::fill_value() -> void {
    // This is where the magic happens.
    // We have already advanced or EOF'd our message iterator.
    // We need to check tags and fill in our message value.
    
#ifdef debug
    cerr << "Fill Protobuf value" << endl;
#endif
    
    if (message_it.has_next()) {
        // Grab the tag and message
        auto& tag_and_message = *message_it;
        auto& tag = tag_and_message.first;
        auto& message = tag_and_message.second;
        
        // See if the tag is valid for what we want to parse.
        // TODO: Do this in a way where we can check this only per-group!
        if (!Registry::check_protobuf_tag<T>(tag)) {
            // The registry doesn't think this tag is legit for what we are parsing.
            // Abort!
            throw runtime_error("[stream::ProtobufIterator] tag \"" + tag +
                "\" for Protobuf that should be \"" + Registry::get_protobuf_tag<T>() + "\"");
                
            // TODO: Skip over these instead of aborting to allow for multiplexing (i.e. VG in with XG)
        }
        
        // Parse the value
        value.ParseFromString(message);
        
#ifdef debug   
        cerr << "Got " << pb2json(value) << " from " << message.size() << " bytes" << endl;
#endif
    } else {
        // Don't waste space. Clean up any old value.
        value.Clear();
    }
}


}

}

#endif
