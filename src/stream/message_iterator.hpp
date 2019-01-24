#ifndef VG_STREAM_MESSAGE_ITERATOR_HPP_INCLUDED
#define VG_STREAM_MESSAGE_ITERATOR_HPP_INCLUDED

/**
 * \file message_iterator.hpp
 * Defines a cursor for reading type-tagged, grouped binary messages from files.
 */

#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <functional>
#include <vector>

#include <google/protobuf/io/coded_stream.h>

#include "stream.hpp"
#include "blocked_gzip_input_stream.hpp"

namespace vg {

namespace stream {

using namespace std;


/**
 * Iterator over messages in VG-format files.
 * Also supports seeking and telling at the group level in bgzip files.
 * Cannot be copied, but can be moved.
 */
class MessageIterator {
public:
    /// Constructor to wrap a stream.
    MessageIterator(std::istream& in);
    
    /// Default constructor for an end iterator.
    MessageIterator() = default;
    
    /// Return true if dereferencing the iterator will produce a valid value, and false otherwise.
    bool has_next() const;
    
    /// Advance the iterator to the next message, or the end if this was the last message.
    void get_next();
    
    /// Get the current item
    const string& operator*() const;
    
    /// Take the current item, which must exist, and advance the iterator to the next one.
    string take();
    
    /// How many serialized, uncompressed bytes did the currently loaded item take.
    /// Will be 0 if no current item is available, but can also be 0 for valid, all-default items.
    size_t get_item_size() const;
    
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
    
    /// Return the raw virtual offset that the cursor is at in the file, or -1
    /// for an unseekable/untellable stream. Not necessarily at a group
    /// boundary, so cannot be used for seeking. Useful for getting the final
    /// virtual offset when the cursor hits the end of a file.
    /// This is NOT the virtual offset at which the currently loaded item occurs!
    int64_t tell_raw() const;
    
    /// Return the virtual offset of the currently loaded item, or tell_raw() if at end.
    /// Returns -1 for an unseekable/untellable stream.
    int64_t tell_item() const;
    
    /// Seek to the given virtual offset and read a single item from there. The
    /// next value produced will be the item that is there, and then the
    /// iterator will end (unless the user seeks somewhere else). Return false
    /// if seeking is unsupported or the seek fails.
    /// Will not work right if EOF is sought.
    bool seek_item_and_stop(int64_t virtual_offset);
    
    /// Returns iterators that act like begin() and end() for a stream containing messages
    static std::pair<MessageIterator, MessageIterator> range(std::istream& in);
    
private:
    
    string value;
    
    /// This holds the number of messages that exist in the current group.
    size_t group_count;
    /// This holds the number of messages read in the current group.
    size_t group_idx;
    /// This holds the virtual offset of the current group's start, or the
    /// number of the current group if seeking is not available.
    /// If the iterator is the end iterator, this is -1.
    int64_t group_vo;
    /// This holds the virtual offset of the current item, or -1 if seeking is not possible.
    /// Useful for seeking back to the item later, although you will have to seek to a group to iterate, after that.
    int64_t item_vo;
    /// This holds the number of serialized bytes that the currently loaded item took up in the stream.
    /// Doesn't count group size or item length overhead.
    size_t item_bytes;
    /// This is a flag for whether we should hit the end on the next get_next()
    /// It is set when we seek to a message individually, and unset when we seek to a group.
    bool end_next;
    
    /// Since these streams can't be copied or moved, we wrap ours in a uniqueptr_t so we can be moved.
    unique_ptr<BlockedGzipInputStream> bgzip_in;
    
    /// Make sure the given Protobuf-library bool return value is true, and fail otherwise with a message.
    void handle(bool ok);
};

}

}

#endif
