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

#include "blocked_gzip_input_stream.hpp"

namespace vg {

namespace stream {

using namespace std;


/**
 * Iterator over messages in VG-format files. Yields pairs of string tag and
 * message data. Also supports seeking and telling at the group level in bgzip
 * files. Cannot be copied, but can be moved.
 *
 * TODO: Right now we always copy all messages into an internal buffer. We
 * should only do that if the message contents are actually accessed, with some
 * kind of lazy fake string that evaluates on string conversion. Then we could
 * more efficiently skip stuff with the wrong tag.
 */
class MessageIterator {
public:

    /// We refuse to serialize individual messages longer than this size.
    const size_t MAX_MESSAGE_SIZE = 1000000000;

    /// Constructor to wrap a stream.
    MessageIterator(istream& in);
    
    /// Represents a pair of a tag value and some message data.
    /// If there is no valid tag for a group, as given in the Registry, the tag will be "".
    using TaggedMessage = pair<string, string>;
    
    ///////////
    // C++ Iterator Interface
    ///////////
    
    /// Default constructor for an end iterator.
    MessageIterator() = default;
    
    /// Get the current item. Caller may move it away.
    /// Only legal to call if we are not an end iterator.
    TaggedMessage& operator*();
    
    /// Get the current item when we are const.
    /// Only legal to call if we are not an end iterator.
    const TaggedMessage& operator*() const;
    
    /// In-place pre-increment to advance the iterator.
    const MessageIterator& operator++();
    
    /// Check if two iterators are equal. Since you can only have one on
    /// a stream, this only has two equality classes: iterators that have hit
    /// the end, and iterators that haven't.
    bool operator==(const MessageIterator& other) const;
    
    /// Check if two iterators are not equal. Since you can only have one on
    /// a stream, this only has two equality classes: iterators that have hit
    /// the end, and iterators that haven't.
    bool operator!=(const MessageIterator& other) const;
    
    /// Returns iterators that act like begin() and end() for a stream containing messages
    static std::pair<MessageIterator, MessageIterator> range(istream& in);
    
    ///////////
    // has_next()/take() interface
    ///////////
    
    /// Return true if dereferencing the iterator will produce a valid value, and false otherwise.
    bool has_next() const;
    
    /// Advance the iterator to the next message, or the end if this was the last message.
    /// Basically the same as ++.
    void get_next();
    
    /// Take the current item, which must exist, and advance the iterator to the next one.
    TaggedMessage take();
    
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
    
    /// Holds the most recently pulled-out message tag and value.
    /// May get moved away.
    TaggedMessage value;
    
    /// Because the whole value pair may get moved away, we keep a backup copy of the tag and replace it.
    /// TODO: This is a bit of a hack.
    string backup_tag;
    
    /// This holds the number of messages that exist in the current group.
    /// Counts the tag, if present.
    size_t group_count;
    /// This holds the number of messages read in the current group.
    /// Counts the tag, if present.
    size_t group_idx;
    /// This holds the virtual offset of the current group's start, or the
    /// number of the current group if seeking is not available.
    /// If the iterator is the end iterator, this is -1.
    int64_t group_vo;
    /// This holds the virtual offset of the current item, or counts up through the group if seeking is not possible.
    /// Useful for seeking back to the item later, although you will have to seek to a group to iterate, after that.
    int64_t item_vo;
    
    /// Since these streams can't be copied or moved, we wrap ours in a uniqueptr_t so we can be moved.
    unique_ptr<BlockedGzipInputStream> bgzip_in;
    
    /// Make sure the given Protobuf-library bool return value is true, and fail otherwise with a message.
    void handle(bool ok);
};

}

}

#endif
