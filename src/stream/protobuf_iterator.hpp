#ifndef VG_STREAM_PROTOBUF_ITERATOR_HPP_INCLUDED
#define VG_STREAM_PROTOBUF_ITERATOR_HPP_INCLUDED

/**
 * \file protobuf_iterator.hpp
 * Defines a cursor for reading Protobuf messages from files.
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
 * Refactored stream::for_each function that follows the unidirectional iterator interface.
 * Also supports seeking and telling at the group level in bgzip files.
 * Cannot be copied, but can be moved.
 */
template <typename T>
class ProtobufIterator {
public:
    /// Constructor
    ProtobufIterator(std::istream& in) :
        group_count(0),
        group_idx(0),
        group_vo(-1),
        item_vo(-1),
        item_bytes(0),
        end_next(false),
        bgzip_in(new BlockedGzipInputStream(in))
    {
        get_next();
    }
    
    /// Return true if dereferencing the iterator will produce a valid value, and false otherwise.
    inline bool has_next() const {
        return item_vo != -1;
    }
    
    /// Advance the iterator to the next message, or the end if this was the last message.
    void get_next() {
        if (end_next || group_count == group_idx) {
            // We have made it to the end of the group we are reading. We will
            // start a new group now.
            
            // Determine exactly where we are positioned, if possible, before
            // creating the CodedInputStream to read the group's item count
            auto virtual_offset = bgzip_in->Tell();
            
            if (virtual_offset == -1) {
                // We don't have seek capability, so we just count up the groups we read.
                // On construction this is -1; bump it up to 0 for the first group.
                group_vo++;
            } else {
                // We can seek. We need to know what offset we are at
                group_vo = virtual_offset;
            }
            
            // Start at the start of the new group
            group_idx = 0;
            
            // Make a CodedInputStream to read the group length
            ::google::protobuf::io::CodedInputStream coded_in(bgzip_in.get());
            // Alot space for group's length (generously)
            coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
            
            // Try and read the group's length
            if (end_next || !coded_in.ReadVarint64((::google::protobuf::uint64*) &group_count)) {
                // We didn't get a length (or we want to end the iteration)
                
                // This is the end of the input stream, switch to state that
                // will match the end constructor
                group_vo = -1;
                item_vo = -1;
                item_bytes = 0;
                value = T();
                return;
            }
            
        }
        
        // Now we know we're in a group.
        
        // Get the item's virtual offset, if available
        auto virtual_offset = bgzip_in->Tell();
        
        // We need a fresh CodedInputStream every time, because of the total byte limit
        ::google::protobuf::io::CodedInputStream coded_in(bgzip_in.get());
        // Alot space for size and item (generously)
        coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
        
        // A message starts here
        if (virtual_offset == -1) {
            // Just track the counter.
            item_vo++;
        } else {
            // We know where here is
            item_vo = virtual_offset;
        }
        
        // The messages are prefixed by their size
        uint32_t msgSize = 0;
        handle(coded_in.ReadVarint32(&msgSize));
        item_bytes = msgSize;
        
        if (msgSize > MAX_PROTOBUF_SIZE) {
            throw std::runtime_error("[stream::ProtobufIterator::get_next] protobuf message of " +
                                     std::to_string(msgSize) + " bytes is too long");
        }
        
        
        // We have a message.
        value.Clear();
        if (msgSize) {
            // It has non-default field values. Parse them.
            std::string s;
            handle(coded_in.ReadString(&s, msgSize));
            handle(value.ParseFromString(s));
        } 
        
        // Move on to the next message in the group
        group_idx++;
    }
    
//    inline ProtobufIterator<T> operator++( int ) {
//        ProtobufIterator<T> temp = *this;
//        get_next();
//        return temp;
//    }
    
    inline const T& operator*() const {
        return value;
    }
    
    /// Take the current item, which must exist, and advance the iterator to the next one.
    inline T take() {
        T temp = std::move(value);
        get_next();
        // Return by value, which gets moved.
        return temp;
    }
    
    /// How many serialized, uncompressed bytes did the currently loaded item take.
    /// Will be 0 if no current item is available, but can also be 0 for valid, all-default items.
    inline size_t get_item_size() const {
        return item_bytes;
    }
    
    /// Return the virtual offset of the group being currently read (i.e. the
    /// group to which the current message belongs), to seek back to. You can't
    /// seek back to the current message, just to the start of the group.
    /// Returns -1 instead if the underlying file doesn't support seek/tell.
    /// Returns the past-the-end virtual offset of the file if EOF is reached.
    inline int64_t tell_group() const {
        if (bgzip_in->Tell() != -1) {
            // The backing file supports seek/tell (which we ascertain by attempting it).
            if (group_vo == -1) {
                // We hit EOF and have no loaded message
                return tell_raw();
            } else {
                // Return the *group's* virtual offset (not the current one)
                return group_vo;
            }
        } else {
            // group_vo holds a count. But we need to say we can't seek.
            return -1;
        }
    }
    
    /// Seek to the given virtual offset and start reading the group that is there.
    /// The next value produced will be the first value in that group.
    /// If already at the start of the group at the given virtual offset, does nothing.
    /// Return false if seeking is unsupported or the seek fails.
    inline bool seek_group(int64_t virtual_offset) {
        if (virtual_offset < 0) {
            // That's not allowed
            return false;
        }
        
        if (group_idx == 0 && group_vo == virtual_offset && !end_next) {
            // We are there already
            return true;
        }
        
        // Try and do the seek
        bool sought = bgzip_in->Seek(virtual_offset);
        
        if (!sought) {
            // We can't seek
            return false;
        }
        
        // Get ready to read the group that's here
        group_count = 0;
        group_idx = 0;
        end_next = false;
        
        // Read it (or detect EOF)
        get_next();
        
        // It worked!
        return true;
    }
    
    /// Return the raw virtual offset that the cursor is at in the file, or -1
    /// for an unseekable/untellable stream. Not necessarily at a group
    /// boundary, so cannot be used for seeking. Useful for getting the final
    /// virtual offset when the cursor hits the end of a file.
    /// This is NOT the virtual offset at which the currently loaded item occurs!
    inline int64_t tell_raw() const {
        return bgzip_in->Tell();
    }
    
    /// Return the virtual offset of the currently loaded item, or tell_raw() if at end.
    /// Returns -1 for an unseekable/untellable stream.
    inline int64_t tell_item() const {
        if (bgzip_in->Tell() != -1) {
            // The backing file supports seek/tell (which we ascertain by attempting it).
            if (item_vo == -1) {
                // We hit EOF and have no loaded message
                return tell_raw();
            } else {
                // Return the item's virtual offset
                return item_vo;
            }
        } else {
            // item_vo holds a count. But we need to say we can't seek.
            return -1;
        }
    }
    
    /// Seek to the given virtual offset and read a single item from there. The
    /// next value produced will be the item that is there, and then the
    /// iterator will end (unless the user seeks somewhere else). Return false
    /// if seeking is unsupported or the seek fails.
    /// Will not work right if EOF is sought.
    inline bool seek_item_and_stop(int64_t virtual_offset) {
        if (virtual_offset < 0) {
            // That's not allowed
            return false;
        }
        
        // Try and do the seek
        bool sought = bgzip_in->Seek(virtual_offset);
        
        if (!sought) {
            // We can't seek
            return false;
        }
        
        // Pretend to be in a 1-element group, which will be the last one until someone seeks again.
        group_count = 1;
        group_idx = 0;
        
        // Allow an item to be read
        end_next = false;
        
        // Read the element
        get_next();
        
        // Stop after this one
        end_next = true;
        
        // It worked!
        return true;
    }
    
    /// Returns iterators that act like begin() and end() for a stream containing protobuf data
    static std::pair<ProtobufIterator<T>, ProtobufIterator<T>> range(std::istream& in) {
        return std::make_pair(ProtobufIterator<T>(in), ProtobufIterator<T>());
    }
    
private:
    
    T value;
    
    // This holds the number of messages that exist in the current group.
    size_t group_count;
    // This holds the number of messages read in the current group.
    size_t group_idx;
    // This holds the virtual offset of the current group's start, or the
    // number of the current group if seeking is not available.
    // If the iterator is the end iterator, this is -1.
    int64_t group_vo;
    // This holds the virtual offset of the current item, or -1 if seeking is not possible.
    // Useful for seeking back to the item later, although you will have to seek to a group to iterate, after that.
    int64_t item_vo;
    // This holds the number of serialized bytes that the currently loaded item took up in the stream.
    // Doesn't count group size or item length overhead.
    size_t item_bytes;
    // This is a flag for whether we should hit the end on the next get_next()
    // It is set when we seek to a message individually, and unset when we seek to a group.
    bool end_next;
    
    // Since Protobuf streams can't be copied or moved, we wrap ours in a uniqueptr_t so we can be moved.
    unique_ptr<BlockedGzipInputStream> bgzip_in;
    
    void handle(bool ok) {
        if (!ok) {
            throw std::runtime_error("[stream::ProtobufIterator] obsolete, invalid, or corrupt protobuf input");
        }
    }
};

}

}

#endif
