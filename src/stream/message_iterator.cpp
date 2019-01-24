/**
 * \file message_iterator.cpp
 * Implementations for the MessageIterator for reading type-tagged grouped message files
 */

#include "message_iterator.hpp"

namespace vg {

namespace stream {

using namespace std;

MessageIterator::MessageIterator(std::istream& in) :
    group_count(0),
    group_idx(0),
    group_vo(-1),
    item_vo(-1),
    end_next(false),
    bgzip_in(new BlockedGzipInputStream(in))
{
    get_next();
}

const string& MessageIterator::operator*() const {
    return value;
}

string& MessageIterator::operator*() {
    return value;
}


const MessageIterator& MessageIterator::operator++() {
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
            value = string();
            return *this;
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
    
    if (msgSize > MAX_PROTOBUF_SIZE) {
        throw std::runtime_error("[stream::MessageIterator::get_next] protobuf message of " +
                                 std::to_string(msgSize) + " bytes is too long");
    }
    
    
    // We have a message.
    value.clear();
    if (msgSize) {
        // It has non-default field values. Parse them.
        handle(coded_in.ReadString(&value, msgSize));
    } 
    
    // Move on to the next message in the group
    group_idx++;
    
    // Return ourselves, after increment
    return *this;
}

bool MessageIterator::operator==(const MessageIterator& other) const {
    // Just ask if we both agree on whether we hit the end.
    return has_next() == other.has_next();
}
    
bool MessageIterator::operator!=(const MessageIterator& other) const {
    // Just ask if we disagree on whether we hit the end.
    return has_next() != other.has_next();
}

bool MessageIterator::has_next() const {
    return item_vo != -1;
}

void MessageIterator::get_next() {
    // Run increment but don't return anything.
    ++(*this);
}

string MessageIterator::take() {
    string temp = std::move(value);
    get_next();
    // Return by value, which gets moved.
    return temp;
}

int64_t MessageIterator::tell_group() const {
    if (bgzip_in->Tell() != -1) {
        // The backing file supports seek/tell (which we ascertain by attempting it).
        if (group_vo == -1) {
            // We hit EOF and have no loaded message
            return bgzip_in->Tell();
        } else {
            // Return the *group's* virtual offset (not the current one)
            return group_vo;
        }
    } else {
        // group_vo holds a count. But we need to say we can't seek.
        return -1;
    }
}

bool MessageIterator::seek_group(int64_t virtual_offset) {
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

pair<MessageIterator, MessageIterator> MessageIterator::range(istream& in) {
    return make_pair(MessageIterator(in), MessageIterator());
}

void MessageIterator::handle(bool ok) {
    if (!ok) {
        throw std::runtime_error("[stream::MessageIterator] obsolete, invalid, or corrupt protobuf input");
    }
}

}

}
