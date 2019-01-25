/**
 * \file message_iterator.cpp
 * Implementations for the MessageIterator for reading type-tagged grouped message files
 */

#include "message_iterator.hpp"
#include "registry.hpp"

namespace vg {

namespace stream {

using namespace std;

MessageIterator::MessageIterator(istream& in) :
    value(),
    backup_tag(),
    group_count(0),
    group_idx(0),
    group_vo(-1),
    item_vo(-1),
    bgzip_in(new BlockedGzipInputStream(in))
{
    get_next();
}

auto MessageIterator::operator*() const -> const TaggedMessage& {
    return value;
}

auto MessageIterator::operator*() -> TaggedMessage& {
    return value;
}


auto MessageIterator::operator++() -> const MessageIterator& {
    if (group_count == group_idx) {
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
        // Alot space for group's length, tag's length, and tag (generously)
        coded_in.SetTotalBytesLimit(MAX_MESSAGE_SIZE * 2, MAX_MESSAGE_SIZE * 2);
        
        // Try and read the group's length
        if (!coded_in.ReadVarint64((::google::protobuf::uint64*) &group_count)) {
            // We didn't get a length
            
            // This is the end of the input stream, switch to state that
            // will match the end constructor
            group_vo = -1;
            item_vo = -1;
            value.first.clear();
            value.second.clear();
            return *this;
        }
        
        // Now we have to grab the tag, which is pretending to be the first item.
        // It could also be the first item, if it isn't a known tag string.
        
        // Get the tag's virtual offset, if available
        virtual_offset = bgzip_in->Tell();
        
        // The tag is prefixed by its size
        uint32_t tagSize = 0;
        handle(coded_in.ReadVarint32(&tagSize));
        
        if (tagSize > MAX_MESSAGE_SIZE) {
            throw runtime_error("[stream::MessageIterator::get_next] tag of " +
                                to_string(tagSize) + " bytes is too long");
        }
        
        // Read it into the tag field of our value
        value.first.clear();
        if (tagSize) {
            handle(coded_in.ReadString(&value.first, tagSize));
        }
        
#ifdef debug
        cerr << "Read what should be the tag of " << tagSize << " bytes" << endl;
#endif
        
        // Update the counters for the tag, which the counters treat as a message.
        if (virtual_offset == -1) {
            // Just track the counter.
            item_vo++;
        } else {
            // We know where here is
            item_vo = virtual_offset;
        }
        
        // Move on to the next message in the group
        group_idx++;
    
        // Work out if this really is a tag.
        bool is_tag = false;
        
        if (!backup_tag.empty() && backup_tag == value.first) {
#ifdef debug
            cerr << "Tag is the same as the last tag of \"" << backup_tag << "\"" << endl;
#endif
            is_tag = true;
        } else {
#ifdef debug
            cerr << "Tag does not match cached backup tag or there is no cached backup tag" << endl;
#endif
        }
    
        if (!is_tag && Registry::is_valid_tag(value.first)) {
#ifdef debug
            cerr << "Tag \"" << value.first << "\" is OK with the registry" << endl;
#endif
            is_tag = true;
        } else if (!is_tag) {
#ifdef debug
            cerr << "Tag is not approved by the registry" << endl;
#endif
        }
    
        if (!is_tag) {
            // If we get here, the registry doesn't think it's a tag.
            // Assume it is actually a message, and make the group's tag ""
            swap(value.first, value.second);
            value.first.clear();
            backup_tag.clear();
            
#ifdef debug
            cerr << "Tag is actually a message probably." << endl;
#endif

#ifdef debug
            cerr << "Found message with tag \"" << value.first << "\"" << endl;
#endif
            
            // Return ourselves, after increment
            return *this;
        }
        
        // Otherwise this is a real tag.
        // Back up its value in case our pair gets moved away.
        backup_tag = value.first;
        
        // We continue and read the real first message into the message half of our pair.
    }
    
    // Now we know we're in a group, and we know the tag, if any.
    
    // Get the item's virtual offset, if available
    auto virtual_offset = bgzip_in->Tell();
    
    // We need a fresh CodedInputStream every time, because of the total byte limit
    ::google::protobuf::io::CodedInputStream coded_in(bgzip_in.get());
    // Alot space for size and item (generously)
    coded_in.SetTotalBytesLimit(MAX_MESSAGE_SIZE * 2, MAX_MESSAGE_SIZE * 2);
    
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
    
    if (msgSize > MAX_MESSAGE_SIZE) {
        throw runtime_error("[stream::MessageIterator::get_next] message of " +
                            to_string(msgSize) + " bytes is too long");
    }
    
    
    // We have a message.
    value.second.clear();
    if (msgSize) {
        handle(coded_in.ReadString(&value.second, msgSize));
    }
    
    // Fill in the tag from the backup to make sure our value pair actually has it.
    // It may have been moved away.
    value.first = backup_tag;
    
    // Move on to the next message in the group
    group_idx++;
    
#ifdef debug
    cerr << "Found message with tag \"" << value.first << "\"" << endl;
#endif
    
    // Return ourselves, after increment
    return *this;
}

auto MessageIterator::operator==(const MessageIterator& other) const -> bool {
    // Just ask if we both agree on whether we hit the end.
    return has_next() == other.has_next();
}
    
auto MessageIterator::operator!=(const MessageIterator& other) const -> bool {
    // Just ask if we disagree on whether we hit the end.
    return has_next() != other.has_next();
}

auto MessageIterator::has_next() const -> bool {
    return item_vo != -1;
}

auto MessageIterator::get_next() -> void {
    // Run increment but don't return anything.
    ++(*this);
}

auto MessageIterator::take() -> TaggedMessage {
    auto temp = std::move(value);
    get_next();
    // Return by value, which gets moved.
    return temp;
}

auto MessageIterator::tell_group() const -> int64_t {
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

auto MessageIterator::seek_group(int64_t virtual_offset) -> bool {
    if (virtual_offset < 0) {
        // That's not allowed
#ifdef debug
        cerr << "Can't seek to negative position" << endl;
#endif
        return false;
    }
    
    if (group_idx == 0 && group_vo == virtual_offset) {
        // We are there already
#ifdef debug
        cerr << "Already at seek position" << endl;
#endif
        return true;
    }
    
    // Try and do the seek
    bool sought = bgzip_in->Seek(virtual_offset);
    
    if (!sought) {
        // We can't seek
#ifdef debug
        cerr << "bgzip_in could not seek" << endl;
#endif
        return false;
    }
    
    // Get ready to read the group that's here
    group_count = 0;
    group_idx = 0;
    
#ifdef debug
    cerr << "Successfully sought" << endl;
#endif
    
    // Read it (or detect EOF)
    get_next();
    
    // It worked!
    return true;
}

auto MessageIterator::range(istream& in) -> pair<MessageIterator, MessageIterator> {
    return make_pair(MessageIterator(in), MessageIterator());
}

auto MessageIterator::handle(bool ok) -> void {
    if (!ok) {
        throw runtime_error("[stream::MessageIterator] obsolete, invalid, or corrupt protobuf input");
    }
}

}

}
