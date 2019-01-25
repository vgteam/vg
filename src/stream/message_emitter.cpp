/**
 * \file message_emitter.cpp
 * Implementations for the MessageEmitter.
 */

#include "message_emitter.hpp"

namespace vg {

namespace stream {

using namespace std;

MessageEmitter::MessageEmitter(ostream& out, size_t max_group_size) :
    group(),
    max_group_size(max_group_size),
    bgzip_out(new BlockedGzipOutputStream(out))
{
#ifdef debug
    cerr << "Creating MessageEmitter" << endl;
#endif
    if (bgzip_out->Tell() == -1) {
        // Say we are starting at the beginnign of the stream, if we don't know where we are.
        bgzip_out->StartFile();
    }
}

MessageEmitter::~MessageEmitter() {
#ifdef debug
    cerr << "Destroying MessageEmitter" << endl;
#endif
    if (bgzip_out.get() != nullptr) {
#ifdef debug
        cerr << "MessageEmitter emitting final group" << endl;
#endif
    
        // Before we are destroyed, write stuff out.
        emit_group();
        
#ifdef debug
        cerr << "MessageEmitter ending file" << endl;
#endif
        
        // Tell our stream to finish the file (since it hasn't been moved away)
        bgzip_out->EndFile();
    }
    
#ifdef debug
    cerr << "MessageEmitter destroyed" << endl;
#endif
}

void MessageEmitter::write(const string& tag, string&& message) {
    if (group.size() >= max_group_size || tag != group_tag) {
        // We have run out of buffer space or changed type
        emit_group();
    }
    if (tag != group_tag) {
        // Adopt the new tag
        group_tag = tag;
    }
    group.emplace_back(std::move(message));
    
    if (group.back().size() > MAX_MESSAGE_SIZE) {
        throw std::runtime_error("stream::MessageEmitter::write: message too large");
    }
}

void MessageEmitter::write_copy(const string& tag, const string& message) {
    if (group.size() >= max_group_size || tag != group_tag) {
        // We have run out of buffer space or changed type
        emit_group();
    }
    if (tag != group_tag) {
        // Adopt the new tag
        group_tag = tag;
    }
    group.push_back(message);
    
    if (group.back().size() > MAX_MESSAGE_SIZE) {
        throw std::runtime_error("stream::MessageEmitter::write_copy: message too large");
    }
}

void MessageEmitter::on_group(group_listener_t&& listener) {
    group_handlers.emplace_back(std::move(listener));
}

void MessageEmitter::emit_group() {
    if (group.empty()) {
        // Nothing to do
        return;
    }
    
    // We can't write a non-empty buffer if our stream is gone/moved away
    assert(bgzip_out.get() != nullptr);
    
    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("stream::MessageEmitter::emit_group: I/O error writing protobuf");
        }
    };

    // Work out where the group we emit will start
    int64_t virtual_offset = bgzip_out->Tell();

    ::google::protobuf::io::CodedOutputStream coded_out(bgzip_out.get());

    // Prefix the group with the number of objects, plus 1 for the tag header
    coded_out.WriteVarint64(group.size() + 1);
    handle(!coded_out.HadError());
   
    // Write the tag length and the tag
    coded_out.WriteVarint32(group_tag.size());
    handle(!coded_out.HadError());
    coded_out.WriteRaw(group_tag.data(), group_tag.size());
    handle(!coded_out.HadError());

    for (auto& message : group) {
        
#ifdef debug
        cerr << "Writing message of " << message.size() << " bytes in group of \""
            << group_tag << "\" @ " << virtual_offset << endl;
#endif
        
        // And prefix each object with its size
        coded_out.WriteVarint32(message.size());
        handle(!coded_out.HadError());
        coded_out.WriteString(message);
        handle(!coded_out.HadError());
    }
    
    // Work out where we ended
    coded_out.Trim();
    int64_t next_virtual_offset = bgzip_out->Tell();
    
    for (auto& handler : group_handlers) {
        // Report the group to each group handler that is listening.
        // We don't report the individual messages. They need to be observed separately.
        handler(group_tag, virtual_offset, next_virtual_offset);
    }
    
    // Empty the buffer because everything in it is written
    group.clear();
}

}

}
