#ifndef STREAM_H
#define STREAM_H

// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <fstream>
#include <google/protobuf/stubs/common.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>

namespace stream {

// serialize a container into the output stream
// count should be equal to the number of objects to write
// but if it is 0, it is not written
// if not all objects are written, return false, otherwise true
template <typename C>
bool save_messages(ostream& out, C& objects, uint64_t count) {

    ::google::protobuf::io::ZeroCopyOutputStream *raw_out =
          new ::google::protobuf::io::OstreamOutputStream(&out);
    ::google::protobuf::io::CodedOutputStream *coded_out =
          new ::google::protobuf::io::CodedOutputStream(raw_out);

    // save the number of the messages to be serialized into the output file
    if (count) {
        coded_out->WriteVarint64(count);
    }

    std::string s;

    uint64_t written = 0;
    for (typename C::iterator o = objects.begin(); o != objects.end(); ++o) {
        o->SerializeToString(&s);
        coded_out->WriteVarint32(s.size());
        coded_out->WriteRaw(s.data(), s.size()); // ->WriteString(s)
        ++written;
    }

    delete coded_out;
    delete raw_out;

    return !count || written == count;
}

// deserialize the input stream into the objects
// count containts the count read
// takes a callback function to be called on the objects
/*
template <typename T>
struct Callback {
    void operator()(T t) { return; }
};
*/

template <typename T, typename C>
bool load_messages(istream& in, C* object, uint64_t& count) {

    ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(&in);
    ::google::protobuf::io::CodedInputStream *coded_in =
          new ::google::protobuf::io::CodedInputStream(raw_in);

    coded_in->ReadVarint64(&count);
    coded_in->SetTotalBytesLimit(1000000000, -1);

    std::string s;

    for (uint64_t i = 0; i < count; ++i) {

        uint32_t msgSize = 0;
        coded_in->ReadVarint32(&msgSize);

        if ((msgSize > 0) &&
            (coded_in->ReadString(&s, msgSize))) {

            T o;
            o.ParseFromString(s);
            // TODO
            // it would be so nice to have a callback here
            // how can this happen?
            //callback(o);
            object->extend(o);

        }
    }

    delete coded_in;
    delete raw_in;

    return true;
}

}

#endif
