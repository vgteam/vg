#ifndef STREAM_H
#define STREAM_H

// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <fstream>
#include <functional>
#include <google/protobuf/stubs/common.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/coded_stream.h>

namespace stream {

// write objects
// count should be equal to the number of objects to write
// but if it is 0, it is not written
// if not all objects are written, return false, otherwise true
template <typename T>
bool write(std::ostream& out, uint64_t count, std::function<T(uint64_t)>& lambda) {

    ::google::protobuf::io::ZeroCopyOutputStream *raw_out =
          new ::google::protobuf::io::OstreamOutputStream(&out);
    ::google::protobuf::io::GzipOutputStream *gzip_out =
          new ::google::protobuf::io::GzipOutputStream(raw_out);
    ::google::protobuf::io::CodedOutputStream *coded_out =
          new ::google::protobuf::io::CodedOutputStream(gzip_out);

    coded_out->WriteVarint64(count);

    std::string s;
    uint64_t written = 0;
    for (uint64_t n = 0; n < count; ++n, ++written) {
        lambda(n).SerializeToString(&s);
        coded_out->WriteVarint32(s.size());
        coded_out->WriteRaw(s.data(), s.size());
    }

    delete coded_out;
    delete gzip_out;
    delete raw_out;

    return !count || written == count;
}


// deserialize the input stream into the objects
// count containts the count read
// takes a callback function to be called on the objects

template <typename T>
bool for_each(std::istream& in,
              std::function<void(T&)>& lambda,
              std::function<void(uint64_t)>& handle_count) {

    ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(&in);
    ::google::protobuf::io::GzipInputStream *gzip_in =
          new ::google::protobuf::io::GzipInputStream(raw_in);

    uint64_t count;

    // this loop handles a chunked file with many pieces
    while (in.good()) {

        ::google::protobuf::io::CodedInputStream *coded_in =
            new ::google::protobuf::io::CodedInputStream(gzip_in);

        coded_in->ReadVarint64(&count);
        handle_count(count);
        delete coded_in;

        std::string s;

        for (uint64_t i = 0; i < count; ++i) {

            ::google::protobuf::io::CodedInputStream *coded_in =
                new ::google::protobuf::io::CodedInputStream(gzip_in);

            uint32_t msgSize = 0;
            coded_in->ReadVarint32(&msgSize);

            if ((msgSize > 0) &&
                (coded_in->ReadString(&s, msgSize))) {
                T object;
                object.ParseFromString(s);
                lambda(object);
            }

            delete coded_in;
        }

    }

    delete gzip_in;
    delete raw_in;

    return !count;
}

template <typename T>
bool for_each(std::istream& in,
              std::function<void(T&)>& lambda) {
    std::function<void(uint64_t)> noop = [](uint64_t) { };
    for_each(in, lambda, noop);
}

}

#endif
