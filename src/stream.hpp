#ifndef STREAM_H
#define STREAM_H

// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <list>
#include <atomic>
#include "google/protobuf/stubs/common.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/io/gzip_stream.h"
#include "google/protobuf/io/coded_stream.h"

namespace stream {

// write objects
// count should be equal to the number of objects to write
// count is written before the objects, but if it is 0, it is not written
// if not all objects are written, return false, otherwise true
template <typename T>
bool write(std::ostream& out, uint64_t count, const std::function<T(uint64_t)>& lambda) {

    ::google::protobuf::io::ZeroCopyOutputStream *raw_out =
          new ::google::protobuf::io::OstreamOutputStream(&out);
    ::google::protobuf::io::GzipOutputStream *gzip_out =
          new ::google::protobuf::io::GzipOutputStream(raw_out);
    ::google::protobuf::io::CodedOutputStream *coded_out =
          new ::google::protobuf::io::CodedOutputStream(gzip_out);

    auto handle = [](bool ok) {
        if (!ok) throw std::runtime_error("stream::write: I/O error writing protobuf");
    };

    // prefix the chunk with the number of objects, if any objects are to be written
    if(count > 0) {
        coded_out->WriteVarint64(count);
        handle(!coded_out->HadError());
    }

    std::string s;
    uint64_t written = 0;
    for (uint64_t n = 0; n < count; ++n, ++written) {
        handle(lambda(n).SerializeToString(&s));
        // and prefix each object with its size
        coded_out->WriteVarint32(s.size());
        handle(!coded_out->HadError());
        coded_out->WriteRaw(s.data(), s.size());
        handle(!coded_out->HadError());
    }

    delete coded_out;
    delete gzip_out;
    delete raw_out;

    return !count || written == count;
}

template <typename T>
bool write_buffered(std::ostream& out, std::vector<T>& buffer, uint64_t buffer_limit) {
    bool wrote = false;
    if (buffer.size() >= buffer_limit) {
        std::function<T(uint64_t)> lambda = [&buffer](uint64_t n) { return buffer.at(n); };
#pragma omp critical (stream_out)
        wrote = write(out, buffer.size(), lambda);
        buffer.clear();
    }
    return wrote;
}

// deserialize the input stream into the objects
// skips over groups of objects with count 0
// takes a callback function to be called on the objects, and another to be called per object group.

template <typename T>
void for_each(std::istream& in,
              const std::function<void(T&)>& lambda,
              const std::function<void(uint64_t)>& handle_count) {

    ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(&in);
    ::google::protobuf::io::GzipInputStream *gzip_in =
          new ::google::protobuf::io::GzipInputStream(raw_in);
    ::google::protobuf::io::CodedInputStream *coded_in =
          new ::google::protobuf::io::CodedInputStream(gzip_in);

    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("[stream::for_each] obsolete, invalid, or corrupt protobuf input");
        }
    };

    uint64_t count;
    // this loop handles a chunked file with many pieces
    // such as we might write in a multithreaded process
    while (coded_in->ReadVarint64((::google::protobuf::uint64*) &count)) {

        handle_count(count);

        std::string s;
        for (uint64_t i = 0; i < count; ++i) {
            uint32_t msgSize = 0;
            delete coded_in;
            coded_in = new ::google::protobuf::io::CodedInputStream(gzip_in);
            // the messages are prefixed by their size
            handle(coded_in->ReadVarint32(&msgSize));
            if (msgSize) {
                handle(coded_in->ReadString(&s, msgSize));
                T object;
                handle(object.ParseFromString(s));
                lambda(object);
            }
        }
    }

    delete coded_in;
    delete gzip_in;
    delete raw_in;
}

template <typename T>
void for_each(std::istream& in,
              const std::function<void(T&)>& lambda) {
    std::function<void(uint64_t)> noop = [](uint64_t) { };
    for_each(in, lambda, noop);
}

template <typename T>
void for_each_parallel(std::istream& in,
                       const std::function<void(T&)>& lambda,
                       const std::function<void(uint64_t)>& handle_count) {

    ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(&in);
    ::google::protobuf::io::GzipInputStream *gzip_in =
          new ::google::protobuf::io::GzipInputStream(raw_in);
    ::google::protobuf::io::CodedInputStream *coded_in =
          new ::google::protobuf::io::CodedInputStream(gzip_in);

    // objects will be handed off to worker threads in batches of this many
    const uint64_t batch_size = 1024;
    // max # of such batches to be holding in memory
    const uint64_t max_batches_outstanding = 256;

    auto handle = [](bool retval) -> void {
        if (!retval) throw std::runtime_error("obsolete, invalid, or corrupt protobuf input");
    };

    // this loop handles a chunked file with many pieces
    // such as we might write in a multithreaded process
    #pragma omp parallel
    #pragma omp single nowait
    {
        std::vector<std::string> *batch = nullptr;
        std::atomic<uint64_t> batches_outstanding;
        batches_outstanding = 0;

        // process chunks prefixed by message count
        uint64_t count;
        while (coded_in->ReadVarint64((::google::protobuf::uint64*) &count)) {
            handle_count(count);
            for (uint64_t i = 0; i < count; ++i) {
                if (!batch) {
                     batch = new std::vector<std::string>();
                     batch->reserve(batch_size);
                }
                uint32_t msgSize = 0;
                // the messages are prefixed by their size
                handle(coded_in->ReadVarint32(&msgSize));
                if (msgSize) {
                    // pick off the message (serialized protobuf object)
                    std::string s;
                    handle(coded_in->ReadString(&s, msgSize));
                    batch->push_back(std::move(s));
                }

                if (batch->size() >= batch_size) {
                    // time to enqueue this batch for processing. first, block if
                    // we've hit max_batches_outstanding.
                    batches_outstanding++;
                    while (batches_outstanding >= max_batches_outstanding) {
                        usleep(1000);
                    }
                    // spawn task to process this batch
                    #pragma omp task firstprivate(batch) shared(batches_outstanding)
                    {
                        {
                            T object;
                            for (const std::string& s_j : *batch) {
                                // parse protobuf object and invoke lambda on it
                                handle(object.ParseFromString(s_j));
                                lambda(object);
                            }
                        } // scope object
                        delete batch;
                        batches_outstanding--;
                    }

                    batch = nullptr;
                }

                // recycle the CodedInputStream in order to avoid its byte limit
                delete coded_in;
                coded_in = new ::google::protobuf::io::CodedInputStream(gzip_in);
            }
        }

        // process final batch
        if (batch) {
            {
                T object;
                for (const std::string& s_j : *batch) {
                    handle(object.ParseFromString(s_j));
                    lambda(object);
                }
            } // scope object
            delete batch;
        }
    }

    delete coded_in;
    delete gzip_in;
    delete raw_in;
}

template <typename T>
void for_each_parallel(std::istream& in,
              const std::function<void(T&)>& lambda) {
    std::function<void(uint64_t)> noop = [](uint64_t) { };
    for_each_parallel(in, lambda, noop);
}

}

#endif
