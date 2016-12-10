#ifndef STREAM_H
#define STREAM_H

// de/serialization of protobuf objects from/to a length-prefixed, gzipped binary stream
// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <list>
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

// Parallelized versions of for_each

// First, an internal implementation underlying several variants below.
// lambda2 is invoked on interleaved pairs of elements from the stream. The
// elements of each pair are in order, but the overall order in which lambda2
// is invoked on pairs is undefined (concurrent). lambda1 is invoked on an odd
// last element of the stream, if any.
template <typename T>
void __for_each_parallel_impl(std::istream& in,
                              const std::function<void(T&,T&)>& lambda2,
                              const std::function<void(T&)>& lambda1,
                              const std::function<void(uint64_t)>& handle_count) {

    // objects will be handed off to worker threads in batches of this many
    const uint64_t batch_size = 256;
    static_assert(batch_size % 2 == 0, "stream::for_each_parallel::batch_size must be even");
    // max # of such batches to be holding in memory
    const uint64_t max_batches_outstanding = 256;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;

    // this loop handles a chunked file with many pieces
    // such as we might write in a multithreaded process
    #pragma omp parallel default(none) shared(in, lambda1, lambda2, handle_count, batches_outstanding)
    #pragma omp single
    {
        auto handle = [](bool retval) -> void {
            if (!retval) throw std::runtime_error("obsolete, invalid, or corrupt protobuf input");
        };

        ::google::protobuf::io::ZeroCopyInputStream *raw_in =
            new ::google::protobuf::io::IstreamInputStream(&in);
        ::google::protobuf::io::GzipInputStream *gzip_in =
            new ::google::protobuf::io::GzipInputStream(raw_in);
        ::google::protobuf::io::CodedInputStream *coded_in =
            new ::google::protobuf::io::CodedInputStream(gzip_in);

        std::vector<std::string> *batch = nullptr;

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

                if (batch->size() == batch_size) {
                    // time to enqueue this batch for processing. first, block if
                    // we've hit max_batches_outstanding.
                    uint64_t b;
                    #pragma omp atomic capture
                    b = ++batches_outstanding;
                    while (b >= max_batches_outstanding) {
                        usleep(1000);
                        #pragma omp atomic read
                        b = batches_outstanding;
                    }
                    // spawn task to process this batch
                    #pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda2, handle)
                    {
                        {
                            T obj1, obj2;
                            for (int i = 0; i<batch_size; i+=2) {
                                // parse protobuf objects and invoke lambda on the pair
                                handle(obj1.ParseFromString(batch->at(i)));
                                handle(obj2.ParseFromString(batch->at(i+1)));
                                lambda2(obj1,obj2);
                            }
                        } // scope obj1 & obj2
                        delete batch;
                        #pragma omp atomic update
                        batches_outstanding--;
                    }

                    batch = nullptr;
                }

                // recycle the CodedInputStream in order to avoid its cumulative byte limit
                delete coded_in;
                coded_in = new ::google::protobuf::io::CodedInputStream(gzip_in);
            }
        }

        #pragma omp taskwait
        // process final batch
        if (batch) {
            {
                T obj1, obj2;
                int i = 0;
                for (; i < batch->size()-1; i+=2) {
                    handle(obj1.ParseFromString(batch->at(i)));
                    handle(obj2.ParseFromString(batch->at(i+1)));
                    lambda2(obj1, obj2);
                }
                if (i == batch->size()-1) { // odd last object
                    handle(obj1.ParseFromString(batch->at(i)));
                    lambda1(obj1);
                }
            } // scope obj1 & obj2
            delete batch;
        }

        delete coded_in;
        delete gzip_in;
        delete raw_in;
    }
}

// parallel iteration over interleaved pairs of elements; error out if there's an odd number of elements
template <typename T>
void for_each_interleaved_pair_parallel(std::istream& in,
                                        const std::function<void(T&,T&)>& lambda2) {
    std::function<void(T&)> err1 = [](T&){
        throw std::runtime_error("stream::for_each_interleaved_pair_parallel: expected input stream of interleaved pairs, but it had odd number of elements");
    };
    __for_each_parallel_impl(in, lambda2, err1, [](uint64_t) { });
}

// parallelized for each individual element
template <typename T>
void for_each_parallel(std::istream& in,
                       const std::function<void(T&)>& lambda1,
                       const std::function<void(uint64_t)>& handle_count) {
    std::function<void(T&,T&)> lambda2 = [&lambda1](T& o1, T& o2) { lambda1(o1); lambda1(o2); };
    __for_each_parallel_impl(in, lambda2, lambda1, handle_count);
}

template <typename T>
void for_each_parallel(std::istream& in,
              const std::function<void(T&)>& lambda) {
    std::function<void(uint64_t)> noop = [](uint64_t) { };
    for_each_parallel(in, lambda, noop);
}

}

#endif
