#ifndef VG_STREAM_HPP_INCLUDED
#define VG_STREAM_HPP_INCLUDED

// de/serialization of protobuf objects from/to a length-prefixed, gzipped binary stream
// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <istream>
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

/// Protobuf will refuse to read messages longer than this size.
const size_t MAX_PROTOBUF_SIZE = 67108864;
/// We aim to generate messages that are this size
const size_t TARGET_PROTOBUF_SIZE = MAX_PROTOBUF_SIZE/2;

/// Write objects using adaptive chunking. Takes a stream to write to, a total
/// element count to write, a guess at how manye elements should be in a chunk,
/// and a function that, given a start element and a length, returns a Protobuf
/// object representing that range of elements.
///
/// Adaptively sets the chunk size, in elements, so that no too-large Protobuf
/// records are serialized.
template <typename T>
bool write(std::ostream& out, uint64_t element_count, uint64_t chunk_elements,
    const std::function<T(uint64_t, uint64_t)>& lambda) {

    // How many elements have we serialized so far
    size_t serialized = 0;
    
    ::google::protobuf::io::OstreamOutputStream raw_out(&out);
    ::google::protobuf::io::GzipOutputStream gzip_out(&raw_out);
    ::google::protobuf::io::CodedOutputStream coded_out(&gzip_out);

    auto handle = [](bool ok) {
        if (!ok) throw std::runtime_error("stream::write: I/O error writing protobuf");
    };
    
    while (serialized < element_count) {
    
        // Work out how many elements can go in this chunk, accounting for the total element count
        chunk_elements = std::min(chunk_elements, element_count - serialized);
    
        // Serialize a chunk
        std::string chunk_data;
        handle(lambda(serialized, chunk_elements).SerializeToString(&chunk_data));
    
        if (chunk_data.size() > MAX_PROTOBUF_SIZE) {
            // This is too big!
            
            if (chunk_elements > 1) {
                // But we can make it smaller. Try again at half this size.
                chunk_elements = chunk_elements / 2;
                continue;
            } else {
                // This single element is too large
                throw std::runtime_error("stream::write: message for element " +
                    std::to_string(serialized) + " too large error writing protobuf");
            }
        } else {
            // We can send this message
            
            // Say we have a group of a single message
            coded_out.WriteVarint64(1);
            handle(!coded_out.HadError());
            // and prefix each object with its size
            coded_out.WriteVarint32(chunk_data.size());
            handle(!coded_out.HadError());
            coded_out.WriteRaw(chunk_data.data(), chunk_data.size());
            handle(!coded_out.HadError());
            
            // Remember how far we've serialized now
            serialized += chunk_elements;
            
            if (chunk_data.size() < TARGET_PROTOBUF_SIZE/2) {
                // We were less than half the target size, so try being twice as
                // big next time.
                chunk_elements *= 2;
            } else if (chunk_data.size() > TARGET_PROTOBUF_SIZE && chunk_elements > 1) {
                // We were larger than the target size and we can be smaller
                chunk_elements /= 2;
            }
        }
    }
    
    

}

// write objects
// count should be equal to the number of objects to write
// count is written before the objects, but if it is 0, it is not written
// if not all objects are written, return false, otherwise true
template <typename T>
bool write(std::ostream& out, uint64_t count, const std::function<T(uint64_t)>& lambda) {

    // Make all our streams on the stack, in case of error.
    ::google::protobuf::io::OstreamOutputStream raw_out(&out);
    ::google::protobuf::io::GzipOutputStream gzip_out(&raw_out);
    ::google::protobuf::io::CodedOutputStream coded_out(&gzip_out);

    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("stream::write: I/O error writing protobuf");
        }
    };

    // prefix the chunk with the number of objects, if any objects are to be written
    if(count > 0) {
        coded_out.WriteVarint64(count);
        handle(!coded_out.HadError());
    }

    std::string s;
    uint64_t written = 0;
    for (uint64_t n = 0; n < count; ++n, ++written) {
        handle(lambda(n).SerializeToString(&s));
        if (s.size() > MAX_PROTOBUF_SIZE) {
            throw std::runtime_error("stream::write: message too large error writing protobuf");
        }
        // and prefix each object with its size
        coded_out.WriteVarint32(s.size());
        handle(!coded_out.HadError());
        coded_out.WriteRaw(s.data(), s.size());
        handle(!coded_out.HadError());
    }

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

    ::google::protobuf::io::IstreamInputStream raw_in(&in);
    ::google::protobuf::io::GzipInputStream gzip_in(&raw_in);
    ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);

    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("[stream::for_each] obsolete, invalid, or corrupt protobuf input");
        }
    };

    uint64_t count;
    // this loop handles a chunked file with many pieces
    // such as we might write in a multithreaded process
    while (coded_in.ReadVarint64((::google::protobuf::uint64*) &count)) {

        handle_count(count);

        std::string s;
        for (uint64_t i = 0; i < count; ++i) {
            uint32_t msgSize = 0;
            // Reconstruct the CodedInputStream in place to reset its maximum-
            // bytes-ever-read counter, because it thinks it's reading a single
            // message.
            coded_in.~CodedInputStream();
            new (&coded_in) ::google::protobuf::io::CodedInputStream(&gzip_in);
            // Alot space for size, and for reading next chunk's length
            coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
            
            // the messages are prefixed by their size
            handle(coded_in.ReadVarint32(&msgSize));
            
            if (msgSize > MAX_PROTOBUF_SIZE) {
                throw std::runtime_error("[stream::for_each] protobuf message of " +
                    std::to_string(msgSize) + " bytes is too long");
            }
            
            if (msgSize) {
                handle(coded_in.ReadString(&s, msgSize));
                T object;
                handle(object.ParseFromString(s));
                lambda(object);
            }
        }
    }
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
void for_each_parallel_impl(std::istream& in,
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

        ::google::protobuf::io::IstreamInputStream raw_in(&in);
        ::google::protobuf::io::GzipInputStream gzip_in(&raw_in);
        ::google::protobuf::io::CodedInputStream coded_in(&gzip_in);

        std::vector<std::string> *batch = nullptr;

        // process chunks prefixed by message count
        uint64_t count;
        while (coded_in.ReadVarint64((::google::protobuf::uint64*) &count)) {
            handle_count(count);
            for (uint64_t i = 0; i < count; ++i) {
                if (!batch) {
                     batch = new std::vector<std::string>();
                     batch->reserve(batch_size);
                }
                
                // Reconstruct the CodedInputStream in place to reset its maximum-
                // bytes-ever-read counter, because it thinks it's reading a single
                // message.
                coded_in.~CodedInputStream();
                new (&coded_in) ::google::protobuf::io::CodedInputStream(&gzip_in);
                // Alot space for size, and for reading next chunk's length
                coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
                
                uint32_t msgSize = 0;
                // the messages are prefixed by their size
                handle(coded_in.ReadVarint32(&msgSize));
                
                if (msgSize > MAX_PROTOBUF_SIZE) {
                    throw std::runtime_error("[stream::for_each] protobuf message of " +
                        std::to_string(msgSize) + " bytes is too long");
                }
                
                if (msgSize) {
                    // pick off the message (serialized protobuf object)
                    std::string s;
                    handle(coded_in.ReadString(&s, msgSize));
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
    }
}

// parallel iteration over interleaved pairs of elements; error out if there's an odd number of elements
template <typename T>
void for_each_interleaved_pair_parallel(std::istream& in,
                                        const std::function<void(T&,T&)>& lambda2) {
    std::function<void(T&)> err1 = [](T&){
        throw std::runtime_error("stream::for_each_interleaved_pair_parallel: expected input stream of interleaved pairs, but it had odd number of elements");
    };
    for_each_parallel_impl(in, lambda2, err1, [](uint64_t) { });
}

// parallelized for each individual element
template <typename T>
void for_each_parallel(std::istream& in,
                       const std::function<void(T&)>& lambda1,
                       const std::function<void(uint64_t)>& handle_count) {
    std::function<void(T&,T&)> lambda2 = [&lambda1](T& o1, T& o2) { lambda1(o1); lambda1(o2); };
    for_each_parallel_impl(in, lambda2, lambda1, handle_count);
}

template <typename T>
void for_each_parallel(std::istream& in,
              const std::function<void(T&)>& lambda) {
    std::function<void(uint64_t)> noop = [](uint64_t) { };
    for_each_parallel(in, lambda, noop);
}

    
/*
 * Refactored stream::for_each function that follows the unidirectional iterator interface
 */
template <typename T>
class ProtobufIterator {
public:
    /// Constructor
    ProtobufIterator(std::istream& in) :
        where(0),
        chunk_count(0),
        chunk_idx(0),
        raw_in(&in),
        gzip_in(&raw_in),
        coded_in(&gzip_in)
    {
        get_next();
    }
    
    
    
//    inline ProtobufIterator<T>& operator=(const ProtobufIterator<T>& other) {
//        value = other.value;
//        where = other.where;
//        chunk_count = other.chunk_count;
//        chunk_idx = other.chunk_idx;
//        raw_in = other.raw_in;
//        gzip_in = other.gzip_in;
//        coded_in = other.coded_in;
//    }

    
    inline bool has_next() {
        return where != 0;
    }
    
    void get_next() {
        if (chunk_count == chunk_idx) {
            chunk_idx = 0;
            if (!coded_in.ReadVarint64((::google::protobuf::uint64*) &chunk_count)) {
                // This is the end of the input stream, switch to state that
                // will match the end constructor
                where = 0;
                value = T();
                return;
            }
        }
        
        std::string s;
        uint32_t msgSize = 0;
        // Reconstruct the CodedInputStream in place to reset its maximum-
        // bytes-ever-read counter, because it thinks it's reading a single
        // message.
        coded_in.~CodedInputStream();
        new (&coded_in) ::google::protobuf::io::CodedInputStream(&gzip_in);
        // Alot space for size, and for reading next chunk's length
        coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
        
        // the messages are prefixed by their size
        handle(coded_in.ReadVarint32(&msgSize));
        
        if (msgSize > MAX_PROTOBUF_SIZE) {
            throw std::runtime_error("[stream::ProtobufIterator::get_next] protobuf message of " +
                                     std::to_string(msgSize) + " bytes is too long");
        }
        
        if (msgSize) {
            value.Clear();
            handle(coded_in.ReadString(&s, msgSize));
            handle(value.ParseFromString(s));
        }
        else {
            get_next();
        }
        
        where++;
        chunk_idx++;
    }
    
//    inline ProtobufIterator<T> operator++( int ) {
//        ProtobufIterator<T> temp = *this;
//        get_next();
//        return temp;
//    }
    
    inline T operator*(){
        return value;
    }
    
private:
    
    T value;
    
    // For testing identity with another iterator
    uint64_t where;
    
    uint64_t chunk_count;
    uint64_t chunk_idx;
    
    ::google::protobuf::io::IstreamInputStream raw_in;
    ::google::protobuf::io::GzipInputStream gzip_in;
    ::google::protobuf::io::CodedInputStream coded_in;
    
    void handle(bool ok) {
        if (!ok) {
            throw std::runtime_error("[stream::for_each] obsolete, invalid, or corrupt protobuf input");
        }
    }
};

/// Returns iterators that act like begin() and end() for a stream containing protobuf data
template <typename T>
std::pair<ProtobufIterator<T>, ProtobufIterator<T>> protobuf_iterator_range(std::istream& in) {
    return std::make_pair(ProtobufIterator<T>(in), ProtobufIterator<T>());
}

}

#endif
