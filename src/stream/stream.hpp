#ifndef VG_STREAM_STREAM_HPP_INCLUDED
#define VG_STREAM_STREAM_HPP_INCLUDED

// de/serialization of protobuf objects from/to a length-prefixed, gzipped binary stream
// from http://www.mail-archive.com/protobuf@googlegroups.com/msg03417.html

#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <functional>
#include <vector>
#include <list>

#include "registry.hpp"
#include "message_iterator.hpp"
#include "protobuf_iterator.hpp"
#include "protobuf_emitter.hpp"

namespace vg {

namespace stream {

using namespace std;

/// Write the EOF marker to the given stream, so that readers won't complain that it might be truncated when they read it in.
/// Internal EOF markers MAY exist, but a file SHOULD have exactly one EOF marker at its end.
void finish(std::ostream& out);

/// Write objects. count should be equal to the number of objects to write.
/// count is written before the objects, but if it is 0, it is not written. To
/// get the objects, calls lambda with the index of the object to retrieve. If
/// not all objects are written, return false, otherwise true.
template <typename T>
bool write(std::ostream& out, size_t count, const std::function<T&(size_t)>& lambda) {

    // Wrap stream in an emitter
    ProtobufEmitter<T> emitter(out);
    
    for (size_t i = 0; i < count; i++) {
        // Write each item.
        emitter.write_copy(lambda(i));
    }
    
    return true;
}

/// Write objects. count should be equal to the number of objects to write.
/// count is written before the objects, but if it is 0, it is not written. To
/// get the objects, calls lambda with the index of the object to retrieve. If
/// not all objects are written, return false, otherwise true.
/// This implementation takes a function that returns actual objects and not references.
template <typename T>
bool write(std::ostream& out, size_t count, const std::function<T(size_t)>& lambda) {

    static_assert(!std::is_reference<T>::value);

    // Wrap stream in an emitter
    ProtobufEmitter<T> emitter(out);
    
    for (size_t i = 0; i < count; i++) {
        // Write each item.
        emitter.write_copy(lambda(i));
    }
    
    return true;
}

/// Start, continue, or finish a buffered stream of objects.
/// If the length of the buffer is greater than the limit, writes the buffer out.
/// Otherwise, leaves the objects in the buffer.
/// Must be called with a buffer limit of 0 after all the objects have been produced, to flush the buffer.
/// When called with a buffer limit of 0, automatically appends an EOF marker.
/// Returns true unless an error occurs.
template <typename T>
bool write_buffered(std::ostream& out, std::vector<T>& buffer, size_t buffer_limit) {
    bool wrote = false;
    if (buffer.size() >= buffer_limit) {
        std::function<T(size_t)> lambda = [&buffer](size_t n) { return buffer.at(n); };
#pragma omp critical (stream_out)
        wrote = write(out, buffer.size(), lambda);
        buffer.clear();
    }
    if (buffer_limit == 0) {
        // The session is over. Append the EOF marker.
        finish(out);
    }
    return wrote;
}

/// Write a single message to a file.
template <typename T>
void write_to_file(const T& item, const string& filename) {
    ofstream out(filename);
    vector<T> items = { item };
    write_buffered(out, items, 1);
    out.close();
}

template <typename T>
void for_each(std::istream& in,
              const std::function<void(int64_t, T&)>& lambda) {
    
    for(ProtobufIterator<T> it(in); it.has_next(); ++it) {
        // For each message in the file, parse and process it with its group VO (or -1)
        lambda(it.tell_group(), *it);
    }
}

template <typename T>
void for_each(std::istream& in,
              const std::function<void(T&)>& lambda) {
    for_each(in, static_cast<const typename std::function<void(int64_t, T&)>&>([&lambda](int64_t virtual_offset, T& item) {
        lambda(item);
    }));
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
                            const std::function<bool(void)>& single_threaded_until_true) {

    // objects will be handed off to worker threads in batches of this many
    // Must be divisible by 2.
    const size_t batch_size = 256;
    static_assert(batch_size % 2 == 0, "stream::for_each_parallel::batch_size must be even");
    // max # of such batches to be holding in memory
    size_t max_batches_outstanding = 256;
    // max # we will ever increase the batch buffer to
    const size_t max_max_batches_outstanding = 1 << 13; // 8192
    // number of batches currently being processed
    size_t batches_outstanding = 0;

    // this loop handles a chunked file with many pieces
    // such as we might write in a multithreaded process
    #pragma omp parallel default(none) shared(in, lambda1, lambda2, batches_outstanding, max_batches_outstanding, single_threaded_until_true)
    #pragma omp single
    {
        auto handle = [](bool retval) -> void {
            if (!retval) throw std::runtime_error("obsolete, invalid, or corrupt protobuf input");
        };
        
        // We do our own multi-threaded Protobuf decoding, but we batch up our strings by pulling them from this iterator.
        MessageIterator message_it(in);

        BlockedGzipInputStream bgzip_in(in);
        ::google::protobuf::io::CodedInputStream coded_in(&bgzip_in);

        std::vector<std::string> *batch = nullptr;
        
        while (message_it.has_next()) {
            // Until we run out of messages, grab them with their tags
            auto tag_and_data = std::move(message_it.take());
            
            // Check the tag.
            // TODO: we should only do this when it changes!
            handle(Registry::check_protobuf_tag<T>(tag_and_data.first));
            
            // If the tag checks out
            
            // Make sure we have a batch
            if (batch == nullptr) {
                batch = new vector<string>();
            }
            
            // Add the message to the batch
            batch->push_back(std::move(tag_and_data.second));
            
            if (batch->size() == batch_size) {
                // time to enqueue this batch for processing. first, block if
                // we've hit max_batches_outstanding.
                size_t b;
#pragma omp atomic capture
                b = ++batches_outstanding;
                
                bool do_single_threaded = !single_threaded_until_true();
                if (b >= max_batches_outstanding || do_single_threaded) {
                    
                    // process this batch in the current thread
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
#pragma omp atomic capture
                    b = --batches_outstanding;
                    
                    if (4 * b / 3 < max_batches_outstanding
                        && max_batches_outstanding < max_max_batches_outstanding
                        && !do_single_threaded) {
                        // we went through at least 1/4 of the batch buffer while we were doing this thread's batch
                        // this looks risky, since we want the batch buffer to stay populated the entire time we're
                        // occupying this thread on compute, so let's increase the batch buffer size
                        // (skip this adjustment if you're in single-threaded mode and thus expect the buffer to be
                        // empty)
                        max_batches_outstanding *= 2;
                    }
                }
                else {
                    // spawn a task in another thread to process this batch
#pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda2, handle, single_threaded_until_true)
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
                }

                batch = nullptr;
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
    std::function<bool(void)> no_wait = [](void) {return true;};
    for_each_parallel_impl(in, lambda2, err1, no_wait);
}
    
template <typename T>
void for_each_interleaved_pair_parallel_after_wait(std::istream& in,
                                                   const std::function<void(T&,T&)>& lambda2,
                                                   const std::function<bool(void)>& single_threaded_until_true) {
    std::function<void(T&)> err1 = [](T&){
        throw std::runtime_error("stream::for_each_interleaved_pair_parallel: expected input stream of interleaved pairs, but it had odd number of elements");
    };
    for_each_parallel_impl(in, lambda2, err1, single_threaded_until_true);
}

// parallelized for each individual element
template <typename T>
void for_each_parallel(std::istream& in,
                       const std::function<void(T&)>& lambda1) {
    std::function<void(T&,T&)> lambda2 = [&lambda1](T& o1, T& o2) { lambda1(o1); lambda1(o2); };
    std::function<bool(void)> no_wait = [](void) {return true;};
    for_each_parallel_impl(in, lambda2, lambda1, no_wait);
}

}

}

#endif
