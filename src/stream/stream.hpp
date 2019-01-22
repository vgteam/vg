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

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/coded_stream.h>

#include "blocked_gzip_output_stream.hpp"
#include "blocked_gzip_input_stream.hpp"

namespace vg {

namespace stream {

using namespace std;

/// Protobuf will refuse to read messages longer than this size.
const size_t MAX_PROTOBUF_SIZE = 1000000000;
/// We aim to generate messages that are this size
const size_t TARGET_PROTOBUF_SIZE = MAX_PROTOBUF_SIZE/2;

/// Write the EOF marker to the given stream, so that readers won't complain that it might be truncated when they read it in.
/// Internal EOF markers MAY exist, but a file SHOULD have exactly one EOF marker at its end.
void finish(std::ostream& out);

/// Write objects using adaptive chunking. Takes a stream to write to, a total
/// element count to write, a guess at how many elements should be in a chunk,
/// and a function that, given a destination virtual offset in the output
/// stream (or -1), a start element, and a length, returns a Protobuf object
/// representing that range of elements.
///
/// Adaptively sets the chunk size, in elements, so that no too-large Protobuf
/// records are serialized.
///
/// Returns true on success, but throws errors on failure.
template <typename T>
bool write(std::ostream& out, size_t element_count, size_t chunk_elements,
    const std::function<T(int64_t, size_t, size_t)>& lambda) {

    // How many elements have we serialized so far
    size_t serialized = 0;
    
    BlockedGzipOutputStream bgzip_out(out);
    ::google::protobuf::io::CodedOutputStream coded_out(&bgzip_out);

    auto handle = [](bool ok) {
        if (!ok) throw std::runtime_error("stream::write: I/O error writing protobuf");
    };
    
    while (serialized < element_count) {
    
        // Work out how many elements can go in this chunk, accounting for the total element count
        chunk_elements = std::min(chunk_elements, element_count - serialized);
    
        // Work out where the chunk is going.
        // TODO: we need to back up the coded output stream after every chunk,
        // and push the partial buffer into BGZF, and get a new buffer, which 
        // wastes time.
#ifdef debug
        cerr << "Trim stream and determine offset" << endl;
#endif
        coded_out.Trim(); 
        int64_t virtual_offset = bgzip_out.Tell();
#ifdef debug
        cerr << "Offset is " << virtual_offset << endl;
#endif
    
        // Serialize a chunk
#ifdef debug
        cerr << "Go get " << chunk_elements << " elements" << endl;
#endif
        std::string chunk_data;
        handle(lambda(virtual_offset, serialized, chunk_elements).SerializeToString(&chunk_data));
    
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
#ifdef debug
            cerr << "Writing message/group of " << chunk_data.size() << " bytes and elements "
                << serialized << "-" << (serialized + chunk_elements) << endl;
#endif
            
            // Say we have a group of a single message
#ifdef debug
            cerr << "\tWrite group length" << endl;
#endif
            coded_out.WriteVarint64(1);
            handle(!coded_out.HadError());
            // and prefix each object with its size
#ifdef debug
            cerr << "\tWrite message length" << endl;
#endif
            coded_out.WriteVarint32(chunk_data.size());
            handle(!coded_out.HadError());
#ifdef debug
            cerr << "\tWrite message data" << endl;
#endif
            coded_out.WriteRaw(chunk_data.data(), chunk_data.size());
            handle(!coded_out.HadError());
#ifdef debug
            cerr << "\tMessage/group written" << endl;
#endif
            
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
    
    return true;

}

/// Write objects using adaptive chunking. Takes a stream to write to, a total
/// element count to write, a guess at how many elements should be in a chunk,
/// and a function that, given a start element and a length, returns a Protobuf
/// object representing that range of elements.
///
/// Adaptively sets the chunk size, in elements, so that no too-large Protobuf
/// records are serialized.
///
/// Returns true on success, but throws errors on failure.
template <typename T>
bool write(std::ostream& out, size_t element_count, size_t chunk_elements,
    const std::function<T(size_t, size_t)>& lambda) {
    
    return write(out, element_count, chunk_elements,
        static_cast<const typename std::function<T(int64_t, size_t, size_t)>&>(
        [&lambda](int64_t virtual_offset, size_t chunk_start, size_t chunk_length) -> T {
        
        // Ignore the virtual offset
        return lambda(chunk_start, chunk_length);
    }));
}

/// Write objects. count should be equal to the number of objects to write.
/// count is written before the objects, but if it is 0, it is not written. To
/// get the objects, calls lambda with the highest virtual offset that can be
/// seek'd to in order to read the object (or -1 if the stream is not
/// tellable), and the index of the object to retrieve. If not all objects are
/// written, return false, otherwise true.
template <typename T>
bool write(std::ostream& out, size_t count, const std::function<T(int64_t, size_t)>& lambda) {

    // Make all our streams on the stack, in case of error.
    BlockedGzipOutputStream bgzip_out(out);
    ::google::protobuf::io::CodedOutputStream coded_out(&bgzip_out);

    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("stream::write: I/O error writing protobuf");
        }
    };
    
    // We can't seek directly to individual messages, because we can only read
    // count-prefixed groups. So the highest seek offset is going to be where
    // we are now, where the group count is being written.
    coded_out.Trim();
    int64_t virtual_offset = bgzip_out.Tell();

    // prefix the chunk with the number of objects, if any objects are to be written
    if(count > 0) {
        coded_out.WriteVarint64(count);
        handle(!coded_out.HadError());
    }

    std::string s;
    size_t written = 0;
    for (size_t n = 0; n < count; ++n, ++written) {
        handle(lambda(virtual_offset, n).SerializeToString(&s));
        if (s.size() > MAX_PROTOBUF_SIZE) {
            throw std::runtime_error("stream::write: message too large error writing protobuf");
        }
        
#ifdef debug
        cerr << "Writing message of " << s.size() << " bytes at " << n << "/" << count << " in group @ " << virtual_offset << endl;
#endif
        
        // and prefix each object with its size
        coded_out.WriteVarint32(s.size());
        handle(!coded_out.HadError());
        coded_out.WriteRaw(s.data(), s.size());
        handle(!coded_out.HadError());
    }

    return !count || written == count;
}

/// Write objects. count should be equal to the number of objects to write.
/// count is written before the objects, but if it is 0, it is not written. To
/// get the objects, calls lambda with the index of the object to retrieve. If
/// not all objects are written, return false, otherwise true.
template <typename T>
bool write(std::ostream& out, size_t count, const std::function<T(size_t)>& lambda) {
    return write(out, count,
        static_cast<const typename std::function<T(int64_t, size_t)>&>(
        [&lambda](int64_t virtual_offset, size_t object_number) -> T {
        // Discard the virtual offset
        return lambda(object_number);
    }));
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

/// Deserialize the input stream into the objects. Skips over groups of objects
/// with count 0. Takes a callback function to be called on the objects, with
/// the object and the blocked gzip virtual offset of its group (or -1 if the
/// input is not blocked gzipped), and another to be called per object group
/// with the group size.
template <typename T>
void for_each_with_group_length(std::istream& in,
              const std::function<void(int64_t, T&)>& lambda,
              const std::function<void(size_t)>& handle_count) {

    BlockedGzipInputStream bgzip_in(in);

    // Have a function to complain if any protobuf things report failure
    auto handle = [](bool ok) {
        if (!ok) {
            throw std::runtime_error("[stream::for_each] obsolete, invalid, or corrupt protobuf input");
        }
    };

    while (true) {
        // For each count-prefixed group
    
        // Get the offset we're at, or -1 if we can't seek/tell
        int64_t virtual_offset = bgzip_in.Tell();
        
#ifdef debug
        cerr << "At virtual offset " << virtual_offset << endl;
#endif
        
        // Read the count
        size_t count;
        {
            
#ifdef debug
            // Peek ahead at the group header
            char* data;
            int size;
            bool worked = bgzip_in.Next((const void**)&data, &size);
            if (worked) {
                cerr << "Next data is " << size << " bytes: " << endl;
                for (size_t i = 0; size > 0 && i < std::min(size, 10); i++) {
                    cerr << (unsigned int)data[i] << " ";
                }
                cerr << endl;
                bgzip_in.BackUp(size);
            } else {
                cerr << "Peek failed!" << endl;
            }
#endif
            
            ::google::protobuf::io::CodedInputStream coded_in(&bgzip_in);
            bool saw_count = coded_in.ReadVarint64((::google::protobuf::uint64*) &count);
            if (!saw_count) {
        
#ifdef debug
                cerr << "Could not read count. Stopping read." << endl;
#endif
                // EOF (probably)
                return;
            }
            
#ifdef debug 
                cerr << "Found group with count " << count << endl;
#endif
        }
        
        // Call the count callback
        handle_count(count);
        
        // Make a shared buffer string to hold message data for each message.
        std::string s;
        for (size_t i = 0; i < count; ++i) {
            uint32_t msgSize = 0;
            
            // Make sure to use a new CodedInputStream every time, because each
            // one limits all input size on the assumption it is reading a
            // single message.
            ::google::protobuf::io::CodedInputStream coded_in(&bgzip_in);
            
            // Alot space for size, and for reading next chunk's length
            coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
            
            // the messages are prefixed by their size. Insist on reading it.
            handle(coded_in.ReadVarint32(&msgSize));
            
#ifdef debug
            cerr << "Found message size of " << msgSize << endl;
#endif
            
            if (msgSize > MAX_PROTOBUF_SIZE) {
                throw std::runtime_error("[stream::for_each] protobuf message of " +
                    std::to_string(msgSize) + " bytes is too long");
            }
            
#ifdef debug
            cerr << "Reading message of " << msgSize << " bytes at " << i << "/" << count << " in group @ " << virtual_offset << endl; 
#endif
            
            // Note that the message may be 0-size, which is a valid (all default values) Protobuf message.
            T object;
            if (msgSize > 0) {
                // Actually need to parse the nonempty message
                handle(coded_in.ReadString(&s, msgSize));
                handle(object.ParseFromString(s));
            }
            // Process the message, passing along the virtual offset of the group, if available
            lambda(virtual_offset, object);
        }
    }
}
    
template <typename T>
void for_each(std::istream& in,
              const std::function<void(int64_t, T&)>& lambda) {
    std::function<void(size_t)> noop = [](size_t) { };
    for_each_with_group_length(in, lambda, noop);
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
                            const std::function<void(size_t)>& handle_count,
                            const std::function<bool(void)>& single_threaded_until_true) {

    // objects will be handed off to worker threads in batches of this many
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
    #pragma omp parallel default(none) shared(in, lambda1, lambda2, handle_count, batches_outstanding, max_batches_outstanding, single_threaded_until_true)
    #pragma omp single
    {
        auto handle = [](bool retval) -> void {
            if (!retval) throw std::runtime_error("obsolete, invalid, or corrupt protobuf input");
        };

        BlockedGzipInputStream bgzip_in(in);
        ::google::protobuf::io::CodedInputStream coded_in(&bgzip_in);

        std::vector<std::string> *batch = nullptr;
        
        // process chunks prefixed by message count
        size_t count;
        while (coded_in.ReadVarint64((::google::protobuf::uint64*) &count)) {
            handle_count(count);
            for (size_t i = 0; i < count; ++i) {
                if (!batch) {
                     batch = new std::vector<std::string>();
                     batch->reserve(batch_size);
                }
                
                // Reconstruct the CodedInputStream in place to reset its maximum-
                // bytes-ever-read counter, because it thinks it's reading a single
                // message.
                coded_in.~CodedInputStream();
                new (&coded_in) ::google::protobuf::io::CodedInputStream(&bgzip_in);
                // Allot space for size, and for reading next chunk's length
                coded_in.SetTotalBytesLimit(MAX_PROTOBUF_SIZE * 2, MAX_PROTOBUF_SIZE * 2);
                
                uint32_t msgSize = 0;
                // the messages are prefixed by their size
                handle(coded_in.ReadVarint32(&msgSize));
                
                if (msgSize > MAX_PROTOBUF_SIZE) {
                    throw std::runtime_error("[stream::for_each] protobuf message of " +
                        std::to_string(msgSize) + " bytes is too long");
                }
                
               
                {
                    std::string s;
                    if (msgSize > 0) {
                        // pick off the message (serialized protobuf object)
                        handle(coded_in.ReadString(&s, msgSize));
                    }
                    // Even empty messages need to be handled; they are all-default Protobuf objects.
                    batch->push_back(std::move(s));
                }

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
    std::function<void(size_t)> no_count = [](size_t i) {};
    std::function<bool(void)> no_wait = [](void) {return true;};
    for_each_parallel_impl(in, lambda2, err1, no_count, no_wait);
}
    
template <typename T>
void for_each_interleaved_pair_parallel_after_wait(std::istream& in,
                                                   const std::function<void(T&,T&)>& lambda2,
                                                   const std::function<bool(void)>& single_threaded_until_true) {
    std::function<void(T&)> err1 = [](T&){
        throw std::runtime_error("stream::for_each_interleaved_pair_parallel: expected input stream of interleaved pairs, but it had odd number of elements");
    };
    std::function<void(size_t)> no_count = [](size_t i) {};
    for_each_parallel_impl(in, lambda2, err1, no_count, single_threaded_until_true);
}

// parallelized for each individual element
template <typename T>
void for_each_parallel(std::istream& in,
                       const std::function<void(T&)>& lambda1,
                       const std::function<void(size_t)>& handle_count) {
    std::function<void(T&,T&)> lambda2 = [&lambda1](T& o1, T& o2) { lambda1(o1); lambda1(o2); };
    std::function<bool(void)> no_wait = [](void) {return true;};
    for_each_parallel_impl(in, lambda2, lambda1, handle_count, no_wait);
}

template <typename T>
void for_each_parallel(std::istream& in,
              const std::function<void(T&)>& lambda) {
    std::function<void(size_t)> noop = [](size_t) { };
    for_each_parallel(in, lambda, noop);
}

}

}

#endif
