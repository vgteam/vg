#ifndef VG_STREAM_FD_STREAMS_HPP_INCLUDED
#define VG_STREAM_FD_STREAMS_HPP_INCLUDED

/**
 * \file fd_streams.hpp
 * Contains istream and ostream implementations that operate on file descriptors.
 * Together with POSIX pipe(), this provides an easy way of talking to yourself over streams.
 * This functionality really should be part of C++, but it isn't. See: https://stackoverflow.com/q/2746168
 */
 
#include <streambuf>
#include <iostream>

namespace vg {

namespace stream {

using namespace std;

/**
 * Streambuf implementation that reads from/writes to a file descriptor. The
 * base streambuf class can use an internal buffer ("controlled" sequence), and
 * can handle the reading/writing from that. It can also have an internal
 * buffer where all the pointers are null, and all reads/writes immediately
 * underflow/overflow. There still *is* a controlled sequence, it just
 * immediately commits to the associated sequence.
 *
 * We are only responsible for overriding some virtual functions that fill
 * in/clear out the buffer from/to the backing data source ("associated"
 * sequence).
 *
 * See: https://en.cppreference.com/w/cpp/io/basic_streambuf and
 * http://www.cplusplus.com/reference/streambuf/streambuf/
 *
 * We MUST implement:
 * - overflow() to set the current controlled output sequence character, and
 *   optionally provide some more write buffer space in the controlled output
 *   sequence and commit what was there. 
 * - underflow() to return the current controlled input sequence character, and
 *   optionally buffer some more in the controlled input sequence and rewrite
 *   pointers.
 * - sync() to commit everything written to the controlled output sequence to
 *   the associated output sequence, IF we don't do that automatically (such as
 *   by having an always-empty buffer.
 *
 * We SHOULD implement:
 * - xsputn() to put a bunch of characters into the controlled output sequence
 *   at once, possibly sending them to the associated output sequence.
 * - xsgetn() to get a bunch of characters fromn the controlled input sequence
 *   at once, reading from the associated input sequence if necessary.
 * - pbackfail() to back up the controlled input sequence by one character,
 *   when we are out of cached input sequnce in the buffer, if possible.
 *   Necessary for peek to work.
 *
 * Our current strategy is to not buffer or do putback at all, and just work on
 * individual characters immediately read/written. The next TODO is batch
 * reads/writes.
 */
class FDStreambuf : public streambuf {
public:
    /**
     * Make a new FDStreambuf wrapping the given file descriptor. Undefined
     * behavior will happen if the FD is closed while the streambuf is alive,
     * or if we try to read/write when the file descriptor doesn't support that
     * direction of IO.
     */
    FDStreambuf(int fd);
    
protected:
    
    /// Put the given character in the controlled sequence, if it is not EOF.
    /// Returns something other than EOF on success and EOF on failure.
    int overflow(int c);
    
    /// Get a character and return it, or EOF on failure.
    int underflow();
    
    /// Stores the actual backing file descriptor
    int fd;
};

// Now we define istream and ostream implementations that create, own, and associate an FDStreambuf.

/**
 * C++ istream that reads from a file descriptor.
 */
class FDIstream : public istream {
public:
    /// Wrap a file descriptor in an istream.
    FDIstream(int fd);
protected:
    /// The streambuf implementation we use as a backend.
    FDStreambuf backend;
};

/**
 * C++ ostream that writes to a file descriptor.
 */
class FDOstream : public ostream {
public:
    /// Wrap a file descriptor in an ostream.
    FDOstream(int fd);
protected:
    /// The streambuf implementation we use as a backend.
    FDStreambuf backend;
};

}

}

#endif
