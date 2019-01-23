/**
 * \file fd_streams.cpp
 * Implementations for file descriptor streams.
 */

#include "fd_streams.hpp"

#include <unistd.h>

namespace vg {

namespace stream {

using namespace std;

fdstreambuf::fdstreambuf(int fd) : fd(fd) {
    // Nothing to do!
    // Leave all the buffers empty/nonexistent
}

auto fdstreambuf::overflow(int c) -> int {

    if (c != traits_type::eof()) {
        // It's not EOF, so write the character
        char byte = (char) c;
        // Do the write and return EOF if it failed.
        return (write(fd, &byte, 1) == 1) ? c : traits_type::eof();
    } else {
        // We got passed EOF. Just EOF right back.
        return c;
    }

}

auto fdstreambuf::underflow() -> int {
    char byte;
    // Read a byte, and return it if successful or EOF if failed.
    return (read(fd, &byte, 1) == 1) ? (int) byte : traits_type::eof();
}

fdistream::fdistream(int fd) : istream(nullptr), backend(fd) {
    // Now that buf is constructed, associate ourselves with it
    rdbuf(&backend);
}

fdostream::fdostream(int fd) : ostream(nullptr), backend(fd) {
    // Now that buf is constructed, associate ourselves with it
    rdbuf(&backend);
}

}

}
