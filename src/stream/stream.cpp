#include "stream.hpp"
#include "blocked_gzip_output_stream.hpp"

namespace vg {

namespace stream {

using namespace std;

void finish(std::ostream& out) {
    // Put an EOF on the stream by making a writer, marking it as EOF, and letting it clean up.
    BlockedGzipOutputStream bgzip_out(out);
    bgzip_out.EndFile();
}

}

}

