/// \file blocked_gzip_output_stream.cpp
///  
/// Unit tests for BlockedGzipOutputStream 

#include <vg/io/blocked_gzip_output_stream.hpp>
#include <vg/io/hfile_cppstream.hpp>
#include "catch.hpp"

#include <htslib/hfile.h>

#include <sstream>

namespace vg {
namespace unittest {
using namespace std;
using namespace vg::io;

// We have a tiny function to get virtual offsets, based on the block's start
// offset in the file, and the offset in the block
static int64_t vo(size_t block_start, size_t offset) {
    return (block_start << 16) | (0xFFFF & offset);
}

TEST_CASE("a BlockedGzipOutputStream can write to a stringstream", "[bgzip]") {
    stringstream datastream;
    
    SECTION("when wrapped manually") {
        // Wrap the stream as an hFILE
        hFILE* wrapped = hfile_wrap((std::ostream&)datastream);
        REQUIRE(wrapped != nullptr);
        // Give ownership of it to a BGZF that writes
        BGZF* bgzf_handle = bgzf_hopen(wrapped, "w");
        REQUIRE(bgzf_handle != nullptr);
        // Give ownership of the BGZF to a BlockedGzipOutputStream
        BlockedGzipOutputStream bgzip_out(bgzf_handle);
        
        REQUIRE(bgzip_out.Tell() == vo(0, 0));
        
        char* buffer;
        int buffer_size;
        REQUIRE(bgzip_out.Next((void**)&buffer, &buffer_size));
        
        REQUIRE(buffer_size >= 5);
        buffer[0] = 'C';
        buffer[1] = 'a';
        buffer[2] = 'n';
        buffer[3] = 'd';
        buffer[4] = 'y';
        
        bgzip_out.BackUp(buffer_size - 5);
        
        REQUIRE(bgzip_out.Tell() == vo(0, 5));
    }
    
    SECTION("when wrapped automatically") {
        BlockedGzipOutputStream bgzip_out(datastream);
        
        REQUIRE(bgzip_out.Tell() == vo(0, 0));
        
        char* buffer;
        int buffer_size;
        REQUIRE(bgzip_out.Next((void**)&buffer, &buffer_size));
        
        REQUIRE(buffer_size >= 5);
        buffer[0] = 'C';
        buffer[1] = 'a';
        buffer[2] = 'n';
        buffer[3] = 'd';
        buffer[4] = 'y';
        
        bgzip_out.BackUp(buffer_size - 5);
        
        REQUIRE(bgzip_out.Tell() == vo(0, 5));
    }
    
    // We must have written something (compressed data) to the stream
    REQUIRE(datastream.str().size() != 0);
    
}

TEST_CASE("BlockedGzipOutputStream can write to the same stream multiple times", "[bgzip]") {

    stringstream datastream;
    
    {
        // Write the first block
        BlockedGzipOutputStream bgzip_out(datastream);
         
        REQUIRE(bgzip_out.Tell() == vo(0, 0));
         
        char* buffer;
        int buffer_size;
        REQUIRE(bgzip_out.Next((void**)&buffer, &buffer_size));
        
        REQUIRE(buffer_size >= 5);
        buffer[0] = 'C';
        buffer[1] = 'a';
        buffer[2] = 'n';
        buffer[3] = 'd';
        buffer[4] = 'y';
        
        bgzip_out.BackUp(buffer_size - 5);
        
        REQUIRE(bgzip_out.Tell() == vo(0, 5));
    }
    
    // Now the stream goes away. The BGZF must be gone, so all the data must be done.
    
    size_t data_size = datastream.str().size();
    
    {
        // We have to be creating a new block now
        BlockedGzipOutputStream bgzip_out(datastream);
        
        // Make sure the stream reports that it is starting a block at the position the file is at, and not at 0.
        REQUIRE(bgzip_out.Tell() == vo(data_size, 0));
         
        char* buffer;
        int buffer_size;
        REQUIRE(bgzip_out.Next((void**)&buffer, &buffer_size));
        
        REQUIRE(buffer_size >= 5);
        buffer[0] = 'D';
        buffer[1] = 'a';
        buffer[2] = 'n';
        buffer[3] = 'd';
        buffer[4] = 'y';
        
        bgzip_out.BackUp(buffer_size - 5);
        
        REQUIRE(bgzip_out.Tell() == vo(data_size, 5));
    }
    
    REQUIRE(datastream.str().size() == data_size * 2);
    

}

TEST_CASE("BlockedGzipOutputStream produces EOF markers on demand", "[bgzip]") {
    stringstream s1;
    stringstream s2;
    
    {
        // Make some streams
        BlockedGzipOutputStream bs1(s1);
        BlockedGzipOutputStream bs2(s2);
        
        // Fake some writes
        char* data;
        int size;
        bs1.Next((void**)&data, &size);
        for (size_t i = 0; i < size; i++) {
            data[i] = 0;
        }
        
        bs2.Next((void**)&data, &size);
        for (size_t i = 0; i < size; i++) {
            data[i] = 0;
        }
        
        // Only give one an EOF
        bs2.EndFile();
    }
    
    // Make sure the extra 28-byte EOF marker is present
    REQUIRE(s2.str().length() == s1.str().length() + 28);
}

}

}
