/// \file blocked_gzip_input_stream.cpp
///  
/// Unit tests for BlockedGzipInputStream

#include "../blocked_gzip_input_stream.hpp"
#include "../blocked_gzip_output_stream.hpp"
#include "../hfile_cppstream.hpp"
#include "catch.hpp"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>


#include <htslib/hfile.h>

#include <sstream>

namespace vg {
namespace unittest {
using namespace std;
using namespace vg::stream;

// We have a tiny function to get virtual offsets, based on the block's start
// offset in the file, and the offset in the block
static int64_t vo(size_t block_start, size_t offset) {
    return (block_start << 16) | (0xFFFF & offset);
}

TEST_CASE("a BlockedGzipInputStream can read from a stringstream", "[bgzip]") {
    stringstream datastream;
    
    string TO_COMPRESS = "It's cheap and it's ethical... well, it's ethical... well, it's magical really";
    
    {
        // Write some data in
        BlockedGzipOutputStream bgzip_out(datastream);
        ::google::protobuf::io::CodedOutputStream coded_out(&bgzip_out);
        coded_out.WriteString(TO_COMPRESS);
    }
    
    SECTION("data can be read the first time through, from the start") {
    
        // Now try and read it back
        BlockedGzipInputStream bgzip_in(datastream);
    
        REQUIRE(bgzip_in.Tell() == vo(0, 0));
        
        const char* buffer;
        int buffer_size;
        REQUIRE(bgzip_in.Next((const void**)&buffer, &buffer_size));
        
        int block = 0;
        int good_through = 0;
        
        do {
        
            // Check each block we read out of the stream
            
            // We know that the stream ought to put us at the end of whatever it read.
            // We also know it ought to read one block per Next if we don't back up.
            REQUIRE(bgzip_in.Tell() == vo(block, buffer_size));
            
            for (size_t i = 0; i < buffer_size; i++) {
                // Check all the characters
                REQUIRE(buffer[i] == TO_COMPRESS[good_through + i]);
            }
            
            good_through += buffer_size;
            block++;
            
        } while (bgzip_in.Next((const void**)&buffer, &buffer_size));
        
        REQUIRE(good_through == TO_COMPRESS.size());
    }
    
    SECTION("data can be seeked into") {
        BlockedGzipInputStream bgzip_in(datastream);
    
        // Make sure we started at the start
        REQUIRE(bgzip_in.Tell() == vo(0, 0));
        
        // Go somewhere else and make sure we get there
        REQUIRE(bgzip_in.Seek(vo(0, 10)));
        REQUIRE(bgzip_in.Tell() == vo(0, 10));
        
        const char* buffer;
        int buffer_size = 0;
        while(buffer_size == 0) {
            // Fish for data until we get some
            REQUIRE(bgzip_in.Next((const void**)&buffer, &buffer_size));
        }
        
        // Make sure we got the right data
        REQUIRE(buffer[0] == TO_COMPRESS[10]);
        
        SECTION("seek can seek back") {
        
            REQUIRE(bgzip_in.Seek(vo(0, 0)));
            REQUIRE(bgzip_in.Tell() == vo(0, 0));
            
            buffer_size = 0;
            while(buffer_size == 0) {
                // Fish for data until we get some
                REQUIRE(bgzip_in.Next((const void**)&buffer, &buffer_size));
            }
            
            // Make sure we got the right data
            REQUIRE(buffer[0] == TO_COMPRESS[0]);
            
        }
        
        
    }
        
}


}

}
