/// \file blocked_gzip_output_stream.cpp
///  
/// Unit tests for BlockedGzipOutputStream and the many adapters that power it.
///

#include "../blocked_gzip_output_stream.hpp"
#include "catch.hpp"

#include <htslib/hfile.h>

#include <sstream>

namespace vg {
namespace unittest {
using namespace std;
using namespace vg::stream;

TEST_CASE("input streams can be wrapped as hFILEs", "[bgzip]") {
    string data = "Candy";
    
    stringstream datastream(data);
    
    hFILE* wrapped = hfile_wrap((std::istream&)datastream);
    
    REQUIRE(wrapped != nullptr);
    
    SECTION("reading works") {
        REQUIRE(hgetc(wrapped) == 'C');
        REQUIRE(hgetc(wrapped) == 'a');
        REQUIRE(hgetc(wrapped) == 'n');
        REQUIRE(hgetc(wrapped) == 'd');
        REQUIRE(hgetc(wrapped) == 'y');
        REQUIRE(hgetc(wrapped) == EOF);
    }
    
    SECTION("seeking and telling works") {
        REQUIRE(htell(wrapped) == 0);
        REQUIRE(hgetc(wrapped) == 'C');
        REQUIRE(htell(wrapped) == 1);
        
        REQUIRE(hseek(wrapped, 1, SEEK_CUR) == 2);
        REQUIRE(htell(wrapped) == 2);
        REQUIRE(hgetc(wrapped) == 'n');
        REQUIRE(htell(wrapped) == 3);
        
        REQUIRE(hseek(wrapped, -1, SEEK_END) == 4);
        REQUIRE(htell(wrapped) == 4);
        REQUIRE(hgetc(wrapped) == 'y');
        REQUIRE(htell(wrapped) == 5);
        
        REQUIRE(hseek(wrapped, 1, SEEK_SET) == 1);
        REQUIRE(hgetc(wrapped) == 'a');
        REQUIRE(htell(wrapped) == 2);
    }
    
    auto close_result = hclose(wrapped);
    REQUIRE(close_result == 0);
    
}

TEST_CASE("output streams can be wrapped as hFILEs", "[bgzip]") {
    stringstream datastream;
    
    hFILE* wrapped = hfile_wrap((std::ostream&)datastream);
    
    REQUIRE(wrapped != nullptr);
    
    SECTION("writing works") {
        REQUIRE(hputc('C', wrapped) == 'C');
        REQUIRE(hputc('a', wrapped) == 'a');
        REQUIRE(hputc('n', wrapped) == 'n');
        REQUIRE(hputc('d', wrapped) == 'd');
        REQUIRE(hputc('y', wrapped) == 'y');
    }
    
    SECTION("seeking and telling works") {
        REQUIRE(htell(wrapped) == 0);
        REQUIRE(hputs("Caddy", wrapped) == 0);
        REQUIRE(htell(wrapped) == 5);
        REQUIRE(hseek(wrapped, -3, SEEK_CUR) == 2);
        REQUIRE(hputc('n', wrapped) == 'n');
    }
    
    auto close_result = hclose(wrapped);
    REQUIRE(close_result == 0);
    
    REQUIRE(datastream.str() == "Candy");
    
}

TEST_CASE("hFILE stream position is preserved across opens and closes", "[bgzip]") {
    // Make a stream that we are not at the beginning of
    stringstream datastream;
    datastream << "datadatadata";
    
    hFILE* wrapped = hfile_wrap((std::ostream&)datastream);
    REQUIRE(wrapped != nullptr);
    REQUIRE(htell(wrapped) == 12);
    
    REQUIRE(hputc('\n', wrapped) == '\n');
    REQUIRE(htell(wrapped) == 13);
    
    REQUIRE(hclose(wrapped) == 0);
    
    REQUIRE(datastream.str() == "datadatadata\n");
    
    wrapped = hfile_wrap((std::ostream&)datastream);
    REQUIRE(wrapped != nullptr);
    REQUIRE(htell(wrapped) == 13);
    
    REQUIRE(hclose(wrapped) == 0);
    
}

// We have a tiny function to get virtual offsets, based on the block's start
// offset in the file, and the offset in the block
int64_t vo(size_t block_start, size_t offset) {
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


}

}
