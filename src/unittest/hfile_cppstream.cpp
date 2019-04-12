/// \file hfile_cppstream.cpp
///  
/// Unit tests for hFILE* C++ stream wrapper
///

#include <vg/io/hfile_cppstream.hpp>
#include "catch.hpp"

#include <htslib/hfile.h>

#include <sstream>

namespace vg {
namespace unittest {
using namespace std;
using namespace vg::io;

TEST_CASE("input streams can be wrapped as hFILEs", "[bgzip][hfile]") {
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

TEST_CASE("output streams can be wrapped as hFILEs", "[bgzip][hfile]") {
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

TEST_CASE("hFILE stream position is preserved across opens and closes", "[bgzip][hfile]") {
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

}

}
