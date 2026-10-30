/// \file site_scheduler.cpp
///
/// Unit tests for the site-calling scheduler's progress tally.

#include "../../src/unittest/catch.hpp"
#include "../site_scheduler.hpp"

#include <sstream>

namespace vg {
namespace unittest {

TEST_CASE("SiteCallTally reports the nested count in nested progress lines", "[site_scheduler]") {
    std::stringstream captured;
    std::streambuf* old_buffer = std::cerr.rdbuf(captured.rdbuf());

    SiteCallTally tally;
    tally.count();
    tally.finish_top_level();
    for (size_t i = 0; i < 100000; i++) {
        tally.count();
    }
    tally.finish_nested();

    std::cerr.rdbuf(old_buffer);

    std::string log = captured.str();
    REQUIRE(log.find("Processed 100000 nested snarls") != std::string::npos);
    REQUIRE(log.find("Finished processing 100000 nested snarls") != std::string::npos);
}

}
}
