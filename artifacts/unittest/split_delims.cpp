/** \file
 *
 * Unit tests for the string helpers in utility.cpp.
 */

#include <string>
#include <vector>

#include "../../src/utility.hpp"

#include "../../src/unittest/catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("split_delims_keep_empty keeps every field in place", "[utility]") {

    SECTION("An empty field between two delimiters is kept") {
        const vector<string> fields = split_delims_keep_empty("a,,b", ",");
        REQUIRE(fields == vector<string>{"a", "", "b"});
    }

    SECTION("An empty string is one empty field") {
        const vector<string> fields = split_delims_keep_empty("", ",");
        REQUIRE(fields == vector<string>{""});
    }

    SECTION("A lone delimiter gives two empty fields") {
        const vector<string> fields = split_delims_keep_empty(",", ",");
        REQUIRE(fields == vector<string>{"", ""});
    }

    SECTION("Any character of the delimiter string splits") {
        const vector<string> fields = split_delims_keep_empty("0|1/.", "/|");
        REQUIRE(fields == vector<string>{"0", "1", "."});
    }

    SECTION("join_delim puts back what split_delims_keep_empty took apart") {
        for (const string& line : {string("a\t\tb\t"), string(""), string("\t"), string("GT:GQ")}) {
            const vector<string> fields = split_delims_keep_empty(line, "\t");
            REQUIRE(join_delim(fields, '\t') == line);
        }
    }
}

}
}
