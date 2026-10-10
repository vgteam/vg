/** \file
 *
 * Unit tests for the string helpers in utility.cpp.
 */

#include <string>
#include <vector>

#include "../utility.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("split_delims_keep_empty keeps every field in place", "[utility]") {

    SECTION("An empty field between two delimiters is kept") {
        vector<string> fields;
        split_delims_keep_empty("a,,b", ",", fields);
        REQUIRE(fields == vector<string>{"a", "", "b"});
    }

    SECTION("An empty string is one empty field") {
        vector<string> fields;
        split_delims_keep_empty("", ",", fields);
        REQUIRE(fields == vector<string>{""});
    }

    SECTION("A lone delimiter gives two empty fields") {
        vector<string> fields;
        split_delims_keep_empty(",", ",", fields);
        REQUIRE(fields == vector<string>{"", ""});
    }

    SECTION("Any character of the delimiter string splits") {
        vector<string> fields;
        split_delims_keep_empty("0|1/.", "/|", fields);
        REQUIRE(fields == vector<string>{"0", "1", "."});
    }

    SECTION("Fields are appended to what the vector already holds") {
        vector<string> fields{"x"};
        split_delims_keep_empty("a:b", ":", fields);
        REQUIRE(fields == vector<string>{"x", "a", "b"});
    }

    SECTION("join_delim puts back what split_delims_keep_empty took apart") {
        for (const string& line : {string("a\t\tb\t"), string(""), string("\t"), string("GT:GQ")}) {
            vector<string> fields;
            split_delims_keep_empty(line, "\t", fields);
            REQUIRE(join_delim(fields, '\t') == line);
        }
    }
}

}
}
