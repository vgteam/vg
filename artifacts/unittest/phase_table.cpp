/// \file phase_table.cpp
/// Unit tests for writing a record's GT from the phase call of its site.

#include <map>
#include <string>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../phase_table.hpp"

namespace vg {
namespace unittest {

using namespace std;
using multipass::phased_gt;

/// Make a phase call that puts traversal `first` on strand 0 and `second` on strand 1.
static LinkageCollector::PhaseCall phase_call(int first, int second, size_t ploidy = 2) {
    LinkageCollector::PhaseCall call;
    call.trav_first = first;
    call.trav_second = second;
    call.ploidy = ploidy;
    return call;
}

TEST_CASE("A diploid record is phased in the order of the strands of its phase call",
          "[phase_table]") {
    // Traversal 5 is allele 0 and traversal 7 is allele 1.
    map<int, int> trav_to_allele{{5, 0}, {7, 1}};
    string gt = "0/1";
    REQUIRE(phased_gt(phase_call(7, 5), {0, 1}, trav_to_allele, gt));
    REQUIRE(gt == "1|0");
}

TEST_CASE("A phase call whose alleles differ from the record's genotype does not phase it",
          "[phase_table]") {
    map<int, int> trav_to_allele{{5, 0}, {7, 1}, {9, 2}};

    SECTION("an allele the record does not call") {
        string gt = "0/1";
        REQUIRE_FALSE(phased_gt(phase_call(9, 5), {0, 1}, trav_to_allele, gt));
        REQUIRE(gt == "0/1");
    }

    SECTION("a traversal with no allele in the record") {
        string gt = "0/1";
        REQUIRE_FALSE(phased_gt(phase_call(11, 5), {0, 1}, trav_to_allele, gt));
        REQUIRE(gt == "0/1");
    }

    SECTION("a ploidy other than the record's") {
        string gt = "1";
        REQUIRE_FALSE(phased_gt(phase_call(7, 5), {1}, trav_to_allele, gt));
        REQUIRE(gt == "1");
    }
}

TEST_CASE("A one-strand record is phased onto its strand", "[phase_table]") {
    map<int, int> trav_to_allele{{5, 0}, {7, 1}};

    SECTION("a site of a nested chain carried by strand 1 names that strand") {
        LinkageCollector::PhaseCall call = phase_call(7, -1, 1);
        call.nested_strand = 1;
        string gt = "1";
        REQUIRE(phased_gt(call, {1}, trav_to_allele, gt));
        REQUIRE(gt == ".|1");
    }

    SECTION("a haploid site keeps its one allele") {
        string gt = "1";
        REQUIRE(phased_gt(phase_call(7, -1, 1), {1}, trav_to_allele, gt));
        REQUIRE(gt == "1");
    }
}

}
}
