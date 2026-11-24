/// \file unittest/ploidy_regions.cpp
///
/// Unit tests for PloidyRegions and the ploidy presets built into vg.

#include <sstream>
#include <string>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../ploidy_regions.hpp"

namespace vg {
namespace unittest {

TEST_CASE("A ploidy BED read from a stream gives each region its ploidy", "[ploidy_regions]") {
    istringstream bed("# a comment\nx\t10\t20\t1\nx\t30\t40\t0\n");
    PloidyRegions table(bed, "test BED");

    REQUIRE(table.region_ploidy("x", 9, 2) == 2);
    REQUIRE(table.region_ploidy("x", 10, 2) == 1);
    REQUIRE(table.region_ploidy("x", 19, 2) == 1);
    REQUIRE(table.region_ploidy("x", 20, 2) == 2);
    REQUIRE(table.region_ploidy("x", 30, 2) == 0);
    REQUIRE(table.region_ploidy("y", 15, 2) == 2);
}

TEST_CASE("A site is in a ploidy BED interval when the first base of its record is", "[ploidy_regions]") {
    // [0, 50) in BED's 0-based coordinates holds VCF POS 1 to 50.
    istringstream bed("x\t0\t50\t1\n");
    PloidyRegions table(bed, "test BED");

    SECTION("A record at POS 50 is inside the interval, and one at POS 51 is not") {
        REQUIRE(table.ploidy_at("x", 49, 0, 2) == 1);
        REQUIRE(table.ploidy_at("x", 50, 0, 2) == 2);
    }
    SECTION("The reference offset moves the site along the contig") {
        REQUIRE(table.ploidy_at("x", 40, 9, 2) == 1);
        REQUIRE(table.ploidy_at("x", 40, 10, 2) == 2);
    }
    SECTION("A subrange path is placed at its start on the contig") {
        REQUIRE(table.ploidy_at("x[30]", 19, 0, 2) == 1);
        REQUIRE(table.ploidy_at("x[30]", 20, 0, 2) == 2);
    }
    SECTION("A site before the start of the contig gets the fallback") {
        REQUIRE(table.ploidy_at("x", 0, -5, 2) == 2);
    }
}

TEST_CASE("Every ploidy preset parses and makes the mitochondrion haploid", "[ploidy_regions]") {
    vector<string> names = PloidyRegions::preset_names();
    REQUIRE(names == vector<string>{"chm13v2.XX", "chm13v2.XY", "grch37.b37.XX", "grch37.b37.XY",
                                    "grch37.ucsc.XX", "grch37.ucsc.XY", "grch38.XX", "grch38.XY"});
    for (const string& name : names) {
        PloidyRegions table = PloidyRegions::preset(name);
        // The b37 names drop the "chr" prefix, and call the mitochondrion MT.
        string mitochondrion = name.starts_with("grch37.b37.") ? "MT" : "chrM";
        INFO(name);
        REQUIRE(table.region_ploidy(mitochondrion, 0, 2) == 1);
        REQUIRE(table.region_ploidy("chr1", 1000000, 2) == 2);
    }
}

TEST_CASE("The chm13v2 presets give chrX and chrY their ploidy", "[ploidy_regions]") {
    SECTION("XY makes chrX haploid between its pseudoautosomal regions, and chrY haploid") {
        PloidyRegions table = PloidyRegions::preset("chm13v2.XY");
        REQUIRE(table.region_ploidy("chrX", 2394369, 2) == 2);
        REQUIRE(table.region_ploidy("chrX", 2394370, 2) == 1);
        REQUIRE(table.region_ploidy("chrX", 153926002, 2) == 1);
        REQUIRE(table.region_ploidy("chrX", 153926003, 2) == 2);
        REQUIRE(table.region_ploidy("chrY", 0, 2) == 1);
        REQUIRE(table.region_ploidy("chrY", 62460028, 2) == 1);
        // A PanSN path name is looked up by its locus.
        REQUIRE(table.region_ploidy("CHM13#0#chrX", 2394370, 2) == 1);
    }
    SECTION("XX leaves chrX diploid and gives chrY ploidy 0") {
        PloidyRegions table = PloidyRegions::preset("chm13v2.XX");
        REQUIRE(table.region_ploidy("chrX", 2394370, 2) == 2);
        REQUIRE(table.region_ploidy("chrY", 0, 2) == 0);
        REQUIRE(table.region_ploidy("chrY", 62460028, 2) == 0);
    }
}

TEST_CASE("The grch37.b37 presets name chromosomes X, Y and MT", "[ploidy_regions]") {
    PloidyRegions xy = PloidyRegions::preset("grch37.b37.XY");
    REQUIRE(xy.region_ploidy("X", 2699520, 2) == 1);
    REQUIRE(xy.region_ploidy("chrX", 2699520, 2) == 2);
    REQUIRE(xy.region_ploidy("Y", 0, 2) == 1);
    REQUIRE(xy.region_ploidy("MT", 0, 2) == 1);

    PloidyRegions xx = PloidyRegions::preset("grch37.b37.XX");
    REQUIRE(xx.region_ploidy("X", 2699520, 2) == 2);
    REQUIRE(xx.region_ploidy("Y", 0, 2) == 0);
    REQUIRE(xx.region_ploidy("chrY", 0, 2) == 2);
}

}
}
