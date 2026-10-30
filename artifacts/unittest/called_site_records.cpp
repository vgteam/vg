/// \file called_site_records.cpp
/// Unit tests for building the VCF record of a called site with a SiteEmitter.

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>

#include "../../src/unittest/catch.hpp"
#include "../../src/site_emitter.hpp"
#include "../../src/utility.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A bubble: node 1 (A), then node 2 (C) or node 3 (G), then node 4 (T), with the reference path
/// "ref" through node 2.
struct Bubble {
    bdsg::HashGraph base;
    handle_t n1, n2, n3, n4;
    Bubble() {
        n1 = base.create_handle("A", 1);
        n2 = base.create_handle("C", 2);
        n3 = base.create_handle("G", 3);
        n4 = base.create_handle("T", 4);
        base.create_edge(n1, n2);
        base.create_edge(n1, n3);
        base.create_edge(n2, n4);
        base.create_edge(n3, n4);
        path_handle_t ref = base.create_path_handle("ref");
        for (handle_t h : {n1, n2, n4}) {
            base.append_step(ref, h);
        }
    }
};

/// Return the data lines `output` writes, split into columns.
vector<vector<string>> written_lines(SortedVCFWriter& output) {
    stringstream out;
    output.write_variants(out);
    vector<vector<string>> lines;
    string line;
    while (getline(out, line)) {
        if (!line.starts_with("#")) {
            lines.push_back(split_delims(line, "\t"));
        }
    }
    return lines;
}

}

TEST_CASE("A SiteEmitter writes a called site as one record", "[site_emitter]") {
    Bubble bubble;
    bdsg::PositionOverlay graph(&bubble.base);
    SortedVCFWriter output("sample");
    SiteEmitter emitter(output);

    // The two traversals, and how often the emitter asks for each.
    const vector<pair<string, vector<NodeVisit>>> spellings{
        {"ACT", {{1, false}, {2, false}, {4, false}}},
        {"AGT", {{1, false}, {3, false}, {4, false}}},
    };
    map<size_t, int> asked;
    allele_source_t alleles = [&](size_t allele) {
        ++asked[allele];
        return spellings.at(allele);
    };
    // The traversal of each VCF allele that the field writer is given.
    vector<int> written_traversals;
    field_writer_t write_fields = [&](const vector<int>& allele_traversals, const vector<int>&,
                                      vcflib::Variant&) {
        written_traversals = allele_traversals;
    };
    const string ref_path = "ref";

    SECTION("a heterozygous call gets a record trimmed to the bases that differ") {
        const vector<int> genotype{0, 1};
        CalledSite site{bubble.n1, bubble.n4, ref_path, 0, genotype, 0, 2, 2};
        REQUIRE(emitter.emit(graph, site, alleles, write_fields));
        REQUIRE(written_traversals == vector<int>({0, 1}));
        // Each traversal is spelled once, however many steps read its spelling.
        REQUIRE(asked == map<size_t, int>({{0, 1}, {1, 1}}));
        vector<vector<string>> lines = written_lines(output);
        REQUIRE(lines.size() == 1);
        REQUIRE(lines[0][0] == "ref");
        REQUIRE(lines[0][1] == "2");
        REQUIRE(lines[0][2] == ">1>4");
        REQUIRE(lines[0][3] == "C");
        REQUIRE(lines[0][4] == "G");
        REQUIRE(lines[0][9] == "0/1");
    }

    SECTION("a reference call is hidden") {
        const vector<int> genotype{0, 0};
        CalledSite site{bubble.n1, bubble.n4, ref_path, 0, genotype, 0, 2, 2};
        REQUIRE(emitter.emit(graph, site, alleles, write_fields));
        REQUIRE(written_lines(output).empty());
    }

    SECTION("a reference call is written with its uncalled alleles when genotyping snarls") {
        emitter.genotype_snarls = true;
        const vector<int> genotype{0, 0};
        CalledSite site{bubble.n1, bubble.n4, ref_path, 0, genotype, 0, 2, 2};
        REQUIRE(emitter.emit(graph, site, alleles, write_fields));
        REQUIRE(written_traversals == vector<int>({0, 1}));
        vector<vector<string>> lines = written_lines(output);
        REQUIRE(lines.size() == 1);
        REQUIRE(lines[0][3] == "C");
        REQUIRE(lines[0][4] == "G");
        REQUIRE(lines[0][9] == "0/0");
    }

    SECTION("a site with an ID of its own gets that ID") {
        const vector<int> genotype{1, 1};
        const string id = "op1";
        CalledSite site{bubble.n1, bubble.n4, ref_path, 0, genotype, 0, 2, 2, &id};
        REQUIRE(emitter.emit(graph, site, alleles, write_fields));
        vector<vector<string>> lines = written_lines(output);
        REQUIRE(lines.size() == 1);
        REQUIRE(lines[0][2] == "op1");
        REQUIRE(lines[0][9] == "1/1");
    }
}

}
}
