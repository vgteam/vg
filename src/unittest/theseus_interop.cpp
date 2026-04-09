/// \file unittest/theseus_interop.cpp
///
/// Unit tests for vg_alignment_from_theseus_alignment in theseus_interop.hpp.

#include "catch.hpp"
#include "../theseus_interop.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {

using namespace std;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a theseus::Alignment directly from components.
static theseus::Alignment make_aln(vector<int> path, vector<char> edit_op,
                                   int start_offset, int end_offset) {
    theseus::Alignment aln;
    aln.path         = std::move(path);
    aln.edit_op      = std::move(edit_op);
    aln.start_offset = start_offset;
    aln.end_offset   = end_offset;
    return aln;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
TEST_CASE("HandleGraphTheseusAdapter: theseus::Graph lines up with HandleGraph",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    vector<bdsg::handle_t> handles = {
        hg.create_handle("GAT", 1),
        hg.create_handle("CAGT", 2),
        hg.create_handle("T", 3),
        hg.create_handle("AAAAA", 4)
    };
    hg.create_edge(handles[0], handles[1]);
    hg.create_edge(handles[0], handles[2]);
    hg.create_edge(handles[1], handles[3]);
    hg.create_edge(handles[2], handles[3]);

    HandleGraphTheseusAdapter adapter(hg);
    theseus::Graph theseus_graph = adapter.take_graph();

    REQUIRE(theseus_graph._vertices.size() == 8); // 4 handles * 2 orientations

    SECTION("theseus vertex names are correct") {
        vector<string> expected_names = {
            "1+", "1-", "2+", "2-", "3+", "3-", "4+", "4-"
        };
        vector<string> theseus_vertex_names; 
    
        for (theseus::Graph::vertex& v: theseus_graph.vertices()) {
            theseus_vertex_names.push_back(v.name);
            cerr << "Vertex " << v.name << ": " << v.value << endl;
        }
        sort(theseus_vertex_names.begin(), theseus_vertex_names.end());
        REQUIRE(theseus_vertex_names == expected_names);
    }

    SECTION("theseus vertex sequences are correct") {
        vector<string> expected_sequences = {
            "GAT", "ATC", "CAGT", "ACTG", "T", "A", "AAAAA", "TTTTT"
        };
        sort(expected_sequences.begin(), expected_sequences.end());
        vector<string> theseus_vertex_sequences; 
    
        for (theseus::Graph::vertex& v: theseus_graph.vertices()) {
            theseus_vertex_sequences.push_back(v.value);
            cerr << "Vertex " << v.name << ": " << v.value << endl;
        }
        sort(theseus_vertex_sequences.begin(), theseus_vertex_sequences.end());
        REQUIRE(theseus_vertex_sequences == expected_sequences);
    }

    SECTION("theseus edges are correct") {
        vector<pair<string, string>> expected_edges = {
            {"1+", "2+"}, {"1+", "3+"}, {"2+", "4+"}, {"3+", "4+"},
            {"2-", "1-"}, {"3-", "1-"}, {"4-", "2-"}, {"4-", "3-"}
        };
        sort(expected_edges.begin(), expected_edges.end());
        vector<pair<string, string>> theseus_edges;
    
        for (theseus::Graph::vertex& v: theseus_graph.vertices()) {
            for (theseus::Graph::edge& e: v.out_edges) {
                string from_name = v.name;
                string to_name   = theseus_graph._vertices[e.to_vertex].name;
                theseus_edges.emplace_back(from_name, to_name);
                cerr << "Edge from " << from_name << " to " << to_name << endl;
            }
        }
        sort(theseus_edges.begin(), theseus_edges.end());
        REQUIRE(theseus_edges == expected_edges);
    }
}

/*
TEST_CASE("vg_alignment_from_theseus_alignment: empty alignment yields empty path",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    theseus::Alignment aln = make_aln({}, {}, 0, 0);
    string sequence = "GATTACA";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.sequence() == sequence);
    REQUIRE(result.path().mapping_size() == 0);
}

TEST_CASE("vg_alignment_from_theseus_alignment: empty path with non-empty edit_op yields empty path",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    // The function guards on path.empty() || edit_op.empty()
    theseus::Alignment aln = make_aln({}, {'M', 'M'}, 0, 2);
    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, "GA", adapter, hg);

    REQUIRE(result.path().mapping_size() == 0);
}

TEST_CASE("vg_alignment_from_theseus_alignment: single-node all-match alignment",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);
    theseus::Alignment aln = make_aln({vtx},
        {'M','M','M','M','M','M','M'}, 0, 7);
    string sequence = "GATTACA";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.sequence() == sequence);
    REQUIRE(result.path().mapping_size() == 1);

    const auto& m = result.path().mapping(0);
    REQUIRE(m.position().node_id() == 1);
    REQUIRE(m.position().offset() == 0);
    REQUIRE(m.position().is_reverse() == false);
    REQUIRE(m.rank() == 1);

    REQUIRE(m.edit_size() == 1);
    const auto& e = m.edit(0);
    REQUIRE(e.from_length() == 7);
    REQUIRE(e.to_length() == 7);
    REQUIRE(e.sequence() == "");   // empty sequence field = match
}

TEST_CASE("vg_alignment_from_theseus_alignment: consecutive same-type ops are merged",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);

    SECTION("all matches produce one edit") {
        // 7 separate M characters → should collapse to 1 edit
        theseus::Alignment aln = make_aln({vtx},
            {'M','M','M','M','M','M','M'}, 0, 7);
        vg::Alignment result =
            vg_alignment_from_theseus_alignment(aln, "GATTACA", adapter, hg);
        REQUIRE(result.path().mapping(0).edit_size() == 1);
        REQUIRE(result.path().mapping(0).edit(0).from_length() == 7);
    }

    SECTION("M run then X run produces two edits") {
        // GATTACA: first 3 match, last 4 are mismatches
        theseus::Alignment aln = make_aln({vtx},
            {'M','M','M','X','X','X','X'}, 0, 7);
        vg::Alignment result =
            vg_alignment_from_theseus_alignment(aln, "GATCCCC", adapter, hg);
        const auto& m = result.path().mapping(0);
        REQUIRE(m.edit_size() == 2);
        REQUIRE(m.edit(0).from_length() == 3);
        REQUIRE(m.edit(0).to_length() == 3);
        REQUIRE(m.edit(0).sequence() == "");   // match
        REQUIRE(m.edit(1).from_length() == 4);
        REQUIRE(m.edit(1).to_length() == 4);
        REQUIRE(m.edit(1).sequence() == "CCCC");
    }
}

TEST_CASE("vg_alignment_from_theseus_alignment: start_offset shifts position on first node",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);
    // Align only bases [2,5) of the node ("TTA" within "GATTACA")
    theseus::Alignment aln = make_aln({vtx}, {'M','M','M'}, 2, 5);

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, "TTA", adapter, hg);

    REQUIRE(result.path().mapping_size() == 1);
    const auto& m = result.path().mapping(0);
    REQUIRE(m.position().node_id() == 1);
    REQUIRE(m.position().offset() == 2);
    REQUIRE(m.edit_size() == 1);
    REQUIRE(m.edit(0).from_length() == 3);
    REQUIRE(m.edit(0).to_length() == 3);
}

TEST_CASE("vg_alignment_from_theseus_alignment: mismatch (X) edit has non-empty sequence",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);
    // Single-base mismatch at position 0
    theseus::Alignment aln = make_aln({vtx}, {'X'}, 0, 1);

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, "C", adapter, hg);

    REQUIRE(result.path().mapping_size() == 1);
    const auto& e = result.path().mapping(0).edit(0);
    REQUIRE(e.from_length() == 1);
    REQUIRE(e.to_length() == 1);
    REQUIRE(e.sequence() == "C");
}

TEST_CASE("vg_alignment_from_theseus_alignment: insertion edit has from_length=0",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);
    // Two inserted bases followed by three matches
    theseus::Alignment aln = make_aln({vtx}, {'I','I','M','M','M'}, 0, 3);
    string sequence = "XXGAT";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.path().mapping_size() == 1);
    const auto& m = result.path().mapping(0);
    REQUIRE(m.edit_size() == 2);

    // First edit: insertion
    const auto& ins = m.edit(0);
    REQUIRE(ins.from_length() == 0);
    REQUIRE(ins.to_length() == 2);
    REQUIRE(ins.sequence() == "XX");

    // Second edit: match
    const auto& mat = m.edit(1);
    REQUIRE(mat.from_length() == 3);
    REQUIRE(mat.to_length() == 3);
    REQUIRE(mat.sequence() == "");
}

TEST_CASE("vg_alignment_from_theseus_alignment: deletion edit has to_length=0",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("GATTACA", 1);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx = adapter.handle_to_vertex(h);
    // Two matches then two deletions
    theseus::Alignment aln = make_aln({vtx}, {'M','M','D','D'}, 0, 4);
    string sequence = "GA";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.path().mapping_size() == 1);
    const auto& m = result.path().mapping(0);
    REQUIRE(m.edit_size() == 2);

    const auto& mat = m.edit(0);
    REQUIRE(mat.from_length() == 2);
    REQUIRE(mat.to_length() == 2);

    const auto& del = m.edit(1);
    REQUIRE(del.from_length() == 2);
    REQUIRE(del.to_length() == 0);
}

TEST_CASE("vg_alignment_from_theseus_alignment: multi-node alignment produces one mapping per node",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h1 = hg.create_handle("ACG", 1);
    bdsg::handle_t h2 = hg.create_handle("TTT", 2);
    hg.create_edge(h1, h2);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx1 = adapter.handle_to_vertex(h1);
    int vtx2 = adapter.handle_to_vertex(h2);

    theseus::Alignment aln = make_aln({vtx1, vtx2},
        {'M','M','M','M','M','M'}, 0, 3);
    string sequence = "ACGTTT";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.path().mapping_size() == 2);

    const auto& m1 = result.path().mapping(0);
    REQUIRE(m1.position().node_id() == 1);
    REQUIRE(m1.position().offset() == 0);
    REQUIRE(m1.rank() == 1);
    REQUIRE(m1.edit_size() == 1);
    REQUIRE(m1.edit(0).from_length() == 3);
    REQUIRE(m1.edit(0).to_length() == 3);

    const auto& m2 = result.path().mapping(1);
    REQUIRE(m2.position().node_id() == 2);
    REQUIRE(m2.position().offset() == 0);
    REQUIRE(m2.rank() == 2);
    REQUIRE(m2.edit_size() == 1);
    REQUIRE(m2.edit(0).from_length() == 3);
    REQUIRE(m2.edit(0).to_length() == 3);
}

TEST_CASE("vg_alignment_from_theseus_alignment: edits are split at node boundaries",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h1 = hg.create_handle("AC", 1);  // len 2
    bdsg::handle_t h2 = hg.create_handle("GT", 2);  // len 2
    hg.create_edge(h1, h2);
    HandleGraphTheseusAdapter adapter(hg);

    int vtx1 = adapter.handle_to_vertex(h1);
    int vtx2 = adapter.handle_to_vertex(h2);

    // ops: M X M M — first two ops belong to node 1, last two to node 2
    // read: A mismatches C (wait: first is M so A matches A; second is X so T mismatches C)
    theseus::Alignment aln = make_aln({vtx1, vtx2},
        {'M','X','M','M'}, 0, 2);
    string sequence = "ATGT";

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.path().mapping_size() == 2);

    const auto& m1 = result.path().mapping(0);
    REQUIRE(m1.edit_size() == 2);
    REQUIRE(m1.edit(0).from_length() == 1);  // M
    REQUIRE(m1.edit(0).to_length() == 1);
    REQUIRE(m1.edit(0).sequence() == "");
    REQUIRE(m1.edit(1).from_length() == 1);  // X
    REQUIRE(m1.edit(1).to_length() == 1);
    REQUIRE(m1.edit(1).sequence() == "T");

    const auto& m2 = result.path().mapping(1);
    REQUIRE(m2.edit_size() == 1);
    REQUIRE(m2.edit(0).from_length() == 2);  // MM merged
    REQUIRE(m2.edit(0).to_length() == 2);
    REQUIRE(m2.edit(0).sequence() == "");
}

TEST_CASE("vg_alignment_from_theseus_alignment: reverse-strand handle sets is_reverse",
          "[theseus_interop]") {
    bdsg::HashGraph hg;
    bdsg::handle_t h = hg.create_handle("ACGT", 1);
    HandleGraphTheseusAdapter adapter(hg);

    bdsg::handle_t h_rev = hg.flip(h);
    int vtx_rev = adapter.handle_to_vertex(h_rev);

    theseus::Alignment aln = make_aln({vtx_rev}, {'M','M','M','M'}, 0, 4);
    string sequence = "ACGT";  // rev-comp of "ACGT" is "ACGT" (self-complementary)

    vg::Alignment result =
        vg_alignment_from_theseus_alignment(aln, sequence, adapter, hg);

    REQUIRE(result.path().mapping_size() == 1);
    const auto& m = result.path().mapping(0);
    REQUIRE(m.position().node_id() == 1);
    REQUIRE(m.position().is_reverse() == true);
}*/

} // namespace unittest
} // namespace vg
