/// \file off_panel_svs.cpp
///
/// Unit tests for finding large off-panel indels. An off-panel indel is an insertion or deletion
/// that recurs in the reads but that no walk of the graph spells. A large one is 50 bases or more
/// long.
///
/// Reads show a large indel in two ways. A long read holds the whole indel, but sequencing errors
/// shift its position and length a little from read to read. A short read that crosses one end of
/// a deletion is aligned up to that end, and its remaining bases are left unaligned as a clip. The
/// clip's bases match the graph past the other end of the deletion.
///
/// The tests check three things. Long reads' slightly different insertions are counted as one
/// candidate. A deletion is found from short reads' clips at both its ends, and not from clips at
/// one end alone. A single read is classed against a known large deletion as carrying it, as
/// carrying the reference, as carrying another indel, or as saying nothing about it.

#include <string>
#include <tuple>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../../src/unittest/catch.hpp"
#include "../off_panel.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A fixed sequence with no long repeats.
static string sv_sequence(size_t length, uint32_t seed) {
    string out;
    for (size_t i = 0; i < length; ++i) {
        seed = seed * 1103515245u + 12345u;
        out.push_back("ACGT"[(seed >> 16) & 3]);
    }
    return out;
}

/// An edit: from length, to length, and the read's bases where they differ.
using SvEdit = tuple<int, int, string>;

/// A read on node 1 from `offset`, with `edits`, its bases those of `node` where it matches.
static Alignment sv_read(const string& name, const string& node, size_t offset,
                         const vector<SvEdit>& edits) {
    Alignment aln;
    aln.set_name(name);
    aln.set_mapping_quality(60);
    Mapping* m = aln.mutable_path()->add_mapping();
    m->mutable_position()->set_node_id(1);
    m->mutable_position()->set_offset(offset);
    string sequence;
    size_t at = offset;
    for (const SvEdit& spec : edits) {
        Edit* e = m->add_edit();
        e->set_from_length(get<0>(spec));
        e->set_to_length(get<1>(spec));
        if (get<0>(spec) == get<1>(spec)) {
            sequence += node.substr(at, get<0>(spec));
        } else {
            e->set_sequence(get<2>(spec));
            sequence += get<2>(spec);
        }
        at += get<0>(spec);
    }
    aln.set_sequence(sequence);
    return aln;
}

static vector<EditCandidate> count_reads(const bdsg::HashGraph& graph,
                                         const vector<Alignment>& reads) {
    EditCountParams params;
    params.min_base_quality = 0;
    params.large_indels = true;
    return count_edits(
        graph,
        [&](const function<void(const Alignment&)>& visit) {
            for (const Alignment& aln : reads) {
                visit(aln);
            }
        },
        [](nid_t) { return true; }, params);
}

TEST_CASE("Long reads' loosely matching large insertions are one candidate", "[off_panel]") {
    const string node = sv_sequence(400, 7);
    bdsg::HashGraph graph;
    graph.create_handle(node, 1);
    const string inserted = sv_sequence(60, 99);
    vector<Alignment> reads;
    // Eight reads with the insertion, moved up to 3 bases and sized 57 to 63.
    for (int r = 0; r < 8; ++r) {
        const int shift = r % 4 - 1;
        const int length = 57 + (r * 5) % 7;
        string ins = inserted.substr(0, std::min(60, length));
        while ((int)ins.size() < length) {
            ins.push_back('A');
        }
        reads.push_back(sv_read("alt" + to_string(r), node, 20,
                                {{180 + shift, 180 + shift, ""},
                                 {0, length, ins},
                                 {150 - shift, 150 - shift, ""}}));
    }
    for (int r = 0; r < 4; ++r) {
        reads.push_back(sv_read("ref" + to_string(r), node, 30, {{320, 320, ""}}));
    }
    vector<EditCandidate> found = count_reads(graph, reads);
    REQUIRE(found.size() == 1);
    const EditCandidate& c = found.front();
    REQUIRE(c.kind == 'I');
    REQUIRE(c.is_large());
    REQUIRE(c.alt_fragments == 8);
    REQUIRE(c.ref_fragments == 4);
    REQUIRE(c.offset >= 197);
    REQUIRE(c.offset <= 201);
    REQUIRE(c.sequence.size() >= 57);
    REQUIRE(c.sequence.size() <= 63);
    REQUIRE(c.walk == vector<nid_t>{1});
}

TEST_CASE("A deletion is found from short reads' clips at both its ends", "[off_panel]") {
    const string node = sv_sequence(400, 11);
    bdsg::HashGraph graph;
    graph.create_handle(node, 1);
    // An 80-base deletion of bases 150 to 229.
    vector<Alignment> reads;
    for (int r = 0; r < 5; ++r) {
        // Aligned up to base 149, then the bases after the deletion, clipped.
        const size_t start = 40 + r * 10;
        reads.push_back(sv_read("left" + to_string(r), node, start,
                                {{(int)(150 - start), (int)(150 - start), ""},
                                 {0, 30, node.substr(230, 30)}}));
        // The bases before the deletion, clipped, then aligned from base 230.
        reads.push_back(sv_read("right" + to_string(r), node, 230,
                                {{0, 30, node.substr(120, 30)},
                                 {60 + r * 10, 60 + r * 10, ""}}));
    }
    for (int r = 0; r < 3; ++r) {
        reads.push_back(sv_read("ref" + to_string(r), node, 100, {{200, 200, ""}}));
    }
    // A clip that matches nothing along the graph finds nothing.
    for (int r = 0; r < 4; ++r) {
        reads.push_back(sv_read("junk" + to_string(r), node, 300,
                                {{50, 50, ""}, {0, 30, string(30, 'T')}}));
    }
    vector<EditCandidate> found = count_reads(graph, reads);
    REQUIRE(found.size() == 1);
    const EditCandidate& c = found.front();
    REQUIRE(c.kind == 'D');
    REQUIRE(c.offset == 149);
    REQUIRE(c.sequence == node.substr(150, 80));
    REQUIRE(c.alt_fragments == 10);
    REQUIRE(c.ref_fragments == 3);

    // Clips at one end alone are not enough.
    vector<Alignment> one_end;
    for (const Alignment& aln : reads) {
        if (aln.name().compare(0, 5, "right") != 0) {
            one_end.push_back(aln);
        }
    }
    REQUIRE(count_reads(graph, one_end).empty());
}

TEST_CASE("A large deletion's reads fall into their classes", "[off_panel]") {
    WalkIndel deletion;
    deletion.at = 100;
    deletion.deleted = 60;
    const int64_t after = 160;
    LargeIndelParams params;
    auto alignment = [](int64_t begin, int64_t end, vector<WalkIndel> indels, size_t clip_begin,
                        size_t clip_end) {
        WalkAlignment aln;
        aln.begin = begin;
        aln.end = end;
        aln.indels = std::move(indels);
        aln.clip_begin = clip_begin;
        aln.clip_end = clip_end;
        return aln;
    };
    WalkIndel shifted;
    shifted.at = 104;
    shifted.deleted = 55;
    int clip = 0;
    WalkAlignment aln = alignment(0, 300, {shifted}, 0, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end, &clip)
            == IndelObservation::ALT);
    REQUIRE(clip == 0);
    aln = alignment(20, 101, {}, 0, 40);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end, &clip)
            == IndelObservation::ALT);
    REQUIRE(clip == 1);
    aln = alignment(158, 300, {}, 40, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end, &clip)
            == IndelObservation::ALT);
    REQUIRE(clip == 2);
    aln = alignment(50, 140, {}, 0, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end)
            == IndelObservation::REF);
    // Past the walk's start, where the read carries on.
    aln = alignment(0, 140, {}, 0, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, -1000, aln.end)
            == IndelObservation::REF);
    aln = alignment(95, 110, {}, 0, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end)
            == IndelObservation::NONE);
    WalkIndel other;
    other.at = 110;
    other.deleted = 25;
    aln = alignment(0, 300, {other}, 0, 0);
    REQUIRE(observe_large_indel(aln, deletion, after, params, aln.begin, aln.end)
            == IndelObservation::OTHER);
}

}
}
