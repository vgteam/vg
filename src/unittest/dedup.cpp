/// \file dedup.cpp
///
/// Unit tests for content-aware chopping and packing of piece sequences.

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../../src/handle.hpp"
#include "../../src/unittest/catch.hpp"
#include "../dedup_internal.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// Return a random sequence of bases.
string random_bases(mt19937_64& rng, size_t length) {
    const char bases[] = "ACGT";
    uniform_int_distribution<int> pick(0, 3);
    string result(length, 'A');
    for (char& c : result) {
        c = bases[pick(rng)];
    }
    return result;
}

/// Return the sequence with the base at the offset replaced by another base.
string with_snv(string sequence, size_t offset) {
    sequence[offset] = (sequence[offset] == 'A' ? 'C' : 'A');
    return sequence;
}

/// Return the anchors by their definition: k-mers of bases whose hash is smaller than that of
/// every other such k-mer starting at most `radius` bases away. This checks the window logic of
/// find_anchors() by brute force, and takes the order of two hashes from find_anchors() on a
/// sequence that holds just those two k-mers.
vector<size_t> anchors_by_definition(const string& sequence, size_t k, size_t radius) {
    auto is_kmer = [&](size_t i) {
        return sequence.find_first_not_of("ACGT", i) >= i + k;
    };
    // With a radius covering both k-mers, find_anchors() reports the one with the smaller hash,
    // or neither if they are equal.
    auto smaller = [&](size_t i, size_t j) {
        string pair = sequence.substr(i, k) + "N" + sequence.substr(j, k);
        vector<KmerAnchor> found = find_anchors(pair, k, 2 * k + 1);
        return found.size() == 1 && found[0].offset == 0;
    };
    vector<size_t> result;
    for (size_t i = 0; i + k <= sequence.size(); i++) {
        if (!is_kmer(i)) {
            continue;
        }
        bool anchor = true;
        size_t first = (i > radius ? i - radius : 0);
        for (size_t j = first; j <= i + radius && j + k <= sequence.size() && anchor; j++) {
            if (j != i && is_kmer(j) && !smaller(i, j)) {
                anchor = false;
            }
        }
        if (anchor) {
            result.push_back(i);
        }
    }
    return result;
}

/// Return the offsets of the anchors.
vector<size_t> anchor_offsets(const vector<KmerAnchor>& anchors) {
    vector<size_t> result;
    for (const KmerAnchor& anchor : anchors) {
        result.push_back(anchor.offset);
    }
    return result;
}

/// Check that the plan is a valid chopping of the graph, that the location of each piece holds
/// its sequence, and that the plan has no pieces of nodes that are not chopped.
void check_plan(const bdsg::HashGraph& graph, const ChoppingPlan& plan, const ContentChopParameters& parameters) {
    // No duplicate nodes.
    map<nid_t, vector<size_t>> cuts_by_node(plan.chop_points.begin(), plan.chop_points.end());
    REQUIRE(cuts_by_node.size() == plan.chop_points.size());

    size_t piece_count = 0;
    graph.for_each_handle([&](const handle_t& handle) {
        nid_t node = graph.get_id(handle);
        size_t length = graph.get_length(handle);
        auto iter = cuts_by_node.find(node);
        if (iter == cuts_by_node.end()) {
            // Cut points and pieces are consistent for short nodes.
            REQUIRE(plan.find_piece(node, 0) == nullptr);
            REQUIRE(length <= parameters.max_node_length);
            return;
        }
        REQUIRE(length > parameters.max_node_length);
        REQUIRE(!iter->second.empty());
        vector<size_t> starts { 0 };
        starts.insert(starts.end(), iter->second.begin(), iter->second.end());
        for (size_t i = 0; i < starts.size(); i++) {
            size_t end = (i + 1 < starts.size() ? starts[i + 1] : length);
            // Each piece exists and the information is correct and consistent.
            REQUIRE(starts[i] < end);
            const PieceLocation* piece = plan.find_piece(node, starts[i]);
            REQUIRE(piece != nullptr);
            REQUIRE(piece->length == end - starts[i]);
            REQUIRE(piece->string_id < plan.strings.size());
            REQUIRE(piece->string_offset + piece->length <= plan.strings[piece->string_id].size());
            REQUIRE(string(plan.piece_sequence(*piece)) == graph.get_subsequence(handle, starts[i], piece->length));
            REQUIRE(piece->length <= parameters.max_node_length);
            piece_count++;
        }
    });
    REQUIRE(plan.pieces.size() == piece_count);
}

}

//------------------------------------------------------------------------------

TEST_CASE("Anchors are the k-mers whose hash is the strict minimum within the radius", "[dedup]") {
    mt19937_64 rng(1);
    size_t k = 8;
    size_t radius = 6;

    SECTION("in random bases") {
        string sequence = random_bases(rng, 400);
        REQUIRE(anchor_offsets(find_anchors(sequence, k, radius)) == anchors_by_definition(sequence, k, radius));
    }

    SECTION("with characters that are not bases") {
        string sequence = random_bases(rng, 150) + "NNN" + random_bases(rng, 5) + "N" + random_bases(rng, 150);
        vector<size_t> found = anchor_offsets(find_anchors(sequence, k, radius));
        REQUIRE(found == anchors_by_definition(sequence, k, radius));
        for (size_t offset : found) {
            REQUIRE(sequence.find('N', offset) >= offset + k);
        }
    }

    SECTION("more than the radius apart") {
        string sequence = random_bases(rng, 2000);
        vector<size_t> found = anchor_offsets(find_anchors(sequence, k, radius));
        REQUIRE(found.size() > 10);
        for (size_t i = 1; i < found.size(); i++) {
            REQUIRE(found[i] - found[i - 1] > radius);
        }
    }

    SECTION("lower case bases are the same as upper case bases") {
        string sequence = random_bases(rng, 300);
        string lower = sequence;
        for (char& c : lower) {
            c = tolower(c);
        }
        vector<KmerAnchor> upper_anchors = find_anchors(sequence, k, radius);
        vector<KmerAnchor> lower_anchors = find_anchors(lower, k, radius);
        REQUIRE(anchor_offsets(upper_anchors) == anchor_offsets(lower_anchors));
    }

    SECTION("but none in a homopolymer run or a short tandem repeat") {
        REQUIRE(find_anchors(string(500, 'A'), k, radius).empty());
        string tandem;
        for (size_t i = 0; i < 100; i++) {
            tandem += "CATTC";
        }
        REQUIRE(find_anchors(tandem, k, radius).empty());
    }

    SECTION("with k-mers of 32 bases") {
        string sequence = random_bases(rng, 300);
        REQUIRE(anchor_offsets(find_anchors(sequence, 32, radius)) == anchors_by_definition(sequence, 32, radius));
    }
}

TEST_CASE("Anchors depend only on the nearby bases", "[dedup]") {
    mt19937_64 rng(2);
    size_t k = 16;
    size_t radius = 16;
    string middle = random_bases(rng, 1000);
    string first = random_bases(rng, 300) + middle + random_bases(rng, 200);
    string second = random_bases(rng, 50) + middle + random_bases(rng, 400);

    // The anchors that are at least the radius from both ends of the shared stretch.
    auto anchors_inside = [&](const string& sequence, size_t start) {
        vector<size_t> result;
        for (const KmerAnchor& anchor : find_anchors(sequence, k, radius)) {
            if (anchor.offset >= start + radius && anchor.offset + k + radius <= start + middle.size()) {
                result.push_back(anchor.offset - start);
            }
        }
        return result;
    };
    vector<size_t> in_first = anchors_inside(first, 300);
    REQUIRE(!in_first.empty());
    REQUIRE(in_first == anchors_inside(second, 50));
    REQUIRE(in_first == anchors_inside(middle, 0));
}

TEST_CASE("Substring fingerprints equal the fingerprints of the substrings", "[dedup]") {
    mt19937_64 rng(3);
    string sequence = random_bases(rng, 500) + "NNacgt" + random_bases(rng, 500);
    size_t max_length = 300;
    SubstringFingerprinter fingerprinter(sequence, max_length);
    uniform_int_distribution<size_t> pick_start(0, sequence.size());
    for (size_t trial = 0; trial < 200; trial++) {
        size_t begin = pick_start(rng);
        size_t end = min(sequence.size(), begin + pick_start(rng) % (max_length + 1));
        REQUIRE(fingerprinter.fingerprint(begin, end) == sequence_fingerprint(string_view(sequence).substr(begin, end - begin)));
    }

    SECTION("and distinguish different sequences") {
        string other = with_snv(sequence, 100);
        SubstringFingerprinter other_fingerprinter(other, max_length);
        REQUIRE(fingerprinter.fingerprint(50, 150) != other_fingerprinter.fingerprint(50, 150));
        REQUIRE(fingerprinter.fingerprint(150, 250) == other_fingerprinter.fingerprint(150, 250));
        REQUIRE(sequence_fingerprint("A") != sequence_fingerprint("AA"));
        REQUIRE(sequence_fingerprint("acgt") != sequence_fingerprint("ACGT"));
    }
}

TEST_CASE("Cut candidates cover the sequence with gaps of at most the maximum node length", "[dedup]") {
    mt19937_64 rng(4);
    ContentChopParameters parameters;
    parameters.max_node_length = 100;

    SECTION("in random bases") {
        string sequence = random_bases(rng, 1000);
        CutCandidateSet candidates(sequence, parameters);
        const vector<size_t>& offsets = candidates.offsets();
        REQUIRE(offsets.front() == 0);
        REQUIRE(offsets.back() == sequence.size());
        for (size_t i = 1; i < offsets.size(); i++) {
            REQUIRE(offsets[i] > offsets[i - 1]);
            REQUIRE(offsets[i] - offsets[i - 1] <= parameters.max_node_length);
        }
        // There is about one anchor per 2 * anchor_radius + 1 bases, and fillers break up the
        // longer gaps.
        REQUIRE(offsets.size() > 15);
    }

    SECTION("with fillers where there are no anchors") {
        string sequence = string(250, 'N');
        CutCandidateSet candidates(sequence, parameters);
        REQUIRE(candidates.offsets() == vector<size_t>({ 0, 100, 200, 250 }));
    }

    SECTION("in an empty sequence") {
        CutCandidateSet candidates("", parameters);
        REQUIRE(candidates.offsets() == vector<size_t>({ 0 }));
        size_t segments = 0;
        candidates.for_each_segment([&](size_t, size_t) {
            segments++;
        });
        REQUIRE(segments == 0);
    }

    SECTION("with a segment between every two candidates within the maximum node length") {
        string sequence = random_bases(rng, 500);
        CutCandidateSet candidates(sequence, parameters);
        const vector<size_t>& offsets = candidates.offsets();
        size_t expected = 0;
        for (size_t from = 0; from < offsets.size(); from++) {
            for (size_t to = from + 1; to < offsets.size(); to++) {
                if (offsets[to] - offsets[from] <= parameters.max_node_length) {
                    expected++;
                }
            }
        }
        size_t segments = 0;
        candidates.for_each_segment([&](size_t from, size_t to) {
            REQUIRE(from < to);
            REQUIRE(offsets[to] - offsets[from] <= parameters.max_node_length);
            REQUIRE(candidates.segment_fingerprint(from, to)
                    == sequence_fingerprint(string_view(sequence).substr(offsets[from], offsets[to] - offsets[from])));
            segments++;
        });
        REQUIRE(segments == expected);
    }
}

TEST_CASE("A segment count table counts each fingerprint", "[dedup]") {
    std::vector<std::uint64_t> fingerprints{ 7, 3, 7, 1000, 7, 3 };
    SegmentCountTable counts(fingerprints);
    REQUIRE(counts.size() == 3);
    REQUIRE(counts.count(3) == 2);
    REQUIRE(counts.count(7) == 3);
    REQUIRE(counts.count(1000) == 1);
    REQUIRE(counts.count(5) == 0);

    std::vector<std::uint64_t> empty_fingerprints;
    SegmentCountTable empty(empty_fingerprints);
    REQUIRE(empty.count(5) == 0);
}

TEST_CASE("Chosen cuts minimize the number of pieces in unique sequence", "[dedup]") {
    mt19937_64 rng(5);
    ContentChopParameters parameters;
    string sequence = random_bases(rng, 2500);
    CutCandidateSet candidates(sequence, parameters);
    vector<uint64_t> fingerprints;
    candidates.for_each_segment([&](size_t from, size_t to) {
        fingerprints.push_back(candidates.segment_fingerprint(from, to));
    });
    SegmentCountTable counts(fingerprints);

    vector<size_t> cuts = choose_cuts(candidates, counts, parameters.piece_weight);
    REQUIRE(cuts.size() == 2);
    REQUIRE(cuts[0] > 0);
    REQUIRE(cuts[1] - cuts[0] <= parameters.max_node_length);
    REQUIRE(sequence.size() - cuts[1] <= parameters.max_node_length);
    for (size_t cut : cuts) {
        REQUIRE(std::binary_search(candidates.offsets().begin(), candidates.offsets().end(), cut));
    }

    SECTION("and leave a sequence of at most the maximum node length whole") {
        string short_sequence = sequence.substr(0, parameters.max_node_length);
        CutCandidateSet short_candidates(short_sequence, parameters);
        REQUIRE(choose_cuts(short_candidates, counts, parameters.piece_weight).empty());
    }
}

TEST_CASE("Content-aware chopping leaves short nodes whole and out of the plan", "[dedup]") {
    mt19937_64 rng(6);
    bdsg::HashGraph graph;
    string first = random_bases(rng, 1024);
    graph.create_handle(first, 1);
    graph.create_handle(random_bases(rng, 30), 2);
    graph.create_handle(first, 3);
    graph.create_handle("", 4);
    graph.create_handle(first + random_bases(rng, 500), 5);

    ContentChopParameters parameters;
    ChoppingPlan plan = plan_content_chopping(graph, parameters);
    check_plan(graph, plan, parameters);
    REQUIRE(plan.chop_points.size() == 1);
    REQUIRE(plan.chop_points[0].first == 5);
    for (const PieceLocation& piece : plan.pieces) {
        REQUIRE(piece.node == 5);
    }
    REQUIRE(plan.stored_length() == first.size() + 500);

    SECTION("with no plan if no node is long") {
        bdsg::HashGraph short_graph;
        short_graph.create_handle(first, 1);
        plan = plan_content_chopping(short_graph, parameters);
        REQUIRE(plan.chop_points.empty());
        REQUIRE(plan.strings.empty());
        REQUIRE(plan.pieces.empty());
    }
}

TEST_CASE("Content-aware chopping chops long nodes into pieces of at most the maximum length", "[dedup]") {
    mt19937_64 rng(7);
    bdsg::HashGraph graph;
    graph.create_handle(random_bases(rng, 4500), 1);
    graph.create_handle(random_bases(rng, 1500) + string(3000, 'N') + random_bases(rng, 500), 2);
    graph.create_handle(string(4000, 'A'), 3);

    ContentChopParameters parameters;
    ChoppingPlan plan = plan_content_chopping(graph, parameters);
    check_plan(graph, plan, parameters);
    REQUIRE(plan.chop_points.size() == 3);

    // Unique sequence is chopped into as few pieces as possible.
    REQUIRE(plan.chop_points[0].first == 1);
    REQUIRE(plan.chop_points[0].second.size() == 4);

    // A homopolymer has no anchors, so it is chopped at fillers.
    REQUIRE(plan.chop_points[2].first == 3);
    REQUIRE(plan.chop_points[2].second.size() == 3);
}

TEST_CASE("Content-aware chopping stores repeated sequence once", "[dedup]") {
    mt19937_64 rng(8);
    ContentChopParameters parameters;
    string repeat = random_bases(rng, 3000);

    SECTION("when the copies are identical") {
        bdsg::HashGraph graph;
        graph.create_handle(random_bases(rng, 700) + repeat + random_bases(rng, 900), 1);
        graph.create_handle(random_bases(rng, 1200) + repeat + random_bases(rng, 300), 2);
        graph.create_handle(repeat, 3);

        ChoppingPlan plan = plan_content_chopping(graph, parameters);
        check_plan(graph, plan, parameters);
        // The unique flanks and one copy of the repeat. Pieces that contain the ends of the
        // repeat and a flank can store some of the repeat again.
        size_t unique_length = 700 + 900 + 1200 + 300;
        REQUIRE(plan.stored_length() >= unique_length + repeat.size());
        REQUIRE(plan.stored_length() < unique_length + repeat.size() + 300);
    }

    SECTION("when each copy has a variant of its own") {
        bdsg::HashGraph graph;
        graph.create_handle(repeat, 1);
        graph.create_handle(with_snv(repeat, 800), 2);
        graph.create_handle(with_snv(repeat, 1500), 3);
        graph.create_handle(with_snv(repeat, 2300), 4);

        ChoppingPlan plan = plan_content_chopping(graph, parameters);
        check_plan(graph, plan, parameters);
        // One copy of the repeat and a short piece for each variant.
        REQUIRE(plan.stored_length() < repeat.size() + 3 * 150);

        // Each variant is in a short piece.
        vector<pair<nid_t, size_t>> variants { { 2, 800 }, { 3, 1500 }, { 4, 2300 } };
        for (auto [node, offset] : variants) {
            for (const PieceLocation& piece : plan.pieces) {
                if (piece.node == node && piece.node_offset <= offset && offset < piece.node_offset + piece.length) {
                    REQUIRE(piece.length < 300);
                }
            }
        }
    }
}

TEST_CASE("Content-aware chopping finds short pieces inside stored long pieces", "[dedup]") {
    mt19937_64 rng(9);
    bdsg::HashGraph graph;
    string long_sequence = random_bases(rng, 2000);
    graph.create_handle(long_sequence, 1);
    // Node 1 is stored whole. The sequence of this node is in the middle of it.
    graph.create_handle(long_sequence.substr(500, 1100), 2);

    ContentChopParameters parameters;
    ChoppingPlan plan = plan_content_chopping(graph, parameters);
    check_plan(graph, plan, parameters);
    const PieceLocation* first_piece = plan.find_piece(1, 0);
    REQUIRE(first_piece != nullptr);
    REQUIRE(plan.strings[first_piece->string_id] == long_sequence);

    // A piece with an anchor at least the anchor radius from both of its ends can be found in
    // the stored sequence. A shorter piece is stored again.
    size_t shortest_found = parameters.kmer_length + 2 * parameters.anchor_radius;
    size_t stored_again = 0;
    for (const PieceLocation& piece : plan.pieces) {
        if (piece.node != 2) {
            continue;
        }
        if (piece.length >= shortest_found) {
            REQUIRE(piece.string_id == first_piece->string_id);
            REQUIRE(piece.string_offset == 500 + piece.node_offset);
        } else if (piece.string_id != first_piece->string_id) {
            stored_again += piece.length;
        }
    }
    REQUIRE(plan.stored_length() == long_sequence.size() + stored_again);
}

TEST_CASE("Chopping by a plan makes each piece a node and finds its location", "[dedup]") {
    mt19937_64 rng(10);
    string repeat = random_bases(rng, 2500);
    bdsg::HashGraph graph;
    handle_t a = graph.create_handle(random_bases(rng, 100), 1);
    handle_t b = graph.create_handle(repeat, 2);
    handle_t c = graph.create_handle(with_snv(repeat, 1000), 3);
    handle_t d = graph.create_handle(random_bases(rng, 40), 4);
    graph.create_edge(a, b);
    graph.create_edge(a, c);
    graph.create_edge(b, d);
    graph.create_edge(c, d);
    path_handle_t path = graph.create_path_handle("path");
    for (handle_t handle : { a, b, d }) {
        graph.append_step(path, handle);
    }
    string path_sequence;
    for (handle_t handle : { a, b, d }) {
        path_sequence += graph.get_sequence(handle);
    }
    // Nodes 1 and 4 are short, so they stay whole.
    vector<string> unchopped_sequences { graph.get_sequence(a), graph.get_sequence(d) };

    ContentChopParameters parameters;
    ChoppingPlan plan = plan_content_chopping(graph, parameters);
    check_plan(graph, plan, parameters);
    REQUIRE(!plan.chop_points.empty());

    PieceDirectory directory = chop_by_plan(graph, plan);
    REQUIRE(directory.size() == plan.pieces.size());

    // Each node of the chopped graph is an unchopped node or a piece of the plan. The directory
    // finds each piece by its new node ID.
    set<pair<nid_t, size_t>> pieces_found;
    vector<string> sequences_not_found;
    graph.for_each_handle([&](const handle_t& handle) {
        nid_t id = graph.get_id(handle);
        const PieceLocation* piece = directory.find(id);
        if (piece == nullptr) {
            sequences_not_found.push_back(graph.get_sequence(handle));
            return;
        }
        // The directory holds the location from the plan.
        const PieceLocation* planned = plan.find_piece(piece->node, piece->node_offset);
        REQUIRE(planned != nullptr);
        REQUIRE(piece->string_id == planned->string_id);
        REQUIRE(piece->string_offset == planned->string_offset);
        REQUIRE(piece->length == planned->length);
        // The location describes the node.
        REQUIRE(piece->length == graph.get_length(handle));
        REQUIRE(string(plan.piece_sequence(*piece)) == graph.get_sequence(handle));
        REQUIRE(graph.get_length(handle) <= parameters.max_node_length);
        // Each piece is a single node.
        REQUIRE(pieces_found.emplace(piece->node, piece->node_offset).second);
    });
    REQUIRE(pieces_found.size() == plan.pieces.size());
    sort(sequences_not_found.begin(), sequences_not_found.end());
    sort(unchopped_sequences.begin(), unchopped_sequences.end());
    REQUIRE(sequences_not_found == unchopped_sequences);

    // The directory has nothing for an ID beyond the largest in the graph.
    REQUIRE(directory.find(graph.max_node_id() + 1) == nullptr);

    // Chopping keeps the path.
    string chopped_path_sequence;
    graph.for_each_step_in_path(graph.get_path_handle("path"), [&](const step_handle_t& step) {
        chopped_path_sequence += graph.get_sequence(graph.get_handle_of_step(step));
    });
    REQUIRE(chopped_path_sequence == path_sequence);

    SECTION("but leaves the graph as it is with an empty plan") {
        bdsg::HashGraph short_graph;
        short_graph.create_handle(random_bases(rng, 100), 1);
        ChoppingPlan empty_plan = plan_content_chopping(short_graph, parameters);
        REQUIRE(empty_plan.chop_points.empty());
        PieceDirectory empty_directory = chop_by_plan(short_graph, empty_plan);
        REQUIRE(empty_directory.size() == 0);
        REQUIRE(short_graph.get_node_count() == 1);
        REQUIRE(short_graph.get_length(short_graph.get_handle(1)) == 100);
    }
}

TEST_CASE("Content-aware chopping rejects unusable settings", "[dedup]") {
    bdsg::HashGraph graph;
    graph.create_handle("ACGT", 1);

    SECTION("a maximum node length of 0") {
        ContentChopParameters parameters;
        parameters.max_node_length = 0;
        REQUIRE_THROWS_AS(plan_content_chopping(graph, parameters), std::invalid_argument);
    }

    SECTION("a negative piece weight") {
        ContentChopParameters parameters;
        parameters.piece_weight = -1.0;
        REQUIRE_THROWS_AS(plan_content_chopping(graph, parameters), std::invalid_argument);
    }

    SECTION("k-mers that do not fit in 64 bits") {
        ContentChopParameters parameters;
        parameters.kmer_length = 33;
        REQUIRE_THROWS_AS(plan_content_chopping(graph, parameters), std::invalid_argument);
    }

    SECTION("k-mers of length 0") {
        ContentChopParameters parameters;
        parameters.kmer_length = 0;
        REQUIRE_THROWS_AS(plan_content_chopping(graph, parameters), std::invalid_argument);
    }
}

}
}
