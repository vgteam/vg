/** \file
 *
 * Unit tests for recombinator.cpp, specifically for clipping gref fragments to a
 * sampled graph and for partitioning top-level chains whose ends only reference
 * paths reach.
 */

#include "../recombinator.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../utility.hpp"

#include "catch.hpp"

#include <gbwtgraph/gfa.h>

#include <fstream>
#include <vector>

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// gbwt::vector_type may store nodes in a narrower type than gbwt::node_type.
typedef gbwt::vector_type::value_type node_t;
node_t fwd(nid_t id) { return gbwt::Node::encode(id, false); }
node_t rev(nid_t id) { return gbwt::Node::encode(id, true); }

// Node length is the node id, so offsets are easy to predict.
size_t id_as_length(gbwt::node_type node) { return gbwt::Node::id(node); }

void check_piece(const PathPiece& piece, size_t start, size_t nodes, size_t bp_offset, size_t bp_length) {
    REQUIRE(piece.start == start);
    REQUIRE(piece.nodes == nodes);
    REQUIRE(piece.bp_offset == bp_offset);
    REQUIRE(piece.bp_length == bp_length);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Paths are clipped to the nodes and edges of a GBWT index", "[gref][recombinator]") {
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    // Paths 1-2-3, 3-4, 5-6 and 8. Node 7 is inside the alphabet range but unused,
    // and nodes 4 and 5 are both present with no edge between them.
    gbwt::GBWTBuilder builder(gbwt::bit_length(rev(8)), 1024);
    std::vector<gbwt::vector_type> paths {
        { fwd(1), fwd(2), fwd(3) }, { fwd(3), fwd(4) }, { fwd(5), fwd(6) }, { fwd(8) }
    };
    for (const gbwt::vector_type& path : paths) {
        builder.insert(path, true);
    }
    builder.finish();
    const gbwt::DynamicGBWT& index = builder.index;
    REQUIRE(index.contains(fwd(7)));

    SECTION("a path that is entirely in the index is one piece") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(3) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("the reverse orientation is in the index too") {
        gbwt::vector_type path { rev(3), rev(2), rev(1) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("a missing node splits the path and still counts toward the offset") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(7), fwd(3), fwd(4) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 3);
        check_piece(pieces[1], 3, 2, 10, 7);
    }

    SECTION("present nodes without an edge between them are split") {
        gbwt::vector_type path { fwd(3), fwd(4), fwd(5), fwd(6) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 7);
        check_piece(pieces[1], 2, 2, 7, 11);
    }

    SECTION("a node outside the alphabet range is missing") {
        gbwt::vector_type path { fwd(2), fwd(100) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 1, 0, 2);
    }

    SECTION("nothing survives in an empty index") {
        gbwt::DynamicGBWT empty;
        gbwt::vector_type path { fwd(1), fwd(2) };
        REQUIRE(clip_path_to_index(path, empty, id_as_length).empty());
    }
}

//------------------------------------------------------------------------------

namespace {

// The indexes vg haplotypes -H uses, built from GFA lines. As with vg index -P REF, the snarl
// decomposition is rooted at the ends of the first path whose name starts with REF.
struct PartitionerInput {
    gbwtgraph::GBZ gbz;
    gbwt::FastLocate r_index;
    SnarlDistanceIndex distance_index;
    HaplotypePartitioner::minimizer_index_type minimizer_index;

    explicit PartitionerInput(const std::vector<std::string>& gfa_lines) : minimizer_index(7, 3) {
        std::string gfa_file = temp_file::create("recombinator");
        {
            std::ofstream out(gfa_file);
            for (const std::string& line : gfa_lines) {
                out << line << "\n";
            }
        }
        auto parsed = gbwtgraph::gfa_to_gbwt(gfa_file);
        temp_file::remove(gfa_file);
        this->gbz = gbwtgraph::GBZ(parsed.first, parsed.second);
        this->r_index = gbwt::FastLocate(this->gbz.index);

        std::unordered_map<nid_t, size_t> extra_node_weight;
        this->gbz.graph.for_each_path_handle([&](const path_handle_t& path) {
            if (extra_node_weight.empty() && this->gbz.graph.get_path_name(path).compare(0, 3, "REF") == 0) {
                extra_node_weight[this->gbz.graph.get_id(this->gbz.graph.get_handle_of_step(this->gbz.graph.path_begin(path)))] += 10000000000;
                extra_node_weight[this->gbz.graph.get_id(this->gbz.graph.get_handle_of_step(this->gbz.graph.path_back(path)))] += 10000000000;
            }
        });
        IntegratedSnarlFinder finder(this->gbz.graph, extra_node_weight);
        fill_in_distance_index(&this->distance_index, &this->gbz.graph, &finder);

        MinimizerIndexParameters params;
        params.k = 7;
        params.w_or_s = 3;
        this->minimizer_index = build_minimizer_index(this->gbz, nullptr, nullptr, params);
    }

    Haplotypes partition(size_t subchain_length) const {
        HaplotypePartitioner partitioner(this->gbz, this->r_index, this->distance_index, this->minimizer_index, Haplotypes::verbosity_silent);
        HaplotypePartitioner::Parameters parameters;
        parameters.subchain_length = subchain_length;
        return partitioner.partition_haplotypes(parameters);
    }
};

// "type start end sequences" for each subchain of the only top-level chain, in the orientation of the
// reference path, with node ids as boundaries (0 for none). Each toy reference path visits its nodes
// forward, so the chain runs the other way if its boundary nodes are reverse.
std::string summarize(const Haplotypes& haplotypes) {
    typedef std::tuple<char, nid_t, nid_t, size_t> subchain_summary;
    REQUIRE(haplotypes.components() == 1);
    std::vector<subchain_summary> result;
    bool reverse = false;
    for (const Haplotypes::Subchain& subchain : haplotypes.chains.front().subchains) {
        char type = "NPSF"[subchain.type];
        nid_t start = (subchain.start == gbwt::ENDMARKER ? 0 : gbwt::Node::id(subchain.start));
        nid_t end = (subchain.end == gbwt::ENDMARKER ? 0 : gbwt::Node::id(subchain.end));
        reverse = (subchain.start != gbwt::ENDMARKER ? gbwt::Node::is_reverse(subchain.start) :
            (subchain.end != gbwt::ENDMARKER ? gbwt::Node::is_reverse(subchain.end) : reverse));
        result.emplace_back(type, start, end, subchain.num_sequences());
    }
    if (reverse) {
        // Prefixes and suffixes swap places.
        std::reverse(result.begin(), result.end());
        for (auto& summary : result) {
            char type = std::get<0>(summary);
            std::get<0>(summary) = (type == 'P' ? 'S' : (type == 'S' ? 'P' : type));
            std::swap(std::get<1>(summary), std::get<2>(summary));
        }
    }
    std::string text;
    for (auto& summary : result) {
        text += std::string(text.empty() ? "" : ", ") + std::get<0>(summary) + " " + std::to_string(std::get<1>(summary)) + " "
            + std::to_string(std::get<2>(summary)) + " " + std::to_string(std::get<3>(summary));
    }
    return text;
}

// A reference path REF that ends in tips of its own, as vg clip -Sk leaves it: 1 and 11 are its copies of
// the first and last bases of its end nodes, and the haplotypes visit 2 and 10 instead. 12 and 13 hang off
// the haplotypes' copies. In between, 3 6 9 are the chain's nodes, with bubbles 4/5 and 7/8.
std::vector<std::string> ref_tips_gfa(const std::vector<std::string>& extra_walks) {
    std::vector<std::string> result = {
        "H\tVN:Z:1.1\tRS:Z:REF",
        "S\t1\tG", "S\t2\tG", "S\t3\tATTACAGATTACA", "S\t4\tC", "S\t5\tT", "S\t6\tTTAGGTTAGGTTAG",
        "S\t7\tA", "S\t8\tC", "S\t9\tGGGATTTGGGATT", "S\t10\tT", "S\t11\tT", "S\t12\tCCGCCG", "S\t13\tAATAAT",
        "L\t1\t+\t3\t+\t0M", "L\t2\t+\t3\t+\t0M", "L\t12\t+\t2\t+\t0M", "L\t3\t+\t4\t+\t0M", "L\t3\t+\t5\t+\t0M",
        "L\t4\t+\t6\t+\t0M", "L\t5\t+\t6\t+\t0M", "L\t6\t+\t7\t+\t0M", "L\t6\t+\t8\t+\t0M", "L\t7\t+\t9\t+\t0M",
        "L\t8\t+\t9\t+\t0M", "L\t9\t+\t10\t+\t0M", "L\t9\t+\t11\t+\t0M", "L\t10\t+\t13\t+\t0M",
        "W\tREF\t0\tchr1\t0\t42\t>1>3>4>6>7>9>11",
        "W\tH1\t1\tchr1\t0\t53\t>12>2>3>4>6>7>9>10>13",
        "W\tH2\t1\tchr1\t0\t41\t>2>3>5>6>8>9>10",
        "W\tH3\t1\tchr1\t0\t39\t>3>5>6>7>9"
    };
    result.insert(result.end(), extra_walks.begin(), extra_walks.end());
    return result;
}

} // anonymous namespace

TEST_CASE("Chain ends that only reference paths reach become prefixes and suffixes", "[recombinator]") {

    SECTION("a chain whose ends only the reference reaches") {
        PartitionerInput input(ref_tips_gfa({}));
        // The outermost subchains would have the reference as their only candidate. Their candidates are now
        // all sequences at the inner boundaries 3 and 9, which a sampled haplotype follows to its own end.
        std::string expected = "P 0 3 4, N 3 6 4, N 6 9 4, S 9 0 4";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("a chain whose ends a haplotype reaches is unchanged") {
        PartitionerInput input(ref_tips_gfa({ "W\tH4\t1\tchr1\t0\t42\t>1>3>4>6>7>9>11" }));
        std::string expected = "N 1 3 2, N 3 6 5, N 6 9 5, N 9 11 2";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("only the end that no haplotype reaches is opened") {
        PartitionerInput input(ref_tips_gfa({ "W\tH4\t1\tchr1\t0\t29\t>6>7>9>11" }));
        std::string expected = "P 0 3 4, N 3 6 4, N 6 9 5, N 9 11 2";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("generic paths and gref covers count as reference paths") {
        std::vector<std::string> gfa = ref_tips_gfa({ "W\tgref_REF\t0\tchr1\t0\t42\t>1>3>4>6>7>9>11", "P\tx\t11-,9-,7-,6-,4-,3-,1-\t*" });
        PartitionerInput input(gfa);
        std::string expected = "P 0 3 6, N 3 6 6, N 6 9 6, S 9 0 6";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("a single subchain with both ends open is split at the chain node most haplotypes visit") {
        // H4 visits 6 and 9, but not 3, so 6 and 9 have the most haplotype visits, and 6 comes first.
        PartitionerInput input(ref_tips_gfa({ "W\tH4\t1\tchr1\t0\t27\t>6>8>9" }));
        std::string expected = "P 0 6 5, S 6 0 5";
        REQUIRE(summarize(input.partition(10000)) == expected);
    }

    SECTION("a single subchain with one end open becomes a prefix") {
        PartitionerInput input(ref_tips_gfa({ "W\tH4\t1\tchr1\t0\t29\t>6>7>9>11" }));
        std::string expected = "P 0 11 2";
        REQUIRE(summarize(input.partition(10000)) == expected);
    }

    SECTION("an open suffix does not take a sequence that folds back out of it") {
        // H5 turns back at 10 and runs back over the chain. From 9, both of its orientations leave the suffix, so
        // neither is a candidate there, while both cross the normal subchains.
        PartitionerInput input(ref_tips_gfa({ "L\t10\t+\t9\t-\t0M", "W\tH5\t1\tchr1\t0\t86\t>2>3>4>6>7>9>10<9<7<6<4<3" }));
        std::string expected = "P 0 3 6, N 3 6 6, N 6 9 6, S 9 0 4";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("an open prefix does not take a sequence that folds back into it") {
        // H6 comes from 6 back to 2, turns, and goes forward through 3 again.
        PartitionerInput input(ref_tips_gfa({ "L\t2\t-\t3\t+\t0M", "W\tH6\t1\tchr1\t0\t57\t<6<4<3<2>3>5>6" }));
        std::string expected = "P 0 3 4, N 3 6 6, N 6 9 4, S 9 0 4";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    // REF starts with a unary path 1 2 3 that GRC (another reference) doesn't visit; haplotypes start at 6.
    std::vector<std::string> unary_start = {
        "H\tVN:Z:1.1\tRS:Z:REF GRC",
        "S\t1\tGATTACA", "S\t2\tCATCAT", "S\t3\tATTACAGATTACA", "S\t4\tC", "S\t5\tT", "S\t6\tTTAGGTTAGGTTAG",
        "S\t7\tA", "S\t8\tC", "S\t9\tGGGATTTGGGATT",
        "L\t1\t+\t2\t+\t0M", "L\t2\t+\t3\t+\t0M", "L\t3\t+\t4\t+\t0M", "L\t3\t+\t5\t+\t0M", "L\t4\t+\t6\t+\t0M",
        "L\t5\t+\t6\t+\t0M", "L\t6\t+\t7\t+\t0M", "L\t6\t+\t8\t+\t0M", "L\t7\t+\t9\t+\t0M", "L\t8\t+\t9\t+\t0M",
        "W\tREF\t0\tchr1\t0\t55\t>1>2>3>4>6>7>9",
        "W\tGRC\t0\tchr1\t0\t42\t>3>5>6>8>9",
        "W\tH1\t1\tchr1\t0\t28\t>6>8>9",
        "W\tH2\t1\tchr1\t0\t28\t>6>7>9"
    };

    SECTION("an open chain end keeps the reference over the unary path the chain starts with") {
        PartitionerInput input(unary_start);
        std::string expected = "P 0 6 4, N 6 9 4";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("a chain end is not opened when a haplotype reaches the first node of the chain") {
        // HX visits 1 and 2, but not 3, where the first subchain starts.
        std::vector<std::string> gfa = unary_start;
        gfa.push_back("W\tHX\t1\tchr1\t0\t13\t>1>2");
        PartitionerInput input(gfa);
        std::string expected = "N 3 6 2, N 6 9 4";
        REQUIRE(summarize(input.partition(1)) == expected);
    }

    SECTION("a single subchain with both ends open is split at the chain node most haplotype paths visit") {
        // H4 and H7 visit 3, so five haplotype paths visit it, against four for 6. H5 visits 6 three times,
        // turning back at 7 and at 4, which makes more visits to 6 than to 3, but it counts once.
        PartitionerInput input(ref_tips_gfa({
            "L\t7\t+\t8\t-\t0M", "L\t4\t-\t5\t+\t0M",
            "W\tH4\t1\tchr1\t0\t14\t>3>4", "W\tH5\t1\tchr1\t0\t46\t>6>7<8<6<4>5>6", "W\tH7\t1\tchr1\t0\t14\t>3>5"
        }));
        std::string expected = "P 0 3 6, S 3 0 6";
        REQUIRE(summarize(input.partition(10000)) == expected);
    }

    SECTION("a single snarl with both ends open has no chain node to split at") {
        // REF and the generic path x cross the one snarl from 1 to 11, and the haplotype only visits 3 inside it.
        PartitionerInput input({
            "H\tVN:Z:1.1\tRS:Z:REF",
            "S\t1\tG", "S\t2\tATTACAGATTACA", "S\t3\tCCGCCGCCG", "S\t11\tT",
            "L\t1\t+\t2\t+\t0M", "L\t1\t+\t3\t+\t0M", "L\t2\t+\t11\t+\t0M", "L\t3\t+\t11\t+\t0M",
            "W\tREF\t0\tchr1\t0\t15\t>1>2>11",
            "P\tx\t1+,3+,11+\t*",
            "W\tH1\t1\tchr1\t0\t9\t>3"
        });
        std::string expected = "N 1 11 2";
        REQUIRE(summarize(input.partition(10000)) == expected);
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
