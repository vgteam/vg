/// \file multipass_caller.cpp
///
/// Unit tests for MultiPassCaller, run on a SnarlManagerDecomposition and a SnarlDistanceIndex
/// built from the same decomposition of one graph: the caller must write the same records from
/// either.

#include "../../src/unittest/catch.hpp"
#include "../../src/alignment_scorer.hpp"
#include "../allele_likelihood.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../multipass_caller.hpp"
#include "../read_likelihood_caller.hpp"
#include "../site_read_source.hpp"
#include "../../src/traversal_support.hpp"
#include "../../src/sorted_vcf_writer.hpp"
#include "support/decomposition_pair.hpp"

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>

#include <limits>
#include <random>
#include <sstream>

namespace vg {
namespace unittest {

using namespace vg::multipass;

/// An all-match alignment along the given nodes, read forward.
static Alignment read_along(const HandleGraph& graph, const string& name,
                            const vector<nid_t>& nodes) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (nid_t n : nodes) {
        string node_seq = graph.get_sequence(graph.get_handle(n, false));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(n);
        m->mutable_position()->set_is_reverse(false);
        m->mutable_position()->set_offset(0);
        Edit* e = m->add_edit();
        e->set_from_length(node_seq.size());
        e->set_to_length(node_seq.size());
        seq += node_seq;
    }
    aln.set_sequence(seq);
    aln.set_quality(string(seq.size(), (char)30));
    aln.set_mapping_quality(60);
    return aln;
}

/// An all-match alignment of the `length` bases of `walk` from base `from`, read forward.
static Alignment read_on_walk(const HandleGraph& graph, const string& name,
                              const vector<nid_t>& walk, size_t from, size_t length) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    size_t node_start = 0;
    for (nid_t n : walk) {
        const string node_seq = graph.get_sequence(graph.get_handle(n, false));
        const size_t begin = std::max(from, node_start);
        if (const size_t end = std::min(from + length, node_start + node_seq.size());
            begin < end) {
            Mapping* m = aln.mutable_path()->add_mapping();
            m->mutable_position()->set_node_id(n);
            m->mutable_position()->set_offset((int32_t)(begin - node_start));
            Edit* e = m->add_edit();
            e->set_from_length((int32_t)(end - begin));
            e->set_to_length((int32_t)(end - begin));
            seq += node_seq.substr(begin - node_start, end - begin);
        }
        node_start += node_seq.size();
    }
    aln.set_sequence(seq);
    aln.set_quality(string(seq.size(), (char)30));
    aln.set_mapping_quality(60);
    return aln;
}

/// Call `graph`, whose sites `decomposition` gives, from `reads`, with the alleles of its paths,
/// "ref" being the reference, and return the VCF.
static string call_vcf(const PathPositionHandleGraph& graph,
                       const SnarlDecomposition& decomposition, const vector<Alignment>& reads,
                       bool group_candidates = false) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(graph, source, qual_scorer, plain_scorer);
    ReadLikelihoodSnarlCaller genotyper(graph, calculator);
    PathTraversalFinder traversal_finder(graph);
    SortedVCFWriter output("sample");
    MultiPassCaller caller(graph, genotyper, decomposition, output, "sample", traversal_finder,
                           ReferencePathTable({"ref"}, {0}, {2}), false,
                           make_pair((size_t)0, numeric_limits<size_t>::max()));
    caller.set_nested_calling(true);
    caller.set_group_candidates(group_candidates);
    caller.add_header_lines();
    genotyper.add_header_lines(output);
    caller.call(GraphCaller::RecurseOnFail, []() {});
    stringstream vcf;
    output.set_contigs({{"ref", graph.get_path_length(graph.get_path_handle("ref"))}});
    output.write_variants(vcf);
    return vcf.str();
}

/// Return the records of the VCF `text`, without its header.
static vector<string> records_of(const string& text) {
    vector<string> records;
    stringstream lines(text);
    for (string line; getline(lines, line);) {
        if (!line.empty() && line[0] != '#') {
            records.push_back(line);
        }
    }
    return records;
}

TEST_CASE("MultiPassCaller writes the same records on both implementations",
          "[multi_pass_caller]") {
    // Three nested sites, 1-8 holding 2-7 holding 3-5, where the alternative path deletes node
    // 4 inside the innermost.
    bdsg::HashGraph graph;
    const vector<string> sequences{"ACGTACGTACGTACGTACGT", "GATTACAGATTACA", "CCCCAGCCCC", "T",
                                   "GGGGTCGGGG", "A", "TTTTGATTTT", "ACGTTGCAACACGTTGCAAC"};
    for (size_t i = 0; i < sequences.size(); ++i) {
        graph.create_handle(sequences[i], (nid_t)i + 1);
    }
    for (const pair<nid_t, nid_t>& edge : vector<pair<nid_t, nid_t>>{
             {1, 2}, {1, 8}, {2, 3}, {2, 6}, {3, 4}, {3, 5}, {4, 5}, {5, 7}, {6, 7}, {7, 8}}) {
        graph.create_edge(graph.get_handle(edge.first), graph.get_handle(edge.second));
    }
    const vector<nid_t> ref_nodes{1, 2, 3, 4, 5, 7, 8};
    const vector<nid_t> alt_nodes{1, 2, 3, 5, 7, 8};
    for (const auto& path : vector<pair<string, vector<nid_t>>>{{"ref", ref_nodes},
                                                                {"alt", alt_nodes}}) {
        path_handle_t handle = graph.create_path_handle(path.first);
        for (nid_t n : path.second) {
            graph.append_step(handle, graph.get_handle(n));
        }
    }
    bdsg::PositionOverlay positioned(&graph);

    IntegratedSnarlFinder finder(graph);
    DecompositionPair pair(graph, finder);
    REQUIRE(!pair.has_known_difference());

    // A heterozygous deletion of node 4.
    vector<Alignment> reads;
    for (size_t i = 0; i < 12; ++i) {
        reads.push_back(read_along(graph, "ref" + std::to_string(i), ref_nodes));
        reads.push_back(read_along(graph, "alt" + std::to_string(i), alt_nodes));
    }

    const string on_adapter = call_vcf(positioned, pair.adapter, reads);
    const string on_index = call_vcf(positioned, pair.distance_index, reads);
    REQUIRE(on_adapter == on_index);

    // The deletion is called at the innermost site, and its enclosing sites, whose alleles differ
    // only inside it, are written as the reference.
    vector<string> records;
    stringstream lines(on_adapter);
    for (string line; getline(lines, line);) {
        if (!line.empty() && line[0] != '#') {
            records.push_back(line);
        }
    }
    REQUIRE(records.size() == 1);
    REQUIRE(records[0].find("\t>3>5\t") != string::npos);
    REQUIRE(records[0].find("\t0/1:") != string::npos);
}

/**
 * Call site 1-12, which holds the chain 2-11 with the SNP children 2-5 (3 or 4), 5-8 (6 or 7) and
 * 8-11 (9 or 10), while the edge 1-12 deletes the chain. The walks "x1", "x2" and "x3" each carry
 * one alternative base, so no one walk spells a strand that carries all three. One strand of the
 * sample crosses the chain with all three, and the other deletes it.
 *
 * Reads of `read_length` bp, with MAPQ `mapq`, start every `step` bp along each strand.
 *
 * @return The records.
 */
static vector<string> call_mosaic_beside_deletion(size_t read_length, size_t step, int mapq,
                                                  bool group_candidates) {
    bdsg::HashGraph graph;
    std::mt19937 random(7);
    auto bases = [&](size_t length) {
        string seq;
        for (size_t i = 0; i < length; ++i) {
            seq += "ACGT"[random() % 4];
        }
        return seq;
    };
    // The long nodes 1 and 12 let reads cover the deletion's junction as often as a base of the
    // chain.
    const vector<string> sequences{bases(300), bases(25), "T", "A", bases(25), "C", "G",
                                   bases(25), "C", "G", bases(25), bases(300)};
    for (size_t i = 0; i < sequences.size(); ++i) {
        graph.create_handle(sequences[i], (nid_t)i + 1);
    }
    for (const auto& [from, to] : vector<pair<nid_t, nid_t>>{
             {1, 2}, {2, 3}, {2, 4}, {3, 5}, {4, 5}, {5, 6}, {5, 7}, {6, 8}, {7, 8}, {8, 9},
             {8, 10}, {9, 11}, {10, 11}, {11, 12}, {1, 12}}) {
        graph.create_edge(graph.get_handle(from), graph.get_handle(to));
    }
    auto crossing = [](nid_t a, nid_t b, nid_t c) {
        return vector<nid_t>{1, 2, a, 5, b, 8, c, 11, 12};
    };
    const vector<nid_t> deletion{1, 12};
    for (const auto& [name, steps] : vector<pair<string, vector<nid_t>>>{
             {"ref", crossing(3, 6, 9)}, {"x1", crossing(4, 6, 9)}, {"x2", crossing(3, 7, 9)},
             {"x3", crossing(3, 6, 10)}, {"del", deletion}}) {
        path_handle_t handle = graph.create_path_handle(name);
        for (nid_t n : steps) {
            graph.append_step(handle, graph.get_handle(n));
        }
    }
    bdsg::PositionOverlay positioned(&graph);
    IntegratedSnarlFinder finder(graph);
    DecompositionPair pair(graph, finder);

    vector<Alignment> reads;
    for (const vector<nid_t>& strand : {crossing(4, 7, 10), deletion}) {
        size_t length = 0;
        for (nid_t n : strand) {
            length += graph.get_length(graph.get_handle(n));
        }
        for (size_t from = 0; from + read_length <= length; from += step) {
            reads.push_back(read_on_walk(graph, std::to_string(reads.size()), strand, from,
                                         read_length));
            reads.back().set_mapping_quality(mapq);
        }
    }
    const string called = call_vcf(positioned, pair.distance_index, reads, group_candidates);
    INFO(called);
    return records_of(called);
}

/// Require `records` to be the call of `call_mosaic_beside_deletion`: the parent heterozygous for
/// the deletion, and each child, which one strand crosses, written at ploidy 1.
static void require_deletion_and_haploid_children(const vector<string>& records) {
    REQUIRE(records.size() == 4);
    REQUIRE(records[0].find("\t>1>12\t") != string::npos);
    REQUIRE(records[0].find("\t0/1:") != string::npos);
    const vector<string> children{"\t>2>5\t", "\t>5>8\t", "\t>8>11\t"};
    for (size_t i = 0; i < children.size(); ++i) {
        REQUIRE(records[i + 1].find(children[i]) != string::npos);
        REQUIRE(records[i + 1].find("\t1:") != string::npos);
    }
}

TEST_CASE("A child chain that one strand of a grouped parent crosses is written at ploidy 1",
          "[multi_pass_caller]") {
    // The walks through the chain form a group, which explains the crossing strand, so the parent
    // is heterozygous for the deletion. The walks alone call both strands through the chain.
    require_deletion_and_haploid_children(call_mosaic_beside_deletion(60, 3, 60, true));
}

}
}
