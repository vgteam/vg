/// \file multipass_caller.cpp
///
/// Unit tests for MultiPassCaller, run on a SnarlManagerDecomposition and a SnarlDistanceIndex
/// built from the same decomposition of one graph: the caller must write the same records from
/// either.

#include "catch.hpp"
#include "../alignment_scorer.hpp"
#include "../allele_likelihood.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../multipass_caller.hpp"
#include "../read_likelihood_caller.hpp"
#include "../site_read_source.hpp"
#include "../traversal_support.hpp"
#include "../vcf_output_caller.hpp"
#include "support/decomposition_pair.hpp"

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>

#include <limits>
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

/// Call `graph`, whose sites `decomposition` gives, from `reads`, with the alleles of its paths,
/// "ref" being the reference, and return the VCF.
static string call_vcf(const PathPositionHandleGraph& graph, SnarlManager& manager,
                       const SnarlDecomposition& decomposition, const vector<Alignment>& reads) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(graph, source, qual_scorer, plain_scorer);
    NullTraversalSupportFinder support(graph, manager);
    ReadLikelihoodSnarlCaller genotyper(graph, manager, support, calculator);
    PathTraversalFinder traversal_finder(graph);
    VCFOutputCaller output("sample");
    MultiPassCaller caller(graph, genotyper, decomposition, output, "sample", traversal_finder,
                           {"ref"}, {0}, {2}, false,
                           make_pair((size_t)0, numeric_limits<size_t>::max()), false, false);
    caller.set_nested_calling(true);
    stringstream vcf;
    vcf << caller.vcf_header(graph, {"ref"}, {});
    caller.call(GraphCaller::RecurseOnFail, []() {});
    output.write_variants(vcf);
    return vcf.str();
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

    const string on_adapter = call_vcf(positioned, pair.manager, pair.adapter, reads);
    const string on_index = call_vcf(positioned, pair.manager, pair.distance_index, reads);
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

}
}
