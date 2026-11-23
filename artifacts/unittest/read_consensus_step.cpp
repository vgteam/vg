/// \file read_consensus_step.cpp
///
/// Unit tests for the read-consensus step of MultiPassCaller: a snarl gains the consensus of the
/// reads of a strand that disagree with what is called on it, and is called with it.

#include "../../src/unittest/catch.hpp"
#include "../../src/alignment_scorer.hpp"
#include "../allele_likelihood.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../linkage_model.hpp"
#include "../multipass_caller.hpp"
#include "../read_likelihood_caller.hpp"
#include "../site_read_source.hpp"
#include "../../src/traversal_support.hpp"
#include "../../src/sorted_vcf_writer.hpp"
#include "support/decomposition_pair.hpp"

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>
#include <gbwt/dynamic_gbwt.h>

#include <limits>
#include <random>
#include <sstream>

namespace vg {
namespace unittest {

using namespace vg::multipass;

namespace {

/// The bases strand 1 of the sample has between nodes 4 and 7, which no walk spells.
const string NOVEL = "TGCATGCAAGTCGATCGGATTACAGCTAGCTAGGATCCATGCAATGCCGTAAGCTTGCAT";

/// A graph of three snarls in a chain, 1-4, 4-7 and 7-10, each a choice of two single bases.
void make_chain(bdsg::HashGraph& graph) {
    std::mt19937 random(11);
    auto bases = [&](size_t length) {
        string seq;
        for (size_t i = 0; i < length; ++i) {
            seq += "ACGT"[random() % 4];
        }
        return seq;
    };
    const vector<string> sequences{bases(40), "A", "T", bases(30), "C", "G", bases(30), "A",
                                   "C", bases(40)};
    for (size_t i = 0; i < sequences.size(); ++i) {
        graph.create_handle(sequences[i], (nid_t)i + 1);
    }
    for (const auto& [from, to] : vector<pair<nid_t, nid_t>>{
             {1, 2}, {1, 3}, {2, 4}, {3, 4}, {4, 5}, {4, 6}, {5, 7}, {6, 7}, {7, 8}, {7, 9},
             {8, 10}, {9, 10}}) {
        graph.create_edge(graph.get_handle(from), graph.get_handle(to));
    }
}

/// The walks of the panel's four haplotypes.
const vector<vector<nid_t>> PANEL{{1, 2, 4, 5, 7, 8, 10},
                                  {1, 3, 4, 6, 7, 9, 10},
                                  {1, 2, 4, 5, 7, 8, 10},
                                  {1, 3, 4, 5, 7, 9, 10}};

/// A read of the whole of `walk`, matching every node, named `name`. With `inserted` not empty,
/// the read has those bases after node 4 and visits neither node 5 nor node 6.
Alignment read_of(const HandleGraph& graph, const string& name, const vector<nid_t>& walk,
                  const string& inserted = "") {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (nid_t n : walk) {
        if (!inserted.empty() && (n == 5 || n == 6)) {
            continue;
        }
        const string node_seq = graph.get_sequence(graph.get_handle(n, false));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(n);
        Edit* e = m->add_edit();
        e->set_from_length((int32_t)node_seq.size());
        e->set_to_length((int32_t)node_seq.size());
        seq += node_seq;
        if (!inserted.empty() && n == 4) {
            Edit* insertion = m->add_edit();
            insertion->set_to_length((int32_t)inserted.size());
            insertion->set_sequence(inserted);
            seq += inserted;
        }
    }
    aln.set_sequence(seq);
    aln.set_mapping_quality(60);
    return aln;
}

/// What a call of the chain gave.
struct Called {
    /// The records of the VCF, and its header.
    vector<string> records;
    string header;
    ReadConsensusCounts counts;
};

/// Call the chain from `reads` with the linkage model, read phasing and the read-consensus
/// step, which is on where `consensus` is.
Called call_chain(const vector<Alignment>& reads, bool consensus, bool flatten = false) {
    bdsg::HashGraph graph;
    make_chain(graph);
    for (size_t h = 0; h < PANEL.size(); ++h) {
        path_handle_t handle = graph.create_path_handle(h == 0 ? "ref" : "h" + std::to_string(h));
        for (nid_t n : PANEL[h]) {
            graph.append_step(handle, graph.get_handle(n));
        }
    }
    bdsg::PositionOverlay positioned(&graph);
    IntegratedSnarlFinder finder(graph);
    DecompositionPair pair(graph, finder);

    // Each haplotype is stored in both orientations, so GBWT sequences 2h and 2h + 1 are
    // haplotype h.
    gbwt::vector_type text;
    vector<size_t> haplotype_of_sequence;
    for (size_t h = 0; h < PANEL.size(); ++h) {
        for (nid_t id : PANEL[h]) {
            text.push_back(gbwt::Node::encode(id, false));
        }
        text.push_back(gbwt::ENDMARKER);
        for (auto it = PANEL[h].rbegin(); it != PANEL[h].rend(); ++it) {
            text.push_back(gbwt::Node::encode(*it, true));
        }
        text.push_back(gbwt::ENDMARKER);
        haplotype_of_sequence.push_back(h);
        haplotype_of_sequence.push_back(h);
    }
    gbwt::DynamicGBWT dynamic_gbwt;
    dynamic_gbwt.insert(text, true);
    const gbwt::GBWT index(dynamic_gbwt);

    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    ReadConsensusOptions options;
    options.enabled = consensus;
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    AlleleLikelihoodParams scoring;
    scoring.collect_read_phasing = true;
    scoring.span_rank = consensus ? options.min_reads : 0;
    GraphAlignedAlleleLikelihoodCalculator calculator(positioned, source, qual_scorer,
                                                      plain_scorer, scoring);
    ReadLikelihoodSnarlCaller genotyper(positioned, calculator);
    PathTraversalFinder traversal_finder(positioned);
    SortedVCFWriter output("sample");
    MultiPassCaller caller(positioned, genotyper, pair.distance_index, output, "sample",
                           traversal_finder, ReferencePathTable({"ref"}, {0}, {2}), false,
                           make_pair((size_t)0, numeric_limits<size_t>::max()));
    LinkageCollector collector(LinkageModel::Params(), PANEL.size());
    caller.set_linkage(&collector, &index, &haplotype_of_sequence);
    caller.set_emit_phasing(true);
    caller.set_nested_calling(true);
    caller.set_flatten_records(flatten);
    caller.set_read_phasing(true, ReadPhasingParams());
    // The chain has too few reads to fit a temper, so it is given.
    RegenotypeParams regenotyping;
    regenotyping.temper = 1.0;
    caller.set_regenotype(true, regenotyping, 2);
    caller.set_read_consensus(options, scoring.sequence_allele_max_bases,
                              scoring.sequence_allele_max_alignments);
    caller.add_header_lines();
    genotyper.add_header_lines(output);
    // The caller reports to stderr.
    stringstream reports;
    std::streambuf* errors = cerr.rdbuf(reports.rdbuf());
    caller.call(GraphCaller::RecurseOnFail, []() {});
    cerr.rdbuf(errors);
    stringstream vcf;
    output.set_contigs({{"ref", positioned.get_path_length(positioned.get_path_handle("ref"))}});
    output.write_variants(vcf);

    Called called;
    called.counts = caller.read_consensus_counts();
    for (string line; getline(vcf, line);) {
        if (!line.empty() && line[0] != '#') {
            called.records.push_back(line);
        } else {
            called.header += line + "\n";
        }
    }
    return called;
}

/// The reads of a sample whose strand 0 is haplotype 0 and whose strand 1 is haplotype 1, with
/// `inserted` in place of the base between nodes 4 and 7 on strand 1 if it is not empty.
vector<Alignment> sample_reads(const string& inserted) {
    bdsg::HashGraph graph;
    make_chain(graph);
    vector<Alignment> reads;
    for (size_t i = 0; i < 12; ++i) {
        reads.push_back(read_of(graph, "x" + std::to_string(i), PANEL[0]));
        reads.push_back(read_of(graph, "y" + std::to_string(i), PANEL[1], inserted));
    }
    return reads;
}

/// Return field `column` (from 0) of the tab-separated `record`.
string field(const string& record, size_t column) {
    stringstream fields(record);
    string value;
    for (size_t i = 0; i <= column; ++i) {
        getline(fields, value, '\t');
    }
    return value;
}

/// Return the record of `records` whose ID is `id` or starts with `id` and `_`, or "".
string record_of(const vector<string>& records, const string& id) {
    for (const string& record : records) {
        if (field(record, 2) == id || field(record, 2).rfind(id + "_", 0) == 0) {
            return record;
        }
    }
    return "";
}

}

TEST_CASE("A snarl whose reads on one strand carry 60 bases no walk spells gains them as an allele",
          "[read_consensus_step]") {
    const vector<Alignment> reads = sample_reads(NOVEL);

    SECTION("with the step on, the snarl is called with the consensus and marked") {
        const Called called = call_chain(reads, true);
        REQUIRE(called.counts.added == 1);
        REQUIRE(called.counts.sites_with_added == 1);
        REQUIRE(called.counts.sites_called == 1);
        REQUIRE(called.counts.candidate_strands == 1);
        REQUIRE(called.header.find("##INFO=<ID=RC,") != string::npos);
        const string record = record_of(called.records, ">4>7");
        REQUIRE(!record.empty());
        // The 60 bases replace the one base between the bounds on one strand.
        REQUIRE(field(record, 3) == "C");
        REQUIRE(field(record, 4) == NOVEL);
        const string genotype = field(record, 9).substr(0, 3);
        REQUIRE((genotype == "0|1" || genotype == "1|0"));
        REQUIRE((";" + field(record, 7) + ";").find(";RC;") != string::npos);
        // The strand with the 60 bases is the one that takes T at the first snarl.
        const string first = field(record_of(called.records, ">1>4"), 9).substr(0, 3);
        REQUIRE(first == genotype);
        // The other records are not marked.
        REQUIRE(field(record_of(called.records, ">1>4"), 7).find("RC") == string::npos);
        REQUIRE(field(record_of(called.records, ">7>10"), 7).find("RC") == string::npos);
    }

    SECTION("with the step off, the snarl has no such allele and no record is marked") {
        const Called called = call_chain(reads, false);
        REQUIRE(called.counts.sites == 0);
        REQUIRE(called.header.find("ID=RC") == string::npos);
        for (const string& record : called.records) {
            REQUIRE(record.find(NOVEL) == string::npos);
            REQUIRE(field(record, 7).find("RC") == string::npos);
        }
    }

    SECTION("flat records spell the consensus and are marked") {
        const Called called = call_chain(reads, true, true);
        REQUIRE(called.counts.sites_called == 1);
        const string record = record_of(called.records, ">4>7");
        REQUIRE(!record.empty());
        REQUIRE(field(record, 4) == NOVEL);
        REQUIRE((";" + field(record, 7) + ";").find(";RC;") != string::npos);
        REQUIRE((";" + field(record, 7) + ";").find(";FL;") != string::npos);
    }
}

TEST_CASE("A snarl whose reads agree with what is called gains no allele",
          "[read_consensus_step]") {
    const vector<Alignment> reads = sample_reads("");
    const Called with_step = call_chain(reads, true);
    REQUIRE(with_step.counts.sites == 3);
    REQUIRE(with_step.counts.candidate_strands == 0);
    REQUIRE(with_step.counts.added == 0);
    REQUIRE(with_step.counts.sites_called == 0);
    // No read has 15 bases between the bounds of a snarl, and no called allele has, so no
    // snarl's reads are fetched again.
    REQUIRE(with_step.counts.prefiltered == 0);
    REQUIRE(with_step.counts.tested == 0);
    // The records are those of a call without the step.
    REQUIRE(with_step.records == call_chain(reads, false).records);
}

TEST_CASE("The read-consensus prefilter leaves out no snarl the step would give an allele",
          "[read_consensus_step]") {
    // The same sample, with every snarl tested.
    bdsg::HashGraph graph;
    make_chain(graph);
    const vector<Alignment> reads = sample_reads(NOVEL);
    const Called filtered = call_chain(reads, true);
    REQUIRE(filtered.counts.prefiltered == 1);
    REQUIRE(filtered.counts.tested == 1);
    REQUIRE(filtered.counts.candidates_outside_prefilter == 0);
}

}
}
