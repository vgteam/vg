/// \file read_consensus_step.cpp
///
/// Unit tests for the read-consensus step of MultiPassCaller: a snarl, or a run of neighbouring
/// snarls, gains the consensus of the reads of a strand that disagree with what is called on
/// it, and is called with it.

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

#include <functional>
#include <limits>
#include <map>
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

/// Call the chain that `make` builds, whose panel is `panel`, from `reads` with the linkage
/// model, read phasing and the read-consensus step, which has the settings `options`.
Called call_graph(const function<void(bdsg::HashGraph&)>& make,
                  const vector<vector<nid_t>>& panel, const vector<Alignment>& reads,
                  const ReadConsensusOptions& options, bool flatten = false) {
    const vector<vector<nid_t>>& PANEL = panel;
    const bool consensus = options.enabled;
    bdsg::HashGraph graph;
    make(graph);
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
    caller.set_read_consensus(options);
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

/// Call the chain of `make_chain` from `reads`, with the read-consensus step on where
/// `consensus` is.
Called call_chain(const vector<Alignment>& reads, bool consensus, bool flatten = false) {
    ReadConsensusOptions options;
    options.enabled = consensus;
    return call_graph(make_chain, PANEL, reads, options, flatten);
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
        REQUIRE(called.counts.sites_changed == 1);
        REQUIRE(called.counts.sites_called == 1);
        REQUIRE(called.counts.candidate_strands == 1);
        REQUIRE(called.counts.runs_treated == 0);
        // Only the snarl with the 60 bases has its reads fetched again.
        REQUIRE(called.counts.tested == 1);
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
    REQUIRE(with_step.counts.sites_changed == 0);
    REQUIRE(with_step.counts.sites_called == 0);
    // No read has 15 bases between the bounds of a snarl, and no called allele has, so no
    // snarl's reads are fetched again.
    REQUIRE(with_step.counts.tested == 0);
    // The records are those of a call without the step.
    REQUIRE(with_step.records == call_chain(reads, false).records);
}

namespace {

/// One unit of a tandem repeat. Its fourth base is the one that varies.
const string UNIT = "CTGATCAGGCTTAGCCATGA";
/// The bases of a unit after the one that varies, and those of the next unit before it.
const string REST = UNIT.substr(4) + UNIT.substr(0, 3);

/// Return `length` bases drawn with `random`.
string drawn(std::mt19937& random, size_t length) {
    string seq;
    for (size_t i = 0; i < length; ++i) {
        seq += "ACGT"[random() % 4];
    }
    return seq;
}

/// Build a chain of three snarls, 1-4, 4-7 and 7-10, each a choice of A or G, from the ten
/// node sequences that `flank`, `between` and the two bases give.
void make_three(bdsg::HashGraph& graph, const vector<string>& sequences) {
    for (size_t i = 0; i < sequences.size(); ++i) {
        graph.create_handle(sequences[i], (nid_t)i + 1);
    }
    for (const auto& [from, to] : vector<pair<nid_t, nid_t>>{
             {1, 2}, {1, 3}, {2, 4}, {3, 4}, {4, 5}, {4, 6}, {5, 7}, {6, 7}, {7, 8}, {7, 9},
             {8, 10}, {9, 10}}) {
        graph.create_edge(graph.get_handle(from), graph.get_handle(to));
    }
}

/// A graph of three units of a tandem repeat between two flanks. The base that varies in each
/// unit is a snarl, so the two nodes between the snarls are `REST`, 19 bases of the repeat.
void make_repeat(bdsg::HashGraph& graph) {
    std::mt19937 random(5);
    const string left = drawn(random, 40);
    const string right = drawn(random, 40);
    make_three(graph, {left + UNIT.substr(0, 3), "A", "G", REST, "A", "G", REST, "A", "G",
                       UNIT.substr(4) + right});
}

/// A graph of three snarls with bases of no repeat between them.
void make_unique(bdsg::HashGraph& graph) {
    std::mt19937 random(7);
    vector<string> sequences{drawn(random, 40), "A", "G", drawn(random, 19), "A", "G",
                             drawn(random, 19), "A", "G", drawn(random, 40)};
    make_three(graph, sequences);
}

/// The walks of the panel's four haplotypes: the first takes G at each snarl, the second A.
const vector<vector<nid_t>> THREE_PANEL{{1, 3, 4, 6, 7, 9, 10},
                                        {1, 2, 4, 5, 7, 8, 10},
                                        {1, 3, 4, 6, 7, 9, 10},
                                        {1, 2, 4, 5, 7, 9, 10}};

/// A read of the whole of `walk`, named `name`, with the bases `inserted` gives after the last
/// base of a node, a base that is not the node's at each of the offsets `changed` gives, and
/// the bases `inside` gives before the base at an offset of a node.
Alignment read_with(const HandleGraph& graph, const string& name, const vector<nid_t>& walk,
                    const map<nid_t, string>& inserted,
                    const map<nid_t, vector<size_t>>& changed = {},
                    const map<nid_t, pair<size_t, string>>& inside = {}) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (nid_t n : walk) {
        const string node_seq = graph.get_sequence(graph.get_handle(n, false));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(n);
        const vector<size_t> offsets = changed.contains(n) ? changed.at(n) : vector<size_t>();
        size_t from = 0;
        auto match_to = [&](size_t to) {
            if (to > from) {
                Edit* e = m->add_edit();
                e->set_from_length((int32_t)(to - from));
                e->set_to_length((int32_t)(to - from));
                seq += node_seq.substr(from, to - from);
            }
            from = to;
        };
        if (inside.contains(n)) {
            match_to(inside.at(n).first);
            Edit* insertion = m->add_edit();
            insertion->set_to_length((int32_t)inside.at(n).second.size());
            insertion->set_sequence(inside.at(n).second);
            seq += inside.at(n).second;
        }
        for (size_t offset : offsets) {
            match_to(offset);
            Edit* e = m->add_edit();
            e->set_from_length(1);
            e->set_to_length(1);
            e->set_sequence(node_seq[offset] == 'A' ? "C" : "A");
            seq += e->sequence();
            from = offset + 1;
        }
        match_to(node_seq.size());
        if (inserted.contains(n)) {
            Edit* insertion = m->add_edit();
            insertion->set_to_length((int32_t)inserted.at(n).size());
            insertion->set_sequence(inserted.at(n));
            seq += inserted.at(n);
        }
    }
    aln.set_sequence(seq);
    aln.set_mapping_quality(60);
    return aln;
}

/// Return the sequence strand `strand` of the sample has from node 1 through node 10 of the
/// graph `make` builds, as the records of `records` spell it: those of `sites`, each the ID of
/// a site and the node of its end bound.
string spelled_strand(const function<void(bdsg::HashGraph&)>& make, const vector<string>& records,
                      size_t strand,
                      const vector<pair<string, nid_t>>& sites = {{">1>4", 4}, {">4>7", 7},
                                                                  {">7>10", 10}}) {
    bdsg::HashGraph graph;
    make(graph);
    auto node = [&](nid_t id) { return graph.get_sequence(graph.get_handle(id)); };
    string spelled = node(1);
    for (const auto& [id, after] : sites) {
        const string record = record_of(records, id);
        REQUIRE(!record.empty());
        vector<string> alleles{field(record, 3)};
        stringstream alts(field(record, 4));
        for (string alt; getline(alts, alt, ',');) {
            alleles.push_back(alt);
        }
        const string genotype = field(record, 9).substr(0, 3);
        REQUIRE(genotype[1] == '|');
        spelled += alleles.at((size_t)(genotype[2 * strand] - '0')) + node(after);
    }
    return spelled;
}

}

TEST_CASE("A run of snarls in a tandem repeat gains an expansion its reads place in different "
          "snarls", "[read_consensus_step]") {
    bdsg::HashGraph graph;
    make_repeat(graph);
    // Strand 0 takes G at each snarl. Strand 1 takes A, and has three more units of the
    // repeat, which the mapper puts after node 1, node 4 or node 7: the read is the same.
    string expansion;
    for (size_t i = 0; i < 3; ++i) {
        expansion += "A" + REST;
    }
    vector<Alignment> reads;
    for (size_t i = 0; i < 12; ++i) {
        reads.push_back(read_with(graph, "x" + std::to_string(i), THREE_PANEL[0], {}));
        reads.push_back(read_with(graph, "y" + std::to_string(i), THREE_PANEL[1],
                                  {{(nid_t)(1 + 3 * (i % 3)), expansion}}));
    }
    REQUIRE(reads[1].sequence() == reads[3].sequence());
    REQUIRE(reads[1].sequence() == reads[5].sequence());
    ReadConsensusOptions options;
    options.enabled = true;

    SECTION("no snarl gains an allele on its own, and the run does") {
        const Called called = call_graph(make_repeat, THREE_PANEL, reads, options);
        // Each snarl has four reads of the strand with the expansion, and eight without, so
        // the consensus of its reads is what is called.
        REQUIRE(called.counts.candidate_strands == 3);
        // The run of the first two snarls is treated, and the third snarl then gives it the
        // reads whose expansion the mapper put there. The three become one site.
        REQUIRE(called.counts.runs_treated == 1);
        REQUIRE(called.counts.runs_merged == 1);
        REQUIRE(called.counts.sites_changed == 1);
        REQUIRE(called.counts.sites_called == 1);
        const vector<pair<string, nid_t>> merged{{">1>10", 10}};
        const string first = spelled_strand(make_repeat, called.records, 0, merged);
        const string second = spelled_strand(make_repeat, called.records, 1, merged);
        REQUIRE((first == reads[1].sequence() || second == reads[1].sequence()));
        REQUIRE((first == reads[0].sequence() || second == reads[0].sequence()));
    }
}

/// The sites of the three-snarl chain once its first two snarls are one site.
const vector<pair<string, nid_t>> FIRST_TWO_MERGED{{">1>7", 7}, {">7>10", 10}};

/// Return the sequence strand `strand` of the sample has, as the flat records of `records`
/// spell it: `reference`, with each record's allele for the strand in place of its REF.
string spelled_flat(const string& reference, const vector<string>& records, size_t strand) {
    string spelled;
    size_t copied = 0;
    for (const string& record : records) {
        vector<string> alleles{field(record, 3)};
        stringstream alts(field(record, 4));
        for (string alt; getline(alts, alt, ',');) {
            alleles.push_back(alt);
        }
        const string genotype = field(record, 9).substr(0, 3);
        REQUIRE(genotype[1] == '|');
        const size_t at = std::stoul(field(record, 1)) - 1;
        const string& allele = alleles.at((size_t)(genotype[2 * strand] - '0'));
        if (&allele == &alleles[0]) {
            continue;
        }
        REQUIRE(at >= copied);
        spelled += reference.substr(copied, at - copied) + allele;
        copied = at + alleles[0].size();
    }
    return spelled + reference.substr(copied);
}

/// Require that `called` made one site of the first two snarls of the chain `make` builds, and
/// that its records spell the sequences of the reads `plain` and `changed`, one on each strand.
/// With `reference`, the sequence of the reference path, the records are flat ones.
void require_first_two_merged(const function<void(bdsg::HashGraph&)>& make, const Called& called,
                              const Alignment& plain, const Alignment& changed,
                              const string& reference = "") {
    REQUIRE(called.counts.runs_treated == 1);
    REQUIRE(called.counts.runs_merged == 1);
    // The second snarl is removed, and the one site takes the place of the first.
    REQUIRE(called.counts.sites_changed == 1);
    REQUIRE(called.counts.sites_called == 1);
    REQUIRE(record_of(called.records, ">1>4").empty());
    REQUIRE(record_of(called.records, ">4>7").empty());
    const string merged = record_of(called.records, ">1>7");
    REQUIRE(!merged.empty());
    REQUIRE((";" + field(merged, 7) + ";").find(";RC;") != string::npos);
    const string first = reference.empty()
                             ? spelled_strand(make, called.records, 0, FIRST_TWO_MERGED)
                             : spelled_flat(reference, called.records, 0);
    const string second = reference.empty()
                              ? spelled_strand(make, called.records, 1, FIRST_TWO_MERGED)
                              : spelled_flat(reference, called.records, 1);
    REQUIRE((first == changed.sequence() || second == changed.sequence()));
    REQUIRE((first == plain.sequence() || second == plain.sequence()));
}

TEST_CASE("A run whose consensus changes the node between two snarls becomes one site",
          "[read_consensus_step]") {
    bdsg::HashGraph graph;
    make_unique(graph);
    // Strand 1 has ten more bases in each of the first two snarls, and two bases of node 4,
    // which is between them, changed.
    vector<Alignment> reads;
    for (size_t i = 0; i < 12; ++i) {
        reads.push_back(read_with(graph, "x" + std::to_string(i), THREE_PANEL[0], {}));
        reads.push_back(read_with(graph, "y" + std::to_string(i), THREE_PANEL[1],
                                  {{1, "TTGACCATGC"}, {4, "GGATCTTACA"}}, {{4, {6, 12}}}));
    }
    ReadConsensusOptions options;
    options.enabled = true;
    // The reference path takes the first haplotype's walk.
    for (const bool flatten : {false, true}) {
        const Called called = call_graph(make_unique, THREE_PANEL, reads, options, flatten);
        // No snarl has 15 bases of the difference.
        REQUIRE(called.counts.tested == 0);
        require_first_two_merged(make_unique, called, reads[0], reads[1],
                                 flatten ? reads[0].sequence() : "");
    }
}

TEST_CASE("A run whose reads have bases inside the node between two snarls becomes one site",
          "[read_consensus_step]") {
    const string novel = "TTGACCATGCGGATCTTACATTGACC";
    ReadConsensusOptions options;
    options.enabled = true;

    SECTION("across a node of few bases") {
        bdsg::HashGraph graph;
        make_unique(graph);
        vector<Alignment> reads;
        for (size_t i = 0; i < 12; ++i) {
            reads.push_back(read_with(graph, "x" + std::to_string(i), THREE_PANEL[0], {}));
            reads.push_back(read_with(graph, "y" + std::to_string(i), THREE_PANEL[1], {}, {},
                                      {{4, {9, novel}}}));
        }
        const Called called = call_graph(make_unique, THREE_PANEL, reads, options);
        // No read has a base more between the bounds of any snarl.
        REQUIRE(called.counts.tested == 0);
        require_first_two_merged(make_unique, called, reads[0], reads[1]);
    }

    SECTION("across a node of more bases than two snarls of a run have between them") {
        // The node between the first two snarls has 80 bases.
        auto make_long = [](bdsg::HashGraph& graph) {
            std::mt19937 random(11);
            vector<string> sequences{drawn(random, 40), "A", "G", drawn(random, 80), "A", "G",
                                     drawn(random, 19), "A", "G", drawn(random, 40)};
            make_three(graph, sequences);
        };
        bdsg::HashGraph graph;
        make_long(graph);
        REQUIRE(graph.get_length(graph.get_handle(4)) > options.run_max_junction);
        vector<Alignment> reads;
        for (size_t i = 0; i < 12; ++i) {
            reads.push_back(read_with(graph, "x" + std::to_string(i), THREE_PANEL[0], {}));
            reads.push_back(read_with(graph, "y" + std::to_string(i), THREE_PANEL[1], {}, {},
                                      {{4, {37, novel}}}));
        }
        const Called called = call_graph(make_long, THREE_PANEL, reads, options);
        require_first_two_merged(make_long, called, reads[0], reads[1]);
    }
}

}
}
