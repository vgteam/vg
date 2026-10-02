/// \file graph_caller.cpp
/// Tests for the VCF output buffer's record ordering.
///
/// Records sharing a position are rare in real output, so the order of ties is checked here.

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "catch.hpp"
#include "../graph_caller.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("The buffered record order is total, so two runs cannot disagree", "[graph_caller]") {
    // Two records at the same position, distinguished only by the snarl they came from. This is the
    // case (contig, POS) cannot order, and it is not rare: a nested site sits at or near its
    // parent's position, and flattening can put two snarls on one anchor base.
    const BufferedRecordKey a{"chr20", 1000, ">1>5"};
    const BufferedRecordKey b{"chr20", 1000, ">6>9"};

    // Antisymmetry: without the ID in the key, the comparator would return false in both directions
    // for this pair, so `std::sort` could leave them in input order.
    REQUIRE(buffered_record_key_less(a, b));
    REQUIRE_FALSE(buffered_record_key_less(b, a));

    // Irreflexivity, on every field, since a strict weak ordering needs it and a `<=` typo here
    // would make std::sort's behaviour undefined rather than merely wrong.
    REQUIRE_FALSE(buffered_record_key_less(a, a));

    SECTION("the earlier field always dominates the later one") {
        // Position beats ID: a later ID at an earlier position still sorts first, so the tie-break
        // cannot reorder the file.
        REQUIRE(buffered_record_key_less(BufferedRecordKey{"chr20", 999, ">9>9"},
                                         BufferedRecordKey{"chr20", 1000, ">1>1"}));
        // Contig beats both.
        REQUIRE(buffered_record_key_less(BufferedRecordKey{"chr1", 5000, ">9>9"},
                                         BufferedRecordKey{"chr20", 1, ">1>1"}));
    }

    SECTION("every input permutation sorts to one output") {
    // Sorting each permutation of a set containing a tie must give the same sequence every time.
        vector<BufferedRecordKey> keys{
            {"chr20", 1000, ">6>9"},
            {"chr20", 1000, ">1>5"},
            {"chr20", 900, ">2>3"},
            {"chr21", 10, ">4>7"},
        };
        sort(keys.begin(), keys.end(), buffered_record_key_less);
        const vector<string> want{">2>3", ">1>5", ">6>9", ">4>7"};

        // Start from the first permutation in id order, so next_permutation walks all of them
        // rather than reporting exhaustion on its first call.
        vector<BufferedRecordKey> permuted = keys;
        sort(permuted.begin(), permuted.end(),
             [](const BufferedRecordKey& x, const BufferedRecordKey& y) { return x.id < y.id; });
        size_t checked = 0;
        do {
            vector<BufferedRecordKey> copy = permuted;
            sort(copy.begin(), copy.end(), buffered_record_key_less);
            for (size_t i = 0; i < want.size(); ++i) {
                REQUIRE(copy[i].id == want[i]);
            }
            ++checked;
        } while (next_permutation(permuted.begin(), permuted.end(),
                                  [](const BufferedRecordKey& x, const BufferedRecordKey& y) {
                                      return x.id < y.id;
                                  }));
        // All 24 permutations of four distinct records, so the claim is exhaustive rather than
        // sampled.
        REQUIRE(checked == 24);
    }
}

TEST_CASE("The GL fold uses the layout its writer actually used", "[graph_caller]") {
    // Two GL orders are in use, and they differ from three alleles up. i-major, which
    // PoissonSupportSnarlCaller writes, is (0,0)(0,1)(0,2)(1,1)(1,2)(2,2). Colexicographic, the VCF
    // specification's order and what ReadLikelihoodSnarlCaller writes, is (0,0)(0,1)(1,1)(0,2)(1,2)
    // (2,2). Indices 2 and 3 swap: (0,2) against (1,1). A fold that assumed the wrong one would
    // transpose two likelihoods, which need not change the best genotype.
    SECTION("the two layouts index the same genotypes differently at n=3") {
        REQUIRE(gl_genotype_index(0, 2, 3, GLLayout::IMajor) == 2);
        REQUIRE(gl_genotype_index(1, 1, 3, GLLayout::IMajor) == 3);
        REQUIRE(gl_genotype_index(1, 1, 3, GLLayout::Colexicographic) == 2);
        REQUIRE(gl_genotype_index(0, 2, 3, GLLayout::Colexicographic) == 3);
        // And agree everywhere else, which is why this went unnoticed.
        for (auto g : {std::make_pair(0u, 0u), std::make_pair(0u, 1u), std::make_pair(1u, 2u),
                       std::make_pair(2u, 2u)}) {
            REQUIRE(gl_genotype_index(g.first, g.second, 3, GLLayout::IMajor)
                    == gl_genotype_index(g.first, g.second, 3, GLLayout::Colexicographic));
        }
    }

    SECTION("folding allele 2 into allele 1 gives a different answer under each layout") {
        // Values chosen so every genotype is distinguishable and the two layouts cannot agree by
        // accident. Read as colexicographic these are:
        //   (0,0)=-9  (0,1)=-8  (1,1)=-7  (0,2)=-2  (1,2)=-1  (2,2)=-6
        const vector<double> gl{-9.0, -8.0, -7.0, -2.0, -1.0, -6.0};
        // Allele 2 is absorbed into allele 1, so the surviving alleles are {0, 1}.
        const vector<int> new_index{0, 1, 1};

        const vector<double> colex = fold_genotype_likelihoods(gl, new_index, 2,
                                                              GLLayout::Colexicographic);
        // New (0,0) takes old (0,0) = -9. New (0,1) takes max of old (0,1) = -8 and (0,2) = -2, so
        // -2. New (1,1) takes max of old (1,1) = -7, (1,2) = -1, (2,2) = -6, so -1.
        REQUIRE(colex.size() == 3);
        REQUIRE(colex[0] == -9.0);
        REQUIRE(colex[1] == -2.0);
        REQUIRE(colex[2] == -1.0);

        // The same input read as i-major puts -7 at (0,2) and -2 at (1,1), so the folded het and hom
        // genotypes take different values, which is wrong for a read-likelihood record.
        const vector<double> imajor = fold_genotype_likelihoods(gl, new_index, 2, GLLayout::IMajor);
        REQUIRE(imajor.size() == 3);
        REQUIRE(imajor[1] != colex[1]);
        // Specifically: the het class loses the -2 it should have absorbed.
        REQUIRE(imajor[1] == -7.0);
    }

    SECTION("a fold that merges nothing is the identity") {
        // The guard against a fold that quietly reorders a record it was not asked to change.
        const vector<double> gl{-5.0, -4.0, -3.0, -2.0, -1.0, -6.0};
        const vector<int> new_index{0, 1, 2};
        for (auto layout : {GLLayout::IMajor, GLLayout::Colexicographic}) {
            const vector<double> same = fold_genotype_likelihoods(gl, new_index, 3, layout);
            REQUIRE(same == gl);
        }
    }
}


/// A traversal over the given node ids, as plain node visits.
static SnarlTraversal make_trav(const vector<nid_t>& nodes) {
    SnarlTraversal t;
    for (nid_t n : nodes) {
        t.add_visit()->set_node_id(n);
    }
    return t;
}

/// A child chain with the given boundary nodes.
static Snarl make_child(nid_t start, nid_t end) {
    Snarl s;
    s.mutable_start()->set_node_id(start);
    s.mutable_end()->set_node_id(end);
    return s;
}

TEST_CASE("offset_of_child reports where a traversal enters a chain", "[graph_caller]") {
    const SnarlTraversal t = make_trav({1, 2, 3, 4, 5});
    SECTION("the entry index, not the exit") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(2, 4)) == 1);
    }
    SECTION("entering from either boundary is the same crossing") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(4, 2)) == 1);
    }
    SECTION("a chain the traversal does not cross has no offset") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(7, 9)) == -1);
    }
    SECTION("touching one boundary only is not a crossing") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(3, 99)) == -1);
    }
}


TEST_CASE("A moved record's quality fields come from the stored direct call", "[graph_caller]") {
    // Settled 1/1, so the settled genotype's GL entry is -3 against a best other of -1: a margin of
    // -20 phred. The printed GQI is capped and the printed GQN is 0, so a divisor recovered from
    // them would be wrong or missing; the stored achievable gap is 100 phred.
    const string line = "chr1\t100\t>1>4\tA\tG\t30\tPASS\t.\tGT:GL:GQ:GQI:GQN\t"
                        "1/1:-5.000000,-1.000000,-3.000000:40:256:0.000";
    auto split = [](const string& text, char delim, vector<string>& out) {
        out.clear();
        size_t start = 0;
        while (true) {
            size_t end = text.find(delim, start);
            out.push_back(text.substr(start, end == string::npos ? string::npos : end - start));
            if (end == string::npos) {
                return;
            }
            start = end + 1;
        }
    };
    auto field = [&](const string& l, const string& key) {
        vector<string> cols, keys, values;
        split(l, '\t', cols);
        split(cols[8], ':', keys);
        split(cols[9], ':', values);
        for (size_t i = 0; i < keys.size(); ++i) {
            if (keys[i] == key) {
                return values[i];
            }
        }
        return string();
    };
    auto filter = [&](const string& l) {
        vector<string> cols;
        split(l, '\t', cols);
        return cols[6];
    };
    LinkageCollector::MovedQuality moved;
    moved.posterior = 0.995;   // -10 log10(0.005) = 23.01
    moved.direct.explained_share = 0.5;
    moved.direct.achievable_gap = 100.0 * log(10.0) / 10.0;

    SECTION("GQ takes the direct call's GQ factor, whatever it is made of") {
        string shared = line, unshared = line;
        moved.direct.gq_factor = 0.5;   // the share, as by default
        REQUIRE(apply_linkage_quality(shared, moved, 0.0));
        REQUIRE(field(shared, "GQ") == "11");
        moved.direct.gq_factor = 1.0;   // --no-share-quality, no depth discount
        REQUIRE(apply_linkage_quality(unshared, moved, 0.0));
        REQUIRE(field(unshared, "GQ") == "23");
    }
    SECTION("GQN divides by the stored achievable gap and multiplies by the share") {
        string rewritten = line;
        REQUIRE(apply_linkage_quality(rewritten, moved, 0.05));
        REQUIRE(field(rewritten, "GQN") == "-0.100");
        REQUIRE(filter(rewritten) == "lowconf");
    }
    SECTION("GQN is missing where the direct call had no achievable gap") {
        string rewritten = line;
        moved.direct.achievable_gap = 0.0;
        REQUIRE(apply_linkage_quality(rewritten, moved, 0.05));
        REQUIRE(field(rewritten, "GQN") == ".");
        REQUIRE(filter(rewritten) == "PASS");
    }
}

/// Exposes the anchor path's strand lookup and the tables it reads.
class StrandLookup : public VCFOutputCaller {
public:
    StrandLookup() : VCFOutputCaller("S") {
        render_lambda_temper = 1.0;
        render_lambda_ceiling = 1.0;
    }
    void add_read(const string& name, size_t phase_set, bool multi) {
        ReadLambda read;
        read.lambda = 2.0;
        read.sites = 1;
        read.phase_set = phase_set;
        read.multi_phase_set = multi;
        render_lambda[(uint64_t)std::hash<string>{}(name)] = read;
    }
    void set_site_phase_set(size_t record_key, size_t phase_set) {
        render_lambda_phase_set[record_key] = phase_set;
    }
    using VCFOutputCaller::read_strand_log_odds;
};

TEST_CASE("A strand from another phase set is NaN in the anchor path, and no strand is 0",
          "[graph_caller]") {
    // --anchors-hom-split drops a NaN read from a split site and places a 0 read by a coin, so the
    // two must stay distinct here; re-genotyping, which does not use this lookup, gives both 0.
    StrandLookup caller;
    caller.add_read("here", 7, false);
    caller.add_read("elsewhere", 9, false);
    caller.add_read("both", 7, true);
    caller.set_site_phase_set(1, 7);

    const double here = caller.read_strand_log_odds(1, "here");
    REQUIRE(std::isfinite(here));
    REQUIRE(here > 0.0);
    REQUIRE(std::isnan(caller.read_strand_log_odds(1, "elsewhere")));
    REQUIRE(std::isnan(caller.read_strand_log_odds(1, "both")));
    REQUIRE(caller.read_strand_log_odds(1, "unseen") == 0.0);

    // A site with no phase set takes any read found in one phase set, but not one found in two.
    REQUIRE(caller.read_strand_log_odds(2, "elsewhere") > 0.0);
    REQUIRE(std::isnan(caller.read_strand_log_odds(2, "both")));
}

TEST_CASE("A parent's phase swap carries its nested strands and their haplotypes with it",
          "[graph_caller]") {
    const size_t W = LinkageModel::WILDCARD;
    // Diploid parent 1, swapped by read phasing. Its ploidy-1 child 2 is on strand 1 with
    // haplotype 7, and 2's own ploidy-1 child 3 is on strand 0 with haplotype 5. Diploid child 4 is
    // not swapped, so its ploidy-1 child 5 keeps its strand.
    auto call = [](size_t key, size_t ploidy, int strand, size_t first, size_t second) {
        LinkageCollector::PhaseCall pc;
        pc.record_key = key;
        pc.ploidy = ploidy;
        pc.nested_strand = (int8_t)strand;
        pc.hap_first = first;
        pc.hap_second = second;
        return pc;
    };
    vector<LinkageCollector::PhaseCall> phased = {
        call(1, 2, -1, 3, 4), call(2, 1, 1, W, 7), call(3, 1, 0, 5, W),
        call(4, 2, -1, 3, 4), call(5, 1, 0, 6, W)};
    std::unordered_map<size_t, size_t> index;
    for (size_t i = 0; i < phased.size(); ++i) {
        index[phased[i].record_key] = i;
    }
    // Out of generation order, as the staged records can be.
    vector<FlowCaller::NestedLink> links = {
        {3, 2, 2}, {5, 4, 2}, {2, 1, 1}, {4, 1, 1}, {1, 0, 0}};
    const unordered_set<size_t> flips = {1};

    SECTION("each strand moves, and its haplotype moves to the slot the strand names") {
        REQUIRE(FlowCaller::cascade_nested_strands(phased, index, links, flips) == 2);
        REQUIRE(phased[1].nested_strand == 0);
        REQUIRE(phased[1].hap_first == 7);
        REQUIRE(phased[1].hap_second == W);
        REQUIRE(phased[2].nested_strand == 1);
        REQUIRE(phased[2].hap_first == W);
        REQUIRE(phased[2].hap_second == 5);
        REQUIRE(phased[4].nested_strand == 0);
        REQUIRE(phased[4].hap_first == 6);
    }
    SECTION("a chain left out of the links stops the swap reaching its children") {
        links.erase(links.begin() + 2);   // chain 2
        REQUIRE(FlowCaller::cascade_nested_strands(phased, index, links, flips) == 0);
        REQUIRE(phased[2].nested_strand == 0);
    }
}

}
}
