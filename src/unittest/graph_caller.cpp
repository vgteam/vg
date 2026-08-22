/// \file graph_caller.cpp
/// Tests for the VCF output buffer's record ordering.
///
/// The property here is only checkable in a unit test. `test/t/18_vg_call.t` already asserts that
/// output is independent of thread count and already passed before this ordering existed, which is
/// direct evidence that the in-tree fixtures produce no records sharing a position -- so the TAP
/// suite cannot see the defect at all. On chr20 it showed up as 72 record pairs swapping between two
/// runs of one binary.

#include <algorithm>
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

    // Antisymmetry, which is what the old comparator lacked: it returned false in both directions
    // for this pair, so `std::sort` was free to leave them in whichever order the per-thread buffers
    // happened to concatenate in -- a function of the input permutation rather than of the data.
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
        // The actual guarantee the gates need. Sorting each permutation of a set containing a tie
        // must give the same sequence every time; under the old comparator the two tied records came
        // out in input order, so this produced two distinct answers.
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
    // Two GL orders are live in this tree and they diverge from three alleles up. i-major, which
    // PoissonSupportSnarlCaller writes, is (0,0)(0,1)(0,2)(1,1)(1,2)(2,2). Colexicographic, the VCF
    // spec's order and what ReadLikelihoodSnarlCaller writes, is (0,0)(0,1)(1,1)(0,2)(1,2)(2,2).
    // Indices 2 and 3 swap: (0,2) against (1,1).
    //
    // The allele-merge fold assumed i-major, in a comment naming the Poisson caller as GL's only
    // writer. It is not. So under -L --read-likelihood a merged record had two of its six
    // likelihoods transposed, folding the het-with-allele-2 class into the hom-allele-1 class and
    // back. Invisible in output: on chr20, 0 of 57 merged records violate the "no genotype is
    // strictly more likely than the called one" invariant, because a transposition need not move
    // the argmax.
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
        // classes take different values. This is the wrong answer for a read-likelihood record, and
        // it is what the code produced before the layout was passed in.
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

}
}
