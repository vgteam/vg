/** \file
 *
 * Unit tests for vcf_genotype_likelihoods.cpp: the two orders of a `Number=G` field, and folding
 * one onto fewer alleles.
 */

#include <utility>
#include <vector>

#include "../vcf_genotype_likelihoods.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("The GL fold uses the layout its writer actually used", "[vcf_genotype_likelihoods]") {
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

}
}
