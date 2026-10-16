#ifndef VG_VCF_GENOTYPE_LIKELIHOODS_HPP_INCLUDED
#define VG_VCF_GENOTYPE_LIKELIHOODS_HPP_INCLUDED

/** \file vcf_genotype_likelihoods.hpp
 *
 * The order of the genotypes in a VCF `Number=G` field, such as GL, and folding such a field onto
 * fewer alleles. Two orders are in use (see `GLLayout`).
 */

#include <cstddef>
#include <vector>

namespace vg {

/// The order of the genotypes in a diploid `Number=G` vector such as GL, over n alleles. The two
/// orders agree up to two alleles and differ from three up: at n=3, index 2 is (1,1) in
/// colexicographic order and (0,2) in i-major order, so code that reindexes a GL vector must know
/// which it holds.
enum class GLLayout {
    /// For i = 0..n-1, for j = i..n-1: every genotype whose smaller allele is 0, then 1, and so on.
    IMajor,
    /// The VCF specification's order: by the larger allele, then the smaller.
    Colexicographic,
};

/// Index of genotype (i, j), i <= j, in a `Number=G` vector of the given layout.
size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout);

/// Max-marginal fold of a diploid GL vector onto fewer alleles. `new_index[a]` is the allele `a`
/// becomes. Genotypes that map to one new genotype merge, taking the best likelihood among them.
/// Returns n_new(n_new+1)/2 entries in the same layout; a new genotype that no old one maps to
/// gets -infinity.
std::vector<double> fold_genotype_likelihoods(const std::vector<double>& old_gl,
                                              const std::vector<int>& new_index,
                                              size_t n_new, GLLayout layout);

}

#endif
