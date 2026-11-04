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

/// The order of the genotypes in a diploid `Number=G` field, such as GL.
///
/// The two orders agree for up to two alleles. At three alleles, index 2 is (1,1) in
/// colexicographic order and (0,2) in i-major order.
enum class GLLayout {
    /// Every genotype whose smaller allele is 0, then every one whose smaller allele is 1, and so
    /// on, each run in order of the larger allele.
    IMajor,
    /// The VCF specification's order: by the larger allele, then the smaller.
    Colexicographic,
};

/// Return the index of genotype (i, j), with i <= j, in a `Number=G` field of the given layout.
size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout);

/**
 * Fold a diploid GL vector onto fewer alleles, keeping the best likelihood of genotypes that merge.
 *
 * `new_index[a]` is the new allele that allele `a` becomes. The result has n_new(n_new+1)/2
 * entries, in the same layout. A new genotype that no old genotype maps to gets -infinity.
 */
std::vector<double> fold_genotype_likelihoods(const std::vector<double>& old_gl,
                                              const std::vector<int>& new_index,
                                              size_t n_new, GLLayout layout);

}

#endif
