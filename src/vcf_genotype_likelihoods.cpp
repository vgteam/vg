#include "vcf_genotype_likelihoods.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace vg {

using namespace std;

size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout) {
    if (layout == GLLayout::Colexicographic) {
        // The VCF specification's order: genotypes sorted by their larger allele, then their
        // smaller.
        return j * (j + 1) / 2 + i;
    }
    // i-major: all genotypes with smaller allele 0, then all with 1, and so on.
    return i * n_alleles - (i * (i - 1)) / 2 + (j - i);
}

vector<double> fold_genotype_likelihoods(const vector<double>& old_gl,
                                         const vector<int>& new_index,
                                         size_t n_new, GLLayout layout) {
    const size_t n_old = new_index.size();
    vector<double> folded(n_new * (n_new + 1) / 2, -std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < n_old; ++i) {
        for (size_t j = i; j < n_old; ++j) {
            const size_t from = gl_genotype_index(i, j, n_old, layout);
            if (from >= old_gl.size()) {
                continue;
            }
            size_t ni = (size_t)new_index[i], nj = (size_t)new_index[j];
            if (ni > nj) {
                std::swap(ni, nj);
            }
            double& slot = folded[gl_genotype_index(ni, nj, n_new, layout)];
            slot = std::max(slot, old_gl[from]);
        }
    }
    return folded;
}

}
