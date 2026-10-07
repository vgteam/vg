#ifndef VG_VCF_RECORD_HPP_INCLUDED
#define VG_VCF_RECORD_HPP_INCLUDED

#include <atomic>
#include <iostream>
#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include <gbwt/cached_gbwt.h>
#include "handle.hpp"
#include "linkage_model.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "anchor.hpp"
#include "read_phasing.hpp"
#include "regenotype.hpp"
#include "snarl_caller.hpp"
#include "symbolic_allele.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "vcf_genotype_likelihoods.hpp"

namespace vg {

using namespace std;

/// Special marker value for star alleles in genotype vectors.
/// A star allele (*) represents a haplotype that spans a nested site in the
/// parent but doesn't have a defined traversal at the child level.
constexpr int STAR_ALLELE_MARKER = -2;

/// Special marker value for missing alleles in genotype vectors.
/// Used when a parent allele doesn't traverse a child snarl and star_allele
/// mode is disabled. Outputs as '.' in VCF to maintain consistent ploidy.
constexpr int MISSING_ALLELE_MARKER = -1;

/// The 1-based position on the base path of the base `along_path` bases into `ref_path_name`, a
/// path that may name a subrange of its base path.
static int64_t base_path_position(const string& ref_path_name, int64_t along_path) {
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    const int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    return along_path + 1 + basepath_offset;
}

}

#endif
