#ifndef VG_PLOIDY_PRESETS_HPP_INCLUDED
#define VG_PLOIDY_PRESETS_HPP_INCLUDED

/** \file
 * Ploidy tables for common reference assemblies and sample sexes, built into vg.
 */

#include <map>
#include <string>

namespace vg {

using namespace std;

/**
 * Get the built-in ploidy tables, by name.
 *
 * Each table is a BED of `CHROM START END PLOIDY`, for `PloidyRegions` to read. A name is an
 * assembly and a sample's sex chromosomes, such as `chm13v2.XY`. In an XY table, chrX is haploid
 * outside its two pseudoautosomal regions (PARs), which are diploid, and chrY and the mitochondrion
 * are haploid. In an XX table, chrX is diploid, chrY has ploidy 0, and the mitochondrion is haploid.
 * The autosomes have no lines, so they keep the ploidy they would otherwise have.
 */
const map<string, string>& ploidy_presets();

}

#endif
