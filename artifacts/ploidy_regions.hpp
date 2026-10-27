#ifndef VG_PLOIDY_REGIONS_HPP_INCLUDED
#define VG_PLOIDY_REGIONS_HPP_INCLUDED

#include <cstdint>
#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

namespace vg {

using namespace std;

/**
 * A table of ploidies for regions of reference contigs, read from a BED of `CHROM START END PLOIDY`.
 *
 * It describes a contig whose copy number changes along it, such as the chrX of a male sample,
 * which is haploid except in the pseudoautosomal regions. A position no interval covers keeps the
 * ploidy of its contig.
 *
 * CHROM is the contig as the output VCF names it, which is the locus of a PanSN path name: `chrX`,
 * not `CHM13#0#chrX`. Intervals are 0-based and half-open, as BED has them.
 *
 * A ploidy is 0, 1 or 2. Ploidy 0 means the sample has no copy of the region, such as chrY in a
 * female sample, so a site there is not genotyped and gets no record.
 *
 * No two intervals overlap.
 */
class PloidyRegions {
public:
    /// Make a table with no regions, in which every lookup returns its fallback.
    PloidyRegions() = default;

    /// Read the BED at `bed_path`.
    ///
    /// An unreadable file, a malformed line, a ploidy other than 0, 1 or 2, or two overlapping
    /// intervals end the run with an error.
    explicit PloidyRegions(const string& bed_path);

    /// Read a BED from `bed`, named `source` in error messages.
    ///
    /// A malformed line, a ploidy other than 0, 1 or 2, or two overlapping intervals end the run
    /// with an error.
    PloidyRegions(istream& bed, const string& source);

    /// Return the ploidy at a position on a reference path, or `fallback` where no interval covers
    /// it.
    ///
    /// `position` is a 0-based offset along the contig, as in the BED.
    int region_ploidy(const string& ref_path_name, size_t position, int fallback) const;

    /// Return the ploidy of a site, or `fallback` where no interval covers it.
    ///
    /// `interval_start` is where the first boundary node of the site begins along `ref_path_name`,
    /// as `get_ref_interval` gives it, and `ref_offset` is added to it. That is the POS of the
    /// site's record, less 1, before its alleles are trimmed.
    int ploidy_at(const string& ref_path_name, int64_t interval_start, int64_t ref_offset,
                  int fallback) const;

private:
    /// One BED interval and its ploidy.
    struct Region {
        size_t start;   ///< 0-based, inclusive
        size_t end;     ///< 0-based, exclusive
        int ploidy;
    };
    /// For each contig, its regions, sorted by start.
    unordered_map<string, vector<Region>> by_contig;
};

}

#endif
