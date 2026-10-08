#ifndef VG_PLOIDY_REGIONS_HPP_INCLUDED
#define VG_PLOIDY_REGIONS_HPP_INCLUDED

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace vg {

using namespace std;

/**
 * Per-region ploidy overrides, from a BED of `CHROM START END PLOIDY`.
 *
 * `-d` and `--ploidy-regex` set ploidy per contig, which cannot express a contig whose copy
 * number changes along it, such as a male sample's chrX, which is haploid except in the
 * pseudoautosomal regions.
 *
 * The CHROM column matches the contig name as it appears in the output VCF -- the locus part
 * of a PanSN path name, so `chrX` rather than `CHM13#0#chrX`. Intervals are BED half-open and
 * 0-based, and a position no interval covers keeps the contig's ploidy from `-d` or
 * `--ploidy-regex`.
 *
 * Overlapping intervals are an error, since a BED that says two things about one base has no
 * correct reading.
 */
class PloidyRegions {
public:
    /// No regions, so that every lookup returns its fallback.
    PloidyRegions() = default;

    /// Read the BED at `bed_path`. An unreadable file, a malformed line, a ploidy other than 1
    /// or 2, or two overlapping intervals end the run with an error.
    explicit PloidyRegions(const string& bed_path);

    /// Ploidy at this reference position, or `fallback` where no interval covers it. `position` is
    /// a 0-based offset along the contig, as in the BED.
    int region_ploidy(const string& ref_path_name, size_t position, int fallback) const;

    /// region_ploidy for a snarl whose reference interval begins at `interval_start`: the first
    /// base of its first boundary node, which is the record's POS less 1 before the record's
    /// alleles are trimmed. Returns `fallback` when no BED is loaded.
    int ploidy_at(const string& ref_path_name, int64_t interval_start, int64_t ref_offset,
                  int fallback) const;

private:
    struct Region {
        size_t start;   ///< 0-based, inclusive
        size_t end;     ///< 0-based, exclusive
        int ploidy;
    };
    /// Contig name -> its regions, sorted by start and not overlapping.
    unordered_map<string, vector<Region>> by_contig;
};

}

#endif
