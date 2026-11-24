#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "ploidy_regions.hpp"
#include "../src/path.hpp"

namespace vg {

namespace {

/**
 * The BED of each preset, by name.
 *
 * A preset gives the ploidy of samples of one karyotype on one reference assembly. Its name is the
 * assembly, a dot, and the karyotype, `XY` or `XX`. The assemblies are:
 *
 * - `chm13v2`: T2T-CHM13v2.0, with contigs `chrX`, `chrY` and `chrM`.
 * - `grch38`: GRCh38, with contigs `chrX`, `chrY` and `chrM`.
 * - `grch37.ucsc`: GRCh37 with UCSC hg19 names, `chrX`, `chrY` and `chrM`.
 * - `grch37.b37`: GRCh37 with b37 and hs37d5 names, `X`, `Y` and `MT`.
 *
 * In an XY preset, chrX is haploid outside its two pseudoautosomal regions (PARs), which are
 * diploid, and chrY and the mitochondrion are haploid. chrY is haploid along its whole length, PARs
 * included, because the second PAR copy of the sample is counted on chrX. In an XX preset, chrX is
 * diploid, chrY has ploidy 0, and the mitochondrion is haploid.
 *
 * A preset has no lines for the autosomes, so a position on an autosome takes the ploidy set for
 * the whole contig, which is 2 by default. One preset therefore serves a whole-genome run and a run
 * on one contig alike.
 *
 * Each BED starts with comment lines that cite the sources of its contig lengths and PAR
 * coordinates.
 */
const map<string, string> PRESET_BEDS = {
    {"chm13v2.XY", R"bed(# Ploidy of an XY sample on T2T-CHM13v2.0: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: an XY sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the CHM13v2.0 assembly (GCA_009914755.4); chrX also matches the CHM13#0#chrX path of
# the HPRC graph. The chrX PARs are as the T2T-Q100 HG002 truth set draws them: it is diploid in
# exactly these two blocks. T2T's own chm13v2.0_PAR.bed puts them at 0-2394410 and
# 153925834-154259566, which differs by 40 bp at PAR1's end and 169 bp at PAR2's start.
chrX	2394370	153926003	1
chrY	0	62460029	1
chrM	0	16569	1
)bed"},
    {"chm13v2.XX", R"bed(# Ploidy of an XX sample on T2T-CHM13v2.0: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the CHM13v2.0 assembly (GCA_009914755.4); chrX also matches the CHM13#0#chrX path of
# the HPRC graph.
chrY	0	62460029	0
chrM	0	16569	1
)bed"},
    {"grch38.XY", R"bed(# Ploidy of an XY sample on GRCh38: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: an XY sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the GRCh38 assembly report (GCA_000001405.15). PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# chrX:10,001-2,781,479 and chrX:155,701,383-156,030,895 (1-based, inclusive).
chrX	0	10000	1
chrX	2781479	155701382	1
chrX	156030895	156040895	1
chrY	0	57227415	1
chrM	0	16569	1
)bed"},
    {"grch38.XX", R"bed(# Ploidy of an XX sample on GRCh38: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh38 assembly report (GCA_000001405.15).
chrY	0	57227415	0
chrM	0	16569	1
)bed"},
    {"grch37.ucsc.XY", R"bed(# Ploidy of an XY sample on GRCh37 with UCSC hg19 names: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: an XY sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1), and UCSC hg19.chrom.sizes for chrM,
# which in hg19 is the 16,571 bp NC_001807 rather than the 16,569 bp rCRS. PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# chrX:60,001-2,699,520 and chrX:154,931,044-155,260,560 (1-based, inclusive).
chrX	0	60000	1
chrX	2699520	154931043	1
chrX	155260560	155270560	1
chrY	0	59373566	1
chrM	0	16571	1
)bed"},
    {"grch37.ucsc.XX", R"bed(# Ploidy of an XX sample on GRCh37 with UCSC hg19 names: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1), and UCSC hg19.chrom.sizes for chrM,
# which in hg19 is the 16,571 bp NC_001807 rather than the 16,569 bp rCRS.
chrY	0	59373566	0
chrM	0	16571	1
)bed"},
    {"grch37.b37.XY", R"bed(# Ploidy of an XY sample on GRCh37 with b37/hs37d5 names: X haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and Y and MT haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Y is haploid along its whole length, PARs included: an XY sample's second PAR copy is
# counted on X, where the reads usually go because analysis sets mask Y's PARs.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1); MT is the 16,569 bp rCRS. PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# X:60,001-2,699,520 and X:154,931,044-155,260,560 (1-based, inclusive).
X	0	60000	1
X	2699520	154931043	1
X	155260560	155270560	1
Y	0	59373566	1
MT	0	16569	1
)bed"},
    {"grch37.b37.XX", R"bed(# Ploidy of an XX sample on GRCh37 with b37/hs37d5 names: X diploid, so it has no line,
# Y absent (ploidy 0, no calls), and MT haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1); MT is the 16,569 bp rCRS.
Y	0	59373566	0
MT	0	16569	1
)bed"},
};

} // namespace

PloidyRegions::PloidyRegions(const string& bed_path) {
    ifstream in(bed_path);
    if (!in) {
        cerr << "error [vg call]: could not open ploidy BED " << bed_path << endl;
        exit(1);
    }
    *this = PloidyRegions(in, "ploidy BED " + bed_path);
}

PloidyRegions::PloidyRegions(istream& bed, const string& source) {
    string line;
    size_t line_number = 0;
    while (getline(bed, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0
            || line.compare(0, 7, "browser") == 0) {
            continue;
        }
        istringstream ss(line);
        string chrom;
        long long start = -1;
        long long end = -1;
        int ploidy = -1;
        if (!(ss >> chrom >> start >> end >> ploidy)) {
            cerr << "error [vg call]: " << source << " line " << line_number
                 << " is not CHROM START END PLOIDY: " << line << endl;
            exit(1);
        }
        if (start < 0 || end < start) {
            cerr << "error [vg call]: " << source << " line " << line_number
                 << " has a negative or reversed interval: " << line << endl;
            exit(1);
        }
        // The callers genotype at ploidy 1 and 2 only, and skip a site at ploidy 0, so reject
        // anything else here.
        if (ploidy < 0 || ploidy > 2) {
            cerr << "error [vg call]: " << source << " line " << line_number
                 << " has ploidy " << ploidy << ", which must be 0, 1 or 2" << endl;
            exit(1);
        }
        if (start == end) {
            // Covers nothing. Keeping it would leave a region no lookup can ever hit.
            continue;
        }
        by_contig[chrom].emplace_back((size_t)start, (size_t)end, ploidy);
    }

    // Sort the regions so lookups can binary-search, and check for overlap while they are in order.
    for (auto& [contig, regions] : by_contig) {
        std::ranges::sort(regions,
                          [](const Region& a, const Region& b) { return a.start < b.start; });
        for (size_t i = 1; i < regions.size(); ++i) {
            if (regions[i].start < regions[i - 1].end) {
                cerr << "error [vg call]: " << source << " has overlapping "
                     << "intervals on " << contig << ": [" << regions[i - 1].start << ","
                     << regions[i - 1].end << ") and [" << regions[i].start << ","
                     << regions[i].end << "). Two ploidies for one base has no correct reading, "
                     << "so this is not resolved by precedence." << endl;
                exit(1);
            }
        }
    }
}

PloidyRegions PloidyRegions::preset(const string& name) {
    auto found = PRESET_BEDS.find(name);
    if (found == PRESET_BEDS.end()) {
        cerr << "error [vg call]: there is no ploidy preset named \"" << name
             << "\". The presets are ";
        vector<string> names = preset_names();
        for (size_t i = 0; i < names.size(); ++i) {
            cerr << (i == 0 ? "" : ", ") << names[i];
        }
        cerr << "." << endl;
        exit(1);
    }
    istringstream bed(found->second);
    return PloidyRegions(bed, "ploidy preset " + name);
}

vector<string> PloidyRegions::preset_names() {
    vector<string> names;
    for (const auto& [name, bed] : PRESET_BEDS) {
        names.push_back(name);
    }
    return names;
}

int PloidyRegions::region_ploidy(const string& ref_path_name, size_t position,
                                 int fallback) const {
    if (by_contig.empty()) {
        return fallback;
    }
    // Match on the contig as the VCF names it, so that a BED written against the output works.
    string contig = Paths::strip_subrange(ref_path_name);
    if (string locus = PathMetadata::parse_locus_name(contig); locus != PathMetadata::NO_LOCUS_NAME) {
        contig = locus;
    }
    auto found = by_contig.find(contig);
    if (found == by_contig.end()) {
        return fallback;
    }
    const vector<Region>& regions = found->second;
    // Find the first region that starts after the position. The regions do not overlap, so only
    // the region before it can cover the position.
    auto it = std::ranges::upper_bound(regions, position, std::less<>{}, &Region::start);
    if (it == regions.begin()) {
        return fallback;
    }
    --it;
    return (position >= it->start && position < it->end) ? it->ploidy : fallback;
}

int PloidyRegions::ploidy_at(const string& ref_path_name, int64_t interval_start,
                             int64_t ref_offset, int fallback) const {
    if (by_contig.empty()) {
        return fallback;
    }
    // Compute where the first boundary node of the site starts along the contig, as a 0-based
    // offset to match the BED. Add the start of any subrange that the path name gives. This is
    // the POS of the site's record less 1, before its alleles are trimmed.
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    int64_t position = interval_start + ref_offset + basepath_offset;
    if (position < 0) {
        return fallback;
    }
    return region_ploidy(ref_path_name, (size_t)position, fallback);
}

}
