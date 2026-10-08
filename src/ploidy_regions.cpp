#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "ploidy_regions.hpp"
#include "path.hpp"

namespace vg {

PloidyRegions::PloidyRegions(const string& bed_path) {
    ifstream in(bed_path);
    if (!in) {
        cerr << "error [vg call]: could not open --ploidy-bed file " << bed_path << endl;
        exit(1);
    }
    string line;
    size_t line_number = 0;
    while (getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0
            || line.compare(0, 7, "browser") == 0) {
            continue;
        }
        istringstream ss(line);
        string chrom;
        long long start = -1, end = -1;
        int ploidy = -1;
        if (!(ss >> chrom >> start >> end >> ploidy)) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " is not CHROM START END PLOIDY: " << line << endl;
            exit(1);
        }
        if (start < 0 || end < start) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " has a negative or reversed interval: " << line << endl;
            exit(1);
        }
        // The callers support ploidy 1 and 2 only, so reject anything else here.
        if (ploidy != 1 && ploidy != 2) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " has ploidy " << ploidy << ", which must be 1 or 2" << endl;
            exit(1);
        }
        if (start == end) {
            // Covers nothing. Keeping it would leave a region no lookup can ever hit.
            continue;
        }
        by_contig[chrom].push_back({(size_t)start, (size_t)end, ploidy});
    }

    // Sorted so lookups can binary-search, and checked for overlap while they are in order.
    for (auto& entry : by_contig) {
        auto& regions = entry.second;
        sort(regions.begin(), regions.end(),
             [](const Region& a, const Region& b) { return a.start < b.start; });
        for (size_t i = 1; i < regions.size(); ++i) {
            if (regions[i].start < regions[i - 1].end) {
                cerr << "error [vg call]: --ploidy-bed " << bed_path << " has overlapping "
                     << "intervals on " << entry.first << ": [" << regions[i - 1].start << ","
                     << regions[i - 1].end << ") and [" << regions[i].start << ","
                     << regions[i].end << "). Two ploidies for one base has no correct reading, "
                     << "so this is not resolved by precedence." << endl;
                exit(1);
            }
        }
    }
}

int PloidyRegions::region_ploidy(const string& ref_path_name, size_t position,
                                 int fallback) const {
    if (by_contig.empty()) {
        return fallback;
    }
    // Match on the contig as the VCF spells it, so a BED written against the output works.
    // Same reduction emit_variant applies when it sets sequenceName.
    string contig = Paths::strip_subrange(ref_path_name);
    string locus = PathMetadata::parse_locus_name(contig);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        contig = locus;
    }
    auto found = by_contig.find(contig);
    if (found == by_contig.end()) {
        return fallback;
    }
    const vector<Region>& regions = found->second;
    // First region starting after the position; its predecessor is the only one that can cover,
    // since the regions are non-overlapping.
    auto it = upper_bound(regions.begin(), regions.end(), position,
                          [](size_t p, const Region& r) { return p < r.start; });
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
    // Same arithmetic emit_variant uses for POS, minus the +1 that makes VCF 1-based: the BED is
    // 0-based, so the two agree on which base an interval boundary falls on.
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
