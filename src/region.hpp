#ifndef VG_REGION_HPP_INCLUDED
#define VG_REGION_HPP_INCLUDED

#include <string>
#include <vector>
#include <sstream>

namespace vg {

using namespace std;

// Represent a parsed genomic region.
// A -1 for start or end indicates that that coordinate is not used.
// Generally regions parsed form user input will be 1-based.
struct Region {
    string seq;
    int64_t start = -1;
    int64_t end = -1;
};

// Parse a genomic contig[:start-end] region. Outputs -1 for missing start or end.
void parse_region(const string& target, string& name, int64_t& start, int64_t& end);


// Parse a genomic contig[:start-end] region. Outputs -1 for missing start or end.
inline void parse_region(string& region,
                         Region& out_region) {
    parse_region(region,
                 out_region.seq,
                 out_region.start,
                 out_region.end);
}
    
// parse a bed file and return a list of regions (like above)
// IMPORTANT: expects usual 0-based BED format.
// So bedline "chr1   5   10" will return start=5 stop=9
void parse_bed_regions(
    const string& bed_path,
    vector<Region>& out_regions,
    vector<string>* out_names = nullptr);
    
}    

#endif
