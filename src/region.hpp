#ifndef VG_REGION_HPP_INCLUDED
#define VG_REGION_HPP_INCLUDED

#include <string>
#include <vector>
#include <sstream>
#include "xg.hpp"

namespace vg {

using namespace std;

// wrap up region parse output
struct Region {
    string seq;
    int64_t start;
    int64_t end;
};


inline void parse_region(string& region,
                         Region& out_region) {
    xg::parse_region(region,
                     out_region.seq,
                     out_region.start,
                     out_region.end);
}
    
// parse a bed file and return a list of regions (like above)
// IMPORTANT: expects usual 0-based BED format.
// So bedline "chr1   5   10" will return start=6 stop=10
void parse_bed_regions(
    const string& bed_path,
    vector<Region>& out_regions);
    
}    

#endif
