#ifndef VG_UTILITY_H
#define VG_UTILITY_H

#include <string>
#include <vector>
#include <sstream>

namespace vg {

using namespace std;

// wrap up region parse output
struct Region {
    string seq;
    int start;
    int end;
};

// parse a region in the form: chrom:start-stop
// (if only chrom given, start, stop will be set to 0,-1)
void parse_region(
    string& region,
    string& startSeq,
    int& startPos,
    int& stopPos);


inline void parse_region(string& region,
                         Region& out_region) {
    parse_region(region,
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
