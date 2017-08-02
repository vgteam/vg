#include <iostream>
#include <fstream>
#include <cassert>
#include "region.hpp"

namespace vg {

void parse_bed_regions(const string& bed_path,
                       vector<Region>& out_regions) {
    out_regions.clear();
    ifstream bedstream(bed_path);
    if (!bedstream) {
        cerr << "Unable to open bed file: " << bed_path << endl;
        return;
    }
    string row;
    string sbuf;
    string ebuf;
    for (int line = 1; getline(bedstream, row); ++line) {
        Region region;
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        if (!getline(ss, region.seq, '\t') ||
            !getline(ss, sbuf, '\t') ||
            !getline(ss, ebuf, '\t')) {
            cerr << "Error parsing bed line " << line << ": " << row << endl;
        } else {
            region.start = std::stoi(sbuf);
            region.end = std::stoi(ebuf);
            assert(region.end > region.start);
            
            // convert from BED-style to 0-based inclusive coordinates
            region.end -= 1;

            out_regions.push_back(region);
        }
    }
}


}
