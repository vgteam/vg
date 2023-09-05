#include <iostream>
#include <fstream>
#include <cassert>
#include "region.hpp"

namespace vg {

void parse_region(const string& target, string& name, int64_t& start, int64_t& end) {
    start = -1;
    end = -1;
    size_t foundLastColon = target.rfind(":");
    // we only have a single string, use the whole sequence as the target
    if (foundLastColon == string::npos) {
        name = target;
    } else {
        name = target.substr(0, foundLastColon);
	    size_t foundRangeDash = target.find("-", foundLastColon);
        if (foundRangeDash == string::npos) {
            start = atoi(target.substr(foundLastColon + 1).c_str());
            end = start;
        } else {
            start = atoi(target.substr(foundLastColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
            end = atoi(target.substr(foundRangeDash + 1).c_str());
        }
    }
}

void parse_bed_regions(const string& bed_path,
                       vector<Region>& out_regions,
                       vector<string>* out_names) {
    out_regions.clear();
    ifstream bedstream(bed_path);
    if (!bedstream) {
        cerr << "Unable to open bed file: " << bed_path << endl;
        return;
    }
    string row;
    string sbuf;
    string ebuf;
    string nbuf;
    for (int line = 1; getline(bedstream, row); ++line) {
        Region region;
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        if (!getline(ss, region.seq, '\t') ||
            !getline(ss, sbuf, '\t') ||
            !getline(ss, ebuf, '\t') ||
            (out_names != nullptr && !getline(ss, nbuf, '\t'))) {
            cerr << "Error parsing bed line " << line << ": " << row << endl;
        } else {
            region.start = std::stoi(sbuf);
            region.end = std::stoi(ebuf);
            assert(region.end > region.start);
            
            // convert from BED-style to 0-based inclusive coordinates
            region.end -= 1;

            out_regions.push_back(region);
            
            if (out_names != nullptr) {
                out_names->push_back(nbuf);
            }
        }
    }
}


}
