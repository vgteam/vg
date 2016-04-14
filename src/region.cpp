#include <iostream>
#include <fstream>
#include <cassert>
#include "region.hpp"

namespace vg {

void parse_region(
    string& region,
    string& startSeq,
    int& startPos,
    int& stopPos) {

    size_t foundFirstColon = region.find(":");

    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        startSeq = region;
        startPos = 0;
        stopPos = -1;
    } else {
        startSeq = region.substr(0, foundFirstColon);
        string sep = "..";
        size_t foundRangeSep = region.find(sep, foundFirstColon);
        if (foundRangeSep == string::npos) {
            sep = "-";
            foundRangeSep = region.find("-", foundFirstColon);
        }
        if (foundRangeSep == string::npos) {
            startPos = atoi(region.substr(foundFirstColon + 1).c_str());
            // differ from bamtools in this regard, in that we process only
            // the specified position if a range isn't given
            stopPos = startPos + 1;
        } else {
            startPos = atoi(region.substr(foundFirstColon + 1, foundRangeSep - foundFirstColon).c_str());
            // if we have range sep specified, but no second number, read to the end of sequence
            if (foundRangeSep + sep.size() != region.size()) {
                stopPos = atoi(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
            } else {
                //stopPos = reference.sequenceLength(startSeq);
                stopPos = -1;
            }
        }
    }
}

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
            
            // convert from BED-style to VCF-style coordinates
            region.start += 1;

            out_regions.push_back(region);
        }
    }
}


}
