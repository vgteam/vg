#ifndef VG_READFILTER_HPP
#define VG_READFILTER_HPP

#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include "vg.hpp"
#include "xg.hpp"
#include "vg.pb.h"

/**
 * Provides a way to filter and transform reads, implementing the bulk of the
 * `vg filter` command.
 *
 */
namespace vg{

using namespace std;

class ReadFilter{
public:
    
    // Filtering parameters
    double min_secondary = 0.;
    double min_primary = 0.;
    double min_sec_delta = 0.;
    double min_pri_delta = 0.;
    bool frac_score = false;
    bool frac_delta = false;
    bool sub_score = false;
    int max_overhang = 99999;
    int context_size = 0;
    bool verbose = false;
    double min_mapq = 0.;
    int repeat_size = 0;
    
    // Extra filename things we need for chunking. TODO: refactor that somehow
    // to maybe be a different class?
    string xg_name;
    string regions_file;
    string outbase;
    
    /**
     * Filter the alignments available from the given stream, placing them on
     * standard output or in the appropriate file. Returns 0 on success, exit
     * code to use on error.
     *
     * TODO: Refactor to be less CLI-aware and more modular-y.
     */
    int filter(istream* alignment_stream);

};
}

#endif
