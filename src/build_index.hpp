#ifndef VG_BUILD_INDEX_HPP_INCLUDED
#define VG_BUILD_INDEX_HPP_INCLUDED

#include <vg/vg.pb.h>
#include <iostream>
#include "handle.hpp"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "kmer.hpp"

/** \file 
 * Functions for building GCSA2/LCP indexes from HandleGraphs
 */

namespace vg {

using namespace std;

void build_gcsa_lcp(const HandleGraph& graph,
                    gcsa::GCSA*& gcsa,
                    gcsa::LCPArray*& lcp,
                    int kmer_size,
                    size_t doubling_steps = 3,
                    size_t size_limit = 500,
                    const string& base_file_name = "vg-kmers-tmp-");

}

#endif
