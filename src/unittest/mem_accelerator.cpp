/// \file path_index.cpp
///  
/// Unit tests for the PathIndex class, which indexes paths for random access.
///

#include <iostream>
#include <string>
#include <sstream>
#include "../mem_accelerator.hpp"
#include "../utility.hpp"
#include "../build_index.hpp"
#include "catch.hpp"
#include "random_graph.hpp"

#include <bdsg/hash_graph.hpp>
#include <gcsa/utils.h>

namespace std {

// For the output operator for ranges
using gcsa::operator<<;

/// Allow ranges to be dumped by the test asserts if they don't match properly.
static std::string to_string(const gcsa::range_type& range) {
    std::stringstream ss;
    ss << range;
    return ss.str();
}

}

namespace vg {
namespace unittest {
using namespace std;

// For the output operator for ranges
using gcsa::operator<<;


TEST_CASE("MEMAccelerator returns same ranges as direct LF queries",
          "[mem][mapping][memaccelerator]" ) {
    
    int num_graphs = 5;
    int seq_size = 100;
    int var_count = 4;
    int var_length = 3;
    int memo_length = 5;
    for (int g = 0; g < num_graphs; ++g) {
        
        bdsg::HashGraph graph;
        random_graph(seq_size, var_length, var_count, &graph);
        
        // Make GCSA quiet
        gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        
        // Make pointers to fill in
        gcsa::GCSA* gcsaidx = nullptr;
        gcsa::LCPArray* lcpidx = nullptr;
        
        // Build the GCSA index
        build_gcsa_lcp(graph, gcsaidx, lcpidx, 8, 2);
        
        MEMAccelerator accelerator(*gcsaidx, memo_length);
        
        // iterate over all k-mers
        for (int k = 0; k < (1 << (2 * memo_length)); ++k) {
            
            string seq(memo_length, 'N');
            for (int i = 0; i < memo_length; ++i) {
                seq[i] = "ACGT"[(k >> i) & 3];
            }
            
            // check that the memoized and direct ranges are equivalent
            auto memo_range = accelerator.memoized_LF(seq.end() - 1);
            auto direct_range = gcsa::range_type(0, gcsaidx->size() - 1);
            auto cursor = seq.end() - 1;
            while (cursor >= seq.begin() && !gcsa::Range::empty(direct_range)) {
                direct_range = gcsaidx->LF(direct_range,
                                           gcsaidx->alpha.char2comp[*cursor]);
                --cursor;
            }
            
            std::cerr << "Kmer " << seq << " memo " << memo_range << " direct " << direct_range << std::endl;
            REQUIRE(memo_range == direct_range);
        }
        
        delete gcsaidx;
        delete lcpidx;
    }
}
   
}
}
        
