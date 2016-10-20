//
//  multipath_alignment.hpp
//
// utility functions for the MultipathAlignment protobuf object
//

#ifndef multipath_alignment_hpp
#define multipath_alignment_hpp

#include <stdio.h>
#include <vector>
#include <list>
#include <unordered_map>
#include "vg.pb.h"
#include "path.hpp"
#include "alignment.hpp"
#include "utility.hpp"

#endif /* multipath_alignment_hpp */

namespace vg {
    
    // finds the start subpaths (i.e. the source nodes of the multipath DAG) and stores
    // them in the 'start' field of the MultipathAlignment
    void identify_start_subpaths(MultipathAlignment& multipath_aln);
    
    // stores the highest scoring alignment contained in the MultipathAlignment in an Alignment
    //
    // note: assumes that each subpath's Path object uses one Mapping per node and that
    // subpaths have been identified
    //
    //  Args:
    //    multipath_aln     multipath alignment to find optimal path through
    //    aln_out           empty alignment to store optimal alignment in (data will be
    //                      overwritten if not empty)
    //
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out);
    
    // stores the reverse complement of a MultipathAlignment in another MultipathAlignment
    //
    //  Args:
    //    multipath_aln     multipath alignment to reverse complement
    //    node_length       a function that returns the length of a node sequence from its node ID
    //    rev_comp_out      empty multipath alignment to store reverse complement in (some data may
    //                      be overwritten if not empty)
    //
    void rev_comp_multipath_alignment(const MultipathAlignment& multipath_aln,
                                      const function<int64_t(int64_t)>& node_length,
                                      MultipathAlignment& rev_comp_out);
    
    // converts a Alignment into a Multipath alignment with one Subpath and stores it in an object
    //
    //  Args:
    //    aln               alignment to convert
    //    multipath_aln     empty multipath alignment to store converted alignment in (data may be
    //                      be overwritten if not empty)
    void to_multipath_alignment(const Alignment& aln, MultipathAlignment& multipath_aln_out);
}





