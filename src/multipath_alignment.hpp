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
#include <algorithm>
#include "vg.pb.h"
#include "path.hpp"
#include "alignment.hpp"
#include "utility.hpp"

namespace vg {
    
    /// Finds the start subpaths (i.e. the source nodes of the multipath DAG) and stores
    /// them in the 'start' field of the MultipathAlignment
    void identify_start_subpaths(MultipathAlignment& multipath_aln);
    
    /// Stores the highest scoring alignment contained in the MultipathAlignment in an Alignment
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal path through
    ///    aln_out           empty alignment to store optimal alignment in (data will be
    ///                      overwritten if not empty)
    ///
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out);
    
    /// Returns the score of the highest scoring alignment contained in the MultipathAlignment
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal score in
    ///
    int32_t optimal_alignment_score(const MultipathAlignment& multipath_aln);
    
    /// Stores the reverse complement of a MultipathAlignment in another MultipathAlignment
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to reverse complement
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///    rev_comp_out      empty multipath alignment to store reverse complement in (some data may
    ///                      be overwritten if not empty)
    ///
    void rev_comp_multipath_alignment(const MultipathAlignment& multipath_aln,
                                      const function<int64_t(int64_t)>& node_length,
                                      MultipathAlignment& rev_comp_out);
    
    /// Stores the reverse complement of a MultipathAlignment in another MultipathAlignment
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to reverse complement in place
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///
    void rev_comp_multipath_alignment_in_place(MultipathAlignment* multipath_aln,
                                               const function<int64_t(int64_t)>& node_length);
    
    /// Converts a Alignment into a Multipath alignment with one Subpath and stores it in an object
    ///
    ///  Args:
    ///    aln               alignment to convert
    ///    multipath_aln     empty multipath alignment to store converted alignment in (data may be
    ///                      be overwritten if not empty)
    ///
    void to_multipath_alignment(const Alignment& aln, MultipathAlignment& multipath_aln_out);
    
    /// Copies metadata from an Alignment object and transfers it to a MultipathAlignment
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const Alignment& from, MultipathAlignment& to);
    
    /// Copies metadata from an MultipathAlignment object and transfers it to a Alignment
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const MultipathAlignment& from, Alignment& to);
    
    // TODO: function for adding a graph augmentation to an existing multipath alignment
}


#endif /* multipath_alignment_hpp */




