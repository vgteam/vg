/// \file multipath_alignment.hpp
///
/// utility functions for the MultipathAlignment protobuf object
///

#ifndef multipath_alignment_hpp
#define multipath_alignment_hpp

#include <stdio.h>
#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>

#include <vg/vg.pb.h>
#include "path.hpp"
#include "alignment.hpp"
#include "utility.hpp"
#include "handle.hpp"

// Declare the haplo::ScoreProvider we use for haplotype-aware traceback generation.
namespace haplo {
    class ScoreProvider;
}

namespace vg {
    
    /// Put subpaths in topological order (assumed to be true for other algorithms)
    void topologically_order_subpaths(MultipathAlignment& multipath_aln);
    
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
    ///    subpath_global    if true, only allows alignments that source subpath to sink subpath
    ///                      in the multipath DAG, else allows any start and end subpath
    ///
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out,
                           bool subpath_global = false);
    
    /// Returns the score of the highest scoring alignment contained in the MultipathAlignment
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal score in
    ///    subpath_global    if true, only allows alignments that source subpath to sink subpath
    ///                      in the multipath DAG, else allows any start and end subpath
    ///
    int32_t optimal_alignment_score(const MultipathAlignment& multipath_aln,
                                    bool subpath_global = false);
    
    /// Returns the top k highest-scoring alignments contained in the MultipathAlignment.
    /// Note that some or all of these may be duplicate Alignments, which were spelled out
    /// by tracebacks through different sequences of subpaths that shared alignment material.
    ///
    /// If the best alignment is no alignment (i.e. the read is unmapped), returns an empty vector.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    count             maximum number of top alignments to return
    ///
    vector<Alignment> optimal_alignments(const MultipathAlignment& multipath_aln, size_t count);
    
    /// Finds k or fewer top-scoring alignments using only distinct subpaths.
    /// Asymmetrical: the optimal alignment for each end subpath is found, greedily, subject to the constraint,
    /// but the other subpaths are first-come first-serve. Also, distinct subpaths may not guarantee distinct
    /// actual alignments, so alignments may need deduplication.
    ///
    /// If the best alignment is no alignment (i.e. the read is unmapped), returns an empty vector.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    count             maximum number of top alignments to return
    ///
    vector<Alignment> optimal_alignments_with_disjoint_subpaths(const MultipathAlignment& multipath_aln, size_t count);
    
    /// Finds all alignments consistent with haplotypes available by incremental search with the given haplotype
    /// score provider. Pads to a certain count with haplotype-inconsistent alignments that are population-scorable
    /// (i.e. use only edges used by some haplotype in the index), and then with unscorable alignments if scorable
    /// ones are unavailable. This may result in an empty vector.
    ///
    /// Output Alignments may not be unique. The input MultipathAlignment may have exponentially many ways to
    /// spell the same Alignment, and we will look at all of them. We also may have duplicates of the optimal
    /// alignment if we are asked to produce it unconsitionally.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    score_provider    a haplo::ScoreProvider that supports incremental search over its haplotype database (such as a GBWTScoreProvider)
    ///    soft_count        maximum number of haplotype-inconsistent alignments to pad to
    ///    hard_count        maximum number of alignments, including haplotype-consistent (0 if no limit)
    ///    optimal_first     always compute and return first the optimal alignment, even if not haplotype-consistent
    ///
    vector<Alignment> haplotype_consistent_alignments(const MultipathAlignment& multipath_aln, const haplo::ScoreProvider& score_provider,
        size_t soft_count, size_t hard_count, bool optimal_first = false);
    
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
    
    /// Copies metadata from an MultipathAlignment object and transfers it to another MultipathAlignment
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const MultipathAlignment& from, MultipathAlignment& to);
    
    /// Merges non-branching paths in a multipath alignment in place
    void merge_non_branching_subpaths(MultipathAlignment& multipath_aln);
    
    /// Returns a vector whose elements are vectors with the indexes of the Subpaths in
    /// each connected component. An unmapped MultipathAlignment with no subpaths produces an empty vector.
    vector<vector<int64_t>> connected_components(const MultipathAlignment& multipath_aln);
    
    /// Extract the MultipathAlignment consisting of the Subpaths with the given indexes
    /// into a new MultipathAlignment object
    void extract_sub_multipath_alignment(const MultipathAlignment& multipath_aln,
                                         const vector<int64_t>& subpath_indexes,
                                         MultipathAlignment& sub_multipath_aln);
    
    /// Debugging function to check that multipath alignment meets the formalism's basic
    /// invariants. Returns true if multipath alignment is valid, else false. Does not
    /// validate alignment score.
    bool validate_multipath_alignment(const MultipathAlignment& multipath_aln, const HandleGraph& handle_graph);
    
    /// Send a formatted string representation of the MultipathAlignment into the ostream
    void view_multipath_alignment(ostream& out, const MultipathAlignment& multipath_aln, const HandleGraph& handle_graph);
    
    /// Converts a MultipathAlignment to a GraphViz Dot representation, output to the given ostream.
    void view_multipath_alignment_as_dot(ostream& out, const MultipathAlignment& multipath_aln, bool show_graph = false);
    
    // TODO: function for adding a graph augmentation to an existing multipath alignment
}


#endif /* multipath_alignment_hpp */




