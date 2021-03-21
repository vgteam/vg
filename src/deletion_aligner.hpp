/**
 * \file deletion_aligner.hpp
 *
 * Defines an aligner for global deletions
 *
 */
#ifndef VG_DELETION_ALIGNER_GRAPH_HPP_INCLUDED
#define VG_DELETION_ALIGNER_GRAPH_HPP_INCLUDED

#include <cstdint>
#include <vg/vg.pb.h>
#include <structures/min_max_heap.hpp>

#include "handle.hpp"

namespace vg {

using namespace std;

/*
 * An aligner for global deletions. Can only produce alignments
 * for empty sequences
 */
class DeletionAligner {
public:
    DeletionAligner(int8_t gap_open, int8_t gap_extension);
    DeletionAligner() = default;
    ~DeletionAligner() = default;
    
    // store a global alignment of an empty sequence to this graph in the alignment
    // crashes if alignment sequence is not empty
    void align(Alignment& aln, const HandleGraph& graph) const;
    
    // store a global alignment of an empty sequence to this graph in the alignment
    // also store the next highest-scoring alignments in the vector
    // crashes if alignment sequence is not empty
    void align_multi(Alignment& aln, vector<Alignment>& alt_alns,
                     const HandleGraph& graph, int32_t max_alt_alns) const;
    
    
private:
    
    // do the dynamic programming and return the trace backs
    vector<vector<handle_t>> run_dp(const HandleGraph& graph,
                                    int32_t max_tracebacks) const;
    
    // calculate min distance with DP. returns the DP table and
    // the sink nodes in the graph, with their full DP value
    pair<vector<size_t>, vector<pair<size_t, size_t>>> min_dists(const vector<handle_t>& order,
                                                                 const unordered_map<handle_t, size_t>& index_of,
                                                                 const HandleGraph& graph) const;
    
    // compute tracebacks using the DP table
    vector<vector<handle_t>> traceback(const vector<handle_t>& order,
                                       const unordered_map<handle_t, size_t>& index_of,
                                       const HandleGraph& graph,
                                       const vector<size_t>& dists,
                                       const vector<pair<size_t, size_t>>& sinks,
                                       size_t max_tracebacks) const;
    
    // convert a trace into an alignment path
    void trace_to_alignment(Alignment& aln, const vector<handle_t>& trace,
                            const HandleGraph& graph) const;
    
    int32_t gap_open;
    int32_t gap_extension;
    
};

}

#endif
