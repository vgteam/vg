/**
 * \file deletion_aligner.hpp
 *
 * Defines an aligner for global deletions
 *
 */
#ifndef VG_DINUCLEOTIDE_MACHINE_GRAPH_HPP_INCLUDED
#define VG_DINUCLEOTIDE_MACHINE_GRAPH_HPP_INCLUDED

#include <cstdint>
#include <vg/vg.pb.h>
#include <structures/min_max_heap.hpp>

#include "handle.hpp"
#include "algorithms/topological_sort.hpp"

namespace vg {

using namespace std;

/*
 * An aligner for global deletions. Can only produce alignments
 * for empty sequences
 */
class DeletionAligner {
public:
    DeletionAligner(int8_t gap_open, int8_t gap_extension);
    ~DeletionAligner() = default;
    
    void align(Alignment& aln, const HandleGraph& graph) const;
    
    void align_multi(Alignment& aln, vector<Alignment>& alt_alns,
                     const HandleGraph& graph, int32_t max_alt_alns) const;
    
    
private:
    
    vector<vector<handle_t>> run_dp(const HandleGraph& graph,
                                    int32_t max_tracebacks) const;
    
    pair<vector<size_t>, vector<pair<size_t, size_t>>> min_dists(const vector<handle_t>& order,
                                                                 const unordered_map<handle_t, size_t>& index_of,
                                                                 const HandleGraph& graph) const;
    
    vector<vector<handle_t>> traceback(const vector<handle_t>& order,
                                       const unordered_map<handle_t, size_t>& index_of,
                                       const HandleGraph& graph,
                                       const vector<size_t>& dists,
                                       const vector<pair<size_t, size_t>>& sinks,
                                       size_t max_tracebacks) const;
    
    void trace_to_alignment(Alignment& aln, const vector<handle_t>& trace,
                            const HandleGraph& graph) const;
    
    int32_t gap_open;
    int32_t gap_extension;
    
};

}

#endif
