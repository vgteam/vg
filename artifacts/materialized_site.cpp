/** \file
 * Implements the spelling of alleles, CandidateAlleles and MaterializedSite.
 */

#include "materialized_site.hpp"

namespace vg {

const StepSequences* own_sequences(const AlleleSequences* sequences, size_t allele) {
    if (sequences == nullptr || allele >= sequences->size() || (*sequences)[allele].empty()) {
        return nullptr;
    }
    return &(*sequences)[allele];
}

string allele_sequence(const HandleGraph& graph, const Traversal& walk, const StepSequences* own) {
    string sequence;
    if (own != nullptr) {
        for (const string& step : *own) {
            sequence += step;
        }
    } else {
        for (const handle_t& handle : walk) {
            sequence += graph.get_sequence(handle);
        }
    }
    return sequence;
}

const StepSequences* CandidateAlleles::own_sequences(size_t allele) const {
    return vg::own_sequences(&sequences, allele);
}

bool MaterializedSite::bounds_inside_nodes(const HandleGraph& graph) const {
    return start_offset != graph.get_length(start) || end_offset != 0;
}

}
