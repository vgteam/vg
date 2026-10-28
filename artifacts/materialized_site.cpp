/** \file
 * Implements CandidateAlleles and MaterializedSite.
 */

#include "materialized_site.hpp"

namespace vg {

const StepSequences* CandidateAlleles::own_sequences(size_t allele) const {
    return vg::own_sequences(&sequences, allele);
}

bool MaterializedSite::bounds_inside_nodes(const HandleGraph& graph) const {
    return start_offset != graph.get_length(start) || end_offset != 0;
}

}
