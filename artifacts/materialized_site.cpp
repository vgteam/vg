/** \file
 * Implements the spelling of alleles, CandidateAlleles and MaterializedSite.
 */

#include "materialized_site.hpp"

#include <algorithm>
#include <tuple>

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

bool CandidateAlleles::spelled_by_walk(size_t allele) const {
    return allele >= edits.size() || edits[allele].empty();
}

StepSequences CandidateAlleles::step_sequences(const HandleGraph& graph, size_t allele) const {
    StepSequences steps;
    if (spelled_by_walk(allele)) {
        return steps;
    }
    const Traversal& walk = walks[allele];
    // Order the changes along the walk, insertions first among those at one offset.
    vector<const AlleleEdit*> changes;
    for (const AlleleEdit& change : edits[allele]) {
        changes.push_back(&change);
    }
    std::ranges::sort(changes, [](const AlleleEdit* a, const AlleleEdit* b) {
        return std::tie(a->step, a->at, a->removed) < std::tie(b->step, b->at, b->removed);
    });
    // Copy each step's own sequence between the changes, and each change's bases in place of the
    // bases it removes.
    steps.reserve(walk.size());
    auto change = changes.begin();
    for (size_t k = 0; k < walk.size(); ++k) {
        const string own = graph.get_sequence(walk[k]);
        string& spelled = steps.emplace_back();
        size_t copied = 0;
        for (; change != changes.end() && (*change)->step == k; ++change) {
            spelled.append(own, copied, (*change)->at - copied);
            spelled += (*change)->added;
            copied = (*change)->at + (*change)->removed;
        }
        spelled.append(own, copied);
    }
    return steps;
}

AlleleSequences CandidateAlleles::all_step_sequences(const HandleGraph& graph) const {
    AlleleSequences all;
    if (edits.empty()) {
        return all;
    }
    all.reserve(walks.size());
    for (size_t a = 0; a < walks.size(); ++a) {
        all.push_back(step_sequences(graph, a));
    }
    return all;
}

string CandidateAlleles::sequence(const HandleGraph& graph, size_t allele) const {
    if (spelled_by_walk(allele)) {
        return allele_sequence(graph, walks[allele], nullptr);
    }
    const StepSequences steps = step_sequences(graph, allele);
    return allele_sequence(graph, walks[allele], &steps);
}

bool MaterializedSite::bounds_inside_nodes(const HandleGraph& graph) const {
    return start_offset != graph.get_length(start) || end_offset != 0;
}

}
