/**
 * \file normalize.cpp
 *
 * Defines an algorithm to normalize a graph.
 */

#include "normalize.hpp"
#include "simplify_siblings.hpp"

#include <unordered_set>
#include <list>
#include <set>
#include <iostream>
#include <sstream>

namespace vg {
namespace algorithms {

using namespace std;

void normalize(handlegraph::MutablePathDeletableHandleGraph* graph, int max_iter, bool debug,
               function<bool(const handle_t&, const handle_t&)> can_merge) {

    size_t last_len = 0;
    if (max_iter > 1) {
        last_len = graph->get_total_length();
    }
    int iter = 0;
    do {
        // Ignore doubly reversing edges; that's not really a coherent concept
        // for all handle graphs, or an obstacle to normality.
        
        // combine diced/chopped nodes (subpaths with no branching)
        handlealgs::unchop(*graph);
        // Resolve forks that shouldn't be
        simplify_siblings(graph, can_merge);
        
        if (max_iter > 1) {
            size_t curr_len = graph->get_total_length();
            if (debug) cerr << "[algorithms::normalize] iteration " << iter+1 << " current length " << curr_len << endl;
            if (curr_len == last_len) break;
            last_len = curr_len;
        }
    } while (++iter < max_iter);
    if (max_iter > 1) {
        if (debug) cerr << "[algorithms::normalize] normalized in " << iter << " steps" << endl;
    }
    
    // there may now be some cut nodes that can be simplified
    // This won't change the length.
    handlealgs::unchop(*graph);
}


}
}

