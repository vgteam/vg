/**
 * \file reverse_complement_graph.cpp: contains the implementation of ReverseComplementGraph
 */


#include "reverse_complement_graph.hpp"
#include "utility.hpp"


namespace vg {

using namespace std;

    ReverseComplementGraph::ReverseComplementGraph(const HandleGraph* forward_graph) : BackwardsGraph(forward_graph) {
        // nothing to do
    }
    
    string ReverseComplementGraph::get_sequence(const handle_t& handle) const {
        // reverse and complement the sequence
        string sequence = forward_graph->get_sequence(handle);
        return reverse_complement(sequence);
    }
}

