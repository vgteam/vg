/**
 * \file find_gbwtgraph.cpp
 */

#include "find_gbwtgraph.hpp"

#include "../gbzgraph.hpp"

namespace vg {
namespace algorithms {

const gbwtgraph::GBWTGraph* find_gbwtgraph(const HandleGraph* graph) {
    if (!graph) {
        // No graph means no translation.
        return nullptr;
    }
    if (dynamic_cast<const gbwt::GBWTGraph*>(graph)) {
        // If it already is one, return it
        return &dynamic_cast<const gbwt::GBWTGraph*>(graph);
    }
    if (dynamic_cast<const GBZGraph*>(graph)) {
        // If it's a GBZGraph, go get the GBWTGraph and return that.
        return &dynamic_cast<const GBZGraph*>(graph)->gbz.graph;
    }
    // Otherwise there's no applicable GBWTGraph
    return nullptr;
}

}
}
