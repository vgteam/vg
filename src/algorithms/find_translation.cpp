/**
 * \file find_translation.cpp
 */

#include "find_translation.hpp"

#include "../io/save_handle_graph.hpp"
#include "../gbzgraph.hpp"

namespace vg {
namespace algorithms {

const NamedNodeBackTranslation* find_translation(const HandleGraph* graph) {
    if (!graph) {
        // No graph means no translation.
        return nullptr;
    }
    if (dynamic_cast<const NamedNodeBackTranslation*>(graph)) {
        // Some graph implementations (like GBWTGraph) just are a NamedNodeBackTranslation already.
        return dynamic_cast<const NamedNodeBackTranslation*>(graph);
    }
    if (dynamic_cast<const GFAHandleGraph*>(graph)) {
        // If we loaded the graph from a GFA we would have attached this translation.
        return &dynamic_cast<const GFAHandleGraph*>(graph)->gfa_id_space;
    }
    if (dynamic_cast<const GBZGraph*>(graph)) {
        // If we loaded the graph from a GBZ we can use the GBWTGraph even though the back-translation stuff isn't proxied.
        // TODO: proxy it?
        return &dynamic_cast<const GBZGraph*>(graph)->gbz.graph;
    }
    // Otherwise there's no applicable translation
    return nullptr;
}

}
}
