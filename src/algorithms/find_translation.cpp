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
    if (const gbwtgraph::GBWTGraph* gg = dynamic_cast<const gbwtgraph::GBWTGraph*>(graph)) {
        // A GBWTGraph is a NamedNodeBackTranslation, but the translation may not be present.
        if (gg->has_segment_names()) { return gg; }
        else { return nullptr; }
    }
    if (const GBZGraph* gg = dynamic_cast<const GBZGraph*>(graph)) {
        // The same goes for a GBWTGraph contained in a GBZ graph.
        if (gg->gbz.graph.has_segment_names()) { return &(gg->gbz.graph); }
        else { return nullptr; }
    }
    if (dynamic_cast<const NamedNodeBackTranslation*>(graph)) {
        // Some graph implementations just are a NamedNodeBackTranslation already.
        return dynamic_cast<const NamedNodeBackTranslation*>(graph);
    }
    if (dynamic_cast<const GFAHandleGraph*>(graph)) {
        // If we loaded the graph from a GFA we would have attached this translation.
        return &dynamic_cast<const GFAHandleGraph*>(graph)->gfa_id_space;
    }
    // Otherwise there's no applicable translation
    return nullptr;
}

}
}
