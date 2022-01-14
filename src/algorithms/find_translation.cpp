/**
 * \file find_translation.cpp
 */

#include "find_translation.hpp"

#include "../io/save_handle_graph.hpp"

namespace vg {
namespace algorithms {

const NamedNodeBackTranslation* find_translation(const HandleGraph* graph) {
    if (!graph) {
        // No graph means no translation.
        return nullptr;
    }
    if (dynamic_cast<const GFAHandleGraph*>(graph)) {
        // If we loaded the graph fropm a GFA we would have attached this translation.
        return &dynamic_cast<const GFAHandleGraph*>(graph)->gfa_id_space;
    }
    // Otherwise there's no applicable translation
    return nullptr;
}

}
}
