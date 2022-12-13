/**
 * \file find_gbwt.cpp
 */

#include "find_gbwt.hpp"
#include "find_gbwtgraph.hpp"
#include <vg/io/vpkg.hpp>

#include "../gbzgraph.hpp"

namespace vg {
namespace algorithms {

const gbwt::GBWT* find_gbwt(const HandleGraph* graph) {
    const gbwtgraph::GBWTGraph* typed_graph = find_gbwtgraph(graph);
    if (!typed_graph) {
        return nullptr;
    }
    return typed_graph->index;
}

const gbwt::GBWT* find_gbwt(const HandleGraph* graph, std::unique_ptr<gbwt::GBWT>& holder, const std::string& filename) {
    if (!filename.empty()) {
        holder = vg::io::VPKG::load_one<gbwt::GBWT>(filename);

        if (holder.get() == nullptr) {
          // Complain if we couldn't get it but were supposed to.
          cerr << "error:[vg::algorithms::find_gbwt] unable to load gbwt index file " << filename << endl;
          exit(1);
        }
        
        return holder.get();
    }
    // If we don't need to load it, try and get it from the graph.
    return find_gbwt(graph);
}

}
}
