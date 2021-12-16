#ifndef VG_GBZGRAPH_HPP_INCLUDED
#define VG_GBZGRAPH_HPP_INCLUDED

/**
 * \file gbzgraph.hpp
 * Defines a GBZGraph that owns a GBZ object and implements PathHandleGraph.
 */
 
#include "handle.hpp"

#include <bdsg/graph_proxy.hpp>

#include <gbwtgraph/gbz.h>

namespace vg {

/**
 * A PathHandleGraph that owns and is backed by a GBZ.
 * Necessary because GBWTGraph implements PathHandleGraph but GBZ doesn't.
 *
 * Should be removed if/when GBZ implements PathHandleGraph.
 */
class GBZGraph : bdsg::PathHandleGraphProxy<gbwtgraph::GBWTGraph> {
public:
    /// This is the GBZ object we own that actually holds the graph and GBWT
    /// data.
    gbwtgraph::GBZ gbz;
    
protected:
    /**
     * Get the object that actually provides the graph methods.
     */
    inline gbwtgraph::GBWTGraph* get() {
        return &gbz.graph;
    }
    
    /**
     * Get the object that actually provides the graph methods.
     */
    inline const gbwtgraph::GBWTGraph* get() const {
        return &gbz.graph;
    }
    
};

}

#endif

