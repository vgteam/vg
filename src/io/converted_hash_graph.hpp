#ifndef VG_IO_CONVERTED_HASH_GRAPH_HPP_INCLUDED
#define VG_IO_CONVERTED_HASH_GRAPH_HPP_INCLUDED

#include <bdsg/hash_graph.hpp>

namespace vg {

namespace io {

/**
 * Define a type that inherits HashGraph so we can tell the difference between
 * a real HashGraph and a HashGraph converted at load time. We care about this
 * so that vg stats can tell you the original input file format.
 */
class ConvertedHashGraph : public bdsg::HashGraph {
    using bdsg::HashGraph::HashGraph;
};

}

}

#endif
