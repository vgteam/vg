#ifndef VG_DISJOINT_COMPONENTS_HPP_INCLUDED
#define VG_DISJOINT_COMPONENTS_HPP_INCLUDED

#include <list>

#include <bdsg/hash_graph.hpp>

#include "../handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Return a list of graphs, one for each connected component in the original graph.
/// Node IDs are preserved from the original graph.
list<bdsg::HashGraph> disjoint_components(const HandleGraph& graph);

}
}

#endif
