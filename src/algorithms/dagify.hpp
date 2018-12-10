#ifndef VG_DAGIFY_HPP_INCLUDED
#define VG_DAGIFY_HPP_INCLUDED


#include "../handle.hpp"
#include "../subgraph.hpp"
#include "shortest_cycle.hpp"
#include "eades_algorithm.hpp"
#include "strongly_connected_components.hpp"
#include "is_single_stranded.hpp"
#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

    // must have single stranded orientation
     unordered_map<id_t, pair<id_t, bool>> dagify(const HandleGraph* graph, MutableHandleGraph* into,
                                                  size_t min_preserved_path_length);
}
}

#endif
