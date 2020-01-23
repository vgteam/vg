#ifndef VG_ALGORITHMS_STRONGLY_CONNECTED_COMPONENTS_HPP_INCLUDED
#define VG_ALGORITHMS_STRONGLY_CONNECTED_COMPONENTS_HPP_INCLUDED

#include <unordered_set>
#include "../handle.hpp"
#include "dfs.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// TODO: is this the return value we really want?
/// Identify strongly connected components
vector<unordered_set<id_t>> strongly_connected_components(const HandleGraph* g);
    
}
}

#endif
