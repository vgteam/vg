#ifndef VG_GRAPH_HPP_INCLUDED
#define VG_GRAPH_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "handle.hpp"
#include <set>
#include <algorithm>

namespace vg {

using namespace std;

// transfer data from a HandleGraph into an empty Graph
void from_handle_graph(const HandleGraph& from, Graph& to);

// transfer data from a PathHandleGraph into an empty Graph
void from_path_handle_graph(const PathHandleGraph& from, Graph& to);
    
}

#endif
