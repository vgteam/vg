#pragma once

#include "../handle.hpp"
#include <vg/vg.pb.h>
#include <string>
#include <functional>
#include "nearest_offsets_in_paths.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// use the embedded paths to get an estimated minimum distance between the positions
size_t min_approx_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos1, const pos_t& pos2, uint64_t max_search);

}

}
