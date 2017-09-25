//
//  vg_algorithms.hpp
//  

// A collection of utility algorithms to manipulating graphs that can be use by both VG and XG objects
//

#ifndef VG_VG_ALGORITHMS_HPP
#define VG_VG_ALGORITHMS_HPP

// Actual algorithms are split up into different files.
#include "extract_connecting_graph.hpp"
#include "extract_containing_graph.hpp"
#include "extract_extending_graph.hpp"

#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <vector>

#include "../position.hpp"
#include "../cached_position.hpp"
#include "../xg.hpp"
#include "../vg.hpp"
#include "../vg.pb.h"
#include "../hash_map.hpp"
#include "../json2pb.h"

namespace vg {
namespace algorithms {
}
}

#endif /* vg_algorithms_hpp */
