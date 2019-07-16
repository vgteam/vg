#include "path_string.hpp"

namespace vg {
namespace algorithms {

using namespace std;

size_t min_approx_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos1, const pos_t& pos2, uint64_t max_search) {
    auto nearest1 = nearest_offsets_in_paths(graph, pos1);
    auto nearest2 = nearest_offsets_in_paths(graph, pos2);
    size_t min_distance = numeric_limits<size_t>::max();
    for (auto& p : nearest1) {
        auto q = nearest2.find(p.first);
        if (q != nearest2.end()) {
            // note, doesn't respect orientation
            for (auto& o1 : p.second) {
                for (auto& o2 : q.second) {
                    min_distance = std::min(min_distance, abs(o1-o2));
                }
            }
        }
    }
    return min_distance;
}
    
}
}
