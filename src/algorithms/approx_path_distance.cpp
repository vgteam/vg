#include "approx_path_distance.hpp"

namespace vg {
namespace algorithms {

using namespace std;

size_t min_approx_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos1, const pos_t& pos2, uint64_t max_search) {
    auto nearest1 = nearest_offsets_in_paths(graph, pos1, max_search);
    auto nearest2 = nearest_offsets_in_paths(graph, pos2, max_search);
    uint64_t min_distance = numeric_limits<uint64_t>::max();
    for (auto& p : nearest1) {
        auto q = nearest2.find(p.first);
        if (q != nearest2.end()) {
            // note, doesn't respect orientation
            for (auto& o1 : p.second) {
                for (auto& o2 : q->second) {
                    uint64_t x = (o1.first > o2.first ? o1.first - o2.first : o2.first - o1.first);
                    min_distance = std::min(min_distance, x);
                }
            }
        }
    }
    return (size_t)min_distance;
}
    
}
}
