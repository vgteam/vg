#include "a_star.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    
//    class AStarTree {
//    private:
//        class AStarNode;
//
//    public:
//        AStarTree(const handle_t& start_node, const handle_t& end_node) {
//            head = new AStarNode(start_node, 0);
//        }
//
//        ~AStarTree {
//            vector<AStarNode*> stack(1, head);
//            while (!stack.empty()) {
//                AStarNode* node = stack.pop_back();
//                for
//            }
//        }
//
//    private:
//        AStarNode* head;
//
//        class AStarNode {
//        public:
//            AStarNode(const handle_t& handle, int64_t distance,
//                      AStarNode* parent = nullptr) :
//            handle(handle), distance(distance), parent(parent) { }
//
//            ~AStarNode() = default;
//
//            const handle_t handle;
//            const int64_t distance;
//
//            AStarNode* parent;
//            vector<AStarNode*> children;
//        };
//    };
    
    template<class DistHeuristic>
    vector<handle_t> a_star(const HandleGraph* graph,
                            const pos_t& pos_1, const pos_t& pos_2,
                            const DistHeuristic& dist_heuristic,
                            bool find_min,
                            int64_t extremal_distance) {
        
        struct AStarPathEnd {
            AStarPathEnd(handle_t handle, int64_t distance, int64_t heuristic_estimate, handle_t predecessor)
                : handle(handle), distance(distance), heuristic_estimate(heuristic_estimate), predecessor(predecessor) {
            }
            
            const handle_t handle;
            const int64_t distance;
            const int64_t heuristic_estimate;
            const handle_t predecessor;
            
            inline bool operator<(const AStarPathEnd& other) {
                return ((find_min && heuristic_estimate > other.heuristic_estimate)
                        || (!find_min && other.heuristic_estimate < heuristic_estimate));
            }
        };
        
        // TODO: this function won't handle cycles as written
        // it needs to be able to move past the end node and circle back around to the start
        // node without breaking the predecessor (i.e. have multiple predecessors)
        
        unordered_map<handle_t, handle_t> best_predecessor;
        
        priority_queue<AStarPathEnd> queue;
        
        handle_t start = graph->get_handle(id(pos_1), is_rev(pos_1));
        handle_t end = graph->get_handle(id(pos_2), is_rev(pos_2));
        
        // set negative distance so the search starts from the beginning of the node
        queue.emplace(start, -offset(pos_1), dist_heuristic(make_pos_t(id(pos_1), is_rev(pos_1), 0), pos_2),
                      as_handle(int64_t(0)));
        
        while (!queue.empty()) {
            auto path_end = queue.top();
            queue.pop();
            
            if ((find_min && path_end.distance + offset(pos_2) > extremal_distance)
                || (!find_min && path_end.distance + offset(pos_2) < extremal_distance)) {
                // we've crossed over the distance beyond which we're not interested
                break;
            }
            
            if (path_end.handle == end) {
                // we've found the end, reconstruct the path and return it
                vector<handle_t> path(1, path_end.handle);
                while (best_predecessor.count(path.back())) {
                    path.emplace_back(best_predecessor[path.back()]);
                }
                // put the path in forward order
                reverse(path.begin(), path.end());
                return path;
            }
            
            if (best_predecessor.count(path_end.handle) && path_end.handle != start) {
                continue;
            }
            
            best_predecessor[path_end.handle] = path_end.predecessor;
            
            int64_t distance = graph->get_length(path_end.handle) + path_end.distance;
            
            auto enqueue_next = [&](const handle_t& next) {
                int64_t remaining = dist_heuristic(make_pos_t(graph->get_id(next), graph->get_is_reverse(next), 0),
                                                   pos_2);
                queue.emplace(next, distance, distance + remaining, path_end.handle);
            };
            
            graph->follow_edges(path_end.handle, enqueue_next);
        }
        
        // we made it through the loop without finding a
        return vector<handle_t>();
//        if (find_min) {
//            return a_star_internal<DistHeuristic, less<pair<int64_t, handle_t>>>(graph,
//                                                                                 pos_1, pos_2,
//                                                                                 dist_heuristic,
//                                                                                 extremal_distance);
//        }
//        else {
//            return a_star_internal<DistHeuristic, greater<pair<int64_t, handle_t>>>(graph,
//                                                                                    pos_1, pos_2,
//                                                                                    dist_heuristic,
//                                                                                    extremal_distance);
//        }
    }
    
//    template<class DistHeuristic, class Compare>
//    vector<handle_t> a_star_internal(const HandleGraph* graph,
//                                     const pos_t& pos_1, const pos_t& pos_2
//                                     const DistHeuristic& dist_heuristic,
//                                     int64_t extremal_distance) {
//
//        unordered_set<handle_t> evaluated_nodes;
//
//        unordered_map<handle_t, handle_t> best_predecessor;
//
//        priority_queue<pair<int64_t,  handle_t>,
//                       vector<pair<int64_t, handle_t>>,
//                       Compare> queue;
//
//
//        handle_t start = graph->get_handle(id(pos_1), is_rev(pos_1));
//        handle_t end = graph->get_handle(id(pos_2), is_rev(pos_2));
//
//        queue.emplace(dist_heuristic(start, end), start);
//
//        while (!queue.empty()) {
//            handle_t node;
//            auto next = queue.top();
//            queue.pop();
//
//            evaluated_nodes.emplace(next.second)
//
//            if (next.second == end) {
//                break;
//            }
//        }
//    }

}
}
