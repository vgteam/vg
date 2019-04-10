#include "are_equivalent.hpp"

namespace vg {
namespace algorithms {

    bool are_equivalent(const HandleGraph* graph_1,
                        const HandleGraph* graph_2, bool verbose) {
        
        
        if (graph_1->node_size() != graph_2->node_size()) {
            if (verbose) {
                cerr << "graphs have different numbers of nodes: " << graph_1->node_size() << " and " << graph_2->node_size() << endl;
            }
            return false;
        }
        
        bool equivalent = true;
        graph_1->for_each_handle([&](const handle_t& handle_1) {
            
            if (!graph_2->has_node(graph_1->get_id(handle_1))) {
                if (verbose) {
                    cerr << "node " << graph_1->get_id(handle_1) << " from graph 1 was not found in graph 2" << endl;
                }
                equivalent = false;
                return false;
            }
            
            handle_t handle_2 = graph_2->get_handle(graph_1->get_id(handle_1),
                                                    graph_1->get_is_reverse(handle_1));
            
            if (graph_1->get_sequence(handle_1) != graph_2->get_sequence(handle_2)) {
                if (verbose) {
                    cerr << "node " << graph_1->get_id(handle_1) << " has different sequence in the two graphs: " << graph_1->get_sequence(handle_1) << " and " << graph_2->get_sequence(handle_2) << endl;
                }
                equivalent = false;
                return false;
            }
            
            for (bool direction : {true, false}) {
                vector<handle_t> nexts_1, nexts_2;
                graph_1->follow_edges(handle_1, direction, [&](const handle_t& next) {
                    nexts_1.push_back(next);
                });
                graph_2->follow_edges(handle_2, direction, [&](const handle_t& next) {
                    nexts_2.push_back(next);
                });
                
                if (nexts_1.size() != nexts_2.size()) {
                    if (verbose) {
                        cerr << "node " << graph_1->get_id(handle_1) << " has a different number of edges on the " << (direction ? "left" : "right") << " side in the two graphs: " << nexts_1.size() << " and " << nexts_2.size() << endl;
                    }
                    equivalent = false;
                    return false;
                }
                
                sort(nexts_1.begin(), nexts_1.end(), [&](const handle_t& a, const handle_t& b) {
                    return (graph_1->get_id(a) < graph_1->get_id(b) ||
                            (graph_1->get_id(a) == graph_1->get_id(b) &&
                             graph_1->get_is_reverse(a) < graph_1->get_is_reverse(b)));
                });
                sort(nexts_2.begin(), nexts_2.end(), [&](const handle_t& a, const handle_t& b) {
                    return (graph_2->get_id(a) < graph_2->get_id(b) ||
                            (graph_2->get_id(a) == graph_2->get_id(b) &&
                             graph_2->get_is_reverse(a) < graph_2->get_is_reverse(b)));
                });
                
                for (size_t i = 0; i < nexts_1.size(); i++) {
                    if (graph_1->get_id(nexts_1[i]) != graph_2->get_id(nexts_2[i]) ||
                        graph_1->get_is_reverse(nexts_1[i]) != graph_2->get_is_reverse(nexts_2[i])) {
                        if (verbose) {
                            cerr << "node " << graph_1->get_id(handle_1) << " has edges to different nodes on the " << (direction ? "left" : "right") << " side in the two graphs" << endl;
                        }
                        equivalent = false;
                        return false;
                    }
                }
            }
            return true;
        });
        
        return equivalent;
    }
    
    bool are_equivalent_with_paths(const PathHandleGraph* graph_1,
                                   const PathHandleGraph* graph_2, bool verbose) {
        
        if (!are_equivalent(graph_1, graph_2, verbose)) {
            return false;
        }
        
        if (graph_1->get_path_count() != graph_2->get_path_count()) {
            if (verbose) {
                cerr << "graphs have different numbers of paths: " << graph_1->get_path_count() << " and " << graph_2->get_path_count() << endl;
            }
            return false;
        }
        
        bool equivalent = true;
        graph_1->for_each_path_handle([&](const path_handle_t& path_handle_1){
            if (!graph_2->has_path(graph_1->get_path_name(path_handle_1))) {
                if (verbose) {
                    cerr << "path " << graph_1->get_path_name(path_handle_1) << " from graph 1 was not found in graph 2" << endl;
                }
                equivalent = false;
                return false;
            }
            
            path_handle_t path_handle_2 = graph_2->get_path_handle(graph_1->get_path_name(path_handle_1));
            
            if (graph_1->get_occurrence_count(path_handle_1) != graph_2->get_occurrence_count(path_handle_2)) {
                equivalent = false;
                return false;
            }
            
            // TODO: if paths are circular, we shouldn't enforce that the start at the same place
            if (!graph_1->is_empty(path_handle_1)) {
                
                auto check_equiv = [&](const occurrence_handle_t& occ_1,
                                       const occurrence_handle_t& occ_2) {
                    handle_t handle_1 = graph_1->get_occurrence(occ_1);
                    handle_t handle_2 = graph_2->get_occurrence(occ_2);
                    bool match = (graph_1->get_id(handle_1) == graph_2->get_id(handle_2) &&
                                  graph_1->get_is_reverse(handle_1) == graph_2->get_is_reverse(handle_2));
                    if (verbose && !match) {
                        cerr << "path " << graph_1->get_path_name(path_handle_1) << " has mismatching occurrences " << graph_1->get_id(handle_1) << (graph_1->get_is_reverse(handle_1) ? "-" : "+") << " and " << graph_2->get_id(handle_2) << (graph_2->get_is_reverse(handle_2) ? "-" : "+") << endl;
                    }
                    
                    return match;
                };
                
                occurrence_handle_t occ_1 = graph_1->get_first_occurrence(path_handle_1);
                occurrence_handle_t occ_2 = graph_2->get_first_occurrence(path_handle_2);
                
                if (!check_equiv(occ_1, occ_2)) {
                    equivalent = false;
                    return false;
                }
                
                while (graph_1->has_next_occurrence(occ_1)) {
                    occ_1 = graph_1->get_next_occurrence(occ_1);
                    occ_2 = graph_1->get_next_occurrence(occ_2);
                    
                    if (!check_equiv(occ_1, occ_2)) {
                        equivalent = false;
                        return false;
                    }
                }
            }
            
            return true;
        });
        
        
        return equivalent;
    }


}
}
