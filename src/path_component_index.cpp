#include "path_component_index.hpp"

#include <queue>
#include "sdsl/bit_vectors.hpp"
#include "algorithms/component_paths.hpp"

//#define debug_component_index

namespace vg {

    PathComponentIndex::PathComponentIndex() {
        // Nothing to do
    }
    
    PathComponentIndex::PathComponentIndex(const PathHandleGraph* graph) {
        
        component_path_sets = algorithms::component_paths(*graph);
        
        // make it so we can index into this with the path rank directly
        component_path_set_of_path.reserve(graph->get_path_count());
        
        // index from the paths to their component set
        for (size_t i = 0; i < component_path_sets.size(); i++) {
            for (const path_handle_t& path : component_path_sets[i]) {
                if (component_path_set_of_path.count(path)) {
                    cerr << "warning:[PathComponentIndex] Graph contains path " << graph->get_path_name(path) << " that spans multiple connected components. This path must follow edges that are not included in the graph. The PathComponentIndex may not be semantically meaningful for this graph." << endl;
                    continue;
                }
                component_path_set_of_path[path] = i;
            }
        }
    }
    
    bool PathComponentIndex::paths_on_same_component(const path_handle_t& path_1,
                                                     const path_handle_t& path_2) const {
        
        return component_path_sets.at(component_path_set_of_path.at(path_1)).count(path_2);
    }
}
