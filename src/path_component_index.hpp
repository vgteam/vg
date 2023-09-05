#ifndef VG_PATH_COMPONENT_INDEX_HPP_INCLUDED
#define VG_PATH_COMPONENT_INDEX_HPP_INCLUDED

/** \file
 *
 * Contains an index that maps embedded paths to the connected components of a graph
 */

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "handle.hpp"

namespace vg {
    
    using namespace std;
    
    /*
     * A class that can keep track of which embedded paths are on which
     * component of the graph.
     */
    class PathComponentIndex {
    public:
        
        /// Constructor
        PathComponentIndex(const PathHandleGraph* graph);
        
        /// Teturns true if the paths are on the same connected component of the graph
        bool paths_on_same_component(const path_handle_t& path_1,
                                     const path_handle_t& path_2) const;
        
        
        
    private:
        
        /// We make the default constructor private so that it can be used
        /// in move's, etc. but isn't exposed
        PathComponentIndex();
        
        /// Memoized sets of the paths that co-occur on a connected component
        vector<unordered_set<path_handle_t>> component_path_sets;
        
        /// An index from a path  to the set of paths that occur on the same
        /// connected component as it
        unordered_map<path_handle_t, size_t> component_path_set_of_path;
    };
}

#endif
