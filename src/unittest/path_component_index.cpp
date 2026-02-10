///
/// \file path_component_index.cpp
///  
/// Unit tests for PathComponentIndex
///

#include "catch.hpp"
#include "path_component_index.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "../io/json2graph.hpp"
#include <bdsg/hash_graph.hpp>
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {
    
    TEST_CASE("Path component memoization produces expected results", "[pathcomponent]") {
        
        string graph_json = R"({"node": [{"sequence": "AAACCC", "id": 1}, {"sequence": "CACACA", "id": 2}, {"sequence": "CACACA", "id": 3}, {"sequence": "TTTTGG", "id": 4}, {"sequence": "ACGTAC", "id": 5}], "path": [{"name": "one", "mapping": [{"position": {"node_id": 1}, "rank": 1}, {"position": {"node_id": 2}, "rank": 2}]}, {"name": "three", "mapping": [{"position": {"node_id": 2}, "rank": 1}, {"position": {"node_id": 3}, "rank": 2}]}, {"name": "two", "mapping": [{"position": {"node_id": 4}, "rank": 1}, {"position": {"node_id": 5}, "rank": 2}]}], "edge": [{"from": 1, "to": 2}, {"from": 2, "to": 3}, {"from": 4, "to": 5}]})";

        // Load the JSON
        bdsg::HashGraph graph;
        vg::io::json2graph(graph_json, &graph);

        // Build the xg index
        xg::XG xg_index;
        xg_index.from_path_handle_graph(graph);
        
        
        unordered_set<path_handle_t> comp_1;
        comp_1.insert(xg_index.get_path_handle("one"));
        comp_1.insert(xg_index.get_path_handle("three"));
        
        PathComponentIndex pc_index(&xg_index);
        
        xg_index.for_each_path_handle([&](const path_handle_t& path_1) {
            xg_index.for_each_path_handle([&](const path_handle_t& path_2) {
                REQUIRE(pc_index.paths_on_same_component(path_1, path_2)
                        == (comp_1.count(path_1) == comp_1.count(path_2)));
            });
        });
    }

}
}
