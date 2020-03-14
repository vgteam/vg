//
//  snarls.cpp
//
//  Unit tests for TraversalSupportFinder
//

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <set>
#include "json2pb.h"
#include <vg/vg.pb.h>
#include "catch.hpp"
#include "traversal_support.hpp"
#include "traversal_finder.hpp"
#include "../handle.hpp"
#include "../json2pb.h"
#include "../proto_handle_graph.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>

//#define debug

namespace vg {
namespace unittest {


/**
 * Get the read support from some maps
 */ 
class TestTraversalSupportFinder : public TraversalSupportFinder {
public:
    TestTraversalSupportFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                               const unordered_map<nid_t, double>& node_supports,
                               const unordered_map<edge_t, double>& edge_supports) :
        TraversalSupportFinder(graph, snarl_manager),
        node_supports(node_supports),
        edge_supports(edge_supports) {
    }    
    ~TestTraversalSupportFinder() = default;

    Support get_edge_support(const edge_t& edge) const {
        Support s;
        s.set_forward(edge_supports.at(edge));
        return s;
    }
    Support get_edge_support(id_t from, bool from_reverse, id_t to, bool to_reverse) const {
        return get_edge_support(graph.edge_handle(graph.get_handle(from, from_reverse),
                                                  graph.get_handle(to, to_reverse)));
    }
    virtual Support get_min_node_support(id_t node) const {
        Support s;
        s.set_forward(node_supports.at(node));
        return s;
    }
    virtual Support get_avg_node_support(id_t node) const {
        return get_min_node_support(node);
    }
    
protected:    
    const unordered_map<nid_t, double> node_supports;
    const unordered_map<edge_t, double> edge_supports;
};

TEST_CASE( "Deletion allele supports found correctly",
           "[traversal_support]" ) {

        string graph_json = R"(
{"edge": [{"from": "31041", "to": "31042"}, {"from": "31040", "to": "31041"}, {"from": "31040", "to": "31043"}, {"from": "134035", "to": "148994"}, {"from": "31042", "to": "134035"}, {"from": "31043", "from_start": true, "to": "134035", "to_end": true}, {"from": "31043", "from_start": true, "to": "148994", "to_end": true}], "node": [{"id": "31041", "sequence": "TATTTCCTAATGGGGTAGTGTCAGAGAGAGTA"}, {"id": "31040", "sequence": "GGCCCTGGAATATC"}, {"id": "134035", "sequence": "ATC"}, {"id": "31042", "sequence": "ATAACGCAGTATTTGTGA"}, {"id": "148994", "sequence": "A"}, {"id": "31043", "sequence": "GATCCCCTCTCCTTTACGAACTGGTAGAAGTG"}]}
    )";
    
    Graph g;
    json2pb(g, graph_json);
    
    // Wrap the graph in a HandleGraph
    ProtoHandleGraph graph(&g);

    unordered_map<nid_t, double> node_supports = {
        {31040, 17.5},
        {31041, 17.3438},
        {31042, 21.2778},
        {31043, 24.3125},
        {134035, 23.6667},
        {148994, 2}
    };

    unordered_map<edge_t, double> edge_supports = {
        {graph.edge_handle(graph.get_handle(31040, false), graph.get_handle(31041, false)), 17},
        {graph.edge_handle(graph.get_handle(31040, false), graph.get_handle(31043, false)), 1},
        {graph.edge_handle(graph.get_handle(31041, false), graph.get_handle(31042, false)), 18},
        {graph.edge_handle(graph.get_handle(31042, false), graph.get_handle(134035, false)), 24},
        {graph.edge_handle(graph.get_handle(134035, false), graph.get_handle(148994, false)), 2},
        {graph.edge_handle(graph.get_handle(134035, false), graph.get_handle(31043, false)), 23},
        {graph.edge_handle(graph.get_handle(148994, false), graph.get_handle(31043, false)), 2},
    };

    Snarl snarl;
    snarl.mutable_start()->set_node_id(31040);
    snarl.mutable_start()->set_backward(false);
    snarl.mutable_end()->set_node_id(31043);
    snarl.mutable_end()->set_backward(false);

    SnarlManager snarl_manager;
    TestTraversalSupportFinder support_finder(graph, snarl_manager, node_supports, edge_supports);
    // work with minimum supports (like call_main does when *not* genotyping)
    support_finder.set_support_switch_threshold(numeric_limits<size_t>::max(), 50);

    // test out our traversal finder while we're at it
    SnarlTraversal trav_0;
    json2pb(trav_0, R"({"visit": [{"node_id": "31040"}, {"node_id": "31041"}, {"node_id": "31042"}, {"node_id": "134035"}, {"node_id": "31043"}]})");
    SnarlTraversal trav_1;
    json2pb(trav_1, R"({"visit": [{"node_id": "31040"}, {"node_id": "31041"}, {"node_id": "31042"}, {"node_id": "134035"}, {"node_id": "148994"}, {"node_id": "31043"}]})");
    SnarlTraversal trav_2;
    json2pb(trav_2, R"({"visit": [{"node_id": "31040"}, {"node_id": "31043"}]})");
    
    function<double(handle_t)> node_weight = [&](handle_t h) {
        return node_supports.at(graph.get_id(h));
    };
    function<double(edge_t)> edge_weight = [&](edge_t e) {
        return edge_supports.at(e);
    };
    FlowTraversalFinder trav_finder(graph, snarl_manager, 100, node_weight, edge_weight);
    pair<vector<SnarlTraversal>, vector<double>> weighted_travs = trav_finder.find_weighted_traversals(snarl);
    REQUIRE(weighted_travs.first.size() == 3);
    REQUIRE(weighted_travs.first[0] == trav_0);
    REQUIRE(weighted_travs.first[1] == trav_1);
    REQUIRE(weighted_travs.first[2] == trav_2);

    REQUIRE(weighted_travs.second[0] == 17);
    REQUIRE(weighted_travs.second[1] == 2);
    REQUIRE(weighted_travs.second[2] == 1);

    const vector<SnarlTraversal>& travs = weighted_travs.first;

    // get the indidivual supports
    vector<Support> ind_supports = support_finder.get_traversal_set_support(travs, {}, {}, {}, false, {}, {}, 0);

    REQUIRE(TraversalSupportFinder::support_val(ind_supports[0]) == weighted_travs.second[0]);
    REQUIRE(TraversalSupportFinder::support_val(ind_supports[1]) == weighted_travs.second[1]);
    REQUIRE(TraversalSupportFinder::support_val(ind_supports[2]) == weighted_travs.second[2]);

    // get the exclusive support of trav 1.
    vector<Support> exc_supports1 = support_finder.get_traversal_set_support(travs, {0}, {}, {1}, true, {}, {}, 0);

    REQUIRE(TraversalSupportFinder::support_val(ind_supports[1]) == 2);

    // get the shared support of trav 0 and 1
    vector<Support> shared_supports = support_finder.get_traversal_set_support(travs, {0,1}, {}, {}, false, {}, {}, 0);

    REQUIRE(TraversalSupportFinder::support_val(shared_supports[0]) == 0.5 * weighted_travs.second[0]);
    REQUIRE(TraversalSupportFinder::support_val(shared_supports[1]) == weighted_travs.second[1]);
    REQUIRE(TraversalSupportFinder::support_val(shared_supports[2]) == weighted_travs.second[2]);

    // get the shared support of trav 0 and 1 but weighting by their relative scores
    vector<Support> rel_shared_supports = support_finder.get_traversal_set_support(travs, {0,1}, ind_supports, {}, false, {}, {}, 0);

    REQUIRE(TraversalSupportFinder::support_val(rel_shared_supports[0]) == (17./19.) * weighted_travs.second[0]);
    REQUIRE(TraversalSupportFinder::support_val(rel_shared_supports[1]) == (2./19.) * weighted_travs.second[0]);
    REQUIRE(TraversalSupportFinder::support_val(rel_shared_supports[2]) == weighted_travs.second[2]);

    
}

}
}
