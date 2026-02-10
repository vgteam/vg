/**
 * unittest/vg.cpp: test cases for vg::VG methods
 */

#include "catch.hpp"
#include "../vg.hpp"
#include "../augment.hpp"
#include "../utility.hpp"
#include "../algorithms/normalize.hpp"
#include "../algorithms/disjoint_components.hpp"
#include "handle.hpp"

namespace vg {
namespace unittest {

using namespace std;

// Turn a JSON string into a VG graph
VG string_to_graph(const string& json) {
    VG graph;
    Graph chunk;
    json2pb(chunk, json.c_str(), json.size());
    graph.merge(chunk);
    
    return graph;
}

TEST_CASE("dagify() should render the graph acyclic", "[vg][cycles][dagify]") {
   
    unordered_map<nid_t, pair<nid_t, bool> > node_translation;
   
    SECTION("a tiny DAG should remain unmodified") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        VG dag = graph.dagify(5, node_translation, 5, 0);
        
        REQUIRE(handlealgs::is_acyclic(&dag));
        REQUIRE(dag.get_node_count() == 2);
        REQUIRE(dag.edge_count() == 1);
    }
    
    SECTION("a tiny cyclic graph should become acyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 2, "to": 1}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        VG dag = graph.dagify(5, node_translation, 5, 0);
        
        REQUIRE(handlealgs::is_acyclic(&dag));
        REQUIRE(dag.get_node_count() >= 2);
    }
    
    SECTION("a tiny cyclic graph with doubly reversing edges should become acyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "from_start": true, "to": 2, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        VG dag = graph.dagify(5, node_translation, 5, 0);
        
        REQUIRE(handlealgs::is_acyclic(&dag));
        REQUIRE(dag.get_node_count() >= 2);
    }
    
}

TEST_CASE("unfold() should properly unfold a graph out to the requested length", "[vg][unfold]") {

    SECTION("Unfolding a graph with no reversing edges should create an isomorphic graph") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2},
                     {"sequence": "TGA","id": 3},
                     {"sequence": "GGC","id": 4}
                     ],
            "edge": [
                     {"from": 1,"to": 2},
                     {"from": 1,"to": 3},
                     {"from": 2,"to": 4},
                     {"from": 3,"to": 4}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(10000, node_translation);
        
        Graph& g = unfolded.graph;
        
        REQUIRE(g.node_size() == 4);
        REQUIRE(g.edge_size() == 4);
        
        bool found_node_1 = false;
        bool found_node_2 = false;
        bool found_node_3 = false;
        bool found_node_4 = false;
        
        for (int i = 0; i < g.node_size(); i++) {
            const Node& n = g.node(i);
            nid_t orig_id = node_translation[n.id()].first;
            if (orig_id == 1) {
                found_node_1 = true;
            }
            else if (orig_id == 2) {
                found_node_2 = true;
            }
            else if (orig_id == 3) {
                found_node_3 = true;
            }
            else if (orig_id == 4) {
                found_node_4 = true;
            }
        }
        
        REQUIRE(found_node_1);
        REQUIRE(found_node_2);
        REQUIRE(found_node_3);
        REQUIRE(found_node_4);
        
        bool found_edge_1 = false;
        bool found_edge_2 = false;
        bool found_edge_3 = false;
        bool found_edge_4 = false;
        
        for (int i = 0; i < g.edge_size(); i++) {
            const Edge& e = g.edge(i);
            nid_t from = node_translation[e.from()].first;
            nid_t to = node_translation[e.to()].first;
            Node& orig_from_node = *graph.get_node(from);
            Node& orig_to_node = *graph.get_node(to);
            Node& unfold_from_node = *unfolded.get_node(e.from());
            Node& unfold_to_node = *unfolded.get_node(e.to());
            if ((from == 1 && to == 2 && unfold_from_node.sequence() == orig_from_node.sequence()
                 && orig_to_node.sequence() == unfold_to_node.sequence()
                 && (!e.from_start() && !e.to_end())) ||
                (from == 2 && to == 1 && unfold_from_node.sequence() == orig_from_node.sequence()
                 && orig_to_node.sequence() == unfold_to_node.sequence()
                 && (e.from_start() && e.to_end())) ||
                (from == 2 && to == 1 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                 && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                 && (!e.from_start() && !e.to_end())) ||
                (from == 1 && to == 2 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                 && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                 && (e.from_start() && e.to_end()))) {
                found_edge_1 = true;
            }
            else if ((from == 1 && to == 3 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 3 && to == 1 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (e.from_start() && e.to_end())) ||
                     (from == 3 && to == 1 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 1 && to == 3 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (e.from_start() && e.to_end()))) {
                found_edge_2 = true;
            }
            else if ((from == 2 && to == 4 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 4 && to == 2 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (e.from_start() && e.to_end())) ||
                     (from == 4 && to == 2 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 2 && to == 4 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (e.from_start() && e.to_end()))) {
                found_edge_3 = true;
            }
            else if ((from == 3 && to == 4 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 4 && to == 3 && unfold_from_node.sequence() == orig_from_node.sequence()
                      && orig_to_node.sequence() == unfold_to_node.sequence()
                      && (e.from_start() && e.to_end())) ||
                     (from == 4 && to == 3 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (!e.from_start() && !e.to_end())) ||
                     (from == 3 && to == 4 && unfold_from_node.sequence() == reverse_complement(orig_from_node.sequence())
                      && orig_to_node.sequence() == reverse_complement(unfold_to_node.sequence())
                      && (e.from_start() && e.to_end()))) {
                found_edge_4 = true;
            }
        }
        
        REQUIRE(found_edge_1);
        REQUIRE(found_edge_2);
        REQUIRE(found_edge_3);
        REQUIRE(found_edge_4);
    }
    
    SECTION("Unfolding can flip the reversed portion of a non-branching path with reversing edges") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2},
                     {"sequence": "TGA","id": 3}
                     ],
            "edge": [
                     {"from": 1,"to": 2,"to_end": true},
                     {"from": 2,"to": 3,"from_start": true}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(10000, node_translation);
        
        Graph& g = unfolded.graph;
        
        REQUIRE(g.node_size() == 3);
        REQUIRE(g.edge_size() == 2);
        
        bool in_orientation_1 = true;
        bool in_orientation_2 = true;
        for (auto record : node_translation) {
            if (record.second.first == 1) {
                in_orientation_1 = in_orientation_1 && !record.second.second;
                in_orientation_2 = in_orientation_2 && record.second.second;
            }
            else if (record.second.first == 2) {
                in_orientation_1 = in_orientation_1 && record.second.second;
                in_orientation_2 = in_orientation_2 && !record.second.second;
            }
            else if (record.second.first == 3) {
                in_orientation_1 = in_orientation_1 && !record.second.second;
                in_orientation_2 = in_orientation_2 && record.second.second;
            }
        }
        REQUIRE(in_orientation_1 != in_orientation_2);
        
        if (in_orientation_1) {
            for (int i = 0; i < g.node_size(); i++) {
                const Node& unfold_node = g.node(i);
                const Node& orig_node = *graph.get_node(node_translation[unfold_node.id()].first);
                if (orig_node.id() == 1) {
                    REQUIRE(unfold_node.sequence() == orig_node.sequence());
                }
                else if (orig_node.id() == 2) {
                    REQUIRE(unfold_node.sequence() == reverse_complement(orig_node.sequence()));
                }
                else if (orig_node.id() == 3) {
                    REQUIRE(unfold_node.sequence() == orig_node.sequence());
                }
            }
        }
        else {
            for (int i = 0; i < g.node_size(); i++) {
                const Node& unfold_node = g.node(i);
                const Node& orig_node = *graph.get_node(node_translation[unfold_node.id()].first);
                if (orig_node.id() == 1) {
                    REQUIRE(unfold_node.sequence() == reverse_complement(orig_node.sequence()));
                }
                else if (orig_node.id() == 2) {
                    REQUIRE(unfold_node.sequence() == orig_node.sequence());
                }
                else if (orig_node.id() == 3) {
                    REQUIRE(unfold_node.sequence() == reverse_complement(orig_node.sequence()));
                }
            }
        }
    }
    
    SECTION("Unfolding can turn a reversing cycle into a directed cycle") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2}
                     ],
            "edge": [
                     {"from": 1,"to": 2},
                     {"from": 2,"to": 2,"to_end": true},
                     {"from": 1,"to": 1,"from_start": true}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(10000, node_translation);
        
        Graph& g = unfolded.graph;
        
        REQUIRE(g.node_size() == 4);
        REQUIRE(g.edge_size() == 4);
        
        nid_t node_1 = 0;
        nid_t node_2 = 0;
        nid_t node_3 = 0;
        nid_t node_4 = 0;
        
        for (int i = 0; i < g.node_size(); i++) {
            const Node& n = g.node(i);
            nid_t orig_id = node_translation[n.id()].first;
            bool flipped =  node_translation[n.id()].second;
            if (orig_id == 1 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_1 = n.id();
            }
            else if (orig_id == 1 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_2 = n.id();
            }
            else if (orig_id == 2 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_3 = n.id();
            }
            else if (orig_id == 2 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_4 = n.id();
            }
        }
        
        REQUIRE(node_1 != 0);
        REQUIRE(node_2 != 0);
        REQUIRE(node_3 != 0);
        REQUIRE(node_4 != 0);
        
        bool found_edge_1 = false;
        bool found_edge_2 = false;
        bool found_edge_3 = false;
        bool found_edge_4 = false;
        
        for (int i = 0; i < g.edge_size(); i++) {
            const Edge& e = g.edge(i);
            if ((e.from() == node_1 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                (e.from() == node_3 && e.to() == node_1 && e.from_start() && e.to_end())) {
                    found_edge_1 = true;
                }
            else if ((e.from() == node_3 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_4 && e.to() == node_3 && e.from_start() && e.to_end())) {
                         found_edge_2 = true;
                     }
            else if ((e.from() == node_4 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                    (e.from() == node_2 && e.to() == node_4 && e.from_start() && e.to_end())) {
                         found_edge_3 = true;
                     }
            else if ((e.from() == node_2 && e.to() == node_1 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_1 && e.to() == node_2 && e.from_start() && e.to_end())) {
                found_edge_4 = true;
            }
        }
        
        REQUIRE(found_edge_1);
        REQUIRE(found_edge_2);
        REQUIRE(found_edge_3);
        REQUIRE(found_edge_4);
    }
    
    SECTION("Unfolding can find reverse strand nodes that require traversing the same node in opposite directions") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2},
                     {"sequence": "GG","id": 3},
                     {"sequence": "TGC","id": 4},
                     {"sequence": "T","id": 5}
                     ],
            "edge": [
                     {"from": 1,"to": 3},
                     {"from": 2,"to": 3},
                     {"from": 3,"to": 4},
                     {"from": 3,"to": 5},
                     {"from": 2,"to": 2, "from_start": true},
                     {"from": 4,"to": 4, "to_end": true}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(10000, node_translation);
        
        Graph& g = unfolded.graph;
        
        REQUIRE(g.node_size() == 10);
        REQUIRE(g.edge_size() == 10);
        
        nid_t node_1 = 0;
        nid_t node_2 = 0;
        nid_t node_3 = 0;
        nid_t node_4 = 0;
        nid_t node_5 = 0;
        nid_t node_6 = 0;
        nid_t node_7 = 0;
        nid_t node_8 = 0;
        nid_t node_9 = 0;
        nid_t node_10 = 0;
        
        for (int i = 0; i < g.node_size(); i++) {
            const Node& n = g.node(i);
            nid_t orig_id = node_translation[n.id()].first;
            bool flipped =  node_translation[n.id()].second;
            if (orig_id == 1 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_1 = n.id();
            }
            else if (orig_id == 1 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_2 = n.id();
            }
            else if (orig_id == 2 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_3 = n.id();
            }
            else if (orig_id == 2 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_4 = n.id();
            }
            else if (orig_id == 3 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_5 = n.id();
            }
            else if (orig_id == 3 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_6 = n.id();
            }
            else if (orig_id == 4 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_7 = n.id();
            }
            else if (orig_id == 4 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_8 = n.id();
            }
            else if (orig_id == 5 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_9 = n.id();
            }
            else if (orig_id == 5 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_10 = n.id();
            }
        }
        
        REQUIRE(node_1 != 0);
        REQUIRE(node_2 != 0);
        REQUIRE(node_3 != 0);
        REQUIRE(node_4 != 0);
        REQUIRE(node_5 != 0);
        REQUIRE(node_6 != 0);
        REQUIRE(node_7 != 0);
        REQUIRE(node_8 != 0);
        REQUIRE(node_9 != 0);
        REQUIRE(node_10 != 0);
        
        
        bool found_edge_1 = false;
        bool found_edge_2 = false;
        bool found_edge_3 = false;
        bool found_edge_4 = false;
        bool found_edge_5 = false;
        bool found_edge_6 = false;
        bool found_edge_7 = false;
        bool found_edge_8 = false;
        bool found_edge_9 = false;
        bool found_edge_10 = false;
        
        for (int i = 0; i < g.edge_size(); i++) {
            const Edge& e = g.edge(i);
            if ((e.from() == node_1 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                (e.from() == node_5 && e.to() == node_1 && e.from_start() && e.to_end())) {
                found_edge_1 = true;
            }
            else if ((e.from() == node_3 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_5 && e.to() == node_3 && e.from_start() && e.to_end())) {
                found_edge_2 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_7 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_7 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_3 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_9 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_9 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_4 = true;
            }
            else if ((e.from() == node_7 && e.to() == node_8 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_8 && e.to() == node_7 && e.from_start() && e.to_end())) {
                found_edge_5 = true;
            }
            else if ((e.from() == node_8 && e.to() == node_6 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_6 && e.to() == node_8 && e.from_start() && e.to_end())) {
                found_edge_6 = true;
            }
            else if ((e.from() == node_10 && e.to() == node_6 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_6 && e.to() == node_10 && e.from_start() && e.to_end())) {
                found_edge_7 = true;
            }
            else if ((e.from() == node_6 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_2 && e.to() == node_6 && e.from_start() && e.to_end())) {
                found_edge_8 = true;
            }
            else if ((e.from() == node_6 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_4 && e.to() == node_6 && e.from_start() && e.to_end())) {
                found_edge_9 = true;
            }
            else if ((e.from() == node_4 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_3 && e.to() == node_4 && e.from_start() && e.to_end())) {
                found_edge_10 = true;
            }
        }
        
        REQUIRE(found_edge_1);
        REQUIRE(found_edge_2);
        REQUIRE(found_edge_3);
        REQUIRE(found_edge_4);
        REQUIRE(found_edge_5);
        REQUIRE(found_edge_6);
        REQUIRE(found_edge_7);
        REQUIRE(found_edge_8);
        REQUIRE(found_edge_9);
        REQUIRE(found_edge_10);
    }
    
    SECTION("Unfolding correctly handles a reversing path along a reverse-oriented path") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2},
                     {"sequence": "GG","id": 3},
                     {"sequence": "TGC","id": 4},
                     {"sequence": "T","id": 5}
                     ],
            "edge": [
                     {"from": 1,"to": 2, "to_end": true},
                     {"from": 2,"to": 3, "from_start": true, "to_end": true},
                     {"from": 3,"to": 4, "from_start": true, "to_end": true},
                     {"from": 4,"to": 5, "from_start": true},
                     {"from": 3,"to": 2, "from_start": true},
                     {"from": 4,"to": 3, "to_end": true}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(10000, node_translation);
        
        Graph& g = unfolded.graph;
        
        REQUIRE(g.node_size() == 10);
        REQUIRE(g.edge_size() == 12);
        
        nid_t node_1 = 0;
        nid_t node_2 = 0;
        nid_t node_3 = 0;
        nid_t node_4 = 0;
        nid_t node_5 = 0;
        nid_t node_6 = 0;
        nid_t node_7 = 0;
        nid_t node_8 = 0;
        nid_t node_9 = 0;
        nid_t node_10 = 0;
        
        for (int i = 0; i < g.node_size(); i++) {
            const Node& n = g.node(i);
            nid_t orig_id = node_translation[n.id()].first;
            bool flipped =  node_translation[n.id()].second;
            if (orig_id == 1 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_1 = n.id();
            }
            else if (orig_id == 1 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_2 = n.id();
            }
            else if (orig_id == 2 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_3 = n.id();
            }
            else if (orig_id == 2 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_4 = n.id();
            }
            else if (orig_id == 3 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_5 = n.id();
            }
            else if (orig_id == 3 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_10 = n.id();
            }
            else if (orig_id == 4 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_6 = n.id();
            }
            else if (orig_id == 4 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_7 = n.id();
            }
            else if (orig_id == 5 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_8 = n.id();
            }
            else if (orig_id == 5 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_9 = n.id();
            }
        }
        
        REQUIRE(node_1 != 0);
        REQUIRE(node_2 != 0);
        REQUIRE(node_3 != 0);
        REQUIRE(node_4 != 0);
        REQUIRE(node_5 != 0);
        REQUIRE(node_6 != 0);
        REQUIRE(node_7 != 0);
        REQUIRE(node_8 != 0);
        REQUIRE(node_9 != 0);
        REQUIRE(node_10 != 0);
        
        bool found_edge_1 = false;
        bool found_edge_2 = false;
        bool found_edge_3 = false;
        bool found_edge_4 = false;
        bool found_edge_5 = false;
        bool found_edge_6 = false;
        bool found_edge_7 = false;
        bool found_edge_8 = false;
        bool found_edge_9 = false;
        bool found_edge_10 = false;
        bool found_edge_11 = false;
        bool found_edge_12 = false;
        
        for (int i = 0; i < g.edge_size(); i++) {
            const Edge& e = g.edge(i);
            if ((e.from() == node_1 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                (e.from() == node_4 && e.to() == node_1 && e.from_start() && e.to_end())) {
                found_edge_1 = true;
            }
            else if ((e.from() == node_4 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_5 && e.to() == node_4 && e.from_start() && e.to_end())) {
                found_edge_2 = true;
            }
            else if ((e.from() == node_9 && e.to() == node_6 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_6 && e.to() == node_9 && e.from_start() && e.to_end())) {
                found_edge_3 = true;
            }
            else if ((e.from() == node_6 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_5 && e.to() == node_6 && e.from_start() && e.to_end())) {
                found_edge_4 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_3 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_5 = true;
            }
            else if ((e.from() == node_3 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_2 && e.to() == node_3 && e.from_start() && e.to_end())) {
                found_edge_6 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_7 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_7 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_7 = true;
            }
            else if ((e.from() == node_7 && e.to() == node_8 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_8 && e.to() == node_7 && e.from_start() && e.to_end())) {
                found_edge_8 = true;
            }
            else if ((e.from() == node_4 && e.to() == node_10 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_10 && e.to() == node_4 && e.from_start() && e.to_end())) {
                found_edge_9 = true;
            }
            else if ((e.from() == node_6 && e.to() == node_10 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_10 && e.to() == node_6 && e.from_start() && e.to_end())) {
                found_edge_10 = true;
            }
            else if ((e.from() == node_10 && e.to() == node_7 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_7 && e.to() == node_10 && e.from_start() && e.to_end())) {
                found_edge_11 = true;
            }
            else if ((e.from() == node_10 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_3 && e.to() == node_10 && e.from_start() && e.to_end())) {
                found_edge_12 = true;
            }
        }
        
        REQUIRE(found_edge_1);
        REQUIRE(found_edge_2);
        REQUIRE(found_edge_3);
        REQUIRE(found_edge_4);
        REQUIRE(found_edge_5);
        REQUIRE(found_edge_6);
        REQUIRE(found_edge_7);
        REQUIRE(found_edge_8);
        REQUIRE(found_edge_9);
        REQUIRE(found_edge_10);
        REQUIRE(found_edge_11);
        REQUIRE(found_edge_12);
    }
    
    SECTION("Unfolding does not duplicate past the length limit") {
        const string graph_json = R"(
        {
            "node": [
                     {"sequence": "ATA","id": 1},
                     {"sequence": "CT","id": 2},
                     {"sequence": "GG","id": 3},
                     {"sequence": "TA","id": 4},
                     {"sequence": "ACT","id": 5}
                     ],
            "edge": [
                     {"from": 1,"to": 2},
                     {"from": 2,"to": 3},
                     {"from": 2,"to": 3, "to_end": true},
                     {"from": 3,"to": 4},
                     {"from": 3,"to": 4, "from_start": true},
                     {"from": 4,"to": 5}
                     ]
        }
        )";
        
        VG graph = string_to_graph(graph_json);
        
        unordered_map<nid_t, pair<nid_t, bool> > node_translation;
        VG unfolded = graph.unfold(2, node_translation);
        
        Graph& g = unfolded.graph;
                
        REQUIRE(g.node_size() == 8);
        REQUIRE(g.edge_size() == 8);
        
        nid_t node_1 = 0;
        nid_t node_2 = 0;
        nid_t node_3 = 0;
        nid_t node_4 = 0;
        nid_t node_5 = 0;
        nid_t node_6 = 0;
        nid_t node_7 = 0;
        nid_t node_8 = 0;
        
        for (int i = 0; i < g.node_size(); i++) {
            const Node& n = g.node(i);
            nid_t orig_id = node_translation[n.id()].first;
            bool flipped =  node_translation[n.id()].second;
            if (orig_id == 1 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_1 = n.id();
            }
            else if (orig_id == 2 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_2 = n.id();
            }
            else if (orig_id == 2 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_3 = n.id();
            }
            else if (orig_id == 3 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_4 = n.id();
            }
            else if (orig_id == 3 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_5 = n.id();
            }
            else if (orig_id == 4 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_6 = n.id();
            }
            else if (orig_id == 4 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                node_7 = n.id();
            }
            else if (orig_id == 5 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                node_8 = n.id();
            }
        }
        
        REQUIRE(node_1 != 0);
        REQUIRE(node_2 != 0);
        REQUIRE(node_3 != 0);
        REQUIRE(node_4 != 0);
        REQUIRE(node_5 != 0);
        REQUIRE(node_6 != 0);
        REQUIRE(node_7 != 0);
        REQUIRE(node_8 != 0);
        
        bool found_edge_1 = false;
        bool found_edge_2 = false;
        bool found_edge_3 = false;
        bool found_edge_4 = false;
        bool found_edge_5 = false;
        bool found_edge_6 = false;
        bool found_edge_7 = false;
        bool found_edge_8 = false;
        bool found_edge_9 = false;
        bool found_edge_10 = false;
        
        for (int i = 0; i < g.edge_size(); i++) {
            const Edge& e = g.edge(i);
            if ((e.from() == node_1 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                (e.from() == node_2 && e.to() == node_1 && e.from_start() && e.to_end())) {
                found_edge_1 = true;
            }
            else if ((e.from() == node_2 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_4 && e.to() == node_2 && e.from_start() && e.to_end())) {
                found_edge_2 = true;
            }
            else if ((e.from() == node_4 && e.to() == node_6 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_6 && e.to() == node_4 && e.from_start() && e.to_end())) {
                found_edge_3 = true;
            }
            else if ((e.from() == node_6 && e.to() == node_8 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_8 && e.to() == node_6 && e.from_start() && e.to_end())) {
                found_edge_4 = true;
            }
            else if ((e.from() == node_2 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_5 && e.to() == node_2 && e.from_start() && e.to_end())) {
                found_edge_5 = true;
            }
            else if ((e.from() == node_4 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_3 && e.to() == node_4 && e.from_start() && e.to_end())) {
                found_edge_6 = true;
            }
            else if ((e.from() == node_7 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_4 && e.to() == node_7 && e.from_start() && e.to_end())) {
                found_edge_7 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_6 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_6 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_8 = true;
            }
            else if ((e.from() == node_7 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_5 && e.to() == node_7 && e.from_start() && e.to_end())) {
                found_edge_9 = true;
            }
            else if ((e.from() == node_5 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                     (e.from() == node_3 && e.to() == node_5 && e.from_start() && e.to_end())) {
                found_edge_10 = true;
            }
        }
        
        REQUIRE(found_edge_1);
        REQUIRE(found_edge_2 != found_edge_9);
        REQUIRE(found_edge_3 != found_edge_10);
        REQUIRE(found_edge_4);
        REQUIRE(found_edge_5);
        REQUIRE(found_edge_6);
        REQUIRE(found_edge_7);
        REQUIRE(found_edge_8);
    }
}

TEST_CASE("expand_context_by_length() should respect barriers", "[vg][context]") {

    const string graph_json = R"(
    {
      "node": [
        {"sequence": "CCATTTGTCCAAAGT","id": 1},
        {"sequence": "AAGCAAACACTG","id": 2},
        {"sequence": "C","id": 3},
        {"sequence": "T","id": 4},
        {"sequence": "TACACTCTTGGAGGGAA","id": 5},
        {"sequence": "T","id": 6},
        {"sequence": "C","id": 7},
        {"sequence": "AAAAACTAG","id": 8},
        {"sequence": "AGTTGCAT","id": 9},
        {"sequence": "TTCTCTGATGATGAG","id": 10},
        {"sequence": "TGATGTTGAGGGTTTTTTTTGTCT","id": 11},
        {"sequence": "ATTGGTCACTTGTACATCTTATTTTTACAA","id": 12},
        {"sequence":"GAACGTTT", "id": 13}
      ],
      "edge": [
        {"from": 1,"to": 9,"from_start": true},
        {"from": 1,"to": 2},
        {"from": 2,"to": 3},
        {"from": 2,"to": 4},
        {"from": 3, "to": 5},
        {"from": 4,"to": 5},
        {"from": 5,"to": 6},
        {"from": 5,"to": 7},
        {"from": 6,"to": 8},
        {"from": 7,"to": 8},
        {"from": 9,"to": 10},
        {"from": 10,"to": 11},
        {"from": 11,"to": 12},
        {"from": 12,"to": 13}
      ]
    }
    )";
    
    VG graph = string_to_graph(graph_json);

    SECTION("barriers on either end of the seed node should stop anything being extracted") {

        VG context;
        context.add_node(*graph.get_node(3));
        graph.expand_context_by_length(context, 1000, false, true, {NodeSide(3, false), NodeSide(3, true)});
        
        REQUIRE(context.size() == 1);
    }
    
    SECTION("barriers should stop edges being formed") {

        VG context;
        context.add_node(*graph.get_node(3));
        context.add_node(*graph.get_node(4));
        // Note that we wouldn't get any edges between 3 and 4, if there were
        // any, because context expansion sees no edges between seed nodes.
        graph.expand_context_by_length(context, 1000, false, true, {NodeSide(3, false), NodeSide(3, true)});
        
        SECTION("node 4 should have both attached edges") {
            REQUIRE(context.has_edge(NodeSide(4, false), NodeSide(2, true)) == true);
            REQUIRE(context.has_edge(NodeSide(4, true), NodeSide(5, false)) == true);
        }
        
        SECTION("node 3 should have no atatched edges") {
            REQUIRE(context.has_edge(NodeSide(3, false), NodeSide(2, true)) == false);
            REQUIRE(context.has_edge(NodeSide(3, true), NodeSide(5, false)) == false);
        }
        
        
    }

}


TEST_CASE("add_nodes_and_edges() should connect all nodes", "[vg][edit]") {
    
    const string graph_json = R"(
        
    {
        "node": [
            {"id": 1, "sequence": "GATT"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "C"},
            {"id": 4, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 3, "to": 4}
        ]
    }

    )";
    
    // Define a graph
    VG graph = string_to_graph(graph_json);
    
    const string path_json = R"(
    {
        "mapping": [
            {
                "position": {
                    "node_id": 1
                },
                "edit": [
                    {
                        "from_length": 4,
                        "to_length": 4
                    }, 
                    {
                        "from_length": 0, 
                        "to_length": 10,
                        "sequence": "AAAAAAAAAA"
                    }
                ]
            },
            {
                "position": {
                    "node_id": 4
                },
                "edit": [
                    {
                        "from_length": 1, 
                        "to_length": 1
                    }
                ]
            }
        ]
    }
    )";
    
    // And a big insert
    Path path;
    json2pb(path, path_json.c_str(), path_json.size());
       
    // First prepare the various state things we need to pass
    
    // This can be empty if no changes have been made yet
    map<pos_t, nid_t> node_translation;
    // As can this
    unordered_map<pair<pos_t, string>, vector<nid_t>> added_seqs;
    // And this
    unordered_map<nid_t, Path> added_nodes;
    
    // This actually needs to be filled in
    unordered_map<nid_t, size_t> orig_node_sizes;
    graph.for_each_node([&](Node* node) {
        orig_node_sizes[node->id()] = node->sequence().size();
    });
    
    // And this can be empty if nothing is dangling in.
    set<NodeSide> dangling;
    
    // Do the addition, but limit node size.
    add_nodes_and_edges(&graph, path, node_translation, added_seqs, added_nodes, orig_node_sizes, dangling, 1);
    
    // Make sure it's still connected
    list<bdsg::HashGraph> subgraphs = algorithms::disjoint_components(graph);
    REQUIRE(subgraphs.size() == 1);
    
}

TEST_CASE("edit() should not get confused even under very confusing circumstances", "[vg][edit]") {
    
    const string graph_json = R"(
        
    {
        "node": [
            {"id": 1, "sequence": "GATT"},
            {"id": 2, "sequence": "T"},
            {"id": 3, "sequence": "C"},
            {"id": 4, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true},
            {"from": 1, "to": 3},
            {"from": 2, "to": 4, "from_start": true},
            {"from": 3, "to": 4}
        ]
    }

    )";
    
    // Define a graph
    VG graph = string_to_graph(graph_json);
    
    // And a path that doubles back on itself through an edge that isn't in the graph yet
    const string path_json = R"(
    {
        "mapping": [
            {
                "position": {
                    "node_id": 1,
                    "offset": 1
                },
                "edit": [
                    {
                        "from_length": 3,
                        "to_length": 3
                    }, 
                    {
                        "from_length": 0, 
                        "to_length": 3,
                        "sequence": "CCC"
                    }
                ]
            },
            {
                "position": {
                    "node_id": 2,
                    "is_reverse": true
                },
                "edit": [
                    {
                        "from_length": 1, 
                        "to_length": 1
                    }
                ]
            },
            {
                "position": {
                    "node_id": 2
                },
                "edit": [
                    {
                        "from_length": 1, 
                        "to_length": 1
                    }
                ]
            },
            {
                "position": {
                    "node_id": 1,
                    "is_reverse": true
                },
                "edit": [
                    {
                        "from_length": 1, 
                        "to_length": 1
                    }
                ]
            }
        ]
    }
    )";
    
    REQUIRE(graph.node_count() == 4);
    REQUIRE(graph.edge_count() == 4);
    
    Path path;
    json2pb(path, path_json.c_str(), path_json.size());
    
    // Needs to be in a vector to apply it
    vector<Path> paths{path};
       
    SECTION("edit() can add the path without modifying it") {
        graph.edit(paths, nullptr, false, false, false);
        
        REQUIRE(pb2json(paths.front()) == pb2json(path));
        
        // The graph should end up with 1 more node and 3 more edges.
        REQUIRE(graph.node_count() == 5);
        REQUIRE(graph.edge_count() == 7);
    }
    
    SECTION("edit() can add the path with modification only") {
        graph.edit(paths, nullptr, false, true, false);
        
        for (const auto& mapping : paths.front().mapping()) {
            // Make sure all the mappings are perfect matches
            REQUIRE(mapping_is_match(mapping));
        }
        
        // The graph should end up with 1 more node and 3 more edges.
        REQUIRE(graph.node_count() == 5);
        REQUIRE(graph.edge_count() == 7);
    }
    
    SECTION("edit() can add the path with end breaking only") {
        graph.edit(paths, nullptr, false, false, true);
        
        REQUIRE(pb2json(paths.front()) == pb2json(path));
        
        // The graph should end up with 3 more nodes (the insert plus 2 new
        // pieces of the original node 1) and 5 more edges.
        REQUIRE(graph.node_count() == 7);
        REQUIRE(graph.edge_count() == 9);
    }
    
}

TEST_CASE("reverse_complement_graph() produces expected results", "[vg]") {
    
    SECTION("reverse_complement_graph() works") {
        
        VG vg;
        
        Node* n0 = vg.create_node("AA");
        Node* n1 = vg.create_node("AC");
        Node* n2 = vg.create_node("AG");
        Node* n3 = vg.create_node("AT");
        Node* n4 = vg.create_node("CA");
        Node* n5 = vg.create_node("CC");
        Node* n6 = vg.create_node("CG");
        
        vg.create_edge(n1, n0, true, true);
        vg.create_edge(n0, n2, false, true);
        vg.create_edge(n2, n3, true, false);
        vg.create_edge(n3, n1, true, true);
        vg.create_edge(n4, n3, true, true);
        vg.create_edge(n0, n4);
        vg.create_edge(n4, n5, false, true);
        vg.create_edge(n0, n5, false, true);
        vg.create_edge(n6, n4, true, true);
        vg.create_edge(n3, n6);
      
        unordered_map<nid_t, pair<nid_t, bool>> trans;
        VG rev = vg.reverse_complement_graph(trans);
        
        REQUIRE(trans.size() == rev.graph.node_size());
        REQUIRE(rev.graph.node_size() == vg.graph.node_size());
        
        for (size_t i = 0; i < rev.graph.node_size(); i++) {
            const Node& node = rev.graph.node(i);
            Node* orig_node = vg.get_node(trans[node.id()].first);
            REQUIRE(reverse_complement(node.sequence()) == orig_node->sequence());
            
            vector<pair<nid_t, bool>> start_edges = vg.edges_start(orig_node);
            vector<pair<nid_t, bool>> end_edges = vg.edges_end(orig_node);
            
            vector<pair<nid_t, bool>> rev_start_edges = rev.edges_start(node.id());
            vector<pair<nid_t, bool>> rev_end_edges = rev.edges_end(node.id());
            
            REQUIRE(start_edges.size() == rev_end_edges.size());
            REQUIRE(end_edges.size() == rev_start_edges.size());
            
            for (auto re : rev_start_edges) {
                bool found = false;
                for (auto e : end_edges) {
                    found = found || re.first == e.first;
                }
                REQUIRE(found);
            }
            for (auto re : rev_end_edges) {
                bool found = false;
                for (auto e : start_edges) {
                    found = found || re.first == e.first;
                }
                REQUIRE(found);
            }
            
        }
    }
    
}

TEST_CASE("find_breakpoints() should determine where the graph needs to break to add a path", "[vg][edit]") {
    
    
    VG vg;
        
    Node* n1 = vg.create_node("GATT");
    Node* n2 = vg.create_node("AAAA");
    Node* n3 = vg.create_node("CA");
    
    vg.create_edge(n1, n2);
    vg.create_edge(n2, n3);
    
    // We will find breakpoints for a path
    Path path;
    // And store them here.
    unordered_map<nid_t, set<pos_t>> breakpoints;
    
    SECTION("find_breakpoints() works on a single edit perfect match") {
        
        // Set the path to a perfect match
        const string path_string = R"(
            {"mapping": [{"position": {"node_id": 1, "offset": 1}, "edit": [{"from_length": 2, "to_length": 2}]}]}
        )";
        json2pb(path, path_string.c_str(), path_string.size());
        
        SECTION("asking for breakpoints at the end gets us the two end breakpoints") {
            find_breakpoints(path, breakpoints);
            
            REQUIRE(breakpoints.size() == 1);
            REQUIRE(breakpoints.count(n1->id()) == 1);
            REQUIRE(breakpoints[n1->id()].size() == 2);
            REQUIRE(breakpoints[n1->id()].count(make_pos_t(n1->id(), false, 1)) == 1);
            REQUIRE(breakpoints[n1->id()].count(make_pos_t(n1->id(), false, 3)) == 1);
            
        
        }
        
        SECTION("asking for no breakpoints at the end gets us no breakpoints") {
            find_breakpoints(path, breakpoints, false);
            
            REQUIRE(breakpoints.empty());
        }
        
    }
    
}
TEST_CASE("create_handle() correctly creates handles using given sequence and id", "[vg]") {
    VG vg;
        
    handle_t h1 = vg.create_handle("GATT", 2);
    handle_t h2 = vg.create_handle("AAAA", 4);
    handle_t h3 = vg.create_handle("CA", 6);
    
    REQUIRE(vg.get_id(h1) == 2);
    REQUIRE(vg.get_id(h2) == 4);
    REQUIRE(vg.get_id(h3) == 6);
    REQUIRE(vg.get_sequence(h1) == "GATT");
    REQUIRE(vg.get_sequence(h2) == "AAAA");
    REQUIRE(vg.get_sequence(h3) == "CA");
    
}

TEST_CASE("normalize() can join nodes and merge siblings", "[vg][normalize]") {

    SECTION("Two redundant SNP values should be merged") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "T"},
                {"id": 3, "sequence": "T"},
                {"id": 4, "sequence": "G"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 1, "to": 4},
                {"from": 2, "to": 5},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // One of the two alternative Ts should have been eliminated
        REQUIRE(graph.get_node_count() == 4);
        
    }
    
    SECTION("Leading identical sequences should be condensed") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "TTA"},
                {"id": 3, "sequence": "TTC"},
                {"id": 4, "sequence": "GGG"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 1, "to": 4},
                {"from": 2, "to": 5},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // Those duplicate Ts should be eliminated
        REQUIRE(graph.length() == 13);
        
    }
    
    SECTION("Multiple families of leading identical sequences should be condensed") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "TTA"},
                {"id": 3, "sequence": "TTC"},
                {"id": 4, "sequence": "GGG"},
                {"id": 6, "sequence": "GGT"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 1, "to": 4},
                {"from": 1, "to": 6},
                {"from": 2, "to": 5},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5},
                {"from": 6, "to": 5}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // Those duplicate Ts and Gs should be eliminated
        REQUIRE(graph.length() == 14);
        
    }
    
    SECTION("Multiple families of trailing identical sequences should be condensed") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "ATT"},
                {"id": 3, "sequence": "CTT"},
                {"id": 4, "sequence": "GGG"},
                {"id": 6, "sequence": "TGG"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 1, "to": 4},
                {"from": 1, "to": 6},
                {"from": 2, "to": 5},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5},
                {"from": 6, "to": 5}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // Those duplicate Ts and Gs should be eliminated
        REQUIRE(graph.length() == 14);
        
    }
}

// TODO: This test case won't pass because VG::simplify_siblings() still thinks
// about "from" and "to" siblings, and doesn't really understand reversing
// edges. It doesn't see any siblings in these cases right now.
#ifdef test_reversing_siblings
TEST_CASE("normalize() can join nodes and merge siblings when nodes are backward", "[vg][normalize]") {

    SECTION("Leading identical sequences should be condensed even when nodes are backward") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "TAA"},
                {"id": 3, "sequence": "GAA"},
                {"id": 4, "sequence": "CCC"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2, "to_end": true},
                {"from": 1, "to": 3, "to_end": true},
                {"from": 1, "to": 4, "to_end": true},
                {"from": 2, "to": 5, "from_start": true},
                {"from": 3, "to": 5, "from_start": true},
                {"from": 4, "to": 5, "from_start": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // Those duplicate Ts (actually As) should be eliminated
        REQUIRE(graph.length() == 13);
        
        bool found = false;
        graph.for_each_node([&](const Node* n) {
            if (n->sequence() == "AA") {
                found = true;
            }
        });
        // They ought to have been combined into an AA node.
        REQUIRE(found);
        
    }

    SECTION("Leading and trailing identical sequences should be condensed correctly when nodes are backward and all siblings match") {
        // TODO: We don't support mixed orientations. But we should support all-backward.
    
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "GAT"},
                {"id": 2, "sequence": "CTAA"},
                {"id": 3, "sequence": "CGAA"},
                {"id": 5, "sequence": "ACA"}
            ],
            "edge": [
                {"from": 1, "to": 2, "to_end": true},
                {"from": 1, "to": 3, "to_end": true},
                {"from": 2, "to": 5, "from_start": true},
                {"from": 3, "to": 5, "from_start": true},
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        algorithms::normalize(&graph);
        
        // Those duplicate Ts (actually As) and Gs (actually Cs) should be eliminated
        REQUIRE(graph.length() == 11);
        
        bool found = false;
        graph.for_each_node([&](const Node* n) {
            if (n->sequence() == "AA") {
                found = true;
            }
        });
        // The Ts ought to have been combined into an AA node.
        REQUIRE(found);
        
    }

}
#endif

}
}
