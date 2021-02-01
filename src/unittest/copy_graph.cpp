#include "catch.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "xg.hpp"

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"

namespace vg {
    namespace unittest {
        
        using namespace std;
        
        TEST_CASE("copy_handle_graph converter works on empty graph, xg to vg", "[handle][vg][xg]") {
            xg::XG xg;
            VG vg;
           handlealgs::copy_handle_graph(&xg, &vg);
            REQUIRE(vg.node_count() == 0);
            REQUIRE(vg.edge_count() == 0);
        }
        TEST_CASE( "copy_handle_graph converter works on empty graph, xg to pg", "[handle][pg][xg]") {
            xg::XG xg;
            bdsg::PackedGraph pg;
            handlealgs::copy_handle_graph(&xg, &pg);
            REQUIRE(pg.get_node_count() == 0);
            
            int edge_count = 0;
            pg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 0);
        }
        TEST_CASE( "copy_handle_graph converter works on empty graph, xg to hg", "[handle][hg][xg]") {
            xg::XG xg;
            bdsg::HashGraph hg;
            handlealgs::copy_handle_graph(&xg, &hg);
            REQUIRE(hg.get_node_count() == 0);
            
            int edge_count = 0;
            hg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 0);
        }
        
        TEST_CASE( "copy_handle_graph converter works on graphs with one node, xg to vg", "[handle][vg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"}
                         ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            VG vg;
            handlealgs::copy_handle_graph(&xg, &vg);
            
            REQUIRE(xg.get_node_count() == 1);
            REQUIRE(vg.get_node_count() == 1);
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with one node, xg to pg", "[handle][pg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"}
                         ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::PackedGraph pg;
            handlealgs::copy_handle_graph(&xg, &pg);
            
            REQUIRE(xg.get_node_count() == 1);
            REQUIRE(pg.get_node_count() == 1);
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with one node, xg to hg", "[handle][hg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"}
                         ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::HashGraph hg;
            handlealgs::copy_handle_graph(&xg, &hg);
            
            REQUIRE(xg.get_node_count() == 1);
            REQUIRE(hg.get_node_count() == 1);
        }

        TEST_CASE( "copy_handle_graph converter works on graphs with one reversing edge, xg to vg", "[handle][vg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":false, "to_end": false},
                        {"from":4, "to":3, "from_start":false, "to_end": true}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            VG vg;
            handlealgs::copy_handle_graph(&xg, &vg);
            
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(vg.get_node_count() == 4);
            REQUIRE(vg.edge_count() == 4);
            REQUIRE(vg.length() == 16);
            
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with one reversing edge, xg to pg", "[handle][pg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":false, "to_end": false},
                        {"from":4, "to":3, "from_start":false, "to_end": true}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::PackedGraph pg;
            handlealgs::copy_handle_graph(&xg, &pg);
            
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(pg.get_node_count() == 4);

            int length = 0;
            pg.for_each_handle([&](const handle_t& here) {
                length += pg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            pg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 4);
        
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with one reversing edge, xg to hg", "[handle][hg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":false, "to_end": false},
                        {"from":4, "to":3, "from_start":false, "to_end": true}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::HashGraph hg;
            handlealgs::copy_handle_graph(&xg, &hg);
            
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(hg.get_node_count() == 4);
            int length = 0;
            hg.for_each_handle([&](const handle_t& here) {
                length += hg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            hg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 4);
            
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with reversing edges and loops", "[handle][vg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            VG vg;
            handlealgs::copy_handle_graph(&xg, &vg);
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(vg.get_node_count() == 4);
            REQUIRE(vg.edge_count() == 7);
            REQUIRE(vg.length() == 16);
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with reversing edges and loops, xg to pg", "[handle][pg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::PackedGraph pg;
            handlealgs::copy_handle_graph(&xg, &pg);
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(pg.get_node_count() == 4);
            
            int length = 0;
            pg.for_each_handle([&](const handle_t& here) {
                length += pg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            pg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 7);
        }
        TEST_CASE( "copy_handle_graph converter works on graphs with reversing edges and loops, xg to hg", "[handle][hg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::HashGraph hg;
            handlealgs::copy_handle_graph(&xg, &hg);
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(hg.get_node_count() == 4);
            
            int length = 0;
            hg.for_each_handle([&](const handle_t& here) {
                length += hg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            hg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 7);
        }
        
        TEST_CASE( "copy_handle_graph converter works on paths, xg to vg", "[handle][vg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ],
                "path":[
                        {"name":"path1","mapping":[{"position":{"node_id":1},"rank":1},
                                                   {"position":{"node_id":2},"rank":2},
                                                   {"position":{"node_id":3},"rank":3}
                                                   ]
                        },
                        {"name":"path2","mapping":[{"position":{"node_id":3},"rank":11},
                                                   {"position":{"node_id":4},"rank":12},
                                                   {"position":{"node_id":1},"rank":13},
                                                   {"position":{"node_id":1},"rank":14}
                                                   ]
                        }
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            VG vg;
           handlealgs::copy_path_handle_graph(&xg, &vg);
            
            
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(vg.get_node_count() == 4);
            REQUIRE(vg.edge_count() == 7);
            REQUIRE(vg.length() == 16);
            REQUIRE(vg.has_path("path1") == true);
            REQUIRE(vg.has_path("path2") == true);
            REQUIRE(vg.get_path_count() == 2);
            vg.for_each_path_handle([&](const path_handle_t& path) {
                string path_name = vg.get_path_name(path);
                if (path_name == "path1"){
                    REQUIRE(vg.get_step_count(path) == 3);
                }
                if (path_name == "path2"){
                    REQUIRE(vg.get_step_count(path) == 4);
                }
            });
        }
        TEST_CASE( "copy_handle_graph converter works on paths, xg to pg", "[handle][pg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ],
                "path":[
                        {"name":"path1","mapping":[{"position":{"node_id":1},"rank":1},
                                                   {"position":{"node_id":2},"rank":2},
                                                   {"position":{"node_id":3},"rank":3}
                                                   ]
                        },
                        {"name":"path2","mapping":[{"position":{"node_id":3},"rank":11},
                                                   {"position":{"node_id":4},"rank":12},
                                                   {"position":{"node_id":1},"rank":13},
                                                   {"position":{"node_id":1},"rank":14}
                                                   ]
                        }
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::PackedGraph pg;
           handlealgs::copy_path_handle_graph(&xg, &pg);
            
            
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(pg.get_node_count() == 4);

            
            int length = 0;
            pg.for_each_handle([&](const handle_t& here) {
                length += pg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            pg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 7);
            
            
            REQUIRE(pg.has_path("path1") == true);
            REQUIRE(pg.has_path("path2") == true);
            REQUIRE(pg.get_path_count() == 2);
            pg.for_each_path_handle([&](const path_handle_t& path) {
                string path_name = pg.get_path_name(path);
                if (path_name == "path1"){
                    REQUIRE(pg.get_step_count(path) == 3);
                }
                if (path_name == "path2"){
                    REQUIRE(pg.get_step_count(path) == 4);
                }
            });
        }
        TEST_CASE( "copy_handle_graph converter works on paths, xg to hg", "[handle][hg][xg]") {
            string graph_json = R"(
            {
                "node": [
                         {"id":1, "sequence":"GATT"},
                         {"id":2,"sequence":"ACA"},
                         {"id":3,"sequence":"CGAT"},
                         {"id":4,"sequence":"TCGAA"}
                         ],
                "edge":[
                        {"from":1, "to":2, "from_start":false, "to_end": false},
                        {"from":2, "to":3, "from_start":false, "to_end": false},
                        {"from":3, "to":4, "from_start":true, "to_end": true},
                        {"from":4, "to":1, "from_start":true, "to_end": true},
                        {"from":1, "to":1, "from_start":true, "to_end": true},
                        {"from":2, "to":3, "from_start":false, "to_end": true},
                        {"from":4, "to":4, "from_start":true, "to_end": false}
                        ],
                "path":[
                        {"name":"path1","mapping":[{"position":{"node_id":1},"rank":1},
                                                   {"position":{"node_id":2},"rank":2},
                                                   {"position":{"node_id":3},"rank":3}
                                                   ]
                        },
                        {"name":"path2","mapping":[{"position":{"node_id":3},"rank":11},
                                                   {"position":{"node_id":4},"rank":12},
                                                   {"position":{"node_id":1},"rank":13},
                                                   {"position":{"node_id":1},"rank":14}
                                                   ]
                        }
                        ]
            }
            )";
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            xg::XG xg;
            xg.from_path_handle_graph(VG(proto_graph));
            bdsg::HashGraph hg;
           handlealgs::copy_path_handle_graph(&xg, &hg);
            
            
            
            REQUIRE(xg.get_sequence(xg.get_handle(1)) == "GATT");
            REQUIRE(xg.get_sequence(xg.get_handle(3)) == "CGAT");
            REQUIRE(xg.get_node_count() == 4);
            REQUIRE(hg.get_node_count() == 4);
            
            
            int length = 0;
            hg.for_each_handle([&](const handle_t& here) {
                length += hg.get_length(here);
                return true;
            });
            REQUIRE(length == 16);
            
            int edge_count = 0;
            hg.for_each_edge([&](const edge_t& edge) {
                edge_count += 1;
                return true;
            });
            REQUIRE(edge_count == 7);
            
            
            REQUIRE(hg.has_path("path1") == true);
            REQUIRE(hg.has_path("path2") == true);
            REQUIRE(hg.get_path_count() == 2);
            hg.for_each_path_handle([&](const path_handle_t& path) {
                string path_name = hg.get_path_name(path);
                if (path_name == "path1"){
                    REQUIRE(hg.get_step_count(path) == 3);
                }
                if (path_name == "path2"){
                    REQUIRE(hg.get_step_count(path) == 4);
                }
            });
        }
    }
}
