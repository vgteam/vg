#include "catch.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../gfa.hpp"
#include "../algorithms/gfa_to_handle.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("Can import a simple graph from GFA", "[gfa]") {
    const string graph_gfa = R"(H	VN:Z:1.0
S	1	G
L	1	+	2	+	0M
L	1	+	4	+	0M
S	2	T
L	2	+	3	+	0M
S	3	G
S	4	C
L	4	+	5	+	0M
S	5	C
L	5	+	2	+	0M
L	5	+	6	+	0M
S	6	T
L	6	+	3	+	0M
P	ref	1+,2+,3+	1M,1M,1M
P	path1	1+,4+,5+,2+	1M,1M,1M,1M
P	path2	1+,4+,5+,6+,3+	1M,1M,1M,1M,1M)";
        
    VG vg;
    stringstream in(graph_gfa);
    algorithms::gfa_to_path_handle_graph(in, &vg);
    REQUIRE(vg.is_valid());
    REQUIRE(vg.length() == 6);
}

TEST_CASE("Can import a simple graph from GFA with out of order lines and elided fields", "[gfa]") {
    const string graph_gfa = R"(H	VN:Z:1.0
S	1	G
L	5	+	6	+	0M
L	1	+	2	+	0M
P	path1	1+,4+,5+,2+	1M,1M,1M,1M
L	1	+	4	+	0M
L	2	+	3	+	0M
S	3	G
S	4	C
L	4	+	5	+	*
S	5	C
L	5	+	2	+	0M
S	6	T
L	6	+	3	+	0M
P	ref	1+,2+,3+	*
S	2	T
P	path2	1+,4+,5+,6+,3+	1M,1M,1M,1M,1M)";
        
    VG vg;
    stringstream in(graph_gfa);
    algorithms::gfa_to_path_handle_graph(in, &vg);
    REQUIRE(vg.is_valid());
    REQUIRE(vg.length() == 6);
}
   
TEST_CASE("Can import a slightly larger graph from GFA", "[gfa]") {
    const string graph_gfa = R"(H	VN:Z:0.1
S	1	A
L	1	+	2	+	0M
S	2	A
L	2	+	3	+	0M
L	2	+	11	+	0M
S	3	A
S	4	G
L	4	+	2	+	0M
L	4	+	8	+	0M
S	5	T
L	5	+	4	+	0M
L	5	+	7	+	0M
S	6	G
L	6	+	5	+	0M
S	7	A
L	7	+	8	+	0M
S	8	C
S	9	G
L	9	+	4	+	0M
S	10	C
L	10	+	5	+	0M
L	10	+	9	+	0M
S	11	A
L	11	+	12	+	0M
L	11	+	15	+	0M
S	12	G
L	12	+	13	+	0M
L	12	+	14	+	0M
S	13	A
S	14	A
S	15	A
L	15	+	14	+	0M
S	17	C
L	17	+	12	+	0M
S	16	G
L	16	+	11	+	0M
L	16	+	17	+	0M
P	ref	1+,2+,3+	1M,1M,1M
P	path1	6+,5+,4+,2+	1M,1M,1M,1M
P	path2	6+,5+,4+,2+,3+	1M,1M,1M,1M,1M
P	path3	2+,11+,12+,13+	1M,1M,1M,1M)";
        
    VG vg;
    stringstream in(graph_gfa);
    algorithms::gfa_to_path_handle_graph(in, &vg);
    REQUIRE(vg.is_valid());
    REQUIRE(vg.length() == 17);
}

TEST_CASE("Can import an even larger graph from GFA", "[gfa]") {
    const string graph_gfa = R"(H	VN:Z:0.1
S	1	A
L	1	+	2	+	0M
S	2	T
L	2	+	3	+	0M
L	2	+	20	+	0M
S	3	T
S	4	C
L	4	+	2	+	0M
L	4	+	13	+	0M
S	5	A
L	5	+	4	+	0M
S	6	C
L	6	+	7	+	0M
S	7	C
L	7	+	8	+	0M
L	7	+	9	+	0M
S	8	C
L	8	+	4	+	0M
L	8	+	10	+	0M
S	9	G
L	9	+	10	+	0M
S	10	G
S	11	A
L	11	+	8	+	0M
S	12	T
L	12	+	7	+	0M
L	12	+	11	+	0M
S	13	T
L	13	+	14	+	0M
L	13	+	16	+	0M
S	14	T
L	14	+	15	+	0M
L	14	+	17	+	0M
S	15	A
S	17	C
S	16	G
L	16	+	17	+	0M
S	19	T
L	19	+	13	+	0M
L	19	+	18	+	0M
S	18	C
L	18	+	14	+	0M
S	21	A
S	20	T
L	20	+	21	+	0M
L	20	+	29	+	0M
S	23	C
L	23	+	24	+	0M
L	23	+	25	+	0M
S	22	A
L	22	+	23	+	0M
S	25	G
L	25	+	26	+	0M
S	24	A
L	24	+	20	+	0M
L	24	+	26	+	0M
S	27	A
L	27	+	24	+	0M
S	26	A
S	29	T
L	29	+	30	+	0M
L	29	+	32	+	0M
S	28	G
L	28	+	23	+	0M
L	28	+	27	+	0M
S	31	A
S	30	C
L	30	+	31	+	0M
L	30	+	33	+	0M
S	34	G
L	34	+	30	+	0M
S	35	G
L	35	+	29	+	0M
L	35	+	34	+	0M
S	32	T
L	32	+	33	+	0M
S	33	G
P	ref	1+,2+,3+	1M,1M,1M)";
        
    VG vg;
    stringstream in(graph_gfa);
    algorithms::gfa_to_path_handle_graph(in, &vg);
    REQUIRE(vg.is_valid());
    REQUIRE(vg.length() == 35);
}


TEST_CASE("Can import paths in GFA 1.0 but rejects GFA 0.1 multiline paths", "[gfa]") {

    // Both these GFA files say the same thing.

    SECTION("Single-line paths work") {

    const string graph_gfa = R"(H	VN:Z:1.0
S	1	CAAATAAG
S	2	A
S	3	G
P	x	1+,3+	8M,1M
L	1	+	2	+	0M
L	1	+	3	+	0M)";


        VG vg;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &vg);

        REQUIRE(vg.is_valid());

        REQUIRE(vg.node_count() == 3);
        REQUIRE(vg.edge_count() == 2);

        REQUIRE(vg.paths.has_path("x"));
        auto path = vg.paths.path("x");
        REQUIRE(path.mapping_size() == 2);

        REQUIRE(mapping_from_length(path.mapping(0)) == 8);
        REQUIRE(mapping_is_match(path.mapping(0)));
        REQUIRE(path.mapping(0).position().is_reverse() == false);

        REQUIRE(mapping_from_length(path.mapping(1)) == 1);
        REQUIRE(mapping_is_match(path.mapping(1)));
        REQUIRE(path.mapping(1).position().is_reverse() == false);

    }

    SECTION("Multi-line paths are rejected") {

    const string graph_gfa = R"(H	VN:Z:0.1
S	1	CAAATAAG
S	2	A
S	3	G
P	1	x	1	+	8M
P	3	x	2	+	1M
L	1	+	2	+	0M
L	1	+	3	+	0M)";

        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        REQUIRE_THROWS_AS(algorithms::gfa_to_path_handle_graph(in, &graph), algorithms::GFAFormatError);
    }
   
}

TEST_CASE("Can reject graphs that are missing nodes", "[gfa]") {

    const string graph_gfa = R"(H	VN:Z:1.0
S	1	CAAATAAG
S	3	G
L	1	+	2	+	0M
L	1	+	3	+	0M)";

    bdsg::HashGraph graph;
    stringstream in(graph_gfa);
    REQUIRE_THROWS_AS(algorithms::gfa_to_path_handle_graph(in, &graph), algorithms::GFAFormatError);

}

TEST_CASE("Can reject graphs that have duplicate nodes", "[gfa]") {

    const string graph_gfa = R"(H	VN:Z:1.0
S	1	CAAATAAG
S	2	A
S	3	G
L	1	+	2	+	0M
L	1	+	3	+	0M
S	1	GATTACA)";

    bdsg::HashGraph graph;
    stringstream in(graph_gfa);
    REQUIRE_THROWS_AS(algorithms::gfa_to_path_handle_graph(in, &graph), algorithms::GFAFormatError);

}

TEST_CASE("Can reject graphs that have duplicate paths", "[gfa]") {

    const string graph_gfa = R"(H	VN:Z:1.0
P	x	1+,3+	8M,1M
S	1	CAAATAAG
S	2	A
S	3	G
P	x	1+,3+	8M,1M
L	1	+	2	+	0M
L	1	+	3	+	0M)";

    bdsg::HashGraph graph;
    stringstream in(graph_gfa);
    REQUIRE_THROWS_AS(algorithms::gfa_to_path_handle_graph(in, &graph), algorithms::GFAFormatError);

}

TEST_CASE("Can tolerate interesting GFA node IDs", "[gfa]") {

    SECTION("An inoffensive graph is accepted") {

        const string graph_gfa = R"(H	VN:Z:1.0
S	1	GATT
S	2	ACA
L	1	+	2	+	0M)";
        
        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &graph);
        
        REQUIRE(graph.get_node_count() == 2);
        // We assume IDs are assigned in file order
        handle_t h1 = graph.get_handle(1);
        handle_t h2 = graph.get_handle(2);
        REQUIRE(graph.get_sequence(h1) == "GATT");
        REQUIRE(graph.get_sequence(h2) == "ACA");
        REQUIRE(graph.has_edge(h1, h2));
    }

    SECTION("A graph that merges the ends of two nodes is rejected with GFAFormatError") {

        const string graph_gfa = R"(H	VN:Z:1.0
S	1	GATTAC
S	2	ATTACA
L	1	+	2	+	5M)";
        
        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        REQUIRE_THROWS_AS(algorithms::gfa_to_path_handle_graph(in, &graph), algorithms::GFAFormatError);
    }
    
    SECTION("A graph that uses a non-numerical identifier is OK") {

        const string graph_gfa = R"(H	VN:Z:1.0
S	1	GATT
S	Chana	ACA
L	1	+	Chana	+	0M)";
        
        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &graph);
        REQUIRE(graph.get_node_count() == 2);
        // We assume IDs are assigned in file order
        handle_t h1 = graph.get_handle(1);
        handle_t h2 = graph.get_handle(2);
        REQUIRE(graph.get_sequence(h1) == "GATT");
        REQUIRE(graph.get_sequence(h2) == "ACA");
        REQUIRE(graph.has_edge(h1, h2));
    }
    
    SECTION("A graph that uses a negative identifier is OK") {

        const string graph_gfa = R"(H	VN:Z:1.0
S	1	GATT
S	-2	ACA
L	1	+	-2	+	0M)";
        
        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &graph);
        REQUIRE(graph.get_node_count() == 2);
        // We assume IDs are assigned in file order
        handle_t h1 = graph.get_handle(1);
        handle_t h2 = graph.get_handle(2);
        REQUIRE(graph.get_sequence(h1) == "GATT");
        REQUIRE(graph.get_sequence(h2) == "ACA");
        REQUIRE(graph.has_edge(h1, h2));
    }
    
    SECTION("A graph that uses a zero identifier is OK") {

        const string graph_gfa = R"(H	VN:Z:1.0
S	1	GATT
S	0	ACA
L	1	+	0	+	0M)";
        
        bdsg::HashGraph graph;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &graph);
        REQUIRE(graph.get_node_count() == 2);
        // We assume IDs are assigned in file order
        handle_t h1 = graph.get_handle(1);
        handle_t h2 = graph.get_handle(2);
        REQUIRE(graph.get_sequence(h1) == "GATT");
        REQUIRE(graph.get_sequence(h2) == "ACA");
        REQUIRE(graph.has_edge(h1, h2));
    }

}

TEST_CASE("Can import W lines from a GFA", "[gfa]") {

    const string graph_gfa = R"(H	VN:Z:1.0
W	GRCh38	0	chr2	*	*	>3<2<1
S	1	CAAATAAG
S	2	ATTACA
S	3	G
L	1	+	2	+	0M
L	1	+	3	-	0M
W	GRCh38	0	chr1	*	*	>1>2<3
W	GRCh38	0	chr3	0	6	>1>2<3
W	GRCh38	0	chr4	5	6	>3
W	GRCh38	0	chr5	99	*	>1>2<3
W	NA19239	1	chr1	0	8	>1
W	NA19239	1	chr1	14	15	<3
W	NA19239	2	chr1	*	*	>1>2<3
)";

    bdsg::HashGraph graph;
    stringstream in(graph_gfa);
    algorithms::gfa_to_path_handle_graph(in, &graph);
    
    REQUIRE(graph.has_path("GRCh38#chr1"));
	REQUIRE(graph.has_path("GRCh38#chr2"));
	REQUIRE(graph.has_path("GRCh38#chr3[0-6]"));
	REQUIRE(graph.has_path("GRCh38#chr4[5-6]"));
	REQUIRE(graph.has_path("GRCh38#chr5[99]"));
	REQUIRE(graph.has_path("NA19239#1#chr1#0"));
	REQUIRE(graph.has_path("NA19239#1#chr1#14"));
	REQUIRE(graph.has_path("NA19239#2#chr1"));

}


        
}
}
