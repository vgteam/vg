#include "catch.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../flow_sort.hpp"
#include "../algorithms/gfa_to_handle.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("sorts input graph using max flow approach", "[flow_sort]") {
    SECTION("Sort simple graph") {
        const string graph_gfa = R"(H	VN:Z:0.1
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
P	1	ref	1	+	1M
P	2	ref	2	+	1M
P	3	ref	3	+	1M
P	1	path1	1	+	1M
P	4	path1	2	+	1M
P	5	path1	3	+	1M
P	2	path1	4	+	1M
P	1	path2	1	+	1M
P	4	path2	2	+	1M
P	5	path2	3	+	1M
P	6	path2	4	+	1M
P	3	path2	5	+	1M)";
        
        VG vg;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &vg);
        REQUIRE(vg.length() == 6);
        
        FlowSort flow_sort(vg);
        string reference_name = "ref";
        
        flow_sort.max_flow_sort(reference_name);
        
        stringstream res;
        for (int i =0; i < vg.graph.node_size(); ++i ) {
            if (i > 0)
               res << " ";
            res << vg.graph.mutable_node(i)->id();
        }
        REQUIRE(res.str().compare("1 4 5 2 6 3") == 0);
    }
    
    SECTION("Sort simple graph test1") {
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
P	1	ref	1	+	1M
P	2	ref	2	+	1M
P	3	ref	3	+	1M
P	6	path1	1	+	1M
P	5	path1	2	+	1M
P	4	path1	3	+	1M
P	2	path1	4	+	1M
P	6	path2	1	+	1M
P	5	path2	2	+	1M
P	4	path2	3	+	1M
P	2	path2	4	+	1M
P	3	path2	5	+	1M
P	2	path3	1	+	1M
P	11	path3	2	+	1M
P	12	path3	3	+	1M
P	13	path3	4	+	1M
P	1	path4	1	+	1M
P	2	path4	2	+	1M
P	11	path4	3	+	1M
P	12	path4	4	+	1M
P	13	path4	5	+	1M)";
        
        VG vg;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &vg);
        REQUIRE(vg.length() == 17);
        
        FlowSort flow_sort(vg);
        string reference_name = "ref";
        flow_sort.max_flow_sort(reference_name);
        
        stringstream res;
        for (int i =0; i < vg.graph.node_size(); ++i ) {
            if (i > 0)
               res << " ";
            res << vg.graph.mutable_node(i)->id();
        }
        REQUIRE(res.str().compare("1 6 10 5 7 9 4 8 2 16 11 15 17 12 14 13 3") == 0);
    }
    
    SECTION("Sort simple graph test2") {
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
P	1	ref	1	+	1M
P	2	ref	2	+	1M
P	3	ref	3	+	1M
P	1	path1	1	+	1M
P	2	path1	2	+	1M
P	20	path1	3	+	1M
P	21	path1	4	+	1M
P	5	path2	1	+	1M
P	4	path2	2	+	1M
P	13	path2	3	+	1M
P	14	path2	4	+	1M
P	15	path2	5	+	1M
P	20	path3	1	+	1M
P	29	path3	2	+	1M
P	30	path3	3	+	1M
P	31	path3	4	+	1M
P	4	path4	1	+	1M
P	13	path4	2	+	1M
P	14	path4	3	+	1M
P	15	path4	4	+	1M
P	5	path5	1	+	1M
P	4	path5	2	+	1M
P	2	path5	3	+	1M
P	3	path5	4	+	1M
P	6	path6	1	+	1M
P	7	path6	2	+	1M
P	8	path6	3	+	1M
P	4	path6	4	+	1M
P	2	path6	5	+	1M
P	2	path7	1	+	1M
P	20	path7	2	+	1M
P	29	path7	3	+	1M
P	30	path7	4	+	1M
P	31	path7	5	+	1M
P	22	path8	1	+	1M
P	23	path8	2	+	1M
P	24	path8	3	+	1M
P	20	path8	4	+	1M
P	21	path8	5	+	1M
P	6	path9	1	+	1M
P	7	path9	2	+	1M
P	8	path9	3	+	1M
P	4	path9	4	+	1M
P	2	path9	5	+	1M
P	3	path9	6	+	1M
P	22	path10	1	+	1M
P	23	path10	2	+	1M
P	24	path10	3	+	1M
P	20	path10	4	+	1M
P	5	path11	1	+	1M
P	4	path11	2	+	1M
P	2	path11	3	+	1M
P	2	path12	1	+	1M
P	20	path12	2	+	1M
P	21	path12	3	+	1M)";
        
        VG vg;
        stringstream in(graph_gfa);
        algorithms::gfa_to_path_handle_graph(in, &vg); 
        REQUIRE(vg.length() == 35);
        
        FlowSort flow_sort(vg);
        
        string reference_name = "ref";
        flow_sort.max_flow_sort(reference_name);
        
        stringstream res;
        for (int i =0; i < vg.graph.node_size(); ++i ) {
            if (i > 0)
               res << " ";
            res << vg.graph.mutable_node(i)->id();
        }
        REQUIRE(res.str().compare("1 5 6 12 7 9 11 8 10 4 19 13 16 18 14 17 15 2 22 28 23 25 27 24 26 20 35 29 32 34 30 33 31 21 3") == 0);
    }
}
}
}
