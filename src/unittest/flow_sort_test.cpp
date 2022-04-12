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
P	ref	1+,2+,3+	*
P	path1	1+,4+,5+,2+	*
P	path2	1+,4+,5+,6+,3+	*)";
        
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
        const string graph_gfa = R"(
H	VN:Z:1.0
S	1	A
S	2	A
S	3	A
S	4	G
S	5	T
S	6	G
S	7	A
S	8	C
S	9	G
S	10	C
S	11	A
S	12	G
S	13	A
S	14	A
S	15	A
S	16	G
S	17	C
P	path1	6+,5+,4+,2+	*
P	path2	6+,5+,4+,2+,3+	*
P	path3	2+,11+,12+,13+	*
P	path4	1+,2+,11+,12+,13+	*
P	ref	1+,2+,3+	*
L	1	+	2	+	*
L	2	+	11	+	*
L	2	+	3	+	*
L	4	+	2	+	*
L	4	+	8	+	*
L	9	+	4	+	*
L	5	+	4	+	*
L	5	+	7	+	*
L	6	+	5	+	*
L	10	+	5	+	*
L	7	+	8	+	*
L	10	+	9	+	*
L	11	+	15	+	*
L	11	+	12	+	*
L	16	+	11	+	*
L	12	+	14	+	*
L	12	+	13	+	*
L	17	+	12	+	*
L	15	+	14	+	*
L	16	+	17	+	*)";
        
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
        const string graph_gfa = R"(
H	VN:Z:1.0
S	1	A
S	2	T
S	3	T
S	4	C
S	5	A
S	6	C
S	7	C
S	8	C
S	9	G
S	10	G
S	11	A
S	12	T
S	13	T
S	14	T
S	15	A
S	16	G
S	17	C
S	18	C
S	19	T
S	20	T
S	21	A
S	22	A
S	23	C
S	24	A
S	25	G
S	26	A
S	27	A
S	28	G
S	29	T
S	30	C
S	31	A
S	32	T
S	33	G
S	34	G
S	35	G
P	path1	1+,2+,20+,21+	*
P	path10	22+,23+,24+,20+	*
P	path11	5+,4+,2+	*
P	path12	2+,20+,21+	*
P	path2	5+,4+,13+,14+,15+	*
P	path3	20+,29+,30+,31+	*
P	path4	4+,13+,14+,15+	*
P	path5	5+,4+,2+,3+	*
P	path6	6+,7+,8+,4+,2+	*
P	path7	2+,20+,29+,30+,31+	*
P	path8	22+,23+,24+,20+,21+	*
P	path9	6+,7+,8+,4+,2+,3+	*
P	ref	1+,2+,3+	*
L	1	+	2	+	*
L	2	+	20	+	*
L	2	+	3	+	*
L	4	+	2	+	*
L	4	+	13	+	*
L	8	+	4	+	*
L	5	+	4	+	*
L	6	+	7	+	*
L	7	+	9	+	*
L	7	+	8	+	*
L	12	+	7	+	*
L	8	+	10	+	*
L	11	+	8	+	*
L	9	+	10	+	*
L	12	+	11	+	*
L	13	+	16	+	*
L	13	+	14	+	*
L	19	+	13	+	*
L	14	+	17	+	*
L	14	+	15	+	*
L	18	+	14	+	*
L	16	+	17	+	*
L	19	+	18	+	*
L	20	+	29	+	*
L	20	+	21	+	*
L	24	+	20	+	*
L	22	+	23	+	*
L	23	+	25	+	*
L	23	+	24	+	*
L	28	+	23	+	*
L	24	+	26	+	*
L	27	+	24	+	*
L	25	+	26	+	*
L	28	+	27	+	*
L	29	+	32	+	*
L	29	+	30	+	*
L	35	+	29	+	*
L	30	+	33	+	*
L	30	+	31	+	*
L	34	+	30	+	*
L	32	+	33	+	*
L	35	+	34	+	*)";
        
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
