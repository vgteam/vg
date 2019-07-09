#ifndef VG_FLOW_SORT_HPP_INCLUDED
#define VG_FLOW_SORT_HPP_INCLUDED

#include <vg/vg.pb.h>

namespace vg {
    
typedef std::map<id_t, std::vector<Edge*>> EdgeMapping;
static const map<char, char> COMPLEMENTARY_NUCLEOTIDES = {
        {'A', 'T'},
        {'T', 'A'},
        {'C', 'G'},
        {'G', 'C'},
        {'t', 'a'},
        {'a', 't'},
        {'c', 'g'},
        {'g', 'c'}
    };


class FlowSort {
public:
    /*
     * Value for nodes on ref path
     */
    static const size_t DEFAULT_PATH_WEIGHT = 5;
     
    FlowSort(VG& vg);
    /*
     * sorts input graph using max-flow algorithm. Returns sorted list
     */
    std::unique_ptr< list<NodeTraversal> > max_flow_sort(const string& ref_name, bool isGrooming = true);
    /*
     * Fast linear sort
     */
    void fast_linear_sort(const string& ref_name, bool isGrooming = true);
    

    //Structure for holding weighted edges of the graph
    struct WeightedGraph {
        EdgeMapping edges_out_nodes;
        EdgeMapping edges_in_nodes;
        map<Edge*, int> edge_weight;
        
        void construct(FlowSort& fs, const string& ref_name, bool isGrooming = true);
    };
    //Structure contains current path and nodes to growth
    struct Growth {
        set<id_t> nodes;
        set<id_t> backbone;
        list<id_t> ref_path;
        Growth(){}
    };
    
    void flow_sort_nodes(list<NodeTraversal>& sorted_nodes, 
        const string& ref_name, bool isGrooming);
    
    int get_node_degree(WeightedGraph &wg, id_t node_id);

    void update_in_out_edges(EdgeMapping& edges_in, EdgeMapping& edges_out, Edge* e);
    void erase_in_out_edges(EdgeMapping& edges_in, EdgeMapping& edges_out, Edge* e);
    void reverse_edge(Edge* &e);
    void reverse_from_start_to_end_edge(Edge* &e);
    // a(from_start ==true) -> b        =>        not a (from_start == false)  -> b
    id_t from_simple_reverse(Edge* &e);
    // b(from_start ==true) -> a        =>        not a (from_start == false)  -> b
    id_t from_simple_reverse_orientation(Edge* &e);
    // a -> b (to_end ==true)       =>        a -> not b(to_end ==false)
    id_t to_simple_reverse(Edge* &e);
    // b -> a (to_end ==true)       =>        a -> not b(to_end ==false)
    id_t to_simple_reverse_orientation(Edge* &e);
    vector<set<id_t>> get_cc_in_wg(EdgeMapping& edges_in, EdgeMapping& edges_out,
                                   const set<id_t>& all_nodes, id_t start_ref_node);
    void groom_components(EdgeMapping& edges_in, EdgeMapping& edges_out, set<id_t>& isolated_nodes, set<id_t>& main_nodes,
                          map<id_t, set<Edge*>> &minus_start, map<id_t, set<Edge*>> &minus_end);
      
    /* Iterate all edges adjacent to node, recalc degrees of related nodes.
     * If node has no incoming edges we add it to the sources and return as next node.
     * */
    id_t get_next_node_recalc_degrees(WeightedGraph& wg, std::vector<std::set<id_t>>& degrees,std::set<id_t> &sources,
                                         id_t node);
    id_t find_max_node(std::vector<std::set<id_t>> nodes_degree);

    bool bfs(set<id_t>& nodes, map<id_t, map<id_t, int>>& edge_weight, id_t s, id_t t, map<id_t, id_t>& parent);
    void dfs(set<id_t>& nodes, id_t s, set<id_t>& visited, map<id_t, map<id_t, int>>& edge_weight);
    void find_in_out_web(   list<NodeTraversal>& sorted_nodes, 
                            Growth& in_out_growth,
                            WeightedGraph& weighted_graph,
                            set<id_t>& unsorted_nodes, 
                            id_t start_node,
                            bool in_out, int count);
    void process_in_out_growth( EdgeMapping& edges_out_nodes, id_t current_id,
                                Growth& in_out_growth,
                                WeightedGraph& weighted_graph,
                                set<id_t>& visited,
                                list<NodeTraversal>& sorted_nodes, 
                                bool reverse, 
                                set<id_t>& unsorted_nodes,
                                bool in_out, int count);
    void mark_dfs(EdgeMapping& graph_matrix, id_t s, set<id_t>& new_nodes, 
                set<id_t>& visited, bool reverse, set<id_t>& nodes, set<id_t>& backbone);
    vector<pair<id_t,id_t>> min_cut(map<id_t, map<id_t, int>>& graph_weight, set<id_t>& nodes, id_t s, id_t t, 
                EdgeMapping& edges_out_nodes, set<Edge*>& in_joins);
    void remove_edge(EdgeMapping& nodes_to_edges, id_t node, id_t to, bool reverse);
    
private:
    VG& vg;

};
}

#endif /* VG_FLOW_SORT_HPP_INCLUDED*/

