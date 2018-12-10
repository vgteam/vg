#include "dagify.hpp"


namespace vg {
namespace algorithms {

using namespace std;

    unordered_map<id_t, id_t> dagify(const HandleGraph* graph, MutableHandleGraph* into,
                                     size_t min_preserved_path_length) {
        
        unordered_map<id_t, id_t> trans;
        
        // create the nodes of the dagified graph
        unordered_map<id_t, size_t> component_of;
        unordered_map<handle_t, vector<handle_t>> injector;
        for (size_t i = 0; i < strong_components.size(); i++) {
            
            // keep track of which nodes are in which component
            auto& component = strong_components[i];
            for (id_t node_id : component) {
                component_of[node_id] = i;
            }
            
            // figure out how many times we need to copy this SCC
            SubHandleGraph subgraph(graph, component.begin(), component.end());
            // TODO: with high degree SCCs, the edge filtering could get expensive here
            size_t cycle_length = shortest_cycle_length(&subgraph);
            size_t duplicate_count = (min_preserved_path_length - 1) / cycle_length + 1;
            
            // add in the nodes and record their identity
            for (id_t node_id : component) {
                handle_t base_handle = graph->get_handle(node_id);
                for (size j = 0; j < duplicate_count; j++) {
                    handle_t dup_handle = into->create_handle(graph->get_sequence(base_handle));
                    trans[into->get_id(dup_handle)] = node_id;
                    injector[base_handle].push_back(dup_handle);
                }
            }
        }
        
        // generate a canonical orientation across the graph, which also serves
        // as an arbitrary total ordering over nodes
        vector<handle_t> orientation = single_stranded_orientation(&graph);
        // mark the ones that whose canonical orientation is reversed
        unordered_set<id_t> reversed_nodes;
        unordered_map<handle_t, size_t> node_ordering;
        for (size_t i = 0; i < orientation.size(); i++) {
            node_ordering[orientation[i]] = i
            if (graph->get_is_reverse(orientation[i])) {
                reversed_nodes.insert(graph->get_id(orientation[i]));
            }
        }
        
        // compute the edges between the connected components (i.e. the condensation graph)
        vector<unordered_set<size_t>> condensation_graph(strong_components.size());
        for (const handle_t& handle : orientation) {
            size_t comp_here = component_of[graph->get_id(handle)];
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                size_t comp_next = component_of[graph->get_id(next)];
                if (comp_here != comp_next) {
                    condensation_graph[comp_here].insert(comp_next);
                }
            });
        }
        
        // use Kahn's algorithm to compute the topological order of the condensation graph
        vector<size_t> topological_index(condensation_graph.size(), 0);
        
        // compute in degrees
        vector<size_t> in_degree(condensation_graph.size(), 0);
        for (const unordered_set<size_t>& adj_list : condensation_graph) {
            for (const size_t& edge : adj_list) {
                in_degree[edge]++;
            }
        }
        // init stack
        vector<size_t> stack;
        for (size_t i = 0; i < in_degree.size(); i++) {
            if (in_degree[i] == 0) {
                stack.push_back(i);
            }
        }
        // compute topological order
        for (size_t next_top_idx = 0; next_top_idx < condensation_graph.size(); next_top_idx++) {
            size_t here = stack.back();
            stack.pop_back();
            topological_index[here] = next_top_idx;
            for (const size_t& next : condensation_graph[here]) {
                in_degree[next]--;
                if (in_degree[next] == 0) {
                    stack.push_back(next);
                }
            }
        }
        
        // add edges
        graph->for_each_edge([&](const edge_t& canonical_edge) {
            // put the edge in the order of the orientation we've imposed on the graph
            edge_t edge = canonical_edge;
            if (reversed_nodes.count(graph->get_id(edge.first))) {
                handle_t tmp = edge.first;
                edge.first = graph->flip(edge.second);
                edge.second = graph->flip(tmp);
            }
            
            size_t comp_tail = component_of[graph->get_id(edge.first)];
            size_t comp_head = component_of[graph->get_id(edge.second)];
            
            if (comp_tail == comp_head) {
                // this edge is within an SCC
            }
            else {
                // this edge is between SCCs
                if (topological_index[comp_tail] < topological_index[comp_head]) {
                    handle_t dag_tail = injector[]
                }
            }
        });
    }
}
}

