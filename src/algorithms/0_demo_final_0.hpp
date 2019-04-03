/*
Robin Rounthwaite
Find function call in ./subcommand/main.cpp
*/
#include <string>
#include "../vg.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "count_walks.hpp"

namespace vg {
    void test_gbwt(MutablePathDeletableHandleGraph& graph);

    void clean_all_snarls(MutablePathDeletableHandleGraph& graph, ifstream& snarl_stream);

    void clean_snarl(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id);

    SubHandleGraph extract_subgraph(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id);

    vector<string> graph_to_strings(MutablePathDeletableHandleGraph& graph, id_t start_id, id_t end_id);
    
    VG strings_to_graph(const vector<string>& walks);

    void integrate_snarl(MutablePathDeletableHandleGraph& graph, HandleGraph& new_snarl, const id_t& start_id, const id_t& end_id);

    pair<vector<handle_t>, vector<handle_t>> get_sources_and_sinks(HandleGraph& graph);
}
