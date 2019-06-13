/*
Robin Rounthwaite
Find function call in ./subcommand/main.cpp
*/
#include <string>
#include "../vg.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "count_walks.hpp"
#include "../gbwt_helper.hpp"

namespace vg {
    void align_haplotypes(const GBWTGraph& haploGraph, const pair< vector< vector<handle_t> >, vector< vector<handle_t> > >& haplotypes);

    pair< vector< vector<handle_t> >, vector< vector<handle_t> > > extract_haplotypes(const GBWTGraph& graph, const id_t& source_id, const id_t& sink_id);

    vector<vector<handle_t>> find_haplotypes_not_at_source(const GBWTGraph& haploGraph, unordered_set<handle_t>& touched_handles, const id_t& sink_id);

    vector< string > format_handle_haplotypes_to_strings(const GBWTGraph& haploGraph, const vector< vector< handle_t > >& haplotype_handle_vectors);

    vector<string> get_embedded_paths(const PathHandleGraph& graph, const handle_t& source_handle, const handle_t& sink_handle);

    VG align_haplotypes(const vector<string>& walks);

    vector<string> debug_graph_to_strings(MutablePathDeletableHandleGraph& graph, id_t start_id, id_t end_id);

    SubHandleGraph debug_extract_subgraph(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id);
}
