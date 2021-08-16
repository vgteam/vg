#include "0_snarl_analyzer.hpp"
#include "../snarls.hpp"

namespace vg {
namespace algorithms{

void print_handles_in_snarl(const HandleGraph& graph, const id_t& source, const id_t& sink)
{
    //Note: extract_subgraph here assumes that the source is leftmost node of the region.
    //todo: somehow avoid that assumption? Can I look up the snarl in snarl_roots, for example?
    SubHandleGraph snarl = extract_subgraph(graph, source, sink, false);
    // cout << "node ids of snarl with source " << source << " and sink " << sink << endl; 
    snarl.for_each_handle([&](const handle_t handle) 
    {
        cout << graph.get_id(handle) << " ";
    });
}


SnarlAnalyzer::SnarlAnalyzer(const HandleGraph& graph, ifstream &snarl_stream, bool skip_source_sink)
    :_graph(graph) 
{
    // extract the stats on each snarl in the _graph.
    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);

    vector<const Snarl *> snarl_roots = snarl_manager->top_level_snarls();
    for (auto roots : snarl_roots) {
        SubHandleGraph snarl = extract_subgraph(graph, roots->start().node_id(), roots->end().node_id(), roots->start().backward());

        int snarl_size = 0;
        if (skip_source_sink) // skip the source and sink to avoid double-counting of seq.
        {
            snarl.for_each_handle([&](const handle_t handle) {
                if (_graph.get_id(handle) == roots->start().node_id() || _graph.get_id(handle) == roots->end().node_id())
                {
                    // count the number of bases in the snarl, except if it's the source or sink.
                    // this is an imperfect solution to avoid the double-counting problem for 
                    // snarls that share a handle at their borders.
                    snarl_size += snarl.get_sequence(handle).size();
                }
            });

        }
        else // allow potential double-counting of sequence in source and sink.
        {
            snarl.for_each_handle([&](const handle_t handle) {
                // count the number of bases in the snarl.
                snarl_size += snarl.get_sequence(handle).size();
            });

        }

        // save the stats on the snarl for later use.
        _snarl_sources.push_back(roots->start().node_id());
        _snarl_sinks.push_back(roots->end().node_id());
        _snarl_sizes.push_back(snarl_size);
    }

    // TODO: if desired, include this count of all snarls, including non-top-level.
    // int general_count = 0;
    // snarl_manager->for_each_snarl_preorder([&](const vg::Snarl * ignored){
    //     general_count++;
    // });
    // cerr << "number of total snarls in graph: " << general_count << endl;

}

void SnarlAnalyzer::output_snarl_sizes(string& file_name)
{
    std::ofstream outfile;
    outfile.open(file_name);

    for (int i=0; i != _snarl_sources.size(); i++)
    {
        outfile << _snarl_sources[i] << "\t" << _snarl_sinks[i] << "\t" << _snarl_sizes[i] << endl;
    }
}



// TODO: Undo this terrible copy-paste from SnarlNormalizer method, and somehow use 
// TODO:   SnarlNormalizer's code instead. (make this a non-object available function?) 
// Given a start and end node id, construct an extract subgraph between the two nodes
// (inclusive). Arguments:
//      _graph: a pathhandlegraph containing the snarl with embedded paths.
//      source_id: the source of the snarl of interest.
//      sink_id: the sink of the snarl of interest.
// Returns:
//      a SubHandleGraph containing only the handles in _graph that are between start_id
//      and sink_id.
SubHandleGraph extract_subgraph(const HandleGraph &graph,
                                                 id_t source_id,
                                                 id_t sink_id,
                                                 const bool backwards) 
{
    // cerr << "extract_subgraph has source and sink: " << source_id << " " << sink_id << endl; 
    // because algorithm moves left to right, determine leftmost and rightmost nodes.
    id_t leftmost_id;
    id_t rightmost_id;
    // if snarl's "backwards," source is rightmost node, sink is leftmost.
    if (backwards) 
    {
        leftmost_id = sink_id;
        rightmost_id = source_id;
    }
    else 
    {
        leftmost_id = source_id;
        rightmost_id = sink_id;
    }
    // cerr << "extract_subgraph" << endl;
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    unordered_set<id_t> visited;  // to avoid counting the same node twice.
    unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

    // initialize with leftmost_handle (because we move only to the right of leftmost_handle):
    handle_t leftmost_handle = graph.get_handle(leftmost_id);
    subgraph.add_handle(leftmost_handle);
    visited.insert(graph.get_id(leftmost_handle));

    // look only to the right of leftmost_handle
    graph.follow_edges(leftmost_handle, false, [&](const handle_t &handle) {
        // mark the nodes to come as to_visit
        if (visited.find(graph.get_id(handle)) == visited.end()) {
            to_visit.insert(graph.get_id(handle));
        }
    });

    /// explore the rest of the snarl:
    while (to_visit.size() != 0) {
        // remove cur_handle from to_visit
        unordered_set<id_t>::iterator cur_index = to_visit.begin();
        handle_t cur_handle = graph.get_handle(*cur_index);

        to_visit.erase(cur_index);

        /// visit cur_handle
        visited.insert(graph.get_id(cur_handle));

        subgraph.add_handle(cur_handle);

        if (graph.get_id(cur_handle) != rightmost_id) { // don't iterate past rightmost node!
            // look for all nodes connected to cur_handle that need to be added
            // looking to the left,
            graph.follow_edges(cur_handle, true, [&](const handle_t &handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
            // looking to the right,
            graph.follow_edges(cur_handle, false, [&](const handle_t &handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
        }
    }
    return subgraph;
}

}//algorithms
}//vg