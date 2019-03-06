#pragma once //TODO: remove this, to avoid warnings + maybe bad coding practice?
#include "0_demo_final_0.hpp"
#include <string>
#include "../vg.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "count_walks.hpp"
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/graph_align.h>
#include "../msa_converter.hpp"
#include "../snarls.hpp"

namespace vg {

void clean_all_snarls(MutablePathDeletableHandleGraph& graph, ifstream& snarl_stream){
    SnarlManager* snarl_manager = new SnarlManager(snarl_stream);
    vector<const Snarl*> snarl_roots = snarl_manager->top_level_snarls();
    for (auto roots : snarl_roots){
        clean_snarl(graph, roots->start().node_id(), roots->end().node_id());
    }
    delete snarl_manager;

    
}

// Given a graph and a start_id and end_id representing the beginning and end of the snarl,
// replaces the nodes between start_id and end_id (inclusive) with the sequence of interest.
void clean_snarl(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id){
    //Convert subgraph of graph, defined by start_id and end_id, into a vector of strings
    //representing all possible walks through the snarl:
    vector<string> walks = graph_to_strings(graph, start_id, end_id);

    //Make a new snarl from walks:
    VG new_snarl = strings_to_graph(walks);

    integrate_snarl(graph, new_snarl, start_id, end_id);
    
}

void integrate_snarl(MutablePathDeletableHandleGraph& graph, HandleGraph& new_snarl, const id_t& start_id, const id_t& end_id){
    //Get old graph snarl
    SubHandleGraph graph_snarl = extract_subgraph(graph, start_id, end_id);

    //Identify old and new snarl start and sink
    pair<vector<handle_t>, vector<handle_t>> graph_snarl_defining_handles = get_sources_and_sinks(graph_snarl);
    pair<vector<handle_t>, vector<handle_t>> new_snarl_defining_handles = get_sources_and_sinks(new_snarl);

    //Check to make sure that newly made snarl has only one start and end.
    if(new_snarl_defining_handles.first.size() > 1 || new_snarl_defining_handles.second.size() > 1){
        cerr << "newly made snarl with more than one start or end. # of starts: " << new_snarl_defining_handles.first.size() << " # of ends: " << new_snarl_defining_handles.second.size() << endl;
        return;
    }
    //extract old and new snarl start and sink: 
    handle_t new_snarl_start = new_snarl_defining_handles.first[0];
    handle_t new_snarl_end = new_snarl_defining_handles.second[0];

    handle_t graph_snarl_start = graph_snarl_defining_handles.first[0];
    handle_t graph_snarl_end = graph_snarl_defining_handles.second[0];

    ///Replace start and end handles of old graph snarl with new_snarl start and end, and delete
    ///rest of old graph snarl.
    
    //Get everything needed to replace graph start and sink.
    string new_start_seq = new_snarl.get_sequence(new_snarl_start);
    string new_end_seq = new_snarl.get_sequence(new_snarl_end);
    id_t new_start_id = graph.get_id(graph_snarl_start);
    id_t new_end_id = graph.get_id(graph_snarl_end);
    vector<handle_t> left_of_start;
    graph.follow_edges(graph_snarl_start, true, [&](const handle_t& handle){
        left_of_start.emplace_back(handle);
    });
    vector<handle_t> right_of_end;
    graph.follow_edges(graph_snarl_end, false, [&](const handle_t& handle){
        right_of_end.emplace_back(handle);
    });

    //Delete all handles in graph_snarl
    graph_snarl.for_each_handle([&](const handle_t& handle){
        graph.destroy_handle(handle);
    }, false);
    
    //Make start and end handles for snarl in graph:
    handle_t new_start_handle = graph.create_handle(new_start_seq, new_start_id);
    handle_t new_end_handle = graph.create_handle(new_end_seq, new_end_id);

    //Insert start and end handles:
    for (handle_t handle : left_of_start) {
        graph.create_edge(handle, new_start_handle);
    }
    for (handle_t handle : right_of_end) {
        graph.create_edge(new_end_handle, handle);
    }

    ///Reintegrate rest of new_snarl. 
    //topologically ordered new_snarl. As I progress through each node in topo_order,
    //I can add all the nodes to the right of the snarl. The final node will be the 
    //end node, which, instead of adding as a new node to graph, I'll re-connect 
    //to the modified end_node, above.
    vector<handle_t> new_snarl_topo_order = algorithms::lazier_topological_order(&new_snarl); 

    //Construct a parallel graph_snarl_topo_order to identify
    //paralogous nodes between new_snarl and graph.
    vector<handle_t> graph_snarl_topo_order = {new_start_handle}; 
    
    for (auto it = ++new_snarl_topo_order.begin(); it != --new_snarl_topo_order.end(); it++){
        //For every handle in new_snarl, make an (unconnected) handle in graph.
        string handle_seq = new_snarl.get_sequence(*it);
        handle_t graph_handle = graph.create_handle(handle_seq);
        graph_snarl_topo_order.push_back(graph_handle);
    }

    graph_snarl_topo_order.push_back(new_end_handle);

    //Connect the rest of the nodes:
    for (int i = 0; i < new_snarl_topo_order.size(); i++){
        // cerr << new_snarl.get_id(new_snarl_topo_order[i]) << endl;

        new_snarl.follow_edges(new_snarl_topo_order[i], false, [&](const handle_t& snarl_handle){
            //get topo_index of nodes to be connected to graph start handle
            auto it = find(new_snarl_topo_order.begin(), new_snarl_topo_order.end(), snarl_handle);
            int topo_index = it - new_snarl_topo_order.begin();
            // cerr << "topo_index" << topo_index << endl;
            // cerr << "i" << i << endl;

            //connect graph start handle
            graph.create_edge(graph_snarl_topo_order[i], graph_snarl_topo_order[topo_index]);
        });
    }
  
}

//Returns tuple of two handles, first being start and second being sink.
pair<vector<handle_t>, vector<handle_t>> get_sources_and_sinks(HandleGraph& graph){
    vector<handle_t> sink;
    vector<handle_t> source;
    
    // identify sources and sinks
    graph.for_each_handle([&](const handle_t& handle) {
        bool is_source = true, is_sink = true;
        graph.follow_edges(handle, true, [&](const handle_t& prev) {
            is_source = false;
            return false;
        });
        graph.follow_edges(handle, false, [&](const handle_t& next) {
            is_sink = false;
            return false;
        });
        
        // base case for dynamic programming
        if (is_source) {
            source.push_back(handle);
        }
        if (is_sink) {
            sink.emplace_back(handle);
        }
    });

    return pair<vector<handle_t>, vector<handle_t>>(source, sink);

}


VG strings_to_graph(const vector<string>& walks){
    seqan::Align<seqan::DnaString> align; // create multiple_sequence_alignment object
    
    seqan::resize(rows(align), walks.size());
    for (int i = 0; i < walks.size(); ++i){
        assignSource(row(align, i), walks[i].c_str());
    }
    

    globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));
    // std::cout << align << "\n";

    stringstream ss;
    ss << align;
    MSAConverter myMSAConverter = MSAConverter();
    myMSAConverter.load_alignments(ss, "seqan");
    VG snarl = myMSAConverter.make_graph();
    snarl.clear_paths();


    // snarl.serialize_to_ostream(cout);
    return snarl;
}




vector<string> graph_to_strings(MutablePathDeletableHandleGraph& graph, id_t start_id, id_t end_id){
    // id_t start_id = 220;
    // id_t end_id = 218;
    SubHandleGraph snarl = extract_subgraph(graph, start_id, end_id);

    unordered_map<handle_t, vector<string>> sequences;
    vector<handle_t> sinks;
    unordered_map<handle_t, size_t> count;
    count.reserve(snarl.node_size()); // resize count to contain enough buckets for size of snarl
    sequences.reserve(snarl.node_size()); // resize sequences to contain enough buckets for size of snarl
    
    // identify sources and sinks
    snarl.for_each_handle([&](const handle_t& handle) {
        bool is_source = true, is_sink = true;
        snarl.follow_edges(handle, true, [&](const handle_t& prev) {
            is_source = false;
            return false;
        });
        snarl.follow_edges(handle, false, [&](const handle_t& next) {
            is_sink = false;
            return false;
        });
        
        // base case for dynamic programming
        if (is_source) {
            count[handle] = 1;
            sequences[handle].push_back(snarl.get_sequence(handle)); //TODO: presented in the handle's local forward orientation. An issue?
        }
        if (is_sink) {
            sinks.emplace_back(handle);
        }
    });

    
    // count walks by dynamic programming
    bool overflowed = false;
    for (const handle_t& handle : algorithms::lazier_topological_order(&snarl)) {
        size_t count_here = count[handle];
        vector<string> seqs_here = sequences[handle];

        snarl.follow_edges(handle, false, [&](const handle_t& next) {
            
            size_t& count_next = count[next];
            string seq_next = snarl.get_sequence(next);
            
            if (numeric_limits<size_t>::max() - count_here < count_next) {
                overflowed = true;
            }

            else {
                count_next += count_here;
                // for (auto it = seqs_here.begin(); it == seqs_here.end(); it++){
                for (string seq : seqs_here){
                    sequences[next].push_back(seq + seq_next);
                }
                // cout << "next_seqs: ";
                // for (string seq : sequences[next]){
                //     cout << seq << endl;
                // }
            }
        });
        ///TODO: figure out how to deal with overflow.         
        // if (overflowed) {
        //     return numeric_limits<size_t>::max();
        // }
    }
    
    // total up the walks at the sinks
    size_t total_count = 0;
    for (handle_t& sink : sinks) {
        total_count += count[sink];
    }
    
    // all the sequences at the sinks will be all the sequences in the snarl.
    vector<string> walks;
    for (handle_t& sink : sinks) {
        for (string seq : sequences[sink]){
            walks.push_back(seq);
        }
    }

    return walks;
}


// given a start and end node id, construct an extract subgraph between the two nodes (inclusive).
// TODO: change the arguments to handles, which contain orientation within themselves. 
// That way, iteration to extract the subgraph will have direction contained within themselves.
// This may actually end up looking like simply parsing an input text file with the handles 
// described from the find_snarl output.
SubHandleGraph extract_subgraph(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id){
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    unordered_set<id_t> visited; // to avoid counting the same node twice.
    unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

    // TODO: how to ensure that "to the right" of start_handle is the correct direction?
    // initialize with start_handle (because we move only to the right of start_handle):
    handle_t start_handle = graph.get_handle(start_id);
    subgraph.add_handle(start_handle);
    visited.insert(graph.get_id(start_handle));

    // look only to the right of start_handle
    graph.follow_edges(start_handle, false, [&](const handle_t& handle){
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

        if (graph.get_id(cur_handle) != end_id){ // don't iterate past end node! 
            // look for all nodes connected to cur_handle that need to be added
            // looking to the left,
            graph.follow_edges(cur_handle, true, [&](const handle_t& handle){
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
            // looking to the right,
            graph.follow_edges(cur_handle, false, [&](const handle_t& handle){
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
        }
    }
    return subgraph;
}
}