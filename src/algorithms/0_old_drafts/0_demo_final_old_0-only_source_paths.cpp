// /*
// In this phase of the code, I've only extracted paths in the depth_first_search 
// that are overlapping the source node. Next step is to integrate paths that start
// in the middle of the snarl.
// */

// #pragma once //TODO: remove this, to avoid warnings + maybe bad coding practice?
// #include "0_demo_final_0.hpp"
// #include <string>
// #include "../vg.hpp"
// #include "../handle.hpp"
// #include "../subgraph.hpp"
// #include "count_walks.hpp"
// #include <seqan/align.h>
// #include <seqan/graph_msa.h>
// #include <seqan/graph_align.h>
// #include "../msa_converter.hpp"
// #include "../snarls.hpp"
// #include "../gbwt_helper.hpp"
// #include "../stream/vpkg.hpp"
// #include "../../include/handlegraph/path_handle_graph.hpp" //TODO: Do I need this?

// namespace vg {

// // void print_kmer(const std::vector<std::pair<pos_t, size_t>>&, const std::string& string){
// //     cerr << string << endl;
// // }

// // vector<string> get_path_strings(PathHandleGraph& graph, handle_t& source_handle, handle_t& sink_handle) {
// //     unordered_map<string, vector<occurrence_handle_t> > handle_paths get_paths(graph, source_handle, sink_handle);
// //     for (auto path : handle_paths) {
// //         for (occuhandle :  
// //     }
// // }

// // Given a graph and a start_id and end_id representing the beginning and end of the snarl,
// // replaces the nodes between start_id and end_id (inclusive) with the sequence of interest.
// void clean_snarl_from_haplotypes(MutablePathDeletableHandleGraph& graph, const id_t& source_id, const id_t& sink_id){
//     //Convert subgraph of graph, defined by start_id and end_id, into a vector of strings
//     //representing all possible walks through the snarl:
//     vg::handle_t source_handle = graph.get_handle(source_id);
//     vg::handle_t sink_handle = graph.get_handle(sink_id);

//     vector<string> haplotypes = depth_first_haplotypes_to_strings(graph, source_id, sink_id);
//     cerr << "finished depth_first, now on to reference." << endl;
//     vector<string> reference = get_paths(graph, source_handle, sink_handle);

//     haplotypes.insert(end(haplotypes), begin(reference), end(reference));
    
//     //Make a new snarl from walks:
//     VG new_snarl = strings_to_graph(haplotypes);

//     integrate_snarl(graph, new_snarl, source_id, sink_id);
    
// }

// // TODO: test/debug this!
// // Given a snarl in graph defined by source_handle and sink_handle, return all walks associated with an embedded path. 
// // Only walks along embedded paths. Returns a map with string keys and values of vectors of handles, 
// // where each vector of handles represents one path from source to sink.
// // alternative function return:
// //unordered_map<string, vector<occurrence_handle_t> > get_paths(PathHandleGraph& graph, handle_t& source_handle, handle_t& sink_handle){
// vector<string> get_paths(const PathHandleGraph& graph, const handle_t& source_handle, const handle_t& sink_handle){
//     unordered_map<string, vector<occurrence_handle_t> > paths;
//     unordered_map<string, int> multiple_occurrences;

//     // TODO: figure out how to ensure that the occurrence handle is in the correct orientation, i.e. towards the sink.
//     graph.for_each_occurrence_on_handle(source_handle, [&] (const occurrence_handle_t& occurrence) {
//         // Each occurrence represents an embedded path 
//         // (note - in the case of a looped path, there will be multiple occurrences for one path.)
//         // For each path represented by an occurrence, we need to walk along the path until we reach 
//         // the sink node. That series of handles represents the the sequence of the path.
        
//         string path = graph.get_path_name(graph.get_path_handle_of_occurrence(occurrence));
//         if (paths.find(path) != paths.end()){ // if there are multiple occurrences on the same path for source_handle (i.e. a loop)
            
//             //record this in multiple_occurrences, and get the number of times we've seen this occurrence.
//             int occ_num;
//             if (multiple_occurrences.find(path) == multiple_occurrences.end()){
//                 occ_num = 1; // counting from 0, where the first ("zeroeth") occurrence doesn't get a special key name in paths.
//                 multiple_occurrences[path] = occ_num;
//             } else {
//                 occ_num = multiple_occurrences[path]++; // also increments multiple_occurrences.
//             }
            
//             //record the other occurrences with an added identifier to differentiate between paths.
//             paths["occurrence_" + to_string(occ_num) + ":::" + path].emplace_back(occurrence);
//         }
//         else{ // this is the first time we've encountered this occurrence.
//             paths[path].emplace_back(occurrence);
//         }
//     });

//     //Now, for every occurrence, walk along the path until we reach the sink.    
//     for (pair<string, vector<occurrence_handle_t> > path : paths){
//         // cerr << "my name" << path.first << endl;
//         // cerr << "my occurences:" << endl;
//         // for (auto occ : path.second) {
//         //     cerr << "occurrence " << graph.get_sequence(graph.get_occurrence(occ)) << endl;
//         // }
//         // cerr << "testing get_next_occurrence:" << endl;
//         // id_t cur_id = graph.get_id(graph.get_occurrence(path.second));
//         // cerr << cur_id;

//         // cur_occurence is the current handle while walking along the path 
//         occurrence_handle_t cur_occurrence = path.second.back(); 
//         id_t cur_id = graph.get_id(graph.get_occurrence(cur_occurrence));
//         // store the path in paths, in the occurrence_handle_t vector.
//         while (cur_id != graph.get_id(sink_handle)){
//             paths[path.first].push_back(graph.get_next_occurrence(cur_occurrence));
//             // path.second.emplace_back(graph.get_next_occurrence(cur_occurrence));
//             cur_occurrence = paths[path.first].back();
//             // cur_occurrence = path.second.back();
//             cur_id = graph.get_id(graph.get_occurrence(cur_occurrence));
//             cerr << "cur id " << cur_id << " sink id " << graph.get_id(sink_handle) << endl;
//         }
//         path.second.emplace_back(graph.get_next_occurrence(cur_occurrence));
//         cerr << path.second.size() << endl;
//         for (auto handle : path.second) {
//             cerr << graph.get_sequence(graph.get_occurrence(handle));
//         }
//     }
//     cerr << "havin' issues here?" << endl;
//     for (auto path : paths) {
//         for (auto handle : path.second) {
//             cerr << graph.get_sequence(graph.get_occurrence(handle));
//         }
//     }
//     // Resolve multiple_occurrences by identifying which entry in paths 
//     // (of those part of the same path) is longest - that will
//     // represent the full breadth of the path through the snarl.
//     for (pair<string, int> element : multiple_occurrences){
//         // A vector of all the path entries in paths:
//         vector<string> same_path_names = {element.first};
        
//         int max_len = paths[element.first].size();
//         string max_path = element.first;

//         for (int occ_num : range_vector(element.second)){
//             occ_num++; // we actually need range_vector[1, ..., end()]
//             string cur_path = "occurrence_" + to_string(occ_num) + ":::" + element.first; 
//             int cur_len = paths[cur_path].size();
//             same_path_names.push_back(cur_path);

//             if (cur_len > max_len){
//                 max_len = cur_len;
//                 max_path = cur_path;
//             }
//         }

//         // get rid of the smaller fragments of path:
//         for (string name : same_path_names) {
//             if (name != max_path){
//                 paths.erase(name);
//             }
//         }
//     }
//     vector<string> path_strings;
//     // get just the strings from the unordered_map<string, vector<occurrence_handle_t> > paths object:
//     for (auto path : paths) {
//         string path_string;
//         for (auto handle : path.second) {
//             path_string += graph.get_sequence(graph.get_occurrence(handle));
//         } 
//         path_strings.push_back(path_string);
//     }
//     return path_strings;
// }

// //TODO: does GBWTgraphs have names associated with haplotypes? 
// //TODO:     If so, I should change return value to an unordered map with key haplotype name 
// //TODO:     and value vector<handle_t> of all handles in haplotype (also, rename fxn).

// //Depth first search here is based on get_exon_haplotypes from transcriptome.cpp.
// //However, is modified to include all haplotypes inside the source/sink handles, 
// //even ones that don't include the source or sink handles.
// vector<string> depth_first_haplotypes_to_strings(const HandleGraph& graph, const id_t& source_id, const id_t& sink_id){
//     ///GBWT graph construction stuff that belongs in mod_main:
//     ifstream gbwt_stream;
//     string gbwt_name = "test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.gbwt";
//     gbwt_stream.open(gbwt_name);

//     unique_ptr<gbwt::GBWT> gbwt;
//     // Load the GBWT from its container
//     gbwt = stream::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
//     GBWTGraph haploGraph = GBWTGraph(*gbwt, graph);
// // -----------------------------------------------------------------------------------------
//     /// Perform depth first search, where whenever the search reaches sink_handle, convert
//     /// vector of handles to string (should be equivalent to haplotype).
//     //TODO: somehow note/account for how recording of haplotype will be terminated the first time it touches the sink_handle - 
//     //TODO:       this function currently doesn't account for if it loops back.
    
//     //haplotype_queue contains all started exon_haplotypes not completed yet.
//     //Every time we encounter a branch in the paths, the next node down the path
//     //Is stored here, along with the vector of handles that represents the path up
//     //to the SearchState.
//     vector< pair< vector<handle_t>, gbwt::SearchState> > haplotype_queue;

//     // source and sink handle for haploGraph:
//     handle_t source_handle = haploGraph.get_handle(source_id);
//     handle_t sink_handle = haploGraph.get_handle(sink_id);

//     //place source in haplotype_queue.
//     vector<handle_t> source_handle_vec(1, source_handle);
//     gbwt::SearchState source_state = haploGraph.get_state(source_handle);
//     // pair< vector<handle_t>, gbwt::SearchState> source = make_pair(source_handle_vec, source_state);
//     haplotype_queue.push_back( make_pair( source_handle_vec, source_state ) );
// /*
//     cerr << "node id from original graph " << source_id << endl;
//     cerr << "node id from handle (immediately after construction) " << haploGraph.get_id(haploGraph.get_handle(source_id)) << endl;
//     cerr << "node id from handle " << haploGraph.get_id(haplotype_queue.back().first.back()) << endl;
//     cerr << "node id from search state " << haplotype_queue.back().second.node << endl;


//     cerr << "here's the code I want to run" << endl;
//     cerr << "here's the handle I care about " << haploGraph.get_sequence(haplotype_queue.back().first.back()) << endl;
//     haploGraph.node_to_handle(haplotype_queue.back().second.node);
//     cerr<<" made a handle "<< endl;
//     cerr << haploGraph.get_sequence(haploGraph.node_to_handle(haplotype_queue.back().second.node)) << endl;
//     cerr << "finished running code. " << endl;
// */
//     //haplotypes contains all "finished" haplotypes - those that were either walked
//     //to their conclusion, or until they reached the sink.
//     vector< vector<handle_t> >  final_haplotypes;

// /*
//     What if I made a list of start_search_states, which would keep track of places I'd like my backtrace to go
//     (maybe make this a tuple to ensure that end_search_states find the start_search_state I intended for them.)
//     Actually, right now I'm only worrying about paths that start at the source_node. The new paths I encounter
//     partway through the traversal will be dealt with a seperate function call that will make use of bidirected 
//     search to ensure I find all the relevant components of the path.

//     Okay. Keep track of the start_search_state, which is the source_node.
//     Then, have a vector of all bidirected SearchStates (variable expanding_search_states) 
//     that are at "the border" of my current search of the snarl.

//     On each iteration of the while loop, get and pop the cur_search_state 
//     out of the back of expanding_search_states. Using haploGraph.follow_paths, find every search_state that
//     proceeds cur_search_state, and add them to either expanding_search_states or end_search_states.
//         - add cur_search_state to end_search_states if next_search_state.empty() = true; or
//         - add next_search_state to end_search_states if haploGraph.get_handle(next_search_state.node) == sink_handle.
//         - 
        
//         (as a bidirectional) if there are no 
// */



//     while (haplotype_queue.size() != 0) {
//         cerr << "iteration! with haplotype_queue:" << endl;
//         // for (auto hap : haplotype_queue){
//         //     cerr << hap.first.size();
//         // }

//         pair< vector<handle_t>, gbwt::SearchState> cur_haplotype = haplotype_queue.back(); // Tuple of (handles_traversed_so_far, last_touched_SearchState)

//         haplotype_queue.pop_back();

//         vector<gbwt::SearchState> next_searches;

//         haploGraph.follow_paths(cur_haplotype.second, [&](const gbwt::SearchState& next_search) -> bool {
//             next_searches.push_back(next_search);
//             return true;
//         });

//         if (next_searches.size()>1){
//             for (gbwt::SearchState next_search : next_searches){
//                 // copy over the vector<handle_t> of cur_haplotype:
//                 vector<handle_t> next_handle_vec(cur_haplotype.first);
//                 // add the new handle to the vec:
//                 next_handle_vec.push_back(haploGraph.node_to_handle(next_search.node)); //TODO: next_search.node is of type node_type, not node_id. Is that okay?
//                 cerr << haploGraph.get_sequence(haploGraph.node_to_handle(next_search.node)) << endl;
//                 pair< vector<handle_t>, gbwt::SearchState> next_haplotype = make_pair(next_handle_vec, next_search);
//                 haplotype_queue.push_back(next_haplotype);
//             }
//         } else if (haploGraph.get_id(haploGraph.node_to_handle(next_searches.back().node)) == sink_id ) { //TODO: once again, is SearchState.node acceptable here?
//             // Then we need to add cur_haplotype + next_search to final_haplotypes.
//             cur_haplotype.first.push_back(haploGraph.node_to_handle(next_searches.back().node));
//             final_haplotypes.push_back(cur_haplotype.first);

//         } else if ( next_searches.back().empty()) {
//             // Then we have reached the end of the path, and need to add cur_haplotype to final_haplotypes.
//             final_haplotypes.push_back(cur_haplotype.first);
//         } else {
//             // Then there is just one next_search, and it's not the end of the path. 
//             // add (cur_haplotype + next_search to haplotype_queue
//             cur_haplotype.first.push_back(haploGraph.node_to_handle(next_searches.back().node));
//             cur_haplotype.second = next_searches.back();
//             haplotype_queue.push_back(cur_haplotype);
//         }
        
//     }

//     //Now, transform the each vector of handles in final_haplotypes into a string, and return as a vector of strings

//     vector<string> string_haplotypes; 
//     for (vector<handle_t> vector_hap : final_haplotypes){
//         string hap;
//         for (handle_t& handle : vector_hap){
//             hap += haploGraph.get_sequence(handle);
//         }
//         string_haplotypes.push_back(hap);
//     }
    
//     return string_haplotypes;
// }
















// //TODO: delete this function once I've decided I don't want it anymore. Should be replaced with (renamed) depth_first_haplotypes_to_strings.
// // Pull out each haplotype passing through a snarl (defined by source_id and sink_id) as a string.
// vector<string> haplotypes_to_strings(MutablePathDeletableHandleGraph& graph, id_t& source_id, id_t& sink_id){

//     ///stuff that will go in mod_main:
//     ifstream gbwt_stream;
//     string gbwt_name = "test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.gbwt";
//     gbwt_stream.open(gbwt_name);

//     unique_ptr<gbwt::GBWT> gbwt;
//     // Load the GBWT from its container
//     gbwt = stream::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
    
// // -----------------------------------------------------------------
//     /// make subgraph for the snarl:

//     // graph.for_each_handle([&] (const handle_t& handle)-> bool{
//     //     cerr << "test for graph ";
//     //     cerr << graph.get_id(handle) << endl;
//     //     return true;
//     // });

//     SubHandleGraph snarl = extract_subgraph(graph, source_id, sink_id);
    
//     // snarl.for_each_handle_impl([&] (const handle_t& handle)-> bool{
//     //     cerr << "test for snarl ";
//     //     cerr << snarl.get_id(handle) << endl;
//     //     return true;
//     // });
//     // cerr << "before 1 \n";

//     // GBWTGraph haploGraph = GBWTGraph(*gbwt, snarl); //TODO: figure out how to prevent error msg here.
//     GBWTGraph haploGraph = GBWTGraph(*gbwt, graph);
//     // cerr << "after 1 \n";

//     // cerr << "before \n";
//     // haploGraph.for_each_handle([&] (const handle_t& handle)-> bool{
//     //     cerr << "test for haploGraph ";
//     //     cerr << haploGraph.get_id(handle) << endl;
//     //     return true;
//     // });
//     // cerr << "after \n";


//     //TODO:identify source and sinks for troubleshooting!
//     unordered_map<handle_t, vector<string>> sequences; // will contain all haplotype walks through snarl
//     handle_t source_handle = haploGraph.get_handle(source_id);
//     sequences[source_handle].push_back(haploGraph.get_sequence(source_handle));

//     for (const handle_t& handle : algorithms::lazier_topological_order(&haploGraph)) {
        
//         vector<string> seqs_here = sequences[handle];
//         gbwt::SearchState cur_state = haploGraph.get_state(handle);
        
//         // id_t cur_id = haploGraph.get_id(handle);
//         // cerr << "cur_id" << cur_id << endl;

//         haploGraph.follow_paths(cur_state, [&](const gbwt::SearchState& next_search) -> bool {
//             handle_t next_handle = GBWTGraph::node_to_handle(next_search.node);

//             id_t next_id = haploGraph.get_id(next_handle);
//             cerr << "next_id" << next_id << endl;

//             string next_seq = haploGraph.get_sequence(next_handle);
//             // transfer the sequences for the preceding handle to next_handle's sequences,
//             // plus the new handle's sequence.
//             for (string seq : seqs_here){
//                 sequences[next_handle].push_back(seq + next_seq);
//             }
//             return true;
            

//         });
//     }

//     // all the sequences at the sinks will be all the sequences in the snarl.
//     handle_t sink_handle = haploGraph.get_handle(sink_id);
//     return sequences[sink_handle];
//     // vector<string> testVec;
//     // return testVec;
// }

// //Iterate over all snarls in a graph, and run clean_snarl on it.
// void clean_all_snarls(MutablePathDeletableHandleGraph& graph, ifstream& snarl_stream){
//     SnarlManager* snarl_manager = new SnarlManager(snarl_stream);

// /* Use this code to count number of snarls in graph.
// *    int top_count = 0;
// *    for (const Snarl* snarl : snarl_manager->top_level_snarls()){
// *        top_count++;
// *    }
// *    cerr << "number of top_level snarls in graph: " << top_count << endl;
// *
// *    int general_count = 0;
// *    snarl_manager->for_each_snarl_preorder([&](const vg::Snarl * ignored){
// *        general_count++;
// *    });
// *    cerr << "number of total snarls in graph: " << general_count << endl;
// */


//     vector<const Snarl*> snarl_roots = snarl_manager->top_level_snarls();
//     for (auto roots : snarl_roots){
//         clean_snarl(graph, roots->start().node_id(), roots->end().node_id());
//     }
    
//     delete snarl_manager;

    
// }

// // Given a graph and a start_id and end_id representing the beginning and end of the snarl,
// // replaces the nodes between start_id and end_id (inclusive) with the sequence of interest.
// void clean_snarl(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id){
//     //Convert subgraph of graph, defined by start_id and end_id, into a vector of strings
//     //representing all possible walks through the snarl:
//     vector<string> walks = graph_to_strings(graph, start_id, end_id);
    
//     //Make a new snarl from walks:
//     VG new_snarl = strings_to_graph(walks);

//     integrate_snarl(graph, new_snarl, start_id, end_id);
    
// }

// // Given a larger graph and a (usually cleaned snarl) subgraph, integrate new_snarl into the graph at start_id and end_id.
// void integrate_snarl(MutablePathDeletableHandleGraph& graph, HandleGraph& new_snarl, const id_t& start_id, const id_t& end_id){
//     //Get old graph snarl
//     SubHandleGraph graph_snarl = extract_subgraph(graph, start_id, end_id);

//     //Identify old and new snarl start and sink
//     pair<vector<handle_t>, vector<handle_t>> graph_snarl_defining_handles = get_sources_and_sinks(graph_snarl);
//     pair<vector<handle_t>, vector<handle_t>> new_snarl_defining_handles = get_sources_and_sinks(new_snarl);

//     //Check to make sure that newly made snarl has only one start and end.
//     if(new_snarl_defining_handles.first.size() > 1 || new_snarl_defining_handles.second.size() > 1){
//         cerr << "newly made snarl with more than one start or end. # of starts: " << new_snarl_defining_handles.first.size() << " # of ends: " << new_snarl_defining_handles.second.size() << endl;
//         return;
//     }
//     //extract old and new snarl start and sink: 
//     handle_t new_snarl_start = new_snarl_defining_handles.first[0];
//     handle_t new_snarl_end = new_snarl_defining_handles.second[0];

//     handle_t graph_snarl_start = graph_snarl_defining_handles.first[0];
//     handle_t graph_snarl_end = graph_snarl_defining_handles.second[0];

//     ///Replace start and end handles of old graph snarl with new_snarl start and end, and delete
//     ///rest of old graph snarl.
    
//     //Get everything needed to replace graph start and sink.
//     string new_start_seq = new_snarl.get_sequence(new_snarl_start);
//     string new_end_seq = new_snarl.get_sequence(new_snarl_end);
//     id_t new_start_id = graph.get_id(graph_snarl_start);
//     id_t new_end_id = graph.get_id(graph_snarl_end);
//     vector<handle_t> left_of_start;
//     graph.follow_edges(graph_snarl_start, true, [&](const handle_t& handle){
//         left_of_start.emplace_back(handle);
//     });
//     vector<handle_t> right_of_end;
//     graph.follow_edges(graph_snarl_end, false, [&](const handle_t& handle){
//         right_of_end.emplace_back(handle);
//     });

//     //Delete all handles in graph_snarl
//     graph_snarl.for_each_handle([&](const handle_t& handle){
//         graph.destroy_handle(handle);
//     }, false);
    
//     //Make start and end handles for snarl in graph:
//     handle_t new_start_handle = graph.create_handle(new_start_seq, new_start_id);
//     handle_t new_end_handle = graph.create_handle(new_end_seq, new_end_id);

//     //Insert start and end handles:
//     for (handle_t handle : left_of_start) {
//         graph.create_edge(handle, new_start_handle);
//     }
//     for (handle_t handle : right_of_end) {
//         graph.create_edge(new_end_handle, handle);
//     }

//     ///Reintegrate rest of new_snarl. 
//     //topologically ordered new_snarl. As I progress through each node in topo_order,
//     //I can add all the nodes to the right of the snarl. The final node will be the 
//     //end node, which, instead of adding as a new node to graph, I'll re-connect 
//     //to the modified end_node, above.
//     vector<handle_t> new_snarl_topo_order = algorithms::lazier_topological_order(&new_snarl); 

//     //Construct a parallel graph_snarl_topo_order to identify
//     //paralogous nodes between new_snarl and graph.
//     vector<handle_t> graph_snarl_topo_order = {new_start_handle}; 
    
//     for (auto it = ++new_snarl_topo_order.begin(); it != --new_snarl_topo_order.end(); it++){
//         //For every handle in new_snarl, make an (unconnected) handle in graph.
//         string handle_seq = new_snarl.get_sequence(*it);
//         handle_t graph_handle = graph.create_handle(handle_seq);
//         graph_snarl_topo_order.push_back(graph_handle);
//     }

//     graph_snarl_topo_order.push_back(new_end_handle);

//     //Connect the rest of the nodes:
//     for (int i = 0; i < new_snarl_topo_order.size(); i++){
//         // cerr << new_snarl.get_id(new_snarl_topo_order[i]) << endl;

//         new_snarl.follow_edges(new_snarl_topo_order[i], false, [&](const handle_t& snarl_handle){
//             //get topo_index of nodes to be connected to graph start handle
//             auto it = find(new_snarl_topo_order.begin(), new_snarl_topo_order.end(), snarl_handle);
//             int topo_index = it - new_snarl_topo_order.begin();
//             // cerr << "topo_index" << topo_index << endl;
//             // cerr << "i" << i << endl;

//             //connect graph start handle
//             graph.create_edge(graph_snarl_topo_order[i], graph_snarl_topo_order[topo_index]);
//         });
//     }
  
// }

// //Returns tuple of two handles, first being start and second being sink.
// pair<vector<handle_t>, vector<handle_t>> get_sources_and_sinks(HandleGraph& graph){
//     vector<handle_t> sink;
//     vector<handle_t> source;
    
//     // identify sources and sinks
//     graph.for_each_handle([&](const handle_t& handle) {
//         bool is_source = true, is_sink = true;
//         graph.follow_edges(handle, true, [&](const handle_t& prev) {
//             is_source = false;
//             return false;
//         });
//         graph.follow_edges(handle, false, [&](const handle_t& next) {
//             is_sink = false;
//             return false;
//         });
        
//         // base case for dynamic programming
//         if (is_source) {
//             source.push_back(handle);
//         }
//         if (is_sink) {
//             sink.emplace_back(handle);
//         }
//     });

//     return pair<vector<handle_t>, vector<handle_t>>(source, sink);

// }


// VG strings_to_graph(const vector<string>& walks){
//     seqan::Align<seqan::DnaString> align; // create multiple_sequence_alignment object
    
//     seqan::resize(rows(align), walks.size());
//     for (int i = 0; i < walks.size(); ++i){
//         assignSource(row(align, i), walks[i].c_str());
//     }
    

//     globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));

//     stringstream ss;
//     ss << align;
//     MSAConverter myMSAConverter = MSAConverter();
//     myMSAConverter.load_alignments(ss, "seqan");
//     VG snarl = myMSAConverter.make_graph();
//     snarl.clear_paths();


//     // snarl.serialize_to_ostream(cerr);
//     return snarl;
// }




// vector<string> graph_to_strings(MutablePathDeletableHandleGraph& graph, id_t start_id, id_t end_id){
//     SubHandleGraph snarl = extract_subgraph(graph, start_id, end_id);

//     unordered_map<handle_t, vector<string>> sequences;
//     vector<handle_t> sinks;
//     unordered_map<handle_t, size_t> count;
//     count.reserve(snarl.node_size()); // resize count to contain enough buckets for size of snarl
//     sequences.reserve(snarl.node_size()); // resize sequences to contain enough buckets for size of snarl
    
//     // identify sources and sinks //TODO: once we've established that this fxn works, we can just use start_id and end_id. 
//     snarl.for_each_handle([&](const handle_t& handle) {
//         bool is_source = true, is_sink = true;
//         snarl.follow_edges(handle, true, [&](const handle_t& prev) {
//             is_source = false;
//             return false;
//         });
//         snarl.follow_edges(handle, false, [&](const handle_t& next) {
//             is_sink = false;
//             return false;
//         });
        
//         // base case for dynamic programming
//         if (is_source) {
//             count[handle] = 1;
//             sequences[handle].push_back(snarl.get_sequence(handle)); //TODO: presented in the handle's local forward orientation. An issue?
//         }
//         if (is_sink) {
//             sinks.emplace_back(handle);
//         }
//     });

    
//     // count walks by dynamic programming
//     bool overflowed = false;
//     for (const handle_t& handle : algorithms::lazier_topological_order(&snarl)) {
//         size_t count_here = count[handle];
//         vector<string> seqs_here = sequences[handle];

//         snarl.follow_edges(handle, false, [&](const handle_t& next) {
            
//             size_t& count_next = count[next];
//             string seq_next = snarl.get_sequence(next);
            
//             if (numeric_limits<size_t>::max() - count_here < count_next) {
//                 overflowed = true;
//             }

//             else {
//                 count_next += count_here;
//                 // for (auto it = seqs_here.begin(); it == seqs_here.end(); it++){
//                 for (string seq : seqs_here){
//                     sequences[next].push_back(seq + seq_next);
//                 }
//                 // cerr << "next_seqs: ";
//                 // for (string seq : sequences[next]){
//                 //     cerr << seq << endl;
//                 // }
//             }
//         });
//         ///TODO: figure out how to deal with overflow.         
//         // if (overflowed) {
//         //     return numeric_limits<size_t>::max();
//         // }
//     }
    
//     // total up the walks at the sinks
//     size_t total_count = 0;
//     for (handle_t& sink : sinks) {
//         total_count += count[sink];
//     }
    
//     // all the sequences at the sinks will be all the sequences in the snarl.
//     vector<string> walks;
//     for (handle_t& sink : sinks) {
//         for (string seq : sequences[sink]){
//             walks.push_back(seq);
//         }
//     }

//     return walks;
// }


// // given a start and end node id, construct an extract subgraph between the two nodes (inclusive).
// // TODO: change the arguments to handles, which contain orientation within themselves. 
// // That way, iteration to extract the subgraph will have direction contained within themselves.
// // This may actually end up looking like simply parsing an input text file with the handles 
// // described from the find_snarl output.
// SubHandleGraph extract_subgraph(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id){
//     /// make a subgraph containing only nodes of interest. (e.g. a snarl)
//     // make empty subgraph
//     SubHandleGraph subgraph = SubHandleGraph(&graph);

//     unordered_set<id_t> visited; // to avoid counting the same node twice.
//     unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

//     // TODO: how to ensure that "to the right" of start_handle is the correct direction?
//     // initialize with start_handle (because we move only to the right of start_handle):
//     handle_t start_handle = graph.get_handle(start_id);
//     subgraph.add_handle(start_handle);
//     visited.insert(graph.get_id(start_handle));

//     // look only to the right of start_handle
//     graph.follow_edges(start_handle, false, [&](const handle_t& handle){
//         // mark the nodes to come as to_visit
//         if (visited.find(graph.get_id(handle)) == visited.end()) {
//             to_visit.insert(graph.get_id(handle));
//         }
//     });

//     /// explore the rest of the snarl:
//     while (to_visit.size() != 0) {
//         // remove cur_handle from to_visit
//         unordered_set<id_t>::iterator cur_index = to_visit.begin();
//         handle_t cur_handle = graph.get_handle(*cur_index);

//         to_visit.erase(cur_index);

//         /// visit cur_handle
//         visited.insert(graph.get_id(cur_handle));

//         subgraph.add_handle(cur_handle);

//         if (graph.get_id(cur_handle) != end_id){ // don't iterate past end node! 
//             // look for all nodes connected to cur_handle that need to be added
//             // looking to the left,
//             graph.follow_edges(cur_handle, true, [&](const handle_t& handle){
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//             // looking to the right,
//             graph.follow_edges(cur_handle, false, [&](const handle_t& handle){
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//         }
//     }
//     return subgraph;
// }
// }