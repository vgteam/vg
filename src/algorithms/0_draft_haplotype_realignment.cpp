#pragma once //TODO: remove this, to avoid warnings + maybe bad coding practice?
#include "0_draft_haplotype_realignment.hpp"

#include <string>
#include <algorithm>

#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/graph_align.h>

#include "../vg.hpp"
#include "../gbwt_helper.hpp"
#include "../stream/vpkg.hpp"
#include "../../include/handlegraph/path_handle_graph.hpp" //TODO: Do I need this?
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/graph_align.h>
#include "../msa_converter.hpp"

//TODO: Tomorrow's goal: edit haplotypes using Jordan's technique to re-integrate your snarl.

namespace vg {

/// Given the haplotypes extracted from the graph in extract_haplotypes,
//         creates a new subgraph made from the realignment of the extracted
//         haplotypes.
void align_haplotypes(const GBWTGraph& haploGraph, const pair< vector< vector<handle_t> >, vector< vector<handle_t> > >& haplotypes){
    vector< string >  haplotypes_from_source_to_sink = format_handle_haplotypes_to_strings(haploGraph, haplotypes.first);
    vector< string > other_haplotypes = format_handle_haplotypes_to_strings(haploGraph, haplotypes.second);
    //TODO: Debug: disamiguate beign/ending regions of nodes by adding leading/trailing AAA seq (essential a special character).
    for (string& hap : haplotypes_from_source_to_sink){
        hap = "AAAAAAAA" + hap + "AAAAAAAA";
    }


    VG new_snarl = align_haplotypes(haplotypes_from_source_to_sink);
    //TODO: Debug workaround to avoid hassle of overwriting inputGraph.
    new_snarl.serialize_to_ostream(cout);
    vector<string> walks = debug_graph_to_strings(new_snarl, 2, 12);

    //TODO: Debug print statements
    // cerr << "source_to_sink haplotypes" << endl;
    // for (string hap : haplotypes_from_source_to_sink){
    //     cerr << hap << endl << endl;
    // }
    // cerr << "source_to_sink_walks" << endl;
    // for (string walk : walks){
    //     cerr << walk << endl << endl;
    // }
    // cerr << "are there any walks that aren't haplotypes?" << endl;
    // for (string walk : walks){
    //     if (find(haplotypes_from_source_to_sink.begin(), haplotypes_from_source_to_sink.end(), walk) != haplotypes_from_source_to_sink.end()){
    //         cerr << "good" << endl;
    //     } else {
    //         cerr << "bad walk" << endl;
    //         cerr << walk << endl;
    //     }
    // }
    cerr << "are there any haps that aren't walks?" << endl;
    for (string hap : haplotypes_from_source_to_sink){
        if (find(walks.begin(), walks.end(), hap) != walks.end()){
            cerr << "good" << endl;
        } else {
            cerr << "bad hap" << endl;
            cerr << hap << endl;
        }
    }

    // cerr << "other haplotypes" << endl;
    // for (string hap : other_haplotypes){
    //     cerr << hap << endl << endl;
    // }
    // vector<string> actually_source_to_sink;
    // vector<string> to_print_other_haps;
    // cerr << "other haplotypes sorted" << endl;
    // for (string hap : other_haplotypes){
    //     if (find(haplotypes_from_source_to_sink.begin(), haplotypes_from_source_to_sink.end(), hap) != haplotypes_from_source_to_sink.end()){
    //         actually_source_to_sink.emplace_back(hap);
    //     } else {
    //         to_print_other_haps.emplace_back(hap);
    //     }

    // }
    // sort(actually_source_to_sink.begin(), actually_source_to_sink.end());
    // cerr << "actually source to sink" << actually_source_to_sink.size() << endl;
    // for (string hap : actually_source_to_sink){
    //     cerr << hap << endl << endl;
    // }
    // cerr << endl << endl << "to_print_other_haps" << to_print_other_haps.size() << endl;
    // sort(to_print_other_haps.begin(), to_print_other_haps.end());
    // for (string hap : to_print_other_haps){
    //     cerr << hap << endl << endl;
    // }

    

}


//Returns: a pair containting two sets of paths (each represented by a vector<handle_t>). The first
//          in the pair represents all paths reaching from source to sink in the snarl, and the
//          second representing all other paths in the snarl (e.g. any that don't reach both
//          source and sink in the graph.)
pair< vector< vector<handle_t> >, vector< vector<handle_t> > > extract_haplotypes(const GBWTGraph& haploGraph,
                                                                                    const id_t& source_id, 
                                                                                    const id_t& sink_id){
    cerr << "depth first begins!" << endl;
    //touched_handles contains all handles that have been touched by the depth_first_search,
    //for later use in other_haplotypes_to_strings, which identifies paths that didn't stretch
    //from source to sink in the snarl.
    unordered_set<handle_t> touched_handles;

    //haplotype_queue contains all started exon_haplotypes not completed yet.
    //Every time we encounter a branch in the paths, the next node down the path
    //Is stored here, along with the vector of handles that represents the path up
    //to the SearchState.
    vector< pair< vector<handle_t>, gbwt::SearchState> > haplotype_queue;

    // source and sink handle for haploGraph:
    handle_t source_handle = haploGraph.get_handle(source_id);
    handle_t sink_handle = haploGraph.get_handle(sink_id);

    //place source in haplotype_queue.
    vector<handle_t> source_handle_vec(1, source_handle);
    gbwt::SearchState source_state = haploGraph.get_state(source_handle);
    haplotype_queue.push_back( make_pair( source_handle_vec, source_state ) );
    touched_handles.emplace(source_handle);

    //haplotypes contains all "finished" haplotypes - those that were either walked
    //to their conclusion, or until they reached the sink.
    vector< vector<handle_t> >  haplotypes_from_source_to_sink;
    vector< vector<handle_t> > other_haplotypes;

    // for every partly-extracted thread, extend the thread until it either reaches
    // the sink of the snarl or the end of the thread.
    while (!haplotype_queue.empty()) {

        // get a haplotype out of haplotype_queue to extend - 
        // a tuple of (handles_traversed_so_far, last_touched_SearchState)
        pair< vector<handle_t>, gbwt::SearchState> cur_haplotype = haplotype_queue.back();
        haplotype_queue.pop_back();

        // get all the subsequent search_states that immediately follow the searchstate from cur_haplotype.
        vector<gbwt::SearchState> next_searches;
        haploGraph.follow_paths(cur_haplotype.second, [&](const gbwt::SearchState next_search) -> bool {
            // cerr << "this node immediately follows cur_haplotypes current search_state." << haploGraph.get_sequence(haploGraph.node_to_handle(next_search.node)) << haploGraph.get_id(haploGraph.node_to_handle(next_search.node)) << endl;
            next_searches.push_back(next_search);
            return true;
        });

        // if next_searches > 1, then we need to make multiple new haplotypes to be recorded in haplotype_queue 
        // or one of the finished haplotype_handle_vectors.
        if (next_searches.size() > 1){

            // for every next_search in next_searches, either create a new, extended cur_haplotype to push into haplotype queue,
            // or place in the haplotypes_from_source_to_sink if haplotype extends to sink,
            // or place in the other_haplotypes if haplotype ends before reaching sink.
            for (gbwt::SearchState next_search : next_searches){
                handle_t next_handle = haploGraph.node_to_handle(next_search.node);

                // copy over the vector<handle_t> of cur_haplotype:
                vector<handle_t> next_handle_vec(cur_haplotype.first);

                // add the new handle to the vec:
                next_handle_vec.push_back(next_handle);

                // if new_handle is the sink, put in haplotypes_from_source_to_sink
                if (haploGraph.get_id(next_handle) == sink_id){
                    haplotypes_from_source_to_sink.push_back(next_handle_vec);
                    
                } else { // keep extending the haplotype!
                
                    pair< vector<handle_t>, gbwt::SearchState> next_haplotype = make_pair(next_handle_vec, next_search);
                    haplotype_queue.push_back(next_haplotype);
                    
                }

                //next_handle will be touched.
                touched_handles.emplace(next_handle);
            }

        } // if next_searches is empty, the path has ended but not reached sink.
        else if ( next_searches.empty() ) {
            //TODO: debug
            // cerr << "next_searches is empty" << endl;

            // We have reached the end of the path, but it doesn't reach the sink.
            // we need to add cur_haplotype to other_haplotypes.
            other_haplotypes.push_back(cur_haplotype.first);

        } // if new_handle is the sink, put in haplotypes_from_source_to_sink
        else if (haploGraph.get_id(haploGraph.node_to_handle(next_searches.back().node)) == sink_id ) {
            // TODO: debug:
            // cerr << "next_searches is sink" << endl;

            // Then we need to add cur_haplotype + next_search to haplotypes_from_source_to_sink.
            handle_t next_handle = haploGraph.node_to_handle(next_searches.back().node);
            cur_haplotype.first.push_back(next_handle);
            haplotypes_from_source_to_sink.push_back(cur_haplotype.first);

            //touched next_search's handle
            touched_handles.emplace(next_handle);

        } //else, there is just one next_search, and it's not the end of the path.
          //just extend the search by adding (cur_haplotype + next_search to haplotype_queue.
        else {

            // get the next_handle from the one next_search.
            handle_t next_handle = haploGraph.node_to_handle(next_searches.back().node);
            // TODO: debug:
            // cerr << "normal extend" << endl;
            // cerr << "this is next_handle" << haploGraph.get_id(next_handle) << endl;


            // modify cur_haplotype with next_handle and next_search. 
            cur_haplotype.first.push_back(next_handle);
            cur_haplotype.second = next_searches.back(); // there's only one next_search in next_searches.
            
            // put cur_haplotype back in haplotype_queue.
            haplotype_queue.push_back(cur_haplotype);
            touched_handles.emplace(next_handle);

        }

    }

    //Find any haplotypes starting from handles not starting at the source, but which
    //still start somewhere inside the snarl.
    vector<vector<handle_t>> haplotypes_not_starting_at_source = find_haplotypes_not_at_source(haploGraph, touched_handles, sink_id);

    // move haplotypes_not_starting_at_source into other_haplotypes:
    other_haplotypes.reserve(other_haplotypes.size() + haplotypes_not_starting_at_source.size());
    move(haplotypes_not_starting_at_source.begin(), haplotypes_not_starting_at_source.end(), back_inserter(other_haplotypes));

    return make_pair(haplotypes_from_source_to_sink, other_haplotypes);
}

vector< string > format_handle_haplotypes_to_strings(const GBWTGraph& haploGraph, const vector< vector< handle_t > >& haplotype_handle_vectors){
    vector< string > haplotype_strings; 
    for (vector<handle_t> haplotype_handles : haplotype_handle_vectors){
        string hap;
        for (handle_t& handle : haplotype_handles){
            hap += haploGraph.get_sequence(handle);
        }
        haplotype_strings.push_back(hap);
    }
    return haplotype_strings;
}

vector<vector<handle_t>> find_haplotypes_not_at_source(const GBWTGraph& haploGraph, unordered_set<handle_t>& touched_handles, const id_t& sink_id){
    //TODO: debug: source handle size?
    // cerr << '\n\n\n\n' << endl;
    // for (id_t node_id = 23493; node_id <= 23505; node_id ++){
    //     handle_t trial_handle = haploGraph.get_handle(node_id);
    //     gbwt::SearchState normal_search = haploGraph.get_state(trial_handle);
    //     cerr << "is normal searchstate at handle " << haploGraph.get_id(trial_handle) << " empty? " << normal_search.empty() << " size: " << normal_search.size() << endl;
    //     gbwt::SearchState new_search = haploGraph.index.prefix(haploGraph.handle_to_node(trial_handle));
    //     cerr << "is the prefix searchstate empty? " << new_search.empty() << " size: " << new_search.size() << endl;
    // }
    
    
    

    
    cerr << "finding haplotypes not at source!" << endl;
    /// Search every handle in touched handles for haplotypes starting at that point.
    // Any new haplotypes will be added to haplotype_queue.     
    vector<pair<vector<handle_t>, gbwt::SearchState>> haplotype_queue;

    // Fully extended haplotypes (or haplotypes extended to the snarl's sink)
    // will be added to finished_haplotypes.
    vector<vector<handle_t>> finished_haplotypes;

    // In addition, we need to put the new handle into to_search, because a path may have
    // started on the new handle (which means we need to start a searchstate there.)
    unordered_set<handle_t> to_search;

    // We don't need to ever check the sink handle, since paths from the sink handle
    // extend beyond snarl.
    handle_t sink_handle = haploGraph.get_handle(sink_id);
    touched_handles.erase(sink_handle);

    // Create nested function for making a new_search:
    auto make_new_search = [&](handle_t handle) { 
        cerr << "lambda" << endl;

        // Are there any new threads starting at this handle?
        gbwt::SearchState new_search = haploGraph.index.prefix(haploGraph.handle_to_node(handle));
        // if (new_search != gbwt::SearchState()){
        if (!new_search.empty()){
            //TODO: Debug code: are searchstates empty?
            cerr << "apparently new thread starts at node: " << haploGraph.get_id(handle) << endl;
            cerr << "is the searchstate empty? " << new_search.empty() << " size: " << new_search.size() << endl;
            // Then add them to haplotype_queue.
            haploGraph.follow_paths(new_search, [&](const gbwt::SearchState& next_search) -> bool {
                
                handle_t next_handle = haploGraph.node_to_handle(next_search.node);

                /// check to make sure that the thread isn't already finished:
                // if next_handle is the sink, or if this thread is only one handle long,
                // then there isn't any useful string to extract from this. 
                if (next_handle != sink_handle || next_search == gbwt::SearchState()){
                    // establish a new thread to walk along.
                    vector<handle_t> new_path;
                    new_path.push_back(handle);
                    new_path.push_back(next_handle);

                    pair<vector<handle_t>, gbwt::SearchState > mypair = make_pair(new_path, next_search);


                    // add the new path to haplotype_queue to be extended.
                    haplotype_queue.push_back(make_pair(new_path, next_search));
                    
                    // if next_handle hasn't been checked for starting threads, add to to_search.
                    if (touched_handles.find(next_handle) == touched_handles.end()){
                        to_search.emplace(next_handle);
                    }
                }
                return true;
            });
        }
    };

    // TODO: Debug code: Search every handle in touched handles for haplotypes starting at that point.
    // for (handle_t handle : touched_handles){
    //     cerr << "isn't a source handle: " << haploGraph.get_sequence(handle) << endl;
    //     make_new_search(handle);
    // }

    /// Extend any paths in haplotype_queue, and add any newly found handles to to_search.
    /// Then, check to see if there are any new threads on handles in to_search.
    /// Extend those threads, and add any newly found handles to to_search, 
    /// then search for threads again in to_search again... repeat until to_search remains
    /// emptied of new handles.

    // for tracking whether the haplotype thread is still extending: 
    bool still_extending;

    // TODO: Debug code: did we find any haplotypes that need extending?
    // cerr << "haps need extending below:" << endl;
    // for (auto handle : to_search){
    //     cerr << "hap needs extending: " << haploGraph.get_id(handle) << " " << haploGraph.get_sequence(handle) << endl;
    // }
    // cerr << "haps queue:" << endl;
    // for (auto hap : haplotype_queue){
    //     handle_t handle = haploGraph.node_to_handle(hap.second.node);
    //     cerr << "need to search hap: " << haploGraph.get_id(handle) << " " << haploGraph.get_sequence(handle) << endl;
    // }
    // extend haplotypes on any nodes found to act as a starting thread.
    while(!to_search.empty() || !haplotype_queue.empty()){
        while (!haplotype_queue.empty()){
            cerr << "extend haplotype_queue" << endl;

            // get a haplotype to extend out of haplotype_queue - a tuple of (handles_traversed_so_far, last_touched_SearchState)
            pair< vector<handle_t>, gbwt::SearchState> cur_haplotype = haplotype_queue.back();
            haplotype_queue.pop_back();

            // get all the subsequent search_states that immediately follow the searchstate from cur_haplotype.
            vector<gbwt::SearchState> next_searches;
            haploGraph.follow_paths(cur_haplotype.second, [&](const gbwt::SearchState& next_search) -> bool {
                next_searches.push_back(next_search);
                return true;
            });

            for (gbwt::SearchState next_search: next_searches){
                handle_t next_handle = haploGraph.node_to_handle(next_search.node);

                // if next_search is empty, then we've fallen off the thread,
                // and cur_haplotype can be placed in finished_haplotypes as is for this thread.
                if (next_search == gbwt::SearchState()){
                    
                    finished_haplotypes.push_back(cur_haplotype.first);

                } 
                // if next_search is on the sink_handle, 
                // then cur_haplotype.first + next_search goes to finished_haplotypes.
                else if (haploGraph.get_id(next_handle) == sink_id){

                    // copy over the vector<handle_t> of cur_haplotype:
                    vector<handle_t> next_handle_vec(cur_haplotype.first);
                    //add next_handle
                    next_handle_vec.push_back(next_handle);
                    //place in finished_haplotypes
                    finished_haplotypes.push_back(next_handle_vec);

                    // also, if next_handle hasn't been checked for new threads, add to to_search.
                    if (touched_handles.find(next_handle) != touched_handles.end()){
                        to_search.emplace(next_handle);
                    }
                    
                } 
                // otherwise, just place an extended cur_haplotype in haplotype_queue.  
                else {

                    // copy over cur_haplotype:
                    pair< vector<handle_t>, gbwt::SearchState> cur_haplotype_copy = cur_haplotype;
                    //modify with next_handle/search
                    cur_haplotype_copy.first.push_back(next_handle);
                    cur_haplotype_copy.second = next_search;
                    // place back in haplotype_queue for further extension.
                    haplotype_queue.push_back(cur_haplotype_copy);

                    // also, if next_handle hasn't been checked for new threads, add to to_search.
                    if (touched_handles.find(next_handle) != touched_handles.end()){
                        to_search.emplace(next_handle);
                    }

                }
            } 



        }
            // Then, make more new_searches from the handles in to_search.
        for (handle_t handle : to_search){
            make_new_search(handle); // will add to haplotype_queue if there's any new_searches to be had.
        }
        to_search.clear();

    }
    return finished_haplotypes;
}


//TODO: make return a vector<vector<handle_t>> instead, then convert using separate fxn.
// Given a snarl in graph defined by source_handle and sink_handle, return all walks associated with an embedded path. 
// Only walks along embedded paths. Returns a map with string keys and values of vectors of handles, 
// where each vector of handles represents one path from source to sink.
// alternative function return:
//unordered_map<string, vector<occurrence_handle_t> > get_paths(PathHandleGraph& graph, handle_t& source_handle, handle_t& sink_handle){
vector<string> get_embedded_paths(const PathHandleGraph& graph, const handle_t& source_handle, const handle_t& sink_handle){
    unordered_map<string, vector<occurrence_handle_t> > paths;
    unordered_map<string, int> multiple_occurrences;

    // TODO: figure out how to ensure that the occurrence handle is in the correct orientation, i.e. towards the sink.
    graph.for_each_occurrence_on_handle(source_handle, [&] (const occurrence_handle_t& occurrence) {
        // Each occurrence represents an embedded path 
        // (note - in the case of a looped path, there will be multiple occurrences for one path.)
        // For each path represented by an occurrence, we need to walk along the path until we reach 
        // the sink node. That series of handles represents the the sequence of the path.
        
        string path = graph.get_path_name(graph.get_path_handle_of_occurrence(occurrence));
        if (paths.find(path) != paths.end()){ // if there are multiple occurrences on the same path for source_handle (i.e. a loop)
            
            //record this in multiple_occurrences, and get the number of times we've seen this occurrence.
            int occ_num;
            if (multiple_occurrences.find(path) == multiple_occurrences.end()){
                occ_num = 1; // counting from 0, where the first ("zeroeth") occurrence doesn't get a special key name in paths.
                multiple_occurrences[path] = occ_num;
            } else {
                occ_num = multiple_occurrences[path]++; // also increments multiple_occurrences.
            }
            
            //record the other occurrences with an added identifier to differentiate between paths.
            paths["occurrence_" + to_string(occ_num) + ":::" + path].emplace_back(occurrence);
        }
        else{ // this is the first time we've encountered this occurrence.
            paths[path].emplace_back(occurrence);
        }
    });

    //Now, for every occurrence, walk along the path until we reach the sink.    
    for (pair<string, vector<occurrence_handle_t> > path : paths){
        // cerr << "my name" << path.first << endl;
        // cerr << "my occurences:" << endl;
        // for (auto occ : path.second) {
        //     cerr << "occurrence " << graph.get_sequence(graph.get_occurrence(occ)) << endl;
        // }
        // cerr << "testing get_next_occurrence:" << endl;
        // id_t cur_id = graph.get_id(graph.get_occurrence(path.second));
        // cerr << cur_id;

        // cur_occurence is the current handle while walking along the path 
        occurrence_handle_t cur_occurrence = path.second.back(); 
        id_t cur_id = graph.get_id(graph.get_occurrence(cur_occurrence));
        // store the path in paths, in the occurrence_handle_t vector.
        while (cur_id != graph.get_id(sink_handle)){
            paths[path.first].push_back(graph.get_next_occurrence(cur_occurrence));
            // path.second.emplace_back(graph.get_next_occurrence(cur_occurrence));
            cur_occurrence = paths[path.first].back();
            // cur_occurrence = path.second.back();
            cur_id = graph.get_id(graph.get_occurrence(cur_occurrence));
            cerr << "cur id " << cur_id << " sink id " << graph.get_id(sink_handle) << endl;
        }
        path.second.emplace_back(graph.get_next_occurrence(cur_occurrence));
        cerr << path.second.size() << endl;
        for (auto handle : path.second) {
            cerr << graph.get_sequence(graph.get_occurrence(handle));
        }
    }
    cerr << "havin' issues here?" << endl;
    for (auto path : paths) {
        for (auto handle : path.second) {
            cerr << graph.get_sequence(graph.get_occurrence(handle));
        }
    }
    // Resolve multiple_occurrences by identifying which entry in paths 
    // (of those part of the same path) is longest - that will
    // represent the full breadth of the path through the snarl.
    for (pair<string, int> element : multiple_occurrences){
        // A vector of all the path entries in paths:
        vector<string> same_path_names = {element.first};
        
        int max_len = paths[element.first].size();
        string max_path = element.first;

        for (int occ_num : range_vector(element.second)){
            occ_num++; // we actually need range_vector[1, ..., end()]
            string cur_path = "occurrence_" + to_string(occ_num) + ":::" + element.first; 
            int cur_len = paths[cur_path].size();
            same_path_names.push_back(cur_path);

            if (cur_len > max_len){
                max_len = cur_len;
                max_path = cur_path;
            }
        }

        // get rid of the smaller fragments of path:
        for (string name : same_path_names) {
            if (name != max_path){
                paths.erase(name);
            }
        }
    }
    vector<string> path_strings;
    // get just the strings from the unordered_map<string, vector<occurrence_handle_t> > paths object:
    for (auto path : paths) {
        string path_string;
        for (auto handle : path.second) {
            path_string += graph.get_sequence(graph.get_occurrence(handle));
        } 
        path_strings.push_back(path_string);
    }
    return path_strings;
}

VG align_haplotypes(const vector<string>& source_to_sink_haplotypes){
    seqan::Align<seqan::DnaString> align; // create multiple_sequence_alignment object
    
    seqan::resize(rows(align), source_to_sink_haplotypes.size());
    for (int i = 0; i < source_to_sink_haplotypes.size(); ++i){
        assignSource(row(align, i), source_to_sink_haplotypes[i].c_str());
    }

    globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));

    stringstream ss;
    ss << align;
    MSAConverter myMSAConverter = MSAConverter();
    myMSAConverter.load_alignments(ss, "seqan");
    VG snarl = myMSAConverter.make_graph();
    snarl.clear_paths();


    // snarl.serialize_to_ostream(cerr);
    return snarl;
}

// ------------------------------ DEBUG CODE BELOW: ------------------------------------------

vector<string> debug_graph_to_strings(MutablePathDeletableHandleGraph& graph, id_t start_id, id_t end_id){
    SubHandleGraph snarl = debug_extract_subgraph(graph, start_id, end_id);

    unordered_map<handle_t, vector<string>> sequences;
    vector<handle_t> sinks;
    unordered_map<handle_t, size_t> count;
    count.reserve(snarl.node_size()); // resize count to contain enough buckets for size of snarl
    sequences.reserve(snarl.node_size()); // resize sequences to contain enough buckets for size of snarl
    
    // identify sources and sinks //TODO: once we've established that this fxn works, we can just use start_id and end_id. 
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
                // cerr << "next_seqs: ";
                // for (string seq : sequences[next]){
                //     cerr << seq << endl;
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
SubHandleGraph debug_extract_subgraph(MutablePathDeletableHandleGraph& graph, const id_t& start_id, const id_t& end_id){
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




/*
Misc. //todo's:
    do I need to fix the fact that find_haplotypes_not_at_source runs forever when given
    a non-snarl? 

    TODO: make it so that gbwt file is customized by user rather than hardcoded.

    TODO: make the demo_0 argument into a better name.

    TODO: make it so that you pass the gbwt file directory to a one-liner function that 
    TODO:       generates gbwt graph, extracts haps, aligns haps, and reintegrates haps.
    TODO:       (eventually will do it for every snarl in the given graph).
*/






/// JUNK:
//TODO: fix the clean_snarl_from_haplotypes fxn to properly combine partial and full alignments. 
//TODO:     make sure that I'm inserting all reference haplotypes in the spot that I wantd
//TODO:    (Now that I've converted depth_first fxn return value to a pair.)
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


//Depth first search here is based on get_exon_haplotypes from transcriptome.cpp.
//However, is modified to include all haplotypes inside the source/sink handles, 
//even ones that don't include the source or sink handles.


