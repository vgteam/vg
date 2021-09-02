#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?

#include "0_oo_normalize_snarls.hpp"
#include "0_snarl_sequence_finder.hpp"
#include <string>

// #include <deps/seqan/include/seqan/align.h>
// #include <deps/seqan/include/seqan/graph_align.h>
// #include <deps/seqan/include/seqan/graph_msa.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"

#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../msa_converter.hpp"
#include "../snarls.hpp"
#include "../vg.hpp"

#include "../types.hpp"
#include "extract_containing_graph.hpp"

#include "multipath_mapper.hpp"

/*
TODO: allow for snarls that have haplotypes that begin or end in the middle of the snarl

TODO: allow normalization of multiple adjacent snarls in one combined realignment.

TODO: test that extract_gbwt haplotypes successfully extracts any haplotypes that start/end in the middle of
TODO:    the snarl.
*/

// todo: add cyclic snarls to the ones to skip, if cyclic snarls turns out to be frequent. Right now, I want to know when I have a cyclic snarl.

int _big_snarl_alignment_job = 200;

namespace vg {
namespace algorithms{
/**
 * To "normalize" a snarl, SnarlNormalizer extracts all the sequences in the snarl as
 * represented in the gbwt, and then realigns them to create a replacement snarl. 
 * This process hopefully results in a snarl with less redundant sequence, and with 
 * duplicate variation combined into a single variant.
*/

SnarlNormalizer::SnarlNormalizer(MutablePathDeletableHandleGraph &graph,
                                 const gbwt::GBWT &gbwt,
                                 const gbwtgraph::GBWTGraph &gbwt_graph,
                                 const int& max_handle_size, 
                                 const int &max_alignment_size /*= MAX_INT*/,
                                 const string &path_finder /*= "GBWT"*/)
    : _graph(graph), _gbwt(gbwt), _max_alignment_size(max_alignment_size),
      _max_handle_size(max_handle_size), _path_finder(path_finder), _gbwt_graph(gbwt_graph){}


/**
 * Iterates over all top-level snarls in _graph, and normalizes them.
 * @param snarl_stream file stream from .snarl.pb output of vg snarls
*/
gbwt::GBWT SnarlNormalizer::normalize_top_level_snarls(ifstream &snarl_stream) {
    // cerr << "disambiguate_top_level_snarls" << endl;
    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);

    int num_snarls_normalized = 0;
    int num_snarls_skipped = 0;
    vector<const Snarl *> snarl_roots = snarl_manager->top_level_snarls();
    
    /**
     * We keep an error record to observe when snarls are skipped because they aren't 
     * normalizable under current restraints. Bools:
     *      0) snarl exceeds max number of threads that can be efficiently aligned,
     *      1) snarl has haplotypes starting/ending in the middle,
     *      2)  some handles in the snarl aren't connected by a thread,
     *      3) snarl is cyclic.
     * There are two additional ints for tracking the snarl size. Ints:
     *      4) number of bases in the snarl before normalization
     *      5) number of bases in the snarl after normalization.
     * Further error records:
     *      6) snarl is trivial (either one or two nodes only), so we skipped normalizing them.
    */ 
    int error_record_size = 7;
    vector<int> one_snarl_error_record(error_record_size, 0);
    vector<int> full_error_record(error_record_size, 0);

    pair<int, int> snarl_sequence_change;

    //todo: debug_code
    int stop_size = 5;
    int num_snarls_touched = 0;

    // int skip_first_few = 1;
    // int skipped = 0;
    
    int snarl_num = 0;
    get_all_gbwt_sequences(1, 15, false);

    for (auto roots : snarl_roots) {
        snarl_num++;

        if (snarl_num==1 || snarl_num==2 || snarl_num==3)
        {
            continue;
            cerr << "not running" << endl;
        }
        
        // if (skipped < skip_first_few){
        //     skipped++;
        //     continue;
        // }
        
        if (num_snarls_touched == stop_size){
            cerr << "breakpoint here" << endl;
            break;
        } else {
            num_snarls_touched++;
        }

        // //todo: debug_print:
        // // if (roots->start().node_id()!=3775521)
        // if (snarl_num!=22530 && snarl_num!=22529)
        // {
        //     continue;
        // }
        // else
        // {
        //     cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        // }
        // cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;

        if (_full_log_print)
        {
            cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        }
        else if (snarl_num%10000 == 0)
        {
            cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        }

        // if (snarl_num%2!=0)
        // {
        if (!roots->start().backward())
        {
            cerr << "not backward" << endl;
            // make_one_edit(1, 15);
            // get_all_gbwt_sequences(roots->start().node_id(), roots->end().node_id(), roots->start().backward());
            make_one_edit(roots->start().node_id(), roots->end().node_id());
        }
        else
        {
            cerr << "backward" << endl;
            make_one_edit(roots->end().node_id(), roots->start().node_id());
        }  

        // }
    }     
        // cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        // cerr << "seq from snarl number " << snarl_num << " with source at: " << _graph.get_sequence(_graph.get_handle(roots->start().node_id())) << " and sink at: " << _graph.get_sequence(_graph.get_handle(roots->end().node_id())) << endl;
        // cerr << "graph range of original graph " << _gbwt_graph.min_node_id() << " " <<  _gbwt_graph.max_node_id() <<endl ;
        // // cerr << "gbwt graph investigation: " << _gbwt_graph.has_node(test) << endl;
        // cerr << _gbwt_graph.get_length(_gbwt_graph.get_handle(roots->start().node_id())) << endl;
        // cerr << "seq from snarl number (using gbwt) " << snarl_num << " with source at: " << _gbwt_graph.get_sequence(_gbwt_graph.get_handle(roots->start().node_id())) << " and sink at: " << _gbwt_graph.get_sequence(_gbwt_graph.get_handle(roots->end().node_id())) << endl;

        // cerr << "backwards value? " << roots->start().backward() << endl;
        // if (roots->start().node_id() == 3881494) {
            // cerr << "root backwards?" << roots->start().backward() << endl;
            // cerr << "disambiguating snarl #"
            //         << (num_snarls_normalized + num_snarls_skipped)
            //         << " source: " << roots->start().node_id()
            //         << " sink: " << roots->end().node_id() << endl;

    //TODO: uncomment the blok of code below, with the "//" part that is indented from the rest.
    //         one_snarl_error_record = normalize_snarl(roots->start().node_id(), roots->end().node_id(), roots->start().backward());
    //         if (!(one_snarl_error_record[0] || one_snarl_error_record[1] ||
    //                 one_snarl_error_record[2] || one_snarl_error_record[3] ||
    //                 one_snarl_error_record[6])) {
    //             // if there are no errors, then we've successfully normalized a snarl.
    //             num_snarls_normalized += 1;
    //             // track the change in size of the snarl.
    //             snarl_sequence_change.first += one_snarl_error_record[4];
    //             snarl_sequence_change.second += one_snarl_error_record[5];
    //             // cerr << "normalized snarl starting at: " << roots->start().node_id() << endl;
    //         } else {
    //             // else, there was an error. Track which errors caused the snarl to not
    //             // normalize.
    //             // note: the ints 4 and 5 are ignored here b/c they're for
    //             // recording the changing size of snarls that are successfully normalized.
    //             for (int i = 0; i < error_record_size; i++) {
    //                 if ( i != 4 && i != 5)
    //                 {
    //                     full_error_record[i] += one_snarl_error_record[i];
    //                 }
    //             }
    //             num_snarls_skipped += 1;
    //         }
            
    //         // //todo: debug_statement for extracting snarl of interest.
    //         // VG outGraph;
    //         // pos_t source_pos = make_pos_t(roots->start().node_id(), false, 0);
    //         // vector<pos_t> pos_vec;
    //         // pos_vec.push_back(source_pos);
    //         // algorithms::extract_containing_graph(&_graph, &outGraph, pos_vec, roots->end().node_id() - roots->start().node_id() + 2);
    //         // outGraph.serialize_to_ostream(cout);
    //         // break;
    //     // }

    // }
    // cerr << endl
    //      << "normalized " << num_snarls_normalized << " snarls, skipped "
    //      << num_snarls_skipped << " snarls because. . .\nthey exceeded the size limit ("
    //      << full_error_record[0] << " snarls),\n"
    //      << "had haplotypes starting/ending in the middle of the snarl ("
    //      << full_error_record[1] << "),\n"
    //     //  << " there were handles not connected by the gbwt info ("
    //     //  << full_error_record[2] << " snarls),\n" 
    //      << "the snarl was cyclic (" << full_error_record[3] << " snarls),\n"
    //      << "or the snarl was trivial - composed of only one or two nodes ("
    //      << full_error_record[6] << " snarls)."
    //      << endl;
    // cerr << "amount of sequence in normalized snarls before normalization: "
    //      << snarl_sequence_change.first << endl;
    // cerr << "amount of sequence in normalized snarls after normalization: "
    //      << snarl_sequence_change.second << endl;

    // cerr << "generating gbwt for normalized graph..." << endl;

    // cerr << "_gbwt_changelog size: " << _gbwt_changelog.size() << endl;

    // for (auto& entry : _gbwt_changelog)
    // {
    //     // if (entry.first.size() - entry.second.size() == 2) //export a good snarl for debugging.
    //     // {
    //     cerr << "old_nodes: " << endl;
    //     for (auto node : entry.first)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    
    //     cerr << "new_nodes: " << endl;
    //     for (auto node : entry.second)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    //     // }
        
    // //     if (false)
    // //     {

    // //         cerr << "hi" << endl;
    // //     }
    // //     for (auto new_node : entry.second)
    // //     {
    // //         cerr << "graph contains the new nodes? " << _graph.contains(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(new_node))) << endl;

    // //     }
    //     // cerr << "size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
    //     // if (_gbwt.find(entry.first.begin(), entry.first.end()).empty())
    //     // {
    //     //     cerr << "found empty path!" << endl;
    //     //     cerr << "size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
            
    //     //     for (auto item : entry.first)
    //     //     {
    //     //         cerr << "seq in node: " << _gbwt_graph.get_sequence(_gbwt_graph.node_to_handle(item)) << endl;
    //     //         entry.second.push_back(item);
    //     //     }
    //     // }
    //     // cerr << _gbwt.find(entry.first.begin(), entry.first.end()) << endl;
    //     // for (auto item : entry.first)
    //     // {
    //     //     cerr << "item" << endl;
    //         // cerr << "seq in node: " << _gbwt_graph.get_sequence(_gbwt_graph.node_to_handle(item)) << endl;
    //         // entry.second.push_back(item);
    //     // }
    // }
    // for (auto& entry : _gbwt_changelog)
    // {
    //     cerr << "2nd size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
    // }
    // cerr << _gbwt.empty() << _gbwt_changelog.empty() << endl;
    // _gbwt_changelog.push_back()
    
    cerr << "contents of gbwt changelog: " << endl;
    for (auto& entry : _gbwt_changelog)
    {
        // if (entry.first.size() - entry.second.size() == 2) //export a good snarl for debugging.
        // {
        cerr << "old_nodes: " << endl;
        for (auto node : entry.first)
        {
            cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
        }
    
        cerr << "new_nodes: " << endl;
        for (auto node : entry.second)
        {
            cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
        }
    }
    
    gbwt::GBWT output_gbwt = rebuild_gbwt(_gbwt, _gbwt_changelog);
    cerr << "finished generating gbwt." << endl;
    
    //todo: debug-code for checking that I can build the gbwt_graph:
    cerr << "making new gbwt graph." << endl;
    gbwtgraph::GBWTGraph output_gbwt_graph = gbwtgraph::GBWTGraph(output_gbwt, _graph);
    cerr << "new gbwt graph created." << endl;

    
    // cerr << "output have second path?" << endl;
    // for (auto& entry : _gbwt_changelog)
    // {
    //     cerr << _gbwt.find(entry.second.begin(), entry.second.end()).empty() << endl;
    // }


    //todo: debug_statement for extracting snarl of interest.
    // VG outGraph;
    // pos_t source_pos = make_pos_t(269695, false, 0);
    // vector<pos_t> pos_vec;
    // pos_vec.push_back(source_pos);
    // algorithms::extract_containing_graph(&_graph, &outGraph, pos_vec, 1000);
    // _graph = outGraph;
    // vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph *>(outGraph.get()), cout);
    // outGraph.serialize_to_ostream(cout);

    delete snarl_manager;
    
    return output_gbwt;


}

void SnarlNormalizer::get_all_gbwt_sequences(id_t source_id, id_t sink_id, bool backwards)
{
    cerr << "in get_all_gbwt_sequences" << endl;
    id_t leftmost_id;
    id_t rightmost_id;
    if (backwards) {
        leftmost_id = sink_id;
        rightmost_id = source_id;
    }
    else {
        leftmost_id = source_id;
        rightmost_id = sink_id;
    }

    
    SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);

    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, source_id, sink_id, backwards);

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();
    // unordered_set<string> hap_strs = format_handle_haplotypes_to_strings(get<0>(gbwt_haplotypes));
    for (auto hap : get<0>(gbwt_haplotypes))
    {
        cerr << "new hap:" << endl;
        for (auto handle : hap)
        {
            cerr << _gbwt_graph.get_id(handle) << endl;
        }
        // cerr << hap_str << endl;
    }
}

void SnarlNormalizer::make_one_edit(id_t leftmost_id, id_t rightmost_id) 
{
    cerr << "this is running. Source: " << leftmost_id << " sink: " << rightmost_id << endl;
    ///first, find a path through the snarl. In gbwt_graph.
    vector<handle_t> path;
    path.push_back(_gbwt_graph.get_handle(leftmost_id));

    gbwt::SearchState first_state = _gbwt_graph.get_state(path.back());

    cerr << "testing out a path example." << endl;
    gbwt::SearchState next_state = first_state;
    cerr << "next state: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) << endl;
    while (_gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) != rightmost_id)
    {
        _gbwt_graph.follow_paths(next_state, [&](const gbwt::SearchState one_next_state) -> bool {
                                     next_state = one_next_state;
                                     return false;
                                 });
        path.push_back(_gbwt_graph.node_to_handle(next_state.node));
        cerr << "next state: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) << endl;
    }
    cerr << "finished a path example." << endl;

    ///store the old path as a gbwt::vector_type. and save the path sequence as a string, for the replacement_handle:
    string replacement_string;
    gbwt::vector_type old_path;
    for (auto old_handle : path)
    {
        if (_gbwt_graph.get_id(old_handle)!= leftmost_id && _gbwt_graph.get_id(old_handle)!= rightmost_id)
        {
           replacement_string.append(_gbwt_graph.get_sequence(old_handle)); 
        }
        old_path.push_back(gbwt::Node::encode(_gbwt_graph.get_id(old_handle), _gbwt_graph.get_is_reverse(old_handle)));
    }


    ///remove that path, except for the source/sink.
    for (handle_t gbwt_handle : path)
    {
        if (_gbwt_graph.get_id(gbwt_handle) != leftmost_id && _gbwt_graph.get_id(gbwt_handle) != rightmost_id )
        {
            handle_t normal_handle = _graph.get_handle(_gbwt_graph.get_id(gbwt_handle), _gbwt_graph.get_is_reverse(gbwt_handle));
            _graph.destroy_handle(normal_handle);
        }
    }

    ///replace it with another
    handle_t replacement_handle = _graph.create_handle(replacement_string);
    _graph.create_edge(_graph.get_handle(leftmost_id), replacement_handle);
    _graph.create_edge(replacement_handle, _graph.get_handle(rightmost_id));
    
    ///log that change in the gbwt_changelog.

    gbwt::vector_type new_path;
    new_path.push_back(gbwt::Node::encode(leftmost_id, _gbwt_graph.get_is_reverse(path.front())));
    new_path.push_back(gbwt::Node::encode(_graph.get_id(replacement_handle), false));
    new_path.push_back(gbwt::Node::encode(rightmost_id, _gbwt_graph.get_is_reverse(path.back())));

    _gbwt_changelog.push_back(make_pair(old_path, new_path));
}


/**
 * Normalize a single snarl defined by a source and sink. Only extracts and realigns 
 * sequences found in the gbwt. 
 * @param source_id the source of the snarl of interest.
 * @param sink_id the sink of the snarl of interest.
 * @param error_record an empty vector of 6 integers.
*/
// Returns: none.
// TODO: allow for snarls that have haplotypes that begin or end in the middle of the
// snarl.
vector<int> SnarlNormalizer::normalize_snarl(id_t source_id, id_t sink_id, const bool backwards) {
    // if (backwards){
    //     // swap the source and sink ids. Essentially, this guarantees I treat the leftmost node in snarl as "source".
    //     // (although some adjustments for paths need be made)
    //     id_t swap_source = sink_id; //temp storage of sink_id value. 
    //     sink_id = source_id;
    //     source_id = swap_source; 
    // }

    id_t leftmost_id;
    id_t rightmost_id;
    if (backwards) {
        leftmost_id = sink_id;
        rightmost_id = source_id;
    }
    else {
        leftmost_id = source_id;
        rightmost_id = sink_id;
    }
    /**
     * We keep an error record to observe when snarls are skipped because they aren't 
     * normalizable under current restraints. Bools:
     *      0) snarl exceeds max number of threads that can be efficiently aligned,
     *      1) snarl has haplotypes starting/ending in the middle,
     *      2)  some handles in the snarl aren't connected by a thread,
     *      3) snarl is cyclic.
     * There are two additional ints for tracking the snarl size. Ints:
     *      4) number of bases in the snarl before normalization
     *      5) number of bases in the snarl after normalization.
     *      6) snarl is trivial (either one or two nodes only), so we skipped normalizing them.
    */ 
    vector<int> error_record(7, 0);
    // //todo: debug_statement: determining whether cyclic problem in yeast graph goes away when I swapo source and sink. 
    // SubHandleGraph snarl = extract_subgraph(_graph, sink_id, source_id);
    SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);

    // //todo: debug_statement: Evaluate connections of all nodes in subgraph.
    // snarl.for_each_handle([&](const handle_t handle){
    //     cerr << "examining left neighbors of handle " << snarl.get_id(handle) << ":" << endl;
    //     snarl.follow_edges(handle, false, [&](const handle_t &next) {
    //         cerr << "     " << snarl.get_id(next) << endl;
    //     });
    // });

    if (!handlealgs::is_acyclic(&snarl)) {
        cerr << "snarl at " << source_id << " is cyclic. Skipping." << endl;
        error_record[3] = true;
        return error_record;
    }

    // only normalize non-trivial snarls (i.e. not composed of just a source and sink.):
    int num_handles_in_snarl = 0;
    snarl.for_each_handle([&](const handle_t handle){
        num_handles_in_snarl++;
        if (num_handles_in_snarl >= 3)
        {
            return;
        }
    });
    if (num_handles_in_snarl <= 2)
    {
        if (_full_log_print)
        {
            cerr << "snarl with source " << source_id << " and sink " << sink_id << " has"
                << " only " << num_handles_in_snarl << " nodes. Skipping normalization of"
                << " trivial snarl." << endl;
        }
        error_record[6] += 1;
        return error_record;
    }


    // extract threads
    // haplotypes is of format:
    // 0: a set of all the haplotypes which stretch from source to sink, in string format.
    //   - it's a set, so doesn't contain duplicates
    // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
    // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
    tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<handle_t>> haplotypes;
    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, source_id, sink_id, backwards);
    
    vector<pair<gbwt::vector_type, string>> source_to_sink_gbwt_paths;
    if (_path_finder == "GBWT") {
        tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();


        // cerr << "naive? gbwt haplotypes extract: " << endl;
        // for (auto hap : get<0>(gbwt_haplotypes)) 
        // {
        //     for (auto handle : hap)
        //     {
        //         cerr << "handle id: " << _graph.get_id(handle) << " seq: " << _graph.get_sequence(handle) << endl;
        //     }
        // }
        // Convert the haplotypes from vector<handle_t> format to string format.
        get<0>(haplotypes) = format_handle_haplotypes_to_strings(get<0>(gbwt_haplotypes));
        //todo: possibly remove the duplicate storage of gbwt info in source_to_sink_gbwt_paths, by finding a way to only pass the gbwt info to the "log_gbwt_changes" function. (currently, get<0>haplotypes will also include any source-to-sink paths embedded in the graph.)
        //deep copy of gbwt_haplotypes.
        for (vector<handle_t> hap_handles : get<0>(gbwt_haplotypes))
        {
            string hap_str;
            gbwt::vector_type hap_ids;
            for (handle_t &handle : hap_handles) 
            {
                // cerr << "id of node: " << _gbwt_graph.get_id(handle) << endl;
                // cerr << "sequence of node: " << _gbwt_graph.get_sequence(handle) << endl;
                hap_ids.emplace_back(_gbwt_graph.handle_to_node(handle));
                // hap_str += _gbwt_graph.get_sequence(handle);
                hap_str += _graph.get_sequence(_graph.get_handle(_gbwt_graph.get_id(handle), _gbwt_graph.get_is_reverse(handle)));
                // cerr << "hap_str: " << hap_str << endl;
            }
            
            pair<gbwt::vector_type, string> hap = make_pair(hap_ids, hap_str);
            source_to_sink_gbwt_paths.emplace_back(hap);
        }
        // for (auto item : source_to_sink_gbwt_paths)
        // {
        //     for (auto nid : item.first)
        //     {
        //         cerr << "is the handle is-reverse? of the handles in 'before': " << _graph.get_is_reverse(_graph.get_handle(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(nid)), _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(nid)))) << endl;

        //     }

        // }
        get<1>(haplotypes) = get<1>(gbwt_haplotypes);
        get<2>(haplotypes) = get<2>(gbwt_haplotypes);
        // cerr << "haplotypes after formatting to strings: " << endl;
        // for (auto hap : get<0>(haplotypes)) 
        // {
        //     cerr << "hap: " << hap << endl;
        // }
        
    } else if (_path_finder == "exhaustive") {
        //todo: to enable support for exhaustive, make tests, run them, and also set up support for when I log changes for the gbwt update.
        cerr << "'exhaustive' path finder currently unsupported. Use 'GBWT'. '" << "'." << endl;
        exit(1);

        // pair<unordered_set<string>, unordered_set<handle_t>> exhaustive_haplotypes =
        //     sequence_finder.find_exhaustive_paths();
        // get<0>(haplotypes) = exhaustive_haplotypes.first;
        // get<2>(haplotypes) = exhaustive_haplotypes.second;
    } else {
        cerr << "path_finder type must be 'GBWT' or 'exhaustive', not '" << _path_finder
             << "'." << endl;
        exit(1);
    }

    // check to make sure that the gbwt _graph has threads connecting all handles:
    // ( needs the unordered_set from extract_gbwt haplotypes to be equal to the number of
    // handles in the snarl).
    unordered_set<handle_t> handles_in_snarl;
    snarl.for_each_handle([&](const handle_t handle) {
        handles_in_snarl.emplace(handle);
        // count the number of bases in the snarl.
        error_record[4] += snarl.get_sequence(handle).size();
    });

    // Print a heads-up about snarls that require an alignment with a greater number of 
    // threads than _big_snarl_alignment_job, so the user knows if they are hung up on a 
    // large job.
    if (get<0>(haplotypes).size() > _big_snarl_alignment_job)
    {
        cerr << "WARNING: aligning a snarl requiring a large number (> " << _big_snarl_alignment_job <<") of threads. Number of threads: " << get<0>(haplotypes).size() << endl;
    }
    // Record start time, for measuring alignment time for a big snarls:
    //todo: also add compute time (vs wall clock) measure.
    auto _big_snarl_time_start = chrono::high_resolution_clock::now();    

    // TODO: this if statement only permits snarls that satsify requirements, i.e.
    // TODO:    there are no haplotype begins/ends in the middle
    // TODO:    of the snarl. Get rid of this once alignment issue is addressed!
    // TODO: also, limits the number of haplotypes to be aligned, since snarl starting at
    // TODO:    2049699 with 258 haplotypes is taking many minutes.
    if (get<1>(haplotypes).empty() && get<0>(haplotypes).size() < _max_alignment_size)
        // the following bool check was to ensure that all the handles in the handlegraph 
        // are touched by the gbwt. Turns out, though, that this isn't necessary. If we 
        // assume that all seq info is in gbwt, the gbwt is all we need to worry about.:
        // && get<2>(haplotypes).size() == handles_in_snarl.size()) {
        {
        // Get the embedded paths in the snarl from _graph, to move them to new_snarl.
        // Any embedded paths not in gbwt are aligned in the new snarl.
        vector<pair<step_handle_t, step_handle_t>> embedded_paths =
            sequence_finder.find_embedded_paths();

        //todo: debug_statement
        // cerr << "strings in path_seq before adding haplotypes: " << endl;
        // for (auto path : get<0>(haplotypes))
        // {
        //     cerr << path << endl;
        // }

        
        // TODO: once haplotypes that begin/end in the middle of the snarl have been
        // TODO:    accounted for in the code, remove next chunk of code that finds 
        // TODO: source-to-sink paths.
        // find the paths that stretch from source to sink:
        // cerr << "~~~~~~~~~~source: " << source_id << "sink: " << sink_id << endl;
        for (auto path : embedded_paths) 
        {

            // cerr << "checking path of name " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << " with source " << _graph.get_id(_graph.get_handle_of_step(path.first)) << " and sink " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << endl;
            // cerr << "SOURCE info: prev step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << "prev prev step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(_graph.get_previous_step(path.second)))) << " source: " << _graph.get_id(_graph.get_handle_of_step(path.second)) << " next step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_next_step(path.second))) << endl;
            // cerr << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << " " << source_id << " source bool: " <<  (_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == source_id) << endl;
            if (_graph.get_id(_graph.get_handle_of_step(path.first)) == source_id &&
                _graph.get_id(_graph.get_handle_of_step(
                    _graph.get_previous_step(path.second))) == sink_id)  {
                // cerr << "path_seq added to haplotypes. " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << endl;

                // cerr << "******************************************\nadding path of name " <<
                // _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) <<
                // endl; 
                // get the sequence of the source to sink path, and add it to the
                // paths to be aligned.
                string path_seq;
                step_handle_t cur_step = path.first;
                while (cur_step != path.second) {
                    // cerr << "while adding path, looking at node " << _graph.get_id(_graph.get_handle_of_step(cur_step)) << " with seq " << _graph.get_sequence(_graph.get_handle_of_step(cur_step)) << endl;
                    path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
                    cur_step = _graph.get_next_step(cur_step);
                }
                // cerr << "path seq:" << path_seq << endl;
                if (backwards) {
                    // cerr << "path seq emplaced (in reverse):" << reverse_complement(path_seq)  << endl;
                    // int init_hap_size = get<0>(haplotypes).size(); // Note: just for debug purposes.
                    get<0>(haplotypes).emplace(reverse_complement(path_seq));
                    // cerr << "was path_seq a new string? " << get<0>(haplotypes).size() - init_hap_size << endl;
                }
                else {
                    // cerr << "path seq emplaced (in forward):" << path_seq  << endl;
                    // int init_hap_size = get<0>(haplotypes).size(); // Note: just for debug purposes.
                    get<0>(haplotypes).emplace(path_seq);
                    // cerr << "was path_seq a copy? " << get<0>(haplotypes).size() - init_hap_size << endl;

                }
            }
        }
        // cerr << "haps in haplotypes: " << endl;
        // for (string hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }
        // Align the new snarl:
        VG new_snarl = align_source_to_sink_haplotypes(get<0>(haplotypes));

        //todo: remove debug:
        
        


        //preprocess new_snarl for log_gbwt_changes:
        bool single_stranded = handlealgs::is_single_stranded(&new_snarl);
        VG* single_stranded_snarl;
        if (!single_stranded) 
        {
            handlealgs::split_strands(&new_snarl, single_stranded_snarl);
            // handlealgs::SplitStrandOverlay(new_snarl)
        }
        else
        {
            single_stranded_snarl=&new_snarl;
        }

        //todo: skipping dagification because I require the input snarl to be a DAG, and I don't think alignments of sequences should produce non-DAGs.
        // bool dag = handlealgs::is_directed_acyclic(&new_snarl);
        // if (!dag)
        // {
        //     handlealgs::dagify(single_stranded_snarl, dagified_snarl, );
        // }
        // else
        // {

        // }

        // count the number of bases in the snarl.
        new_snarl.for_each_handle([&](const handle_t handle) {
            error_record[5] += new_snarl.get_sequence(handle).size();
        });
        force_maximum_handle_size(new_snarl);
        
        // integrate the new_snarl into the _graph, removing the old snarl as you go.
        // //todo: debug_statement
        // integrate_snarl(new_snarl, embedded_paths, sink_id, source_id);
        pair<handle_t, handle_t> new_left_right = integrate_snarl(snarl, new_snarl, embedded_paths, source_id, sink_id, backwards);

        // make a subhandlegraph of the normalized snarl to find the new gbwt paths in the graph.
        SubHandleGraph integrated_snarl = extract_subgraph(_graph, _graph.get_id(new_left_right.first), _graph.get_id(new_left_right.second));

        log_gbwt_changes(source_to_sink_gbwt_paths, integrated_snarl);

        // Print a heads-up about snarls that require an alignment with a greater number of 
        // threads than _big_snarl_alignment_job, so the user knows if they are hung up on a 
        // large job.
        if (get<0>(haplotypes).size() > _big_snarl_alignment_job)
        {
            // Record end time
            auto _big_snarl_time_finish = std::chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = _big_snarl_time_finish - _big_snarl_time_start;
            cerr << "big snarl with " << get<0>(haplotypes).size() << " threads for alignment finished normalization." << endl;
            cerr << "Elapsed time normalizing snarl after sequence extraction: " << elapsed.count() << " s\n";
        }
    } else {
        if (!get<1>(haplotypes).empty()) {
            cerr << "found a snarl with source " << source_id << " and sink "
                 << sink_id
                 << " with haplotypes that start or end in the middle. Skipping." << endl;
            cerr << "There are " << sizeof(get<1>(haplotypes)) << " haplotypes of that description." << endl;
            // vector<string> string_haps = format_handle_haplotypes_to_strings(get<1>(haplotypes).front());
            // cerr << "First example: " << get<1>(haplotypes) << endl;
            error_record[1] = true;
        }
        if (get<0>(haplotypes).size() > _max_alignment_size) {
            cerr << "found a snarl with source " << source_id << " and sink "
                 << sink_id << " with too many haplotypes (" << get<0>(haplotypes).size()
                 << "). Is greater than max, " << _max_alignment_size <<" Skipping." << endl;
            error_record[0] = true;
        }
        // if (get<2>(haplotypes).size() != handles_in_snarl.size()) {
        //     cerr << "some handles in the snarl with source " << source_id
        //          << " and sink " << sink_id
        //          << " aren't accounted for by the gbwt_graph. "
        //             "Skipping."
        //          << endl;
        //     cerr << "handles in snarl:" << handles_in_snarl.size() << "number of handles touched by gbwt graph: " << get<2>(haplotypes).size() << endl;
        //     cerr << "these handles are:" << endl << "\t";
        //     for (auto handle : handles_in_snarl) {
        //         if (get<2>(haplotypes).find(handle) == get<2>(haplotypes).end()) {
        //             cerr << _graph.get_id(handle) << " ";
        //         }
        //     }
        //     cerr << endl;
        //     error_record[2] = true;
        // }
    }
    // todo: decide if we should only normalize snarls that decrease in size.
    if (error_record[5] > error_record[4]) {
        if (_full_log_print)
        {
            cerr << "**************************in UNIT-TEST for normalize_snarl: **************************" << endl;
            cerr << "NOTE: normalized a snarl which *increased* in sequence quantity, "
                    "with source: " << source_id << " and sink: " << sink_id << endl
                << "\tsize before: " << error_record[4] << " size after: " << error_record[5]
                << endl;
        }
    } else if (error_record[5] <= 0) {
        cerr << "normalized snarl size is <= zero: " << error_record[5] << endl;
    }
    return error_record;

}


// Given a vector of haplotypes of format vector< handle_t >, returns a vector of
// haplotypes of
//      format string (which is the concatenated sequences in the handles).
// Arguments:
//      > haplotypes. haplotypte_handle_vectors: a vector of haplotypes in vector<
//      handle_t > format. the handles are from the _gbwt_graph.
// Returns: a vector of haplotypes of format string (which is the concatenated sequences
// in the handles).
unordered_set<string> SnarlNormalizer::format_handle_haplotypes_to_strings(
    const vector<vector<handle_t>> &haplotype_handle_vectors) {
    unordered_set<string> haplotype_strings;
    for (vector<handle_t> haplotype_handles : haplotype_handle_vectors) {
        string hap;
        for (handle_t &handle : haplotype_handles) {
            // hap += _gbwt_graph.get_sequence(handle);
            hap += _graph.get_sequence(_graph.get_handle(_gbwt_graph.get_id(handle), _gbwt_graph.get_is_reverse(handle)));
        }
        haplotype_strings.emplace(hap);
    }
    return haplotype_strings;
}

// TODO: eventually change to deal with haplotypes that start/end in middle of snarl.
// Aligns haplotypes to create a new _graph using MSAConverter's seqan converter.
//      Assumes that each haplotype stretches from source to sink.
// Arguments:
//      source_to_sink_haplotypes: a vector of haplotypes in string format (concat of
//      handle sequences).
// Returns:
//      VG object representing the newly realigned snarl.
VG SnarlNormalizer::align_source_to_sink_haplotypes(
    const unordered_set<string>& source_to_sink_haplotypes) {
    // cerr << "align_source_to_sink_haplotypes" << endl;
    // cerr << " haplotypes in source_to_sink_haplotypes: " << endl;
    // for (string hap : source_to_sink_haplotypes) {
    //     cerr << hap << endl;
    // }
    // cerr << "number of strings to align: " << source_to_sink_haplotypes.size() << endl;
    // TODO: make the following comment true, so that I can normalize haplotypes that
    // TODO:    aren't source_to_sink by adding a similar special character to strings in
    // TODO:    the middle of the snarl.
    // modify source_to_sink_haplotypes to replace the leading and
    // trailing character with a special character. This ensures that the leading char of
    // the haplotype becomes the first character in the newly aligned snarl's source - it
    // maintains the context of the snarl.

    // store the source/sink chars for later reattachment to source and sink.
    string random_element;
    for (auto hap : source_to_sink_haplotypes){
        random_element = hap;
        break;
    }
    string source_char(1, random_element.front());
    string sink_char(1, random_element.back());

    // cerr << "strings in path_seq before replacing final character: " << endl;
    // for (auto path : get<0>(haplotypes))
    // {
    //     cerr << path << endl;f
    // }

    // replace the source and sink chars with X, to force match at source and sink.
    unordered_set<string> edited_source_to_sink_haplotypes;
    // for (auto it = source_to_sink_haplotypes.begin(); it != source_to_sink_haplotypes.end(); it++)
    for (auto hap : source_to_sink_haplotypes)
    {
        // cerr << "hap before replace: " << hap << endl;
        hap.replace(0, 1, "X");
        hap.replace(hap.size() - 1, 1, "X");
        // cerr << "hap after replace: " << hap << endl;
        edited_source_to_sink_haplotypes.emplace(hap);
    }
    // cerr << "source_char: " << source_char << endl;
    // cerr << "sink_char: " << sink_char << endl;

    // //todo: debug_statement
    // source_to_sink_haplotypes.emplace_back("XX");

    // /// make a new scoring matrix with _match=5, _mismatch = -3, _gap_extend = -1, and
    // _gap_open = -3, EXCEPT that Q has to be matched with Q (so match score between Q
    // and Q =len(seq)+1)
    // // 1. Define type and constants.
    // //
    // // Define types for the score value and the scoring scheme.
    // typedef int TValue;
    // typedef seqan::Score<TValue, seqan::ScoreMatrix<seqan::Dna5, seqan::Default> >
    // TScoringScheme;
    // // Define our gap scores in some constants.
    // int const gapOpenScore = -1;
    // int const gapExtendScore = -1;

    // static int const _data[TAB_SIZE] =
    //     {
    //         1, 0, 0, 0, 0,
    //         0, 1, 0, 0, 0,
    //         0, 0, 1, 0, 0,
    //         0, 0, 0, 1, 0,
    //         0, 0, 0, 0, 0
    //     };

    // create seqan multiple_sequence_alignment object
    //// seqan::Align<seqan::DnaString>   align;
    seqan::Align<seqan::CharString> align;

    seqan::resize(rows(align), edited_source_to_sink_haplotypes.size());
    int i = 0;
    for (auto hap : edited_source_to_sink_haplotypes) {
        assignSource(row(align, i), hap.c_str());
        i++;
    }

    globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));

    vector<string> row_strings;
    for (auto &row : rows(align)) {
        string row_string;
        auto it = begin(row);
        auto itEnd = end(row);
        for (; it != itEnd; it++) {
            row_string += *it;
        }
        // todo: debug_statement
        // cerr << "ROW_STRING: " << row_string << endl;
        // edit the row so that the proper source and sink chars are added to the
        // haplotype instead of the special characters added to ensure correct alignment
        // of source and sink.
        // cerr << "row_string before: " << row_string << endl;
        row_string.replace(0, 1, source_char);
        row_string.replace(row_string.size() - 1, 1, sink_char);
        row_strings.push_back(row_string);
        // cerr << "row_string after: " << row_string << endl;
    }

    stringstream ss;
    for (string seq : row_strings) {
        // todo: debug_statement
        // cerr << "seq in alignment:" << seq << endl;
        ss << endl << seq;
    }
    // ss << align;
    MSAConverter myMSAConverter = MSAConverter();
    myMSAConverter.load_alignments(ss, "seqan");
    VG snarl = myMSAConverter.make_graph();
    snarl.clear_paths();

    pair<vector<handle_t>, vector<handle_t>> source_and_sink =
        debug_get_sources_and_sinks(snarl);


    // TODO: throw exception(?) instead of cerr, or remove these messages if I'm confident
    // TODO:    code works.
    if (source_and_sink.first.size() != 1) {
        cerr << "WARNING! Snarl realignment has generated "
             << source_and_sink.first.size() << " source nodes." << endl;
    }

    if (source_and_sink.second.size() != 1) {
        cerr << "WARNING! Snarl realignment has generated "
             << source_and_sink.second.size() << " sink nodes." << endl;
    }

    return snarl;
}

/** For each handle in a given _graph, divides any handles greater than max_size into
 * parts that are equal to or less than the size of max_size.
 *
 * @param  {MutableHandleGraph} _graph : the _graph in which we want to force a maximum
 * handle size for all handles.
 * @param  {size_t} max_size          : the maximum size we want a handle to be.
 */
void SnarlNormalizer::force_maximum_handle_size(MutableHandleGraph &graph) {
    // forcing each handle in the _graph to have a maximum sequence length of max_size:
    graph.for_each_handle([&](handle_t handle) {
        // all the positions we want to make in the handle are in offsets.
        vector<size_t> offsets;

        size_t sequence_len = graph.get_sequence(handle).size();
        int number_of_divisions = floor(sequence_len / _max_handle_size);

        // if the handle divides evenly into subhandles of size _max_handle_size, we don't need to
        // make the last cut (which would be at the very end of the handle - cutting off
        // no sequence).
        if (sequence_len % _max_handle_size == 0) {
            number_of_divisions--;
        }

        // calculate the position of all the divisions we want to make.
        for (int i = 1; i <= number_of_divisions; i++) {
            offsets.push_back(i * _max_handle_size);
        }

        // divide the handle into parts.
        graph.divide_handle(handle, offsets);
    });
}

// TODO: change the arguments to handles, which contain orientation within themselves.
// Given a start and end node id, construct an extract subgraph between the two nodes
// (inclusive). Arguments:
//      _graph: a pathhandlegraph containing the snarl with embedded paths.
//      source_id: the source of the snarl of interest.
//      sink_id: the sink of the snarl of interest.
// Returns:
//      a SubHandleGraph containing only the handles in _graph that are between start_id
//      and sink_id.
SubHandleGraph SnarlNormalizer::extract_subgraph(const HandleGraph &graph,
                                                 const id_t &leftmost_id,
                                                 const id_t &rightmost_id) {
    // cerr << "extract_subgraph has source and sink: " << source_id << " " << sink_id << endl; 
    // because algorithm moves left to right, determine leftmost and rightmost nodes.
    // id_t leftmost_id;
    // id_t rightmost_id;
    //// if snarl's "backwards," source is rightmost node, sink is leftmost.
    // if (backwards) 
    // {
    //     leftmost_id = sink_id;
    //     rightmost_id = source_id;
    // }
    // else 
    // {
    //     leftmost_id = source_id;
    //     rightmost_id = sink_id;
    // }
    // cerr << "extract_subgraph" << endl;
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    unordered_set<id_t> visited;  // to avoid counting the same node twice.
    unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

    // initialize with leftmost_handle (because we move only to the right of leftmost_handle):
    handle_t leftmost_handle = _graph.get_handle(leftmost_id);
    subgraph.add_handle(leftmost_handle);
    visited.insert(graph.get_id(leftmost_handle));

    // look only to the right of leftmost_handle
    _graph.follow_edges(leftmost_handle, false, [&](const handle_t &handle) {
        // mark the nodes to come as to_visit
        if (visited.find(graph.get_id(handle)) == visited.end()) {
            to_visit.insert(graph.get_id(handle));
        }
    });

    /// explore the rest of the snarl:
    while (to_visit.size() != 0) {
        // remove cur_handle from to_visit
        unordered_set<id_t>::iterator cur_index = to_visit.begin();
        handle_t cur_handle = _graph.get_handle(*cur_index);

        to_visit.erase(cur_index);

        /// visit cur_handle
        visited.insert(graph.get_id(cur_handle));

        subgraph.add_handle(cur_handle);

        if (graph.get_id(cur_handle) != rightmost_id) { // don't iterate past rightmost node!
            // look for all nodes connected to cur_handle that need to be added
            // looking to the left,
            _graph.follow_edges(cur_handle, true, [&](const handle_t &handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
            // looking to the right,
            _graph.follow_edges(cur_handle, false, [&](const handle_t &handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
        }
    }
    return subgraph;
}

// Integrates the snarl into the _graph, replacing the snarl occupying the space between
// source_id and sink_id.
//      In the process, transfers any embedded paths traversing the old snarl into the new
//      snarl.
// Arguments:
//      _graph: the _graph in which we want to insert the snarl.
//      to_insert_snarl: a *separate* handle_graph from _graph, often generated from
//      MSAconverter. embedded_paths: a vector of paths, where each is a pair.
//                        pair.first is the first step_handle of interest in the
//                        old_embedded_path, and pair.second is the step_handle *after*
//                        the last step_handle of interest in the old_embedded_path (can
//                        be the null step at the end of the path.)
//                        Note: these paths will be altered to represent the way they
//                        overlap in the new snarl. Otherwise, they would be invalidated.
//      source_id: the source of the old (to be replaced) snarl in _graph
//      sink_id: the sink of the old (to be replaced) snarl in _graph.
// Return: a pair of node ids, representing source and sink of the newly integrated snarl.
pair<handle_t, handle_t> SnarlNormalizer::integrate_snarl(SubHandleGraph &old_snarl, 
    const HandleGraph &to_insert_snarl,
    vector<pair<step_handle_t, step_handle_t>>& embedded_paths, 
    const id_t &source_id, const id_t &sink_id, const bool backwards) {
    // cerr << "integrate_snarl" << endl;

    //todo: debug_statement
    // cerr << "\nhandles in to_insert_snarl:" << endl;
    // to_insert_snarl.for_each_handle([&](const handle_t &handle) {
    //     cerr << to_insert_snarl.get_id(handle) << " "
    //          << to_insert_snarl.get_sequence(handle) << " ";
    //     cerr << "neighbors: ";
    //     to_insert_snarl.follow_edges(handle, false, [&](const handle_t &next) {
    //         cerr << "     " << to_insert_snarl.get_id(next) << endl;
    //     });
    //     cerr << " \n";
    // });
    // cerr << endl;
    // Get old _graph snarl
    // SubHandleGraph old_snarl = extract_subgraph(_graph, source_id, sink_id, backwards);

    // TODO: debug_statement: Check to make sure that newly made snarl has only one start
    // and end.
    // TODO:     (shouldn't be necessary once we've implemented alignment with
    // leading/trailing special chars.) Identify old and new snarl start and sink
    pair<vector<handle_t>, vector<handle_t>> to_insert_snarl_defining_handles =
        debug_get_sources_and_sinks(to_insert_snarl);

    if (to_insert_snarl_defining_handles.first.size() > 1 ||
        to_insert_snarl_defining_handles.second.size() > 1) {
        cerr << "ERROR: newly made snarl from a snarl with source " << source_id
             << " has more than one start or end. # of starts: "
             << to_insert_snarl_defining_handles.first.size()
             << " # of ends: " << to_insert_snarl_defining_handles.second.size() << endl;
        exit(1);
    }

    /// Replace start and end handles of old _graph snarl with to_insert_snarl start and
    /// end, and delete rest of old _graph snarl:

    // add to_insert_snarl into _graph without directly attaching the snarl to the _graph
    // (yet).
    vector<handle_t> to_insert_snarl_topo_order =
        handlealgs::lazier_topological_order(&to_insert_snarl);

    // Construct a parallel new_snarl_topo_order to identify
    // paralogous nodes between to_insert_snarl and the new snarl inserted in _graph.
    vector<handle_t> new_snarl_topo_order;

    // integrate the handles from to_insert_snarl into the _graph, and keep track of their
    // identities by adding them to new_snarl_topo_order.
    for (handle_t to_insert_snarl_handle : to_insert_snarl_topo_order) {
        // //todo: debug_statement:
        // cerr << "About to insert snarl handle from normalized graph of id, seq: "
        //      << to_insert_snarl.get_id(to_insert_snarl_handle) << " "
        //      << to_insert_snarl.get_sequence(to_insert_snarl_handle) << endl;

        handle_t graph_handle =
            _graph.create_handle(to_insert_snarl.get_sequence(to_insert_snarl_handle));
        // cerr << "here is the new snarl handle: " 
        //      << _graph.get_id(graph_handle) << " "
        //      << _graph.get_sequence(graph_handle) << endl;
        new_snarl_topo_order.push_back(graph_handle);
    }
    // cerr << "finished inserting the snarls from to_insert_snarl into normalized graph." << endl;

    // Connect the newly made handles in the _graph together the way they were connected
    // in to_insert_snarl:
    for (int i = 0; i < to_insert_snarl_topo_order.size(); i++) {
        to_insert_snarl.follow_edges(
            to_insert_snarl_topo_order[i], false, [&](const handle_t &snarl_handle) {
                // get topo_index of nodes to be connected to _graph start handle
                auto it = find(to_insert_snarl_topo_order.begin(),
                               to_insert_snarl_topo_order.end(), snarl_handle);
                int topo_index = it - to_insert_snarl_topo_order.begin();

                // connect _graph start handle
                _graph.create_edge(new_snarl_topo_order[i],
                                   new_snarl_topo_order[topo_index]);
            });
    }


    // save the source and sink values of new_snarl_topo_order, since topological order is
    // not necessarily preserved by move_path_to_snarl. Is temporary b/c we need to
    // replace the handles with ones with the right id_t label for source and sink later
    // on.
    id_t temp_snarl_leftmost_id = _graph.get_id(new_snarl_topo_order.front());
    id_t temp_snarl_rightmost_id = _graph.get_id(new_snarl_topo_order.back());
    if (new_snarl_topo_order.size() == 1)
    {
        // in case the normalized snarl is only one handle in size, split it into two.
        // This allows the front to be renamed after the source, and the end after the sink.
        std::pair<handle_t, handle_t> split_handle = _graph.divide_handle(new_snarl_topo_order.back(), 1);
        temp_snarl_leftmost_id = _graph.get_id(split_handle.first);
        temp_snarl_rightmost_id = _graph.get_id(split_handle.second);
    }
    // cerr << "the temp source id: " << temp_snarl_leftmost_id << endl;
    // cerr << "the temp sink id: " << temp_snarl_rightmost_id << endl;

    // Add the neighbors of the source and sink of the original snarl to the new_snarl's
    // source and sink.
    // source integration:
    if (!backwards)
    {
    _graph.follow_edges(
        _graph.get_handle(source_id), true, [&](const handle_t &prev_handle) {
            _graph.create_edge(prev_handle, _graph.get_handle(temp_snarl_leftmost_id));
        });
    _graph.follow_edges(
        _graph.get_handle(sink_id), false, [&](const handle_t &next_handle) {
            _graph.create_edge(_graph.get_handle(temp_snarl_rightmost_id), next_handle);
        });
    }
    else 
    {
        _graph.follow_edges(
        _graph.get_handle(source_id), false, [&](const handle_t &next_handle) {
            _graph.create_edge(_graph.get_handle(temp_snarl_rightmost_id), next_handle);
        });
    _graph.follow_edges(
        _graph.get_handle(sink_id), true, [&](const handle_t &prev_handle) {
            _graph.create_edge(prev_handle, _graph.get_handle(temp_snarl_leftmost_id));
        });
    }
    // For each path of interest, move it onto the new_snarl.
    for (int i = 0; i != embedded_paths.size(); i++)
    {
        // //todo: debug_statement
        // cerr << "the new sink id: " << temp_snarl_rightmost_id << endl;
        // //todo: debug_statement
        // move_path_to_snarl(path, new_snarl_topo_order, temp_snarl_rightmost_id,
        //                    temp_snarl_leftmost_id, sink_id, source_id);
        // move_path_to_snarl(path, new_snarl_topo_order, temp_snarl_leftmost_id,
        //                    temp_snarl_rightmost_id, source_id, sink_id, backwards);
        // cerr << "is path backwards? " << backwards << endl;
        // cerr << "path first: " << _graph.get_id(_graph.get_handle_of_step(path.first)) << " step after path first: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_next_step(path.first))) << " path second: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << endl;
        // cerr << "source: " << source_id << " sink: " << sink_id << endl;
        // pair<bool, bool> path_spans_left_right;
        // path_spans_left_right.first = (!backwards && _graph.get_id(_graph.get_handle_of_step(path.first)) == source_id) || (backwards && _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == source_id);
        // path_spans_left_right.second = (!backwards && _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == sink_id) || (backwards && _graph.get_id(_graph.get_handle_of_step(path.first)) == sink_id);
        // cerr << "first: " << path_spans_left_right.first << "second: " << path_spans_left_right.second << endl;
        pair<bool, bool> path_spans_left_right;
        path_spans_left_right.first = (_graph.get_id(_graph.get_handle_of_step(embedded_paths[i].first)) == source_id);
        path_spans_left_right.second = (_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(embedded_paths[i].second))) == sink_id);

        embedded_paths[i] = move_path_to_new_snarl(embedded_paths[i], temp_snarl_leftmost_id, temp_snarl_rightmost_id, path_spans_left_right, !backwards, make_pair(source_id, sink_id));
    }

    // Destroy the old snarl.
    old_snarl.for_each_handle([&](const handle_t &handle) 
    {
        // //todo: debug_statement these are the handles in old_snarl:
        // cerr << "destroying old_snarl handle: " << old_snarl.get_id(handle) << " with sequence: " << old_snarl.get_sequence(handle) << endl;
        _graph.destroy_handle(handle);
    });

    // Replace the source and sink handles with ones that have the original source/sink id
    // (for compatibility with future iterations on neighboring top-level snarls using the
    // same snarl manager. Couldn't replace it before b/c we needed the old handles to
    // move the paths.
    handle_t new_leftmost_handle;
    handle_t new_rightmost_handle;
    if (!backwards) 
    {
        // cerr << "!backwards" << endl;
        // cerr << "overwriting node id " << temp_snarl_leftmost_id <<  " with " << source_id << " (which is source_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_leftmost_id)) << endl;
        new_leftmost_handle = overwrite_node_id(temp_snarl_leftmost_id, source_id);
        // cerr << "overwriting node id " << temp_snarl_rightmost_id <<  " with " << sink_id << " (which is sink_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_rightmost_id)) << endl;
        new_rightmost_handle = overwrite_node_id(temp_snarl_rightmost_id, sink_id);
    }
    else
    {
        // cerr << "backwards" << endl;
        // cerr << "overwriting node id " << temp_snarl_leftmost_id <<  " with " << sink_id << " (which is sink_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_leftmost_id)) << endl;
        new_leftmost_handle = overwrite_node_id(temp_snarl_leftmost_id, sink_id);
        // cerr << "overwriting node id " << temp_snarl_rightmost_id <<  " with " << source_id << " (which is source_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_rightmost_id)) << endl;
        new_rightmost_handle = overwrite_node_id(temp_snarl_rightmost_id, source_id);
    }    
    pair<handle_t, handle_t> new_left_right = make_pair(new_leftmost_handle, new_rightmost_handle);
    return new_left_right;
}


/**
 * Deletes the given handle's underlying node, and returns a new handle to a new node 
 * with the desired node_id
 * 
 * @param  {id_t} handle     : The old node id, to be replaced with a new node id.
 * @param  {id_t} node_id    : The node id for the new node. Cannot be currently in use in
 *                              the graph.
 * @return {handle_t}        : The new handle, in the same position as the original handle
 *                              in the graph, but with the new node_id.
 */
handle_t SnarlNormalizer::overwrite_node_id(const id_t& old_node_id, const id_t& new_node_id)
{
    handle_t old_handle = _graph.get_handle(old_node_id);
    handle_t new_handle = _graph.create_handle(_graph.get_sequence(old_handle), new_node_id);

    // move the edges:
    _graph.follow_edges(old_handle, true, [&](const handle_t &prev_handle) 
    {
        _graph.create_edge(prev_handle, new_handle);
    });
    _graph.follow_edges(old_handle, false, [&](const handle_t &next_handle)
    {
        _graph.create_edge(new_handle, next_handle);
    });

    // move the paths:
    _graph.for_each_step_on_handle(old_handle, [&](step_handle_t step) 
    {
        handle_t properly_oriented_old_handle = _graph.get_handle_of_step(step); 
        if (_graph.get_is_reverse(properly_oriented_old_handle) != _graph.get_is_reverse(new_handle))
        {
            new_handle = _graph.flip(new_handle);
        }
        _graph.rewrite_segment(step, _graph.get_next_step(step), vector<handle_t>{new_handle});
    });

    // delete the old_handle:
    _graph.destroy_handle(old_handle);
    return new_handle;
}

/**
 * Updates the changes that need making to the gbwt after the graph is finished being
 * normalized, so that an updated gbwt can be made.
 * @param  {list<string>} old_paths : the paths in the gbwt that need moving to the new
 * graph.
 * @param  {HandleGraph} new_snarl  : the normalized portion of the graph. Probably a 
 * subhandlegraph.
 */
void SnarlNormalizer::log_gbwt_changes(const vector<pair<gbwt::vector_type, string>>& source_to_sink_gbwt_paths, const HandleGraph &new_snarl){
    //todo: move Aligner to initialization of object, since I'm not supposed to make a new one each time I do alignments.
    Aligner aligner = Aligner();
    // cerr << "in log_gbwt_changes" << endl;
    // cerr << old_paths.size() << endl;
    for (auto path : source_to_sink_gbwt_paths)
    {
        Alignment alignment;
        alignment.set_sequence(path.second);
        aligner.align_global_banded(alignment, new_snarl,0, false);
        // cerr << "ALIGNMENT PATH FOR " << path << ":" << endl;
        gbwt::vector_type alignment_full_path;
        for (auto mapping : alignment.path().mapping())
        {
            // gbwt::Node::encode(id, is_reverse)
            // mapping.position().
            alignment_full_path.emplace_back(gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse()));
            // cerr << "mapping.position().node_id() " << mapping.position().node_id() << endl;
        }
        _gbwt_changelog.emplace_back(path.first, alignment_full_path);
        // cerr << pb2json(alignment.path()) << endl << alignment.query_position() << endl << alignment.path().mapping().begin() << endl << endl;
        // alignment.path().mapping()
    }
    
    // use banded global aligner. optimizations for finidng one perfect match from source to sink.

}




/** Used to help move_path_to_snarl map paths from an old snarl to its newly
 * normalized counterpart. In particular, ensures that any paths which touch the
 * source and/or sink of the old snarl still do so in the new snarl (which is
 * important to ensure that we don't break any paths partway through the snarl.)
 *
 * @param  {HandleGraph} _graph         : the _graph that contains the old and new snarl
 * nodes.
 * @param  {id_t} new_source_id        : the node id of the newly created source.
 * @param  {id_t} new_sink_id          : the node id of the newly created sink.
 * @param  {bool} touching_source      : true if the path is connected to the old
 * source.
 * @param  {bool} touching_sink        : true if the path is connected to the old
 * sink.
 * @param  {handle_t} path_start : proposed source for the path in the new snarl.
 * @param  {handle_t} path_end   : proposed sink for the path in the new snarl.
 * @return {bool}                      : true if the path satisfies the requirement
 * that, if the original path covered the old source or sink, the new path also covers
 * the same respective nodes in the new snarl.
 */
bool SnarlNormalizer::source_and_sink_handles_map_properly(
    const HandleGraph &graph, const id_t &new_source_id, const id_t &new_sink_id,
    const bool &touching_source, const bool &touching_sink, const handle_t &path_start,
    const handle_t &path_end) {

    bool path_map = false;
    // cerr << "touching source? " << touching_source << "touching_sink" << touching_sink
    //      << "source is source?" << (graph.get_id(path_start) == new_source_id)
    //      << " sink is sink: " << (graph.get_id(path_end) == new_sink_id) << endl;
    if (touching_source && touching_sink) {
        path_map = ((graph.get_id(path_start) == new_source_id) &&
                    (graph.get_id(path_end) == new_sink_id));
    } else if (touching_source) {
        path_map = (graph.get_id(path_start) == new_source_id);
    } else if (touching_sink) {
        path_map = (graph.get_id(path_end) == new_sink_id);
    } else {
        path_map = true;
    }
    // cerr << "path_map " << path_map << endl;
    return path_map;
}

// Determines whether some subsequence in a handle satisfies the condition of being
// the beginning of a path.
//      If the path_seq is longer than the handle_seq, only checks subsequences that
//      reach from the beginning/middle of the handle_seq to the end. If path_seq is
//      shorter than handle_seq, checks for any substring of length path_seq within
//      the handle_seq, as well as substrings smaller than length path_seq that extend
//      beyond the current handle.
// Arguments:
//      handle_seq: the sequence in the handle we're trying to identify as a
//      start_of_path_seq. path_seq: the sequence in the path we're trying to find
//      starting points for in handle_seq
// Return: a vector of all potential starting index of the subsequence in the
// handle_seq.
vector<int> SnarlNormalizer::check_handle_as_start_of_path_seq(const string &handle_seq,
                                                               const string &path_seq) {
    vector<int> possible_start_indices;
    // If the handle_seq.size <= path_seq.size, look for subsequences reaching from
    // beginning/middle of handle_seq to the end - where path_seq may run off the end
    // of this handle to the next in the snarl.
    if (handle_seq.size() <= path_seq.size()) {
        // iterate through all possible starting positions in the handle_seq.
        for (int handle_start_i = 0; handle_start_i < handle_seq.size();
             handle_start_i++) {
            int subseq_size = handle_seq.size() - handle_start_i;
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
    }
    // if handle_seq.size > path_seq.size, look for any subsequence within handle_seq
    // of path_seq.size, as well as any subsequence smaller than path_seq reaching
    // from middle of handle_seq to the end of handle_seq.
    else {
        // first, search through all handle_seq for any comparable subsequence of
        // path_seq.size. Note: only differences between this for loop and above for
        // loop is that handle_start_i stops at (<= path_seq.size() -
        // handle_seq.size()), and subseq.size() = path_seq.size()
        for (int handle_start_i = 0;
             handle_start_i <= (handle_seq.size() - path_seq.size()); handle_start_i++) {
            int subseq_size = path_seq.size();
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
        // second, search through the last few bases of handle_seq for the beginning
        // of path_seq. Note: nearly identical for loop to the one in "if
        // (handle_seq.size()
        // <= path_seq.size())"
        for (int handle_start_i = (handle_seq.size() - path_seq.size() + 1);
             handle_start_i < handle_seq.size(); handle_start_i++) {
            int subseq_size = handle_seq.size() - handle_start_i;
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
    }
    // Note: if we passed through the above check without returning anything, then
    // there isn't any satisfactory subsequence and we'll return an empty vector.
    return possible_start_indices;
}

// ------------------------------ DEBUG CODE BELOW:
// ------------------------------------------

// Returns pair where pair.first is a vector of all sources of the given _graph and
// path.second is all the sinks of the given _graph. If _graph is a subhandlegraph of a
// snarl, there should only be one source and sink each.
pair<vector<handle_t>, vector<handle_t>>
SnarlNormalizer::debug_get_sources_and_sinks(const HandleGraph &graph) {
    // cerr << "debug_get_source_and_sinks" << endl;
    vector<handle_t> sink;
    vector<handle_t> source;

    // identify sources and sinks
    graph.for_each_handle([&](const handle_t &handle) {
        //todo: debug_statements in code below:
        // cerr << "identifying if " << graph.get_id(handle) << "is a source/sink." <<endl;
        bool is_source = true, is_sink = true;
        // cerr << "handles to the left: ";
        graph.follow_edges(handle, true, [&](const handle_t &prev) {
            // cerr << graph.get_id(prev) << endl;
            is_source = false;
            return false;
        });
        // cerr << "handles to the right: ";
        graph.follow_edges(handle, false, [&](const handle_t &next) {
            // cerr << graph.get_id(next) << endl;
            is_sink = false;
            return false;
        });

        if (is_source) {
            // cerr<< "determined is_source" << endl;
            source.push_back(handle);
        }
        if (is_sink) {
            // cerr<< "determined is_sink" << endl;
            sink.emplace_back(handle);
        }
    });
    return pair<vector<handle_t>, vector<handle_t>>(source, sink);
}

}
}