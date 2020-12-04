#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?

#include "0_oo_normalize_snarls.hpp"
#include "0_snarl_sequence_finder.hpp"
// #include <algorithm>
#include <string>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

#include <gbwtgraph/gbwtgraph.h>

#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../msa_converter.hpp"
#include "../snarls.hpp"
#include "../vg.hpp"
#include "is_acyclic.hpp"

#include "../types.hpp"
#include "extract_containing_graph.hpp"

/*
TODO: allow for snarls that have haplotypes that begin or end in the middle of the snarl

TODO: allow normalization of multiple adjacent snarls in one combined realignment.

TODO: test that extract_gbwt haplotypes successfully extracts any haplotypes that start/end in the middle of
TODO:    the snarl.
*/
namespace vg {
namespace algorithms{
/**
 * To "normalize" a snarl, SnarlNormalizer extracts all the sequences in the snarl as
 * represented in the gbwt, and then realigns them to create a replacement snarl. 
 * This process hopefully results in a snarl with less redundant sequence, and with 
 * duplicate variation combined into a single variant.
*/
SnarlNormalizer::SnarlNormalizer(MutablePathDeletableHandleGraph &graph,
                                 const gbwtgraph::GBWTGraph &haploGraph,
                                 const int &max_alignment_size, const string &path_finder)
    : _haploGraph(haploGraph), _graph(graph), _max_alignment_size(max_alignment_size),
      _path_finder(path_finder) {}


/**
 * Iterates over all top-level snarls in _graph, and normalizes them.
 * @param snarl_stream file stream from .snarl.pb output of vg snarls
*/
void SnarlNormalizer::normalize_top_level_snarls(ifstream &snarl_stream) {
    cerr << "disambiguate_top_level_snarls" << endl;
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
    */ 
    int error_record_size = 5;
    vector<int> one_snarl_error_record(error_record_size, 0);
    vector<int> full_error_record(error_record_size, 0);

    pair<int, int> snarl_sequence_change;

    for (auto roots : snarl_roots) {
        // if (roots->start().node_id() > 269600 && roots->start().node_id() < 269700) {
            cerr << "disambiguating snarl #"
                 << (num_snarls_normalized + num_snarls_skipped)
                 << " source: " << roots->start().node_id()
                 << " sink: " << roots->end().node_id() << endl;

            normalize_snarl(roots->start().node_id(), roots->end().node_id());

            if (!((one_snarl_error_record[0]) || (one_snarl_error_record[1]) ||
                  (one_snarl_error_record[2]) || (one_snarl_error_record[3]))) {
                // if there are no errors, then we've successfully normalized a snarl.
                num_snarls_normalized += 1;
                // track the change in size of the snarl.
                snarl_sequence_change.first += one_snarl_error_record[4];
                snarl_sequence_change.second += one_snarl_error_record[5];
            } else {
                // else, there was an error. Track which errors caused the snarl to not
                // normalize.
                // note: the last two ints are ignored here b/c they're for
                // recording the changing size of snarls that are successfully normalized.
                for (int i = 0; i < error_record_size - 2; i++) {
                    full_error_record[i] += one_snarl_error_record[i];
                }
                num_snarls_skipped += 1;
            }
        // }
    }
    cerr << endl
         << "normalized " << num_snarls_normalized << " snarl(s), skipped "
         << num_snarls_skipped << " snarls because. . .\nthey exceeded the size limit ("
         << full_error_record[0]
         << "snarls),\nhad haplotypes starting/ending in the middle of the snarl ("
         << full_error_record[1] << "),\nthe snarl was cyclic (" << full_error_record[3]
         << " snarls),\nor there "
            "were handles not connected by the gbwt info ("
         << full_error_record[2] << " snarls)." << endl;
    cerr << "amount of sequence in normalized snarls before normalization: "
         << snarl_sequence_change.first << endl;
    cerr << "amount of sequence in normalized snarls after normalization: "
         << snarl_sequence_change.second << endl;

    // //todo: debug_statement for extracting snarl of interest.
    // VG outGraph;
    // pos_t source_pos = make_pos_t(269695, false, 0);
    // vector<pos_t> pos_vec;
    // pos_vec.push_back(source_pos);
    // algorithms::extract_containing_graph(&_graph, &outGraph, pos_vec, 1000);
    // outGraph.serialize_to_ostream(cout);

    delete snarl_manager;
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
vector<int> SnarlNormalizer::normalize_snarl(const id_t &source_id, const id_t &sink_id) {
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
    */ 
    vector<int> error_record(6, 0);
    SubHandleGraph snarl = extract_subgraph(_graph, source_id, sink_id);

    if (!algorithms::is_acyclic(&snarl)) {
        cerr << "snarl at " << source_id << " is cyclic. Skipping." << endl;
        error_record[3] = true;
        return error_record;
    }

    // extract threads
    tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<handle_t>> haplotypes;
    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _haploGraph, source_id, sink_id);
    
    if (_path_finder == "GBWT") {
        tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();
        // Convert the haplotypes from vector<handle_t> format to string format.
        get<0>(haplotypes) = format_handle_haplotypes_to_strings(get<0>(gbwt_haplotypes));
        get<1>(haplotypes) = get<1>(gbwt_haplotypes);
        get<2>(haplotypes) = get<2>(gbwt_haplotypes);
    } else if (_path_finder == "exhaustive") {
        pair<unordered_set<string>, unordered_set<handle_t>> exhaustive_haplotypes =
            sequence_finder.find_exhaustive_paths();
        get<0>(haplotypes) = exhaustive_haplotypes.first;
        get<2>(haplotypes) = exhaustive_haplotypes.second;
    } else {
        cerr << "path_finder type must be 'GBWT' or 'exhaustive', not " << _path_finder
             << endl;
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

    // TODO: this if statement removes snarls where a haplotype begins/ends in the middle
    // TODO:    of the snarl. Get rid of this once alignment issue is addressed!
    // TODO: also, limits the number of haplotypes to be aligned, since snarl starting at
    // TODO:    2049699 with 258 haplotypes is taking many minutes.
    if (get<1>(haplotypes).empty() && get<0>(haplotypes).size() < _max_alignment_size &&
        get<2>(haplotypes).size() == handles_in_snarl.size()) {
        // Get the embedded paths in the snarl from _graph, to move them to new_snarl.
        // Any embedded paths not in gbwt are aligned in the new snarl.
        vector<pair<step_handle_t, step_handle_t>> embedded_paths =
            sequence_finder.find_embedded_paths();

        //todo: debug_statement
        cerr << "Let's see what sequences I have before adding embedded paths to seq info:" << endl;
        for (string seq : get<0>(haplotypes)) {
            cerr << seq << endl;
        }

        // TODO: once haplotypes that begin/end in the middle of the snarl have been
        // TODO:    accounted for in the code, remove next chunk of code that finds 
        // TODO: source-to-sink paths.
        // find the paths that stretch from source to sink:
        for (auto path : embedded_paths) {
            // cerr << "checking path of name " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << " with start " << _graph.get_id(_graph.get_handle_of_step(path.first)) << " and sink " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << endl;
            if (_graph.get_id(_graph.get_handle_of_step(path.first)) == source_id &&
                _graph.get_id(_graph.get_handle_of_step(
                    _graph.get_previous_step(path.second))) == sink_id) {
                // cerr << "adding path of name " <<
                // _graph.get_path_name(graph.get_path_handle_of_step(path.first)) <<
                // endl; get the sequence of the source to sink path, and add it to the
                // paths to be aligned.
                string path_seq;
                step_handle_t cur_step = path.first;
                while (cur_step != path.second) {
                    path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
                    cur_step = _graph.get_next_step(cur_step);
                }
                get<0>(haplotypes).emplace(path_seq);
            }
        }
        // Align the new snarl:
        VG new_snarl = align_source_to_sink_haplotypes(get<0>(haplotypes));

        // count the number of bases in the snarl.
        new_snarl.for_each_handle([&](const handle_t handle) {
            error_record[5] += new_snarl.get_sequence(handle).size();
        });

        force_maximum_handle_size(new_snarl, _max_alignment_size);

        // integrate the new_snarl into the _graph, removing the old snarl as you go.
        integrate_snarl(new_snarl, embedded_paths, source_id, sink_id);
    } else {
        if (!get<1>(haplotypes).empty()) {
            cerr << "found a snarl starting at " << source_id << " and ending at "
                 << sink_id
                 << " with haplotypes that start or end in the middle. Skipping." << endl;
            error_record[1] = true;
        }
        if (get<0>(haplotypes).size() > _max_alignment_size) {
            cerr << "found a snarl starting at " << source_id << " and ending at "
                 << sink_id << " with too many haplotypes (" << get<0>(haplotypes).size()
                 << ") to efficiently align. Skipping." << endl;
            error_record[0] = true;
        }
        if (get<2>(haplotypes).size() != handles_in_snarl.size()) {
            cerr << "some handles in the snarl starting at " << source_id
                 << " and ending at " << sink_id
                 << " aren't accounted for by the gbwt_graph. "
                    "Skipping."
                 << endl;
            cerr << "these handles are:" << endl << "\t";
            for (auto handle : handles_in_snarl) {
                if (get<2>(haplotypes).find(handle) == get<2>(haplotypes).end()) {
                    cerr << _graph.get_id(handle) << " ";
                }
            }
            cerr << endl;
            error_record[2] = true;
        }
    }
    // todo: decide if we should only normalize snarls that decrease in size.
    if (error_record[5] > error_record[4]) {
        cerr << "NOTE: normalized a snarl which *increased* in sequence quantity, "
                "starting at "
             << source_id << endl
             << "\tsize before: " << error_record[4] << " size after: " << error_record[5]
             << endl;
    } else if (error_record[5] <= 0) {
        cerr << "normalized snarl size is <= zero: " << error_record[5] << endl;
    }
    return error_record;

}

// Given a vector of haplotypes of format vector< handle_t >, returns a vector of
// haplotypes of
//      format string (which is the concatenated sequences in the handles).
// Arguments:
//      _haploGraph: a gbwtgraph::GBWTGraph which contains the handles in vector< handle_t
//      > haplotypes. haplotypte_handle_vectors: a vector of haplotypes in vector<
//      handle_t > format.
// Returns: a vector of haplotypes of format string (which is the concatenated sequences
// in the handles).
unordered_set<string> SnarlNormalizer::format_handle_haplotypes_to_strings(
    const vector<vector<handle_t>> &haplotype_handle_vectors) {
    unordered_set<string> haplotype_strings;
    for (vector<handle_t> haplotype_handles : haplotype_handle_vectors) {
        string hap;
        for (handle_t &handle : haplotype_handles) {
            hap += _haploGraph.get_sequence(handle);
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
    unordered_set<string> source_to_sink_haplotypes) {
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

    // replace the source and sink chars with X, to force match at source and sink.
    for (auto hap : source_to_sink_haplotypes) {
        hap.replace(0, 1, "X");
        hap.replace(hap.size() - 1, 1, "X");
    }

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

    seqan::resize(rows(align), source_to_sink_haplotypes.size());
    int i = 0;
    for (auto hap : source_to_sink_haplotypes) {
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
        cerr << "ROW_STRING: " << row_string << endl;
        // edit the row so that the proper source and sink chars are added to the
        // haplotype instead of the special characters added to ensure correct alignment
        // of source and sink.
        row_string.replace(0, 1, source_char);
        row_string.replace(row_string.size() - 1, 1, sink_char);
        row_strings.push_back(row_string);
    }

    stringstream ss;
    for (string seq : row_strings) {
        // todo: debug_statement
        cerr << "seq in alignment:" << seq << endl;
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
void SnarlNormalizer::force_maximum_handle_size(MutableHandleGraph &graph,
                                                const size_t &max_size) {
    // forcing each handle in the _graph to have a maximum sequence length of max_size:
    _graph.for_each_handle([&](handle_t handle) {
        // all the positions we want to make in the handle are in offsets.
        vector<size_t> offsets;

        size_t sequence_len = _graph.get_sequence(handle).size();
        int number_of_divisions = floor(sequence_len / max_size);

        // if the handle divides evenly into subhandles of size max_size, we don't need to
        // make the last cut (which would be at the very end of the handle - cutting off
        // no sequence).
        if (sequence_len % max_size == 0) {
            number_of_divisions--;
        }

        // calculate the position of all the divisions we want to make.
        for (int i = 1; i <= number_of_divisions; i++) {
            offsets.push_back(i * max_size);
        }

        // divide the handle into parts.
        _graph.divide_handle(handle, offsets);
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
                                                 const id_t &start_id,
                                                 const id_t &sink_id) {
    // cerr << "extract_subgraph" << endl;
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    unordered_set<id_t> visited;  // to avoid counting the same node twice.
    unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

    // TODO: how to ensure that "to the right" of start_handle is the correct direction?
    // initialize with start_handle (because we move only to the right of start_handle):
    handle_t start_handle = _graph.get_handle(start_id);
    subgraph.add_handle(start_handle);
    visited.insert(graph.get_id(start_handle));

    // look only to the right of start_handle
    _graph.follow_edges(start_handle, false, [&](const handle_t &handle) {
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

        if (graph.get_id(cur_handle) != sink_id) { // don't iterate past end node!
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
//      source_id: the source of the old (to be replaced) snarl in _graph
//      sink_id: the sink of the old (to be replaced) snarl in _graph.
// Return: None.
// TODO: Note: How to ensure that step_handle_t's walk along the snarl in the same
// TODO:     orientation as we expect? i.e. that they don't move backward? I think
// TODO:     we want match_orientation to be = true, but this may cause problems
// TODO:     in some cases given the way we currently construct handles (fixed when we
// TODO:     create snarl-scanning interface).
// TODO:     It may also be that we *don't want match_orientation to be true,
// TODO:     if we're tracking a path that loops backward in the snarl. Hmm... Will think
// about this.
void SnarlNormalizer::integrate_snarl(
    const HandleGraph &to_insert_snarl,
    const vector<pair<step_handle_t, step_handle_t>> embedded_paths, 
    const id_t &source_id, const id_t &sink_id) {
    // cerr << "integrate_snarl" << endl;

    //todo: debug_statement
    cerr << "handles in to_insert_snarl:" << endl;
    to_insert_snarl.for_each_handle([&](const handle_t &handle) {
        cerr << to_insert_snarl.get_id(handle) << " "
             << to_insert_snarl.get_sequence(handle) << " \t";
    });
    cerr << endl;
    // Get old _graph snarl
    SubHandleGraph old_snarl = extract_subgraph(_graph, source_id, sink_id);

    // TODO: debug_statement: Check to make sure that newly made snarl has only one start
    // and end.
    // TODO:     (shouldn't be necessary once we've implemented alignment with
    // leading/trailing special chars.) Identify old and new snarl start and sink
    pair<vector<handle_t>, vector<handle_t>> to_insert_snarl_defining_handles =
        debug_get_sources_and_sinks(to_insert_snarl);

    if (to_insert_snarl_defining_handles.first.size() > 1 ||
        to_insert_snarl_defining_handles.second.size() > 1) {
        cerr << "ERROR: newly made snarl from a snarl starting at " << source_id
             << " has more than one start or end. # of starts: "
             << to_insert_snarl_defining_handles.first.size()
             << " # of ends: " << to_insert_snarl_defining_handles.second.size() << endl;
        return;
    }

    /// Replace start and end handles of old _graph snarl with to_insert_snarl start and
    /// end, and delete rest of old _graph snarl:

    // add to_insert_snarl into _graph without directly attaching the snarl to the _graph
    // (yet).
    vector<handle_t> to_insert_snarl_topo_order =
        algorithms::lazier_topological_order(&to_insert_snarl);

    // Construct a parallel new_snarl_topo_order to identify
    // paralogous nodes between to_insert_snarl and the new snarl inserted in _graph.
    vector<handle_t> new_snarl_topo_order;

    // integrate the handles from to_insert_snarl into the _graph, and keep track of their
    // identities by adding them to new_snarl_topo_order.
    for (handle_t to_insert_snarl_handle : to_insert_snarl_topo_order) {
        // //todo: debug_statement:
        // cerr << " pre-inserted snarl handle: "
        //      << to_insert_snarl.get_id(to_insert_snarl_handle) << " "
        //      << to_insert_snarl.get_sequence(to_insert_snarl_handle) << endl;

        handle_t graph_handle =
            _graph.create_handle(to_insert_snarl.get_sequence(to_insert_snarl_handle));
        new_snarl_topo_order.push_back(graph_handle);
        cerr << "graph handle being inserted into new_snarl_topo_order:" << _graph.get_id(graph_handle) << endl;
    }

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
    id_t temp_snarl_source_id = _graph.get_id(new_snarl_topo_order.front());
    id_t temp_snarl_sink_id = _graph.get_id(new_snarl_topo_order.back());
    cerr << "the temp source id: " << temp_snarl_source_id << endl;
    cerr << "the temp sink id: " << temp_snarl_sink_id << endl;

    // Add the neighbors of the source and sink of the original snarl to the new_snarl's
    // source and sink.
    // source integration:
    _graph.follow_edges(
        _graph.get_handle(source_id), true, [&](const handle_t &prev_handle) {
            _graph.create_edge(prev_handle, _graph.get_handle(temp_snarl_source_id));
        });
    _graph.follow_edges(
        _graph.get_handle(sink_id), false, [&](const handle_t &next_handle) {
            _graph.create_edge(_graph.get_handle(temp_snarl_sink_id), next_handle);
        });

    // For each path of interest, move it onto the new_snarl.
    for (auto path : embedded_paths) {
        // //todo: debug_statement
        // cerr << "the new sink id: " << temp_snarl_sink_id << endl;
        move_path_to_snarl(path, new_snarl_topo_order, temp_snarl_source_id,
                           temp_snarl_sink_id, source_id, sink_id);
    }

    // Destroy the old snarl.
    old_snarl.for_each_handle(

        [&](const handle_t &handle) {
            // //todo: debug_statement these are the handles in old_snarl:
            // cerr << old_snarl.get_id(handle) << old_snarl.get_sequence(handle) << endl;
            _graph.destroy_handle(handle);
        });

    // Replace the source and sink handles with ones that have the original source/sink id
    // (for compatibility with future iterations on neighboring top-level snarls using the
    // same snarl manager. Couldn't replace it before b/c we needed the old handles to
    // move the paths.
    handle_t new_source_handle = _graph.create_handle(
        _graph.get_sequence(_graph.get_handle(temp_snarl_source_id)), source_id);
    handle_t new_sink_handle = _graph.create_handle(
        _graph.get_sequence(new_snarl_topo_order.back()), sink_id);

    // move the source edges:
    // TODO: note the copy/paste. Ask if there's a better way to do this (I totally could
    // in Python!)
    _graph.follow_edges(_graph.get_handle(temp_snarl_source_id), true,
                        [&](const handle_t &prev_handle) {
                            _graph.create_edge(prev_handle, new_source_handle);
                        });
    _graph.follow_edges(_graph.get_handle(temp_snarl_source_id), false,
                        [&](const handle_t &next_handle) {
                            _graph.create_edge(new_source_handle, next_handle);
                        });

    // move the sink edges:
    _graph.follow_edges(_graph.get_handle(temp_snarl_sink_id), true,
                        [&](const handle_t &prev_handle) {
                            _graph.create_edge(prev_handle, new_sink_handle);
                        });
    _graph.follow_edges(_graph.get_handle(temp_snarl_sink_id), false,
                        [&](const handle_t &next_handle) {
                            _graph.create_edge(new_sink_handle, next_handle);
                        });

    // move the paths:
    _graph.for_each_step_on_handle(
        _graph.get_handle(temp_snarl_source_id), [&](step_handle_t step) {
            _graph.rewrite_segment(step, _graph.get_next_step(step),
                                   vector<handle_t>{new_source_handle});
        });
    _graph.for_each_step_on_handle(
        _graph.get_handle(temp_snarl_sink_id), [&](step_handle_t step) {
            _graph.rewrite_segment(step, _graph.get_next_step(step),
                                   vector<handle_t>{new_sink_handle});
        });
    cerr << "the temp source id: " << temp_snarl_source_id << endl;
    cerr << "the temp sink id: " << temp_snarl_sink_id << endl;

    // delete the previously created source and sink:
    for (handle_t handle : {_graph.get_handle(temp_snarl_source_id),
                            _graph.get_handle(temp_snarl_sink_id)}) {
        cerr << "id of handle to delete from tem source/sink: " << _graph.get_id(handle) << endl;
        _graph.destroy_handle(handle);
    }
}

// Moves a path from its original location in the _graph to a new snarl,
//      defined by a vector of interconnected handles.
//      NOTE: the handles in new_snarl_handles may not preserve topological order after
//      being passed to this method, if they were ordered before.
// Arguments: _graph: the _graph containing the old_embedded_path and the handles in
// new_snarl_topo_order
//            old_embedded_path: a pair, where
//                          pair.first is the first step_handle of interest in the
//                          old_embedded_path, and pair.second is the step_handle *after*
//                          the last step_handle of interest in the old_embedded_path (can
//                          be the null step at the end of the path.)
//            new_snarl_topo_order: all the handles in the new snarl, inside the _graph.
// Return: None.
void SnarlNormalizer::move_path_to_snarl(
    const pair<step_handle_t, step_handle_t> &old_embedded_path,
    vector<handle_t> &new_snarl_handles, id_t &new_source_id, id_t &new_sink_id,
    const id_t &old_source_id, const id_t &old_sink_id) {
    // cerr << "\nmove_path_to_snarl" << endl;
    // //TODO: debug_statement:
    // cerr << "path name: "
    //      <<
    //      _graph.get_path_name(_graph.get_path_handle_of_step(old_embedded_path.first))
    //      << endl;
    // cerr << "source: " << new_source_id << " sink: " << new_sink_id << endl;
    // if (_graph.get_path_name(_graph.get_path_handle_of_step(old_embedded_path.first))
    // ==
    //     "chr10") {
    //     cerr << "\t\tstart and end of old embedded path: "
    //          << _graph.get_id(_graph.get_handle_of_step(old_embedded_path.first))
    //          << "end id"
    //          << _graph.get_id(_graph.get_handle_of_step(old_embedded_path.second))
    //          << endl;
    // }
    // cerr << "#### handles in snarl (according to move_path_to_snarl): ####" << endl;
    // for (handle_t handle : new_snarl_handles) {
    //     cerr << "\t" << _graph.get_id(handle) << " " << _graph.get_sequence(handle);
    // }
    // cerr << endl << endl;
    // cerr << "~~~~~ Handles following each handle:" << endl;
    // for (handle_t handle : new_snarl_handles) {
    //     cerr << "neighbors of handle " << _graph.get_id(handle) << " ("
    //          << _graph.get_sequence(handle) << "):" << endl;
    //     _graph.follow_edges(handle, false, [&](const handle_t &next_handle) {
    //         cerr << "\t" << _graph.get_id(next_handle) << " "
    //              << _graph.get_sequence(next_handle) << endl;
    //     });
    // }

    // get the sequence associated with the path
    string path_seq;
    step_handle_t cur_step = old_embedded_path.first;

    // if the old path is touching either/both the source/sink, we want to make sure that
    // the newly moved path also touches those. Otherwise, any paths that extend beyond
    // the source or sink may be cut into pieces when the portion of the path overlapping
    // the snarl is moved to a region inside the snarl.
    bool touching_source =
        (_graph.get_id(_graph.get_handle_of_step(old_embedded_path.first)) ==
         old_source_id);
    bool touching_sink =
        (_graph.get_id(_graph.get_handle_of_step(
             _graph.get_previous_step(old_embedded_path.second))) == old_sink_id);

    // extract the path sequence of the embedded path:
    while (cur_step != old_embedded_path.second) {
        path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
        cur_step = _graph.get_next_step(cur_step);
    }

    // TODO: debug_statement:
    // cerr << "\t\tpath sequence length: " << path_seq.size() << endl;
    // cerr << "path sequence: " << path_seq << endl;

    // for the given path, find every good possible starting handle in the new_snarl
    //      format of pair is < possible_path_handle_vec,
    //      starting_index_in_the_first_handle, current_index_in_path_seq>
    // //todo: debug_statement
    // cerr << "checking handles as start of path-seq" << endl;
    vector<tuple<vector<handle_t>, int, int>> possible_paths;
    for (handle_t handle : new_snarl_handles) {
        string handle_seq = _graph.get_sequence(handle);

        // starting index is where the path would begin in the handle,
        // since it could begin in the middle of the handle.
        vector<int> starting_indices =
            check_handle_as_start_of_path_seq(handle_seq, path_seq);

        // if there is a starting index,
        if (starting_indices.size() != 0) {
            for (int starting_index : starting_indices) {
                if ((handle_seq.size() - starting_index) >= path_seq.size() &&
                    source_and_sink_handles_map_properly(_graph, new_source_id,
                                                         new_sink_id, touching_source,
                                                         touching_sink, handle, handle)) {
                    // if the entire path fits inside the current handle, and if any
                    // paths that touched source and sink in the old snarl would be
                    // touching source and sink in the new snarl, then we've already
                    // found the full mapping location of the path! Move the path, end
                    // the method.
                    vector<handle_t> new_path{handle};
                    _graph.rewrite_segment(old_embedded_path.first,
                                           old_embedded_path.second, new_path);
                    // //todo: debug_statement
                    // cerr << "found a full mapping at " << _graph.get_id(handle)
                    //      << " w/ seq " << _graph.get_sequence(handle) << endl;
                    return;
                } else {
                    // this is a potential starting handle for the path. Add as a
                    // possible_path.
                    vector<handle_t> possible_path_handle_vec{handle};
                    possible_paths.push_back(
                        make_tuple(possible_path_handle_vec, starting_index,
                                   handle_seq.size() - starting_index));
                }
            }
        }
    }

    // //todo: debug_statement:
    // cerr << "done checking handles as start of path seq" << endl;

    // //TODO: debug_statement:
    // cerr << "possible paths so far: " << endl;
    // for (tuple<vector<handle_t>, int, int> path : possible_paths) {
    //     cerr << " possible start: ";
    //     for (handle_t handle : get<0>(path)) {
    //         cerr << _graph.get_id(handle) << " ";
    //     }
    //     cerr << endl;
    // }

    // for every possible path, extend it to determine if it really is the path we're
    // looking for:
    while (!possible_paths.empty()) {
        // take a path off of possible_paths, which will be copied for every iteration
        // through _graph.follow_edges, below:
        tuple<vector<handle_t>, int, int> possible_path_query = possible_paths.back();
        possible_paths.pop_back();

        // //TODO: debug_statement:
        // for (tuple<vector<handle_t>, int, int> path : possible_paths) {
        //     cerr << "*\tpossible path query: ";
        //     for (handle_t handle : get<0>(possible_path_query)) {
        //         cerr << _graph.get_id(handle) << " " << _graph.get_sequence(handle)
        //              << " ";
        //     }
        //     cerr << endl;
        // }

        // extend the path through all right-extending edges to see if any subsequent
        // paths still satisfy the requirements for being a possible_path:
        bool no_path = _graph.follow_edges(
            get<0>(possible_path_query).back(), false, [&](const handle_t &next) {
                // //todo: debug_statement
                // cerr << "cur handle id: "
                //      << _graph.get_id(get<0>(possible_path_query).back()) << endl;

                // cerr << "next handle id and seq: " << _graph.get_id(next) << " "
                //      << _graph.get_sequence(next) << endl;
                // make a copy to be extended for through each possible next handle in
                // follow edges.
                tuple<vector<handle_t>, int, int> possible_path = possible_path_query;

                // extract relevant information to make code more readable.
                string next_seq = _graph.get_sequence(next);
                id_t next_id = _graph.get_id(next);
                int &cur_index_in_path = get<2>(possible_path);
                if (cur_index_in_path <= path_seq.size() &&
                    (find(new_snarl_handles.cbegin(), new_snarl_handles.cend(), next) !=
                     new_snarl_handles.cend())) {
                    // if the next handle would be the ending handle for the path,
                    if (next_seq.size() >= (path_seq.size() - cur_index_in_path)) {
                        // cerr << "next handle would be the ending handle for the path"
                        //      << endl;
                        //     check to see if the sequence in the handle is suitable
                        // for ending the path:
                        int compare_length = path_seq.size() - cur_index_in_path;

                        // //todo: debug_statement
                        // cerr << "about to compare. compare val: "
                        //      << (next_seq.compare(0, compare_length, path_seq,
                        //                           cur_index_in_path, compare_length) ==
                        //                           0)
                        //      << " source_and_sink_handles_map "
                        //      << source_and_sink_handles_map_properly(
                        //             _graph, new_source_id, new_sink_id,
                        //             touching_source, touching_sink,
                        //             get<0>(possible_path).front(), next)
                        //      << endl;
                        // cerr << "arguments of compare: "
                        //      << " " << 0 << " " << compare_length << " " << path_seq
                        //      << " " << cur_index_in_path << " " << compare_length << "
                        //      "
                        //      << endl;
                        if ((next_seq.compare(0, compare_length, path_seq,
                                              cur_index_in_path, compare_length) == 0) &&
                            source_and_sink_handles_map_properly(
                                _graph, new_source_id, new_sink_id, touching_source,
                                touching_sink, get<0>(possible_path).front(), next)) {
                            // todo: debug_statement
                            // cerr << "compared." << endl;

                            // we've found the new path! Move path to the new sequence,
                            // and end the function.

                            if (compare_length < next_seq.size()) {
                                // If the path ends before the end of next_seq, then split
                                // the handle so that the path ends flush with the end of
                                // the first of the two split handles.

                                // divide the handle where the path ends;
                                pair<handle_t, handle_t> divided_next =
                                    _graph.divide_handle(next, compare_length);
                                get<0>(possible_path).push_back(divided_next.first);

                                // Special case if next is the sink or the source, to
                                // preserve the reassignment of source and sink ids in
                                // integrate_snarl.
                                if (next_id == new_sink_id) {
                                    new_sink_id = _graph.get_id(divided_next.second);
                                }

                                // TODO: NOTE: finding the old "next" handle is expensive.
                                // TODO:   Use different container?
                                auto it = find(new_snarl_handles.begin(),
                                               new_snarl_handles.end(), next);

                                // replace the old invalidated handle with one of the new
                                // ones
                                *it = divided_next.first;
                                // stick the other new handle on the end of
                                // new_snarl_handles.
                                new_snarl_handles.push_back(divided_next.second);

                            } else {
                                // otherwise, the end of the path already coincides with
                                // the end of the handle. In that case, just add it to the
                                // path.
                                get<0>(possible_path).push_back(next);
                            }
                            _graph.rewrite_segment(old_embedded_path.first,
                                                   old_embedded_path.second,
                                                   get<0>(possible_path));
                            // //todo: debug_statement:
                            // cerr << "got a full path: ";
                            // for (handle_t handle : get<0>(possible_path)) {
                            //     cerr << _graph.get_id(handle) << " ";
                            // }
                            // cerr << endl;

                            // we've already found the path. No need to keep looking for
                            // more paths.
                            return false;
                        }
                    }
                    // see if the next handle would be the continuation of the path, but
                    // not the end,
                    else {

                        // check to see if the sequence in the handle is suitable for
                        // extending the path:
                        int compare_length = next_seq.size();
                        // //todo: debug_statement
                        // cerr << "compare returned false" << endl;
                        // cerr << "compare in returned false: "
                        //      << " next_seq len " << next_seq.size() << " compare_length
                        //      "
                        //      << compare_length << " path_seq len " << path_seq.size()
                        //      << " cur_index_in_path " << cur_index_in_path << endl;
                        // cerr << "if statement eval: cur_index_in_path <=
                        // next_seq.size() "
                        //      << (cur_index_in_path <= next_seq.size())
                        //      << " next_seq.compare(0, compare_length, path_seq, "
                        //         "cur_index_in_path, compare_length) == 0) "
                        //      << (next_seq.compare(0, compare_length, path_seq,
                        //                           cur_index_in_path, compare_length) ==
                        //                           0)
                        //      << endl;
                        if (next_seq.compare(0, compare_length, path_seq,
                                             cur_index_in_path, compare_length) == 0) {
                            // cerr << "compared in return false" << endl;
                            // extend the path
                            get<0>(possible_path).push_back(next);

                            // update the current index in path_seq.
                            get<2>(possible_path) += next_seq.size();

                            // place back into possible_paths
                            possible_paths.push_back(possible_path);
                            // cerr << "extending the path!" << endl;
                        }
                    }
                }
                // continue to iterate through follow_edges.
                return true;
            });

        // //todo: debug_statement:
        // if
        // (graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first))
        // ==
        //     "_alt_19f9bc9ad2826f58f113965edf36bb93740df46d_0") {
        //     cerr << "mystery node 4214930: "
        //          << _graph.get_sequence(graph.get_handle(4214930)) << endl;
        // }

        // if we've found a complete path in the above follow_edges, then we've
        // already moved the path, and we're done.
        if (!no_path) {
            return;
        }
    }
    // //todo: figure out how to do some better error message instead of cerr.
    // if we failed to find a path, show an error message.
    cerr << "##########################\nWarning! Didn't find a corresponding path of "
            "name "
         << _graph.get_path_name(_graph.get_path_handle_of_step(old_embedded_path.first))
         << " from the old snarl at " << old_source_id
         << " in the newly aligned snarl. This snarl WILL be "
            "normalized, resulting in a probably incorrectly-constructed snarl."
            "\n##########################"
         << endl
         << endl;
    // throw _graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first));
    // assert(true && "Warning! Didn't find a corresponding path of name " +
    //         _graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first))
    //         + " from the old snarl in the newly aligned snarl.");
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
    _graph.for_each_handle([&](const handle_t &handle) {
        bool is_source = true, is_sink = true;
        _graph.follow_edges(handle, true, [&](const handle_t &prev) {
            is_source = false;
            return false;
        });
        _graph.follow_edges(handle, false, [&](const handle_t &next) {
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

}
}