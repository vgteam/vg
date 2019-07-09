/**
 * \file
 * surjector.cpp: implements a class that surjects alignments onto paths
 */

#include "surjector.hpp"

//#define debug_anchored_surject
//#define debug_validate_anchored_multipath_alignment

namespace vg {

using namespace std;
    
    Surjector::Surjector(const XG* xg_index) : xindex(xg_index) {
        if (!xindex) {
            cerr << "error:[Surjector] Failed to provide an XG index to the Surjector" << endl;
        }
    }
    
    Alignment Surjector::surject(const Alignment& source, const set<string>& path_names,
                                 bool allow_negative_scores) const {
    
        // Allocate the annotation info
        string path_name_out;
        int64_t path_pos_out;
        bool path_rev_out;
        
        // Do the surjection
        Alignment surjected = surject(source, path_names, path_name_out, path_pos_out, path_rev_out, allow_negative_scores);
        
        // Pack all the info into the refpos field
        surjected.clear_refpos();
        auto* pos = surjected.add_refpos();
        pos->set_name(path_name_out);
        pos->set_offset(path_pos_out);
        pos->set_is_reverse(path_rev_out);
    
        return surjected;
    }
    
    Alignment Surjector::surject(const Alignment& source, const set<string>& path_names,
                                 string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                 bool allow_negative_scores) const {
        
#ifdef debug_anchored_surject
        cerr << "surjecting alignment: " << pb2json(source) << " onto paths ";
        for (const string& path_name : path_names) {
            cerr << path_name << " ";
        }
        cerr << endl;
#endif
        
        if (source.path().mapping_size() != 0) {
            // The read is mapped. Check the input alignment for basic
            // consistency. If the sequence and the graph path don't agree
            // about the read length, something is very wrong with the input.
            size_t source_to_length = path_to_length(source.path());
            if (source.sequence().size() != source_to_length) {
                cerr << "error[Surjector::surject]: read " << source.name() << " has "
                    << source.sequence().size() << " sequence bases but an input alignment that aligns "
                    << source_to_length << " bases instead. This is invalid and uninterpretable; check your mapper." << endl;
                exit(1);
            }
        }

        // translate the path names into ranks for the XG and path hanles
        unordered_set<path_handle_t> surjection_path_handles;
        for (const string& path_name : path_names) {
            surjection_path_handles.insert(xindex->get_path_handle(path_name));
        }
        
        // make an overlay that will memoize the results of some expensive XG operations
        MemoizingGraph memoizing_graph(xindex);
        
        // get the chunks of the aligned path that overlap the ref path
        auto path_overlapping_anchors = extract_overlapping_paths(&memoizing_graph, source, surjection_path_handles);
        
#ifdef debug_anchored_surject
        cerr << "got path overlapping segments" << endl;
        for (const auto& path_record : path_overlapping_anchors) {
            cerr << "path " << xindex->get_path_name(path_record.first) << endl;
            for (auto& anchor : path_record.second) {
                cerr << "\t read[" << (anchor.first.first - source.sequence().begin()) << ":" << (anchor.first.second - source.sequence().begin()) << "] : ";
                for (auto iter = anchor.first.first; iter != anchor.first.second; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
                cerr << "\t" << pb2json(anchor.second) << endl;
            }
        }
#endif
        
        // the surjected alignment for each path we overlapped
        unordered_map<path_handle_t, Alignment> path_surjections;
        for (pair<const path_handle_t, vector<path_chunk_t>>& path_record : path_overlapping_anchors) {
#ifdef debug_anchored_surject
            cerr << "found overlaps on path " << xindex->get_path_name(path_record.first) << ", performing surjection" << endl;
#endif
            
            // find the interval of the ref path we need to consider
            pair<size_t, size_t> ref_path_interval = compute_path_interval(&memoizing_graph, source, path_record.first,
                                                                           path_record.second);
            
#ifdef debug_anchored_surject
            cerr << "final path interval is " << ref_path_interval.first << ":" << ref_path_interval.second << endl;
#endif
            
            // get the path graph corresponding to this interval
            sglib::HashGraph path_graph;
            unordered_map<id_t, pair<id_t, bool>> path_trans = extract_linearized_path_graph(&memoizing_graph, &path_graph, path_record.first,
                                                                                             ref_path_interval.first, ref_path_interval.second);
            
            // split it into a forward and reverse strand
            sglib::HashGraph split_path_graph;
            unordered_map<id_t, pair<id_t, bool>> split_trans = algorithms::split_strands(&path_graph, &split_path_graph);
            
            auto node_trans = overlay_node_translations(split_trans, path_trans);
            
#ifdef debug_anchored_surject
            cerr << "made split, linearized path graph" << endl;
#endif
            
            // compute the connectivity between the path chunks
            MultipathAlignmentGraph mp_aln_graph(split_path_graph, path_record.second, source, node_trans);
            
            // we don't overlap this reference path at all or we filtered out all of the path chunks, so just make a sentinel
            if (mp_aln_graph.empty()) {
                path_surjections[path_record.first] = make_null_alignment(source);
                continue;
            }
            
            // TODO: is this necessary in a linear graph?
            vector<size_t> topological_order;
            mp_aln_graph.topological_sort(topological_order);
            mp_aln_graph.remove_transitive_edges(topological_order);
            
            // align the intervening segments and store the result in a multipath alignment
            MultipathAlignment mp_aln;
            mp_aln_graph.align(source, split_path_graph, get_aligner(), false, 1, false, 1, mp_aln, allow_negative_scores);
            topologically_order_subpaths(mp_aln);
            
            for (size_t i = 0; i < mp_aln.subpath_size(); i++) {
                // translate back into the original ID space
                translate_oriented_node_ids(*mp_aln.mutable_subpath(i)->mutable_path(), node_trans);
            }
            
            // identify the source subpaths (necessary for subpath-global optimal alignment algorithm)
            identify_start_subpaths(mp_aln);
            
#ifdef debug_anchored_surject
            cerr << "made multipath alignment " << pb2json(mp_aln) << endl;
#endif
            
#ifdef debug_validate_anchored_multipath_alignment
            if (!validate_multipath_alignment(mp_aln, *xindex)) {
                cerr << "WARNING: multipath alignment for surjection of " << source.name() << " failed to validate" << endl;
            }
#endif
            
            // concatenate the subpaths either locally or globally, depending on whether we're
            // allowing negative scores
            Alignment& surjected = path_surjections[path_record.first];
            optimal_alignment(mp_aln, surjected, allow_negative_scores);
            
            // transfer applicable metadata
            surjected.set_mapping_quality(source.mapping_quality());
            if (source.has_fragment_next()) {
                *surjected.mutable_fragment_next() = source.fragment_next();
            }
            if (source.has_fragment_prev()) {
                *surjected.mutable_fragment_prev() = source.fragment_prev();
            }
            
#ifdef debug_anchored_surject
            cerr << "concatenated and translated alignment " << pb2json(surjected) << endl;
#endif
            
        }
        
        // in case we didn't overlap any paths, add a sentinel so the following code still executes correctly
        if (path_surjections.empty()) {
            path_surjections[handlegraph::as_path_handle(0)] = make_null_alignment(source);
        }
        
        // choose which path surjection was best
        path_handle_t best_path_handle;
        int32_t score = numeric_limits<int32_t>::min();
        for (const auto& surjection : path_surjections) {
            if (surjection.second.score() >= score) {
                score = surjection.second.score();
                best_path_handle = surjection.first;
            }
        }
        Alignment& best_surjection = path_surjections[best_path_handle];
        
        // find the position along the path
        set_path_position(&memoizing_graph, best_surjection, best_path_handle,
                          path_name_out, path_pos_out, path_rev_out);
        
#ifdef debug_anchored_surject
        cerr << "chose path " << path_name_out << " at position " << path_pos_out << (path_rev_out ? "-" : "+") << endl;
#endif
        return move(best_surjection);
    }
    
    unordered_map<path_handle_t, vector<Surjector::path_chunk_t>>
    Surjector::extract_overlapping_paths(const PathPositionHandleGraph* graph, const Alignment& source,
                                         const unordered_set<path_handle_t>& surjection_paths) const {
        
        unordered_map<path_handle_t, vector<path_chunk_t>> to_return;
        
        const Path& path = source.path();
        
        // for each path that we're extending, the previous step and the strand we were at on it
        // mapped to the index of that path chunk in the path's vector
        unordered_map<pair<step_handle_t, bool>, size_t> extending_steps;
        int64_t through_to_length = 0;
        
        for (size_t i = 0; i < path.mapping_size(); i++) {
            
            int64_t before_to_length = through_to_length;
            through_to_length += mapping_to_length(path.mapping(i));
            
            const Position& pos = path.mapping(i).position();
            handle_t handle = graph->get_handle(pos.node_id(), pos.is_reverse());
            
#ifdef debug_anchored_surject
            cerr << "looking for paths on mapping " << i << " at position " << make_pos_t(pos) << endl;
#endif
            
            unordered_map<pair<step_handle_t, bool>, size_t> next_extending_steps;
            
            for (const step_handle_t& step : graph->steps_of_handle(handle)) {
                
#ifdef debug_anchored_surject
                cerr << "found a step on " << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
#endif
                
                path_handle_t path_handle = graph->get_path_handle_of_step(step);
                if (!surjection_paths.count(path_handle)) {
                    // we are not surjecting onto this path
#ifdef debug_anchored_surject
                    cerr << "not surjecting to this path, skipping" << endl;
#endif
                    continue;
                }
                
                // the chunks of the alignments along this path
                vector<path_chunk_t>& path_chunks = to_return[path_handle];
                
                bool path_strand = graph->get_is_reverse(handle) != graph->get_is_reverse(graph->get_handle_of_step(step));
                
                step_handle_t prev_step = path_strand ? graph->get_next_step(step) : graph->get_previous_step(step);
                
#ifdef debug_anchored_surject
                cerr << "path strand is " << (path_strand ? "rev" : "fwd") << ", prev step is ";
                if (prev_step == graph->path_end(path_handle)) {
                    cerr << " path end";
                }
                else if (prev_step == graph->path_front_end(path_handle)) {
                    cerr << " path front end";
                }
                else {
                    cerr << graph->get_id(graph->get_handle_of_step(prev_step)) << (graph->get_is_reverse(graph->get_handle_of_step(prev_step)) ? "-" : "+");
                }
                cerr << endl;
                cerr << "possible extensions from: " << endl;
                for (const auto& record : extending_steps) {
                    cerr << "\t" << graph->get_id(graph->get_handle_of_step(record.first.first)) << (graph->get_is_reverse(graph->get_handle_of_step(record.first.first)) ? "-" : "+") << " on " << graph->get_path_name(graph->get_path_handle_of_step(record.first.first)) << " " << (record.first.second ? "rev" : "fwd") << endl;
                }
#endif
                
                if (extending_steps.count(make_pair(prev_step, path_strand))) {
                    // we are extending from the previous step, so we continue with the extension
                    
                    auto& path_chunks = to_return[graph->get_path_handle_of_step(step)];
                    path_chunk_t& chunk = path_chunks[extending_steps[make_pair(prev_step, path_strand)]];
                    
                    // move the end of the sequence out
                    chunk.first.second = source.sequence().begin() + through_to_length;
                    Mapping* mapping = chunk.second.add_mapping();
                    // add this mapping
                    *mapping = path.mapping(i);
                    mapping->set_rank(chunk.second.mapping(chunk.second.mapping_size() - 2).rank() + 1);
                    
                    // in the next iteration, this step should point into the chunk it just extended
                    next_extending_steps[make_pair(step, path_strand)] = extending_steps[make_pair(prev_step, path_strand)];
                }
                else {
                    // this step does not extend a previous step, so we start a new chunk
                    auto& path_chunks = to_return[graph->get_path_handle_of_step(step)];
                    path_chunks.emplace_back();
                    path_chunk_t& chunk = path_chunks.back();
                    
                    // init the new chunk with the sequence interval
                    chunk.first.first = source.sequence().begin() + before_to_length;
                    chunk.first.second = source.sequence().begin() + through_to_length;
                    
                    // and with the first mapping
                    Mapping* mapping = chunk.second.add_mapping();
                    *mapping = path.mapping(i);
                    mapping->set_rank(1);
                    
                    // keep track of where this chunk is in the vector and which step it came from
                    // for the next iteration
                    next_extending_steps[make_pair(step, path_strand)] = path_chunks.size() - 1;
                }
            }
            
            // we've finished extending the steps from the previous mapping, so we replace them
            // with the steps we found in this iteration that we want to extend on the next one
            extending_steps = next_extending_steps;
        }
        
        return to_return;
    }
    
    pair<size_t, size_t>
    Surjector::compute_path_interval(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                                     const vector<path_chunk_t>& path_chunks) const {
        
        pair<size_t, size_t> interval(numeric_limits<size_t>::max(), numeric_limits<size_t>::min());
        
        for (const auto& path_chunk : path_chunks) {
            
            size_t path_length = graph->get_path_length(path_handle);
            
            size_t left_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.first)
                                    + (path_chunk.first.first - source.sequence().begin()));
            
            size_t right_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.second)
                                     + (source.sequence().end() - path_chunk.first.second));
            
            const Position& first_pos = path_chunk.second.mapping(0).position();
            handle_t first_handle = graph->get_handle(first_pos.node_id(), first_pos.is_reverse());
            for (const step_handle_t& step : graph->steps_of_handle(first_handle)) {
                
                if (graph->get_path_handle_of_step(step) != path_handle) {
                    // this step isn't on the path we're considering
                    continue;
                }
                
                if (first_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) {
                    size_t path_offset = graph->get_position_of_step(step) + graph->get_length(first_handle) - first_pos.offset();
                    interval.second = max(interval.second, min(path_offset + left_overhang, path_length));
                }
                else {
                    size_t path_offset = graph->get_position_of_step(step) + first_pos.offset();
                    if (left_overhang > path_offset) {
                        // avoid underflow
                        interval.first = 0;
                    }
                    else {
                        interval.first = min(interval.first, path_offset - left_overhang);
                    }
                }
            }
            
            const Mapping& final_mapping = path_chunk.second.mapping(path_chunk.second.mapping_size() - 1);
            const Position& final_pos = final_mapping.position();
            handle_t final_handle = graph->get_handle(final_pos.node_id(), final_pos.is_reverse());
            for (const step_handle_t& step : graph->steps_of_handle(final_handle)) {
                
                if (graph->get_path_handle_of_step(step) != path_handle) {
                    // this step isn't on the path we're considering
                    continue;
                }
                
                if (final_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) {
                    size_t path_offset = graph->get_position_of_step(step) + graph->get_length(final_handle) - final_pos.offset() - mapping_from_length(final_mapping);
                    if (right_overhang > path_offset) {
                        // avoid underflow
                        interval.first = 0;
                    }
                    else {
                        interval.first = min(interval.first, path_offset - right_overhang);
                    }
                }
                else {
                    size_t path_offset = graph->get_position_of_step(step) + first_pos.offset() + mapping_to_length(final_mapping);
                    interval.second = max(interval.second, min(path_offset + right_overhang, path_length));
                }
            }
        }
        
        return interval;
    }
    
    unordered_map<id_t, pair<id_t, bool>>
    Surjector::extract_linearized_path_graph(const PathPositionHandleGraph* graph, MutableHandleGraph* into,
                                             path_handle_t path_handle, size_t first, size_t last) const {
        
        // TODO: we need better semantics than an unsigned interval for surjecting to circular paths
        
#ifdef debug_anchored_surject
        cerr << "extracting path graph for position interval " << first << ":" << last << " in path of length " << graph->get_path_length(path_handle) << endl;
#endif
        
        unordered_map<id_t, pair<id_t, bool>> node_trans;
        
        step_handle_t begin = graph->get_step_at_position(path_handle, first);
        step_handle_t end = graph->get_step_at_position(path_handle, last);
        
        if (graph->get_position_of_step(end) < last && end != graph->path_end(path_handle)) {
            // we actually want part of this step too, so we use the next one as the end iterator
            end = graph->get_next_step(end);
        }
        
        handle_t prev_node;
        for (step_handle_t step = begin; step != end; step = graph->get_next_step(step)) {
            // copy the node with the local orientation now forward
            handle_t copying = graph->get_handle_of_step(step);
            handle_t node_here = into->create_handle(graph->get_sequence(copying));
            
            if (step != begin) {
                // add an edge from the previous node
                into->create_edge(prev_node, node_here);
            }
            
            // record the translation
            node_trans[into->get_id(node_here)] = pair<id_t, bool>(graph->get_id(copying),
                                                                   graph->get_is_reverse(copying));
            
            prev_node = node_here;
        }
        
        return node_trans;
    }
    
    void Surjector::set_path_position(const PathPositionHandleGraph* graph, const Alignment& surjected,
                                      path_handle_t best_path_handle,
                                      string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const {

        const Path& path = surjected.path();
        
        if (path.mapping_size() == 0){
            // hack: we want a 0 position once we convert this to 1-based indexes
            path_pos_out = -1;
            path_rev_out = false;
            path_name_out = "";
            return;
        }
        
#ifdef debug_anchored_surject
        cerr << "choosing a path position for surjected alignment on path " << graph->get_path_name(best_path_handle) << endl;
#endif
        
        const Position& start_pos = path.mapping(0).position();
        handle_t start_handle = graph->get_handle(start_pos.node_id(), start_pos.is_reverse());
        // TODO: depending on path coverage, it could be inefficient to iterate over all steps on the handle
        for (const step_handle_t& step : graph->steps_of_handle(start_handle)) {
#ifdef debug_anchored_surject
            cerr << "found step on " << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
#endif
            if (graph->get_path_handle_of_step(step) != best_path_handle) {
                // this is not the path we surjected onto
#ifdef debug_anchored_surject
                cerr << "wrong path, skipping" << endl;
#endif
                continue;
            }
            
            bool strand = graph->get_is_reverse(graph->get_handle_of_step(step)) != start_pos.is_reverse();
            
            // does the alignment follow the path here?
            bool match = true;
            step_handle_t path_step = strand ? graph->get_previous_step(step) : graph->get_next_step(step);
            for (size_t i = 1; i < path.mapping_size(); ++i) {
                
                if (path_step == graph->path_end(best_path_handle) ||
                    path_step == graph->path_front_end(best_path_handle)) {
                    // we've gone off the end of the path
#ifdef debug_anchored_surject
                    cerr << "encountered end of path at mapping index " << i << ", skipping" << endl;
#endif
                    match = false;
                    break;
                }
                
                const Position& pos = path.mapping(i).position();
#ifdef debug_anchored_surject
                cerr << "at mapping pos " << make_pos_t(pos) << " and path pos " << graph->get_id(graph->get_handle_of_step(path_step)) << (graph->get_is_reverse(graph->get_handle_of_step(path_step)) ? "-" : "+") << endl;
#endif
                
                if (pos.node_id() != graph->get_id(graph->get_handle_of_step(path_step)) ||
                    (pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(path_step)) != strand)) {
                    // the path here doesn't match the surjected path
#ifdef debug_anchored_surject
                    cerr << "encountered mismatch at mapping index " << i << ", skipping" << endl;
#endif
                    match = false;
                    break;
                }
                
                path_step = strand ? graph->get_previous_step(path_step) : graph->get_next_step(path_step);
            }
            
            
            if (match) {
                // we found a match, so set the position values and return
                path_name_out = graph->get_path_name(best_path_handle);
                path_rev_out = strand;
                size_t start_offset = graph->get_position_of_step(step);
                if (strand) {
                    path_pos_out = start_offset + graph->get_length(start_handle) - start_pos.offset() - path_from_length(path);
                }
                else {
                    path_pos_out = start_offset + start_pos.offset();
                }
                return;
            }
        }
        
        // we ran through all of the occurrences without finding a full match...
        
        cerr << "error:[Surjector] could not identify path position of surjected alignment " << surjected.name() << endl;
        exit(1);
    }
    
    Alignment Surjector::make_null_alignment(const Alignment& source) {
        Alignment null;
        null.set_name(source.name());
        null.set_sequence(source.sequence());
        null.set_quality(source.quality());
        null.set_read_group(source.read_group());
        null.set_sample_name(source.sample_name());
        null.set_is_secondary(source.is_secondary());
        if (source.has_fragment_next()) {
            *null.mutable_fragment_next() = source.fragment_next();
        }
        if (source.has_fragment_prev()) {
            *null.mutable_fragment_prev() = source.fragment_prev();
        }
        return null;
    }
}


