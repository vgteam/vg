/**
 * \file
 * surjector.cpp: implements a class that surjects alignments onto paths
 */

#include "surjector.hpp"

//#define debug_surject
//#define debug_anchored_surject
//#define debug_validate_anchored_multipath_alignment

namespace vg {

using namespace std;
    
    Surjector::Surjector(xg::XG* xg_index) : BaseMapper(xg_index, nullptr, nullptr) {
        if (!xindex) {
            cerr << "error:[Surjector] Failed to provide an XG index to the Surjector" << endl;
        }
    }
    
    Surjector::~Surjector() {
        
    }
    
    Alignment Surjector::surject(const Alignment& source, const set<string>& path_names,
                                 string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                 bool allow_negative_scores) {
        
#ifdef debug_anchored_surject
        cerr << "surjecting alignment: " << pb2json(source) << " onto paths ";
        for (const string& path_name : path_names) {
            cerr << path_name << " ";
        }
        cerr << endl;
#endif
        
        // translate the path names into ranks for the XG
        unordered_map<size_t, string> path_rank_to_name;
        for (const string& path_name : path_names) {
            path_rank_to_name[xindex->path_rank(path_name)] = path_name;
        }
        
        
        // memos for expensive succinct operations that may be repeated
        unordered_map<int64_t, vector<size_t>> paths_of_node_memo;
        unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>> oriented_occurrences_memo;
        
        // get the chunks of the aligned path that overlap the ref path
        auto path_overlapping_anchors = extract_overlapping_paths(source, path_rank_to_name, &paths_of_node_memo, &oriented_occurrences_memo);
        
#ifdef debug_anchored_surject
        cerr << "got path overlapping segments" << endl;
        for (const auto& path_record : path_overlapping_anchors) {
            cerr << "path rank " << path_record.first << endl;
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
        unordered_map<size_t, Alignment> path_surjections;
        for (pair<const size_t, vector<path_chunk_t>>& path_record : path_overlapping_anchors) {
#ifdef debug_anchored_surject
            cerr << "found overlaps on path " << path_record.first << ", performing surjection" << endl;
#endif
            
            const xg::XGPath& xpath = xindex->get_path(path_rank_to_name[path_record.first]);
            
            // find the interval of the ref path we need to consider
            pair<size_t, size_t> ref_path_interval = compute_path_interval(source, path_record.first, xpath, path_record.second,
                                                                           &oriented_occurrences_memo);
            
#ifdef debug_anchored_surject
            cerr << "final path interval is " << ref_path_interval.first << ":" << ref_path_interval.second << endl;
#endif
            
            // get the path graph corresponding to this interval
            unordered_map<id_t, pair<id_t, bool>> path_trans;
            VG path_graph = extract_linearized_path_graph(ref_path_interval.first, ref_path_interval.second, xpath, path_trans);
            
            // split it into a forward and reverse strand
            VG split_path_graph;
            unordered_map<id_t, pair<id_t, bool>> split_trans = algorithms::split_strands(&path_graph, &split_path_graph);
            
            algorithms::lazier_sort(&split_path_graph);
            
            auto node_trans = split_path_graph.overlay_node_translations(split_trans, path_trans);
            
#ifdef debug_anchored_surject
            cerr << "made split, linearized path graph " << pb2json(split_path_graph.graph) << endl;
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
            mp_aln_graph.align(source, split_path_graph, get_aligner(), false, 1, false, 1, mp_aln);
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
            path_surjections[0] = make_null_alignment(source);
        }
        
        // choose which path surjection was best
        size_t best_path_rank;
        int32_t score = numeric_limits<int32_t>::min();
        for (const auto& surjection : path_surjections) {
            if (surjection.second.score() >= score) {
                score = surjection.second.score();
                best_path_rank = surjection.first;
            }
        }
        
        // which path was it?
        path_name_out = path_rank_to_name[best_path_rank];
        
        Alignment& best_surjection = path_surjections[best_path_rank];
        
        // find the position along the path
        const xg::XGPath& best_xpath = xindex->get_path(path_name_out);
        set_path_position(best_surjection, best_path_rank, best_xpath, path_name_out, path_pos_out, path_rev_out, &oriented_occurrences_memo);
        
#ifdef debug_anchored_surject
        cerr << "chose path " << path_name_out << " at position " << path_pos_out << (path_rev_out ? "-" : "+") << endl;
#endif
        return move(best_surjection);
    }
    
    unordered_map<size_t, vector<Surjector::path_chunk_t>>
    Surjector::extract_overlapping_paths(const Alignment& source, const unordered_map<size_t, string>& path_rank_to_name,
                                         unordered_map<int64_t, vector<size_t>>* paths_of_node_memo,
                                         unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo) {
        
        
        unordered_map<size_t, vector<path_chunk_t>> to_return;
        
        const Path& path = source.path();
        
        // for each path rank that we're extending, the offset and relative orientation of the previous node
        unordered_map<size_t, unordered_set<pair<size_t, bool>>> offset_and_orientations_on_paths;
        int64_t through_to_length = 0;
        
        for (size_t i = 0; i < path.mapping_size(); i++) {
            
            int64_t before_to_length = through_to_length;
            through_to_length += mapping_to_length(path.mapping(i));
            
            const Position& pos = path.mapping(i).position();
            vector<size_t> paths_of_node = xindex->memoized_paths_of_node(pos.node_id(), paths_of_node_memo);
            
            unordered_set<size_t> paths_here;
            for (size_t path_rank : paths_of_node) {
                if (path_rank_to_name.count(path_rank)) {
                    paths_here.insert(path_rank);
                }
            }
            
            if (paths_here.empty()) {
                // none of the paths that this node is on are in the list we're surjecting onto
                offset_and_orientations_on_paths.clear();
            }
            else {
                // we're on at least one path that we're surjecting onto
                
                // for each path
                for (size_t path_rank : paths_here) {
                    
                    // we'll need to know where this node occurs on the path
                    auto occurrences = xindex->memoized_oriented_occurrences_on_path(pos.node_id(), path_rank,
                                                                                     oriented_occurrences_memo);
                    
                    // the chunks of the alignments along this path
                    vector<path_chunk_t>& path_chunks = to_return[path_rank];
                    
                    // get the location(s) we were extending along the path in the previous iteration (if any)
                    auto& offset_and_orientations_on_path = offset_and_orientations_on_paths[path_rank];
                    
                    if (offset_and_orientations_on_paths.count(path_rank)) {
                        // we were on the path before too, so we might be extending an existing chunk
                        
                        // do any of these locations match up with where we are now?
                        unordered_set<pair<size_t, bool>> next_offsets_and_orientations;
                        for (pair<size_t, bool>& occurrence : occurrences) {
                            bool on_reverse = (occurrence.second != pos.is_reverse());
                            size_t prev_offset = on_reverse ? occurrence.first + 1 : occurrence.first - 1;
                            
                            if (offset_and_orientations_on_path.count(make_pair(prev_offset, on_reverse))) {
                                next_offsets_and_orientations.emplace(occurrence.first, on_reverse);
                            }
                        }
                        
                        if (next_offsets_and_orientations.empty()) {
                            // we're still on the path, but we didn't follow the next edge in the path so this
                            // is actually the start of a new chunk not a continuation of the current chunk
                            
                            // initialize a path chunk
                            path_chunks.emplace_back();
                            path_chunks.back().first.first = source.sequence().begin() + before_to_length;
                            
                            // mark the locations of this node along the path
                            offset_and_orientations_on_path.clear();
                            for (pair<size_t, bool>& occurrence : occurrences) {
                                offset_and_orientations_on_path.emplace(occurrence.first, occurrence.second != pos.is_reverse());
                            }
                        }
                        else {
                            // keep track of where the next positions should come from
                            offset_and_orientations_on_path = move(next_offsets_and_orientations);
                        }
                    }
                    else {
                        // this is the start of a new chunk, initialize it with a mpapping
                        path_chunks.emplace_back();
                        path_chunks.back().first.first = source.sequence().begin() + before_to_length;
                        
                        // mark the locations of this node along the path
                        offset_and_orientations_on_path.clear();
                        for (pair<size_t, bool>& occurrence : occurrences) {
                            offset_and_orientations_on_path.emplace(occurrence.first, occurrence.second != pos.is_reverse());
                        }
                    }
                    
                    // extend the path chunk by this mapping
                    path_chunks.back().first.second = source.sequence().begin() + through_to_length;
                    *path_chunks.back().second.add_mapping() = path.mapping(i);
                }
                
                // if we've left any paths, we need to remove them from the locations index
                vector<size_t> to_erase;
                for (const auto& path_record : offset_and_orientations_on_paths) {
                    if (!paths_here.count(path_record.first)) {
                        to_erase.push_back(path_record.first);
                    }
                }
                for (size_t path_rank : to_erase) {
                    offset_and_orientations_on_paths.erase(path_rank);
                }
            }
        }
        
        return to_return;
    }
    
    pair<size_t, size_t>
    Surjector::compute_path_interval(const Alignment& source, size_t path_rank, const xg::XGPath& xpath, const vector<path_chunk_t>& path_chunks,
                                     unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo) {
        
        pair<size_t, size_t> interval(numeric_limits<size_t>::max(), numeric_limits<size_t>::min());
        
        for (const auto& path_chunk : path_chunks) {
            
            string::const_iterator read_pos = path_chunk.first.first;
            
            // TODO: do I need to do this at every mapping? it might be enough to just look at the first and last
            for (size_t i = 0; i < path_chunk.second.mapping_size(); i++) {
                const Position& pos = path_chunk.second.mapping(i).position();
                
                // the distance the read could align to the left of this mapping (oriented by the read)
                int64_t left_overhang = get_aligner()->longest_detectable_gap(source, read_pos) + (read_pos - source.sequence().begin());
                
                read_pos += mapping_to_length(path_chunk.second.mapping(i));
                
                // the distance the read could align to the right of this mapping (oriented by the read)
                int64_t right_overhang = get_aligner()->longest_detectable_gap(source, read_pos) + (source.sequence().end() - read_pos);
                
                auto oriented_occurrences = xindex->memoized_oriented_occurrences_on_path(pos.node_id(), path_rank,
                                                                                          oriented_occurrences_memo);
                
                // the length forward along the path that the end of the mapping is
                int64_t mapping_length = mapping_from_length(path_chunk.second.mapping(i));
                
                for (const pair<size_t, bool>& occurrence : oriented_occurrences) {
                    if (occurrence.second == pos.is_reverse()) {
                        int64_t path_offset = xpath.positions[occurrence.first];
                        
                        int64_t left_boundary = max<int64_t>(0, path_offset + pos.offset() - left_overhang);
                        interval.first = min<size_t>(interval.first, left_boundary);
                        
                        int64_t right_boundary = min<int64_t>(path_offset + pos.offset() + mapping_length + right_overhang, xpath.offsets.size() - 1);
                        interval.second = max<size_t>(interval.second, right_boundary);
                        
#ifdef debug_anchored_surject
                        cerr << "path chunk " << pb2json(path_chunk.second) << " can be aligned to forward strand in interval " << left_boundary << ":" << right_boundary << endl;
#endif
                    }
                    else {
                        int64_t path_offset = occurrence.first + 1 < xpath.positions.size() ? xpath.positions[occurrence.first + 1] : xpath.offsets.size();
                        
                        int64_t left_boundary = max<int64_t>(0, path_offset - pos.offset() - mapping_length - right_overhang);
                        interval.first = min<size_t>(interval.first, left_boundary);
                        
                        int64_t right_boundary = min<int64_t>(path_offset - pos.offset() + left_overhang, xpath.offsets.size() - 1);
                        interval.second = max<size_t>(interval.second, right_boundary);
                        
#ifdef debug_anchored_surject
                        cerr << "path chunk " << pb2json(path_chunk.second) << " can be aligned to reverse strand in interval " << left_boundary << ":" << right_boundary << endl;
#endif
                    }
                }
            }
        }
        
        return interval;
    }
    
    VG Surjector::extract_linearized_path_graph(size_t first, size_t last, const xg::XGPath& xpath,
                                                unordered_map<id_t, pair<id_t, bool>>& node_trans) {
        
#ifdef debug_anchored_surject
        cerr << "extracting path graph for position interval " << first << ":" << last << " in path of length " << xpath.positions[xpath.positions.size() - 1] + xindex->node_length(xpath.node(xpath.ids.size() - 1)) << endl;
#endif
        
        VG path_graph;
        
        size_t begin = xpath.offset_at_position(first);
        size_t end = min<size_t>(xpath.positions.size(), xpath.offset_at_position(last) + 1);
        
        Node* prev_node = nullptr;
        for (size_t i = begin; i < end; i++) {
            
            id_t node_id = xpath.node(i);
            string seq = xindex->node_sequence(node_id);
            bool rev = xpath.directions[i];
            
            Node* node;
            if (rev) {
                node = path_graph.create_node(reverse_complement(seq));
            }
            else {
                node = path_graph.create_node(seq);
            }
            
            if (prev_node) {
                path_graph.create_edge(prev_node, node);
            }
            prev_node = node;
            
            node_trans[node->id()] = make_pair(node_id, rev);
        }
        
        return path_graph;
    }
    
    void Surjector::set_path_position(const Alignment& surjected, size_t best_path_rank, const xg::XGPath& xpath,
                                      string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                      unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo) {
        
        const Path& path = surjected.path();
        
        if (path.mapping_size() == 0){
            // hack: we want a 0 position once we convert this to 1-based indexes
            path_pos_out = -1;
            path_rev_out = false;
            path_name_out = "";
            return;
        }
        
        const Position& start_pos = path.mapping(0).position();
        
        vector<pair<size_t, bool>> oriented_occurrences = xindex->memoized_oriented_occurrences_on_path(start_pos.node_id(), best_path_rank,
                                                                                                        oriented_occurrences_memo);

        for (pair<size_t, bool>& occurrence : oriented_occurrences) {
            if (occurrence.second == start_pos.is_reverse()) {
                // the first node in this alignment occurs on the forward strand of the path
                
                if (occurrence.first + path.mapping_size() > xpath.ids.size()) {
                    // but it doesn't fit on the path
                    continue;
                }
                
                // does the alignment follow the path here?
                bool match = true;
                for (size_t i = 1; i < path.mapping_size(); i++) {
                    const Position& pos = path.mapping(i).position();
                    if (pos.node_id() != xpath.node(occurrence.first + i) || pos.is_reverse() != xpath.is_reverse(occurrence.first + i)) {
                        match = false;
                        break;
                    }
                }
                
                // we found where the alignment could be from
                if (match) {
                    path_pos_out = xpath.positions[occurrence.first] + start_pos.offset();
                    path_rev_out = false;
                    return;
                }
            }
            else {
                // the first node in this alignment occurs on the reverse strand of the path
                
                if (occurrence.first + 1 < path.mapping_size()) {
                    // but it doesn't fit on the path
                    continue;
                }
                
                // does the alignment follow the path here?
                bool match = true;
                for (size_t i = 1; i < path.mapping_size(); i++) {
                    const Position& pos = path.mapping(i).position();
                    if (pos.node_id() != xpath.node(occurrence.first - i) || pos.is_reverse() == xpath.is_reverse(occurrence.first - i)) {
                        match = false;
                        break;
                    }
                }
                
                // we found where the alignment could be from
                if (match) {
                    const Mapping& last_mapping = path.mapping(path.mapping_size() - 1);
                    size_t last_offset = occurrence.first + 1 - path.mapping_size();
                    int64_t node_start = last_offset + 1 < xpath.positions.size() ? xpath.positions[last_offset + 1] : xpath.offsets.size();
                    path_pos_out = node_start - last_mapping.position().offset() - mapping_from_length(last_mapping);
                    path_rev_out = true;
                    return;
                }
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


