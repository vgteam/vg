/**
 * \file
 * surjector.cpp: implements a class that surjects alignments onto paths
 */

#include "surjector.hpp"


namespace vg {

using namespace std;
    
    Surjector::Surjector(xg::XG* xg_index) : Mapper(xg_index, nullptr, nullptr) {
        
    }
    
    Surjector::~Surjector() {
        
    }

    Alignment Surjector::surject_classic(const Alignment& source,
                                         const set<string>& path_names,
                                         string& path_name,
                                         int64_t& path_pos,
                                         bool& path_reverse) {
        
        Alignment surjection = source;
        // Leave the original mapping quality in place (because that's the quality
        // on the placement of this read in this region at all)
        surjection.clear_mapping_quality();
        surjection.clear_score();
        surjection.clear_identity();
        surjection.clear_path();
        
        int count_forward=0, count_reverse=0;
        for (auto& mapping : source.path().mapping()) {
            if (mapping.position().is_reverse()) {
                ++count_reverse;
            } else {
                ++count_forward;
            }
        }
        //cerr << "fwd " << count_forward << " rev " << count_reverse << endl;
        
        // here we assume that people will use this on DAGs
        // require that we have an alignment with a score, and that it is on one strand
        if (!source.has_path() || source.path().mapping_size() == 0
            || alignment_from_length(source) == 0
            || count_forward > 0 && count_reverse > 0) {
#ifdef debug_mapper
            
#pragma omp critical (cerr)
            cerr << "Alignment " << source.name() << " is unmapped and cannot be surjected" << endl;
            
#endif
            return surjection;
        }
        
        set<id_t> nodes;
        for (int i = 0; i < source.path().mapping_size(); ++ i) {
            nodes.insert(source.path().mapping(i).position().node_id());
        }
        VG graph;
        for (auto& node : nodes) {
            *graph.graph.add_node() = xindex->node(node);
        }
        xindex->expand_context(graph.graph, 3, true); // get connected edges and path
        graph.paths.append(graph.graph);
        graph.rebuild_indexes();
        VG base_graph = graph;
        
        // non-fiddly approach, rest on augmentation
        // 0) remove softclips from the read
        // 1) augment the graph with the read
        // 2) keep the ref path and the aln path both in the graph
        // 3) detach the nodes on the other sides of the aln path start and end from all other nodes
        // 4) remove the non-path component
        
        Alignment trimmed_source = strip_from_end(strip_from_start(source, softclip_start(source)), softclip_end(source));
        // check if we'd fail
        if (trimmed_source.sequence().size() == 0) {
            return surjection;
        }
        
        vector<Path> source_path;
        source_path.push_back(trimmed_source.path());
        source_path.back().set_name(source.name());
        // Make sure to pass true here to embed the alignment
        auto translation = graph.edit(source_path, true); //, true, true);
        Translator translator(translation);
        Path source_in_graph = graph.paths.path(source.name());
        Position start_pos = make_position(initial_position(source_in_graph));
        Position end_pos = make_position(final_position(source_in_graph));
        //cerr << "start and end pos " << pb2json(start_pos) << " " << pb2json(end_pos) << endl;
        
        //Position start_pos = make_position(initial_position(trimmed_source.path()));
        //Position end_pos = make_position(final_position(trimmed_source.path()));
        
        // find then unlink the next and previous path nodes from the rest of the graph to isolate the path-specific component
        handle_t start = graph.get_handle(start_pos.node_id(), start_pos.is_reverse());
        handle_t end = graph.get_handle(end_pos.node_id(), end_pos.is_reverse());
        handle_t cut;
        bool found = false;
        unordered_set<handle_t> curr;
        unordered_set<handle_t> next;
        auto find_path = [&](const handle_t& h) {
            vector<string> path_intersection;
            set<string> node_paths = graph.paths.of_node(graph.get_id(h));
            //cerr << "Node paths for " << graph.get_id(h) << " " << node_paths.size() << endl;
            if (!node_paths.empty()) {
                std::set_intersection(path_names.begin(), path_names.end(),
                                      node_paths.begin(), node_paths.end(),
                                      std::back_inserter(path_intersection));
            }
            cut = h;
            found = path_intersection.size() > 0;
            //cerr << "path intersection size " << path_intersection.size() << endl;
            next.insert(h);
            return !found;
        };
        found = false;
        curr.insert(start);
        //cerr << "going back" << endl;
        while (!curr.empty()) {
            bool finished = false;
            //cerr << "cur has " << curr.size() << endl;
            for (auto& h : curr) {
                finished |= !graph.follow_edges(h, true, find_path);
                if (finished) break;
            }
            if (finished) {
                curr.clear();
                next.clear();
                break;
            } else {
                curr = next;
                next.clear();
            }
        }
        handle_t cut_before = cut;
        bool found_forward = found;
        curr.insert(end);
        //ncerr << "going forward" << endl;
        while (!curr.empty()) {
            bool finished = false;
            //cerr << "cur has " << curr.size() << endl;
            for (auto& h : curr) {
                finished |= !graph.follow_edges(h, false, find_path);
                if (finished) break;
            }
            if (finished) {
                curr.clear();
                next.clear();
                break;
            } else {
                curr = next;
                next.clear();
            }
        }
        handle_t cut_after = cut;
        //cerr << "cut before " << graph.get_id(cut_before) << endl;
        //cerr << "cut after " << graph.get_id(cut_after) << endl;
        bool found_reverse = found;
        //graph.serialize_to_file("before-" + source.name() + ".vg");
        //graph.serialize_to_file("before-" + graph.hash() + ".vg");
        
        set<string> kept_paths;
        graph.keep_paths(path_names, kept_paths);
        graph.remove_non_path();
        // by definition we have found path
        if (found_forward && found_reverse && cut_before == cut_after) {
            graph.destroy_handle(cut_before);
        } else {
            if (found_forward) graph.destroy_handle(cut_before);
            if (found_reverse) graph.destroy_handle(cut_after);
        }
        //graph.serialize_to_file("after-" + source.name() + ".vg");
        //graph.serialize_to_file("after-" + graph.hash() + ".vg");
        
#ifdef debug_mapper
        cerr << "src " << pb2json(source) << endl;
        cerr << "start " << pb2json(start_pos) << endl;
        cerr << "end " << pb2json(end_pos) << endl;
        cerr << "graph " << pb2json(graph.graph) << endl;
#endif
        //Position end_pos = alignment_end(source);
        // assume DAG
        set<vg::id_t> target_ids;
        for (auto& mapping : source_in_graph.mapping()) {
            target_ids.insert(mapping.position().node_id());
        }
        
        // otherwise, two cuts
        // remove the links in both cases
        // we can clean up by removing
        
        // get only the subgraph that we want to align to
        list<VG> subgraphs;
        graph.disjoint_subgraphs(subgraphs);
        
        // Align the old alignment to the graph in both orientations. Apparently
        // align only does a single oriantation, and we have no idea, even looking
        // at the mappings, which of the orientations will correspond to the one the
        // alignment is actually in.
        
        Graph subgraph;
        for (auto& graph : subgraphs) {
            //cerr << pb2json(graph.graph) << endl;
            bool found = false;
            graph.for_each_handle([&](const handle_t& h) {
                if (!found && target_ids.count(graph.get_id(h))) {
                    found = true;
                }
            });
            if (found) {
                subgraph = graph.graph;
                break;
            }
        }
        if (subgraph.node_size() == 0) {
            // couldn't find subgraph, try the one we've got
            subgraph = graph.graph;
        }
        
        if (subgraph.node_size() == 0) {
            return surjection; //empty graph, avoid further warnings
        }
        
        // DAG assumption
        sort_by_id_dedup_and_clean(subgraph);
#ifdef debug_mapper
        cerr << "sub " << pb2json(subgraph) << endl;
#endif
        
        // Flip the string and its quality around
        Alignment surjection_fwd = surjection;
        Alignment surjection_rev = surjection;
        int start_softclip_length = softclip_start(source);
        int end_softclip_length = softclip_end(source);
        Alignment start_softclip_fwd, end_softclip_fwd;
        Alignment start_softclip_rev, end_softclip_rev;
        if (start_softclip_length) {
            start_softclip_fwd = strip_from_end(surjection_fwd, surjection_fwd.sequence().size() - start_softclip_length);
            end_softclip_rev = reverse_complement_alignment(start_softclip_fwd, [&](id_t id) { return xindex->node_length(id); });
        }
        if (end_softclip_length) {
            end_softclip_fwd = strip_from_start(surjection_fwd, surjection_fwd.sequence().size() - end_softclip_length);
            start_softclip_rev = reverse_complement_alignment(end_softclip_fwd, [&](id_t id) { return xindex->node_length(id); });
        }
        surjection_fwd = strip_from_end(strip_from_start(surjection_fwd, start_softclip_length), end_softclip_length);
        surjection_rev = reverse_complement_alignment(surjection_fwd, [&](id_t id) { return xindex->node_length(id); });
        
        // align to the graph with a big full len, and simplify without removal of internal deletions, as we'll need these for BAM reconstruction
        Alignment surjection_forward, surjection_reverse;
        int fwd_score = 0, rev_score = 0;
        // override the full length bonus
        int8_t full_length_bonus_override = 30;
        int8_t saved_bonus = get_aligner(!surjection.quality().empty())->full_length_bonus;
        get_aligner(!surjection.quality().empty())->full_length_bonus = full_length_bonus_override;
        if (count_forward) {
            surjection_forward = simplify(align_to_graph(surjection_fwd, graph.graph, max_query_graph_ratio, true, false, false, false, false, false), false);
            fwd_score = surjection_forward.score();
        }
        if (count_reverse) {
            surjection_reverse = simplify(align_to_graph(surjection_rev, graph.graph, max_query_graph_ratio, true, false, false, false, false, false), false);
            rev_score = surjection_reverse.score();
        }
        // reset bonus because hacks
        get_aligner(!surjection.quality().empty())->full_length_bonus = saved_bonus;
        
#ifdef debug_mapper
        cerr << "fwd " << pb2json(surjection_forward) << endl;
        cerr << "rev " << pb2json(surjection_reverse) << endl;
#endif
        
        graph = base_graph;
        // We need this for inverting mappings to the correct strand
        function<int64_t(id_t)> node_length = [&graph](id_t node) {
            return graph.get_node(node)->sequence().size();
        };
        
#ifdef debug_mapper
#pragma omp critical (cerr)
        cerr << surjection.name() << " " << surjection_forward.score() << " forward score, " << surjection_reverse.score() << " reverse score" << endl;
#endif
        
        // translate
        try {
            if (count_forward) surjection_forward = translator.translate(surjection_forward);
            if (count_reverse) surjection_reverse = translator.translate(surjection_reverse);
        } catch (...) {
            cerr << "[vg Mapper::surject_alignment] warning: surjection failure with read " << source.name() << endl;
            return surjection;
        }
        
        // patch
        if (count_forward) patch_alignment(surjection_forward, surjection_forward.sequence().size(), false);
        if (count_reverse) patch_alignment(surjection_reverse, surjection_reverse.sequence().size(), false);
        
        // reattach soft clips and set original score (score isn't really used through...)
        if (count_forward) {
            surjection_forward = merge_alignments({start_softclip_fwd, surjection_forward, end_softclip_fwd});
            surjection_forward.set_score(fwd_score);
        }
        if (count_reverse) {
            surjection_reverse = merge_alignments({start_softclip_rev, surjection_reverse, end_softclip_rev});
            surjection_reverse.set_score(rev_score);
        }

        // choose
        if (count_reverse && count_forward) {
            if (surjection_reverse.score() > surjection_forward.score()) {
                surjection = reverse_complement_alignment(surjection_reverse, [&](id_t id) { return xindex->node_length(id); });
            } else {
                surjection = surjection_forward;
            }
        } else {
            if (count_reverse) {
                surjection = reverse_complement_alignment(surjection_reverse, [&](id_t id) { return xindex->node_length(id); });
            } else {
                surjection = surjection_forward;
            }
        }
        surjection = simplify(surjection, false);
        /*
        cerr << "surj " << pb2json(surjection) << endl;
        for (size_t j = 0; j < surjection.path().mapping_size(); ++j) {
            Mapping* m = surjection.mutable_path()->mutable_mapping(j);
            for (size_t i = 0; i < m->edit_size(); ++i) {
                cerr << pb2json(m->edit(i)) << endl;
                assert(!edit_is_empty(m->edit(i)));
            }
        }
        */

#ifdef debug_mapper
        
#pragma omp critical (cerr)
        cerr << surjection.path().mapping_size() << " mappings, " << kept_paths.size() << " paths" << endl;
        
#endif
        //assert(check_alignment(surjection));
        if (surjection.path().mapping_size() > 0 && kept_paths.size() == 1) {
            // determine the paths of the node we mapped into
            //  ... get the id of the first node, get the paths of it
            assert(kept_paths.size() == 1);
            path_name = *kept_paths.begin();
            
            int64_t path_id = xindex->path_rank(path_name);
            auto& first_pos = surjection.path().mapping(0).position();
            int64_t hit_id = surjection.path().mapping(0).position().node_id();
            bool hit_backward = surjection.path().mapping(0).position().is_reverse();
            // we pick up positional information using the index
            
            //cerr << "hit id " << hit_id << endl;
            auto path_posns = xindex->position_in_path(hit_id, path_name);
            if (path_posns.size() > 1) {
                cerr << "[vg map] surject_alignment: warning, multiple positions for node " << hit_id << " in " << path_name << " but will use only first: " << path_posns.front() << endl;
            } else if (path_posns.size() == 0) {
                cerr << "[vg map] surject_alignment: error, no positions for alignment " << source.name() << endl;
                exit(1);
            }
            
            // if we are reversed
            path_pos = path_posns.front();
            bool reversed_path = xindex->mapping_at_path_position(path_name, path_pos).position().is_reverse();
            if (reversed_path) {
                // if we got the start of the node position relative to the path
                // we need to offset to make things right
                // but which direction
                if (hit_backward) {
                    path_pos = path_posns.front() + first_pos.offset();
                } else {
                    auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                    path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
                }
                path_reverse = !hit_backward;
            } else {
                if (!hit_backward) {
                    path_pos = path_posns.front() + first_pos.offset();
                } else {
                    auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                    path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
                }
                path_reverse = hit_backward;
            }
            
#ifdef debug_mapper
            cerr << "path position " << path_name << ":" << path_pos << endl;
#endif
            
        } else {
            
            surjection = source;
#ifdef debug_mapper
            
#pragma omp critical (cerr)
            cerr << "Alignment " << source.name() << " did not align to the surjection subgraph" << endl;
            
#endif
            
        }
        
#ifdef debug_mapper
        
#pragma omp critical (cerr)
        cerr << "Surjection on reverse strand? " << path_reverse << endl;
        cerr << "Surjected alignment: " << pb2json(surjection) << endl;
        
#endif
        //cerr << "final " << pb2json(surjection) << endl;
        Alignment final = source;
        *final.mutable_path() = surjection.path();
        final.set_score(surjection.score());
        final.set_identity(identity(surjection.path()));
        return final;
    }
    
    Alignment Surjector::path_anchored_surject(const Alignment& source, const string& path_name) {
        
        // the path we're surecting onto
        size_t path_rank = xindex->path_rank(path_name);
        const xg::XGPath& xpath = xindex->get_path(path_name);
        
        // memos for expensive succinct operations that may be repeated
        unordered_map<int64_t, vector<size_t>> paths_of_node_memo;
        unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>> oriented_occurrences_memo;
        
        // get the chunks of the aligned path that overlap the ref path
        auto path_overlapping_anchors = extract_overlapping_paths(source, path_rank, &paths_of_node_memo);
        
        // find the interval of the ref path we need to cnsider
        pair<size_t, size_t> ref_path_interval = compute_path_interval(source, path_rank, xpath, path_overlapping_anchors,
                                                                       &oriented_occurrences_memo);
        
        // get the path graph corresponding to this interval
        unordered_map<id_t, pair<id_t, bool>> node_trans;
        VG path_graph = extract_linearized_path_graph(ref_path_interval.first, ref_path_interval.second, xpath, node_trans);
        
        // compute the connectivity between the path chunks
        MultipathAlignmentGraph mp_aln_graph(path_graph, path_overlapping_anchors, node_trans);
        
        // TODO: is this necessary?
        vector<size_t> topological_order;
        mp_aln_graph.topological_sort(topological_order);
        mp_aln_graph.remove_transitive_edges(topological_order);
        
        // align the intervening segments and store the result in a multipath alignment
        MultipathAlignment mp_aln;
        mp_aln_graph.align(source, path_graph, get_aligner(), false, 1, 1, mp_aln);
        
        // concatenate the subpaths
        Alignment surjected;
        optimal_alignment(mp_aln, surjected);
        
        // translate back into the original ID space and return
        translate_oriented_node_ids(*surjected.mutable_path(), node_trans);
        return surjected;
    }
    
    vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>
    Surjector::extract_overlapping_paths(const Alignment& source, size_t path_rank,
                                         unordered_map<int64_t, vector<size_t>>* paths_of_node_memo) {
        
        vector<pair<pair<string::const_iterator, string::const_iterator>, Path>> to_return;
        
        const Path& path = source.path();
        
        bool currently_on_path = false;
        int64_t through_to_length = 0;
        for (size_t i = 0; i < path.mapping_size(); i++) {
            
            int64_t before_to_length = through_to_length;
            through_to_length += mapping_to_length(path.mapping(i));
            
            vector<size_t> paths_of_node = xindex->memoized_paths_of_node(path.mapping(i).position().node_id(),
                                                                          paths_of_node_memo);
            
            size_t j = 0;
            for (; j < paths_of_node.size(); j++) {
                if (paths_of_node[j] == path_rank) {
                    break;
                }
            }
            
            if (j == paths_of_node.size()) {
                currently_on_path = false;
            }
            else {
                if (!currently_on_path) {
                    to_return.emplace_back();
                    to_return.back().first.first = source.sequence().begin() + before_to_length;
                }
                
                to_return.back().first.second = source.sequence().begin() + through_to_length;
                *to_return.back().second.add_mapping() = path.mapping(i);
                currently_on_path = true;
            }
        }
        
        return to_return;
    }
    
    pair<size_t, size_t>
    Surjector::compute_path_interval(const Alignment& source, size_t path_rank, const xg::XGPath& xpath,
                                     const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                     unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo) {
        
        pair<size_t, size_t> interval(numeric_limits<size_t>::max(), numeric_limits<size_t>::min());
        
        for (const auto& path_chunk : path_chunks) {
            
            string::const_iterator read_pos = path_chunk.first.first;
            
            for (size_t i = 0; i < path_chunk.second.mapping_size(); i++) {
                const Position& pos = path_chunk.second.mapping(i).position();
                
                int64_t left_gap = get_aligner()->longest_detectable_gap(source, read_pos);
                int64_t left_overhang = left_gap + (read_pos - source.sequence().begin());
                
                read_pos += mapping_to_length(path_chunk.second.mapping(i));
                
                int64_t right_gap = get_aligner()->longest_detectable_gap(source, read_pos);
                int64_t right_overhang = right_gap + (source.sequence().end() - read_pos);
                
                auto oriented_occurrences = xindex->memoized_oriented_occurrences_on_path(pos.node_id(), path_rank,
                                                                                          oriented_occurrences_memo);
                
                int64_t mapping_length = mapping_from_length(path_chunk.second.mapping(i));
                
                for (const pair<size_t, bool>& occurrence : oriented_occurrences) {
                    if (occurrence.second == pos.is_reverse()) {
                        int64_t path_offset = xpath.positions[occurrence.first];
                        
                        int64_t left_boundary = max<int64_t>(0, path_offset + pos.offset() - left_overhang);
                        interval.first = min<size_t>(interval.first, left_boundary);
                        
                        int64_t right_boundary = min<int64_t>(path_offset + mapping_length + right_overhang, xpath.offsets.size());
                        interval.second = max<size_t>(interval.second, right_boundary);
                    }
                    else {
                        int64_t path_offset = occurrence.first + 1 < xpath.positions.size() ? xpath.positions[occurrence.first + 1] : xpath.offsets.size();
                        
                        int64_t left_boundary = max<int64_t>(0, path_offset - mapping_length - right_overhang);
                        interval.first = min<size_t>(interval.first, left_boundary);
                        
                        int64_t right_boundary = min<int64_t>(path_offset - pos.offset() + left_overhang, xpath.offsets.size());
                        interval.second = max<size_t>(interval.second, right_boundary);
                    }
                }
            }
        }
        
        return interval;
    }
    
    VG Surjector::extract_linearized_path_graph(size_t begin, size_t end, const xg::XGPath& xpath,
                                                unordered_map<id_t, pair<id_t, bool>>& node_trans) {
        
        VG path_graph;
        
        size_t first = xpath.offset_at_position(begin);
        size_t last = min<size_t>(xpath.positions.size(), xpath.offset_at_position(end) + 1);
        
        Node* prev_node = nullptr;
        for (size_t i = first; i < last; i++) {
            
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
}


