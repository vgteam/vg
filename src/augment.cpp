#include "vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/alignment_emitter.hpp>

#include "augment.hpp"
#include "alignment.hpp"
#include "packer.hpp"
#include "annotation.hpp"
//#define debug

using namespace vg::io;

namespace vg {

using namespace std;

// The correct way to edit the graph
void augment(MutablePathMutableHandleGraph* graph,
             const string& gam_path,
             const string& aln_format,
             vector<Translation>* out_translations,
             const string& gam_out_path,
             bool embed_paths,
             bool break_at_ends,
             bool remove_softclips,
             bool filter_out_of_graph_alignments,
             double min_baseq,
             double min_mapq,
             Packer* packer,
             size_t min_bp_coverage,
             double max_frac_n,
             bool edges_only) {

    // memory-wasting hack: we need node lengths from the original graph in order to parse the GAF.  Unlesss we
    // store them, they will be lost in the 2nd pass
    unordered_map<nid_t, int64_t> id_to_length;
    if (aln_format == "GAF") {
        graph->for_each_handle([&](handle_t handle) {
                id_to_length[graph->get_id(handle)] = graph->get_length(handle);
            });
    }

    function<void(function<void(Alignment&)>, bool, bool)> iterate_gam =
        [&gam_path, &aln_format, &graph, &packer, &id_to_length] (function<void(Alignment&)> aln_callback, bool second_pass, bool parallel) {
        if (aln_format == "GAM") {
            get_input_file(gam_path, [&](istream& gam_stream) {
                    if (parallel) {
                        vg::io::for_each_parallel(gam_stream, aln_callback, Packer::estimate_batch_size(get_thread_count()));
                    } else {
                        vg::io::for_each(gam_stream, aln_callback);
                    }
                });
        } else {
            assert(aln_format == "GAF");
            function<size_t(nid_t)> node_to_length;
            function<string(nid_t, bool)> node_to_sequence;
            if (second_pass) {
                // graph has changed, need to fall back on our table we saved from the original graph
                node_to_length = [&id_to_length](nid_t node_id) {
                    return id_to_length[node_id];
                };
                // try to do without sequences
                node_to_sequence = nullptr;
            } else {
                // graph is valid on the first pass
                node_to_length = [&graph](nid_t node_id) {
                    return graph->get_length(graph->get_handle(node_id));
                };
                node_to_sequence = [&graph](nid_t node_id, bool is_reversed) {
                    return graph->get_sequence(graph->get_handle(node_id, is_reversed));
                };
            }
            if (parallel) {
                vg::io::gaf_unpaired_for_each_parallel(node_to_length, node_to_sequence, gam_path, aln_callback);
            } else {
                vg::io::gaf_unpaired_for_each(node_to_length, node_to_sequence, gam_path, aln_callback);
            }
        }
    };

    augment_impl(graph,                 
                 iterate_gam,
                 aln_format,
                 out_translations,
                 gam_out_path,
                 embed_paths,
                 break_at_ends,
                 remove_softclips,
                 filter_out_of_graph_alignments,
                 min_baseq,
                 min_mapq,                 
                 packer,
                 min_bp_coverage,
                 max_frac_n,
                 edges_only);
}

void augment(MutablePathMutableHandleGraph* graph,
             vector<Path>& path_vector,
             const string& aln_format,
             vector<Translation>* out_translations,
             const string& gam_out_path,
             bool embed_paths,
             bool break_at_ends,
             bool remove_softclips,
             bool filter_out_of_graph_alignments,
             double min_baseq,
             double min_mapq,
             Packer* packer,
             size_t min_bp_coverage,
             double max_frac_n,
             bool edges_only) {
    
    function<void(function<void(Alignment&)>, bool, bool)> iterate_gam =
        [&path_vector] (function<void(Alignment&)> aln_callback, bool second_pass, bool parallel) {
        if (parallel) {
#pragma omp parallel for
            for (size_t i = 0; i < path_vector.size(); ++i) {
                Path& path = path_vector[i];
                Alignment aln;
                *aln.mutable_path() = path;
                aln.set_name(path.name());
                aln_callback(aln);
            }
            
        }
        else {
            for (Path& path : path_vector) {
                Alignment aln;
                *aln.mutable_path() = path;
                aln.set_name(path.name());
                aln_callback(aln);
            }
        }
    };

    augment_impl(graph,
                 iterate_gam,
                 aln_format,
                 out_translations,
                 gam_out_path,
                 embed_paths,
                 break_at_ends,
                 remove_softclips,
                 filter_out_of_graph_alignments,
                 min_baseq,
                 min_mapq,
                 packer,
                 min_bp_coverage,
                 max_frac_n,
                 edges_only);
}

// Check if alignment contains node that's not in the graph
static inline bool check_in_graph(const Path& path, HandleGraph* graph) {
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        if (!graph->has_node(path.mapping(i).position().node_id())) {
            return false;
        }
    }
    return true;
}

// Check if alignment contains node that's not in the graph (via node sizes map)
static inline bool check_in_graph(const Path& path, const unordered_map<id_t, size_t>& node_map) {
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        if (!node_map.count(path.mapping(i).position().node_id())) {
            return false;
        }
    }
    return true;
}

void augment_impl(MutablePathMutableHandleGraph* graph,
                  function<void(function<void(Alignment&)>,bool, bool)> iterate_gam,
                  const string& aln_format,                  
                  vector<Translation>* out_translations,
                  const string& gam_out_path,
                  bool embed_paths,
                  bool break_at_ends,
                  bool remove_softclips,
                  bool filter_out_of_graph_alignments,
                  double min_baseq,
                  double min_mapq,
                  Packer* packer,
                  size_t min_bp_coverage,
                  double max_frac_n,
                  bool edges_only) {

    if (edges_only) {
        // just add edges between consecutive mappings that aren't already in the graph
        // note: offsets are completely ignored.  if we want to take them into account,
        // it's probably best to use below logic, but filtering non-deletion edits out
        // of the GAM first.  but keeping just the node-adjacency information as we
        // do here is hopefully sufficient to improve sv genotyping with pack/call
        add_edges_only(graph, iterate_gam, min_mapq, min_bp_coverage);
        return;
    }
    
    // toggle between using Packer to store breakpoints or the STL map
    bool packed_mode = min_bp_coverage > 0 || min_baseq > 0 || max_frac_n < 1.;
    assert(!packed_mode || packer != nullptr);
    
    unordered_map<id_t, set<pos_t>> breakpoints;
        
    // First pass: find the breakpoints
    iterate_gam((function<void(Alignment&)>)[&](Alignment& aln) {
#ifdef debug
            cerr << pb2json(aln.path()) << endl;
#endif
            if (aln.mapping_quality() < min_mapq || (filter_out_of_graph_alignments && !check_in_graph(aln.path(), graph))) {
                return;
            }

            if (aln_format == "GAF" && has_annotation(aln, "from_cg") && get_annotation<bool>(aln, "from_cg")) {
#pragma omp critical (cerr)
                {
                    cerr << "[vg augment] error: GAF with cg cigars contains insufficient information for augmenting: cs cigars required." << endl;
                }
                exit(1);
            }                    

            if (remove_softclips) {
                softclip_trim(aln);
            }

            // Simplify the path, just to eliminate adjacent match Edits in the same
            // Mapping (because we don't have or want a breakpoint there)
            Path simplified_path = simplify(aln.path());

            // Add in breakpoints from each path
            if (packed_mode) {
                find_packed_breakpoints(simplified_path, *packer, break_at_ends, aln.quality(), min_baseq, max_frac_n);
            } else {
                // note: we cannot pass non-zero min_baseq here.  it relies on filter_breakpoints_by_coverage
                // to work correctly, and must be passed in only via find_packed_breakpoints.
                find_breakpoints(simplified_path, breakpoints, break_at_ends, "", 0, 1.);
            }
        }, false, packed_mode);

    if (packed_mode) {
        // Filter the breakpoints by coverage
        breakpoints = filter_breakpoints_by_coverage(*packer, min_bp_coverage);
    } else {
        // Invert the breakpoints that are on the reverse strand
        breakpoints = forwardize_breakpoints(graph, breakpoints);
    }

    // don't need this anymore: free up some memory
    if (packer != nullptr) {
        packer->clear();
    }

    // get the node sizes, for use when making the translation
    unordered_map<id_t, size_t> orig_node_sizes;
    orig_node_sizes.reserve(graph->get_node_count());
    graph->for_each_handle([&](handle_t node) {
            orig_node_sizes[graph->get_id(node)] = graph->get_length(node);
        });

    // Break any nodes that need to be broken. Save the map we need to translate
    // from offsets on old nodes to new nodes. Note that this would mess up the
    // ranks of nodes in their existing paths, which is why we clear and rebuild
    // them.
    auto node_translation = ensure_breakpoints(graph, breakpoints);

    // we remember the sequences of nodes we've added at particular positions on the forward strand
    unordered_map<pair<pos_t, string>, vector<id_t>> added_seqs;
    // we will record the nodes that we add, so we can correctly make the returned translation
    unordered_map<id_t, Path> added_nodes;
    // output alignment emitter and buffer
    unique_ptr<vg::io::AlignmentEmitter> aln_emitter;
    if (!gam_out_path.empty()) {
        aln_emitter = vg::io::get_non_hts_alignment_emitter(gam_out_path, aln_format, {}, get_thread_count(), graph);
    }
    vector<Alignment> aln_buffer;

    // Second pass: add the nodes and edges
    iterate_gam((function<void(Alignment&)>)[&](Alignment& aln) {
            if (aln.mapping_quality() < min_mapq || (filter_out_of_graph_alignments && !check_in_graph(aln.path(), orig_node_sizes))) {
                return;
            }
            
            if (remove_softclips) {
                softclip_trim(aln);
            }

            // Simplify the path, just to eliminate adjacent match Edits in the same
            // Mapping (because we don't have or want a breakpoint there)
            // Note: We're electing to re-simplify in a second pass to avoid storing all
            // the input paths in memory
            Path simplified_path = simplify(aln.path());

            // Filter out edits corresponding to breakpoints that didn't meet our coverage
            // criteria
            if (min_bp_coverage > 0) {
                simplify_filtered_edits(graph, aln, simplified_path, node_translation, orig_node_sizes,
                                        min_baseq, max_frac_n);
            }

            // Create new nodes/wire things up. Get the added version of the path.
            Path added = add_nodes_and_edges(graph, simplified_path, node_translation, added_seqs,
                                             added_nodes, orig_node_sizes);

            // Copy over the name
            *added.mutable_name() = aln.name();

            if (embed_paths) {
                add_path_to_graph(graph, added);
            }

            // something is off about this check.
            // assuming the GAM path is sorted, let's double-check that its edges are here
            for (size_t i = 1; i < added.mapping_size(); ++i) {
                auto& m1 = added.mapping(i-1);
                auto& m2 = added.mapping(i);
                // we're no longer sorting our input paths, so we assume they are sorted
                assert((m1.rank() == 0 && m2.rank() == 0) || (m1.rank() + 1 == m2.rank()));
                //if (!adjacent_mappings(m1, m2)) continue; // the path is completely represented here
                auto s1 = graph->get_handle(m1.position().node_id(), m1.position().is_reverse());
                auto s2 = graph->get_handle(m2.position().node_id(), m2.position().is_reverse());
                // Ensure that we always have an edge between the two nodes in the correct direction
                graph->create_edge(s1, s2);
            }

            // optionally write out the modified path to GAM
            if (!gam_out_path.empty()) {
                *aln.mutable_path() = added;
                aln_buffer.push_back(aln);
                if (aln_buffer.size() >= 100) {
                    aln_emitter->emit_singles(vector<Alignment>(aln_buffer));
                    aln_buffer.clear();
                }
            }
        }, true, false);
    if (!aln_buffer.empty()) {
        // Flush the buffer
        aln_emitter->emit_singles(vector<Alignment>(aln_buffer));
    }

    // perform the same check as above, but on the paths that were already in the graph
    // assuming the graph's paths are sorted, let's double-check that the edges are here
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            step_handle_t prev_handle;
            int i = 0;
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    handle_t handle = graph->get_handle_of_step(step_handle);
                    if (i > 0) {
                        // Ensure the edge that the path follows exists.
                        // TODO: Should this be an error if it doesn't exist instead?
                        graph->create_edge(graph->get_handle_of_step(prev_handle), handle);
                    }
                    prev_handle = step_handle;
                });
        });

    // make the translation
    if (out_translations != nullptr) {
        *out_translations = make_translation(graph, node_translation, added_nodes, orig_node_sizes);
    }

    VG* vg_graph = dynamic_cast<VG*>(graph);
    
    // This code got run after augment in VG::edit, so we make sure it happens here too
    if (vg_graph != nullptr) {
        // Rebuild path ranks, aux mapping, etc. by compacting the path ranks
        // Todo: can we just do this once?
       vg_graph->paths.compact_ranks();
        
       // execute a semi partial order sort on the nodes
       vg_graph->sort();
    }

}

double get_avg_baseq(const Edit& edit, const string& base_quals, size_t position_in_read) {
    double avg_qual = numeric_limits<int>::max();
    if (!base_quals.empty() && !edit.sequence().empty() && (edit_is_sub(edit) || edit_is_insertion(edit))) {
        double tot_qual = 0;
        for (int i = 0; i < edit.sequence().length(); ++i) {
            tot_qual += base_quals[position_in_read + i];
        }
        avg_qual = tot_qual / (double)edit.sequence().length();
    }
    return avg_qual;
}

// returns breakpoints on the forward strand of the nodes
void find_breakpoints(const Path& path, unordered_map<id_t, set<pos_t>>& breakpoints, bool break_ends,
                      const string& base_quals, double min_baseq, double max_frac_n) {
    // We need to work out what offsets we will need to break each node at, if
    // we want to add in all the new material and edges in this path.

#ifdef debug
    cerr << "Processing path..." << endl;
#endif

    // The base position in the edit
    size_t position_in_read = 0;

    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        const Mapping& m = path.mapping(i);

        // What node are we on?
        id_t node_id = m.position().node_id();

        if(node_id == 0) {
            // Skip Mappings that aren't actually to nodes.
            continue;
        }

        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        pos_t edit_first_position = make_pos_t(m.position());
    
#ifdef debug
        cerr << "Processing mapping " << pb2json(m) << endl;
#endif

        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            const Edit& e = m.edit(j);

            // We know where the mapping starts in its node. But where does it
            // end (inclusive)? Note that if the edit has 0 reference length,
            // this may not actually be included in the edit (and
            // edit_first_position will be further along than
            // edit_last_position).
            pos_t edit_last_position = edit_first_position;
            if (e.from_length()) {
                get_offset(edit_last_position) += e.from_length();
            }

#ifdef debug
            cerr << "Edit on " << node_id << " from " << edit_first_position << " to " << edit_last_position << endl;
            cerr << pb2json(e) << endl;
#endif

            // Do the base quality check if applicable.  If it fails we just ignore the edit
            if ((min_baseq == 0 || get_avg_baseq(e, base_quals, position_in_read) >= min_baseq) &&
                (max_frac_n == 1. || get_fraction_of_ns(e.sequence()) <= max_frac_n)) {

                
                if (!edit_is_match(e) || (j == 0 && (i != 0 || break_ends))) {
                    // If this edit is not a perfect match, or if this is the first
                    // edit in this mapping and either we had a previous mapping we
                    // may need to connect to or we want to break at the path's
                    // start, we need to make sure we have a breakpoint at the start
                    // of this edit.

#ifdef debug
                    cerr << "Need to break " << node_id << " at edit lower end " <<
                        edit_first_position << endl;
#endif

                    // We need to snip between edit_first_position and edit_first_position - direction.
                    // Note that it doesn't matter if we put breakpoints at 0 and 1-past-the-end; those will be ignored.
                    breakpoints[node_id].insert(edit_first_position);
                }

                if (!edit_is_match(e) || (j == m.edit_size() - 1 && (i != path.mapping_size() - 1 || break_ends))) {
                    // If this edit is not a perfect match, or if it is the last
                    // edit in a mapping and we have a subsequent mapping we might
                    // need to connect to or we want to break at the path ends, make
                    // sure we have a breakpoint at the end of this edit.

#ifdef debug
                    cerr << "Need to break " << node_id << " at past edit upper end " <<
                        edit_last_position << endl;
#endif

                    // We also need to snip between edit_last_position and edit_last_position + direction.
                    breakpoints[node_id].insert(edit_last_position);
                }
            }
            // TODO: for an insertion or substitution, note that we need a new
            // node and two new edges.

            // TODO: for a deletion, note that we need an edge. TODO: Catch
            // and complain about some things we can't handle (like a path with
            // a leading/trailing deletion)? Or just skip deletions when wiring.

            // Use up the portion of the node taken by this mapping, so we know
            // where the next mapping will start.
            edit_first_position = edit_last_position;

            position_in_read += e.to_length();
        }
    }

}

unordered_map<id_t, set<pos_t>> forwardize_breakpoints(const HandleGraph* graph,
                                                       const unordered_map<id_t, set<pos_t>>& breakpoints) {
    unordered_map<id_t, set<pos_t>> fwd;
    for (auto& p : breakpoints) {
        id_t node_id = p.first;
        if (!graph->has_node(node_id)) {
            throw runtime_error("Node from GAM \"" + std::to_string(node_id) + "\" not found in graph.  If you are sure"
                                " the input graph is a subgraph of that used to create the GAM, you can ignore this error"
                                " with \"vg augment -s\"");
        }
        assert(graph->has_node(node_id));
        size_t node_length = graph->get_length(graph->get_handle(node_id));
        auto bp = p.second;
        for (auto& pos : bp) {
            pos_t x = pos;
            if (offset(pos) == node_length) continue;
            if (offset(pos) > node_length) {
                cerr << "forwardize_breakpoints error: failure, position " << pos << " is not inside node "
                     << node_id << endl;
                assert(false);
            }
            if (is_rev(pos)) {
                fwd[node_id].insert(reverse(pos, node_length));
            } else {
                fwd[node_id].insert(pos);
            }
        }
    }
    return fwd;
}


// returns breakpoints on the forward strand of the nodes
void find_packed_breakpoints(const Path& path, Packer& packed_breakpoints, bool break_ends,
                             const string& base_quals, double min_baseq, double max_frac_n) {
    // use existing methods to find the breakpoints, then copy them into a packer
    // todo: streamline?
    unordered_map<id_t, set<pos_t>> breakpoints;
    find_breakpoints(path, breakpoints, break_ends, base_quals, min_baseq, max_frac_n);
    breakpoints = forwardize_breakpoints(packed_breakpoints.get_graph(), breakpoints);
    const HandleGraph* graph = packed_breakpoints.get_graph();
    for (auto& id_set : breakpoints) {
        size_t node_len = graph->get_length(graph->get_handle(id_set.first));
        Position position;
        position.set_node_id(id_set.first);
        for (auto pos : id_set.second) {
            size_t offset = get_offset(pos);
            if (offset <= node_len - 1) {
                position.set_offset(offset);
                packed_breakpoints.increment_coverage(packed_breakpoints.position_in_basis(position));
            }
        }
    }
}

unordered_map<id_t, set<pos_t>> filter_breakpoints_by_coverage(const Packer& packed_breakpoints, size_t min_bp_coverage) {
    vector<unordered_map<id_t, set<pos_t>>> bp_maps(get_thread_count());
    size_t n = packed_breakpoints.coverage_size();
    const VectorizableHandleGraph* vec_graph = dynamic_cast<const VectorizableHandleGraph*>(packed_breakpoints.get_graph());
    // we assume our position vector is much larger than the number of filtered breakpoints
    // and scan it in parallel in a first pass
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        if (packed_breakpoints.coverage_at_position(i) >= min_bp_coverage) {
            auto& bp_map = bp_maps[omp_get_thread_num()];
            nid_t node_id = vec_graph->node_at_vector_offset(i+1);
            size_t offset = i - vec_graph->node_vector_offset(node_id);
            bp_map[node_id].insert(make_pos_t(node_id, false, offset));
        }
    }
    // then collect up the breakpoints sequentially in a second pass
    for (size_t i = 1; i < bp_maps.size(); ++i) {
        for (auto& kv : bp_maps[i]) {
            bp_maps[0][kv.first].insert(kv.second.begin(), kv.second.end());
        }
    }
    
    return bp_maps[0];
}
    

path_handle_t add_path_to_graph(MutablePathHandleGraph* graph, const Path& path) {
    path_handle_t path_handle = graph->create_path_handle(path.name(), path.is_circular());
    for (int i = 0; i < path.mapping_size(); ++i) {
        graph->append_step(path_handle, graph->get_handle(path.mapping(i).position().node_id(),
                                                          path.mapping(i).position().is_reverse()));
    }
    return path_handle;
}

map<pos_t, id_t> ensure_breakpoints(MutableHandleGraph* graph,
                                    const unordered_map<id_t, set<pos_t>>& breakpoints) {
    // Set up the map we will fill in with the new node start positions in the
    // old nodes.
    map<pos_t, id_t> toReturn;

    for(auto& kv : breakpoints) {
        // Go through all the nodes we need to break up
        auto original_node_id = kv.first;

        // Save the original node length. We don't want to break here (or later)
        // because that would be off the end.
        id_t original_node_length = graph->get_length(graph->get_handle(original_node_id));

        // We are going through the breakpoints left to right, so we need to
        // keep the node pointer for the right part that still needs further
        // dividing.
        handle_t right_part = graph->get_handle(original_node_id);
        handle_t left_part;

        pos_t last_bp = make_pos_t(original_node_id, false, 0);
        // How far into the original node does our right part start?
        id_t current_offset = 0;

        for(auto breakpoint : kv.second) {
            // For every point at which we need to make a new node, in ascending
            // order (due to the way sets store ints)...

            // ensure that we're on the forward strand (should be the case due to forwardize_breakpoints)
            assert(!is_rev(breakpoint));

            // This breakpoint already exists, because the node starts or ends here
            if(offset(breakpoint) == 0
               || offset(breakpoint) == original_node_length) {
                continue;
            }

            // How far in do we need to break the remaining right part? And how
            // many bases will be in this new left part?
            id_t divide_offset = offset(breakpoint) - current_offset;


#ifdef debug
            cerr << "Need to divide original " << original_node_id << " at " << breakpoint << "/" <<

                original_node_length << endl;
            cerr << "Translates to " << graph->get_id(right_part) << " at " << divide_offset << "/" <<
                graph->get_length(right_part) << endl;
            cerr << "divide offset is " << divide_offset << endl;
#endif

            if (offset(breakpoint) <= 0) { cerr << "breakpoint is " << breakpoint << endl; }
            assert(offset(breakpoint) > 0);
            if (offset(breakpoint) >= original_node_length) { cerr << "breakpoint is " << breakpoint << endl; }
            assert(offset(breakpoint) < original_node_length);

            // Make a new left part and right part. This updates all the
            // existing perfect match paths in the graph.
            std::tie(left_part, right_part) = graph->divide_handle(right_part, divide_offset);

#ifdef debug

            cerr << "Produced " << graph->get_id(left_part) << " (" << graph->get_length(left_part) << " bp)" << endl;
            cerr << "Left " << graph->get_id(right_part) << " (" << graph->get_length(right_part) << " bp)" << endl;
#endif

            // The left part is now done. We know it started at current_offset
            // and ended before breakpoint, so record it by start position.

            // record forward and reverse
            toReturn[last_bp] = graph->get_id(left_part);
            toReturn[reverse(breakpoint, original_node_length)] = graph->get_id(left_part);

            // Record that more sequence has been consumed
            current_offset += divide_offset;
            last_bp = breakpoint;

        }

        // Now the right part is done too. It's going to be the part
        // corresponding to the remainder of the original node.
        toReturn[last_bp] = graph->get_id(right_part);
        toReturn[make_pos_t(original_node_id, true, 0)] = graph->get_id(right_part);

        // and record the start and end of the node
        toReturn[make_pos_t(original_node_id, true, original_node_length)] = 0;
        toReturn[make_pos_t(original_node_id, false, original_node_length)] = 0;

    }

    return toReturn;
}

// We use this function to get the id of the node that contains a position on an
// original node.
static nid_t find_new_node(HandleGraph* graph, pos_t old_pos, const map<pos_t, id_t>& node_translation) {
    if(node_translation.find(make_pos_t(id(old_pos), false, 0)) == node_translation.end()) {
        // The node is unchanged
        return id(old_pos);
    }
    // Otherwise, get the first new node starting after that position, and
    // then look left.
    auto found = node_translation.upper_bound(old_pos);
    assert(found != node_translation.end());
    if (id(found->first) != id(old_pos)
        || is_rev(found->first) != is_rev(old_pos)) {
        return id_t(0);
    }
    // Get the thing before that (last key <= the position we want
    --found;
    assert(graph->has_node(found->second));

    // Return the node we found.
    return found->second;
};


bool simplify_filtered_edits(HandleGraph* graph, Alignment& aln, Path& path, const map<pos_t, id_t>& node_translation,
                             const unordered_map<id_t, size_t>& orig_node_sizes,
                             double min_baseq, double max_frac_n) {

    // check if an edit position is chopped at its next or prev position
    auto is_chopped = [&](pos_t edit_position, bool forward) {
        // todo: better coverage support at node ends (problem is pack structure doesn't have that extra bin)
        bool chopped = offset(edit_position) >= orig_node_sizes.find(id(edit_position))->second - 1 || offset(edit_position) <= 0;
        if (!chopped) {
            if (forward) {
                auto edit_next_position = edit_position;
                ++get_offset(edit_next_position);
                chopped = find_new_node(graph, edit_position, node_translation) != find_new_node(graph, edit_next_position, node_translation);
            } else {
                auto edit_prev_position = edit_position;
                --get_offset(edit_prev_position);
                chopped = find_new_node(graph, edit_position, node_translation) != find_new_node(graph, edit_prev_position, node_translation);
            }
        }
        return chopped;
    };

    bool filtered_an_edit = false;
    bool kept_an_edit = false;

    // The base position in the edit
    size_t position_in_read = 0;

    // stuff that's getting cut out of the read, which requires cuts to
    // quality and and the alignment string
    vector<pair<size_t, size_t>> read_deletions;

    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        Mapping& m = *path.mutable_mapping(i);

        // What node are we on? In old node ID space.
        id_t node_id = m.position().node_id();

        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        pos_t edit_first_position = make_pos_t(m.position());

        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            Edit& e = *m.mutable_edit(j);
            size_t orig_to_length = e.to_length(); // remember here, as we may filter an insertion

            // Work out where its end position on the original node is (inclusive)
            // We don't use this on insertions, so 0-from-length edits don't matter.
            pos_t edit_last_position = edit_first_position;
            get_offset(edit_last_position) += (e.from_length()?e.from_length()-1:0);

            // skip edits whose breakpoitns weren't added due to the coverage filter
            // or edits whose avg base quality fails the min_baseq filter
            if (!edit_is_match(e)) {
                bool chopped;
                if (e.from_length() == 0) {
                    // Just need one-side (prev) test when insertion
                    chopped = is_chopped(edit_first_position, false);
                } else {
                    chopped = is_chopped(edit_first_position, false) && is_chopped(edit_last_position, true);
                }
                if (!chopped || 
                    (min_baseq > 0 && get_avg_baseq(e, aln.quality(), position_in_read) < min_baseq) ||
                    (max_frac_n < 1. && get_fraction_of_ns(e.sequence()) > max_frac_n)) {
                    if (e.from_length() == e.to_length() && !aln.sequence().empty()) {
                        // if we're smoothing a match out, patch the alignment sequence right away
                        // todo: actually look up the correct sequence from the translation.
                    } else if (edit_is_insertion(e)) {
                        // we're trimming off the filtered insertion.  so the alignment's sequence and quality
                        // will need to get updated. 
                        read_deletions.push_back(make_pair(position_in_read, e.to_length()));
                    }
                    e.set_to_length(e.from_length());
                    e.set_sequence("");
                    filtered_an_edit = true;
                    
                } else {
                    kept_an_edit = true;
                }
            }

            // Advance in the right direction along the original node for this edit.
            // This way the next one will start at the right place.
            get_offset(edit_first_position) += e.from_length();

            position_in_read += orig_to_length;
        }
    }

    if (filtered_an_edit) {
        // there's something to simplify
        path = simplify(path);

        if (!read_deletions.empty()) {
            // cut out deleted parts of the read from the sequence and quality
            const string& seq = aln.sequence();
            const string& qual = aln.quality();
            string cut_seq;
            string cut_qual;
            int j = 0;
            for (int i = 0; i < seq.length(); ++i) {
                if (j < read_deletions.size() && i == read_deletions[j].first) {
                    // skip a deleted interval
                    i += read_deletions[j].second - 1;
                    ++j;
                } else {
                    // copy a single position that wasn't skipped
                    cut_seq.push_back(seq[i]);
                    if (!qual.empty()) {
                        cut_qual.push_back(qual[i]);
                    }
               } 
            }
            aln.set_sequence(cut_seq);
            aln.set_quality(cut_qual);
        }
    }

    return kept_an_edit;
}

Path add_nodes_and_edges(MutableHandleGraph* graph,
                         const Path& path,
                         const map<pos_t, id_t>& node_translation,
                         unordered_map<pair<pos_t, string>, vector<id_t>>& added_seqs,
                         unordered_map<id_t, Path>& added_nodes,
                         const unordered_map<id_t, size_t>& orig_node_sizes,
                         size_t max_node_size) {
    
    set<NodeSide> dangling;
    return add_nodes_and_edges(graph,
        path,
        node_translation,
        added_seqs,
        added_nodes,
        orig_node_sizes,
        dangling,
        max_node_size);  

}


Path add_nodes_and_edges(MutableHandleGraph* graph,
                         const Path& path,
                         const map<pos_t, id_t>& node_translation,
                         unordered_map<pair<pos_t, string>, vector<id_t>>& added_seqs,
                         unordered_map<id_t, Path>& added_nodes,
                         const unordered_map<id_t, size_t>& orig_node_sizes,
                         set<NodeSide>& dangling,
                         size_t max_node_size) {
    
    // The basic algorithm is to traverse the path edit by edit, keeping track
    // of a NodeSide for the last piece of sequence we were on. If we hit an
    // edit that creates new sequence, we check if it has been added before If
    // it has, we use it. If not, we create that new sequence as a node, and
    // attach it to the dangling NodeSide(s), and leave its end dangling
    // instead. If we hit an edit that corresponds to a match, we know that
    // there's a breakpoint on each end (since it's bordered by a non-perfect-
    // match or the end of a node), so we can attach its start to the dangling
    // NodeSide(s) and leave its end dangling instead.

    // We need node_translation to translate between node ID space, where the
    // paths are articulated, and new node ID space, where the edges are being
    // made.

    // We need orig_node_sizes so we can remember the sizes of nodes that we
    // modified, so we can interpret the paths. It only holds sizes for modified
    // nodes.
    
    // This is where we will keep the version of the path articulated as
    // actually embedded in the graph.
    Path embedded;
    embedded.set_name(path.name());

    auto create_new_mappings = [&](pos_t p1, pos_t p2, bool is_rev) {
        vector<Mapping> mappings;
        vector<id_t> nodes;
        for (pos_t p = p1; p <= p2; ++get_offset(p)) {
            auto n = find_new_node(graph, p, node_translation);
            assert(n != 0);
            nodes.push_back(find_new_node(graph, p, node_translation));
        }
        auto np = nodes.begin();
        while (np != nodes.end()) {
            size_t c = 0;
            auto n1 = np;
            while (np != nodes.end() && *n1 == *np) {
                ++c;
                ++np; // we'll always increment once
            }
            assert(c);
            // set the mapping position
            Mapping m;
            m.mutable_position()->set_node_id(*n1);
            m.mutable_position()->set_is_reverse(is_rev);
            // and the edit that says we match
            Edit* e = m.add_edit();
            e->set_from_length(c);
            e->set_to_length(c);
            mappings.push_back(m);
        }
        return mappings;
    };

    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        const Mapping& m = path.mapping(i);

        // What node are we on? In old node ID space.
        id_t node_id = m.position().node_id();

        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        pos_t edit_first_position = make_pos_t(m.position());

        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            const Edit& e = m.edit(j);

            // Work out where its end position on the original node is (inclusive)
            // We don't use this on insertions, so 0-from-length edits don't matter.
            pos_t edit_last_position = edit_first_position;
            //get_offset(edit_last_position) += (e.from_length()?e.from_length()-1:0);
            get_offset(edit_last_position) += (e.from_length()?e.from_length()-1:0);

//#define debug_edit true
#ifdef debug_edit
            cerr << "Edit on " << node_id << " from " << edit_first_position << " to " << edit_last_position << endl;
            cerr << pb2json(e) << endl;
#endif

            if(edit_is_insertion(e) || edit_is_sub(e)) {
                // This edit introduces new sequence.
#ifdef debug_edit
                cerr << "Handling ins/sub relative to " << node_id << endl;
#endif
                // store the path representing this novel sequence in the translation table
                auto prev_position = edit_first_position;
                Path from_path;
                auto prev_from_mapping = from_path.add_mapping();
                *prev_from_mapping->mutable_position() = make_position(prev_position);
                auto from_edit = prev_from_mapping->add_edit();
                from_edit->set_sequence(e.sequence());
                from_edit->set_to_length(e.to_length());
                from_edit->set_from_length(e.from_length());
                // find the position after the edit
                // if the edit is not the last in a mapping, the position after is from_length of the edit after this
                pos_t next_position;
                if (j + 1 < m.edit_size()) {
                    next_position = prev_position;
                    get_offset(next_position) += e.from_length();
                    auto next_from_mapping = from_path.add_mapping();
                    *next_from_mapping->mutable_position() = make_position(next_position);
                } else { // implicitly (j + 1 == m.edit_size())
                    // if the edit is the last in a mapping, look at the next mapping position
                    if (i + 1 < path.mapping_size()) {
                        auto& next_mapping = path.mapping(i+1);
                        auto next_from_mapping = from_path.add_mapping();
                        *next_from_mapping->mutable_position() = next_mapping.position();
                    } else {
                        // if we are at the end of the path, then this insertion has no end, and we do nothing
                    }
                }
                // TODO what about forward into reverse????
                if (is_rev(prev_position)) {
                    from_path = simplify(
                        reverse_complement_path(from_path, [&](int64_t id) {
                                auto l = orig_node_sizes.find(id);
                                if (l == orig_node_sizes.end()) {
                                    // The node has no entry, so it must not have been broken
                                    return graph->get_length(graph->get_handle(id));
                                } else {
                                    return l->second;
                                }
                            }));
                }

                // Create the new nodes, reversing it if we are reversed
                vector<id_t> new_nodes;
                pos_t start_pos = make_pos_t(from_path.mapping(0).position());
                // We put in the reverse of our sdequence if we are an insert on
                // the revers of a node, to keep the graph pointing mostly the
                // same direction.
                auto fwd_seq = m.position().is_reverse() ?
                    reverse_complement(e.sequence())
                    : e.sequence();
                auto novel_edit_key = make_pair(start_pos, fwd_seq);
                auto added = added_seqs.find(novel_edit_key);
                if (added != added_seqs.end()) {
                    // if we have the node run already, don't make it again, just use the existing one
                    new_nodes = added->second;
#ifdef debug_edit
                    cerr << "Re-using already added nodes: ";
                    for (auto n : new_nodes) {
                        cerr << n << " ";
                    }
                    cerr << endl;
#endif
                } else {
                    // Make a new run of nodes of up to max_node_size each
                    
                    // Make sure that we are trying to make a run of nodes of
                    // the length we're supposed to be.
                    assert(path_to_length(from_path) == fwd_seq.size());
                    
                    size_t cursor = 0;
                    while (cursor < fwd_seq.size()) {
                        // Until we used up all the sequence, make nodes
                        handle_t new_node = graph->create_handle(fwd_seq.substr(cursor, max_node_size));
                        cursor += max_node_size;
                        
#ifdef debug_edit
                        cerr << "Create new node " << pb2json(*new_node) << endl;
#endif
                        if (!new_nodes.empty()) {
                            // Connect each to the previous node in the chain.
                            graph->create_edge(graph->get_handle(new_nodes.back()), new_node);
#ifdef debug_edit
                            cerr << "Create edge " << new_nodes.back() << "," << graph->get_id(new_node) << endl;
#endif
                        }
                        
                        // Remember the new node
                        new_nodes.push_back(graph->get_id(new_node));
                        
                        // Chop the front of the from path off and associate it
                        // with this node. TODO: this is n^2 in number of nodes
                        // we add because we copy the whole path each time.
                        Path front_path;
                        
                        if (path_to_length(from_path) > graph->get_length(new_node)) {
                            // There will still be path left, so we cut the path
                            tie(front_path, from_path) = cut_path(from_path, graph->get_length(new_node));
                        } else {
                            // We consume the rest of the path. Don't bother cutting it.
                            swap(front_path, from_path);
                        }
                        
                        // The front bit of the path belongs to this new node
                        added_nodes[graph->get_id(new_node)] = front_path;
                        
                    }

                    // reverse the order of the nodes if we did a rev-comp
                    if (m.position().is_reverse()) {
                        std::reverse(new_nodes.begin(), new_nodes.end());
                    }

                    // TODO: fwd_seq can't be empty or problems will be happen
                    // because we'll have an empty vector of created nodes. I
                    // think the edit won't be an insert or sub if it is,
                    // though.

                    // Remember that this run belongs to this edit
                    added_seqs[novel_edit_key] = new_nodes;
                    
                }

                for (auto new_node : new_nodes) {
                    // Add a mapping to each newly created node
                    Mapping& nm = *embedded.add_mapping();
                    nm.mutable_position()->set_node_id(new_node);
                    nm.mutable_position()->set_is_reverse(m.position().is_reverse());

                    // Don't set a rank; since we're going through the input
                    // path in order, the auto-generated ranks will put our
                    // newly created mappings in order.

                    Edit* e = nm.add_edit();
                    size_t l = graph->get_length(graph->get_handle(new_node));
                    e->set_from_length(l);
                    e->set_to_length(l);
                }

                for (auto& dangler : dangling) {
                    // This actually referrs to a node.
                    
                    // Attach what was dangling to the early-in-the-alignment side of the newly created run.
                    auto to_attach = NodeSide(m.position().is_reverse() ? new_nodes.back() : new_nodes.front(),
                        m.position().is_reverse());
                    
#ifdef debug_edit
                    cerr << "Connecting " << dangler << " and " << to_attach << endl;
#endif
                    // Add an edge from the dangling NodeSide to the start of this new node
                    auto from_handle = graph->get_handle(dangler.node, !dangler.is_end);
                    auto to_handle = graph->get_handle(to_attach.node, to_attach.is_end);
                    graph->create_edge(from_handle, to_handle);
                }

                // Dangle the late-in-the-alignment end of this run of new nodes
                dangling.clear();
                dangling.insert(NodeSide(m.position().is_reverse() ? new_nodes.front() : new_nodes.back(),
                    !m.position().is_reverse()));

                // save edit into translated path

            } else if(edit_is_match(e)) {
                // We're using existing sequence

                // We know we have breakpoints on both sides, but we also might
                // have additional breakpoints in the middle. So we need the
                // left node, that contains the first base of the match, and the
                // right node, that contains the last base of the match.
                id_t left_node = find_new_node(graph, edit_first_position, node_translation);
                id_t right_node = find_new_node(graph, edit_last_position, node_translation);

                // TODO: we just assume the outer edges of these nodes are in
                // the right places. They should be if we cut the breakpoints
                // right.

                // get the set of new nodes that we map to
                // and use the lengths of each to create new mappings
                // and append them to the path we are including
                for (auto nm : create_new_mappings(edit_first_position,
                                                   edit_last_position,
                                                   m.position().is_reverse())) {

                    *embedded.add_mapping() = nm;

                    // Don't set a rank; since we're going through the input
                    // path in order, the auto-generated ranks will put our
                    // newly created mappings in order.
                }

#ifdef debug_edit
                cerr << "Handling match relative to " << node_id << endl;
#endif

                for (auto& dangler : dangling) {
#ifdef debug_edit
                    cerr << "Connecting " << dangler << " and " << NodeSide(left_node->id(), m.position().is_reverse()) << endl;
#endif

                    // Connect the left end of the left node we matched in the direction we matched it
                    auto from_handle = graph->get_handle(dangler.node, !dangler.is_end);
                    auto to_handle = graph->get_handle(left_node,  m.position().is_reverse());
                    graph->create_edge(from_handle, to_handle);
                }

                // Dangle the right end of the right node in the direction we matched it.
                if (right_node != 0) {
                    dangling.clear();
                    dangling.insert(NodeSide(right_node, !m.position().is_reverse()));
                }
            } else {
                // We don't need to deal with deletions since we'll deal with the actual match/insert edits on either side
#ifdef debug_edit
                cerr << "Skipping other edit relative to " << node_id << endl;
#endif
            }

            // Advance in the right direction along the original node for this edit.
            // This way the next one will start at the right place.
            get_offset(edit_first_position) += e.from_length();

//#undef debug_edut
        }

    }
    
    // Actually return the embedded path.
    return embedded;

}

vector<Translation> make_translation(const HandleGraph* graph,
                                     const map<pos_t, id_t>& node_translation,
                                     const unordered_map<id_t, Path>& added_nodes,
                                     const unordered_map<id_t, size_t>& orig_node_sizes) {
    vector<Translation> translation;
    // invert the translation
    map<id_t, pos_t> inv_node_trans;
    for (auto& t : node_translation) {
        if (!is_rev(t.first)) {
            inv_node_trans[t.second] = t.first;
        }
    }
    // walk the whole graph
    graph->for_each_handle([&](handle_t handle) {
            id_t node = graph->get_id(handle);
            translation.emplace_back();
            auto& trans = translation.back();
            auto f = inv_node_trans.find(node);
            auto added = added_nodes.find(node);
            if (f != inv_node_trans.end()) {
                // if the node is in the inverted translation, use the position to make a mapping
                auto pos = f->second;
                auto from_mapping = trans.mutable_from()->add_mapping();
                auto to_mapping = trans.mutable_to()->add_mapping();
                // Make sure the to mapping is in the same orientation as the
                // from mapping, since we're going to be making translations on
                // both strands and the new node is the same local orientation
                // as the old node.
                *to_mapping->mutable_position() = make_position(node, is_rev(pos), 0);
                *from_mapping->mutable_position() = make_position(pos);
                auto match_length = graph->get_length(handle);
                auto to_edit = to_mapping->add_edit();
                to_edit->set_to_length(match_length);
                to_edit->set_from_length(match_length);
                auto from_edit = from_mapping->add_edit();
                from_edit->set_to_length(match_length);
                from_edit->set_from_length(match_length);
            } else if (added != added_nodes.end()) {
                // the node is novel
                auto to_mapping = trans.mutable_to()->add_mapping();
                *to_mapping->mutable_position() = make_position(node, false, 0);
                auto to_edit = to_mapping->add_edit();
                to_edit->set_to_length(graph->get_length(handle));
                to_edit->set_from_length(graph->get_length(handle));
                auto from_path = trans.mutable_from();
                *trans.mutable_from() = added->second;
            } else {
                // otherwise we assume that the graph is unchanged
                auto from_mapping = trans.mutable_from()->add_mapping();
                auto to_mapping = trans.mutable_to()->add_mapping();
                *to_mapping->mutable_position() = make_position(node, false, 0);
                *from_mapping->mutable_position() = make_position(node, false, 0);
                auto match_length = graph->get_length(handle);
                auto to_edit = to_mapping->add_edit();
                to_edit->set_to_length(match_length);
                to_edit->set_from_length(match_length);
                auto from_edit = from_mapping->add_edit();
                from_edit->set_to_length(match_length);
                from_edit->set_from_length(match_length);
            }
        });
    std::sort(translation.begin(), translation.end(),
              [&](const Translation& t1, const Translation& t2) {
                  if (!t1.from().mapping_size() && !t2.from().mapping_size()) {
                      // warning: this won't work if we don't have to mappings
                      // this guards against the lurking segfault
                      return t1.to().mapping_size() && t2.to().mapping_size()
                          && make_pos_t(t1.to().mapping(0).position())
                          < make_pos_t(t2.to().mapping(0).position());
                  } else if (!t1.from().mapping_size()) {
                      return true;
                  } else if (!t2.from().mapping_size()) {
                      return false;
                  } else {
                      return make_pos_t(t1.from().mapping(0).position())
                          < make_pos_t(t2.from().mapping(0).position());
                  }
              });
    // append the reverse complement of the translation
    translation.reserve(translation.size() * 2);
    auto get_curr_node_length = [&](id_t id) {
        return graph->get_length(graph->get_handle(id));
    };
    auto get_orig_node_length = [&](id_t id) {
        auto f = orig_node_sizes.find(id);
        if (f == orig_node_sizes.end()) {
            // The node has no entry, so it must not have been broken
            return graph->get_length(graph->get_handle(id));
        }
        return f->second;
    };
    for (auto& trans : translation) {
        translation.emplace_back();
        auto& rev_trans = translation.back();
        *rev_trans.mutable_to() = simplify(reverse_complement_path(trans.to(), get_curr_node_length));
        *rev_trans.mutable_from() = simplify(reverse_complement_path(trans.from(), get_orig_node_length));
    }
    return translation;
}


void add_edges_only(MutableHandleGraph* graph,
                    function<void(function<void(Alignment&)>,bool, bool)> iterate_gam,
                    double min_mapq,
                    size_t min_bp_coverage) {
    // occurrence of each non-graph edge
    // todo: is this too big? do we need something more compact?
    //       novel non-graph edges from read-mappings should be pretty local (and not quadratically exploding)
    //       in general, i think.
    vector<unordered_map<edge_t, size_t>> edge_counts(get_thread_count());

    // scan every non-graph edge in the alignment paths.  if we have a coverage threshold,
    // then fill in the edge_counts map, otherwise just add the edges as soon as they're found
    iterate_gam((function<void(Alignment&)>)[&](Alignment& aln) {

            if (aln.mapping_quality() < min_mapq) {
                return;
            }

            handle_t prev_handle;
            for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
                const Mapping& mapping = aln.path().mapping(i);
                const Position& pos = mapping.position();
                handle_t handle = graph->get_handle(pos.node_id(), pos.is_reverse());
                if (i > 0) {
                    edge_t edge = graph->edge_handle(prev_handle, handle);
                    if (!graph->has_edge(edge)) {
                        if (min_bp_coverage > 1) {
                            edge_counts[omp_get_thread_num()][edge]++;
                        } else {
                            graph->create_edge(edge);
                        }
                    }
                }
                prev_handle = handle;
            }
        }, false, min_bp_coverage > 1);

    if (min_bp_coverage > 1) {
        // second pass required to add edges that meet threshold
        
        // start by merging the thread counters into the first
        for (size_t i = 1; i < edge_counts.size(); ++i) {
            for (const auto& ec : edge_counts[i]) {
                edge_counts[0][ec.first] += ec.second;
            }
            edge_counts[i].clear();
        }
        // then add all the edges that meet threshold
        for (const auto& ec : edge_counts[0]) {
            if (ec.second >= min_bp_coverage) {
                graph->create_edge(ec.first);
            }
        }
    }
}

}
