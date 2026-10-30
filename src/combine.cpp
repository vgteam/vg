/**
 * \file combine.cpp
 *
 * Implements GraphCombiner.
 */

#include "combine.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <handlegraph/algorithms/path_length.hpp>

#include "algorithms/fuse_nodes.hpp"

namespace vg {

using namespace std;
using handlegraph::nid_t;
using handlegraph::offset_t;
using handlegraph::PathMetadata;
using handlegraph::PathSense;
using handlegraph::subrange_t;
using handlegraph::algorithms::get_path_length;

offset_t GraphCombiner::path_start_offset(const subrange_t& subrange,
                                          size_t phase_block) const {
    if (subrange != PathMetadata::NO_SUBRANGE) {
        return subrange.first;
    }
    if (phase_block_is_offset && phase_block != PathMetadata::NO_PHASE_BLOCK) {
        return (offset_t)phase_block;
    }
    return 0;
}

int64_t GraphCombiner::shift_ids_above(MutableHandleGraph* graph, int64_t above) {
    int64_t delta = above - (int64_t)graph->min_node_id();
    if (delta < 0) {
        return 0;
    }
    graph->increment_node_ids(delta + 1);
    return delta + 1;
}

GraphCombiner::GroupKey GraphCombiner::group_key_of(const PathHandleGraph& graph,
                                                    const path_handle_t& path) const {
    // If phase blocks are coordinates, the key leaves them out.
    size_t phase_block = graph.get_phase_block(path);

    return GroupKey{graph.get_sense(path),
                    graph.get_sample_name(path),
                    graph.get_locus_name(path),
                    graph.get_haplotype(path),
                    phase_block_is_offset ? PathMetadata::NO_PHASE_BLOCK : phase_block};
}

string GraphCombiner::path_prefix(const PathHandleGraph& graph, const path_handle_t& path,
                                  size_t length) {
    string sequence;
    for (handle_t h : graph.scan_path(path)) {
        if (sequence.size() >= length) {
            break;
        }
        sequence += graph.get_sequence(h);
    }
    // The last node read can overshoot; a short path just yields all of itself.
    sequence.resize(std::min(sequence.size(), length));
    return sequence;
}

string GraphCombiner::path_suffix(const PathHandleGraph& graph, const path_handle_t& path,
                                  size_t length) {
    if (graph.is_empty(path)) {
        return string();
    }
    string sequence;
    step_handle_t step = graph.path_back(path);
    while (sequence.size() < length) {
        sequence.insert(0, graph.get_sequence(graph.get_handle_of_step(step)));
        if (!graph.has_previous_step(step)) {
            break;
        }
        step = graph.get_previous_step(step);
    }
    // Drop whatever the first node read overshot by, keeping the tail.
    if (sequence.size() > length) {
        sequence.erase(0, sequence.size() - length);
    }
    return sequence;
}

GraphCombiner::GraphCombiner(SeamPolicy seam_policy, bool fuse, bool phase_block_is_offset,
                             const string& log_context)
    : seam_policy(seam_policy), fuse(fuse), phase_block_is_offset(phase_block_is_offset),
      logger(log_context) {
}

unique_ptr<MutablePathDeletableHandleGraph> GraphCombiner::combine(vector<Input>&& inputs) {

    // renumber and shared need at least one input graph; trim needs at least two.
    if (inputs.empty()) {
        throw std::invalid_argument("At least one input graph is required.");
    }
    if (seam_policy == SeamPolicy::COORD_TRIM && inputs.size() < 2) {
        throw std::invalid_argument("Trimming measures overlap between inputs, so it needs "
                                    "at least two input graphs.");
    }

    // Check each input: no path may have multiple fragments within one graph, and
    // every subrange end must match the path's sequence length.
    for (const Input& in : inputs) {
        check_no_split_paths(in);
        check_subrange_matches_length(in);
    }

    if (seam_policy != SeamPolicy::COORD_TRIM) {
        // Non-trim modes just ingest each input in turn.
        for (Input& in : inputs) {
            ingest(in);
        }
    } else {
        // Describe each input's reference path, and sort by base path then start offset.
        vector<Chunk> chunks;
        chunks.reserve(inputs.size());
        for (Input& in : inputs) {
            ReferenceDescription ref = describe_reference(in);
            chunks.push_back(Chunk{std::move(in), std::move(ref)});
        }
        sort_by_base_path_offset(chunks);
        check_overlaps_agree(chunks);

        for (size_t i = 0; i < chunks.size(); ++i) {
            Chunk& chunk = chunks[i];

            // The sort made each base path's chunks contiguous, so two bools tell
            // whether this is that base path's first or last chunk.
            auto found = seams.find(chunk.ref.base_path_key);
            bool is_first_of_ref = (found == seams.end());
            bool is_last_of_ref = (i + 1 == chunks.size())
                || chunks[i + 1].ref.base_path_key != chunk.ref.base_path_key;

            if (!is_last_of_ref) {
                check_all_paths_end_at_reference_last_node(chunk);
            }
            if (!is_first_of_ref) {
                check_all_paths_start_at_reference_first_node(chunk);
            }

            // Handles will be invalidated, so record IDs and orientations first.
            nid_t new_start_id = chunk.graph->get_id(chunk.ref.start_handle);
            bool start_is_reverse = chunk.graph->get_is_reverse(chunk.ref.start_handle);
            nid_t ref_end_id = chunk.graph->get_id(chunk.ref.end_handle);
            bool end_is_reverse = chunk.graph->get_is_reverse(chunk.ref.end_handle);

            // offset_t is unsigned, so check for a gap before subtracting; then
            // trim_and_shift cuts the overlap away and gives the new start ID.
            if (!is_first_of_ref) {
                const SeamState& running = found->second;

                if (running.is_circular != chunk.ref.is_circular) {
                    stringstream message;
                    message << "Reference path \"" << chunk.ref.path_name << "\" in \""
                            << chunk.name << "\" is "
                            << (chunk.ref.is_circular ? "circular" : "not circular")
                            << ", but its base path in \"" << running.last_input_name
                            << "\" is " << (running.is_circular ? "circular" : "not circular")
                            << ". Pieces of one path have to agree on whether it is circular.";
                    throw InputConflictError(message.str());
                }
                if (running.end_offset < chunk.ref.start_offset) {
                    stringstream message;
                    message << "Base path offset gap of "
                            << (chunk.ref.start_offset - running.end_offset)
                            << " bp between \"" << running.last_input_name
                            << "\" (its reference path ends at offset " << running.end_offset
                            << ") and \"" << chunk.name << "\" (its reference path starts at offset "
                            << chunk.ref.start_offset
                            << "). Trimming does not fill in gaps.";
                    throw ChunkOverlapError(message.str());
                }
                offset_t overlap = running.end_offset - chunk.ref.start_offset;
                size_t left_end_node_len = dest->get_length(
                    dest->get_handle(running.end_id, running.end_is_reverse));
                new_start_id = trim_and_shift(chunk, overlap, left_end_node_len);
            }

            // ingest() shifts the chunk's node IDs, so look up its end nodes in dest by
            // their shifted IDs, then join it onto the end of the base path built so far.
            int64_t shift = ingest(chunk);

            handle_t chunk_ref_end = dest->get_handle((nid_t)((int64_t)ref_end_id + shift),
                                                     end_is_reverse);
            handle_t new_end = chunk_ref_end;
            if (!is_first_of_ref) {
                handle_t new_start = dest->get_handle((nid_t)((int64_t)new_start_id + shift),
                                                     start_is_reverse);
                new_end = connect_seam(found->second, new_start, chunk_ref_end, chunk.name);
            }

            // Update the SeamState for this chunk.
            seams[chunk.ref.base_path_key] = SeamState{dest->get_id(new_end),
                                                       dest->get_is_reverse(new_end),
                                                       chunk.ref.end_offset,
                                                       chunk.name,
                                                       chunk.ref.is_circular};
        }
    }
    inputs.clear();

    merge_path_fragments();
    return std::move(dest);
}

void GraphCombiner::check_no_split_paths(const Input& input) const {

    // Group by GroupKey, recording each path's start offset and name.
    map<GroupKey, vector<pair<offset_t, string>>> groups;
    input.graph->for_each_path_handle([&](const path_handle_t& path) {
        groups[group_key_of(*input.graph, path)]
            .emplace_back(path_start_offset(input.graph->get_subrange(path),
                                            input.graph->get_phase_block(path)),
                          input.graph->get_path_name(path));
    });

    for (auto& group : groups) {
        vector<pair<offset_t, string>>& pieces = group.second;
        if (pieces.size() <= 1) {
            continue;
        }

        // Sort by offset then name, and build the list of names for the message.
        std::sort(pieces.begin(), pieces.end());
        string names;
        for (const auto& piece : pieces) {
            names += (names.empty() ? "\"" : ", \"") + piece.second + "\"";
        }
        stringstream message;
        message << "\"" << input.name << "\" holds " << pieces.size()
                << " pieces of one path (" << names << "). Only pieces that are in "
                << "separate inputs are joined; a path already split inside one input is "
                << "not reassembled.";
        throw InputConflictError(message.str());
    }
}

void GraphCombiner::check_subrange_matches_length(const Input& input) {
    input.graph->for_each_path_handle([&](const path_handle_t& path) {
        subrange_t subrange = input.graph->get_subrange(path);
        if (subrange == PathMetadata::NO_SUBRANGE
            || subrange.second == PathMetadata::NO_END_POSITION) {
            return;
        }
        offset_t actual_end = subrange.first + (offset_t)get_path_length(input.graph.get(), path);
        if (subrange.second != actual_end) {
            stringstream message;
            message << "Path \"" << input.graph->get_path_name(path) << "\" in \"" << input.name
                    << "\" claims to end at offset " << subrange.second
                    << ", but its sequence ends at offset " << actual_end << ".";
            throw PathLengthMismatchError(message.str());
        }
    });
}

GraphCombiner::ReferenceDescription GraphCombiner::describe_reference(const Input& input) const {

    // Collect every REFERENCE-sense path into refs.
    vector<path_handle_t> refs;
    input.graph->for_each_path_handle([&](const path_handle_t& p) {
        if (input.graph->get_sense(p) == PathSense::REFERENCE) {
            refs.push_back(p);
        }
    });

    // Two checks: there must be exactly one REFERENCE-sense path.
    if (refs.empty()) {
        stringstream message;
        message << "No REFERENCE-sense path found in \"" << input.name
                << "\"; trimming reads each input's position off its REFERENCE-sense path, "
                << "and needs exactly one per graph.";
        throw ReferencePathError(message.str());
    }
    if (refs.size() > 1) {
        stringstream message;
        message << refs.size() << " REFERENCE-sense paths found in \"" << input.name
                << "\"; trimming needs exactly one, so it knows which one gives the "
                << "input its position.";
        throw ReferencePathError(message.str());
    }
    path_handle_t path = refs[0];

    // Record the reference path's name, the base path it is part of, and where
    // it starts and ends.
    ReferenceDescription ref;
    ref.path_name = input.graph->get_path_name(path);
    ref.base_path_key = BasePathKey{input.graph->get_sample_name(path),
                                    input.graph->get_locus_name(path),
                                    input.graph->get_haplotype(path)};
    ref.is_circular = input.graph->get_is_circular(path);
    ref.start_offset = path_start_offset(input.graph->get_subrange(path),
                                         input.graph->get_phase_block(path));
    ref.end_offset = ref.start_offset + (offset_t)get_path_length(input.graph.get(), path);

    ref.start_handle = input.graph->get_handle_of_step(input.graph->path_begin(path));
    ref.end_handle = input.graph->get_handle_of_step(input.graph->path_back(path));

    // Check boundary orientation: the reference path must visit its first and last
    // nodes forward.
    if (input.graph->get_is_reverse(ref.start_handle)
        || input.graph->get_is_reverse(ref.end_handle)) {
        stringstream message;
        message << "Reference path \"" << ref.path_name << "\" in \"" << input.name
                << "\" visits its boundary node in reverse orientation; trimming only "
                << "supports forward-strand boundaries.";
        throw ReferencePathError(message.str());
    }
    return ref;
}

void GraphCombiner::sort_by_base_path_offset(vector<Chunk>& chunks) {

    // Different base paths keep first-seen order; chunks of one base path sort by
    // ascending offset.
    map<BasePathKey, size_t> first_seen;
    for (const Chunk& chunk : chunks) {
        first_seen.emplace(chunk.ref.base_path_key, first_seen.size());
    }
    std::stable_sort(chunks.begin(), chunks.end(),
                     [&](const Chunk& a, const Chunk& b) {
                         size_t ra = first_seen.at(a.ref.base_path_key);
                         size_t rb = first_seen.at(b.ref.base_path_key);
                         if (ra != rb) {
                             return ra < rb;
                         }
                         return a.ref.start_offset < b.ref.start_offset;
                     });
}

void GraphCombiner::check_overlaps_agree(const vector<Chunk>& chunks) const {

    // The sort made each base path's chunks contiguous and offset-ascending, so only
    // neighbours can overlap.
    for (size_t i = 1; i < chunks.size(); ++i) {
        const Chunk& left = chunks[i - 1];
        const Chunk& right = chunks[i];

        // Skip a different base path, and a gap, which combine() reports itself.
        if (right.ref.base_path_key != left.ref.base_path_key
            || left.ref.end_offset <= right.ref.start_offset) {
            continue;
        }
        offset_t overlap = left.ref.end_offset - right.ref.start_offset;

        // Read the shared span off each side. A nested chunk reads short; the checks
        // in trim_and_shift catch that.
        string from_left = path_suffix(*left.graph,
                                       left.graph->get_path_handle(left.ref.path_name),
                                       overlap);
        string from_right = path_prefix(*right.graph,
                                        right.graph->get_path_handle(right.ref.path_name),
                                        overlap);

        // Find the first base they disagree on, if any.
        size_t shared = std::min(from_left.size(), from_right.size());
        size_t at = 0;
        while (at < shared && from_left[at] == from_right[at]) {
            ++at;
        }
        if (at < shared) {
            stringstream message;
            message << "Chunks \"" << left.name << "\" and \"" << right.name
                    << "\" overlap on base path coordinates ["
                    << right.ref.start_offset << ", " << left.ref.end_offset
                    << ") but disagree there: at offset "
                    << (right.ref.start_offset + at) << " \"" << left.name
                    << "\" has " << from_left[at] << " and \"" << right.name
                    << "\" has " << from_right[at]
                    << ". Trimming drops the right chunk's copy of the overlap, "
                    << "so the two have to spell the same sequence.";
            throw ChunkOverlapError(message.str());
        }
    }
}

void GraphCombiner::check_all_paths_start_at_reference_first_node(const Chunk& chunk) const {
    chunk.graph->for_each_path_handle([&](const path_handle_t& p) {
        handle_t first = chunk.graph->get_handle_of_step(chunk.graph->path_begin(p));
        if (first != chunk.ref.start_handle) {
            stringstream message;
            message << "Path \"" << chunk.graph->get_path_name(p) << "\" in \""
                    << chunk.name << "\" does not start at the reference path's start node "
                    << "(id " << chunk.graph->get_id(chunk.ref.start_handle)
                    << "). Trimming cuts and joins at that one node, so every path in a "
                    << "chunk that has a chunk to its left has to begin there.";
            throw ChunkOverlapError(message.str());
        }
    });
}

void GraphCombiner::check_all_paths_end_at_reference_last_node(const Chunk& chunk) const {
    chunk.graph->for_each_path_handle([&](const path_handle_t& p) {
        handle_t last = chunk.graph->get_handle_of_step(chunk.graph->path_back(p));
        if (last != chunk.ref.end_handle) {
            stringstream message;
            message << "Path \"" << chunk.graph->get_path_name(p) << "\" in \""
                    << chunk.name << "\" does not end at the reference path's end node "
                    << "(id " << chunk.graph->get_id(chunk.ref.end_handle)
                    << "). Trimming cuts and joins at that one node, so every path in a "
                    << "chunk that has a chunk to its right has to end there.";
            throw ChunkOverlapError(message.str());
        }
    });
}

nid_t GraphCombiner::trim_and_shift(Chunk& chunk, offset_t overlap,
                                    size_t left_end_node_len) {

    // A raw pointer to this chunk's graph, and the ID of its reference path's start node.
    MutablePathDeletableHandleGraph* graph = chunk.graph.get();
    nid_t new_start_id = graph->get_id(chunk.ref.start_handle);

    // Two checks: the overlap must be strictly shorter than the first node.
    if (overlap > 0) {
        size_t start_len = graph->get_length(chunk.ref.start_handle);
        if (overlap > start_len) {
            stringstream message;
            message << "Overlap of " << overlap << " bp before \""
                    << chunk.name << "\" exceeds the length of its first node ("
                    << start_len << " bp on node " << new_start_id
                    << "), so the overlap is not all on that node. Re-chunk so it is.";
            throw ChunkOverlapError(message.str());
        }
        if (overlap == start_len) {
            stringstream message;
            message << "Overlap of " << overlap
                    << " bp exactly equals the length of the first node of \""
                    << chunk.name << "\"; trimming it would leave a zero-length "
                    << "node. Re-chunk so that first node extends past the overlap.";
            throw ChunkOverlapError(message.str());
        }

        // Split the first node in two: the overlapping piece, and the remainder,
        // which becomes the new start.
        vector<handle_t> pieces = graph->divide_handle(chunk.ref.start_handle,
                                                       vector<size_t>{overlap});
        handle_t overlapping = pieces[0];
        handle_t remainder = pieces[1];
        new_start_id = graph->get_id(remainder);

        // Drop every path's first step, so each now starts at the remainder.
        vector<path_handle_t> all_paths;
        graph->for_each_path_handle([&](const path_handle_t& p) {
            all_paths.push_back(p);
        });
        for (path_handle_t p : all_paths) {
            graph->pop_front_step(p);
        }

        // Check that no edge enters overlapping from the left; if none does,
        // destroy it.
        size_t inbound = 0;
        graph->follow_edges(overlapping, true, [&](const handle_t&) { inbound++; });
        if (inbound > 0) {
            stringstream message;
            message << "The first " << overlap << " bp of \"" << chunk.name
                    << "\" overlap the chunk to its left and have to be trimmed, "
                    << "but " << inbound << " edge(s) attach to node "
                    << graph->get_id(overlapping)
                    << " there and would be lost. Trimming needs the chunk "
                    << "boundary to be a clean left edge.";
            throw ChunkOverlapError(message.str());
        }
        graph->destroy_handle(overlapping);
    }

    // Compute the metadata shift.
    int64_t meta_shift = (int64_t)overlap - (fuse ? (int64_t)left_end_node_len : 0);
    if (meta_shift == 0) {
        return new_start_id;
    }

    // Collect every path's metadata and all steps.
    struct PathSnapshot {
        PathSense sense;         ///< Sense of the path.
        string sample;           ///< Sample name.
        string locus;            ///< Locus name.
        size_t haplotype;        ///< Haplotype number.
        size_t phase_block;      ///< Phase block.
        subrange_t subrange;     ///< Coordinate range.
        bool is_circular;        ///< Whether it is circular.
        vector<handle_t> steps;  ///< Nodes visited, in order.
    };

    vector<PathSnapshot> snapshots;
    vector<path_handle_t> to_destroy;
    graph->for_each_path_handle([&](const path_handle_t& p) {
        PathSnapshot s;
        s.sense = graph->get_sense(p);
        s.sample = graph->get_sample_name(p);
        s.locus = graph->get_locus_name(p);
        s.haplotype = graph->get_haplotype(p);
        s.phase_block = graph->get_phase_block(p);
        s.subrange = graph->get_subrange(p);
        s.is_circular = graph->get_is_circular(p);
        // A path without coordinates starts at 0.
        int64_t start = 0;
        if (s.subrange != PathMetadata::NO_SUBRANGE) {
            start = (int64_t)s.subrange.first;
        } else if (phase_block_is_offset && s.phase_block != PathMetadata::NO_PHASE_BLOCK) {
            start = (int64_t)s.phase_block;
        }
        if (start + meta_shift < 0) {
            stringstream message;
            message << "Path \"" << graph->get_path_name(p) << "\" in \"" << chunk.name
                    << "\" starts at offset " << start << ", but joining it on moves it back by "
                    << -meta_shift << " bp, which would put its start before 0. Its "
                    << "coordinates do not match the chunk to its left.";
            throw ChunkOverlapError(message.str());
        }
        for (handle_t h : graph->scan_path(p)) {
            s.steps.push_back(h);
        }
        snapshots.push_back(std::move(s));
        to_destroy.push_back(p);
    });
    for (path_handle_t p : to_destroy) {
        graph->destroy_path(p);
    }
    for (PathSnapshot& s : snapshots) {
        // Shift each path's coordinates by meta_shift, and rebuild the path.
        if (s.subrange != PathMetadata::NO_SUBRANGE) {
            s.subrange.first = (offset_t)((int64_t)s.subrange.first + meta_shift);
        } else if (phase_block_is_offset && s.phase_block != PathMetadata::NO_PHASE_BLOCK) {
            s.phase_block = (size_t)((int64_t)s.phase_block + meta_shift);
        } else {
            s.subrange = subrange_t{(offset_t)meta_shift, PathMetadata::NO_END_POSITION};
        }
        path_handle_t rebuilt = graph->create_path(s.sense, s.sample, s.locus, s.haplotype,
                                                  s.phase_block, s.subrange, s.is_circular);
        for (handle_t h : s.steps) {
            graph->append_step(rebuilt, h);
        }
    }
    return new_start_id;
}

int64_t GraphCombiner::ingest(Input& input) {

    if (!dest) {
        // The first graph becomes dest: nothing to copy or renumber.
        dest = std::move(input.graph);
        max_node_id = dest->max_node_id();
        return 0;
    }

    int64_t shift = 0;
    if (seam_policy == SeamPolicy::SHARED_IDS) {
        share_nodes_by_id(input);
        max_node_id = std::max(max_node_id, (int64_t)input.graph->max_node_id());
    } else {
        shift = shift_ids_above(input.graph.get(), max_node_id);
        max_node_id = std::max(max_node_id, (int64_t)input.graph->max_node_id());
        handlealgs::copy_handle_graph(input.graph.get(), dest.get());
    }
    copy_paths_checked(input);
    // Free the source rather than holding every input until combine finishes.
    input.graph.reset();
    return shift;
}

void GraphCombiner::share_nodes_by_id(Input& input) {

    input.graph->for_each_handle([&](const handle_t& source) {
        nid_t id = input.graph->get_id(source);
        string sequence = input.graph->get_sequence(source);
        if (!dest->has_node(id)) {
            dest->create_handle(sequence, id);
            return;
        }
        // Under this policy a repeated ID is the same node, so two sequences for
        // one ID is a contradiction with no way to keep both.
        if (dest->get_sequence(dest->get_handle(id)) != sequence) {
            stringstream message;
            message << "Node " << id << " has different sequences in \"" << input.name
                    << "\" and in an earlier input. When node IDs are shared, a repeated "
                    << "node ID means the same node, so the sequences have to match.";
            throw InputConflictError(message.str());
        }
    });
    input.graph->for_each_edge([&](const edge_t& source) {
        handle_t from = dest->get_handle(input.graph->get_id(source.first),
                                        input.graph->get_is_reverse(source.first));
        handle_t to = dest->get_handle(input.graph->get_id(source.second),
                                       input.graph->get_is_reverse(source.second));
        dest->create_edge(from, to);
    });
}

void GraphCombiner::copy_paths_checked(Input& input) {

    input.graph->for_each_path_handle([&](const path_handle_t& source) {
        string path_name = input.graph->get_path_name(source);
        if (dest->has_path(path_name)) {
            stringstream message;
            message << "Path \"" << path_name << "\" is in \"" << input.name
                    << "\" and in an earlier input. Pieces of one path have to be "
                    << "told apart by name, so give them distinct subranges.";
            throw InputConflictError(message.str());
        }
        handlealgs::copy_path(input.graph.get(), source, dest.get());
    });
}

handle_t GraphCombiner::connect_seam(const SeamState& base_path_end, handle_t chunk_start,
                                     handle_t chunk_end, const string& input_name) {

    handle_t left = dest->get_handle(base_path_end.end_id, base_path_end.end_is_reverse);

    // If not fusing nodes, just add an edge.
    if (!fuse) {
        dest->create_edge(left, chunk_start);
        return chunk_end;
    }

    // can_fuse_nodes decides whether the two nodes can be fused.
    if (!algorithms::can_fuse_nodes(*dest, left, chunk_start)) {
        // Fusing removes both facing ends, and something is attached to one.
        stringstream message;
        message << "Cannot fuse node " << dest->get_id(left) << " from \""
                << base_path_end.last_input_name << "\" to node " << dest->get_id(chunk_start)
                << " from \"" << input_name << "\": one of them still has an edge "
                << "on the side being fused, which fusing would have to drop.";
        throw FuseError(message.str());
    }

    // Note the end node before fusing, fuse into a new node, and update the max ID.
    bool ref_ends_at_seam = (dest->get_id(chunk_end) == dest->get_id(chunk_start));
    handle_t fused = algorithms::fuse_nodes(dest.get(), left, chunk_start);
    max_node_id = std::max(max_node_id, (int64_t)dest->get_id(fused));
    return ref_ends_at_seam ? fused : chunk_end;
}

void GraphCombiner::merge_path_fragments() {

    struct Fragment {
        path_handle_t handle;            ///< Handle to the path.
        subrange_t subrange;             ///< Coordinate range.
        size_t phase_block;              ///< Phase block.
        size_t length;                   ///< Sequence length.
        bool positioned_by_phase_block;  ///< Whether the phase block gives its position.
        bool is_circular;                ///< Whether it is circular.
    };
    map<GroupKey, vector<Fragment>> groups;

    // Walk every path in dest, filling in a Fragment for each.
    dest->for_each_path_handle([&](const path_handle_t& path) {
        Fragment f;
        f.handle = path;
        f.subrange = dest->get_subrange(path);
        f.phase_block = dest->get_phase_block(path);
        f.length = get_path_length(dest.get(), path);
        f.positioned_by_phase_block = phase_block_is_offset
            && f.subrange == PathMetadata::NO_SUBRANGE
            && f.phase_block != PathMetadata::NO_PHASE_BLOCK;
        f.is_circular = dest->get_is_circular(path);
        groups[group_key_of(*dest, path)].push_back(std::move(f));
    });

    // Lambdas computing a fragment's start and end coordinates.
    auto fragment_start = [&](const Fragment& f) -> offset_t {
        return path_start_offset(f.subrange, f.phase_block);
    };
    auto fragment_end = [&](const Fragment& f) -> offset_t {
        // Fusing lengthens the left fragment's last node, so its measured length can
        // overshoot; the subrange end, checked against the input, is the real end.
        if (f.subrange != PathMetadata::NO_SUBRANGE
            && f.subrange.second != PathMetadata::NO_END_POSITION) {
            return f.subrange.second;
        }
        return fragment_start(f) + (offset_t)f.length;
    };

    // Merge each group in turn.
    for (auto& group : groups) {
        vector<Fragment>& frags = group.second;
        if (frags.size() <= 1) {
            // Nothing to merge, and nothing to rename: a lone path is untouched.
            continue;
        }
        PathSense sense = std::get<0>(group.first);
        const string& sample = std::get<1>(group.first);
        const string& locus = std::get<2>(group.first);
        size_t haplotype = std::get<3>(group.first);

        // Circularity is not part of the key, so the fragments have to agree on it.
        bool is_circular = frags.front().is_circular;
        for (const Fragment& f : frags) {
            if (f.is_circular != is_circular) {
                stringstream message;
                message << "Fragments \"" << dest->get_path_name(frags.front().handle)
                        << "\" and \"" << dest->get_path_name(f.handle)
                        << "\" are pieces of one path, but only one of them is circular.";
                throw InputConflictError(message.str());
            }
        }

        std::sort(frags.begin(), frags.end(),
                  [&](const Fragment& a, const Fragment& b) {
                      return fragment_start(a) < fragment_start(b);
                  });

        // Detect any gap, updating the furthest offset along the base path reached so far.
        offset_t covered_end = fragment_end(frags.front());
        for (size_t i = 1; i < frags.size(); ++i) {
            offset_t next_start = fragment_start(frags[i]);
            if (covered_end < next_start) {
                stringstream message;
                message << "Gap of " << (next_start - covered_end)
                        << " bp before fragment \""
                        << dest->get_path_name(frags[i].handle)
                        << "\"; the fragments before it only reach offset "
                        << covered_end << ". Gaps between fragments are not filled in.";
                throw FragmentMergeError(message.str());
            }
            covered_end = std::max(covered_end, fragment_end(frags[i]));
        }

        // With phase blocks read as coordinates, check they run continuously.
        for (size_t i = 0; phase_block_is_offset && i + 1 < frags.size(); ++i) {
            if (!frags[i].positioned_by_phase_block
                || !frags[i + 1].positioned_by_phase_block) {
                continue;
            }
            offset_t expected = (offset_t)frags[i].phase_block + (offset_t)frags[i].length;
            if (expected != (offset_t)frags[i + 1].phase_block) {
                stringstream message;
                message << "Phase blocks are read as offsets, but phase blocks of \""
                        << dest->get_path_name(frags[i].handle)
                        << "\" do not line up as coordinates (block "
                        << frags[i].phase_block << " + length " << frags[i].length
                        << " != block " << frags[i + 1].phase_block
                        << "). Are these phase block numbers rather than offsets?";
                throw FragmentMergeError(message.str());
            }
        }

        // Compute the merged subrange, phase block, and name.
        offset_t combined_start = fragment_start(frags.front());
        offset_t combined_end = 0;
        for (const Fragment& f : frags) {
            combined_end = std::max(combined_end, fragment_end(f));
        }
        subrange_t combined_subrange{combined_start, combined_end};

        // HAPLOTYPE requires a phase block. Keep the group's own; treating phase blocks
        // as coordinates added it to the subrange, so 0 is used there.
        size_t group_phase_block = std::get<4>(group.first);
        size_t merged_phase_block = PathMetadata::NO_PHASE_BLOCK;
        if (sense == PathSense::HAPLOTYPE) {
            merged_phase_block = (group_phase_block == PathMetadata::NO_PHASE_BLOCK)
                ? 0 : group_phase_block;
        }
        string merged_name = PathMetadata::create_path_name(sense, sample, locus, haplotype,
                                                           merged_phase_block,
                                                           combined_subrange);

        // Can't lean on create_path(): PackedGraph throws on a duplicate, HashGraph
        // leaves the old path unreachable by name, and VG just appends to it. So check here.
        if (dest->has_path(merged_name)) {
            bool is_own_fragment = false;
            for (const Fragment& f : frags) {
                if (dest->get_path_name(f.handle) == merged_name) {
                    is_own_fragment = true;
                    break;
                }
            }
            if (!is_own_fragment) {
                // A path with this name has the same GroupKey, so it is one of frags.
                throw std::logic_error("merged path name collides with a path outside its group");
            }
        }

        // Save each fragment's node steps into routes.
        vector<vector<handle_t>> routes;
        routes.reserve(frags.size());
        for (const Fragment& f : frags) {
            vector<handle_t> steps;
            for (handle_t h : dest->scan_path(f.handle)) {
                steps.push_back(h);
            }
            routes.push_back(std::move(steps));
        }

        // Destroy the old fragments, and create the merged path.
        for (const Fragment& f : frags) {
            dest->destroy_path(f.handle);
        }

        path_handle_t merged = dest->create_path(sense, sample, locus, haplotype,
                                                merged_phase_block, combined_subrange,
                                                is_circular);
        const vector<handle_t>* previous = nullptr;
        // Walk each fragment's steps in order, appending them to the merged path.
        for (const vector<handle_t>& steps : routes) {
            if (steps.empty()) {
                continue;
            }
            // `skip` is how many steps at the start of `steps` repeat the end of `previous`.
            // Try each length `run`, longest first, down to 1 (a run of 0 would always match);
            // `k` goes through the run, comparing step k of `steps` with step k of `previous`'s last `run` steps.
            size_t skip = 0;
            if (previous != nullptr) {
                size_t most = std::min(previous->size(), steps.size());
                for (size_t run = most; run >= 1; --run) {
                    bool same = true;
                    for (size_t k = 0; k < run; ++k) {
                        if ((*previous)[previous->size() - run + k] != steps[k]) {
                            same = false;
                            break;
                        }
                    }
                    if (same) {
                        skip = run;
                        break;
                    }
                }
                if (skip == 0) {
                    // Nothing shared, so they abut and the graph must say so.
                    dest->create_edge(previous->back(), steps.front());
                }
            }
            for (size_t k = skip; k < steps.size(); ++k) {
                dest->append_step(merged, steps[k]);
            }
            previous = &steps;
        }

        logger.info() << "Merged " << frags.size() << " fragments into \"" << merged_name
                      << "\"" << endl;
    }
}

}
