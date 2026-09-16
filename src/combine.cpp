/**
 * \file combine.cpp
 *
 * Implements GraphCombiner.
 */

#include "combine.hpp"

#include <algorithm>

#include "algorithms/fuse_nodes.hpp"

namespace vg {

using namespace std;
using handlegraph::nid_t;
using handlegraph::offset_t;
using handlegraph::PathMetadata;
using handlegraph::PathSense;
using handlegraph::subrange_t;

/// Compute the length of a path.
size_t GraphCombiner::compute_path_length(const PathHandleGraph* graph,
                                          const path_handle_t& path) {
    // A path-position graph already knows this.
    if (const PathPositionHandleGraph* positional
        = dynamic_cast<const PathPositionHandleGraph*>(graph)) {
        return positional->get_path_length(path);
    }
    size_t length = 0;
    for (handle_t h : graph->scan_path(path)) {
        length += graph->get_length(h);
    }
    return length;
}

/// Return a path's start offset.
offset_t GraphCombiner::path_start_offset(const subrange_t& subrange, size_t phase_block,
                                          bool phase_block_is_offset) {
    if (subrange != PathMetadata::NO_SUBRANGE) {
        return subrange.first;
    }
    if (phase_block_is_offset && phase_block != PathMetadata::NO_PHASE_BLOCK) {
        return (offset_t)phase_block;
    }
    return 0;
}

/// Raise every node ID above the previous graph's max ID; return the shift.
int64_t GraphCombiner::renumber_out_of_way(MutableHandleGraph* graph, int64_t above) {
    int64_t delta = above - (int64_t)graph->min_node_id();
    if (delta < 0) {
        return 0;
    }
    graph->increment_node_ids(delta + 1);
    return delta + 1;
}

GraphCombiner::GroupKey GraphCombiner::group_key_of(const PathHandleGraph& graph,
                                                    const path_handle_t& path,
                                                    bool phase_block_is_offset) {
    // If phase blocks are coordinates, the key folds them out.
    size_t phase_block = graph.get_phase_block(path);

    return GroupKey{graph.get_sense(path),
                    graph.get_sample_name(path),
                    graph.get_locus_name(path),
                    graph.get_haplotype(path),
                    graph.get_is_circular(path),
                    phase_block_is_offset ? PathMetadata::NO_PHASE_BLOCK : phase_block};
}

/// Read up to `length` bp off the front of a path.
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

/// Read up to `length` bp off the back of a path.
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

/// Constructor.
GraphCombiner::GraphCombiner(SeamPolicy seam, bool fuse, bool phase_block_is_offset)
    : seam(seam), fuse(fuse), phase_block_is_offset(phase_block_is_offset),
      logger("vg combine") {
}

/// Combine the input graphs.
unique_ptr<MutablePathDeletableHandleGraph> GraphCombiner::combine(vector<Input>&& inputs) {

    // renumber and shared need at least one input graph; trim needs at least two.
    if (inputs.empty()) {
        logger.error() << "At least one input graph is required." << endl;
    }
    if (seam == SeamPolicy::COORD_TRIM && inputs.size() < 2) {
        logger.error() << "--seam trim measures overlap between chunks, so it needs "
                       << "at least two input graphs." << endl;
    }

    // Turn each public Input into a private Chunk.
    vector<Chunk> chunks;
    chunks.reserve(inputs.size());
    for (Input& in : inputs) {
        Chunk chunk;
        chunk.name = std::move(in.name);
        chunk.graph = std::move(in.graph);
        chunks.push_back(std::move(chunk));
    }
    inputs.clear();

    // Check each input: no path may have multiple fragments within one graph.
    for (const Chunk& chunk : chunks) {
        check_no_split_paths(chunk);
    }

    // For COORD_TRIM, read each chunk's REFERENCE identity and sort by reference
    // then start offset.
    if (seam == SeamPolicy::COORD_TRIM) {
        for (Chunk& chunk : chunks) {
            describe_reference(chunk);
        }
        sort_by_reference_offset(chunks);
        check_overlaps_agree(chunks);
    }

    for (size_t i = 0; i < chunks.size(); ++i) {
        Chunk& chunk = chunks[i];

        // Non-trim modes just ingest each chunk in turn.
        if (seam != SeamPolicy::COORD_TRIM) {
            ingest(chunk);
            continue;
        }

        // The sort made each reference's chunks contiguous, so two bools settle
        // whether this is that reference's first or last chunk.
        auto found = seams.find(chunk.ref_key);
        bool is_first_of_ref = (found == seams.end());
        bool is_last_of_ref = (i + 1 == chunks.size())
            || chunks[i + 1].ref_key != chunk.ref_key;

        if (!is_last_of_ref) {
            check_all_paths_end_at(chunk);
        }
        if (!is_first_of_ref) {
            check_all_paths_start_at(chunk);
        }

        // Handles will be invalidated, so record IDs and orientations first.
        nid_t new_start_id = chunk.graph->get_id(chunk.ref_start_handle);
        bool start_is_reverse = chunk.graph->get_is_reverse(chunk.ref_start_handle);
        nid_t ref_end_id = chunk.graph->get_id(chunk.ref_end_handle);
        bool end_is_reverse = chunk.graph->get_is_reverse(chunk.ref_end_handle);

        // offset_t is unsigned, so check for a gap before subtracting; then
        // trim_and_shift cuts the overlap away and gives the new start ID.
        if (!is_first_of_ref) {
            const SeamState& running = found->second;

            if (running.end_offset < chunk.ref_start_offset) {
                logger.error() << "REFERENCE offset gap of "
                               << (chunk.ref_start_offset - running.end_offset)
                               << " bp between \"" << running.source
                               << "\" (REFERENCE ends at offset " << running.end_offset
                               << ") and \"" << chunk.name << "\" (REFERENCE starts at offset "
                               << chunk.ref_start_offset
                               << "). --seam trim does not insert gaps." << endl;
            }
            offset_t overlap = running.end_offset - chunk.ref_start_offset;
            size_t left_boundary_len = dest->get_length(
                dest->get_handle(running.end_id, running.end_is_reverse));
            new_start_id = trim_and_shift(chunk, overlap, left_boundary_len);
        }

        // Ingesting renumbers IDs, so re-fetch the handles by the shift, then
        // seam onto the chunk to the left.
        int64_t shift = ingest(chunk);

        handle_t chunk_ref_end = dest->get_handle((nid_t)((int64_t)ref_end_id + shift),
                                                 end_is_reverse);
        handle_t new_end = chunk_ref_end;
        if (!is_first_of_ref) {
            handle_t new_start = dest->get_handle((nid_t)((int64_t)new_start_id + shift),
                                                 start_is_reverse);
            new_end = connect_seam(found->second, chunk.name, new_start,
                                   chunk_ref_end);
        }

        // Update the SeamState for this chunk.
        seams[chunk.ref_key] = SeamState{dest->get_id(new_end),
                                         dest->get_is_reverse(new_end),
                                         chunk.ref_end_offset,
                                         chunk.name};
    }

    merge_path_fragments();
    return std::move(dest);
}

/// Check whether a chunk holds multiple fragments of a logical path.
void GraphCombiner::check_no_split_paths(const Chunk& chunk) const {

    // Group by GroupKey, recording each path's start offset and name.
    map<GroupKey, vector<pair<offset_t, string>>> groups;
    chunk.graph->for_each_path_handle([&](const path_handle_t& path) {
        groups[group_key_of(*chunk.graph, path, phase_block_is_offset)]
            .emplace_back(path_start_offset(chunk.graph->get_subrange(path),
                                            chunk.graph->get_phase_block(path),
                                            phase_block_is_offset),
                          chunk.graph->get_path_name(path));
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
        logger.error() << "\"" << chunk.name << "\" holds " << pieces.size()
                       << " pieces of one path (" << names << "). vg combine joins pieces "
                       << "that arrive in separate inputs; a path already split inside one "
                       << "input is not something it will reassemble." << endl;
    }
}

/// Fill in `chunk`'s REFERENCE description and do some checks.
void GraphCombiner::describe_reference(Chunk& chunk) {

    // Collect every REFERENCE-sense path into refs.
    vector<path_handle_t> refs;
    chunk.graph->for_each_path_handle([&](const path_handle_t& p) {
        if (chunk.graph->get_sense(p) == PathSense::REFERENCE) {
            refs.push_back(p);
        }
    });

    // Two checks: there must be exactly one REFERENCE path.
    if (refs.empty()) {
        logger.error() << "No REFERENCE-sense path found in \"" << chunk.name
                       << "\"; --seam trim reads each chunk's position off its "
                       << "REFERENCE path, and needs exactly one per graph." << endl;
    }
    if (refs.size() > 1) {
        logger.error() << refs.size() << " REFERENCE-sense paths found in \"" << chunk.name
                       << "\"; --seam trim needs exactly one, so it knows which one "
                       << "gives the chunk its position." << endl;
    }
    path_handle_t ref = refs[0];

    // Record the two identity fields, ref_name and ref_key, and the REFERENCE
    // path's start and end.
    chunk.ref_name = chunk.graph->get_path_name(ref);
    chunk.ref_key = RefKey{chunk.graph->get_sample_name(ref),
                           chunk.graph->get_locus_name(ref),
                           chunk.graph->get_haplotype(ref),
                           chunk.graph->get_is_circular(ref)};
    chunk.ref_start_offset = path_start_offset(chunk.graph->get_subrange(ref),
                                               chunk.graph->get_phase_block(ref),
                                               phase_block_is_offset);
    chunk.ref_end_offset = chunk.ref_start_offset
        + (offset_t)compute_path_length(chunk.graph.get(), ref);

    chunk.ref_start_handle = chunk.graph->get_handle_of_step(chunk.graph->path_begin(ref));
    chunk.ref_end_handle = chunk.graph->get_handle_of_step(chunk.graph->path_back(ref));

    // Check boundary orientation: the REFERENCE path must visit its first and last
    // nodes forward.
    if (chunk.graph->get_is_reverse(chunk.ref_start_handle)
        || chunk.graph->get_is_reverse(chunk.ref_end_handle)) {
        logger.error() << "REFERENCE path \"" << chunk.ref_name << "\" in \"" << chunk.name
                       << "\" visits its boundary node in reverse orientation; --seam trim "
                       << "only supports forward-strand boundaries." << endl;
    }
}

/// Group each reference's chunks together and sort them by start position.
void GraphCombiner::sort_by_reference_offset(vector<Chunk>& chunks) {

    // Different references keep first-seen order; chunks of one reference sort by
    // ascending offset.
    map<RefKey, size_t> first_seen;
    for (const Chunk& chunk : chunks) {
        first_seen.emplace(chunk.ref_key, first_seen.size());
    }
    std::stable_sort(chunks.begin(), chunks.end(),
                     [&](const Chunk& a, const Chunk& b) {
                         size_t ra = first_seen.at(a.ref_key);
                         size_t rb = first_seen.at(b.ref_key);
                         if (ra != rb) {
                             return ra < rb;
                         }
                         return a.ref_start_offset < b.ref_start_offset;
                     });
}

/// Check that overlapping chunks spell the same sequence over the overlap.
void GraphCombiner::check_overlaps_agree(const vector<Chunk>& chunks) const {

    // The sort made each reference's chunks contiguous and offset-ascending, so only
    // neighbours can overlap.
    for (size_t i = 1; i < chunks.size(); ++i) {
        const Chunk& left = chunks[i - 1];
        const Chunk& right = chunks[i];

        // Skip a different reference, and a gap, which combine() reports itself.
        if (right.ref_key != left.ref_key
            || left.ref_end_offset <= right.ref_start_offset) {
            continue;
        }
        offset_t overlap = left.ref_end_offset - right.ref_start_offset;

        // Read the shared span off each side. A nested chunk reads short; the checks
        // in trim_and_shift catch that.
        string from_left = path_suffix(*left.graph,
                                       left.graph->get_path_handle(left.ref_name),
                                       overlap);
        string from_right = path_prefix(*right.graph,
                                        right.graph->get_path_handle(right.ref_name),
                                        overlap);

        // Find the first base they disagree on, if any.
        size_t shared = std::min(from_left.size(), from_right.size());
        size_t at = 0;
        while (at < shared && from_left[at] == from_right[at]) {
            ++at;
        }
        if (at < shared) {
            logger.error() << "Chunks \"" << left.name << "\" and \"" << right.name
                           << "\" overlap on REFERENCE coordinates ["
                           << right.ref_start_offset << ", " << left.ref_end_offset
                           << ") but disagree there: at offset "
                           << (right.ref_start_offset + at) << " \"" << left.name
                           << "\" has " << from_left[at] << " and \"" << right.name
                           << "\" has " << from_right[at]
                           << ". --seam trim drops the right chunk's copy of the overlap, "
                           << "so the two have to spell the same sequence." << endl;
        }
    }
}

/// Check that every path in the chunk starts at its REFERENCE start node.
void GraphCombiner::check_all_paths_start_at(const Chunk& chunk) const {
    chunk.graph->for_each_path_handle([&](const path_handle_t& p) {
        handle_t first = chunk.graph->get_handle_of_step(chunk.graph->path_begin(p));
        if (first != chunk.ref_start_handle) {
            logger.error() << "Path \"" << chunk.graph->get_path_name(p) << "\" in \""
                           << chunk.name << "\" does not start at the REFERENCE start node "
                           << "(id " << chunk.graph->get_id(chunk.ref_start_handle)
                           << "). --seam trim trims and joins at that one node, so every "
                           << "path in a chunk that has a chunk to its left has to begin "
                           << "there." << endl;
        }
    });
}

/// Check that every path in the chunk ends at its REFERENCE end node.
void GraphCombiner::check_all_paths_end_at(const Chunk& chunk) const {
    chunk.graph->for_each_path_handle([&](const path_handle_t& p) {
        handle_t last = chunk.graph->get_handle_of_step(chunk.graph->path_back(p));
        if (last != chunk.ref_end_handle) {
            logger.error() << "Path \"" << chunk.graph->get_path_name(p) << "\" in \""
                           << chunk.name << "\" does not end at the REFERENCE end node "
                           << "(id " << chunk.graph->get_id(chunk.ref_end_handle)
                           << "). --seam trim trims and joins at that one node, so every "
                           << "path in a chunk that has a chunk to its right has to end "
                           << "there." << endl;
        }
    });
}

/// Trim `chunk`'s overlap with the chunk to its left, and fix up path coordinates.
nid_t GraphCombiner::trim_and_shift(Chunk& chunk, offset_t overlap,
                                    size_t left_boundary_len) {

    // A raw pointer to this chunk's graph, and the ID of its REFERENCE start node.
    MutablePathDeletableHandleGraph* graph = chunk.graph.get();
    nid_t new_start_id = graph->get_id(chunk.ref_start_handle);

    // Two checks: the overlap must be strictly shorter than the first node.
    if (overlap > 0) {
        size_t start_len = graph->get_length(chunk.ref_start_handle);
        if (overlap > start_len) {
            logger.error() << "REFERENCE overlap of " << overlap << " bp before \""
                           << chunk.name << "\" exceeds the length of its first node ("
                           << start_len << " bp on node " << new_start_id
                           << "), so the overlap is not all on that node. Re-chunk so it "
                           << "is." << endl;
        }
        if (overlap == start_len) {
            logger.error() << "REFERENCE overlap of " << overlap
                           << " bp exactly equals the length of the first node of \""
                           << chunk.name << "\"; trimming it would leave a zero-length "
                           << "node. Re-chunk so that first node extends past the overlap."
                           << endl;
        }

        // Split the first node in two: the overlapping piece, and the remainder,
        // which becomes the new start.
        vector<handle_t> pieces = graph->divide_handle(chunk.ref_start_handle,
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
            logger.error() << "The first " << overlap << " bp of \"" << chunk.name
                           << "\" overlap the chunk to its left and have to be trimmed, "
                           << "but " << inbound << " edge(s) arrive at node "
                           << graph->get_id(overlapping)
                           << " there and would be lost. --seam trim needs the chunk "
                           << "boundary to be a clean left edge." << endl;
        }
        graph->destroy_handle(overlapping);
    }

    // Compute the metadata shift.
    int64_t meta_shift = (int64_t)overlap - (fuse ? (int64_t)left_boundary_len : 0);
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
            s.subrange.first =
                (offset_t)max<int64_t>(0, (int64_t)s.subrange.first + meta_shift);
        } else if (phase_block_is_offset && s.phase_block != PathMetadata::NO_PHASE_BLOCK) {
            s.phase_block = (size_t)max<int64_t>(0, (int64_t)s.phase_block + meta_shift);
        } else {
            s.subrange = subrange_t{(offset_t)max<int64_t>(0, meta_shift),
                                    PathMetadata::NO_END_POSITION};
        }
        path_handle_t rebuilt = graph->create_path(s.sense, s.sample, s.locus, s.haplotype,
                                                  s.phase_block, s.subrange, s.is_circular);
        for (handle_t h : s.steps) {
            graph->append_step(rebuilt, h);
        }
    }
    return new_start_id;
}

/// Move `chunk` into the accumulator; return how far its IDs were shifted.
int64_t GraphCombiner::ingest(Chunk& chunk) {

    if (!dest) {
        // The first graph becomes the accumulator: nothing to copy or renumber.
        dest = std::move(chunk.graph);
        max_node_id = dest->max_node_id();
        return 0;
    }

    int64_t shift = 0;
    if (seam == SeamPolicy::SHARED_IDS) {
        share_nodes_by_id(chunk);
        max_node_id = std::max(max_node_id, (int64_t)chunk.graph->max_node_id());
    } else {
        shift = renumber_out_of_way(chunk.graph.get(), max_node_id);
        max_node_id = std::max(max_node_id, (int64_t)chunk.graph->max_node_id());
        handlealgs::copy_handle_graph(chunk.graph.get(), dest.get());
    }
    copy_paths_checked(chunk);
    // Free the source rather than holding every input until combine finishes.
    chunk.graph.reset();
    return shift;
}

/// Copy nodes and edges, keeping IDs as given and distinguishing nodes by ID.
void GraphCombiner::share_nodes_by_id(Chunk& chunk) {

    chunk.graph->for_each_handle([&](const handle_t& source) {
        nid_t id = chunk.graph->get_id(source);
        string sequence = chunk.graph->get_sequence(source);
        if (!dest->has_node(id)) {
            dest->create_handle(sequence, id);
            return;
        }
        // Under this policy a repeated ID is the same node, so two sequences for
        // one ID is a contradiction with no way to keep both.
        if (dest->get_sequence(dest->get_handle(id)) != sequence) {
            logger.error() << "Node " << id << " has different sequences in \"" << chunk.name
                           << "\" and in an earlier input. Under --seam shared a repeated "
                           << "node ID means the same node, so the sequences have to match."
                           << endl;
        }
    });
    chunk.graph->for_each_edge([&](const edge_t& source) {
        handle_t from = dest->get_handle(chunk.graph->get_id(source.first),
                                        chunk.graph->get_is_reverse(source.first));
        handle_t to = dest->get_handle(chunk.graph->get_id(source.second),
                                       chunk.graph->get_is_reverse(source.second));
        if (!dest->has_edge(from, to)) {
            dest->create_edge(from, to);
        }
    });
}

/// Copy paths, rejecting duplicate names.
void GraphCombiner::copy_paths_checked(Chunk& chunk) {

    chunk.graph->for_each_path_handle([&](const path_handle_t& source) {
        string path_name = chunk.graph->get_path_name(source);
        if (dest->has_path(path_name)) {
            logger.error() << "Path \"" << path_name << "\" is in \"" << chunk.name
                           << "\" and in an earlier input. Pieces of one path have to be "
                           << "told apart by name, so give them distinct subranges." << endl;
        }
        handlealgs::copy_path(chunk.graph.get(), source, dest.get());
    });
}

/// Join the left chunk's end to this chunk's start; return the reference's new end.
handle_t GraphCombiner::connect_seam(const SeamState& running, const string& chunk_name,
                                     handle_t new_start, handle_t chunk_ref_end) {

    handle_t left = dest->get_handle(running.end_id, running.end_is_reverse);

    // If not fusing nodes, just add an edge.
    if (!fuse) {
        if (!dest->has_edge(left, new_start)) {
            dest->create_edge(left, new_start);
        }
        return chunk_ref_end;
    }

    // can_fuse_nodes decides whether the two nodes can be fused.
    if (!algorithms::can_fuse_nodes(*dest, left, new_start)) {
        // Welding buries both facing ends, and something is attached to one.
        logger.error() << "--fuse cannot weld node " << dest->get_id(left) << " from \""
                       << running.source << "\" to node " << dest->get_id(new_start)
                       << " from \"" << chunk_name << "\": one of them still has an edge "
                       << "on the side being welded, which the weld would have to drop. "
                       << "--seam trim without --fuse joins them with an edge instead and "
                       << "keeps everything." << endl;
    }

    // Note the end node before fusing, fuse into a new node, and update the max ID.
    bool ref_ends_at_seam = (dest->get_id(chunk_ref_end) == dest->get_id(new_start));
    handle_t fused = algorithms::fuse_nodes(dest.get(), left, new_start);
    max_node_id = std::max(max_node_id, (int64_t)dest->get_id(fused));
    return ref_ends_at_seam ? fused : chunk_ref_end;
}

/// Merge fragments from the same path into a single path spanning their entire
/// coordinate range.
void GraphCombiner::merge_path_fragments() {

    struct Fragment {
        path_handle_t handle;            ///< Handle to the path.
        subrange_t subrange;             ///< Coordinate range.
        size_t phase_block;              ///< Phase block.
        size_t length;                   ///< Sequence length.
        bool positioned_by_phase_block;  ///< Whether the phase block gives its position.
    };
    map<GroupKey, vector<Fragment>> groups;

    // Walk every path in dest, filling in a Fragment for each.
    dest->for_each_path_handle([&](const path_handle_t& path) {
        Fragment f;
        f.handle = path;
        f.subrange = dest->get_subrange(path);
        f.phase_block = dest->get_phase_block(path);
        f.length = compute_path_length(dest.get(), path);
        f.positioned_by_phase_block = phase_block_is_offset
            && f.subrange == PathMetadata::NO_SUBRANGE
            && f.phase_block != PathMetadata::NO_PHASE_BLOCK;
        groups[group_key_of(*dest, path, phase_block_is_offset)].push_back(std::move(f));
    });

    // Lambdas computing a fragment's start and end coordinates.
    auto fragment_start = [&](const Fragment& f) -> offset_t {
        return path_start_offset(f.subrange, f.phase_block, phase_block_is_offset);
    };
    auto fragment_end = [&](const Fragment& f) -> offset_t {
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
        bool is_circular = std::get<4>(group.first);

        std::sort(frags.begin(), frags.end(),
                  [&](const Fragment& a, const Fragment& b) {
                      return fragment_start(a) < fragment_start(b);
                  });

        // Detect any gap, updating the furthest reference offset reached so far.
        offset_t covered_end = fragment_end(frags.front());
        for (size_t i = 1; i < frags.size(); ++i) {
            offset_t next_start = fragment_start(frags[i]);
            if (covered_end < next_start) {
                logger.error() << "Gap of " << (next_start - covered_end)
                               << " bp before fragment \""
                               << dest->get_path_name(frags[i].handle)
                               << "\"; the fragments before it only reach offset "
                               << covered_end << ". vg combine does not insert gaps."
                               << endl;
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
                logger.error() << "--phase-block-is-offset given, but phase blocks of \""
                               << dest->get_path_name(frags[i].handle)
                               << "\" do not line up as coordinates (block "
                               << frags[i].phase_block << " + length " << frags[i].length
                               << " != block " << frags[i + 1].phase_block
                               << "). Are these phase block numbers rather than offsets?"
                               << endl;
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
        // as coordinates folded it into the subrange, so 0 stands in there.
        size_t group_phase_block = std::get<5>(group.first);
        size_t merged_phase_block = PathMetadata::NO_PHASE_BLOCK;
        if (sense == PathSense::HAPLOTYPE) {
            merged_phase_block = (group_phase_block == PathMetadata::NO_PHASE_BLOCK)
                ? 0 : group_phase_block;
        }
        string merged_name = PathMetadata::create_path_name(sense, sample, locus, haplotype,
                                                           merged_phase_block,
                                                           combined_subrange);

        // Can't lean on create_path(): PackedGraph throws on a duplicate, HashGraph
        // quietly orphans the old path, and VG just appends to it. So check here.
        if (dest->has_path(merged_name)) {
            bool is_own_fragment = false;
            for (const Fragment& f : frags) {
                if (dest->get_path_name(f.handle) == merged_name) {
                    is_own_fragment = true;
                    break;
                }
            }
            if (!is_own_fragment) {
                logger.error() << "The merge of " << frags.size() << " fragments would be "
                               << "named \"" << merged_name << "\", which is already an "
                               << "unrelated path in the output." << endl;
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
            // `skip` is the shared run at the seam: `run` is the candidate length, tried
            // longest first down to 1 (a run of 0 would match trivially), and `k` compares
            // `previous`'s tail against `steps`'s head.
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
                if (skip == 0 && !dest->has_edge(previous->back(), steps.front())) {
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
