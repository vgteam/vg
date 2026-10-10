#ifndef VG_COMBINE_HPP_INCLUDED
#define VG_COMBINE_HPP_INCLUDED

/** \file combine.hpp
 *
 * Combining graphs that hold fragments of the same logical paths.
 */

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "handle.hpp"
#include "log.hpp"

namespace vg {

namespace unittest {
/// Reaches into GraphCombiner's internals; see src/unittest/combine.cpp.
class TestGraphCombiner;
}

/// The three modes of vg combine, selected by --seam.
enum class SeamPolicy {
    RENUMBER,    ///< --seam renumber
    SHARED_IDS,  ///< --seam shared
    COORD_TRIM   ///< --seam trim
};

/// Combine fragmented graphs into a whole graph.
class GraphCombiner {
public:
    struct Input {
        std::string name;
        std::unique_ptr<MutablePathDeletableHandleGraph> graph;
    };

    /// Constructor.
    GraphCombiner(SeamPolicy seam, bool fuse, bool phase_block_is_offset);

    /// Combine the input graphs.
    std::unique_ptr<MutablePathDeletableHandleGraph> combine(std::vector<Input>&& inputs);

private:
    /// Identity of a reference: sample, locus, haplotype, circular.
    using RefKey = std::tuple<std::string, std::string, size_t, bool>;

    /// Seaming state.
    struct SeamState {
        handlegraph::nid_t end_id;         ///< Node the reference so far ends on.
        bool end_is_reverse;               ///< Its orientation.
        handlegraph::offset_t end_offset;  ///< REFERENCE offset it ends at.
        std::string source;                ///< Name of chunk/inputs.
    };

    /// An input, plus the REFERENCE description COORD_TRIM positions it by.
    struct Chunk {
        std::string name;                                        ///< Name of the input.
        std::unique_ptr<MutablePathDeletableHandleGraph> graph;  ///< The input graph.

        std::string ref_name;                        ///< Its REFERENCE path.
        handle_t ref_start_handle;                   ///< First node of that path.
        handle_t ref_end_handle;                     ///< Last node of that path.
        handlegraph::offset_t ref_start_offset = 0;  ///< Where that path starts.
        handlegraph::offset_t ref_end_offset = 0;    ///< Where that path ends.
        RefKey ref_key;                              ///< Which reference it belongs to.
    };

    /// A path's identity as a comparable tuple: sense, sample, locus, haplotype, circular, phase block.
    using GroupKey = std::tuple<handlegraph::PathSense, std::string, std::string,
                                size_t, bool, size_t>;

    /// Compute the length of a path.
    static size_t compute_path_length(const PathHandleGraph* graph, const path_handle_t& path);
    /// Return a path's start offset.
    static handlegraph::offset_t path_start_offset(const handlegraph::subrange_t& subrange,
                                                   size_t phase_block, bool phase_block_is_offset);
    /// Raise every node ID above the previous graph's max ID; return the shift.
    static int64_t renumber_out_of_way(MutableHandleGraph* graph, int64_t above);
    /// Build a path's GroupKey.
    static GroupKey group_key_of(const PathHandleGraph& graph, const path_handle_t& path,
                                 bool phase_block_is_offset);
    /// Read up to `length` bp off the front of a path.
    static std::string path_prefix(const PathHandleGraph& graph, const path_handle_t& path,
                                   size_t length);
    /// Read up to `length` bp off the back of a path.
    static std::string path_suffix(const PathHandleGraph& graph, const path_handle_t& path,
                                   size_t length);

    /// Check whether a chunk holds multiple fragments of a logical path.
    void check_no_split_paths(const Chunk& chunk) const;
    /// Fill in `chunk`'s REFERENCE description and do some checks.
    void describe_reference(Chunk& chunk);
    /// Group each reference's chunks together and sort them by start position.
    void sort_by_reference_offset(std::vector<Chunk>& chunks);
    /// Check that overlapping chunks spell the same sequence over the overlap.
    void check_overlaps_agree(const std::vector<Chunk>& chunks) const;

    /// Check that every path in the chunk starts at its REFERENCE start node.
    void check_all_paths_start_at(const Chunk& chunk) const;
    /// Check that every path in the chunk ends at its REFERENCE end node.
    void check_all_paths_end_at(const Chunk& chunk) const;

    /// Trim `chunk`'s overlap with the chunk to its left, and fix up path coordinates.
    handlegraph::nid_t trim_and_shift(Chunk& chunk, handlegraph::offset_t overlap,
                                      size_t left_boundary_len);

    /// Move `chunk` into the accumulator; return how far its IDs were shifted.
    int64_t ingest(Chunk& chunk);
    /// Copy nodes and edges, keeping IDs as given and distinguishing nodes by ID.
    void share_nodes_by_id(Chunk& chunk);
    /// Copy paths, rejecting duplicate names.
    void copy_paths_checked(Chunk& chunk);

    /// Join the left chunk's end to this chunk's start; return the reference's new end.
    handle_t connect_seam(const SeamState& seam, const std::string& chunk_name,
                          handle_t new_start, handle_t chunk_ref_end);

    /// Merge fragments from the same path into a single path spanning their entire
    /// coordinate range.
    void merge_path_fragments();

    SeamPolicy seam;
    bool fuse;
    bool phase_block_is_offset;
    Logger logger;

    std::unique_ptr<MutablePathDeletableHandleGraph> dest;
    int64_t max_node_id = 0;
    std::map<RefKey, SeamState> seams;

    friend class unittest::TestGraphCombiner;
};

}

#endif
