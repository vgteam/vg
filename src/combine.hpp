#ifndef VG_COMBINE_HPP_INCLUDED
#define VG_COMBINE_HPP_INCLUDED

/** \file combine.hpp
 *
 * Combining graphs whose paths are fragments of the same paths, told apart by
 * subrange, and merging the fragments back into whole paths.
 *
 * Fragment: one path in an input, which is a piece of a contig.
 * Base path: the whole contig a reference path is a piece of.
 * Chunk: an input graph, plus where its one reference path lies on the base path.
 * Seam: where two neighbouring chunks are joined, after their overlap is trimmed off.
 * Seam nodes: at a seam, the last node of the left chunk and the first node of the right chunk.
 */

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "handle.hpp"
#include "log.hpp"

namespace vg {

/// How overlap between the inputs is resolved.
enum class SeamPolicy {
    RENUMBER,    ///< Inputs don't overlap; node IDs are shifted to avoid collisions.
    SHARED_IDS,  ///< Nodes with the same ID are the same node.
    COORD_TRIM   ///< Overlap in reference coordinates is trimmed off.
};

/// Thrown when GraphCombiner's inputs cannot be combined.
class CombineError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// An input's reference path is missing, ambiguous, or cannot be trimmed at.
class ReferencePathError : public CombineError {
public:
    using CombineError::CombineError;
};

/// Neighbouring chunks of a base path do not fit together where they meet.
class ChunkOverlapError : public CombineError {
public:
    using CombineError::CombineError;
};

/// Inputs contradict each other, or one input contradicts itself.
class InputConflictError : public CombineError {
public:
    using CombineError::CombineError;
};

/// Two seam nodes cannot be fused into one.
class FuseError : public CombineError {
public:
    using CombineError::CombineError;
};

/// Fragments of one path cannot be merged into a single path.
class FragmentMergeError : public CombineError {
public:
    using CombineError::CombineError;
};

/// A path's subrange claims a different end than its sequence gives it.
class PathLengthMismatchError : public CombineError {
public:
    using CombineError::CombineError;
};

/// Combine fragmented graphs into a whole graph.
class GraphCombiner {
public:
    /// A graph to combine.
    struct Input {
        std::string name;  ///< Identifies the graph in error messages, usually the original filename.
        std::unique_ptr<MutablePathDeletableHandleGraph> graph;  ///< The graph itself.
    };

    /// Make a combiner. seam_policy is the mode for resolving overlap between the inputs;
    /// fuse merges the two seam nodes into one where trimmed inputs are joined,
    /// instead of adding an edge; phase_block_is_offset reads a path's phase block
    /// as its start coordinate; log_context prefixes progress messages.
    GraphCombiner(SeamPolicy seam_policy, bool fuse, bool phase_block_is_offset,
                  const std::string& log_context);

    /// Combine the input graphs into one, which is empty if there are none, or throw a
    /// CombineError if they cannot be combined.
    std::unique_ptr<MutablePathDeletableHandleGraph> combine(std::vector<Input>&& inputs);

protected:
    /// Identity of a base path: sample, locus, haplotype.
    using BasePathKey = std::tuple<std::string, std::string, size_t>;

    /// A path's identity as a comparable tuple: sense, sample, locus, haplotype, phase block.
    using GroupKey = std::tuple<handlegraph::PathSense, std::string, std::string,
                                size_t, size_t>;

    /// A path meeting a seam: where it starts or ends on its own coordinates, and its name for error messages.
    struct PathAtSeam {
        handlegraph::offset_t offset;  ///< Its start or end coordinate.
        std::string name;              ///< Its name.
    };

    /// The current end of a base path in the combined graph, where the next input is seamed on.
    struct SeamState {
        handlegraph::nid_t end_id;         ///< Node the base path so far ends on.
        bool end_is_reverse;               ///< Its orientation.
        handlegraph::offset_t end_offset;  ///< Base path offset it ends at.
        std::string last_input_name;       ///< Name of the input that last extended this base path, for error messages.
        std::map<GroupKey, PathAtSeam> path_ends;  ///< The paths ending on that node, by GroupKey: where each ends, and its name.
    };

    /// Metadata about an input's reference path, as it is in that input's graph.
    struct ReferenceDescription {
        std::string path_name;                   ///< Name of its reference path in this chunk's graph, such as GRCh38#0#chr1[100-200], for looking it up and for error reporting.
        BasePathKey base_path_key;               ///< The base path its reference path is a piece of, such as {"GRCh38", "chr1", 0}.
        handle_t start_handle;                   ///< First node of that path.
        handle_t end_handle;                     ///< Last node of that path.
        handlegraph::offset_t start_offset = 0;  ///< Where that path starts on the base path.
        handlegraph::offset_t end_offset = 0;    ///< Where that path ends on the base path.
        bool is_circular = false;                ///< Whether that path is circular.
    };

    /// An input, plus a description of its reference path.
    struct Chunk : Input {
        ReferenceDescription ref;
    };

    /// Build a path's GroupKey.
    GroupKey group_key_of(const PathHandleGraph& graph, const path_handle_t& path) const;

    /// Return a path's start offset.
    handlegraph::offset_t path_start_offset(const handlegraph::subrange_t& subrange,
                                            size_t phase_block) const;
    /// Increase graph nodes IDs so that the minimum node ID is greater than "above".
    /// Returns the amount IDs were incremented by, which is 0 if no IDs needed to change.
    static int64_t shift_ids_above(MutableHandleGraph* graph, int64_t above);
    /// Read up to `length` bp off the front of a path.
    static std::string path_prefix(const PathHandleGraph& graph, const path_handle_t& path,
                                   size_t length);
    /// Read up to `length` bp off the back of a path.
    static std::string path_suffix(const PathHandleGraph& graph, const path_handle_t& path,
                                   size_t length);

    /// Throw PathLengthMismatchError if a path in the input has a subrange end other than its start plus its sequence length.
    static void validate_input(const Input& input);
    /// Throw InputConflictError or ChunkOverlapError if `chunks` cannot be trimmed and joined end to end.
    /// The chunks must already be sorted by sort_by_base_path_offset().
    void validate_chunks(const std::vector<Chunk>& chunks) const;
    /// Describe `input`'s reference path. Throw ReferencePathError if it does not have
    /// exactly one, or if that path is empty or visits its first or last node in reverse.
    ReferenceDescription describe_reference(const Input& input) const;
    /// Group each base path's chunks together and sort them by start position.
    void sort_by_base_path_offset(std::vector<Chunk>& chunks);

    /// Cut the overlap off the front of `chunk`, and shift the coordinates of the paths that start there to match.
    /// `chunk` is the right chunk, being joined onto the end of the base path built so far.
    /// `overlap` is how many bp of `chunk` the base path already covers; it must be shorter than `chunk`'s first node.
    /// `left_end_node_len` is the length of the base path's current last node, which fusing adds in front of `chunk`.
    /// Returns the ID of `chunk`'s new first node.
    /// Throws ChunkOverlapError if the shift would move a path's start coordinate below 0.
    handlegraph::nid_t trim_and_shift(Chunk& chunk, handlegraph::offset_t overlap,
                                      size_t left_end_node_len);

    /// Return each path on `node` with its start (at_start) or end coordinate, by GroupKey; throw ChunkOverlapError if one path has two pieces there.
    std::map<GroupKey, PathAtSeam> paths_at_seam(const Chunk& chunk, handle_t node,
                                                 bool at_start) const;
    /// Throw ChunkOverlapError if a path's pieces on the two sides of a seam do not end and start at the same coordinate, or, when fusing, if a path is on only one side.
    void check_paths_meet_at_seam(const SeamState& left,
                                  const std::map<GroupKey, PathAtSeam>& right,
                                  size_t left_end_node_len,
                                  const std::string& input_name) const;

    /// Add `input`'s nodes, edges and paths to dest, and free its graph.
    /// Unless nodes are shared by ID, its node IDs are first shifted above every ID already in dest; returns the shift, or 0 if none.
    int64_t ingest(Input& input);
    /// Copy `input`'s nodes and edges into dest, keeping their IDs; a node whose ID is already in dest is the same node and is not copied again.
    /// Throw InputConflictError if such a node has a different sequence.
    void share_nodes_by_id(Input& input);
    /// Copy `input`'s paths into dest; throw InputConflictError if a path name is already in dest.
    void copy_paths_checked(Input& input);

    /// Join a chunk already copied into dest onto the end of the base path built so far, by an edge or, when fusing, by fusing the two nodes.
    /// `base_path_end` is where that base path ends.
    /// `chunk_start` and `chunk_end` are the chunk's reference path's first and last nodes in dest.
    /// `input_name` names the chunk's input in error messages.
    /// Returns the base path's new last node.
    /// Throws FuseError if fusing is on and the two nodes cannot be fused.
    handle_t connect_seam(const SeamState& base_path_end, handle_t chunk_start,
                          handle_t chunk_end, const std::string& input_name);

    /// Merge each path's consecutive fragments with no gap between them into one path.
    /// Throw FragmentMergeError on a gap in a reference or circular path.
    void merge_path_fragments();

    SeamPolicy seam_policy;      ///< The mode for resolving overlap between the inputs.
    bool fuse;                   ///< Merge the two seam nodes into one where trimmed inputs are joined, instead of adding an edge.
    bool phase_block_is_offset;  ///< Read a path's phase block as its start coordinate.
    Logger logger;               ///< Reports progress.

    std::unique_ptr<MutablePathDeletableHandleGraph> dest;  ///< The combined graph built so far.
    int64_t max_node_id = 0;     ///< Largest node ID used in dest so far.
    std::map<BasePathKey, SeamState> seams;  ///< For each base path, the node, orientation, offset and input it ends at so far.
};

}

#endif
