#ifndef VG_FEATURE_LIFT_HPP_INCLUDED
#define VG_FEATURE_LIFT_HPP_INCLUDED

/**
 * \file feature_lift.hpp
 *
 * Lifting features from linear sequences to the graph, the gref and the paths.
 *
 * A feature is a BED record on a linear sequence. The sequence is a cover
 * sequence, which is a gref path or a tig of vg unitigs (see unitigs.hpp), or
 * any path in the GBZ. Each of these sequences spells a walk in the graph, so
 * a feature is a WalkInterval. The interval has a gref placement, and it has
 * places on the paths that visit it.
 *
 * A place is exact when a path visits the whole walk. Exact places suit short
 * features found by scanning a cover. Every window of at most k bases of a
 * haplotype lies within some cover sequence, so a tool that looks at no more
 * than k bases at a time finds each feature on every haplotype. Cover
 * sequences overlap by fewer than k bases, so a window shorter than k can lie
 * within several of them. Its copies lift to the same interval and are kept
 * once.
 *
 * A place is a liftover when a path visits enough of the walk's bases, in
 * order. Liftovers suit long features, such as genes, which haplotypes rarely
 * share whole.
 *
 * Each feature also gets a path coordinate (see PathCoordinate), which names it
 * exactly as an interval on one path that carries it.
 */

#include "gref.hpp"
#include "walk_interval.hpp"

#include <functional>
#include <iosfwd>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace vg {

/// Index of the strand in the fields of a feature. The fields start at the 4th
/// BED column, so this is the 6th.
constexpr size_t BED_STRAND_FIELD = 2;

/**
 * A feature in linear coordinates.
 *
 * It covers bases [start, end) of the named sequence. It reads them in reverse
 * if its strand is "-".
 */
struct BedFeature {
    std::string sequence;
    size_t start = 0;
    size_t end = 0;
    /// BED columns from the 4th on: name, score, strand and any more.
    std::vector<std::string> fields;
};

/**
 * Read the features in a BED file, skipping comments and track and browser lines.
 *
 * @throws std::runtime_error on a malformed line.
 */
std::vector<BedFeature> read_bed_features(std::istream& in);

/**
 * A feature's path coordinate: its interval on the first path, in priority
 * order, that carries it whole.
 *
 * Gref paths come first in the order, so a feature that lies along one is
 * named there. A feature off the gref is named on the first reference or
 * assembly path that carries it (see PathPriority). Paths keep their
 * coordinates when a graph is clipped, filtered or sampled, so these names
 * hold across such graphs while the naming path is in them. No two graph
 * intervals share a path coordinate.
 */
struct PathCoordinate {
    /// The path, named as for_each_path_start() names it.
    std::string path;
    /// Place of the path in the priority order.
    size_t rank = 0;
    /// Coordinates on the full path.
    uint64_t start = 0;
    uint64_t end = 0;
    /// True if the feature reads the path in reverse.
    bool is_reverse = false;

    /// Return true if this coordinate comes first in the priority order.
    bool operator<(const PathCoordinate& other) const;
};

/**
 * A feature in the graph.
 *
 * It is a walk interval and the BED fields that describe it. A stranded
 * feature reads along its interval, and its strand field is "+". A feature
 * without a strand has "." there, or no strand field.
 */
struct GraphFeature {
    WalkInterval interval;
    std::vector<std::string> fields;
    /// Where the interval lies on the gref. Empty text without a gref.
    GrefPlacement placement;
    /// Number of visits to the interval, in either orientation, by paths other
    /// than gref paths.
    size_t support = 0;
    /// The path coordinate. make_graph_feature() gives one to a feature along
    /// a gref path; others get theirs from the paths that carry them.
    std::optional<PathCoordinate> coordinate;

    /// Return true if the feature has a strand.
    bool is_stranded() const;
};

/**
 * Make the graph feature for an interval of a linear sequence and the BED
 * fields of a feature on it.
 *
 * A stranded feature is oriented along its strand, which reverses the interval
 * if the strand is "-". A feature without a strand takes the orientation with
 * the most bases on forward gref strands. On a tie, it takes the orientation
 * whose walk sorts first. The gref may be null. A feature along one gref path
 * gets its path coordinate there.
 */
GraphFeature make_graph_feature(const gbwtgraph::GBWTGraph& graph, const GrefIndex* gref,
                                WalkInterval interval, std::vector<std::string> fields);

/**
 * The priority order of the paths that name features off the gref.
 *
 * Paths whose names start with the given prefixes come first, in the order of
 * the prefixes. Then come reference paths, and then the rest. Each group is in
 * name order. All ranks come after those of the gref paths.
 */
class PathPriority {
public:
    /// Rank the paths of the GBZ other than gref paths, starting at first_rank.
    PathPriority(const gbwtgraph::GBZ& gbz, const std::vector<std::string>& prefixes,
                 const std::string& gref_prefix, size_t first_rank);

    /// Return the rank of the path with this name, as for_each_path_start()
    /// names it. Requires a ranked path.
    size_t rank(const std::string& name) const;

private:
    std::unordered_map<std::string, size_t> ranks;
};

/**
 * A directory of the sequences of a cover.
 *
 * It knows the walk that each gref path and each loaded tig spells.
 */
class Cover {
public:
    /// Make a cover of the gref paths of the index, which may be null, with no tigs yet.
    Cover(const gbwtgraph::GBWTGraph& graph, const GrefIndex* gref);

    /**
     * Load the tigs of a GAF from vg unitigs for which wanted(name) is true.
     *
     * @throws std::runtime_error on a malformed line or a repeated tig name.
     */
    void add_tigs(std::istream& gaf, const std::function<bool(const std::string&)>& wanted);

    /// Return true if the name is that of a gref path.
    bool is_gref_path(const std::string& name) const;
    /// Return the length of the named sequence, or nothing if the cover lacks it.
    std::optional<size_t> sequence_length(const std::string& name) const;
    /// Return bases [start, end) of the named sequence. Requires the sequence
    /// to be in the cover, and start < end <= its length.
    WalkInterval interval(const std::string& name, size_t start, size_t end) const;

private:
    const gbwtgraph::GBWTGraph& graph;
    /// Spells the gref paths. Null without a gref.
    const GrefIndex* gref;
    /// The walk interval of each loaded tig, by name.
    std::unordered_map<std::string, WalkInterval> tigs;
};

/**
 * A directory of paths in a GBZ, by the names that for_each_path_start() gives.
 *
 * A path can be split into fragments, such as the pieces of an assembly contig
 * in a clipped graph. The directory finds the walk of a stretch of the full
 * path in the fragment that holds it.
 */
class PathDirectory {
public:
    /// Load the fragments of each path whose name passes the filter.
    PathDirectory(const gbwtgraph::GBZ& gbz, const std::function<bool(const std::string&)>& wanted);

    /// Return true if the directory has fragments of the named path.
    bool has_path(const std::string& name) const;
    /// Return bases [start, end) of the named full path, or nothing if no
    /// fragment holds them all. Requires start < end.
    std::optional<WalkInterval> interval(const std::string& name, uint64_t start, uint64_t end) const;

private:
    struct Fragment {
        /// Where the fragment starts on the full path.
        uint64_t start;
        PathSteps steps;
    };
    /// The fragments of each path, by name, sorted by start.
    std::unordered_map<std::string, std::vector<Fragment>> fragments;
};

/**
 * Sort the features along the gref, then by interval and fields, and merge the
 * copies of each feature.
 *
 * Copies have the same interval and the same fields.
 */
void sort_and_merge_features(std::vector<GraphFeature>& features);

/// A place where a path visits a feature.
struct PathOccurrence {
    /// Coordinates on the full path of the first base of the place and of the
    /// base after its last.
    uint64_t start;
    uint64_t end;
    /// Index of the feature.
    size_t feature;
    /// True if the path reads the feature in reverse.
    bool is_reverse;
    /// Number of feature bases that the path visits at this place.
    size_t matched;
    /// The stretches of the place that hold feature bases, as [start, end)
    /// coordinates on the full path, in increasing order.
    std::vector<std::pair<uint64_t, uint64_t>> blocks;
};

/**
 * A matcher that finds where paths visit the whole walks of graph features.
 *
 * It indexes the features by the first node of their walks, read either way.
 * At each step of a path, it checks the features indexed at that node against
 * the steps that follow.
 */
class FeatureMatcher {
public:
    /// Index the features, which must outlive the matcher.
    FeatureMatcher(const gbwtgraph::GBWTGraph& graph, const std::vector<GraphFeature>& features);

    /**
     * Call visit(feature, step, is_reverse) for each place where the path
     * visits the walk of a feature.
     *
     * The feature is its index in the features. The step is the path's first
     * step on the walk. The flag is_reverse is true if the path reads the walk
     * in reverse.
     */
    void for_each_match(const gbwt::vector_type& path,
                        const std::function<void(size_t, size_t, bool)>& visit) const;

private:
    const gbwtgraph::GBWTGraph& graph;
    const std::vector<GraphFeature>& features;
    handlegraph::nid_t min_id = 0;
    /// Where the entries of each node start, by node ID minus min_id, followed
    /// by the number of entries.
    std::vector<size_t> starts;
    /// The walks that start on each node, grouped by node. Each is a feature
    /// index shifted left by one, with is_reverse in the low bit. A walk read
    /// in reverse starts on its last node.
    std::vector<size_t> entries;
};

/**
 * An index of the nodes on the walks of graph features.
 *
 * It turns each step of a path onto such a node into an anchor, which is the
 * stretch of the feature on that node at its place on the path.
 */
class FeatureNodeIndex {
public:
    /// A stretch of a feature at a place on a path.
    struct Anchor {
        /// Index of the feature.
        size_t feature;
        /// Bases [from, to) of the feature.
        size_t from;
        size_t to;
        /// Path coordinate of the first base of the stretch on the path.
        uint64_t start;
        /// True if the path reads the stretch in reverse.
        bool is_reverse;
    };

    /// Index the features, which must outlive the index.
    FeatureNodeIndex(const gbwtgraph::GBWTGraph& graph, const std::vector<GraphFeature>& features);

    /// Call visit(anchor) for each step of the path onto a node of a feature
    /// walk, in path order. The path coordinate of each step is in coords.
    void for_each_anchor(const gbwt::vector_type& path, const std::vector<uint64_t>& coords,
                         const std::function<void(const Anchor&)>& visit) const;

private:
    /// A stretch of a feature on one node of its walk.
    struct Entry {
        size_t feature;
        /// Bases [from, to) of the feature.
        size_t from;
        size_t to;
        /// Offsets of those bases on the node, as the walk reads it.
        size_t node_from;
        size_t node_to;
        /// True if the walk reads the node in reverse.
        bool is_reverse;
    };

    const gbwtgraph::GBWTGraph& graph;
    handlegraph::nid_t min_id = 0;
    /// Where the entries of each node start, by node ID minus min_id, followed
    /// by the number of entries.
    std::vector<size_t> starts;
    /// The stretches of the features on each node, grouped by node.
    std::vector<Entry> entries;
};

/**
 * Chain the anchors that one path gives a set of features into liftovers.
 *
 * A chain is a run of anchors of one feature that the path reads in the same
 * orientation, in the order of the feature's bases. Each anchor ends before the
 * next one starts on the path, with a gap of at most the feature's length.
 *
 * Chains are taken best first. The first chain holds the most feature bases,
 * the next one holds the most among the anchors left, and so on. Return each
 * chain that holds at least min_match of its feature's bases, as a place on
 * the path.
 */
std::vector<PathOccurrence> chain_anchors(std::vector<FeatureNodeIndex::Anchor> anchors,
                                          const std::vector<GraphFeature>& features, double min_match);

/**
 * Find where each path visits the whole walks of the features.
 *
 * For each GBWT path whose name passes the filter, call visit(name, occurrences)
 * with the exact places where the path visits the features, sorted by start.
 * The name leaves out any subrange or phase block, and the coordinates are on
 * the full path. For a haplotype path cut from an assembly, that is the
 * contig (see for_each_path_start()).
 *
 * Paths are searched with OpenMP threads. One thread at a time calls visit(),
 * in path order.
 */
void for_each_path_occurrences(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features,
                               const std::function<bool(const std::string&)>& include,
                               const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit);

/**
 * Find where each path visits at least min_match of the bases of each
 * feature, in order.
 *
 * This works like for_each_path_occurrences(), with the places found by
 * chain_anchors() in place of exact ones.
 */
void for_each_path_liftover(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features, double min_match,
                            const std::function<bool(const std::string&)>& include,
                            const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit);

}

#endif
