#ifndef VG_COVER_LIFT_HPP_INCLUDED
#define VG_COVER_LIFT_HPP_INCLUDED

/**
 * \file cover_lift.hpp
 *
 * Lifting features from a cover to the graph, the gref and the paths.
 *
 * A cover (see unitigs.hpp) is a set of sequences that a linear tool can scan
 * in place of the haplotypes of a graph. It is made of the gref paths and the
 * tigs of vg unitigs. Each cover sequence spells a walk in the graph, so a
 * feature that the tool reports on one, as a BED record, is a WalkInterval.
 * The interval has a gref placement. Each path that visits the interval gives
 * the feature a place on that path's assembly.
 *
 * Every window of at most k bases of a haplotype lies within some cover
 * sequence. So a tool that looks at no more than k bases at a time finds each
 * feature on every haplotype. Cover sequences overlap by fewer than k bases,
 * so a window shorter than k can lie within several of them. Its copies lift
 * to the same interval and are kept once.
 */

#include "unitigs.hpp"

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
 * A feature in cover coordinates.
 *
 * It covers bases [start, end) of the named cover sequence. It reads them in
 * reverse if its strand is "-".
 */
struct CoverFeature {
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
std::vector<CoverFeature> read_bed_features(std::istream& in);

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

    /// Return true if the feature has a strand.
    bool is_stranded() const;
};

/**
 * A directory of the sequences of a cover.
 *
 * It knows the walk that each gref path and each loaded tig spells. It
 * translates features on those sequences into graph features.
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

    /**
     * Translate the feature into a graph feature.
     *
     * A stranded feature is oriented along its strand. A feature without a
     * strand takes the orientation with the most bases on forward gref
     * strands. On a tie, it takes the orientation whose walk sorts first.
     *
     * @pre The feature's sequence is in the cover, and its range is within it.
     */
    GraphFeature lift(const CoverFeature& feature) const;

private:
    const gbwtgraph::GBWTGraph& graph;
    /// Spells the gref paths and places intervals on the gref. Null without a gref.
    const GrefIndex* gref;
    /// The walk interval of each loaded tig, by name.
    std::unordered_map<std::string, WalkInterval> tigs;
};

/**
 * Sort the features along the gref, then by interval and fields, and merge the
 * copies of each feature.
 *
 * Copies have the same interval and the same fields.
 */
void sort_and_merge_features(std::vector<GraphFeature>& features);

/**
 * A matcher that finds where paths visit the intervals of graph features.
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

/// A place where a path visits a feature.
struct PathOccurrence {
    /// Coordinate on the full path of the first base of the feature.
    uint64_t start;
    /// Index of the feature.
    size_t feature;
    /// True if the path reads the feature in reverse.
    bool is_reverse;
};

/**
 * Find where each path visits the features.
 *
 * For each GBWT path whose name passes the filter, call visit(name, occurrences)
 * with the places where the path visits the features, sorted by start. The name
 * leaves out any subrange or phase block, and the coordinates are on the full
 * path. For a haplotype path cut from an assembly, that is the contig (see
 * for_each_path_start()).
 *
 * Paths are searched with OpenMP threads. One thread at a time calls visit(),
 * in path order.
 */
void for_each_path_occurrences(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features,
                               const std::function<bool(const std::string&)>& include,
                               const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit);

}

#endif
