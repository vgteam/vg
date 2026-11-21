#ifndef VG_FILL_PINS_HPP_INCLUDED
#define VG_FILL_PINS_HPP_INCLUDED

/** \file fill_pins.hpp
 *
 * Fill anchors: anchors at points that the graph alone chooses, spaced evenly inside runs of nodes
 * with no branch, so that the anchor file pins sequence that lies far from every site's pins.
 *
 * A *unitig* is a longest run of nodes in which each node's right side has one edge, to the next
 * node of the run, and that node's left side has that one edge only. The two nodes of each such
 * join must differ. A single node whose sides both branch, or have no edge, is a unitig of its own.
 * Inside each unitig of length L longer than the spacing K, ceil(L/K) - 1 *fill points* divide the
 * length evenly. A fill point is a pin inside a node: a point between two bases of one node.
 *
 * Each read that crosses a fill point is placed there, and a point with reads becomes a fill
 * anchor of one slot. The slot comes from the called walks of the innermost site that holds the
 * point's node and has a called genotype whose strands are known (see `FillAnchors::slot_of`).
 */

#include <array>
#include <atomic>
#include <cstdint>
#include <iosfwd>
#include <utility>
#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "../src/handle.hpp"
#include "../src/traversal_finder.hpp"
#include "anchor.hpp"
#include "site_read_source.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/// The width, in node IDs, of the windows in which `FillAnchors::resolve_reads` takes the points
/// where the read source has no windows of its own. It sets only how the work is shared among
/// threads.
constexpr size_t FILL_READ_WINDOW = 16384;

/// A fill point: a pin inside one node.
struct FillPoint {
    /// The node the point lies in.
    nid_t node = 0;
    /// The node's bases before the point, along its forward strand. It is at least 1 and at most
    /// the node's length minus 1.
    uint32_t pin = 0;

    /// Order by node, then pin.
    auto operator<=>(const FillPoint& other) const = default;
};

/// A tally of the fill points placed and dropped.
struct FillPointTally {
    /// The unitigs longer than the spacing, which get points.
    size_t unitigs = 0;
    /// The points placed.
    size_t placed = 0;
    /// The points dropped because they fell on a join between two nodes, and moving them forward
    /// off the join reached the end of the unitig.
    size_t dropped_at_end = 0;
    /// The points dropped because the point before them, moved forward off a join between two
    /// nodes, reached or passed them.
    size_t merged = 0;

    /// Add another tally to this one.
    void add(const FillPointTally& other);
};

/**
 * Return a unitig, given as a run of handles each followed by the next, in its canonical
 * direction.
 *
 * A unitig of one node reads the node forward. A longer one reads from the end whose node has the
 * smaller ID, so that the points it gets do not depend on which end it was found from.
 */
vector<handle_t> canonical_unitig(const HandleGraph& graph, vector<handle_t> unitig);

/**
 * Append the fill points of one unitig to `out`, in the order they lie along it.
 *
 * The unitig is read in the direction given. For a unitig of L bases, with L larger than
 * `spacing`, there are n = ceil(L/spacing) - 1 points, the i-th of them with round(i L / (n+1))
 * bases of the unitig before it, rounding half up. A point that falls on a join between two nodes
 * moves forward one base at a time until it lies inside a node. It is dropped if it reaches the end
 * of the unitig, and a later point that would move to the same place, or that the move passed, is
 * dropped too.
 */
void place_unitig_points(const HandleGraph& graph, const vector<handle_t>& unitig, size_t spacing,
                         vector<FillPoint>& out, FillPointTally& tally);

/**
 * Return the fill points of every unitig of `graph`, each read in its canonical direction, sorted
 * by node, then pin.
 *
 * Runs on several threads. The points depend on the graph and `spacing` alone.
 *
 * @param spacing The spacing in bases, or 0 for no points.
 * @param tally Set to what was placed and dropped.
 */
vector<FillPoint> place_fill_points(const HandleGraph& graph, size_t spacing,
                                    FillPointTally& tally);

/// The bits of a fill point's mark, which record which called walks cross the point.
///
/// The walk of strand 0 of a genotype with two different alleles.
constexpr uint8_t FILL_STRAND_0 = 1;
/// The walk of strand 1 of a genotype with two different alleles.
constexpr uint8_t FILL_STRAND_1 = 2;
/// The walk of an allele that every strand of the genotype carries, as in a homozygote, or in a
/// region of ploidy 1.
constexpr uint8_t FILL_EVERY_STRAND = 4;

/**
 * Return the walks of a called genotype that mark fill points, as pairs of an allele and the bits
 * its walk sets.
 *
 * A genotype of one allele at a nested site whose strand is not known gives no walks, since its
 * one allele lies on one strand of the region, and which one is not known. The walks of the sites
 * around it then decide its points.
 *
 * @param genotype The called alleles, in phase order.
 * @param nested_strand For a genotype of one allele at a nested haploid chain, the strand of the
 * parent that carries the chain, 0 or 1, or -1 where it is not known.
 * @param nested Whether the site is nested in another. A site that is not, with a genotype of one
 * allele, is in a region of ploidy 1, and its allele is on every strand of the region.
 */
vector<pair<int, uint8_t>> fill_marks_of_call(const vector<int>& genotype, int nested_strand,
                                              bool nested);

/// The kind and slot of a fill anchor.
struct FillSlot {
    AnchorKind kind = AnchorKind::fill_off_call;
    int slot = 0;
};

/**
 * Return the kind and slot of the fill anchor at a point.
 *
 * A point that no site holds, or that both strands' walks cross, is fill-shared, in slot 0. A
 * point that one strand's walk alone crosses is fill-strand, in that strand's slot. A point that
 * no called walk crosses is fill-off-call, in slot 0.
 *
 * @param held_by_root Whether the point's node lies in no site.
 * @param mark The bits that the called walks of the site deciding the point set on it.
 */
FillSlot classify_fill_point(bool held_by_root, uint8_t mark);

/// How a read's mapping quality gives its score at a fill anchor.
struct FillScoreParams {
    /// Whether the mismapping probability comes from the mapping quality. If not, every read has
    /// `min_mismap`.
    bool use_mapping_quality = true;
    /// The floor on the mismapping probability.
    double min_mismap = 0.0;
    /// The ceiling on the mismapping probability.
    double max_mismap = 1.0;
};

/**
 * Return the score of a read with mapping quality `mapq` at a fill anchor.
 *
 * This is the phred-scaled mismapping probability, clamped between the floor and the ceiling, which
 * is the score a read gets at a site of one allele. A probability of 0 gives 99.
 */
float fill_read_score(int32_t mapq, const FillScoreParams& params);

/// One read crossing one fill point.
struct FillPlacement {
    /// The point, as its index in `FillAnchors::points()`.
    uint32_t point = 0;
    /// The read's row in the point's anchor. Its score is from `fill_read_score`, and its
    /// direction is 0 when the read crosses the point along the node's forward strand, 1 against
    /// it.
    AnchorWriter::ReadRow row;
};

/**
 * The fill anchors of a run: their points, the reads that cross them, and the called walks that
 * cross them.
 *
 * Use it in this order. `place_points` places the points and groups them by window of node IDs.
 * While the sites are genotyped, a read source can call `resolve_window_reads` as its window
 * observer (see `SiteReadSource::set_window_observer`), which places the reads crossing the points
 * of each window the source fetches whole, while that window is in memory. `place` finds the sites
 * that hold each point's node. `resolve_reads` fetches the reads of the windows that
 * `resolve_window_reads` did not take, before the read source is closed, and gathers every
 * window's placements. `mark_walk`, called for each strand of each called site as its record is
 * rendered, records which strands cross each point. `add_anchors` then makes the anchors, and
 * `report` writes the counts.
 */
class FillAnchors {
public:
    FillAnchors() = default;
    FillAnchors(const FillAnchors&) = delete;
    FillAnchors& operator=(const FillAnchors&) = delete;

    /**
     * Set the spacing of the points, and how reads are scored.
     *
     * @param point_spacing The spacing in bases, or 0 for no fill anchors.
     * @param read_scores How a read's mapping quality gives its score.
     * @param min_score The smallest score a read must have to be placed.
     */
    void configure(size_t point_spacing, const FillScoreParams& read_scores, double min_score);

    /// Return whether fill anchors are made.
    bool is_enabled() const { return spacing > 0; }

    /**
     * Place the points in the unitigs of `point_graph`, and group them by window of `window` node
     * IDs.
     *
     * To take the windows of a read source, `window` must be the source's window width, since
     * `resolve_window_reads` takes only whole windows of this width. `point_graph` must outlive this
     * object.
     */
    void place_points(const HandleGraph& point_graph, size_t window);

    /**
     * Find the sites of `decomposition` that hold each point's node.
     *
     * Call it after `place_points`. `decomposition` must outlive this object.
     */
    void place(const SnarlDecomposition& decomposition);

    /**
     * Place the reads of `reads` at the points in the window of node IDs from `first` to `last`,
     * as `resolve_reads` would, unless that window's reads were placed already.
     *
     * It has the signature of a `SiteReadSource::WindowObserver`, so that the caller can have a
     * read source pass it each window the source fetches whole, and `resolve_reads` then fetches
     * only the windows it did not take. A range that is not exactly one window of `place_points`
     * is ignored. Safe to call from many threads.
     */
    void resolve_window_reads(nid_t first, nid_t last, const SiteReadSource::WindowReads& reads);

    /**
     * Find where each read of `reads` crosses each point, in the windows whose reads
     * `resolve_window_reads` did not place.
     *
     * The windows are taken on several threads, and the reads touching each window's points are
     * fetched without counting them (see
     * `SiteReadSource::for_each_read_uncounted`). A read is placed at a point by
     * `resolve_anchor_pin`, as an exit pin in the direction of the node's forward strand, so a
     * read that visits the node twice, or has no base aligned to the node base before the point,
     * is not placed. Nor is a read whose alignment does not reach the node base before the point,
     * or that has no read base after the point outside its soft clips.
     */
    void resolve_reads(const SiteReadSource& reads);

    /**
     * Record that `site` has a called genotype, and that the walk of its allele `walk` crosses the
     * points on the walk's nodes that `site` holds.
     *
     * The walk sets `bits` on each such point. At a point that several sites hold, only the marks
     * of the innermost of them with a called genotype count (see `slot_of`). The marks do not
     * depend on the order of the calls. Safe to call from many threads.
     *
     * @param bits Which strands the walk is, from `fill_marks_of_call`.
     */
    void mark_walk(const net_handle_t& site, const Traversal& walk, uint8_t bits);

    /// Turn each point with reads into a fill anchor, and add it to `anchors`.
    ///
    /// The placements are freed afterwards.
    void add_anchors(AnchorCollector& anchors);

    /// Write what was placed, resolved and classified.
    void report(ostream& out) const;

    /// Return the points, sorted by node, then pin.
    const vector<FillPoint>& points() const { return fill_points; }

    /// Return the reads placed at the points, grouped by point.
    const vector<FillPlacement>& placements() const { return placed_reads; }

    /**
     * Return the kind and slot of the fill anchor at point `point`, from the marks so far.
     *
     * The innermost site that holds the point's node and was marked decides, as
     * `classify_fill_point` says. A site that was not genotyped, or that is a nested chain whose
     * strand is not known, marks nothing, so the site around it decides. A point that sites hold,
     * none of them marked, is fill-off-call.
     */
    FillSlot slot_of(size_t point) const;

private:
    /// A node of the points that a read visits.
    struct VisitedNode {
        /// The node.
        nid_t node = 0;
        /// The index in the read's path of the read's first mapping on the node.
        int mapping = -1;
        /// How many of the read's mappings are on the node.
        size_t mappings = 0;
    };

    /// Set `visited` to the nodes of `nodes`, which is sorted, that `read` visits, once each and
    /// sorted.
    static void nodes_visited(const SiteRead& read, const vector<nid_t>& nodes,
                              vector<VisitedNode>& visited);

    /// The points in one window of node IDs.
    struct PointWindow {
        /// The window's number, its first node ID divided by `window_size`.
        nid_t number = 0;
        /// The points, as [begin, end) in `fill_points`.
        size_t begin = 0;
        size_t end = 0;
    };

    /// Return the placements of the reads of `reads` touching the points `fill_points[begin, end)`,
    /// which lie in one window of node IDs.
    vector<FillPlacement> resolve_window(const SiteReadSource::WindowReads& reads, size_t begin,
                                         size_t end);

    /// Place `read` at each of the points `fill_points[begin, end)` on the nodes `visited`, which
    /// it visits, and append the placements to `out`.
    void place_read(const SiteRead& read, const vector<VisitedNode>& visited, size_t begin,
                    size_t end, vector<FillPlacement>& out);

    /// The spacing in bases, or 0 for no fill anchors.
    size_t spacing = 0;
    /// How reads are scored.
    FillScoreParams scores;
    /// The smallest score a read must have to be placed.
    double min_read_score = 0.0;
    /// The graph the points are in, or null before `place_points`.
    const HandleGraph* graph = nullptr;
    /// The sites, or null before `place`.
    const SnarlDecomposition* sites = nullptr;

    /// The points, sorted by node, then pin.
    vector<FillPoint> fill_points;
    /// The width of the windows of node IDs the points are grouped by, or 0 before
    /// `place_points`.
    size_t window_size = 0;
    /// The windows that hold points, in order.
    vector<PointWindow> windows;
    /// Whether the reads of each window of `windows` have been placed, or are being placed.
    vector<std::atomic<bool>> window_claimed;
    /// The placements in each window of `windows`, until `resolve_reads` gathers them.
    vector<vector<FillPlacement>> window_placements;
    /// The windows whose reads `resolve_window_reads` placed.
    std::atomic<size_t> windows_observed{0};
    /// The windows whose reads `resolve_reads` fetched.
    size_t windows_fetched = 0;
    /// The sites that hold some point's node, each as its canonical net handle, sorted.
    vector<net_handle_t> holding_sites;
    /// Whether each site of `holding_sites` was marked by `mark_walk`.
    vector<std::atomic<uint8_t>> site_marked;
    /// The sites that hold each point's node, innermost first, as indices into `holding_sites`.
    /// Those of point i are entries `holder_start[i]` to `holder_start[i + 1]` of `holders`. A
    /// point with none lies in no site.
    vector<uint32_t> holders;
    /// Where the sites of each point start in `holders`, with one more entry at the end.
    vector<uint32_t> holder_start;
    /// The mark of each point: the bits the called walks set on it, in the low byte, from the
    /// innermost site to mark it, whose depth among the point's sites, counting the innermost as 1,
    /// is the high byte. 0 for a point no walk has marked.
    vector<std::atomic<uint16_t>> marks;
    /// The reads placed at the points, grouped by point once `resolve_reads` is done.
    vector<FillPlacement> placed_reads;

    /// What placing the points did.
    FillPointTally point_tally;
    /// How the reads were resolved at the points.
    AnchorCounters read_counters;
    /// The reads whose alignment on a point's node does not reach the node base before the point.
    std::atomic<size_t> not_reached{0};
    /// The reads that stop at a point, or go on past it only in a soft clip.
    std::atomic<size_t> stopped{0};
    /// The reads resolved at a point.
    std::atomic<size_t> resolved{0};
    /// The resolved reads left out for a score below `min_read_score`.
    std::atomic<size_t> low_score{0};
    /// The points with reads of each fill kind, from `add_anchors`: fill-shared, fill-strand and
    /// fill-off-call.
    std::array<size_t, 3> kind_counts{};
};

}

#endif
