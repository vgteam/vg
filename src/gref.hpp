#ifndef VG_GREF_HPP_INCLUDED
#define VG_GREF_HPP_INCLUDED

/**
 * \file gref.hpp
 *
 * Interface for computing and querying graph reference path covers.
 *
 * An gref cover is a set of path fragments (stored as separate paths) in the graph.
 * They are always relative to an existing reference sample (ie GRCh38 or CHM13).
 * Unlike rGFA paths which use complex metadata embedding, gref paths use a simple naming
 * scheme: the base reference path name moved into the "gref_" namespace, with the
 * fragments numbered off it as {gref_base_path_name}_{N}_alt
 *
 * For example, if the reference path is "CHM13#0#chr1", the cover is written as:
 *   - gref_CHM13#0#chr1        (a copy of the reference path itself)
 *   - gref_CHM13#0#chr1_1_alt  (fragments hanging off it)
 *   - gref_CHM13#0#chr1_2_alt
 *   - etc.
 *
 * The naming is a convention, not an option: because the prefix lands on the sample,
 * a graph carries both references at once, "is this a gref path/sample" is a prefix
 * test, and the base path a gref path came from is recoverable by dropping the prefix.
 *
 * The data structures used in this class are always relative to the original paths in the
 * graph.  The REFERENCE-sense fragments the cover is written as are output only; nothing
 * reads them back.
 */

#include <bdsg/snarl_distance_index.hpp>

#include "handle.hpp"

namespace vg {

using namespace std;

class GrefCover {
public:

    // The prefix every gref path name carries
    static const string gref_prefix;  // "gref_"

    // The suffix used to identify gref fragments
    static const string gref_suffix;  // "_alt"

    // Name of the gref copy of a path: the same name moved into the gref namespace,
    // keeping its subrange so that subpaths of one contig stay distinct and keep their
    // coordinates.  Example: make_gref_copy_name("CHM13#0#chr1[100-200]")
    //                        -> "gref_CHM13#0#chr1[100-200]"
    // For PanSN names the prefix renames the sample, which is what makes "is this path
    // a gref path" and "which base path is it derived from" answerable from the name
    // alone.  Anything matching a base path against its gref copy must use this.
    static string make_gref_copy_name(const string& path_name);

    // Name a fragment's gref base is built from: make_gref_copy_name() with the subrange
    // dropped as well, since make_gref_name() appends "_{N}_alt" and that has to land on
    // the locus.  Fragments off different subpaths of one contig therefore share a base
    // name; the per-base counter keeps their full names distinct.
    // Example: make_gref_base_name("CHM13#0#chr1[100-200]") -> "gref_CHM13#0#chr1"
    static string make_gref_base_name(const string& base_path_name);

    // Test if a name is in the gref namespace. Works on path names and, because
    // the prefix lands on the sample for PanSN names, on sample names too.
    // Example: is_gref_derived("gref_CHM13#0#chr1_3_alt") -> true
    //          is_gref_derived("gref_CHM13") -> true
    static bool is_gref_derived(const string& name);


    // Create an gref path name from a base reference path name and an index.
    // Example: make_gref_name("gref_CHM13#0#chr1", 1) -> "gref_CHM13#0#chr1_1_alt"
    static string make_gref_name(const string& base_path_name, int64_t gref_index);

    // Test if a path name is an gref fragment (contains "_{N}_alt" suffix).
    // Copied base contigs are gref paths but not fragments, so they fail this.
    static bool is_gref_name(const string& path_name);

    // Parse an gref path name to extract the base reference path name.
    // Returns the original base path name, or the input if not an gref path.
    // Example: parse_base_path("gref_CHM13#0#chr1_3_alt") -> "gref_CHM13#0#chr1"
    static string parse_base_path(const string& gref_name);

    // Parse an gref path name to extract the gref index.
    // Returns -1 if the path is not an gref path.
    // Example: parse_gref_index("gref_CHM13#0#chr1_3_alt") -> 3
    static int64_t parse_gref_index(const string& gref_name);

public:
    // Clear out any existing gref paths from the graph. Recommended to run this
    // before compute().
    void clear(MutablePathMutableHandleGraph* graph);

    // Compute the gref cover from the graph, starting with a given set of reference paths.
    //
    // distance_index, when non-null, is a decomposition-only snarl distance index over the
    // same graph.  It is not used to build the cover; it is used afterwards to map each
    // fragment to the top-level snarl containing it (see assign_top_level_snarls()).
    void compute(const PathHandleGraph* graph,
                 const bdsg::SnarlDistanceIndex* distance_index,
                 const unordered_set<path_handle_t>& reference_paths,
                 int64_t minimum_length);


    // Apply the gref cover to a graph (must have been computed first).  Everything
    // lands in the gref namespace:
    // 1. The base reference paths are copied over (CHM13#0#chr1 -> gref_CHM13#0#chr1)
    // 2. The fragments are added alongside them (gref_CHM13#0#chr1_1_alt, etc.)
    // The original paths are left untouched, so the graph carries both the base
    // reference and the gref reference and they can be told apart by name.
    void apply(MutablePathMutableHandleGraph* mutable_graph);

    // Enable verbose output (coverage summary, etc.)
    void set_verbose(bool verbose);


    // Summary of how the final cover sits inside the top-level snarl decomposition.
    // Filled in by assign_top_level_snarls(); all zero when no distance index was supplied.
    //
    // The unit of the decomposition is a child of a top-level chain: either a top-level
    // snarl (a bubble hanging off the spine) or a bare spine node.  A fragment is
    // independently processable exactly when all of its nodes sit in one unit.
    struct TopLevelSnarlStats {
        // Shape of the decomposition
        int64_t top_level_chains = 0;
        int64_t top_level_snarls = 0;
        int64_t spine_nodes = 0;             // node children of top-level chains
        int64_t unexpected_chain_children = 0;  // neither snarl nor node: should be 0

        // Where the fragments landed
        int64_t fragments = 0;
        int64_t fragments_in_one_snarl = 0;     // wholly inside a single top-level snarl
        int64_t fragments_on_one_spine_node = 0;
        int64_t fragments_spanning_snarls = 0;  // touch >= 2 top-level snarls: THE CLAIM
        int64_t fragments_touching_spine = 0;   // contain at least one spine node
        int64_t fragments_spanning_units = 0;   // touch >= 2 units of any kind
        int64_t fragments_with_unmapped_nodes = 0;  // decomposition did not place a node
        int64_t max_snarls_in_one_fragment = 0;
        int64_t max_units_in_one_fragment = 0;
        // How many distinct top-level snarls the fragments landed in.  Equal to
        // fragments_in_one_snarl when no two fragments share a snarl; far below it would
        // mean the decomposition is coarse, and 1 would mean the mapping is degenerate.
        int64_t distinct_snarls_with_fragments = 0;
        // Nodes the descent placed in two different units.  The decomposition is a
        // partition, so this must be 0; a non-zero value would mean the unit assignment is
        // ambiguous and every count above it is meaningless.
        int64_t nodes_claimed_twice = 0;
    };

    // Require every top-level snarl to be anchored on the reference, and take the sequence
    // of any that is not out of the cover entirely.
    //
    // Anchored means both bounds of the snarl are nodes of the SAME reference contig.  Both
    // halves matter.  Reference nodes are pre-assigned before any candidate is considered, so
    // a fragment can never contain one; a snarl with both bounds on one reference contig is
    // therefore uncrossable, and everything inside it is reachable only through those two
    // nodes.  That is what makes each top-level snarl an independent problem, and it is why
    // no fragment can ever need to be merged with one from another snarl.  Requiring the SAME
    // contig, not merely two reference nodes, is what rules out a snarl bridging two contigs,
    // which is exactly the case where "inside" stops being well defined.
    //
    // A snarl that fails is not repaired and is not covered: its nodes are withheld from the
    // cover and reported.  Withholding sequence is a real cost, so it is counted in nodes and
    // bp and warned about unconditionally -- never silently.  On a graph whose reference is
    // full-length nothing is withheld; the failures are clipped or subranged references, where
    // the spine is not continuous and the guarantee genuinely does not hold.
    //
    // Must run after the reference paths are pre-assigned and before any fragment is claimed.
    void enforce_top_level_anchoring(const bdsg::SnarlDistanceIndex& distance_index);

    // Map every fragment onto the top-level snarl decomposition carried by distance_index.
    // Fills interval_snarl_bounds with the boundary nodes of the containing top-level snarl
    // for every fragment that lies wholly inside one, leaves {0,0} for the rest, and records
    // the summary returned by get_top_level_snarl_stats().
    void assign_top_level_snarls(const bdsg::SnarlDistanceIndex& distance_index);


    // Work out, for every fragment, which gref contig supplies the coordinates its VCF records
    // will be reported against, and how many such contigs lie between it and the base
    // reference.  Fills interval_parent and interval_level.
    //
    // The rule follows the caller rather than the graph.  Walk up the snarl tree from the
    // fragment until reaching a snarl whose boundary nodes are owned by some gref interval;
    // that interval is the parent, because a snarl's records are reported against whichever
    // gref path traverses its boundaries.  A fragment whose enclosing snarl is bounded by
    // reference nodes is level 1 -- one coordinate-system change from the base reference --
    // and the levels chain from there.
    //
    // This is the measurement e20e1f277 reverted a column for getting wrong.  That version
    // counted hops over interval adjacency, which is a fact about the graph, while INFO/CH
    // counts coordinate-system changes among ancestors, which is a fact about the caller; they
    // agreed on only 76% of chr22 fragments.  The difference then was that the cover had no
    // snarl decomposition and had to approximate containment.  It has one now, and it is the
    // same decomposition deconstruct recurses, so this computes the caller's answer directly
    // instead of a proxy for it.  Gate any change here on agreement with INFO/CH, not on
    // whether the numbers look plausible.
    void assign_nesting(const bdsg::SnarlDistanceIndex& distance_index);

    // Write a tab-separated table describing gref segments.
    // Each line contains: source_path, source_start, source_end, gref_contig, level, strand,
    //                     ref_contig, ref_start, ref_end, top_level_snarl
    //
    // Columns 1-6 are a valid BED6, with the nesting level in the score slot -- a small
    // non-negative integer is what that field wants, and it avoids a column that would always
    // be ".".  `cut -f1-6` therefore feeds bedtools directly, strand included.
    //
    // The first line is a '#'-prefixed header naming the columns.  Consumers skip it by the
    // usual conventions (grep -v '^#', pandas comment='#').  Cactus's merge_gref_segs() keeps
    // the first header and drops the rest when concatenating the per-chromosome tables.
    //
    // ref_path is derivable from gref_path_name and is emitted anyway.  Recovering it means
    // finding the last '_' before the '_alt' suffix, and reference contig names contain
    // underscores (gref_GRCh38#0#chr4_GL000008v2_random_1_alt), so leaving it to the consumer
    // invites an off-by-one-underscore bug for no saving worth having.
    //
    // top_level_snarl names the site the fragment is an allele of, spelled ">start>end" as vg
    // deconstruct and vg call spell that snarl's VCF record ID, so it joins directly to an ID
    // or to the PS of a record nested inside the fragment.  Both are reference nodes on ref_path, so fragments that
    // are alleles of the same site share the pair, and grouping the table on it groups the
    // cover by site.  0 0 when no decomposition was supplied.
    //
    // level equals INFO/CH of every record on the contig -- NOT INFO/LV, which counts only
    // ancestors on the same CHROM and is a different number.
    //
    // There is deliberately no parent column.  It restated "gref_" + ref_contig on every one
    // of the 201,630 level-1 fragments genome-wide, and on the rest it named a contig without
    // a position on it, so nothing could be located or reproduced from it.  The VCF answers
    // the question properly: a record's PS names the enclosing site, and its CHROM names the
    // contig, with a position attached to both.
    //
    // strand is '-' when apply() reverse-complemented the run, so the sequence at
    // source_path[source_start:source_end] is the reverse complement of the gref contig.  26
    // of 2137 chr22 fragments; without it a samtools faidx on columns 1-3 silently disagrees
    // with the VCF's REF/ALT.
    //
    // Note ref_start/ref_end is a different thing from the site in columns 8-9: it is the
    // reference interval
    // the fragment's own source path brackets it with.  A top-level snarl is much coarser --
    // 77% of chrY's reference nodes are inside one -- so it identifies the site but does not
    // localise the fragment on the reference.
    //
    // Must be called after compute() and knows what gref path names will be used.
    void write_gref_segments(ostream& os);

protected:

    // Record an interval and assign it every node it walks.  Nothing merges or extends:
    // candidates are maximal runs of uncovered nodes, so a new interval can neither abut nor
    // overlap an existing one.  See the body for the measurement that established this.
    void add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_gref_intervals,
                      unordered_map<nid_t, int64_t>& thread_node_to_interval,
                      const pair<step_handle_t, step_handle_t>& new_interval,
                      vector<pair<nid_t, nid_t>>* snarl_bounds_vec = nullptr,
                      pair<nid_t, nid_t> snarl_bounds = {0, 0});


    // Remove non-reference intervals shorter than minimum_length, then defragment.
    void filter_short_intervals(int64_t minimum_length);

    // Build the cover: enumerate the maximal runs of uncovered nodes along every candidate
    // source path, then claim them longest-first.  This is the whole cover.
    void fill_uncovered_nodes(int64_t minimum_length);

    // Search back to the reference and return <distance, node_id> when found.
    // (here distance is the number of intervals crossed, aka rank)
    // "first" toggles returning the first interval found vs all of them.
    vector<pair<int64_t, nid_t>> get_reference_nodes(nid_t node_id, bool first) const;

    // Resolve the gref base path name a fragment interval gets named after: the reference
    // path it traces back to, or its own source path when there is no reference in its
    // component.
    //
    // apply() and write_gref_segments() must agree here.  Both number fragments by
    // incrementing a per-base-name counter over the same interval sequence, and
    // write_gref_segments() runs first, so it cannot read the names apply() picks -- it has
    // to predict them.  A single disagreement does not just mislabel one row: it renumbers
    // every later fragment on both of the names involved.  So the prediction is not a
    // reimplementation, it is this same call.
    //
    // When the interval traces back to a reference, *out_ref_path and *out_ref_node are set
    // to the reference path and the reference node the trace landed on.  When it does not,
    // both are left untouched, so initialise out_ref_node to 0 to detect that case.
    string resolve_base_path_name(int64_t interval_index,
                                  path_handle_t* out_ref_path = nullptr,
                                  nid_t* out_ref_node = nullptr) const;

    // Debug function: verify that every node in the graph is covered by the gref cover.
    // Prints a summary of coverage statistics to stderr.
    void verify_cover(int64_t minimum_length) const;

    // Enforce node-disjointness: no two gref intervals may share a node, and no fragment
    // may claim a node belonging to a reference interval.
    //
    // This is a hard invariant, not a quality measure.  Every gref path becomes its own
    // reference contig with its own coordinate space, so a node owned by two of them is the
    // same sequence at two different coordinates -- downstream that is not a worse cover,
    // it is an incoherent one.  node_to_interval cannot enforce it on its own: it maps each
    // node to ONE interval index, so an interval that overwrites another's entry silently
    // wins the index while the loser keeps the step in its own range and still gets written
    // by apply().  Only walking the intervals themselves can detect that, which is what
    // this does.
    //
    // Called at the end of compute() over the final cover, so it holds whatever route the
    // intervals took to get there.  Exits non-zero on violation rather than warning: a
    // cover that violates this should never reach a downstream tool.
    void verify_disjoint() const;

    const PathHandleGraph* graph = nullptr;

    // Intervals are end-exclusive (like BED).
    vector<pair<step_handle_t, step_handle_t>> gref_intervals;

    // Top-level snarl boundary nodes for each interval, parallel to gref_intervals.
    // (0, 0) sentinel for reference intervals, and for any fragment that is not wholly
    // contained in one top-level snarl (or that was computed without a distance index).
    vector<pair<nid_t, nid_t>> interval_snarl_bounds;

    // Summary from assign_top_level_snarls().
    TopLevelSnarlStats top_level_snarl_stats;

    // Nodes inside a top-level unit that failed the reference-anchoring requirement, and so
    // are withheld from the cover.  Empty unless enforce_top_level_anchoring() ran.
    // Parallel to gref_intervals, filled by assign_nesting().
    // interval_parent: index of the gref interval supplying this one's coordinates, or -1 for
    //   a reference interval and for a fragment with no gref ancestor at all.
    // interval_level: 0 for a reference interval, 1 for a fragment hanging directly off the
    //   base reference, n+1 for one inside a level-n fragment.
    vector<int64_t> interval_parent;
    vector<int64_t> interval_level;

    unordered_set<nid_t> excluded_nodes;
    int64_t excluded_units = 0;
    int64_t excluded_bp = 0;

    // gref_intervals[0, num_ref_intervals-1] are all rank-0 reference intervals.
    int64_t num_ref_intervals = 0;

    // Map from node ID to interval index.
    unordered_map<nid_t, int64_t> node_to_interval;

    // Counter for generating unique gref indices per base path.
    // Using mutable so it can be updated in apply() which is logically const for the cover.
    unordered_map<string, int64_t> base_path_gref_counter;

    // Whether to print verbose output (coverage summary, etc.)
    bool verbose = false;

    // Copy the base reference paths into the gref namespace.
    // Creates new paths like "gref_CHM13#0#chr1" from "CHM13#0#chr1".
    void copy_base_paths_to_gref(MutablePathMutableHandleGraph* mutable_graph,
                                 const unordered_set<path_handle_t>& reference_paths);

};

}

#endif
