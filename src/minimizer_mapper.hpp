#ifndef VG_MINIMIZER_MAPPER_HPP_INCLUDED
#define VG_MINIMIZER_MAPPER_HPP_INCLUDED

/** 
 * \file minimizer_mapper.hpp
 * Defines a mapper that uses the minimizer index and GBWT-based extension.
 */

#include "aligner.hpp"
#include "xg.hpp"
#include "minimizer.hpp"
#include "alignment_emitter.hpp"
#include "gapless_extender.hpp"
#include "snarls.hpp"
#include "min_distance.hpp"
#include "seed_clusterer.hpp"
#include "tree_subgraph.hpp"
#include "algorithms/nearest_offsets_in_paths.hpp"

#include <structures/immutable_list.hpp>

namespace vg {

using namespace std;

class MinimizerMapper : public AlignerClient {
public:

    /**
     * Construct a new MinimizerMapper using the given indexes. The XG index can be nullptr,
     * as we only use it for correctness tracking.
     */

    MinimizerMapper(const GBWTGraph& graph, const MinimizerIndex& minimizer_index,
         MinimumDistanceIndex& distance_index, const XG* xg_index = nullptr);

    /**
     * Map the given read, and send output to the given AlignmentEmitter. May be run from any thread.
     * TODO: Can't be const because the clusterer's cluster_seeds isn't const.
     */
    void map(Alignment& aln, AlignmentEmitter& alignment_emitter);

    // Mapping settings.
    // TODO: document each

    /// Use all minimizers with at most hit_cap hits
    size_t hit_cap = 10;

    /// Ignore all minimizers with more than hard_hit_cap hits
    size_t hard_hit_cap = 300;

    /// Take minimizers between hit_cap and hard_hit_cap hits until this fraction
    /// of total score
    double minimizer_score_fraction = 0.6;

    /// How many clusters should we align?
    size_t max_extensions = 48;

    /// How many extended clusters should we align, max?
    size_t max_alignments = 8;

    //If a cluster's score is smaller than the best score of any cluster by more than
    //this much, then don't extend it
    double cluster_score_threshold = 0;

    //If the read coverage of a cluster is less than the best coverage of any cluster
    //by more than this much, don't extend it
    double cluster_coverage_threshold = 0;

    //If an extension set's score is smaller than the best 
    //extension's score by more than this much, don't align it
    double extension_set_score_threshold = 0;

    //If an extension's score is smaller than the best extension's score by
    //more than this much, don't align it
    int extension_score_threshold = 0;

    size_t max_multimaps = 1;
    size_t distance_limit = 1000;
    bool do_chaining = true;
    bool use_xdrop_for_tails = true;
    bool linear_tails = false;
    /// Use GBWT states from extensions to seed connectivity and tail searches.
    bool reuse_gbwt_states = false;
    string sample_name;
    string read_group;
    
    /// Track which internal work items came from which others during each
    /// stage of the mapping algorithm.
    bool track_provenance = false;

    /// Guess which seed hits are correct by location in the linear reference
    /// and track if/when their descendants make it through stages of the
    /// algorithm. Only works if track_provenance is true.
    bool track_correctness = false;
    
protected:
    // These are our indexes
    const XG* xg_index; // Can be nullptr; only needed for correctness tracking.
    const MinimizerIndex& minimizer_index;
    MinimumDistanceIndex& distance_index;

    /// This is our primary graph.
    const GBWTGraph& gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender;
    
    /// We have a clusterer
    SnarlSeedClusterer clusterer;
    
    /**
     * Estimate the score it may be possible to achieve using the given group of GaplessExtensions.
     * Supports single full-length extensions and groups that need chaining.
     * May reorder the input extended_seeds vector if it is not sorted in read space.
     * Is not always an overestimate of the actual score.
     */
    int estimate_extension_group_score(const Alignment& aln, vector<GaplessExtension>& extended_seeds) const;
    
    /**
     * Determine if a score estimate is significant enough to justify computing the real Alignment.
     * Returns true if it might win or affect mapping quality, and false otherwise.
     */
    bool score_is_significant(int score_estimate, int best_score, int second_best_score) const; 
    
    /**
     * Operating on the given input alignment, chain together the given
     * extended perfect-match seeds and produce an alignment into the given
     * output Alignment object.
     */
    void chain_extended_seeds(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& out) const; 
    
    /**
     * Find for each pair of extended seeds all the haplotype-consistent graph
     * paths against which the intervening read sequence needs to be aligned.
     *
     * Limits walks from each extended seed end to the longest detectable gap
     * plus the remaining to-be-alinged sequence, both computed using the read
     * length.
     *
     * extended_seeds must be sorted by read start position. Any extended seeds
     * that overlap in the read will be precluded from connecting.
     *
     * numeric_limits<size_t>::max() is used to store sufficiently long Paths
     * ending before sources (which cannot be reached from other extended
     * seeds) and starting after sinks (which cannot reach any other extended
     * seeds). Only sources and sinks have these "tail" paths.
     *
     * Tail paths are only calculated if the MinimizerMapper has linear_tails
     * set to true.
     */
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> find_connecting_paths(const vector<GaplessExtension>& extended_seeds,
        size_t read_length) const;
        
        
    /**
     * For gapless extensions that can't reach/be reached by anything in
     * connecting_paths, get all the trees defining tails off the specified
     * side of the gapless extension. Assumes that connecting_paths contains no
     * tail entries itself (i.e. that linear_tails was false on the
     * MinimizerMapper when find_connecting_paths computed it), but uses it to
     * identify the source or sink extensions.
     *
     * If the gapless extension starts or ends at a node boundary, there may be
     * multiple trees produced, each with a distinct root.
     *
     * If the gapless extension abuts the edge of the read, no forests will be
     * produced for it.
     *
     * Each tree is represented as a TreeSubgraph over our gbwt_graph.
     *
     * If left_tails is true, the trees read out of the left sides of the
     * gapless extensions. Otherwise they read out of the right sides.
     *
     * Gapless extensions that are not sources or sinks get no map entries.
     * Gapless extensions with dangling read sequence but no viable paths in
     * the graph will at least get map entries with empty forests. Sources or
     * sinks with no dangling sequence don't necessarily get map entries at
     * all.
     */
    unordered_map<size_t, vector<TreeSubgraph>> get_tail_forests(const vector<GaplessExtension>& extended_seeds,
        size_t read_length, const unordered_map<size_t, unordered_map<size_t, vector<Path>>>& connecting_paths, bool left_tails) const;
        
    /**
     * Find the best alignment of the given sequence against any of the paths
     * defined in paths.
     *
     * If no mapping is possible, produce a pure insert at default_position.
     *
     * If pinned is true, pin the alignment on one end to the start or end of
     * each path.
     *
     * When pinning, if pin_left is true, pin it on the left to the start of
     * each path. Otherwise pin it on the right to the end.
     *
     * Returns alingments in gbwt_graph space.
     */
    pair<Path, size_t> get_best_alignment_against_any_path(const vector<Path>& paths, const string& sequence,
        const Position& default_position, bool pinned, bool pin_left) const;
    
    /**
     * Find the best alignment of the given sequence against any of the trees
     * provided in trees, where each tree is a TreeSubgraph over the GBWT
     * graph. Each tree subgraph is rooted at the left in its own local
     * coordinate space, even if we are pinning on the right.
     *
     * If no mapping is possible (for example, because there are no trees),
     * produce a pure insert at default_position.
     *
     * Alignment is always pinned.
     *
     * If pin_left is true, pin the alignment on the left to the root of each
     * tree. Otherwise pin it on the right to the root of each tree.
     *
     * Returns alingments in gbwt_graph space.
     */
    pair<Path, size_t> get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees, const string& sequence,
        const Position& default_position, bool pin_left) const;
        
    /// We define a type for shared-tail lists of Mappings, to avoid constantly
    /// copying Path objects.
    using ImmutablePath = structures::ImmutableList<Mapping>;
    
    /**
     * Get the from length of an ImmutabelPath.
     *
     * Can't be called path_from_length or it will shadow the one for Paths
     * instead of overloading.
     */
    static size_t immutable_path_from_length(const ImmutablePath& path);
    
    /**
     * Convert an ImmutablePath to a Path.
     */
    static Path to_path(const ImmutablePath& path);

    /**
     * Given a Position, explore the GBWT graph out to the given maximum walk
     * distance.
     *
     * If from_state is not null, uses that starting GBWT search state, which
     * must be on the node that from is on and facking in the same orientation
     * as from.
     *
     * Calls the visit callback with the list of Mappings being extended (in
     * reverse order) and the handle it is being extended with.
     *
     * Only considers paths that visit at least one node after the node the
     * from Position is on. The from Position cuts immediately before the
     * first included base.
     *
     * If the callback returns false, that GBWT search state is not extended
     * further.
     *
     * If the walk_distance limit is exceeded, or a dead end in the graph is
     * hit, calls the limit callback with the list of Mappings (in reverse
     * order) that passed the limit or hit the dead end.
     */
    void explore_gbwt(const Position& from, const gbwt::SearchState* from_state, 
        size_t walk_distance, 
        const function<bool(const ImmutablePath&, const handle_t&)>& visit_callback,
        const function<void(const ImmutablePath&)>& limit_callback) const;
        
    /**
     * The same as explore_gbwt on a Position, but takes a handle in the
     * backing gbwt_graph and an offset from the start of the handle instead.
     */
    void explore_gbwt(handle_t from_handle, size_t from_offset, const gbwt::SearchState* from_state,
        size_t walk_distance,
        const function<bool(const ImmutablePath&, const handle_t&)>& visit_callback,
        const function<void(const ImmutablePath&)>& limit_callback) const;
    
    /**
     * Run a DFS on valid haplotypes in the GBWT starting from the given
     * Position, and continuing up to the given number of bases.
     *
     * Calls enter_handle when the DFS enters a haplotype visit to a particular
     * handle, and exit_handle when it exits a visit. These let the caller
     * maintain a stack and track the traversals.
     *
     * The starting node is only entered if its offset isn't equal to its
     * length (i.e. bases remain to be visited).
     *
     * Stopping early is not permitted.
     */
    void dfs_gbwt(const Position& from, size_t walk_distance,
        const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const;
     
    /**
     * The same as dfs_gbwt on a Position, but takes a handle in the
     * backing gbwt_graph and an offset from the start of the handle instead.
     */ 
    void dfs_gbwt(handle_t from_handle, size_t from_offset, size_t walk_distance,
        const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const;
     
};

}



#endif
