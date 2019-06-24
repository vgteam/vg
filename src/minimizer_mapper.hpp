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

#include <structures/immutable_list.hpp>

namespace vg {

using namespace std;

class MinimizerMapper : public AlignerClient {
public:

    /**
     * Construct a new MinimizerMapper using the given indexes.
     */

    MinimizerMapper(const XG* xg_index, const gbwt::GBWT* gbwt_index, const MinimizerIndex* minimizer_index,
         MinimumDistanceIndex* distance_index);

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

    size_t max_multimaps = 1;
    size_t distance_limit = 1000;
    bool do_chaining = true;
    bool all_tails = false;
    bool use_xdrop_for_tails = false;
    string sample_name;
    string read_group;


protected:
    // These are our indexes
    const XG* xg_index;
    const gbwt::GBWT* gbwt_index;
    const MinimizerIndex* minimizer_index;
    MinimumDistanceIndex* distance_index;

    /// We have a GBWTGraph over the GBWT and the XG
    GBWTGraph gbwt_graph;
    
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
     */
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> find_connecting_paths(const vector<GaplessExtension>& extended_seeds,
        size_t read_length) const;
        
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
    void explore_gbwt(const Position& from, size_t walk_distance, const function<bool(const ImmutablePath&, const handle_t&)>& visit_callback,
        const function<void(const ImmutablePath&)>& limit_callback) const;
     
};

}



#endif
