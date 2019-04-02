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
#include "distance.hpp"
#include "seed_clusterer.hpp"

namespace vg {

using namespace std;

class MinimizerMapper : public AlignerClient {
public:

    /**
     * Construct a new MinimizerMapper using the given indexes.
     */
    MinimizerMapper(const xg::XG* xg_index, const gbwt::GBWT* gbwt_index, const MinimizerIndex* minimizer_index,
        SnarlManager* snarl_manager, DistanceIndex* distance_index);

    /**
     * Map the given read, and send output to the given AlignmentEmitter. May be run from any thread.
     */
    void map(Alignment& aln, AlignmentEmitter& alignment_emitter);

    // Mapping settings.
    // TODO: document each
    size_t max_alignments;
    size_t max_multimaps;
    size_t hit_cap;
    size_t distance_limit;
    string sample_name;
    string read_group;


protected:
    // These are our indexes
    const xg::XG* xg_index;
    const gbwt::GBWT* gbwt_index;
    const MinimizerIndex* minimizer_index;
    SnarlManager* snarl_manager;
    DistanceIndex* distance_index;

    /// We have a GBWTGraph over the GBWT and the XG
    GBWTGraph gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender;
    
    /// We have a clusterer
    SnarlSeedClusterer clusterer;

    /**
     * Find for each pair of extended seeds all the haplotype-consistent graph
     * paths against which the intervening read sequence needs to be aligned.
     *
     * Limits walks to the longest detectable gap plus the remaining
     * to-be-alinged sequence, both computed using the read length.
     *
     * extended_seeds must be sorted by read start position. Any extended seeds
     * that overlap in the read will be precluded from connecting.
     *
     * numeric_limits<size_t>::max() is used to store sufficeintly long Paths
     * ending before sources (which cannot be reached from other extended
     * seeds) and starting after sinks (which cannot reach any other extended
     * seeds).
     *
     * Note that paths from all sinks are included, even if there would be no
     * read sequence to align against the path, because the read sequence
     * length is not passed.
     */
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> find_connecting_paths(const vector<pair<Path, size_t>>& extended_seeds,
        size_t read_length) const;

    /**
     * Given a Position, explore the GBWT graph out to the given maximum walk
     * distance.
     *
     * Calls the visit callback with the Path being extended and the handle it
     * is being extended with.
     *
     * Only considers paths that visit at least one node after the node the
     * from Position is on. The from Position cuts immediately before the
     * first included base.
     *
     * If the callback returns false, that GBWT search state is not extended
     * further.
     *
     * If the walk_distance limit is exceeded, or a dead end in the graph is
     * hit, calls the limit callback with the Path that passed the limit or hit
     * the dead end.
     */
    void explore_gbwt(const Position& from, size_t walk_distance, const function<bool(const Path&, const handle_t&)>& visit_callback,
        const function<void(const Path&)>& limit_callback) const;

};

}



#endif
