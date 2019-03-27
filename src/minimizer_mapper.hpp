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

    /// Find for each pair of extended seeds all the haplotype-consistent graph paths against which the intervening read sequence needs to be aligned.
    /// extended_seeds must be sorted by read start position. Any extended seeds that overlap in the read will be precluded from connecting.
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> find_connecting_paths(const vector<pair<Path, size_t>>& extended_seeds) const;


};

}



#endif
