#ifndef VG_RYMER_SEED_RESCUE_HPP_INCLUDED
#define VG_RYMER_SEED_RESCUE_HPP_INCLUDED

/** \file
 * SAFARI's way of finding seeds that damage would hide from an ordinary minimizer index.
 *
 * For each minimizer of the read under the RYmer index's k and w, the RYmer index gives the graph
 * k-mers with the same RYmer. A graph k-mer is accepted if `damage_posterior()` says the read
 * k-mer is likely enough to be a damaged copy of it. For each minimizer with an accepted graph
 * k-mer, the first one accepted in index order, the rescue makes a stand-in sequence: the length
 * of the read, all 'N' except for the graph k-mer at the read k-mer's offset. The mapper then
 * finds minimizers in each stand-in sequence with its ordinary minimizer index and seeds from them
 * as if they came from the read, so a seed exists even where damage changed the read's own k-mer.
 *
 * See: Rubin J, van Waaij J, Kraft L, Sirén J, Sackett PW, Renaud G. SAFARI: pangenome alignment
 * of ancient DNA using purine/pyrimidine encodings. NAR Genomics and Bioinformatics 7(4), lqaf170
 * (2025). https://doi.org/10.1093/nargab/lqaf170
 */

#include <ostream>
#include <string>
#include <vector>

#include <gbwtgraph/minimizer.h>

#include "deamination_profile.hpp"

namespace vg {

using namespace std;

class RymerSeedRescue {
public:
    /// The rescue keeps references to the index and the profile, which must outlive it. Throws
    /// std::invalid_argument if the index is not an RYmer index or a setting is out of range.
    RymerSeedRescue(const gbwtgraph::DefaultMinimizerIndex& rymer_index,
                    const DeaminationProfile& profile,
                    double posterior_threshold,
                    const DamageModelParameters& model);

    /// Returns the stand-in sequences, as the file comment describes, for a read with this
    /// sequence. Safe to call from several threads at once.
    vector<string> stand_in_sequences(const string& read) const;

private:
    /// A graph k-mer is accepted if its posterior is greater than this, or greater than 1e-10
    /// where this is 0.
    double posterior_threshold;

    /// Settings of the damage model.
    DamageModelParameters model;

    const gbwtgraph::DefaultMinimizerIndex& rymer_index;
    const DeaminationProfile& profile;
};

/// Prints the request to cite the SAFARI paper, with the citation and its BibTeX entry.
void print_safari_citation(ostream& out);

}

#endif
