#include "rymer_seed_rescue.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include "rymer_index.hpp"
#include <gbwtgraph/utils.h>

namespace vg {

using namespace std;

RymerSeedRescue::RymerSeedRescue(const gbwtgraph::DefaultMinimizerIndex& rymer_index,
                                 const DeaminationProfile& profile,
                                 double posterior_threshold,
                                 const DamageModelParameters& model) :
    posterior_threshold(posterior_threshold), model(model), rymer_index(rymer_index), profile(profile) {
    if (!is_rymer_index(rymer_index)) {
        throw invalid_argument("RymerSeedRescue needs an RYmer index, which index_rymers() builds");
    }
    if (rymer_index.uses_syncmers()) {
        throw invalid_argument("RymerSeedRescue cannot use an RYmer index of syncmers");
    }
    if (posterior_threshold < 0.0 || posterior_threshold > 1.0) {
        throw invalid_argument("RymerSeedRescue needs a posterior threshold from 0 to 1");
    }
    if (model.spurious_prior < 0.0 || model.spurious_prior > 1.0) {
        throw invalid_argument("RymerSeedRescue needs a prior probability of a spurious match from 0 to 1");
    }
}

vector<string> RymerSeedRescue::stand_in_sequences(const string& read) const {
    vector<string> result;
    size_t k = rymer_index.k();
    // A threshold of 0 would accept graph k-mers the damage model rules out, so it means "any
    // positive posterior" instead.
    double threshold = posterior_threshold == 0.0 ? 1e-10 : posterior_threshold;

    for (const auto& minimizer : rymer_index.minimizers(read)) {
        if (minimizer.empty()) {
            continue;
        }
        // The offset of a reverse minimizer is its last base in the read.
        size_t read_start = minimizer.is_reverse ? minimizer.offset - (k - 1) : minimizer.offset;
        string read_kmer = read.substr(read_start, k);

        // The same graph k-mer is stored once for each place it occurs, so test each one once,
        // keeping the index order so that the first accepted is the one SAFARI accepts.
        unordered_set<gbwtgraph::Key64::code_type> tested;
        for (gbwtgraph::Key64 graph_key : graph_kmers_for_rymer(rymer_index, to_rymer_minimizer(minimizer, k))) {
            if (!tested.insert(graph_key.get_key()).second) {
                continue;
            }
            // The graph k-mer is stored in the orientation chosen as a minimizer in the graph,
            // and that orientation has the same RYmer as the read's minimizer, so it is in
            // read orientation exactly when the read's minimizer is.
            string graph_kmer = graph_key.decode(k);
            if (minimizer.is_reverse) {
                graph_kmer = gbwtgraph::reverse_complement(graph_kmer);
            }
            if (damage_posterior(graph_kmer, read_kmer, read.size(), read_start, profile, model) > threshold) {
                string stand_in(read.size(), 'N');
                stand_in.replace(read_start, k, graph_kmer);
                result.push_back(std::move(stand_in));
                break;
            }
        }
    }
    return result;
}

void print_safari_citation(ostream& out) {
    out << "If you use SAFARI (vg safari or vg rymer) in published work, please cite:" << endl
        << endl
        << "  Joshua Rubin, Jan van Waaij, Louis Kraft, Jouni Sirén, Peter Wad Sackett, Gabriel Renaud." << endl
        << "  SAFARI: pangenome alignment of ancient DNA using purine/pyrimidine encodings." << endl
        << "  NAR Genomics and Bioinformatics, Volume 7, Issue 4, December 2025, lqaf170." << endl
        << "  https://doi.org/10.1093/nargab/lqaf170" << endl
        << endl
        << "BibTeX:" << endl
        << endl
        << "  @article{rubin2025safari," << endl
        << "    title   = {{SAFARI}: pangenome alignment of ancient {DNA} using purine/pyrimidine encodings}," << endl
        << "    author  = {Rubin, Joshua and van Waaij, Jan and Kraft, Louis and Sir{\\'e}n, Jouni and" << endl
        << "               Sackett, Peter Wad and Renaud, Gabriel}," << endl
        << "    journal = {NAR Genomics and Bioinformatics}," << endl
        << "    volume  = {7}," << endl
        << "    number  = {4}," << endl
        << "    pages   = {lqaf170}," << endl
        << "    year    = {2025}," << endl
        << "    month   = dec," << endl
        << "    doi     = {10.1093/nargab/lqaf170}," << endl
        << "    url     = {https://academic.oup.com/nargab/article/7/4/lqaf170/8376687}" << endl
        << "  }" << endl;
}

}
