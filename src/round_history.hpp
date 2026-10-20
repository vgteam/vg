#ifndef VG_ROUND_HISTORY_HPP_INCLUDED
#define VG_ROUND_HISTORY_HPP_INCLUDED

#include <array>
#include <unordered_map>
#include <vector>

#include "linkage_model.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A log of the states that the re-genotyping rounds reach, for spotting convergence or a cycle.
 *
 * A state is the genotype chosen at every staged site. Rounds can cycle without converging,
 * because dropping and reinstating a subtree is a discrete change and the phase between sites
 * changes with the genotypes, so no single score has to improve each round.
 */
class RoundHistory {
public:
    /// Each staged site's chosen pair and ploidy, `{trav_first, trav_second, ploidy}`, keyed by
    /// record key. A site the model chose nothing for is left out.
    using State = unordered_map<size_t, std::array<int, 3>>;

    /// Return the state `model` has chosen for `sites`, or an empty state if `model` is null.
    static State state(StagedSiteTable& sites, const LinkageCollector* model);

    /// Count the sites whose chosen pair or ploidy differs between `before` and `after`.
    ///
    /// A site chosen in only one of the two states counts as changed.
    static size_t changed(const State& before, const State& after);

    /// Remember `state` as the next round's.
    void remember(const State& state);

    /// Return the first remembered round that reached `state`, counting from 1, or 0 if none did.
    ///
    /// Rounds are compared by digest, so a digest collision reports a match that is not one.
    size_t first_round_with(const State& state) const;

    /// Return how many rounds are remembered.
    size_t rounds() const { return digests.size(); }

private:
    /// Hash a state, visiting its sites in the map's iteration order.
    static size_t digest(const State& state);

    /// Each remembered round's digest, in order.
    vector<size_t> digests;
};

}

#endif
