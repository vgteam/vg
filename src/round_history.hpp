#ifndef VG_ROUND_HISTORY_HPP_INCLUDED
#define VG_ROUND_HISTORY_HPP_INCLUDED

#include <array>
#include <unordered_map>
#include <vector>

#include "linkage_model.hpp"
#include "staged_site.hpp"

namespace vg {

using namespace std;

/**
 * The states the re-genotyping rounds have reached, so that the rounds can tell whether they have
 * converged or are cycling.
 *
 * A state is every staged site's chosen genotype. The rounds can cycle: dropping and reinstating
 * a subtree is a discrete change, and the phase is a chain whose links move with the genotypes,
 * so no single quantity must increase.
 */
class RoundHistory {
public:
    /// Each staged site's chosen pair and ploidy, `{trav_first, trav_second, ploidy}`, keyed by
    /// record key. A site the model chose nothing for is left out.
    using State = unordered_map<size_t, std::array<int, 3>>;

    /// The state `model` has chosen for `sites`. Empty without a linkage model.
    static State state(StagedSiteTable& sites, const LinkageCollector* model);

    /// How many sites have a different chosen pair or ploidy in `after` than in `before`, counting
    /// a site that gained or lost a chosen answer as moved.
    static size_t changed(const State& before, const State& after);

    /// Remember `state` as the next round's.
    void remember(const State& state);

    /// The first remembered round that reached `state`, counting the first remembered state as
    /// round 1, or 0 when none did.
    size_t first_round_with(const State& state) const;

    /// How many rounds are remembered.
    size_t rounds() const { return digests.size(); }

private:
    /// A digest of a state, built by visiting its sites in the map's own order.
    static size_t digest(const State& state);

    /// Each remembered round's digest, in order.
    vector<size_t> digests;
};

}

#endif
