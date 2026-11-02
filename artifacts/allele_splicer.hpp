#ifndef VG_ALLELE_SPLICER_HPP_INCLUDED
#define VG_ALLELE_SPLICER_HPP_INCLUDED

#include <unordered_map>
#include <vector>

#include "../src/handle.hpp"
#include "phase_table.hpp"
#include "site_values.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A splicer that spells the alleles of a site with the chosen alleles of the sites staged under
 * it.
 *
 * A record spells each called ALT on the first strand of the site's phase that carries it. On
 * that strand, it spells every site staged under the site by its own chosen allele, recursively.
 * The splicer reads the chosen alleles from the phase table as the render has it.
 */
class AlleleSplicer {
public:
    using SplicedAlleles = multipass::SplicedAlleles;

    /// One change to the sequence of one step of a spliced walk.
    struct StepEdit {
        /// The index of the step.
        size_t step;
        /// The offset in the step's sequence where the change starts.
        size_t at;
        /// The number of bases the change removes.
        size_t removed;
        /// The bases the change adds.
        string added;
        /// The record key of the child site that makes the change.
        size_t key;
    };

    /// Make a splicer over `graph` that reads the chosen alleles from `phase_table`.
    ///
    /// Neither is owned.
    AlleleSplicer(const HandleGraph& graph, const EditSnarlDecompositionOverlay& sites,
                  const PhaseTable& phase_table);

    /// Index the staged sites under each site of `staged_sites`.
    ///
    /// Index nothing unless `phasing` is on and the phase table has its render phases, since
    /// without a phase there is no strand to splice on. The sites must not move until `clear()`.
    void index(StagedSiteTable& staged_sites, bool phasing);

    /// Forget the index.
    void clear();

    /// Return true if any staged site is indexed under the site with record key `record_key`.
    bool has_children(size_t record_key) const;

    /// Fill `out` with the alleles the record of `site` spells.
    ///
    /// @return False if no allele changed.
    bool splice(const StagedSite& site, SplicedAlleles& out) const;

private:
    /// The chosen allele of one child site on one strand, as the steps that replace the parent's
    /// visits [begin, end).
    struct Replacement {
        int begin;
        int end;
        Traversal walk;
        StepSequences steps;
        /// Whether the child lies inside nodes, so that its allele changes the bases of the
        /// parent's steps rather than replacing them.
        bool inside_node;
        /// The record key of the child.
        size_t key;
    };

    /**
     * Spell allele `trav` of `site` on strand `strand` of the site's phase, with each child staged
     * under it spelled by its chosen allele on that strand, recursively.
     *
     * @param haploid Whether the site has one strand, which all its one-copy children are on.
     * @param walk Set to the spliced walk.
     * @param steps Set to the sequence of each step of `walk`.
     * @return False, leaving `walk` and `steps` empty, if no child changes the allele.
     */
    bool spliced_allele(const StagedSite& site, int trav, int strand, bool haploid,
                        Traversal& walk, StepSequences& steps) const;

    /**
     * Make the replacement of a child site's chosen allele on one strand of its parent's walk.
     *
     * @param first_visit The index in `parent_walk` of the first visit to each node.
     * @return False if the child has no chosen allele on the strand, or the parent's walk does
     * not cross it.
     */
    bool chosen_replacement(const StagedSite& child, int strand, bool haploid,
                            const Traversal& parent_walk,
                            const unordered_map<nid_t, int>& first_visit, Replacement& rep) const;

    /// Spell `parent_walk` into `walk` and `steps`, with each of `routes` in place of the visits
    /// it covers.
    ///
    /// A route that overlaps one before it is skipped.
    ///
    /// @param own The step sequences of the parent's allele, or null if its walk spells it.
    /// @param routes Sorted by where they begin.
    /// @return True if any route changed the allele.
    bool splice_routes(const Traversal& parent_walk, const StepSequences* own,
                       const vector<Replacement>& routes, Traversal& walk,
                       StepSequences& steps) const;

    /// Find the changes that children inside nodes make to the steps of the spliced walk `walk`.
    ///
    /// Each child changes the bases its allele differs from its nodes in, where `walk` first
    /// visits its first node and then the rest of its walk, read the way `walk` reads them.
    vector<StepEdit> step_edits_of(const vector<Replacement>& edits, const Traversal& walk) const;

    /// Make the edits to the sequences of the steps, right to left within each step.
    ///
    /// Going right to left, an edit that changes the length does not move the offsets of those
    /// before it. An edit that overlaps one already made to its right is skipped.
    ///
    /// @return True if any sequence changed.
    static bool apply_step_edits(vector<StepEdit> step_edits, StepSequences& steps);

    const HandleGraph& graph;
    /// The overlay of sites that the staged sites' net handles are in.
    const EditSnarlDecompositionOverlay& sites;
    /// Gives the chosen allele of each site on each strand.
    const PhaseTable& phase_table;
    /// The staged sites under each site, by the record key of the parent.
    unordered_map<size_t, vector<const StagedSite*>> children;
};

}

#endif
