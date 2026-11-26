#ifndef VG_MOSAIC_WRITER_HPP_INCLUDED
#define VG_MOSAIC_WRITER_HPP_INCLUDED

#include <array>
#include <atomic>
#include <cstddef>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include "linkage_model.hpp"
#include "materialized_site.hpp"
#include "panel_lookup.hpp"

namespace vg::multipass {

using namespace std;

/// Where and how to write the mosaic.
struct MosaicParams {
    /// Where to write the mosaic; empty for none.
    string path;
    /// The graph the mosaic is to be read against, for the header.
    string graph_name;
    /// Panel index -> "sample#phase", the unit the linkage model works in: a haplotype stored as
    /// several GBWT paths is one haplotype. With a row's contig, that is enough to find its paths.
    /// The index means nothing outside this run, so the header writes the whole mapping.
    vector<string> haplotype_names;
    /// The full names of the reference paths called against, such as `CHM13#0#chr20`. The rows
    /// give only the contig as the VCF names it, and a graph can hold several references, so the
    /// header lists them.
    vector<string> reference_paths;
    /// Fill a gap across which no panel haplotype can be followed with the reference, so that a
    /// strand stays one walk. The fill is marked `ref` in the file.
    bool patch_gaps = true;
    /// Include nested sites in the runs, so that a switch of haplotype at a nested site starts a
    /// new row. Off leaves nested sites out, so a strand follows its enclosing site's haplotype
    /// through them, and its walk need not spell the nested sites' called alleles. Recorded in the
    /// file's #nested header.
    bool keep_nested = true;
    /// Carry the flanking haplotype through a stretch the panel cannot explain, rather than
    /// writing an unwalkable row and breaking the path. It gives up the called alleles across
    /// those sites for a contiguous path.
    bool connect_unexplained = true;
    /// Write an `A` row wherever a strand, as the records spell it, is not what its rows'
    /// haplotypes spell, so that the strand rebuilt from the file is the strand the records
    /// describe. The `H` rows are the same either way.
    bool exact = false;
};

/**
 * What each strand carries at each top-level site, as the records spell it. The exact mosaic's
 * `A` rows are written from these.
 *
 * A site's alleles are freed once its record is rendered, and the mosaic is written after every
 * record, so the strands are copied here before the render.
 */
class StrandAlleles {
public:
    /// The entry of one step of a walk in an `A` row's bases, where it is not `*`: the step's
    /// bases are not those of its node, or no edge of the graph leads to the next step.
    struct StepBases {
        /// The index of the step in the walk.
        uint32_t step;
        /// The entry of the step in an `A` row's bases.
        string entry;
    };

    /// One strand at a top-level snarl.
    struct Allele {
        /// The strand's walk, from the site's start bound node to its end bound node, as
        /// oriented GBWT nodes.
        std::span<const gbwt::node_type> walk;
        /// The steps of the walk that have an entry other than `*`, in walk order.
        std::span<const StepBases> changed;
    };

    /// A change that one strand makes to the bases of a node at a top-level edit site, along the
    /// forward strand of the node: `removed` bases from offset `at` on become `added`.
    struct NodeChange {
        nid_t node;
        size_t at;
        size_t removed;
        string added;
    };

    /**
     * Keep the walk of strand `strand` at the top-level snarl with `record_key`, and the bases
     * of its steps.
     *
     * @param steps The bases of each step of `walk`, or empty where each has those of its node.
     * @return False, keeping nothing, if the strand changes the bases of a bound node of the
     * site, which no row could say. Bases put after the start bound node are not such a change.
     */
    bool add_site(const HandleGraph& graph, size_t record_key, int strand, const Traversal& walk,
                  const StepSequences& steps);

    /// Keep the changes `edits` that strand `strand` makes to the nodes of `walk` at the
    /// top-level edit site with `record_key`.
    void add_edit_site(const HandleGraph& graph, size_t record_key, int strand,
                       const Traversal& walk, const vector<AlleleEdit>& edits);

    /// Get what strand `strand` carries at the top-level snarl with `record_key`. Return false
    /// if nothing was kept.
    bool find(size_t record_key, int strand, Allele& out) const;

    /// Get the changes strand `strand` makes at the top-level edit site with `record_key`.
    std::span<const NodeChange> changes_of(size_t record_key, int strand) const;

    /// The top-level sites that have sites under them and no phase, so that no row says what a
    /// strand carries at them or under them.
    size_t unphased_parents = 0;
    /// The strands that `add_site` refused.
    size_t bound_changed = 0;

private:
    /// A stretch of one of the vectors below.
    struct Slice {
        size_t at = 0;
        size_t size = 0;
    };
    /// Where the walk and the changed steps of each strand of a snarl are.
    struct SiteEntry {
        Slice walk[2];
        Slice changed[2];
    };
    unordered_map<size_t, SiteEntry> sites;
    /// Where the changes of each strand of an edit site are in `node_changes`.
    unordered_map<size_t, array<Slice, 2>> edit_sites;
    /// The walks of every strand, one after another.
    vector<gbwt::node_type> nodes;
    vector<StepBases> step_bases;
    vector<NodeChange> node_changes;
};

/// Counters for the mosaic writer, reported when it writes the mosaic.
struct MosaicCounters {
    /// Runs with no position to walk from. Panel haplotypes are often fragments, so this is
    /// reported rather than expected to be zero.
    std::atomic<size_t> unwalkable{0};
    /// Of those, the ones that are only a head: the run could be walked from a later site, so the
    /// walkable rest is written separately.
    std::atomic<size_t> head_clipped{0};
    /// Run boundaries closed by following the first run's haplotype to the next run.
    std::atomic<size_t> extended{0};
    /// Run boundaries that no haplotype could cross, left as a gap that starts a new fragment.
    std::atomic<size_t> gap_left{0};
    /// Run boundaries that no panel haplotype could cross, filled with the reference.
    std::atomic<size_t> patched{0};
    /// Rows whose own haplotype does not span them, rewritten as a reference substitution.
    std::atomic<size_t> row_to_ref{0};
    /// Run boundaries where the walk enters a child snarl, stated at the child's start node.
    std::atomic<size_t> nested_enter{0};
    /// Run boundaries where the walk leaves a child snarl, stated at the child's end node.
    std::atomic<size_t> nested_leave{0};
    /// Rows the carried direction could not walk but the other direction could, as at an inversion.
    std::atomic<size_t> direction_broken{0};
    /// Run boundaries closed by following the second run's haplotype back to the first run.
    std::atomic<size_t> extended_left{0};

    // The rest count the `A` rows of an exact mosaic, and are zero otherwise.

    /// `A` rows written.
    std::atomic<size_t> allele_rows{0};
    /// Of those, the rows of one node, which give the changes of edit sites to it.
    std::atomic<size_t> allele_rows_node{0};
    /// Of those, the rows that give bases, since a step of the walk does not spell its node, or
    /// is not followed by an edge of the graph.
    std::atomic<size_t> allele_rows_with_bases{0};
    /// Of those, the rows of sites that have no VCF line.
    std::atomic<size_t> allele_rows_unwritten{0};
    /// Strands of top-level snarls that were compared with the strand's walk.
    std::atomic<size_t> strands_compared{0};
    /// Strands of top-level snarls whose bounds the strand's walk does not cross, so that no
    /// row could put the strand's allele on it.
    std::atomic<size_t> strands_not_crossed{0};
    /// Strands of top-level snarls that the strand's walk crosses from end to start, which get
    /// no row.
    std::atomic<size_t> strands_backward{0};
    /// Nodes that a strand's edit sites change where no row could say so: the strand's walk
    /// does not visit the node, or visits it inside a top-level snarl.
    std::atomic<size_t> node_changes_unwritten{0};
};

/**
 * A writer of the mosaic, a file that describes each strand of the sample as panel haplotypes.
 *
 * Each row of the mosaic is a run of sites over which one strand copies one panel haplotype. A row
 * is located by node IDs and a GBWT position, so that a consumer can rebuild the strand by walking
 * the haplotype through the GBWT.
 */
class MosaicWriter {
public:
    /// Where and how to write the mosaic. No mosaic is written until a path is given.
    void set_params(MosaicParams new_params) { params = std::move(new_params); }

    /// Whether a mosaic is to be written.
    bool is_enabled() const { return !params.path.empty(); }

    /// Whether the mosaic is to be exact, which needs `strand_alleles()` filled before `write`.
    bool is_exact() const { return is_enabled() && params.exact; }

    /// What each strand carries at each top-level site, which an exact mosaic's `A` rows are
    /// written from.
    StrandAlleles& strand_alleles() { return alleles; }

    /// Write the mosaic of `phasing`, the phase calls of the sites that have a VCF line, in the
    /// order the linkage model produced them. The runs are walked through `panel`'s GBWT, and the
    /// header names `sample_name`. Reports the counters on stderr.
    ///
    /// An exact mosaic also needs `genotyped`, the phase calls of every genotyped site in the
    /// same order, with or without a line, and `graph`, for the sequences of nodes.
    void write(const vector<LinkageCollector::PhaseCall>& phasing, const PanelLookup& panel,
               const string& sample_name,
               const vector<LinkageCollector::PhaseCall>* genotyped = nullptr,
               const HandleGraph* graph = nullptr);

private:
    MosaicParams params;
    MosaicCounters counters;
    StrandAlleles alleles;
};

}

#endif
