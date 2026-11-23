#ifndef VG_MOSAIC_WRITER_HPP_INCLUDED
#define VG_MOSAIC_WRITER_HPP_INCLUDED

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
    /// Write an `A` row wherever a strand's called allele is not what its rows' haplotypes spell,
    /// so that the strand rebuilt from the file takes the called allele at every genotyped site.
    /// The `H` rows are the same either way.
    bool exact = false;
};

/**
 * The allele each strand was called with at each genotyped site, as a walk and the changes the
 * allele makes to it. The exact mosaic's `A` rows are written from these.
 *
 * A site's alleles are freed once its record is rendered, and the mosaic is written after every
 * record, so the called alleles are copied here before the render.
 */
class CalledAlleles {
public:
    /// One strand's called allele at a site.
    struct Allele {
        /// The candidate traversal the allele is, as `PhaseCall::trav_first` numbers it.
        int traversal = -1;
        /// The allele's walk, from the site's start bound node to its end bound node, as
        /// oriented GBWT nodes.
        std::span<const gbwt::node_type> walk;
        /// The changes the allele makes to the sequences of its walk's steps, or null when the
        /// walk spells the allele.
        const vector<AlleleEdit>* edits = nullptr;
        /// True for the allele of an edit site that makes no change. It says that the strand
        /// lacks the site's change, and nothing about the walk the strand takes: its walk is
        /// only the route the change would lie on, and can pass through a snarl that the strand
        /// crosses another way.
        bool no_change = false;
    };

    /// Copy the alleles that `phase` puts on the strands of the site with `alleles`, which is an
    /// edit site if `edit_site`.
    void add(const HandleGraph& graph, const LinkageCollector::PhaseCall& phase,
             const CandidateAlleles& alleles, bool edit_site);

    /// Get the allele in slot `slot` of the site with `record_key`: 0 is the allele of
    /// `PhaseCall::trav_first` and 1 that of `trav_second`. Return false if none was copied.
    bool find(size_t record_key, int slot, Allele& out) const;

private:
    /// Where one site's two alleles are in `nodes` and `edit_sets`.
    struct Entry {
        int traversal[2] = {-1, -1};
        size_t offset[2] = {0, 0};
        uint32_t length[2] = {0, 0};
        /// The index in `edit_sets` of each allele's changes, or -1 for none.
        int32_t edits[2] = {-1, -1};
        /// Whether the site is an edit site.
        bool edit_site = false;
    };
    unordered_map<size_t, Entry> entries;
    /// The walks of every allele, one after another.
    vector<gbwt::node_type> nodes;
    vector<vector<AlleleEdit>> edit_sets;
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
    /// Of those, the rows that give a sequence, since the allele's walk does not spell it.
    std::atomic<size_t> allele_rows_with_sequence{0};
    /// Of those, the rows that apply to the walk of an earlier `A` row.
    std::atomic<size_t> allele_rows_nested{0};
    /// Of those, the rows of sites that have no VCF line.
    std::atomic<size_t> allele_rows_unwritten{0};
    /// Of those, the rows whose walk joins no row, since the site's row cannot be walked.
    std::atomic<size_t> allele_rows_standalone{0};
    /// Of those, the rows of sites the strand's walk crosses from end to start.
    std::atomic<size_t> allele_rows_reversed{0};
    /// Strands of genotyped sites whose called allele the strand's walk was compared with.
    std::atomic<size_t> strands_compared{0};
    /// Strands of genotyped sites whose bounds the strand's walk does not cross, so that no row
    /// could put the called allele on it.
    std::atomic<size_t> strands_not_crossed{0};
    /// Strands of genotyped sites whose called allele was not copied before the render.
    std::atomic<size_t> strands_without_allele{0};
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

    /// Whether the mosaic is to be exact, which needs `called_alleles()` filled before `write`.
    bool is_exact() const { return is_enabled() && params.exact; }

    /// The called alleles an exact mosaic's `A` rows are written from.
    CalledAlleles& called_alleles() { return alleles; }

    /// Write the mosaic of `phasing`, the phase calls of the sites that have a VCF line, in the
    /// order the linkage model produced them. The runs are walked through `panel`'s GBWT, and the
    /// header names `sample_name`. Reports the counters on stderr.
    ///
    /// An exact mosaic also needs `genotyped`, the phase calls of every genotyped site in the
    /// same order, with or without a line, and `graph`, for the sequences of changed alleles.
    void write(const vector<LinkageCollector::PhaseCall>& phasing, const PanelLookup& panel,
               const string& sample_name,
               const vector<LinkageCollector::PhaseCall>* genotyped = nullptr,
               const HandleGraph* graph = nullptr);

private:
    MosaicParams params;
    MosaicCounters counters;
    CalledAlleles alleles;
};

}

#endif
