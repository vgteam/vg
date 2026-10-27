#ifndef VG_MOSAIC_WRITER_HPP_INCLUDED
#define VG_MOSAIC_WRITER_HPP_INCLUDED

#include <atomic>
#include <cstddef>
#include <string>
#include <vector>

#include "linkage_model.hpp"
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

    /// Write the mosaic of `phasing`, the phase calls of the sites that have a VCF line, in the
    /// order the linkage model produced them. The runs are walked through `panel`'s GBWT, and the
    /// header names `sample_name`. Reports the counters on stderr.
    void write(const vector<LinkageCollector::PhaseCall>& phasing, const PanelLookup& panel,
               const string& sample_name);

private:
    MosaicParams params;
    MosaicCounters counters;
};

}

#endif
