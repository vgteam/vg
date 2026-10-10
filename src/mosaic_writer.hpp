#ifndef VG_MOSAIC_WRITER_HPP_INCLUDED
#define VG_MOSAIC_WRITER_HPP_INCLUDED

#include <atomic>
#include <cstddef>
#include <string>
#include <vector>

#include "linkage_model.hpp"
#include "panel_lookup.hpp"

namespace vg {
namespace multipass {

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
    /// Run boundaries across which the first run's haplotype could be followed, and those across
    /// which it could not, which a reference fill or a new fragment has to cover.
    std::atomic<size_t> extended{0}, gap_left{0}, patched{0};
    /// Rows whose own haplotype does not span them, rewritten as a reference substitution.
    std::atomic<size_t> row_to_ref{0};
    /// Run boundaries between a parent and a child snarl, at the child's own boundary nodes.
    std::atomic<size_t> nested_enter{0}, nested_leave{0};
    /// Rows the current direction could not walk but the other direction could, at an inversion.
    std::atomic<size_t> direction_broken{0}, extended_left{0};
};

/**
 * Writes the mosaic: the per-site phasing collapsed into runs, stretches of sites over which a
 * strand copies one panel haplotype, each located by node ID so that a consumer can rebuild the
 * strand by walking the haplotype through the GBWT.
 */
class MosaicWriter {
public:
    /// Where and how to write the mosaic. No mosaic is written until a path is given.
    void set_params(MosaicParams params) { this->params = std::move(params); }

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
}

#endif
