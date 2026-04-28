#ifndef VG_GIRAFFE_ENGINE_HPP_INCLUDED
#define VG_GIRAFFE_ENGINE_HPP_INCLUDED

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <bdsg/overlays/reference_path_overlay.hpp>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include <handlegraph/path_position_handle_graph.hpp>

#include "haplotype_assignment.hpp"
#include "minimizer_mapper.hpp"
#include "snarl_distance_index.hpp"
#include "surjector.hpp"
#include "zip_code.hpp"

namespace vg {

// Forward declaration so we can mention Alignment in surjection types
// without forcing every includer to drag in vg.pb.h.
class Alignment;

struct GiraffeEnginePaths {
    std::string gbz_path;
    std::string minimizer_path;
    std::string distance_path;
    std::string zipcode_path;
};

/// Outcome of attempting to surject one alignment onto one named haplotype
/// path. The status field captures every observable failure mode; the rest
/// is only meaningful when status == OK.
struct HaplotypeSurjectionResult {
    enum class Status {
        OK,                  ///< Surjection succeeded against the target path.
        EMPTY_INPUT,         ///< Source alignment had no graph path.
        UNKNOWN_PATH,        ///< Target name is not a known path in the GBZ.
        PATH_NOT_INDEXED,    ///< Target name exists but was not pre-indexed for surjection.
        INCOMPATIBLE,        ///< Target haplotype is not in the assigner's candidate set.
        SURJECTION_FAILED,   ///< Surjector returned an unmapped/empty result.
    };

    Status status = Status::OK;

    /// Path name on success; otherwise empty.
    std::string path_name;
    /// 0-based path position of the surjected alignment; -1 if unset.
    int64_t path_position = -1;
    /// Strand of the alignment relative to the path's forward orientation.
    bool path_reverse = false;
    /// Score of the surjected linear alignment.
    int32_t score = 0;
    /// Mapping quality copied from the source alignment.
    int32_t mapping_quality = 0;
    /// SAM-style CIGAR for the surjected alignment.
    std::string cigar;
};

/// Convert a status enum value into a short stable token suitable for
/// emitting in GAF tags or framed protocol fields.
const char* haplotype_surjection_status_token(HaplotypeSurjectionResult::Status status);

/**
 * Owns a Surjector configured against a PathPositionHandleGraph that already
 * has the desired target paths indexed for positions. Used by GiraffeEngine
 * to project graph alignments onto a specific GBZ-known haplotype path.
 *
 * The position graph and the GBZ must outlive this object. Thread-safe after
 * construction: surject() is const and Surjector::surject is const.
 */
class HaplotypeSurjector {
public:
    HaplotypeSurjector(const gbwtgraph::GBZ& gbz,
                       const handlegraph::PathPositionHandleGraph& position_graph);

    /// True if the given haplotype path name is known to the GBZ AND has
    /// positional indexing in the position graph (i.e. is surjectable).
    bool path_is_indexed(const std::string& path_name) const;

    /// Attempt to surject `aln` onto the path named `target_path_name`.
    /// `candidate_seq_ids` is consulted for the compatibility check; if
    /// neither orientation of the target path's GBWT sequence id is
    /// present, status is INCOMPATIBLE and no surjection is attempted.
    HaplotypeSurjectionResult surject(const Alignment& aln,
                                      const std::string& target_path_name,
                                      const std::vector<gbwt::size_type>& candidate_seq_ids) const;

private:
    const gbwtgraph::GBZ* gbz_;
    const handlegraph::PathPositionHandleGraph* position_graph_;
    Surjector surjector_;
};

struct GiraffeEngineConfig {
    size_t threads = 1;
    size_t max_multimaps = 1;
    bool preload_distance_index = true;

    /// Names of haplotype paths in the GBZ to pre-index for surjection.
    /// If empty, surjection is disabled regardless of per-read requests.
    /// These are passed to bdsg::ReferencePathOverlay as extra_path_names so
    /// haplotype-sense paths become surjectable, not just reference paths.
    std::vector<std::string> surjection_target_paths;
};

struct GiraffeFastqRead {
    std::string name;
    std::string sequence;
    std::string quality;
    /// Optional: name of a pre-indexed haplotype path to surject this read
    /// onto. Empty means "no surjection". If set, the engine will surject
    /// every reported alignment for this read onto the given path and
    /// append the result as GAF tags.
    std::string surjection_target;
};

class GiraffeEngine {
public:
    GiraffeEngine() = default;

    void load(const GiraffeEnginePaths& paths, const GiraffeEngineConfig& config);

    bool is_loaded() const;

    const GiraffeEngineConfig& config() const;

    /// If true, map_reads() post-processes each final Alignment to annotate
    /// it with the set of haplotypes that traverse its path.
    bool assign_haplotypes = true;
    /// If true and assign_haplotypes is true, the haplotype search is
    /// extended outward to the nearest enclosing snarl boundaries.
    bool assign_haplotypes_extend_to_snarls = true;

    std::vector<std::string> gaf_header_lines() const;

    // Returns one vector per input read, preserving order.
    // Each inner vector contains one or more GAF lines (multimaps).
    std::vector<std::vector<std::string>> map_reads(const std::vector<GiraffeFastqRead>& reads);

private:
    void require_loaded() const;

    GiraffeEngineConfig config_{};
    std::unique_ptr<gbwtgraph::GBZ> gbz_;
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index_;
    std::unique_ptr<SnarlDistanceIndex> distance_index_;
    std::unique_ptr<ZipCodeCollection> zipcodes_;
    std::unique_ptr<MinimizerMapper> mapper_;
    std::unique_ptr<HaplotypeAssigner> haplotype_assigner_;

    /// Path-position overlay built once at load time over gbz_->graph,
    /// indexing the configured surjection_target_paths so the surjector
    /// can answer position queries on those haplotypes.
    std::unique_ptr<bdsg::ReferencePathOverlay> position_overlay_;
    std::unique_ptr<HaplotypeSurjector> haplotype_surjector_;
};

} // namespace vg

#endif // VG_GIRAFFE_ENGINE_HPP_INCLUDED
