#ifndef VG_GIRAFFE_ENGINE_HPP_INCLUDED
#define VG_GIRAFFE_ENGINE_HPP_INCLUDED

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <bdsg/overlays/reference_path_overlay.hpp>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include <vg/vg.pb.h>

#include "haplotype_assignment.hpp"
#include "minimizer_mapper.hpp"
#include "snarl_distance_index.hpp"
#include "surjector.hpp"
#include "zip_code.hpp"

namespace vg {

/// Where one alignment lands on one haplotype's path, or why it could not be
/// placed. Fields other than status are valid only when status == OK.
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

    /// Exact edit breakdown of the surjected alignment, counted from the
    /// alignment's edits rather than the CIGAR: a SAM CIGAR writes both matches
    /// and substitutions as 'M', so identity cannot be recovered from `cigar`
    /// alone. `matches` + `mismatches` is the aligned (M) length, which is what
    /// a PSL needs for its matches / misMatches columns.
    uint64_t matches = 0;
    uint64_t mismatches = 0;
    uint64_t inserted = 0;    ///< bases in the read absent from the target
    uint64_t deleted = 0;     ///< bases in the target absent from the read
};

/// Convert a status enum value into a short stable token suitable for
/// emitting in GAF tags or framed protocol fields.
const char* haplotype_surjection_status_token(HaplotypeSurjectionResult::Status status);

/// What GiraffeEngine::load() loads and how the engine then runs.
struct GiraffeEngineConfig {
    /// The graph plus Giraffe's minimizer, distance and zipcode indexes.
    std::string gbz_path;
    std::string minimizer_path;
    std::string distance_path;
    std::string zipcode_path;

    size_t threads = 1;

    /// Sets the mapper's parameters once it is built, given the minimizer index
    /// it uses. giraffe-server fills it from vg giraffe's presets and rules.
    /// Unset: MinimizerMapper's defaults.
    std::function<void(MinimizerMapper&, const gbwtgraph::DefaultMinimizerIndex&)> configure_mapper;

    /// If > 0, replaces max_multimaps after configure_mapper. Reporting is
    /// still capped by max_alignments.
    size_t max_multimaps = 0;

    bool preload_distance_index = true;

    /// Names of haplotype paths in the GBZ to pre-index for surjection.
    /// If empty, surjection is disabled regardless of per-read requests.
    /// These are passed to bdsg::ReferencePathOverlay as extra_path_names so
    /// haplotype-sense paths become surjectable, not just reference paths.
    std::vector<std::string> surjection_target_paths;

    /// Tag each alignment with the haplotypes that carry its path.
    bool assign_haplotypes = true;
    /// When assigning haplotypes, extend the search out to the enclosing snarl's
    /// boundaries, keeping only haplotypes that also agree there.
    bool assign_haplotypes_extend_to_snarls = true;
};

/// One read to map, as giraffe-server receives it.
struct GiraffeServerRequest {
    /// Name, sequence and base qualities (raw Phred, as vg stores them; empty
    /// if unknown, which like FASTA input to vg giraffe leaves MAPQ uncapped).
    Alignment read;
    /// Pre-indexed haplotype path to surject every alignment of the read onto,
    /// reported as GAF tags. Empty for none.
    std::string surjection_target;
};

/// Wire-format payload for one pre-computed anchor sent over the
/// giraffe-server stdin protocol. Mirrors panindexer::PrecomputedAnchor
/// fields that we serialize; step_handle_t is reconstructed locally from
/// the (node, offset_in_record) pair using gbwtgraph's documented
/// step encoding (as_integers(step)[0]=node, [1]=offset).
struct WireAnchor {
    uint64_t step_begin_node = 0;
    uint64_t step_begin_offset = 0;
    uint64_t step_end_node = 0;
    uint64_t step_end_offset = 0;
    size_t path_offset_step_begin = 0;
    size_t path_offset_step_end = 0;
    size_t read_begin_offset = 0;
    size_t read_end_offset = 0;
    size_t source_mapping_begin = 0;
    size_t source_mapping_end = 0;
};

/// Giraffe with its indexes loaded once: maps reads to GAF, and surjects
/// alignments onto haplotypes. Backs vg giraffe-server.
class GiraffeEngine {
public:
    GiraffeEngine() = default;

    void load(const GiraffeEngineConfig& config);

    bool is_loaded() const;

    const GiraffeEngineConfig& config() const;

    std::vector<std::string> gaf_header_lines() const;

    // Returns one vector per input read, preserving order.
    // Each inner vector contains one or more GAF lines (multimaps).
    std::vector<std::vector<std::string>> map_reads(const std::vector<GiraffeServerRequest>& reads);

    /// Surject a graph alignment (`gaf_line`) onto `target_haplotype` using
    /// pre-computed anchors instead of a built-in ReferencePathOverlay.
    /// Returns one GAF line containing the original alignment with appended
    /// surjection tags (sj:, sn:, sp:, sr:, ss:, sm:, sc:), mirroring the
    /// per-read surjection-tag schema emitted by map_reads.
    ///
    /// `target_path_length` is the caller-precomputed total length in bases
    /// of the target path; it's what AnchorBackedPositionGraph reports for
    /// `get_path_length`. Pass 0 to fall back to walking the path here.
    ///
    /// Throws if the engine wasn't loaded, if the GAF can't be parsed, or
    /// if `target_haplotype` is not a known path in the GBZ.
    std::vector<std::string> surject_with_anchors(
        const std::string& read_name,
        const std::string& gaf_line,
        const std::vector<WireAnchor>& anchors,
        const std::string& target_haplotype,
        size_t target_path_length);

private:
    void require_loaded() const;

    /// Surject aln onto the pre-indexed path target_path_name, if one of
    /// candidate_seq_ids (from the haplotype assigner) is that path.
    HaplotypeSurjectionResult surject_to_indexed_path(
        const Alignment& aln, const std::string& target_path_name,
        const std::vector<gbwt::size_type>& candidate_seq_ids) const;

    GiraffeEngineConfig config_{};
    std::unique_ptr<gbwtgraph::GBZ> gbz_;
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index_;
    std::unique_ptr<SnarlDistanceIndex> distance_index_;
    std::unique_ptr<ZipCodeCollection> zipcodes_;
    std::unique_ptr<MinimizerMapper> mapper_;
    std::unique_ptr<HaplotypeAssigner> haplotype_assigner_;

    /// Position index over surjection_target_paths, and a Surjector on it.
    /// Null when no targets are configured.
    std::unique_ptr<bdsg::ReferencePathOverlay> position_overlay_;
    std::unique_ptr<Surjector> overlay_surjector_;
};

} // namespace vg

#endif // VG_GIRAFFE_ENGINE_HPP_INCLUDED
