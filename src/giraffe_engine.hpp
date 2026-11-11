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
#include <vg/vg.pb.h>

#include "haplotype_assignment.hpp"
#include "minimizer_mapper.hpp"
#include "snarl_distance_index.hpp"
#include "surjector.hpp"
#include "zip_code.hpp"

namespace vg {

/// Where one alignment lands on one haplotype path, or why it could not be
/// placed. Fields other than status are set only when status is OK.
struct HaplotypeSurjectionResult {
    enum class Status {
        OK,                  ///< The alignment was surjected onto the target path.
        EMPTY_INPUT,         ///< The alignment has no graph path.
        UNKNOWN_PATH,        ///< The graph has no path with the target's name.
        PATH_NOT_INDEXED,    ///< The graph has the target path, but its positions are not indexed.
        INCOMPATIBLE,        ///< The target path does not carry the alignment's path.
        SURJECTION_FAILED,   ///< The surjector found no alignment on the target path.
    };

    Status status = Status::OK;

    /// Name of the target path.
    std::string path_name;
    /// 0-based start of the alignment on the path's forward strand.
    int64_t path_position = -1;
    /// True if the alignment is on the path's reverse strand.
    bool path_reverse = false;
    /// Score of the surjected alignment.
    int32_t score = 0;
    /// MAPQ of the surjected alignment, which the surjector copies from the
    /// graph alignment.
    int32_t mapping_quality = 0;
    /// CIGAR of the surjected alignment.
    std::string cigar;

    /// Bases in each kind of edit of the surjected alignment. The CIGAR writes
    /// matches and mismatches both as M, so these are counted from the edits.
    uint64_t matches = 0;
    uint64_t mismatches = 0;
    uint64_t inserted = 0;    ///< Read bases missing from the path.
    uint64_t deleted = 0;     ///< Path bases missing from the read.
};

/// Return the token that stands for a status in GAF tags.
const char* haplotype_surjection_status_token(HaplotypeSurjectionResult::Status status);

/// What GiraffeEngine::load() loads, and how the engine maps.
struct GiraffeEngineConfig {
    /// The GBZ graph and Giraffe's minimizer, distance and zipcode indexes.
    std::string gbz_path;
    std::string minimizer_path;
    std::string distance_path;
    std::string zipcode_path;

    /// Number of mapping threads.
    size_t threads = 1;

    /// Page the whole distance index into memory at load time.
    bool preload_distance_index = true;

    /// Names of the paths that map_reads() can surject onto. If empty,
    /// map_reads() surjects nothing.
    std::vector<std::string> surjection_target_paths;

    /// Have the surjectors drop low-complexity anchors. vg giraffe turns this
    /// on with --prune-low-cplx.
    bool prune_suspicious_anchors = false;

    /// Tag each alignment with the haplotypes that carry its path.
    bool assign_haplotypes = true;
    /// When assigning haplotypes, extend the search out to the enclosing
    /// snarl's boundaries, keeping only haplotypes that also agree there.
    bool assign_haplotypes_extend_to_snarls = true;
};

/// One read to map, and the path to surject its alignments onto.
struct GiraffeServerRequest {
    /// The read's name, sequence and base qualities. The qualities are raw
    /// Phred scores. If they are empty, base quality does not cap MAPQ.
    Alignment read;
    /// One of the surjection target paths, or empty to surject nothing.
    std::string surjection_target;
};

/// One precomputed anchor, as giraffe-server reads it: a run of steps of the
/// target path, from step_begin to step_end inclusive, with the path offsets
/// of both steps.
///
/// Each step is a GBWT node and an offset in that node's GBWT record, which
/// is how GBWTGraph encodes a step handle. The engine does not use the read
/// and source mapping ranges.
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

/**
 * Giraffe with its indexes loaded once, for a mapping server.
 *
 * The engine maps batches of reads to GAF lines. It tags each alignment with
 * the haplotypes that carry its path, and can surject it onto one haplotype
 * path. A path indexed at load time can be a surjection target for
 * map_reads(). Any other path can be a target for surject_with_anchors(),
 * given anchors that place some of its steps on the path.
 */
class GiraffeEngine {
public:
    GiraffeEngine() = default;

    /**
     * Load the indexes and build the mapper.
     *
     * This also sets the OpenMP thread count for the whole process.
     *
     * @throws std::runtime_error if an index file is missing, the indexes do
     * not match, or a surjection target path is not in the graph.
     */
    void load(const GiraffeEngineConfig& config);

    /// Return true once load() has succeeded.
    bool is_loaded() const;

    const GiraffeEngineConfig& config() const;

    /// The mapper. The caller sets its parameters after load() and before
    /// mapping. The engine cannot do it: vg giraffe's options are in
    /// src/subcommand, which is linked into the vg binary and not into libvg.
    MinimizerMapper& mapper();

    /// The minimizer index the mapper uses.
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index() const;

    /// Return the GAF header lines that name the graph.
    std::vector<std::string> gaf_header_lines() const;

    /**
     * Map a batch of reads in parallel.
     *
     * Return one vector per read, in input order, with one GAF line per
     * alignment. An unmapped read gets one unmapped line. Each line carries
     * the haplotype tags and, if the read names a surjection target, the
     * surjection tags.
     *
     * @throws std::runtime_error if a read's qualities and sequence differ in
     * length. If mapping a read throws, the first such exception is rethrown
     * after the batch.
     */
    std::vector<std::vector<std::string>> map_reads(const std::vector<GiraffeServerRequest>& reads);

    /**
     * Surject a graph alignment onto a haplotype path, using anchors on the
     * path for positions.
     *
     * The anchors stand in for a position index of the whole path, which the
     * engine has only for the paths indexed at load time.
     *
     * `target_haplotype` is a path name, such as GRCh38#0#chrM, or a sample
     * and haplotype, such as HG002#2. For a name the graph does not have, the
     * target is the path of the first anchor's step.
     *
     * Return one GAF line: `gaf_line` with the surjection tags appended. On
     * success it also has an:i, the number of anchors, and ap:Z, the target
     * path's name. If the surjector throws, se:Z holds the error message.
     *
     * @pre All anchors lie on one path.
     *
     * @param target_path_length Length of the target path in bases, or 0 to
     * have it computed here by walking the path.
     *
     * @throws std::runtime_error if the engine is not loaded or `gaf_line` is
     * empty.
     */
    std::vector<std::string> surject_with_anchors(
        const std::string& read_name,
        const std::string& gaf_line,
        const std::vector<WireAnchor>& anchors,
        const std::string& target_haplotype,
        size_t target_path_length);

private:
    /// Throw if load() has not succeeded.
    void require_loaded() const;

    /// Surject an alignment onto a path whose positions were indexed at load
    /// time. The path must visit the nodes of the alignment's mappings
    /// [begin_mapping, end_mapping) in order, in either orientation, or the
    /// status is INCOMPATIBLE.
    HaplotypeSurjectionResult surject_to_indexed_path(
        const Alignment& aln, const std::string& target_path_name,
        size_t begin_mapping, size_t end_mapping) const;

    GiraffeEngineConfig config_{};
    std::unique_ptr<gbwtgraph::GBZ> gbz_;
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index_;
    std::unique_ptr<SnarlDistanceIndex> distance_index_;
    std::unique_ptr<ZipCodeCollection> zipcodes_;
    std::unique_ptr<MinimizerMapper> mapper_;
    std::unique_ptr<HaplotypeAssigner> haplotype_assigner_;

    /// Position index of the surjection target paths, and a surjector that
    /// uses it. Both are null if there are no target paths.
    std::unique_ptr<bdsg::ReferencePathOverlay> position_overlay_;
    std::unique_ptr<Surjector> overlay_surjector_;
};

} // namespace vg

#endif // VG_GIRAFFE_ENGINE_HPP_INCLUDED
