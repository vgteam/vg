#ifndef VG_GIRAFFE_ENGINE_HPP_INCLUDED
#define VG_GIRAFFE_ENGINE_HPP_INCLUDED

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
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

/// The result of surjecting an alignment onto a haplotype path: where it
/// lands, or why it could not be placed. Fields other than status are set
/// only when status is OK.
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

    /// Find the haplotypes that follow each alignment's path.
    bool assign_haplotypes = true;
    /// When assigning haplotypes, also find the boundary nodes of the
    /// smallest snarl around each alignment's path.
    bool find_snarl_boundaries = true;
};

/// A read to map, and the path to surject its alignments onto.
struct GiraffeServerRequest {
    /// The read's name, sequence and base qualities. The qualities are raw
    /// Phred scores. If they are empty, base quality does not cap MAPQ.
    Alignment read;
    /// One of the surjection target paths, or empty to surject nothing.
    std::string surjection_target;
};

/// An alignment of a read, with the haplotypes that follow its path and its
/// surjection. A read can have several.
struct MappedAlignment {
    Alignment alignment;
    /// The runs of the alignment's path that haplotypes follow, in path
    /// order. Unset if haplotype assignment is off.
    std::optional<std::vector<HaplotypeSegment>> haplotype_segments;
    /// Where the haplotypes that follow the whole path leave the snarl around
    /// it. Both sides are unset if finding the boundaries is off.
    SnarlBoundaries snarl_boundaries;
    /// Where the alignment lands on the read's surjection target. Unset if
    /// the read names no target.
    std::optional<HaplotypeSurjectionResult> surjection;
};

/// The alignments of a read, best first. An unmapped read has one unmapped
/// alignment.
struct MappedRead {
    std::vector<MappedAlignment> alignments;
    /// Why the read could not be mapped. If set, there are no alignments.
    std::optional<std::string> error;
};

/// A run of steps of a path, from step_begin to step_end inclusive, with the
/// path offsets of both steps. surject_with_anchors() gets its positions on
/// the target path from these.
///
/// Each step is a GBWT node and an offset in that node's GBWT record, which
/// is how GBWTGraph encodes a step handle.
struct PathAnchor {
    uint64_t step_begin_node = 0;
    uint64_t step_begin_offset = 0;
    uint64_t step_end_node = 0;
    uint64_t step_end_offset = 0;
    size_t path_offset_step_begin = 0;
    size_t path_offset_step_end = 0;
};

/**
 * Holds Giraffe's indexes and a mapper, and maps and surjects reads with them.
 *
 * The engine maps batches of reads. It can find the haplotypes that carry each
 * alignment's path, and can surject the alignment onto one haplotype path. A
 * path indexed at load time can be a surjection target for map_reads(). Any
 * other path can be a target for surject_with_anchors(), given anchors that
 * place some of its steps on the path.
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

    /// The graph.
    const gbwtgraph::GBZ& gbz() const;

    /**
     * Map a batch of reads.
     *
     * Return one MappedRead per read, in input order. Each alignment has its
     * haplotypes and, if the read names a surjection target, its surjection.
     *
     * Each read is mapped in an OpenMP task. Call this from inside a parallel
     * region to map on all its threads. Outside of one, the reads are mapped
     * one at a time.
     *
     * If mapping a read throws, the read has no name, or its qualities and
     * sequence differ in length, the read gets the message as its error. The
     * other reads are still mapped.
     *
     * @throws std::runtime_error if the engine is not loaded.
     */
    std::vector<MappedRead> map_reads(const std::vector<GiraffeServerRequest>& reads);

    /**
     * Surject a graph alignment onto a haplotype path, using anchors on the
     * path for positions.
     *
     * The anchors stand in for a position index of the whole path, which the
     * engine has only for the paths indexed at load time.
     *
     * `target_haplotype` is a path name, such as GRCh38#0#chrM, or a sample
     * and haplotype, such as HG002#2. For a name the graph does not have, the
     * target is the path of the first anchor's step. The result's path_name
     * is the name of the path that was used.
     *
     * @pre All anchors lie on one path.
     *
     * @param target_path_length Length of the target path in bases, or 0 to
     * have it computed here by walking the path.
     *
     * @throws std::runtime_error if the engine is not loaded. An exception
     * from the surjector is not caught.
     */
    HaplotypeSurjectionResult surject_with_anchors(
        const Alignment& aln,
        const std::vector<PathAnchor>& anchors,
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
