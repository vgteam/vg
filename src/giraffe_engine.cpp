#include "giraffe_engine.hpp"

#include <algorithm>
#include <fstream>
#include <atomic>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <omp.h>

#include <vg/io/alignment_io.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/vg.pb.h>

#include <limits>

#include "alignment.hpp"
#include "annotation.hpp"
#include "anchor_backed_position_graph.hpp"
#include "gbwtgraph_helper.hpp"
#include "surjector.hpp"
#include "utility.hpp"

#include <vg/io/gafkluge.hpp>
#include <handlegraph/util.hpp>

namespace vg {

using namespace std;

const char* haplotype_surjection_status_token(HaplotypeSurjectionResult::Status status) {
    using S = HaplotypeSurjectionResult::Status;
    switch (status) {
        case S::OK:                return "ok";
        case S::EMPTY_INPUT:       return "empty";
        case S::UNKNOWN_PATH:      return "unknown_path";
        case S::PATH_NOT_INDEXED:  return "not_indexed";
        case S::INCOMPATIBLE:      return "incompatible";
        case S::SURJECTION_FAILED: return "failed";
    }
    return "unknown";
}

HaplotypeSurjector::HaplotypeSurjector(const gbwtgraph::GBZ& gbz,
                                       const handlegraph::PathPositionHandleGraph& position_graph)
    : gbz_(&gbz),
      position_graph_(&position_graph),
      surjector_(&position_graph)
{
    // Surjector keeps default vg scoring; the engine can tune it later if
    // we expose hooks. AlignerClient's default ctor builds default aligners.
}

bool HaplotypeSurjector::path_is_indexed(const string& path_name) const {
    // The position graph (e.g. ReferencePathOverlay) reports has_path() only
    // for paths it has actually indexed. That is the surjectable set.
    return position_graph_->has_path(path_name);
}

namespace {

/// Apply the "blat" preset to a MinimizerMapper.
///
/// Mirrors the "blat" entry in giraffe_main.cpp's `presets` map (defined
/// next to "hifi"). Inherits all of "hifi"'s long-read code-path values and
/// then applies the BLAT overrides flagged below. Keep this in sync with
/// the preset definition in giraffe_main.cpp — the engine does not use the
/// CLI parser, so the values are duplicated here for direct field
/// assignment.
void apply_blat_preset(MinimizerMapper& m) {
    // ---- inherited from "hifi" (long-read code path) ----
    m.align_from_chains = true;
    m.use_explored_cap = false;
    m.max_unique_min = 79;
    m.num_bp_per_min = 152;
    m.minimizer_downsampling_window_count = 15;
    m.minimizer_downsampling_max_window_length = 227;
    m.hit_cap = 0;
    m.minimizer_score_fraction = 1.0;
    m.hard_hit_cap = 13614;
    m.gapless_extension_limit = 0;
    m.zipcode_tree_score_threshold = 100.0;
    m.pad_zipcode_tree_score_threshold = 50.0;
    m.zipcode_tree_coverage_threshold = 0.5;
    m.zipcode_tree_scale = 2.0;
    m.min_chaining_problems = 6;
    m.max_chaining_problems = std::numeric_limits<int>::max();
    m.max_graph_lookback_bases = 20000;
    m.max_graph_lookback_bases_per_base = 0.10501002120802233;
    m.max_indel_bases = 5000;
    m.max_indel_bases_per_base = 2.45;
    m.item_bonus = 20;
    m.item_scale = 1.0;
    m.gap_scale = 0.2;
    m.max_min_chain_score = 100;
    m.max_skipped_bases = 1000;
    m.max_chain_connection = 233;
    m.max_tail_length = 68;
    m.max_tail_gap = 150;
    m.max_middle_gap = 500;
    m.max_dp_cells = 8000000000ULL;
    m.wfa_distance = 33;
    m.wfa_distance_per_base = 0.195722;
    m.wfa_max_distance = 240;
    m.wfa_max_mismatches = 2;
    m.wfa_max_mismatches_per_base = 0.05;
    m.wfa_max_max_mismatches = 15;
    // watchdog-timeout / batch-size / prune-low-cplx live on
    // GiraffeMainOptions in vg, not on MinimizerMapper; they're
    // program-level concerns we don't need in the engine.

    // ---- BLAT-specific overrides ----
    // BLAT: keep mapq off (hifi sets 0.001; explicit here so future
    // hifi changes don't silently flip it).
    m.mapq_score_scale = 0.001;
    // BLAT: report many hits per read instead of one best.
    // hifi: max_alignments=3, max_multimaps not set (=1).
    // blat: 100 each so we surface every haplotype/paralog hit.
    // max_alignments must be >= max_multimaps or it caps reporting.
    m.max_multimaps = 100;
    m.max_alignments = 100;
    // BLAT: loosen chain filtering so secondary hits aren't pruned.
    // hifi: chain_score_threshold=200 (chains within 200 of best).
    // blat: 1000 — admit much weaker chains relative to the best.
    m.chain_score_threshold = 1000.0;
    // BLAT: drop the per-base chain score floor.
    // hifi: 0.1 (chain must score >= 0.1 per read base).
    // blat: 0.05 — let noisy/partial matches survive.
    m.min_chain_score_per_base = 0.05;
    // BLAT: always try many chains, even when one dominates.
    // hifi: min_chains=2. blat: 50, so paralog/cross-haplotype hits
    // surface even when one chain has a clear best score.
    m.min_chains = 50;
    // BLAT: allow many chains per zipcode-tree cluster.
    // hifi: max_chains_per_tree=3. blat: 20 — tandem repeats and
    // high-copy gene families need more than 3 to produce all hits.
    m.max_chains_per_tree = 20;
}

/// Convert a SAM-style cigar vector into a string like "55M". Returns "*"
/// for the empty CIGAR so the tag is unambiguous in framed output.
string cigar_to_string(const vector<pair<int, char>>& cigar) {
    if (cigar.empty()) {
        return "*";
    }
    stringstream ss;
    for (const auto& op : cigar) {
        ss << op.first << op.second;
    }
    return ss.str();
}

} // namespace

HaplotypeSurjectionResult HaplotypeSurjector::surject(
    const Alignment& aln,
    const string& target_path_name,
    const vector<gbwt::size_type>& candidate_seq_ids) const
{
    HaplotypeSurjectionResult result;
    using S = HaplotypeSurjectionResult::Status;

    if (aln.path().mapping_size() == 0) {
        result.status = S::EMPTY_INPUT;
        return result;
    }
    if (target_path_name.empty()) {
        result.status = S::UNKNOWN_PATH;
        return result;
    }

    // The position graph only reports has_path() for paths that are actually
    // indexed for surjection. Anything else is "not indexed" even if the
    // underlying GBZ knows about the path under a different sense.
    if (!position_graph_->has_path(target_path_name)) {
        if (gbz_->graph.has_path(target_path_name)) {
            result.status = S::PATH_NOT_INDEXED;
        } else {
            result.status = S::UNKNOWN_PATH;
        }
        return result;
    }

    // Translate the requested path name to its GBWT path id, then to the
    // forward and backward GBWT sequence ids that represent it. The
    // assigner returns sequence ids in either orientation depending on the
    // alignment's traversal of the path, so we accept either as
    // "compatible". Path handles from the GBZ's GBWTGraph are stable across
    // the overlay (the overlay doc explicitly preserves them).
    path_handle_t gbz_path_handle = gbz_->graph.get_path_handle(target_path_name);
    gbwt::size_type gbwt_path_id = gbz_->graph.handle_to_path(gbz_path_handle);
    gbwt::size_type fwd_seq = 2 * gbwt_path_id;
    gbwt::size_type rev_seq = fwd_seq + 1;

    bool compatible = false;
    for (gbwt::size_type seq : candidate_seq_ids) {
        if (seq == fwd_seq || seq == rev_seq) {
            compatible = true;
            break;
        }
    }
    if (!compatible) {
        result.status = S::INCOMPATIBLE;
        return result;
    }

    // Run vg's Surjector against just this one path. Default settings: no
    // softclip suppression, deletions not preserved, no spliced semantics.
    unordered_set<path_handle_t> paths;
    paths.insert(gbz_path_handle);
    vector<tuple<string, int64_t, bool>> positions_out;
    vector<Alignment> surjected = surjector_.surject(aln, paths, positions_out,
                                                     /*allow_negative_scores=*/false,
                                                     /*preserve_deletions=*/false);

    if (surjected.empty() || positions_out.empty()) {
        result.status = S::SURJECTION_FAILED;
        return result;
    }

    // Pick the highest-scoring surjected alignment when there's more than
    // one (multimap_to_all_paths is off by default but be defensive).
    size_t best = 0;
    for (size_t i = 1; i < surjected.size(); ++i) {
        if (surjected[i].score() > surjected[best].score()) {
            best = i;
        }
    }
    const Alignment& chosen = surjected[best];
    const auto& chosen_pos = positions_out[best];

    if (chosen.path().mapping_size() == 0 || std::get<0>(chosen_pos).empty()) {
        result.status = S::SURJECTION_FAILED;
        return result;
    }

    result.status = S::OK;
    result.path_name = std::get<0>(chosen_pos);
    result.path_position = std::get<1>(chosen_pos);
    result.path_reverse = std::get<2>(chosen_pos);
    result.score = chosen.score();
    result.mapping_quality = chosen.mapping_quality();

    // Build the SAM-style CIGAR. cigar_against_path may adjust pos when
    // soft-clip suppression is enabled; we pass 0 to disable that, so the
    // returned position equals the supplied input.
    int64_t pos_adj = result.path_position;
    path_handle_t overlay_path_handle = position_graph_->get_path_handle(result.path_name);
    size_t path_len = position_graph_->get_path_length(overlay_path_handle);
    auto cigar_ops = cigar_against_path(chosen, result.path_reverse, pos_adj,
                                        path_len, /*softclip_suppress=*/0);
    result.cigar = cigar_to_string(cigar_ops);

    return result;
}

namespace {

/// Collapse raw GBWT sequence ids to the distinct, sorted set of haplotype
/// names (sample#phase#contig) via GBWT metadata. The raw ids double-count
/// orientations, repeats, and fragments, so a conserved locus can list
/// hundreds of numbers that are really a few dozen haplotypes. Falls back to
/// the numeric id when metadata is unavailable.
std::vector<string> compute_haplotype_names(const gbwtgraph::GBZ& gbz,
                                            const std::vector<gbwt::size_type>& seq_ids) {
    std::set<string> uniq;
    const gbwt::GBWT& index = gbz.index;
    const bool have_md = index.hasMetadata();
    const gbwt::Metadata& md = index.metadata;
    for (gbwt::size_type s : seq_ids) {
        // GBWT sequence id -> path id (forward/reverse share a path).
        gbwt::size_type path_id = s / 2;
        if (!have_md || path_id >= md.paths()) {
            uniq.insert(std::to_string(s));  // no metadata: fall back to id
            continue;
        }
        gbwt::PathName pn = md.path(path_id);
        // Same name construction as build_translation_tables.cpp's
        // build_base_name: sample#phase#contig, with safe fallbacks.
        string sample = (md.samples() > 0 && pn.sample < md.samples())
                        ? md.sample(pn.sample) : std::to_string(pn.sample);
        string contig = (md.contigs() > 0 && pn.contig < md.contigs())
                        ? md.contig(pn.contig) : std::to_string(pn.contig);
        uniq.insert(sample + "#" + std::to_string(pn.phase) + "#" + contig);
    }
    return std::vector<string>(uniq.begin(), uniq.end());
}

string join_csv(const std::vector<string>& items) {
    string out;
    for (size_t i = 0; i < items.size(); ++i) {
        if (i) out.push_back(',');
        out += items[i];
    }
    return out;
}

void attach_haplotype_annotations(Alignment& aln, const HaplotypeAssignmentResult& r) {
    // Report deduplicated haplotype NAMES (+ their count) rather than the raw,
    // inflated sequence-id list.
    set_annotation(aln, "haplotype_names", join_csv(r.haplotype_names));
    set_annotation(aln, "haplotype_distinct_count", static_cast<double>(r.haplotype_names.size()));
    set_annotation(aln, "haplotype_candidate_count", static_cast<double>(r.candidate_count));
    set_annotation(aln, "haplotype_left_boundary", static_cast<double>(r.left_boundary_node));
    set_annotation(aln, "haplotype_right_boundary", static_cast<double>(r.right_boundary_node));
    set_annotation(aln, "haplotype_truncated", r.truncated);
}

/// Append GAF optional tags describing a surjection result. We always emit
/// the status tag (sj:Z:<token>); the remaining tags are only emitted on
/// success so downstream parsers can short-circuit.
void append_surjection_tags(stringstream& ss, const HaplotypeSurjectionResult& r) {
    ss << "\tsj:Z:" << haplotype_surjection_status_token(r.status);
    if (r.status != HaplotypeSurjectionResult::Status::OK) {
        return;
    }
    ss << "\tsn:Z:" << r.path_name
       << "\tsp:i:" << r.path_position
       << "\tsr:i:" << (r.path_reverse ? 1 : 0)
       << "\tss:i:" << r.score
       << "\tsm:i:" << r.mapping_quality
       << "\tsc:Z:" << (r.cigar.empty() ? string("*") : r.cigar);
}

} // namespace

void GiraffeEngine::load(const GiraffeEnginePaths& paths, const GiraffeEngineConfig& config) {
    if (paths.gbz_path.empty() || paths.minimizer_path.empty() || paths.distance_path.empty() || paths.zipcode_path.empty()) {
        throw runtime_error("All index paths must be provided: GBZ, minimizer, distance, zipcodes");
    }
    if (config.threads == 0) {
        throw runtime_error("Thread count must be >= 1");
    }
    // max_multimaps == 0 is allowed and means "use the preset / mapper
    // default". A positive value is an explicit override (see below).
    if (!file_exists(paths.gbz_path) || !file_exists(paths.minimizer_path)
        || !file_exists(paths.distance_path) || !file_exists(paths.zipcode_path)) {
        throw runtime_error("One or more required index files do not exist");
    }

    auto gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(paths.gbz_path);
    auto minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(paths.minimizer_path);
    auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(paths.distance_path);
    auto zipcodes = make_unique<ZipCodeCollection>();

    ifstream zip_in(paths.zipcode_path);
    if (!zip_in) {
        throw runtime_error("Failed to open zipcode file: " + paths.zipcode_path);
    }
    zipcodes->deserialize(zip_in);
    zip_in.close();

    require_payload(*minimizer_index, MinimizerIndexParameters::PAYLOAD_ZIPCODES);
    require_compatible_graphs(*gbz, "GBZ", *minimizer_index, "Minimizer Index");

    if (config.preload_distance_index) {
        distance_index->preload(true);
    }

    auto mapper = make_unique<MinimizerMapper>(gbz->graph, *minimizer_index, distance_index.get(), zipcodes.get());

    // The engine always uses the long-read/BLAT mapping configuration — per
    // Benedict, the long-read code path is intended to replace the
    // short-read path and works well enough for short reads too. Applied
    // before explicit overrides so config.max_multimaps (if > 0) wins.
    apply_blat_preset(*mapper);
    if (config.max_multimaps > 0) {
        mapper->max_multimaps = config.max_multimaps;
    }

    auto assigner = make_unique<HaplotypeAssigner>(*gbz, *distance_index);

    // Build the path-position overlay over the GBZ graph if the caller
    // configured target haplotype paths to surject onto. We do this once at
    // load time so per-read surjection is just a position-graph query plus
    // the surjector's realigning DP. Validating that every requested path
    // actually exists up front makes server misconfiguration fail loudly.
    unique_ptr<bdsg::ReferencePathOverlay> position_overlay;
    unique_ptr<HaplotypeSurjector> surjector;
    if (!config.surjection_target_paths.empty()) {
        for (const auto& name : config.surjection_target_paths) {
            if (!gbz->graph.has_path(name)) {
                throw runtime_error(
                    "Surjection target path not present in GBZ: " + name);
            }
        }
        std::unordered_set<string> extras(config.surjection_target_paths.begin(),
                                          config.surjection_target_paths.end());
        position_overlay = make_unique<bdsg::ReferencePathOverlay>(
            &gbz->graph, extras, /*all_paths=*/false);
        surjector = make_unique<HaplotypeSurjector>(*gbz, *position_overlay);
    }

    // Keep the same behavior as vg giraffe: global OMP thread count controls mapping parallelism.
    omp_set_num_threads(static_cast<int>(config.threads));

    config_ = config;
    gbz_ = move(gbz);
    minimizer_index_ = move(minimizer_index);
    distance_index_ = move(distance_index);
    zipcodes_ = move(zipcodes);
    mapper_ = move(mapper);
    haplotype_assigner_ = move(assigner);
    position_overlay_ = move(position_overlay);
    haplotype_surjector_ = move(surjector);
}

bool GiraffeEngine::is_loaded() const {
    return static_cast<bool>(mapper_);
}

const GiraffeEngineConfig& GiraffeEngine::config() const {
    return config_;
}

vector<string> GiraffeEngine::gaf_header_lines() const {
    require_loaded();
    gbwtgraph::GraphName graph_name = gbz_->graph_name();
    return graph_name.gaf_header_lines();
}

vector<vector<string>> GiraffeEngine::map_reads(const vector<GiraffeFastqRead>& reads) {
    require_loaded();
    vector<vector<string>> results(reads.size());
    if (reads.empty()) {
        return results;
    }
    for (const auto& read : reads) {
        if (!read.quality.empty() && read.quality.size() != read.sequence.size()) {
            throw runtime_error("Read " + read.name + " has mismatched sequence/quality lengths");
        }
    }
    exception_ptr worker_error = nullptr;
    atomic<bool> has_worker_error(false);
    mutex worker_error_mutex;

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < reads.size(); ++i) {
        if (has_worker_error.load(memory_order_acquire)) {
            continue;
        }
        const auto& read = reads[i];
        try {
            Alignment aln;
            aln.set_name(read.name);
            aln.set_sequence(read.sequence);
            aln.set_quality(read.quality.empty() ? string(read.sequence.size(), 'I') : read.quality);
            check_quality_length(aln);
            toUppercaseInPlace(*aln.mutable_sequence());

            vector<Alignment> mapped = mapper_->map(aln);
            vector<HaplotypeAssignmentResult> hap_results;
            if (this->assign_haplotypes && this->haplotype_assigner_) {
                hap_results.reserve(mapped.size());
                for (auto& m : mapped) {
                    auto hap = this->haplotype_assigner_->assign(
                        m, this->assign_haplotypes_extend_to_snarls);
                    // Resolve raw sequence ids to distinct haplotype names once,
                    // here, so both the annotation and the hp/hn tags reuse it.
                    hap.haplotype_names = compute_haplotype_names(*gbz_, hap.seq_ids);
                    attach_haplotype_annotations(m, hap);
                    hap_results.emplace_back(std::move(hap));
                }
            }

            // If the caller requested surjection for this read, attach a
            // per-alignment result for every reported mapping. We always
            // emit a status when a target was requested, so the client can
            // tell the difference between "engine has no target indexed"
            // and "the requested haplotype isn't compatible with this
            // alignment".
            vector<HaplotypeSurjectionResult> surj_results;
            if (!read.surjection_target.empty()) {
                surj_results.resize(mapped.size());
                if (!this->haplotype_surjector_) {
                    for (auto& r : surj_results) {
                        r.status = HaplotypeSurjectionResult::Status::PATH_NOT_INDEXED;
                    }
                } else if (hap_results.empty()) {
                    // The compatibility check needs the assigner output.
                    // Without it we can't safely surject.
                    for (auto& r : surj_results) {
                        r.status = HaplotypeSurjectionResult::Status::INCOMPATIBLE;
                    }
                } else {
                    for (size_t j = 0; j < mapped.size(); ++j) {
                        surj_results[j] = this->haplotype_surjector_->surject(
                            mapped[j], read.surjection_target, hap_results[j].seq_ids);
                    }
                }
            }

            auto& lines = results[i];
            lines.reserve(mapped.size());
            for (size_t j = 0; j < mapped.size(); ++j) {
                auto& m = mapped[j];
                auto gaf = io::alignment_to_gaf(gbz_->graph, m);
                stringstream ss;
                ss << gaf;
                if (!hap_results.empty()) {
                    const auto& hap = hap_results[j];
                    // hp:Z: = deduplicated haplotype NAMES (was raw seq ids);
                    // hn:i: = number of distinct haplotypes; hc:i: = raw
                    // candidate occurrences (pre-cap, for reference).
                    ss << "\thp:Z:" << join_csv(hap.haplotype_names)
                       << "\thn:i:" << hap.haplotype_names.size()
                       << "\thc:i:" << hap.candidate_count
                       << "\thl:i:" << hap.left_boundary_node
                       << "\thr:i:" << hap.right_boundary_node
                       << "\tht:i:" << (hap.truncated ? 1 : 0);
                }
                if (j < surj_results.size()) {
                    append_surjection_tags(ss, surj_results[j]);
                }
                lines.emplace_back(ss.str());
            }
        } catch (...) {
            lock_guard<mutex> guard(worker_error_mutex);
            if (worker_error == nullptr) {
                worker_error = current_exception();
                has_worker_error.store(true, memory_order_release);
            }
        }
    }
    if (worker_error != nullptr) {
        rethrow_exception(worker_error);
    }

    return results;
}

void GiraffeEngine::require_loaded() const {
    if (!mapper_) {
        throw runtime_error("GiraffeEngine is not loaded");
    }
}

namespace {

/// Construct a libhandlegraph step_handle_t from a (node, offset_in_record)
/// pair using gbwtgraph's documented step encoding (gbwtgraph.cpp:911-1024:
/// as_integers(step)[0] = node, [1] = offset).
handlegraph::step_handle_t make_step_handle(uint64_t node, uint64_t offset) {
    handlegraph::step_handle_t step;
    handlegraph::as_integers(step)[0] = node;
    handlegraph::as_integers(step)[1] = offset;
    return step;
}

/// Walk the target's GBWT path summing node lengths. Used as a fallback when
/// the caller doesn't pass a precomputed target_path_length.
size_t walk_target_path_length(const gbwtgraph::GBZ& gbz, size_t target_gbwt_path_id) {
    gbwt::vector_type nodes = gbz.index.extract(
        gbwt::Path::encode(target_gbwt_path_id, false));
    size_t total = 0;
    for (gbwt::node_type node : nodes) {
        if (node == gbwt::ENDMARKER) break;
        total += gbz.graph.get_length(
            gbz.graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
    }
    return total;
}

} // anonymous namespace

std::vector<std::string> GiraffeEngine::surject_with_anchors(
    const std::string& read_name,
    const std::string& gaf_line,
    const std::vector<WireAnchor>& anchors,
    const std::string& target_haplotype,
    size_t target_path_length)
{
    require_loaded();
    using S = HaplotypeSurjectionResult::Status;

    if (gaf_line.empty()) {
        throw runtime_error("surject_with_anchors: empty GAF input");
    }

    // 1. Parse GAF → vg::Alignment using the GBZ's graph for node lengths.
    gafkluge::GafRecord gaf_record;
    gafkluge::parse_gaf_record(gaf_line, gaf_record);
    Alignment aln;
    vg::io::gaf_to_alignment(gbz_->graph, gaf_record, aln);
    // Preserve the caller-supplied read name on the resulting Alignment.
    if (!read_name.empty()) {
        aln.set_name(read_name);
    } else if (aln.name().empty()) {
        aln.set_name(gaf_record.query_name);
    }

    // 2. Validate target haplotype. The caller may pass either a full GBZ
    //    path name (e.g. "GRCh38#0#chrM") or a haplotype prefix that T1
    //    resolved during build_surject_anchors (e.g. "HG002#2"). For the
    //    second case the GBZ doesn't know "HG002#2" as a path, so we
    //    recover the concrete path name from the first anchor's step_handle
    //    — which was pinned during build to exactly one GBWT path id.
    HaplotypeSurjectionResult surj_result;
    std::string resolved_target = target_haplotype;
    path_handle_t target_path{};
    bool target_resolved = false;

    if (aln.path().mapping_size() == 0) {
        surj_result.status = S::EMPTY_INPUT;
    } else if (gbz_->graph.has_path(target_haplotype)) {
        target_path = gbz_->graph.get_path_handle(target_haplotype);
        target_resolved = true;
    } else if (!anchors.empty()) {
        // Anchors carry step_handles in gbwtgraph's native (node,
        // offset_in_record) encoding; get_path_handle_of_step decodes that
        // back to the path the anchor pinned during the build step.
        handlegraph::step_handle_t s = make_step_handle(
            anchors[0].step_begin_node, anchors[0].step_begin_offset);
        try {
            target_path = gbz_->graph.get_path_handle_of_step(s);
            resolved_target = gbz_->graph.get_path_name(target_path);
            target_resolved = !resolved_target.empty();
        } catch (const std::exception&) {
            target_resolved = false;
        }
        if (!target_resolved) {
            surj_result.status = S::UNKNOWN_PATH;
        }
    } else {
        surj_result.status = S::UNKNOWN_PATH;
    }

    if (!target_resolved || surj_result.status != S::OK) {
        if (surj_result.status == S::OK) {
            surj_result.status = S::UNKNOWN_PATH;
        }
        stringstream ss;
        ss << gaf_line;
        ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status);
        return {ss.str()};
    }

    // 3. Determine the target path length (cached by the caller if possible).
    size_t path_length = target_path_length;
    if (path_length == 0) {
        gbwt::size_type gbwt_path_id = gbz_->graph.handle_to_path(target_path);
        path_length = walk_target_path_length(*gbz_, gbwt_path_id);
    }

    // 4. Translate wire anchors → vg::PrecomputedAnchor records. The
    //    anchor-backed position graph only needs step_handles + base offsets.
    std::vector<PrecomputedAnchor> vg_anchors;
    vg_anchors.reserve(anchors.size());
    for (const auto& a : anchors) {
        PrecomputedAnchor pa;
        pa.step_begin = make_step_handle(a.step_begin_node, a.step_begin_offset);
        pa.step_end   = make_step_handle(a.step_end_node, a.step_end_offset);
        pa.path_offset_step_begin = a.path_offset_step_begin;
        pa.path_offset_step_end   = a.path_offset_step_end;
        vg_anchors.push_back(pa);
    }

    // 5. Build the anchor-backed position graph and a Surjector over it.
    AnchorBackedPositionGraph anchor_graph(&gbz_->graph, vg_anchors,
                                           target_path, path_length);
    Surjector anchor_surjector(&anchor_graph);

    unordered_set<path_handle_t> paths;
    paths.insert(target_path);
    vector<tuple<string, int64_t, bool>> positions_out;
    vector<Alignment> surjected;
    try {
        surjected = anchor_surjector.surject(aln, paths, positions_out,
                                             /*allow_negative_scores=*/false,
                                             /*preserve_deletions=*/false);
    } catch (const std::exception& e) {
        surj_result.status = S::SURJECTION_FAILED;
        stringstream ss;
        ss << gaf_line;
        ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status);
        ss << "\tse:Z:" << e.what();
        return {ss.str()};
    }

    if (surjected.empty() || positions_out.empty()) {
        surj_result.status = S::SURJECTION_FAILED;
        stringstream ss;
        ss << gaf_line;
        ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status);
        return {ss.str()};
    }

    // 6. Pick best (highest-scoring) surjected alignment, mirroring HaplotypeSurjector::surject.
    size_t best = 0;
    for (size_t i = 1; i < surjected.size(); ++i) {
        if (surjected[i].score() > surjected[best].score()) {
            best = i;
        }
    }
    const Alignment& chosen = surjected[best];
    const auto& chosen_pos = positions_out[best];

    if (chosen.path().mapping_size() == 0 || std::get<0>(chosen_pos).empty()) {
        surj_result.status = S::SURJECTION_FAILED;
        stringstream ss;
        ss << gaf_line;
        ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status);
        return {ss.str()};
    }

    surj_result.status = S::OK;
    surj_result.path_name = std::get<0>(chosen_pos);
    surj_result.path_position = std::get<1>(chosen_pos);
    surj_result.path_reverse = std::get<2>(chosen_pos);
    surj_result.score = chosen.score();
    surj_result.mapping_quality = chosen.mapping_quality();

    int64_t pos_adj = surj_result.path_position;
    size_t plen = anchor_graph.get_path_length(target_path);
    auto cigar_ops = cigar_against_path(chosen, surj_result.path_reverse, pos_adj,
                                        plen, /*softclip_suppress=*/0);
    surj_result.cigar = cigar_to_string(cigar_ops);

    // 7. Append surjection tags to the original GAF line and return one line.
    //    Matches the schema map_reads + HaplotypeSurjector::surject emits so
    //    callers don't have to special-case the anchor path.
    stringstream ss;
    ss << gaf_line;
    ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status)
       << "\tsn:Z:" << surj_result.path_name
       << "\tsp:i:" << surj_result.path_position
       << "\tsr:i:" << (surj_result.path_reverse ? 1 : 0)
       << "\tss:i:" << surj_result.score
       << "\tsm:i:" << surj_result.mapping_quality
       << "\tsc:Z:" << (surj_result.cigar.empty() ? string("*") : surj_result.cigar)
       << "\tan:i:" << anchors.size()
       << "\tap:Z:" << resolved_target;
    return {ss.str()};
}

} // namespace vg
