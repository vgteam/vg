#include "giraffe_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
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
#include <gbwtgraph/utils.h>   // sample_name_set, GENERIC_PATH_SAMPLE_NAME

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

namespace {

/// Count matches / substitutions / indels from a surjected alignment's edits.
/// vg edit semantics: from_length == to_length with an empty sequence is a run
/// of matches; equal lengths WITH a sequence is a substitution; to_length only
/// is an insertion in the read; from_length only is a deletion from the read.
/// The SAM CIGAR collapses the first two into 'M', so this is the only place
/// identity can be measured.
void count_edits(const Alignment& aln, uint64_t& matches, uint64_t& mismatches,
                 uint64_t& inserted, uint64_t& deleted) {
    matches = mismatches = inserted = deleted = 0;
    for (const auto& mapping : aln.path().mapping()) {
        for (const auto& edit : mapping.edit()) {
            if (edit.from_length() == edit.to_length()) {
                if (edit.sequence().empty()) matches += edit.from_length();
                else                         mismatches += edit.from_length();
            } else if (edit.from_length() == 0) {
                inserted += edit.to_length();
            } else if (edit.to_length() == 0) {
                deleted += edit.from_length();
            } else {
                // Mixed-length edit: the shared prefix is substituted, the
                // remainder is an indel.
                const uint64_t common = std::min(edit.from_length(), edit.to_length());
                mismatches += common;
                if (edit.to_length() > common) inserted += edit.to_length() - common;
                else                           deleted  += edit.from_length() - common;
            }
        }
    }
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

/// Fill result from the highest-scoring of the Surjector's outputs on path,
/// a path of position_graph. Sets SURJECTION_FAILED if there is none.
void take_best_surjection(const vector<Alignment>& surjected,
                          const vector<tuple<string, int64_t, bool>>& positions,
                          const PathPositionHandleGraph& position_graph,
                          const path_handle_t& path,
                          HaplotypeSurjectionResult& result) {
    using S = HaplotypeSurjectionResult::Status;
    if (surjected.empty() || positions.empty()) {
        result.status = S::SURJECTION_FAILED;
        return;
    }
    size_t best = 0;
    for (size_t i = 1; i < surjected.size(); ++i) {
        if (surjected[i].score() > surjected[best].score()) {
            best = i;
        }
    }
    const Alignment& chosen = surjected[best];
    const auto& chosen_pos = positions[best];
    if (chosen.path().mapping_size() == 0 || std::get<0>(chosen_pos).empty()) {
        result.status = S::SURJECTION_FAILED;
        return;
    }

    result.status = S::OK;
    result.path_name = std::get<0>(chosen_pos);
    result.path_position = std::get<1>(chosen_pos);
    result.path_reverse = std::get<2>(chosen_pos);
    result.score = chosen.score();
    result.mapping_quality = chosen.mapping_quality();

    // softclip_suppress = 0 leaves the position as given.
    int64_t pos_adj = result.path_position;
    size_t path_len = position_graph.get_path_length(path);
    result.cigar = cigar_to_string(cigar_against_path(chosen, result.path_reverse, pos_adj,
                                                      path_len, /*softclip_suppress=*/0));
    count_edits(chosen, result.matches, result.mismatches,
                result.inserted, result.deleted);
}

} // namespace

HaplotypeSurjectionResult GiraffeEngine::surject_to_indexed_path(
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
    if (!position_overlay_ || !position_overlay_->has_path(target_path_name)) {
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
    vector<Alignment> surjected = overlay_surjector_->surject(aln, paths, positions_out,
                                                              /*allow_negative_scores=*/false,
                                                              /*preserve_deletions=*/false);
    take_best_surjection(surjected, positions_out, *position_overlay_,
                         position_overlay_->get_path_handle(target_path_name), result);
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

/// Build the sample#phase#contig name for a GBWT path id (safe fallbacks).
string path_name_for(const gbwt::Metadata& md, gbwt::size_type path_id) {
    gbwt::PathName pn = md.path(path_id);
    string sample = (md.samples() > 0 && pn.sample < md.samples())
                    ? md.sample(pn.sample) : std::to_string(pn.sample);
    string contig = (md.contigs() > 0 && pn.contig < md.contigs())
                    ? md.contig(pn.contig) : std::to_string(pn.contig);
    return sample + "#" + std::to_string(pn.phase) + "#" + contig;
}

/// Break the (often large) tie among haplotypes that carry the read's path
/// identically by choosing one representative. Policy:
///   1. prefer a haplotype whose sample is a reference sample (or the generic
///      backbone), then
///   2. among the eligible class, take the lowest GBWT sequence id
///      (deterministic and reproducible).
/// Returns (representative_name, chosen_because_reference). Empty name if there
/// are no candidates or no metadata.
std::pair<string, bool> choose_representative(const gbwtgraph::GBZ& gbz,
                                              const std::vector<gbwt::size_type>& seq_ids) {
    const gbwt::GBWT& index = gbz.index;
    if (seq_ids.empty() || !index.hasMetadata()) {
        return {string(), false};
    }
    const gbwt::Metadata& md = index.metadata;
    const gbwtgraph::sample_name_set& refs = gbz.get_reference_samples();

    bool any = false, ref_found = false;
    gbwt::size_type best_any = 0, best_ref = 0;
    for (gbwt::size_type s : seq_ids) {
        gbwt::size_type path_id = s / 2;
        if (path_id >= md.paths()) continue;
        gbwt::PathName pn = md.path(path_id);
        string sample = (md.samples() > 0 && pn.sample < md.samples())
                        ? md.sample(pn.sample) : string();
        bool is_ref = !sample.empty() &&
                      (refs.count(sample) > 0 ||
                       sample == gbwtgraph::GENERIC_PATH_SAMPLE_NAME);
        if (!any || s < best_any) { best_any = s; any = true; }
        if (is_ref && (!ref_found || s < best_ref)) { best_ref = s; ref_found = true; }
    }
    if (!any) {
        return {string(), false};
    }
    gbwt::size_type chosen = ref_found ? best_ref : best_any;
    return {path_name_for(md, chosen / 2), ref_found};
}

void attach_haplotype_annotations(Alignment& aln, const HaplotypeAssignmentResult& r) {
    // Report deduplicated haplotype NAMES (+ their count) rather than the raw,
    // inflated sequence-id list.
    set_annotation(aln, "haplotype_names", join_csv(r.haplotype_names));
    set_annotation(aln, "haplotype_distinct_count", static_cast<double>(r.haplotype_names.size()));
    set_annotation(aln, "haplotype_representative", r.representative_name);
    set_annotation(aln, "haplotype_representative_is_reference", r.representative_is_reference);
    set_annotation(aln, "haplotype_candidate_count", static_cast<double>(r.candidate_count));
    if (r.left_boundary_node != 0) {
        set_annotation(aln, "haplotype_left_boundary", static_cast<double>(r.left_boundary_node));
    }
    if (r.right_boundary_node != 0) {
        set_annotation(aln, "haplotype_right_boundary", static_cast<double>(r.right_boundary_node));
    }
    set_annotation(aln, "haplotype_truncated", r.seq_ids.size() < r.candidate_count);
    // Path-coverage / mosaic accounting: whether one haplotype carries the whole
    // read, how much of it the reported (best) block covers, and how many
    // maximal haplotype-consistent segments the read decomposes into.
    set_annotation(aln, "haplotype_fully_covered", r.fully_covered);
    set_annotation(aln, "haplotype_best_covered_bp", static_cast<double>(r.best_covered_bp));
    set_annotation(aln, "haplotype_path_length_bp", static_cast<double>(r.path_length_bp));
    set_annotation(aln, "haplotype_num_segments", static_cast<double>(r.segments.size()));
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
       << "\tsc:Z:" << (r.cigar.empty() ? string("*") : r.cigar)
       << "\tsM:i:" << r.matches
       << "\tsX:i:" << r.mismatches
       << "\tsI:i:" << r.inserted
       << "\tsD:i:" << r.deleted;
}

} // namespace

void GiraffeEngine::load(const GiraffeEngineConfig& config) {
    if (config.gbz_path.empty() || config.minimizer_path.empty() || config.distance_path.empty() || config.zipcode_path.empty()) {
        throw runtime_error("All index paths must be provided: GBZ, minimizer, distance, zipcodes");
    }
    if (config.threads == 0) {
        throw runtime_error("Thread count must be >= 1");
    }
    // max_multimaps == 0 is allowed and means "use the preset / mapper
    // default". A positive value is an explicit override (see below).
    if (!file_exists(config.gbz_path) || !file_exists(config.minimizer_path)
        || !file_exists(config.distance_path) || !file_exists(config.zipcode_path)) {
        throw runtime_error("One or more required index files do not exist");
    }

    auto gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(config.gbz_path);
    auto minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(config.minimizer_path);
    auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(config.distance_path);
    auto zipcodes = make_unique<ZipCodeCollection>();

    ifstream zip_in(config.zipcode_path);
    if (!zip_in) {
        throw runtime_error("Failed to open zipcode file: " + config.zipcode_path);
    }
    zipcodes->deserialize(zip_in);
    zip_in.close();

    if (!has_payload(*minimizer_index, {MinimizerIndexParameters::PAYLOAD_ZIPCODES,
                                        MinimizerIndexParameters::PAYLOAD_ZIPCODES_WITH_PATHS})) {
        throw runtime_error("Minimizer index has no zipcode payload: " + config.minimizer_path);
    }
    require_compatible_graphs(*gbz, "GBZ", *minimizer_index, "Minimizer Index");

    if (config.preload_distance_index) {
        distance_index->preload(true);
    }

    auto mapper = make_unique<MinimizerMapper>(gbz->graph, *minimizer_index, distance_index.get(), zipcodes.get());

    // Mapping parameters come from the caller (giraffe-server takes them from
    // vg giraffe's presets), so the engine keeps no copy of preset values.
    if (config.configure_mapper) {
        config.configure_mapper(*mapper, *minimizer_index);
    }
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
    unique_ptr<Surjector> surjector;
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
            &gbz->graph, /*all_paths=*/false, extras);
        surjector = make_unique<Surjector>(position_overlay.get());
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
    overlay_surjector_ = move(surjector);
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

vector<vector<string>> GiraffeEngine::map_reads(const vector<GiraffeServerRequest>& reads) {
    require_loaded();
    vector<vector<string>> results(reads.size());
    if (reads.empty()) {
        return results;
    }
    for (const auto& req : reads) {
        if (!req.read.quality().empty() && req.read.quality().size() != req.read.sequence().size()) {
            throw runtime_error("Read " + req.read.name() + " has mismatched sequence/quality lengths");
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
        const auto& req = reads[i];
        try {
            Alignment aln = req.read;
            check_quality_length(aln);
            toUppercaseInPlace(*aln.mutable_sequence());

            vector<Alignment> mapped = mapper_->map(aln);
            vector<HaplotypeAssignmentResult> hap_results;
            if (config_.assign_haplotypes && this->haplotype_assigner_) {
                hap_results.reserve(mapped.size());
                for (auto& m : mapped) {
                    auto hap = this->haplotype_assigner_->assign(
                        m, config_.assign_haplotypes_extend_to_snarls);
                    // Name every mosaic segment and pick each one's representative
                    // (reference-first, else lowest id). For a fully-covered read
                    // there is a single segment spanning the whole path.
                    for (auto& seg : hap.segments) {
                        seg.haplotype_names = compute_haplotype_names(*gbz_, seg.seq_ids);
                        seg.representative_name = choose_representative(*gbz_, seg.seq_ids).first;
                    }
                    // Promote the best (longest) segment's names/representative to
                    // the top level so the annotation and hp/hn/hb tags describe
                    // the haplotype(s) covering the most of the read path.
                    if (hap.best_segment < hap.segments.size()) {
                        const auto& best = hap.segments[hap.best_segment];
                        hap.haplotype_names = best.haplotype_names;
                        auto rep = choose_representative(*gbz_, best.seq_ids);
                        hap.representative_name = rep.first;
                        hap.representative_is_reference = rep.second;
                    }
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
            if (!req.surjection_target.empty()) {
                surj_results.resize(mapped.size());
                if (!overlay_surjector_) {
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
                        surj_results[j] = surject_to_indexed_path(
                            mapped[j], req.surjection_target, hap_results[j].seq_ids);
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
                       << "\thb:Z:" << (hap.representative_name.empty()
                                        ? string("*") : hap.representative_name)
                       << "\thf:i:" << (hap.representative_is_reference ? 1 : 0)
                       << "\thc:i:" << hap.candidate_count;
                    // hl:i: / hr:i: snarl-boundary nodes, only when one was found.
                    if (hap.left_boundary_node != 0) {
                        ss << "\thl:i:" << hap.left_boundary_node;
                    }
                    if (hap.right_boundary_node != 0) {
                        ss << "\thr:i:" << hap.right_boundary_node;
                    }
                    ss << "\tht:i:" << (hap.seq_ids.size() < hap.candidate_count ? 1 : 0)
                       // Path-coverage / mosaic tags:
                       //   hq:i: 1 if one haplotype carries the whole read, else 0
                       //   hv:i: % of the read path the reported (best) block covers
                       //   hs:i: number of maximal haplotype-consistent segments
                       << "\thq:i:" << (hap.fully_covered ? 1 : 0)
                       << "\thv:i:" << (hap.path_length_bp
                                        ? (100 * hap.best_covered_bp / hap.path_length_bp)
                                        : 0)
                       << "\ths:i:" << hap.segments.size();
                    // hm:Z: mosaic breakdown, emitted only when the read is NOT
                    // carried by a single haplotype (>1 segment). Each entry is
                    // "<covered_bp>:<distinct_haplotypes>:<representative>", in
                    // path order, separated by ';'.
                    if (hap.segments.size() > 1) {
                        ss << "\thm:Z:";
                        for (size_t s = 0; s < hap.segments.size(); ++s) {
                            const auto& seg = hap.segments[s];
                            ss << (s ? ";" : "")
                               << seg.covered_bp << ":"
                               << seg.haplotype_names.size() << ":"
                               << (seg.representative_name.empty()
                                   ? string("*") : seg.representative_name);
                        }
                    }
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

    // Optional phase timing for surjection performance debugging. Enable with
    // the env var GIRAFFE_SURJECT_TIMING=1; prints a per-read breakdown of
    // where time goes (parse / path-length walk / anchor build / surjector DP /
    // post) to stderr. No overhead when unset.
    static const bool surject_timing =
        (std::getenv("GIRAFFE_SURJECT_TIMING") != nullptr);
    using surj_clock = std::chrono::steady_clock;
    surj_clock::time_point t_begin   = surj_clock::now();
    surj_clock::time_point t_parse   = t_begin, t_pathlen = t_begin,
                           t_anchors = t_begin, t_build   = t_begin,
                           t_surject = t_begin;

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

    t_parse = surj_clock::now();

    // 3. Determine the target path length (cached by the caller if possible).
    size_t path_length = target_path_length;
    if (path_length == 0) {
        gbwt::size_type gbwt_path_id = gbz_->graph.handle_to_path(target_path);
        path_length = walk_target_path_length(*gbz_, gbwt_path_id);
    }
    t_pathlen = surj_clock::now();

    // Per-read timing reporter (no-op unless GIRAFFE_SURJECT_TIMING is set).
    auto report_timing = [&](const char* outcome) {
        if (!surject_timing) return;
        auto ms = [](surj_clock::time_point a, surj_clock::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };
        surj_clock::time_point now = surj_clock::now();
        std::cerr << "[surject-timing]"
                  << " read=" << aln.name()
                  << " read_len=" << aln.sequence().size()
                  << " n_anchors=" << anchors.size()
                  << " path_len=" << path_length
                  << " | parse="   << ms(t_begin,   t_parse)
                  << " pathlen="   << ms(t_parse,   t_pathlen)
                  << " anchors="   << ms(t_pathlen, t_anchors)
                  << " build="     << ms(t_anchors, t_build)
                  << " surject="   << ms(t_build,   t_surject)
                  << " post="      << ms(t_surject, now)
                  << " posq="      << g_anchor_graph_stats.pos_of_step_calls
                     << "/"        << g_anchor_graph_stats.pos_of_step_walk
                  << " stepq="     << g_anchor_graph_stats.step_at_pos_calls
                     << "/"        << g_anchor_graph_stats.step_at_pos_walk
                  << " total="     << ms(t_begin,   now) << "ms"
                  << " outcome="   << outcome << std::endl;
    };

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
    t_anchors = surj_clock::now();

    // 5. Build the anchor-backed position graph and a Surjector over it.
    AnchorBackedPositionGraph anchor_graph(&gbz_->graph, vg_anchors,
                                           target_path, path_length);
    Surjector anchor_surjector(&anchor_graph);
    t_build = surj_clock::now();

    unordered_set<path_handle_t> paths;
    paths.insert(target_path);
    vector<tuple<string, int64_t, bool>> positions_out;
    vector<Alignment> surjected;
    g_anchor_graph_stats = {};  // measure position-query cost during this surjection
    try {
        surjected = anchor_surjector.surject(aln, paths, positions_out,
                                             /*allow_negative_scores=*/false,
                                             /*preserve_deletions=*/false);
        t_surject = surj_clock::now();
    } catch (const std::exception& e) {
        t_surject = surj_clock::now();
        report_timing("surject_exception");
        surj_result.status = S::SURJECTION_FAILED;
        stringstream ss;
        ss << gaf_line;
        ss << "\tsj:Z:" << haplotype_surjection_status_token(surj_result.status);
        ss << "\tse:Z:" << e.what();
        return {ss.str()};
    }

    take_best_surjection(surjected, positions_out, anchor_graph, target_path, surj_result);
    stringstream ss;
    ss << gaf_line;
    append_surjection_tags(ss, surj_result);
    if (surj_result.status != S::OK) {
        report_timing("surject_empty");
        return {ss.str()};
    }
    ss << "\tan:i:" << anchors.size()
       << "\tap:Z:" << resolved_target;
    report_timing("ok");
    return {ss.str()};
}

} // namespace vg
