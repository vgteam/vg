#include "giraffe_engine.hpp"

#include <algorithm>
#include <fstream>
#include <atomic>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <omp.h>

#include <vg/io/alignment_io.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/vg.pb.h>

#include "alignment.hpp"
#include "annotation.hpp"
#include "gbwtgraph_helper.hpp"
#include "utility.hpp"

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

string haplotype_seq_ids_csv(const HaplotypeAssignmentResult& r) {
    string seq_ids_csv;
    seq_ids_csv.reserve(r.seq_ids.size() * 8);
    for (size_t i = 0; i < r.seq_ids.size(); ++i) {
        if (i) seq_ids_csv.push_back(',');
        seq_ids_csv += std::to_string(r.seq_ids[i]);
    }
    return seq_ids_csv;
}

void attach_haplotype_annotations(Alignment& aln, const HaplotypeAssignmentResult& r) {
    string seq_ids_csv = haplotype_seq_ids_csv(r);
    set_annotation(aln, "haplotype_seq_ids", seq_ids_csv);
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
    if (config.max_multimaps == 0) {
        throw runtime_error("max_multimaps must be >= 1");
    }
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
    mapper->max_multimaps = config.max_multimaps;

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
                    ss << "\thp:Z:" << haplotype_seq_ids_csv(hap)
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

} // namespace vg
