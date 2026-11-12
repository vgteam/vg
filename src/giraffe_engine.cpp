#include "giraffe_engine.hpp"

#include <algorithm>
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

/// Count the matched, mismatched, inserted and deleted bases in an alignment's
/// edits.
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
                // The edit changes both lengths. Count the shorter length as
                // mismatches and the rest as an insertion or deletion.
                const uint64_t common = std::min(edit.from_length(), edit.to_length());
                mismatches += common;
                if (edit.to_length() > common) inserted += edit.to_length() - common;
                else                           deleted  += edit.from_length() - common;
            }
        }
    }
}

/// Return a CIGAR as a string such as "55M", or "*" if it is empty.
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

/// Apply the engine's surjector settings. vg giraffe has four more, set by
/// flags the server does not have.
void set_up_surjector(Surjector& surjector, const GiraffeEngineConfig& config) {
    surjector.prune_suspicious_anchors = config.prune_suspicious_anchors;
}

/// Fill `result` from the highest-scoring surjected alignment on `path`.
/// `positions` holds the path name, start and strand of each alignment. Set
/// SURJECTION_FAILED if no alignment landed on the path.
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

    // cigar_against_path moves the position only when it suppresses softclips,
    // and we ask it to suppress none.
    int64_t pos_adj = result.path_position;
    size_t path_len = position_graph.get_path_length(path);
    result.cigar = cigar_to_string(cigar_against_path(chosen, result.path_reverse, pos_adj,
                                                      path_len, /*softclip_suppress=*/0));
    count_edits(chosen, result.matches, result.mismatches,
                result.inserted, result.deleted);
}

/// Return true if path visits the nodes of the mappings [begin, end) of
/// aln_path in order, in either orientation.
///
/// We look only at the path's visits to the first mapping's node, so this is
/// cheap on a graph that lists just the visits of a few indexed paths there.
bool path_carries_mappings(const PathPositionHandleGraph& graph, const path_handle_t& path,
                           const Path& aln_path, size_t begin, size_t end) {
    if (begin >= end) {
        return false;
    }
    auto handle_of = [&](size_t k) {
        const auto& pos = aln_path.mapping(k).position();
        return graph.get_handle(pos.node_id(), pos.is_reverse());
    };
    const handle_t first = handle_of(begin);
    bool carried = false;
    graph.for_each_step_on_handle(first, [&](const step_handle_t& step) {
        if (graph.get_path_handle_of_step(step) != path) {
            return true;
        }
        // A visit in the first mapping's orientation reads the mappings
        // forward along the path. A visit in the other orientation reads them
        // backward, with each node flipped.
        const bool forward = graph.get_handle_of_step(step) == first;
        step_handle_t cursor = step;
        size_t k = begin + 1;
        for (; k < end; ++k) {
            if (forward ? !graph.has_next_step(cursor) : !graph.has_previous_step(cursor)) {
                break;
            }
            cursor = forward ? graph.get_next_step(cursor) : graph.get_previous_step(cursor);
            handle_t visited = graph.get_handle_of_step(cursor);
            if ((forward ? visited : graph.flip(visited)) != handle_of(k)) {
                break;
            }
        }
        carried = (k == end);
        return !carried;
    });
    return carried;
}

} // namespace

HaplotypeSurjectionResult GiraffeEngine::surject_to_indexed_path(
    const Alignment& aln,
    const string& target_path_name,
    size_t begin_mapping, size_t end_mapping) const
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

    // The overlay has only the paths whose positions it indexes.
    if (!position_overlay_ || !position_overlay_->has_path(target_path_name)) {
        if (gbz_->graph.has_path(target_path_name)) {
            result.status = S::PATH_NOT_INDEXED;
        } else {
            result.status = S::UNKNOWN_PATH;
        }
        return result;
    }

    // The overlay lists only its indexed paths' visits on each node, so the
    // check walks just the target path. The overlay keeps the GBZ graph's path
    // handles, so the GBZ's handle works with the overlay's surjector.
    path_handle_t gbz_path_handle = gbz_->graph.get_path_handle(target_path_name);
    if (!path_carries_mappings(*position_overlay_,
                               position_overlay_->get_path_handle(target_path_name),
                               aln.path(), begin_mapping, end_mapping)) {
        result.status = S::INCOMPATIBLE;
        return result;
    }

    unordered_set<path_handle_t> paths;
    paths.insert(gbz_path_handle);
    vector<tuple<string, int64_t, bool>> positions_out;
    vector<Alignment> surjected = overlay_surjector_->surject(aln, paths, positions_out,
                                                              /*allow_negative_scores=*/true,
                                                              /*preserve_deletions=*/false);
    take_best_surjection(surjected, positions_out, *position_overlay_,
                         position_overlay_->get_path_handle(target_path_name), result);
    return result;
}

namespace {

/// Return the distinct PanSN names, sorted, of the haplotypes with these GBWT
/// sequence IDs.
///
/// A haplotype can have many sequence IDs: one per orientation, per visit
/// and per fragment of its path. A sequence ID that GBWT metadata does not
/// name stands for itself, as a number.
std::vector<string> compute_haplotype_names(const gbwtgraph::GBZ& gbz,
                                            const std::vector<gbwt::size_type>& seq_ids) {
    std::set<string> uniq;
    const gbwt::GBWT& index = gbz.index;
    const bool have_md = index.hasMetadata();
    const gbwt::Metadata& md = index.metadata;
    for (gbwt::size_type s : seq_ids) {
        gbwt::size_type path_id = s / 2;
        if (!have_md || path_id >= md.paths()) {
            uniq.insert(std::to_string(s));
            continue;
        }
        gbwt::PathName pn = md.path(path_id);
        string sample = (md.samples() > 0 && pn.sample < md.samples())
                        ? md.sample(pn.sample) : std::to_string(pn.sample);
        string contig = (md.contigs() > 0 && pn.contig < md.contigs())
                        ? md.contig(pn.contig) : std::to_string(pn.contig);
        uniq.insert(sample + "#" + std::to_string(pn.phase) + "#" + contig);
    }
    return std::vector<string>(uniq.begin(), uniq.end());
}

/// Return the PanSN name of a GBWT path. A sample or contig that the metadata
/// does not name appears as its number.
string path_name_for(const gbwt::Metadata& md, gbwt::size_type path_id) {
    gbwt::PathName pn = md.path(path_id);
    string sample = (md.samples() > 0 && pn.sample < md.samples())
                    ? md.sample(pn.sample) : std::to_string(pn.sample);
    string contig = (md.contigs() > 0 && pn.contig < md.contigs())
                    ? md.contig(pn.contig) : std::to_string(pn.contig);
    return sample + "#" + std::to_string(pn.phase) + "#" + contig;
}

/// Choose one haplotype to stand for all the haplotypes with these GBWT
/// sequence IDs.
///
/// Prefer a reference sample or the generic path sample, and then the lowest
/// sequence ID, so that the choice is the same on every run. Return the
/// chosen haplotype's PanSN name, and whether it is a reference. Return an
/// empty name if there are no sequence IDs or no GBWT metadata.
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

} // namespace

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

void GiraffeEngine::load(const GiraffeEngineConfig& config) {
    if (config.gbz_path.empty() || config.minimizer_path.empty() || config.distance_path.empty() || config.zipcode_path.empty()) {
        throw runtime_error("All index paths must be provided: GBZ, minimizer, distance, zipcodes");
    }
    if (config.threads == 0) {
        throw runtime_error("Thread count must be >= 1");
    }
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

    auto assigner = make_unique<HaplotypeAssigner>(*gbz, *distance_index);

    // Index the positions of the surjection target paths once, here. The
    // overlay indexes reference and generic paths, and the target names are
    // passed as extra paths so that haplotype paths are indexed too. A target
    // that is not in the graph fails the load.
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
        set_up_surjector(*surjector, config);
    }

    // map_reads() maps on OpenMP threads, whose count is global to the process.
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

MinimizerMapper& GiraffeEngine::mapper() {
    require_loaded();
    return *mapper_;
}

const gbwtgraph::DefaultMinimizerIndex& GiraffeEngine::minimizer_index() const {
    require_loaded();
    return *minimizer_index_;
}

const gbwtgraph::GBZ& GiraffeEngine::gbz() const {
    require_loaded();
    return *gbz_;
}

vector<string> GiraffeEngine::gaf_header_lines() const {
    require_loaded();
    gbwtgraph::GraphName graph_name = gbz_->graph_name();
    return graph_name.gaf_header_lines();
}

vector<MappedRead> GiraffeEngine::map_reads(const vector<GiraffeServerRequest>& reads) {
    require_loaded();
    vector<MappedRead> results(reads.size());
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
                    // Name the haplotypes of each segment and choose its
                    // representative.
                    for (auto& seg : hap.segments) {
                        seg.haplotype_names = compute_haplotype_names(*gbz_, seg.seq_ids);
                        seg.representative_name = choose_representative(*gbz_, seg.seq_ids).first;
                    }
                    // Report the longest segment's haplotypes for the whole
                    // alignment, so the tags name the haplotypes that cover
                    // most of the read.
                    if (hap.best_segment < hap.segments.size()) {
                        const auto& best = hap.segments[hap.best_segment];
                        hap.haplotype_names = best.haplotype_names;
                        auto rep = choose_representative(*gbz_, best.seq_ids);
                        hap.representative_name = rep.first;
                        hap.representative_is_reference = rep.second;
                    }
                    hap_results.emplace_back(std::move(hap));
                }
            }

            // When the read names a target, every alignment gets a surjection
            // status, so the client can tell why a surjection is missing.
            vector<HaplotypeSurjectionResult> surj_results;
            if (!req.surjection_target.empty()) {
                surj_results.resize(mapped.size());
                if (!overlay_surjector_) {
                    for (auto& r : surj_results) {
                        r.status = HaplotypeSurjectionResult::Status::PATH_NOT_INDEXED;
                    }
                } else {
                    for (size_t j = 0; j < mapped.size(); ++j) {
                        // The target must carry the mappings that the reported
                        // haplotypes carry: the best segment, or the whole
                        // path without a haplotype assignment.
                        size_t begin = 0;
                        size_t end = mapped[j].path().mapping_size();
                        if (!hap_results.empty()) {
                            const auto& hap = hap_results[j];
                            if (hap.best_segment < hap.segments.size()) {
                                begin = hap.segments[hap.best_segment].start_mapping;
                                end = hap.segments[hap.best_segment].end_mapping;
                            } else {
                                // No haplotype visits any node of the path.
                                end = begin;
                            }
                        }
                        surj_results[j] = surject_to_indexed_path(
                            mapped[j], req.surjection_target, begin, end);
                    }
                }
            }

            auto& alignments = results[i].alignments;
            alignments.reserve(mapped.size());
            for (size_t j = 0; j < mapped.size(); ++j) {
                MappedAlignment out;
                out.alignment = std::move(mapped[j]);
                if (!hap_results.empty()) {
                    out.haplotypes = std::move(hap_results[j]);
                }
                if (j < surj_results.size()) {
                    out.surjection = std::move(surj_results[j]);
                }
                alignments.emplace_back(std::move(out));
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

/// Return the GBWTGraph step handle for a GBWT node and an offset in that
/// node's GBWT record.
handlegraph::step_handle_t make_step_handle(uint64_t node, uint64_t offset) {
    handlegraph::step_handle_t step;
    handlegraph::as_integers(step)[0] = node;
    handlegraph::as_integers(step)[1] = offset;
    return step;
}

/// Return the length in bases of a GBWT path, by walking it.
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

    gafkluge::GafRecord gaf_record;
    gafkluge::parse_gaf_record(gaf_line, gaf_record);
    Alignment aln;
    vg::io::gaf_to_alignment(gbz_->graph, gaf_record, aln);
    if (!read_name.empty()) {
        aln.set_name(read_name);
    } else if (aln.name().empty()) {
        aln.set_name(gaf_record.query_name);
    }

    // Find the target path. A sample and haplotype such as HG002#2 is not a
    // path name in the graph. Then the target is the path of the first
    // anchor's step, since all the anchors lie on the target path.
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

    // Walk the target path for its length if the caller did not give it.
    size_t path_length = target_path_length;
    if (path_length == 0) {
        gbwt::size_type gbwt_path_id = gbz_->graph.handle_to_path(target_path);
        path_length = walk_target_path_length(*gbz_, gbwt_path_id);
    }

    // The anchor-backed position graph needs only the steps and their path
    // offsets.
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

    // Surject onto the target path. The anchor-backed position graph answers
    // the surjector's position queries.
    AnchorBackedPositionGraph anchor_graph(&gbz_->graph, vg_anchors,
                                           target_path, path_length);
    Surjector anchor_surjector(&anchor_graph);
    set_up_surjector(anchor_surjector, config_);

    unordered_set<path_handle_t> paths;
    paths.insert(target_path);
    vector<tuple<string, int64_t, bool>> positions_out;
    vector<Alignment> surjected;
    try {
        surjected = anchor_surjector.surject(aln, paths, positions_out,
                                             /*allow_negative_scores=*/true,
                                             /*preserve_deletions=*/false);
    } catch (const std::exception& e) {
        // The error fails only this surjection, and goes back to the client.
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
        return {ss.str()};
    }
    ss << "\tan:i:" << anchors.size()
       << "\tap:Z:" << resolved_target;
    return {ss.str()};
}

} // namespace vg
