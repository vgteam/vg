#include "giraffe_engine.hpp"

#include <fstream>
#include <atomic>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include <omp.h>

#include <vg/io/alignment_io.hpp>
#include <vg/io/vpkg.hpp>

#include "alignment.hpp"
#include "gbwtgraph_helper.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

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

    // Keep the same behavior as vg giraffe: global OMP thread count controls mapping parallelism.
    omp_set_num_threads(static_cast<int>(config.threads));

    config_ = config;
    gbz_ = move(gbz);
    minimizer_index_ = move(minimizer_index);
    distance_index_ = move(distance_index);
    zipcodes_ = move(zipcodes);
    mapper_ = move(mapper);
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
            auto& lines = results[i];
            lines.reserve(mapped.size());
            for (auto& m : mapped) {
                auto gaf = io::alignment_to_gaf(gbz_->graph, m);
                stringstream ss;
                ss << gaf;
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
