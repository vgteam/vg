#ifndef VG_GIRAFFE_ENGINE_HPP_INCLUDED
#define VG_GIRAFFE_ENGINE_HPP_INCLUDED

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include "minimizer_mapper.hpp"
#include "snarl_distance_index.hpp"
#include "zip_code.hpp"

namespace vg {

struct GiraffeEnginePaths {
    std::string gbz_path;
    std::string minimizer_path;
    std::string distance_path;
    std::string zipcode_path;
};

struct GiraffeEngineConfig {
    size_t threads = 1;
    size_t max_multimaps = 1;
    bool preload_distance_index = true;
};

struct GiraffeFastqRead {
    std::string name;
    std::string sequence;
    std::string quality;
};

/// One oriented node traversal on the GBWT graph (from an Alignment path).
struct GiraffeGraphStep {
    int64_t node_id = 0;
    bool is_reverse = false;
    int32_t from_length = 0;
};

/// One alignment hypothesis for a read: GAF line plus the underlying graph path.
struct GiraffeAlignmentRecord {
    std::string gaf_line;
    std::vector<GiraffeGraphStep> graph_path;
};

/// All alignments returned for a single read (primary + multimaps).
struct GiraffeReadMappings {
    std::vector<GiraffeAlignmentRecord> alignments;
};

class GiraffeEngine {
public:
    GiraffeEngine() = default;

    void load(const GiraffeEnginePaths& paths, const GiraffeEngineConfig& config);

    bool is_loaded() const;

    const GiraffeEngineConfig& config() const;

    std::vector<std::string> gaf_header_lines() const;

    // Returns one vector per input read, preserving order.
    // Each inner vector contains one or more GAF lines (multimaps).
    std::vector<std::vector<std::string>> map_reads(const std::vector<GiraffeFastqRead>& reads);

    /// Same as map_reads but includes per-alignment graph path (node, strand, from_length per Mapping).
    std::vector<GiraffeReadMappings> map_reads_detailed(const std::vector<GiraffeFastqRead>& reads);

private:
    void require_loaded() const;

    GiraffeEngineConfig config_{};
    std::unique_ptr<gbwtgraph::GBZ> gbz_;
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index_;
    std::unique_ptr<SnarlDistanceIndex> distance_index_;
    std::unique_ptr<ZipCodeCollection> zipcodes_;
    std::unique_ptr<MinimizerMapper> mapper_;
};

} // namespace vg

#endif // VG_GIRAFFE_ENGINE_HPP_INCLUDED
