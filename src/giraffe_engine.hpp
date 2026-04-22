#ifndef VG_GIRAFFE_ENGINE_HPP_INCLUDED
#define VG_GIRAFFE_ENGINE_HPP_INCLUDED

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include "haplotype_assignment.hpp"
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
};

} // namespace vg

#endif // VG_GIRAFFE_ENGINE_HPP_INCLUDED
