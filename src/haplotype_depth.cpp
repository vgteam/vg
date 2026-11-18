#include "haplotype_depth.hpp"

#include "feature_lift.hpp"
#include "gref.hpp"

#include <gbwtgraph/utils.h>
#include <omp.h>

#include <algorithm>
#include <map>
#include <unordered_set>

namespace vg {

using gbwtgraph::GBWTGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;
using handlegraph::PathSense;

//------------------------------------------------------------------------------

PathRegions::PathRegions(std::istream& bed) {
    for (BedFeature& feature : read_bed_features(bed)) {
        regions[feature.sequence].emplace_back(feature.start, feature.end);
    }
    for (auto& [name, list] : regions) {
        std::sort(list.begin(), list.end());
        size_t kept = 0;
        for (const Region& region : list) {
            if (kept > 0 && region.first <= list[kept - 1].second) {
                list[kept - 1].second = std::max(list[kept - 1].second, region.second);
            } else {
                list[kept++] = region;
            }
        }
        list.resize(kept);
    }
}

const std::vector<PathRegions::Region>* PathRegions::regions_of(const std::string& name) const {
    auto found = regions.find(name);
    return found == regions.end() ? nullptr : &found->second;
}

std::vector<std::string> PathRegions::path_names() const {
    std::vector<std::string> names;
    for (const auto& [name, list] : regions) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

//------------------------------------------------------------------------------

HaplotypeDepth::HaplotypeDepth(const gbwtgraph::GBZ& gbz, const Options& options) {
    const GBWTGraph& graph = gbz.graph;
    if (graph.get_node_count() == 0) {
        return;
    }
    min_id = graph.min_node_id();
    whole.assign(graph.max_node_id() - min_id + 1, 0);

    // Group the counted paths by haplotype. One haplotype is usually several
    // paths: the contigs of an assembly, or their fragments in a clipped graph.
    // Counting by haplotype makes a base that two of its paths read count once.
    HaplotypePaths by_haplotype = counted_paths(gbz, options.count_references);
    haplotypes = by_haplotype.size();
    if (options.regions != nullptr) {
        match_regions(by_haplotype, *options.regions);
    }

    std::vector<const std::vector<CountedPath>*> groups;
    for (const auto& [haplotype, paths] : by_haplotype) {
        groups.push_back(&paths);
    }
    // One thread counts each haplotype, so that it can merge the stretches that
    // the haplotype's paths read before adding them.
    std::vector<NodeStretch> cut_stretches;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < groups.size(); i++) {
        std::vector<NodeStretch> stretches;
        for (const CountedPath& path : *groups[i]) {
            if (options.regions == nullptr) {
                add_path_stretches(gbz, path, stretches);
            } else if (const std::vector<PathRegions::Region>* regions = options.regions->regions_of(path.name)) {
                add_region_stretches(gbz, path, *regions, stretches);
            }
        }
        std::vector<NodeStretch> cut_here;
        add_stretches(std::move(stretches), graph, !options.count_visits, cut_here);
        #pragma omp critical (haplotype_depth_cuts)
        cut_stretches.insert(cut_stretches.end(), cut_here.begin(), cut_here.end());
    }
    add_cuts(std::move(cut_stretches));
}

HaplotypeDepth::HaplotypePaths HaplotypeDepth::counted_paths(const gbwtgraph::GBZ& gbz, bool count_references) {
    HaplotypePaths by_haplotype;
    gbwtgraph::sample_name_set reference_samples = gbwtgraph::parse_reference_samples_tag(gbz.index);
    for_each_path_start(gbz.index, [&](gbwt::size_type path_id, const std::string& name, uint64_t start) {
        if (GrefCover::is_gref_derived(name)) {
            return;
        }
        PathSense sense = gbwtgraph::get_path_sense(gbz.index, path_id, reference_samples);
        if (sense == PathSense::GENERIC || (sense == PathSense::REFERENCE && !count_references)) {
            return;
        }
        std::string haplotype = gbwtgraph::get_path_sample_name(gbz.index, path_id, sense) + "#"
                              + std::to_string(gbwtgraph::get_path_haplotype(gbz.index, path_id, sense));
        by_haplotype[haplotype].emplace_back(path_id, name, start);
    });
    return by_haplotype;
}

void HaplotypeDepth::match_regions(const HaplotypePaths& by_haplotype, const PathRegions& regions) {
    std::unordered_set<std::string> counted_names;
    for (const auto& [haplotype, paths] : by_haplotype) {
        if (std::none_of(paths.begin(), paths.end(), [&](const CountedPath& path) {
                return regions.regions_of(path.name) != nullptr;
            })) {
            without_regions++;
        }
        for (const CountedPath& path : paths) {
            counted_names.insert(path.name);
        }
    }
    for (const std::string& name : regions.path_names()) {
        if (!counted_names.contains(name)) {
            unmatched_regions.push_back(name);
        }
    }
}

void HaplotypeDepth::add_path_stretches(const gbwtgraph::GBZ& gbz, const CountedPath& path,
                                        std::vector<NodeStretch>& stretches) {
    for (gbwt::node_type node : gbz.index.extract(gbwt::Path::encode(path.id, false))) {
        nid_t id = gbwt::Node::id(node);
        auto length = static_cast<uint32_t>(gbz.graph.get_length(GBWTGraph::node_to_handle(node)));
        stretches.push_back({ id, 0, length });
    }
}

void HaplotypeDepth::add_region_stretches(const gbwtgraph::GBZ& gbz, const CountedPath& path,
                                          const std::vector<PathRegions::Region>& regions,
                                          std::vector<NodeStretch>& stretches) {
    uint64_t coord = path.start;
    size_t next_region = 0;
    for (gbwt::node_type node : gbz.index.extract(gbwt::Path::encode(path.id, false))) {
        nid_t id = gbwt::Node::id(node);
        auto length = static_cast<uint32_t>(gbz.graph.get_length(GBWTGraph::node_to_handle(node)));
        // Skip the regions that end before this step. The regions are
        // disjoint and sorted, and later steps lie further along the path, so
        // no later step needs them.
        while (next_region < regions.size() && regions[next_region].second <= coord) {
            next_region++;
        }
        for (size_t r = next_region; r < regions.size() && regions[r].first < coord + length; r++) {
            auto from = static_cast<uint32_t>(std::max(coord, regions[r].first) - coord);
            auto to = static_cast<uint32_t>(std::min(coord + length, regions[r].second) - coord);
            // A reverse step reads the node from its end.
            if (gbwt::Node::is_reverse(node)) {
                std::tie(from, to) = std::make_pair(length - to, length - from);
            }
            stretches.push_back({ id, from, to });
        }
        coord += length;
    }
}

void HaplotypeDepth::add_cuts(std::vector<NodeStretch> cut_stretches) {
    std::sort(cut_stretches.begin(), cut_stretches.end(), [](const NodeStretch& a, const NodeStretch& b) {
        return a.node < b.node;
    });
    for (size_t i = 0; i < cut_stretches.size();) {
        std::vector<std::pair<uint32_t, int32_t>> changes;
        size_t j = i;
        for (; j < cut_stretches.size() && cut_stretches[j].node == cut_stretches[i].node; j++) {
            changes.emplace_back(cut_stretches[j].from, 1);
            changes.emplace_back(cut_stretches[j].to, -1);
        }
        std::sort(changes.begin(), changes.end());
        size_t kept = 0;
        for (const auto& change : changes) {
            if (kept > 0 && changes[kept - 1].first == change.first) {
                changes[kept - 1].second += change.second;
            } else {
                changes[kept++] = change;
            }
        }
        changes.resize(kept);
        cuts.try_emplace(cut_stretches[i].node, std::move(changes));
        i = j;
    }
}

void HaplotypeDepth::add_stretches(std::vector<NodeStretch> stretches, const GBWTGraph& graph, bool merge,
                                   std::vector<NodeStretch>& cut_stretches) {
    if (merge) {
        // A haplotype that reads a base more than once counts once there.
        std::sort(stretches.begin(), stretches.end(), [](const NodeStretch& a, const NodeStretch& b) {
            return std::tie(a.node, a.from, a.to) < std::tie(b.node, b.from, b.to);
        });
        size_t kept = 0;
        for (const NodeStretch& stretch : stretches) {
            NodeStretch* last = kept > 0 ? &stretches[kept - 1] : nullptr;
            if (last != nullptr && last->node == stretch.node && stretch.from <= last->to) {
                last->to = std::max(last->to, stretch.to);
            } else {
                stretches[kept++] = stretch;
            }
        }
        stretches.resize(kept);
    }
    for (const NodeStretch& stretch : stretches) {
        if (stretch.from == 0 && stretch.to == graph.get_length(graph.get_handle(stretch.node))) {
            #pragma omp atomic
            whole[stretch.node - min_id]++;
        } else {
            cut_stretches.push_back(stretch);
        }
    }
}

uint64_t HaplotypeDepth::bases_read(nid_t id, uint32_t length) const {
    if (whole.empty() || id < min_id || size_t(id - min_id) >= whole.size()) {
        return 0;
    }
    uint64_t total = uint64_t(whole[id - min_id]) * length;
    if (auto found = cuts.find(id); found != cuts.end()) {
        // The count of the stretch after each change holds until the next one.
        int64_t depth = 0;
        uint32_t at = 0;
        for (const auto& [offset, change] : found->second) {
            total += uint64_t(depth) * (offset - at);
            depth += change;
            at = offset;
        }
    }
    return total;
}

namespace {

/**
 * A builder of the runs of equal depth along a path.
 *
 * It takes stretches of bases in order along the path. It joins a stretch to
 * the run before it if they touch and have the same depth, and reports each
 * run once it ends.
 */
class RunBuilder {
public:
    /// Report each run with visit(start, end, depth).
    explicit RunBuilder(const std::function<void(uint64_t, uint64_t, uint32_t)>& visit) : visit(visit) {
    }

    /// Add bases [start, end) of the full path, at the given depth.
    void extend(uint64_t start, uint64_t end, uint32_t depth) {
        if (start == end) {
            return;
        }
        if (open && depth == run_depth && start == run_end) {
            run_end = end;
            return;
        }
        finish();
        open = true;
        run_start = start;
        run_end = end;
        run_depth = depth;
    }

    /// Report the run in progress, if any.
    void finish() {
        if (open) {
            visit(run_start, run_end, run_depth);
            open = false;
        }
    }

private:
    const std::function<void(uint64_t, uint64_t, uint32_t)>& visit;
    /// True if a run is in progress.
    bool open = false;
    uint64_t run_start = 0;
    uint64_t run_end = 0;
    uint32_t run_depth = 0;
};

/// A stretch [from, to) of a node, in its forward orientation, at one depth.
struct DepthPiece {
    uint32_t from;
    uint32_t to;
    uint32_t depth;
};

/**
 * Return the stretches of a cut node of the given length between the changes
 * of its depth, in its forward orientation.
 *
 * The depth starts at base_depth and changes at each offset in changes.
 */
std::vector<DepthPiece> depth_pieces(const std::vector<std::pair<uint32_t, int32_t>>& changes, uint32_t base_depth,
                                     uint32_t length) {
    std::vector<DepthPiece> pieces;
    int64_t depth = base_depth;
    uint32_t at = 0;
    for (const auto& [offset, change] : changes) {
        if (offset > at) {
            pieces.push_back({ at, offset, uint32_t(depth) });
        }
        depth += change;
        at = offset;
    }
    if (at < length) {
        pieces.push_back({ at, length, uint32_t(depth) });
    }
    return pieces;
}

}

void HaplotypeDepth::for_each_run(const PathSteps& path, uint64_t offset,
                                  const std::function<void(uint64_t, uint64_t, uint32_t)>& visit) const {
    RunBuilder runs(visit);
    const std::vector<PathSteps::Step>& steps = path.steps();
    for (size_t i = 0; i < steps.size(); i++) {
        uint64_t start = offset + steps[i].start;
        auto length = static_cast<uint32_t>((i + 1 < steps.size() ? steps[i + 1].start : path.length()) - steps[i].start);
        gbwt::node_type node = GBWTGraph::handle_to_node(steps[i].handle);
        nid_t id = gbwt::Node::id(node);
        uint32_t base_depth = whole.empty() ? 0 : whole[id - min_id];
        auto found = cuts.find(id);
        if (found == cuts.end()) {
            runs.extend(start, start + length, base_depth);
            continue;
        }
        std::vector<DepthPiece> pieces = depth_pieces(found->second, base_depth, length);
        if (gbwt::Node::is_reverse(node)) {
            for (auto piece = pieces.rbegin(); piece != pieces.rend(); ++piece) {
                runs.extend(start + length - piece->to, start + length - piece->from, piece->depth);
            }
        } else {
            for (const DepthPiece& piece : pieces) {
                runs.extend(start + piece.from, start + piece.to, piece.depth);
            }
        }
    }
    runs.finish();
}

}
