/**
 * \file feature_lift.cpp
 * Lifting features, see feature_lift.hpp.
 */

#include "feature_lift.hpp"

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <istream>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

namespace {

/**
 * Parse a BED line into a feature.
 *
 * @throws std::runtime_error if the line has fewer than 3 columns, bad
 * coordinates, no bases, or a strand other than +, - and ".".
 */
BedFeature parse_bed_line(const std::string& line, size_t line_number) {
    auto fail = [&](const std::string& problem) {
        throw std::runtime_error("BED line " + std::to_string(line_number) + ": " + problem + ": " + line);
    };
    std::vector<std::string> columns;
    std::stringstream stream(line);
    std::string column;
    while (std::getline(stream, column, '\t')) {
        columns.push_back(column);
    }
    if (columns.size() < 3) {
        fail("fewer than 3 columns");
    }
    BedFeature feature;
    feature.sequence = columns[0];
    try {
        size_t used_start = 0;
        size_t used_end = 0;
        feature.start = std::stoull(columns[1], &used_start);
        feature.end = std::stoull(columns[2], &used_end);
        if (used_start != columns[1].size() || used_end != columns[2].size()) {
            fail("bad coordinates");
        }
    } catch (const std::invalid_argument&) {
        fail("bad coordinates");
    } catch (const std::out_of_range&) {
        fail("bad coordinates");
    }
    if (feature.start >= feature.end) {
        fail("empty feature");
    }
    feature.fields.assign(columns.begin() + 3, columns.end());
    if (feature.fields.size() > BED_STRAND_FIELD) {
        const std::string& strand = feature.fields[BED_STRAND_FIELD];
        if (strand != "+" && strand != "-" && strand != ".") {
            fail("strand is not +, - or .");
        }
    }
    return feature;
}

}

std::vector<BedFeature> read_bed_features(std::istream& in) {
    std::vector<BedFeature> features;
    std::string line;
    size_t line_number = 0;
    while (std::getline(in, line)) {
        line_number++;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty() || line[0] == '#' || line.starts_with("track") || line.starts_with("browser")) {
            continue;
        }
        features.push_back(parse_bed_line(line, line_number));
    }
    return features;
}

bool GraphFeature::is_stranded() const {
    return fields.size() > BED_STRAND_FIELD && fields[BED_STRAND_FIELD] != ".";
}

bool PathCoordinate::operator<(const PathCoordinate& other) const {
    return std::tie(rank, start, end, is_reverse) < std::tie(other.rank, other.start, other.end, other.is_reverse);
}

std::string PathCoordinate::text() const {
    return path + ':' + std::to_string(start + 1) + '-' + std::to_string(end) + ':' + (is_reverse ? '-' : '+');
}

GraphFeature make_graph_feature(const GBWTGraph& graph, const GrefIndex* gref,
                                WalkInterval interval, std::vector<std::string> fields,
                                std::optional<handle_t> before, std::optional<handle_t> after) {
    GraphFeature result;
    result.interval = std::move(interval);
    result.fields = std::move(fields);
    // Read the interval the other way, with the nodes on either side.
    auto reverse = [&]() {
        result.interval = result.interval.reverse(graph);
        std::optional<handle_t> old_before = before;
        before = after ? std::optional<handle_t>(graph.flip(*after)) : std::nullopt;
        after = old_before ? std::optional<handle_t>(graph.flip(*old_before)) : std::nullopt;
    };
    bool stranded = result.is_stranded();
    if (stranded && result.fields[BED_STRAND_FIELD] == "-") {
        reverse();
        result.fields[BED_STRAND_FIELD] = "+";
    }

    if (gref != nullptr) {
        result.placement = gref->place(result.interval, before, after);
    }
    if (!stranded) {
        const GrefPlacement& placement = result.placement;
        bool flip = placement.reverse_bases > placement.forward_bases ||
                    (placement.reverse_bases == placement.forward_bases &&
                     result.interval.reverse(graph).compare(result.interval) < 0);
        if (flip) {
            reverse();
            if (gref != nullptr) {
                result.placement = gref->place(result.interval, before, after);
            }
        }
    }

    // Gref paths come first in the priority order, and their ranks are their
    // places in the gref.
    if (gref != nullptr && result.placement.along_one_path) {
        uint32_t path = result.placement.path;
        uint64_t start = gref->path_shift(path) + result.placement.coord;
        result.coordinate = PathCoordinate { gref->path_name(path), path, start, start + result.interval.length,
                                             result.placement.reverse_bases > 0 };
    }

    // The GBWT holds every path in both orientations, so this counts visits in
    // either orientation. Leave out the gref path that the interval lies
    // along, if any, because a gref path is a copy of part of another path.
    result.support = graph.bd_find(result.interval.walk).size();
    if (result.placement.along_one_path) {
        result.support--;
    }
    return result;
}

//------------------------------------------------------------------------------

Cover::Cover(const GBWTGraph& graph, const GrefIndex* gref) : graph(graph), gref(gref) {
}

void Cover::add_tigs(std::istream& gaf, const std::function<bool(const std::string&)>& wanted) {
    std::string line;
    std::string name;
    while (std::getline(gaf, line)) {
        // GAF header lines start with @.
        if (line.empty() || line[0] == '@' || line[0] == '#') {
            continue;
        }
        name = line.substr(0, line.find('\t'));
        if (!wanted(name)) {
            continue;
        }
        WalkInterval interval = parse_gaf_columns(line, graph, name);
        if (!tigs.try_emplace(name, std::move(interval)).second) {
            throw std::runtime_error("tig " + name + " appears twice in the GAF");
        }
    }
}

bool Cover::is_gref_path(const std::string& name) const {
    return gref != nullptr && gref->find_path(name) != GrefIndex::NO_PATH;
}

std::optional<size_t> Cover::sequence_length(const std::string& name) const {
    if (is_gref_path(name)) {
        return gref->path_length(gref->find_path(name));
    }
    if (auto found = tigs.find(name); found != tigs.end()) {
        return found->second.length;
    }
    return std::nullopt;
}

WalkInterval Cover::interval(const std::string& name, size_t start, size_t end) const {
    if (is_gref_path(name)) {
        return gref->interval(gref->find_path(name), start, end);
    }
    return tigs.at(name).subinterval(graph, start, end);
}

//------------------------------------------------------------------------------

PathPriority::PathPriority(const gbwtgraph::GBZ& gbz, const std::vector<std::string>& prefixes,
                           const std::string& gref_prefix, size_t first_rank) {
    const gbwt::GBWT& index = gbz.index;
    gbwtgraph::sample_name_set reference_samples = gbwtgraph::parse_reference_samples_tag(index);
    // Sort keys: the group (a prefix, reference paths, the rest), then the name.
    std::vector<std::pair<size_t, std::string>> keys;
    std::unordered_map<std::string, size_t> groups;
    for_each_path_start(index, [&](gbwt::size_type path_id, const std::string& name, uint64_t) {
        if (name.compare(0, gref_prefix.size(), gref_prefix) == 0 || groups.contains(name)) {
            return;
        }
        size_t group = prefixes.size() + 1;
        for (size_t i = 0; i < prefixes.size(); i++) {
            if (name.compare(0, prefixes[i].size(), prefixes[i]) == 0) {
                group = i;
                break;
            }
        }
        if (group > prefixes.size() &&
            gbwtgraph::get_path_sense(index, path_id, reference_samples) == handlegraph::PathSense::REFERENCE) {
            group = prefixes.size();
        }
        groups[name] = group;
        keys.emplace_back(group, name);
    });
    std::sort(keys.begin(), keys.end());
    for (size_t i = 0; i < keys.size(); i++) {
        ranks[keys[i].second] = first_rank + i;
    }
}

size_t PathPriority::rank(const std::string& name) const {
    return ranks.at(name);
}

//------------------------------------------------------------------------------

PathDirectory::PathDirectory(const gbwtgraph::GBZ& gbz, const std::function<bool(const std::string&)>& wanted) {
    for_each_path_start(gbz.index, [&](gbwt::size_type path_id, const std::string& name, uint64_t start) {
        if (wanted(name)) {
            fragments[name].emplace_back(start, PathSteps(gbz.graph, path_id));
        }
    });
    for (auto& [name, path_fragments] : fragments) {
        std::sort(path_fragments.begin(), path_fragments.end(), [](const Fragment& a, const Fragment& b) {
            return a.start < b.start;
        });
    }
}

bool PathDirectory::has_path(const std::string& name) const {
    return fragments.contains(name);
}

std::optional<WalkInterval> PathDirectory::interval(const std::string& name, uint64_t start, uint64_t end) const {
    auto found = fragments.find(name);
    if (found == fragments.end()) {
        return std::nullopt;
    }
    for (const Fragment& fragment : found->second) {
        if (fragment.start <= start && end <= fragment.start + fragment.steps.length()) {
            return fragment.steps.interval(start - fragment.start, end - fragment.start);
        }
    }
    return std::nullopt;
}

//------------------------------------------------------------------------------

void sort_and_merge_features(std::vector<GraphFeature>& features) {
    std::sort(features.begin(), features.end(), [](const GraphFeature& a, const GraphFeature& b) {
        if (a.placement.path != b.placement.path) {
            return a.placement.path < b.placement.path;
        }
        if (a.placement.coord != b.placement.coord) {
            return a.placement.coord < b.placement.coord;
        }
        if (int order = a.interval.compare(b.interval); order != 0) {
            return order < 0;
        }
        return a.fields < b.fields;
    });
    auto end = std::unique(features.begin(), features.end(), [](const GraphFeature& a, const GraphFeature& b) {
        return a.interval.compare(b.interval) == 0 && a.fields == b.fields;
    });
    features.erase(end, features.end());
}

//------------------------------------------------------------------------------

FeatureMatcher::FeatureMatcher(const GBWTGraph& graph, const std::vector<GraphFeature>& features) :
    graph(graph), features(features) {
    if (graph.get_node_count() == 0) {
        return;
    }
    min_id = graph.min_node_id();
    starts.assign(graph.max_node_id() - min_id + 2, 0);
    // The walk read in reverse starts on its last node.
    auto first_node = [&](const GraphFeature& feature, bool is_reverse) {
        const std::vector<handle_t>& walk = feature.interval.walk;
        return graph.get_id(is_reverse ? walk.back() : walk.front()) - min_id;
    };
    for (const GraphFeature& feature : features) {
        for (bool is_reverse : { false, true }) {
            starts[first_node(feature, is_reverse) + 1]++;
        }
    }
    for (size_t i = 1; i < starts.size(); i++) {
        starts[i] += starts[i - 1];
    }
    entries.resize(starts.back());
    std::vector<size_t> filled(starts.begin(), starts.end() - 1);
    for (size_t i = 0; i < features.size(); i++) {
        for (bool is_reverse : { false, true }) {
            entries[filled[first_node(features[i], is_reverse)]++] = (i << 1) | size_t(is_reverse);
        }
    }
}

void FeatureMatcher::for_each_match(const gbwt::vector_type& path,
                                    const std::function<void(size_t, size_t, bool)>& visit) const {
    if (starts.empty()) {
        return;
    }
    for (size_t step = 0; step < path.size(); step++) {
        handle_t handle = GBWTGraph::node_to_handle(path[step]);
        size_t node = graph.get_id(handle) - min_id;
        for (size_t i = starts[node]; i < starts[node + 1]; i++) {
            size_t feature = entries[i] >> 1;
            bool is_reverse = entries[i] & 1;
            if (visits_walk(path, step, features[feature].interval.walk, is_reverse)) {
                visit(feature, step, is_reverse);
            }
        }
    }
}

bool FeatureMatcher::visits_walk(const gbwt::vector_type& path, size_t step, const std::vector<handle_t>& walk,
                                 bool is_reverse) const {
    if (step + walk.size() > path.size()) {
        return false;
    }
    for (size_t j = 0; j < walk.size(); j++) {
        handle_t expected = is_reverse ? graph.flip(walk[walk.size() - 1 - j]) : walk[j];
        if (GBWTGraph::node_to_handle(path[step + j]) != expected) {
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------

FeatureNodeIndex::FeatureNodeIndex(const GBWTGraph& graph, const std::vector<GraphFeature>& features) :
    graph(graph) {
    if (graph.get_node_count() == 0) {
        return;
    }
    min_id = graph.min_node_id();
    // Call visit(node, entry) for each node of the feature's walk that holds feature bases.
    auto for_each_entry = [&](size_t feature, const std::function<void(size_t, const Entry&)>& visit) {
        const WalkInterval& interval = features[feature].interval;
        size_t node_start = 0;
        for (const handle_t& handle : interval.walk) {
            size_t node_end = node_start + graph.get_length(handle);
            // Walk coordinates of the feature bases on this node.
            size_t low = std::max(node_start, interval.offset);
            if (size_t high = std::min(node_end, interval.offset + interval.length); low < high) {
                visit(graph.get_id(handle) - min_id, { feature, low - interval.offset, high - interval.offset,
                                                       low - node_start, high - node_start, graph.get_is_reverse(handle) });
            }
            node_start = node_end;
        }
    };
    starts.assign(graph.max_node_id() - min_id + 2, 0);
    for (size_t i = 0; i < features.size(); i++) {
        for_each_entry(i, [&](size_t node, const Entry&) {
            starts[node + 1]++;
        });
    }
    for (size_t i = 1; i < starts.size(); i++) {
        starts[i] += starts[i - 1];
    }
    entries.resize(starts.back());
    std::vector<size_t> filled(starts.begin(), starts.end() - 1);
    for (size_t i = 0; i < features.size(); i++) {
        for_each_entry(i, [&](size_t node, const Entry& entry) {
            entries[filled[node]++] = entry;
        });
    }
}

void FeatureNodeIndex::for_each_anchor(const gbwt::vector_type& path, const std::vector<uint64_t>& coords,
                                       const std::function<void(const Anchor&)>& visit) const {
    if (starts.empty()) {
        return;
    }
    for (size_t step = 0; step < path.size(); step++) {
        handle_t handle = GBWTGraph::node_to_handle(path[step]);
        size_t node = graph.get_id(handle) - min_id;
        size_t node_length = graph.get_length(handle);
        for (size_t i = starts[node]; i < starts[node + 1]; i++) {
            const Entry& entry = entries[i];
            // A path that reads the node the other way finds the bases at the
            // mirrored offsets.
            bool same = (graph.get_is_reverse(handle) == entry.is_reverse);
            uint64_t start = coords[step] + (same ? entry.node_from : node_length - entry.node_to);
            visit({ entry.feature, entry.from, entry.to, start, !same });
        }
    }
}

namespace {

using Anchor = FeatureNodeIndex::Anchor;

/// The index of the anchor before the first anchor of a chain.
constexpr size_t NO_ANCHOR = std::numeric_limits<size_t>::max();

/// Return the path coordinate after the last base of the anchor.
uint64_t end_of(const Anchor& anchor) {
    return anchor.start + (anchor.to - anchor.from);
}

/**
 * Return true if a chain on a feature of the given length can go from anchor
 * before to anchor after.
 *
 * The path must read both anchors the same way. On the path, after must start
 * no more than length bases past the end of before. In the feature, after must
 * come later in the order the path reads it, which is decreasing for a path
 * that reads the feature in reverse.
 */
bool can_follow(const Anchor& before, const Anchor& after, size_t length) {
    return before.is_reverse == after.is_reverse && end_of(before) <= after.start &&
           after.start - end_of(before) <= length &&
           (after.is_reverse ? after.to <= before.from : after.from >= before.to);
}

/// Return true if anchor after starts where anchor before ends, both on the
/// path and in the feature as the path reads it.
bool abuts(const Anchor& before, const Anchor& after) {
    return after.start == end_of(before) && (after.is_reverse ? after.to == before.from : after.from == before.to);
}

/// The best chain that ends at an anchor.
struct ChainEnd {
    /// The number of feature bases in the chain.
    size_t score = 0;
    /// The index of the chain's previous anchor, or NO_ANCHOR.
    size_t previous = NO_ANCHOR;
};

/**
 * Find the chain ending at each anchor that holds the most feature bases.
 *
 * The anchors are those of one feature, of the given length, sorted by their
 * start on the path.
 */
std::vector<ChainEnd> best_chains(std::span<const Anchor> anchors, size_t length) {
    std::vector<ChainEnd> ends(anchors.size());
    for (size_t i = 0; i < anchors.size(); i++) {
        const Anchor& anchor = anchors[i];
        size_t bases = anchor.to - anchor.from;
        ends[i].score = bases;
        for (size_t j = i; j-- > 0;) {
            const Anchor& before = anchors[j];
            if (before.start + 2 * length < anchor.start) {
                // This anchor and all earlier ones end more than the
                // feature's length back.
                break;
            }
            if (can_follow(before, anchor, length) && ends[j].score + bases > ends[i].score) {
                ends[i] = { ends[j].score + bases, j };
            }
        }
    }
    return ends;
}

/**
 * Return the place on the path of a chain of anchors.
 *
 * The chain holds indices into the anchors, in chain order. Anchors that abut
 * share a block.
 */
PathOccurrence chain_place(std::span<const Anchor> anchors, const std::vector<size_t>& chain) {
    const Anchor& first = anchors[chain.front()];
    PathOccurrence place { first.start, first.start, first.feature, first.is_reverse, 0, {} };
    const Anchor* before = nullptr;
    for (size_t i : chain) {
        const Anchor& anchor = anchors[i];
        if (before != nullptr && abuts(*before, anchor)) {
            place.blocks.back().second = end_of(anchor);
        } else {
            place.blocks.emplace_back(anchor.start, end_of(anchor));
        }
        place.end = end_of(anchor);
        place.matched += anchor.to - anchor.from;
        before = &anchor;
    }
    return place;
}

/**
 * Chain the anchors of one feature, of the given length, and add each chain
 * that holds at least min_match of the feature's bases to the places.
 *
 * The anchors are sorted by their start on the path. Chains are taken best first, and each
 * anchor goes in one chain, so a later chain stops where it reaches an anchor
 * already taken.
 */
void take_chains(std::span<const Anchor> anchors, size_t length, double min_match,
                 std::vector<PathOccurrence>& places) {
    std::vector<ChainEnd> ends = best_chains(anchors, length);
    std::vector<size_t> order(anchors.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return ends[a].score > ends[b].score;
    });
    std::vector<bool> used(anchors.size(), false);
    for (size_t last : order) {
        std::vector<size_t> chain;
        for (size_t i = last; i != NO_ANCHOR && !used[i]; i = ends[i].previous) {
            chain.push_back(i);
            used[i] = true;
        }
        if (chain.empty()) {
            continue;
        }
        std::reverse(chain.begin(), chain.end());
        PathOccurrence place = chain_place(anchors, chain);
        if (static_cast<double>(place.matched) >= min_match * static_cast<double>(length)) {
            places.push_back(std::move(place));
        }
    }
}

}

std::vector<PathOccurrence> chain_anchors(std::vector<FeatureNodeIndex::Anchor> anchors,
                                          const std::vector<GraphFeature>& features, double min_match) {
    std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
        return std::tie(a.feature, a.start, a.from) < std::tie(b.feature, b.start, b.from);
    });
    std::vector<PathOccurrence> result;
    auto group = anchors.begin();
    while (group != anchors.end()) {
        size_t feature = group->feature;
        auto group_end = std::find_if(group, anchors.end(), [&](const Anchor& anchor) {
            return anchor.feature != feature;
        });
        take_chains(std::span<const Anchor>(group, group_end), features[feature].interval.length, min_match, result);
        group = group_end;
    }
    return result;
}

//------------------------------------------------------------------------------

namespace {

/**
 * Search each GBWT path whose name passes the filter with find(steps, coords,
 * places), and call visit(name, places) with the places sorted by start.
 *
 * The coords hold the coordinate of each step on the full path. Paths are
 * searched in batches with OpenMP threads. Each batch is reported in path
 * order before the next starts.
 */
void sweep_paths(const gbwtgraph::GBZ& gbz, const std::function<bool(const std::string&)>& include,
                 const std::function<void(const gbwt::vector_type&, const std::vector<uint64_t>&,
                                          std::vector<PathOccurrence>&)>& find,
                 const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit) {
    const GBWTGraph& graph = gbz.graph;
    struct PathStart {
        gbwt::size_type id;
        std::string name;
        uint64_t start;
    };
    std::vector<PathStart> paths;
    for_each_path_start(gbz.index, [&](gbwt::size_type path_id, const std::string& name, uint64_t start) {
        if (include(name)) {
            paths.emplace_back(path_id, name, start);
        }
    });

    size_t batch_size = 4 * omp_get_max_threads();
    std::vector<std::vector<PathOccurrence>> found(batch_size);
    for (size_t batch_start = 0; batch_start < paths.size(); batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, paths.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = batch_start; i < batch_end; i++) {
            const PathStart& path = paths[i];
            gbwt::vector_type steps = gbz.index.extract(gbwt::Path::encode(path.id, false));
            std::vector coords { path.start };
            for (gbwt::node_type node : steps) {
                coords.push_back(coords.back() + graph.get_length(GBWTGraph::node_to_handle(node)));
            }
            std::vector<PathOccurrence>& places = found[i - batch_start];
            places.clear();
            find(steps, coords, places);
            std::sort(places.begin(), places.end(), [](const PathOccurrence& a, const PathOccurrence& b) {
                return std::tie(a.start, a.end, a.feature, a.is_reverse) < std::tie(b.start, b.end, b.feature, b.is_reverse);
            });
        }
        for (size_t i = batch_start; i < batch_end; i++) {
            visit(paths[i].name, found[i - batch_start]);
        }
    }
}

}

void for_each_path_occurrences(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features,
                               const std::function<bool(const std::string&)>& include,
                               const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit) {
    FeatureMatcher matcher(gbz.graph, features);
    sweep_paths(gbz, include, [&](const gbwt::vector_type& steps, const std::vector<uint64_t>& coords,
                                  std::vector<PathOccurrence>& places) {
        matcher.for_each_match(steps, [&](size_t feature, size_t step, bool is_reverse) {
            const WalkInterval& interval = features[feature].interval;
            uint64_t walk_start = coords[step];
            uint64_t walk_end = coords[step + interval.walk.size()];
            uint64_t start = is_reverse ? walk_end - interval.offset - interval.length : walk_start + interval.offset;
            uint64_t end = start + interval.length;
            places.push_back({ start, end, feature, is_reverse, interval.length, { { start, end } } });
        });
    }, visit);
}

void for_each_path_liftover(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features, double min_match,
                            const std::function<bool(const std::string&)>& include,
                            const std::function<void(const std::string&, const std::vector<PathOccurrence>&)>& visit) {
    FeatureNodeIndex index(gbz.graph, features);
    sweep_paths(gbz, include, [&](const gbwt::vector_type& steps, const std::vector<uint64_t>& coords,
                                  std::vector<PathOccurrence>& places) {
        std::vector<FeatureNodeIndex::Anchor> anchors;
        index.for_each_anchor(steps, coords, [&](const FeatureNodeIndex::Anchor& anchor) {
            anchors.push_back(anchor);
        });
        places = chain_anchors(std::move(anchors), features, min_match);
    }, visit);
}

}
