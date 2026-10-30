/**
 * \file feature_lift.cpp
 * Lifting features, see feature_lift.hpp.
 */

#include "feature_lift.hpp"

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <istream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

std::vector<BedFeature> read_bed_features(std::istream& in) {
    std::vector<BedFeature> features;
    std::string line;
    size_t line_number = 0;
    while (std::getline(in, line)) {
        line_number++;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0 || line.compare(0, 7, "browser") == 0) {
            continue;
        }
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
            size_t used_start = 0, used_end = 0;
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
        features.push_back(std::move(feature));
    }
    return features;
}

bool GraphFeature::is_stranded() const {
    return fields.size() > BED_STRAND_FIELD && fields[BED_STRAND_FIELD] != ".";
}

namespace {

/// Return -1, 0 or 1 as a sorts before, with or after b. Intervals sort by
/// their walks as GBWT node numbers, then by offset and length.
int compare_intervals(const WalkInterval& a, const WalkInterval& b) {
    for (size_t i = 0; i < std::min(a.walk.size(), b.walk.size()); i++) {
        gbwt::node_type x = GBWTGraph::handle_to_node(a.walk[i]), y = GBWTGraph::handle_to_node(b.walk[i]);
        if (x != y) {
            return x < y ? -1 : 1;
        }
    }
    if (a.walk.size() != b.walk.size()) {
        return a.walk.size() < b.walk.size() ? -1 : 1;
    }
    if (a.offset != b.offset) {
        return a.offset < b.offset ? -1 : 1;
    }
    if (a.length != b.length) {
        return a.length < b.length ? -1 : 1;
    }
    return 0;
}

}

bool PathCoordinate::operator<(const PathCoordinate& other) const {
    return std::tie(rank, start, end, is_reverse) < std::tie(other.rank, other.start, other.end, other.is_reverse);
}

std::string PathCoordinate::text() const {
    return path + ':' + std::to_string(start + 1) + '-' + std::to_string(end) + ':' + (is_reverse ? '-' : '+');
}

GraphFeature make_graph_feature(const GBWTGraph& graph, const GrefIndex* gref,
                                WalkInterval interval, std::vector<std::string> fields) {
    GraphFeature result;
    result.interval = std::move(interval);
    result.fields = std::move(fields);
    bool stranded = result.is_stranded();
    if (stranded && result.fields[BED_STRAND_FIELD] == "-") {
        result.interval = result.interval.reverse(graph);
        result.fields[BED_STRAND_FIELD] = "+";
    }

    if (gref != nullptr) {
        result.placement = gref->place(result.interval);
    }
    if (!stranded) {
        const GrefPlacement& placement = result.placement;
        WalkInterval reverse = result.interval.reverse(graph);
        bool flip = placement.reverse_bases > placement.forward_bases ||
                    (placement.reverse_bases == placement.forward_bases &&
                     compare_intervals(reverse, result.interval) < 0);
        if (flip) {
            result.interval = std::move(reverse);
            if (gref != nullptr) {
                result.placement = gref->place(result.interval);
            }
        }
    }

    // Gref paths come first in the priority order, and their ranks are their
    // places in the gref.
    if (result.placement.along_one_path) {
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
        if (!tigs.emplace(name, std::move(interval)).second) {
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
    auto found = tigs.find(name);
    if (found != tigs.end()) {
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
        if (name.compare(0, gref_prefix.size(), gref_prefix) == 0 || groups.count(name)) {
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
            fragments[name].push_back({ start, PathSteps(gbz.graph, path_id) });
        }
    });
    for (auto& entry : fragments) {
        std::sort(entry.second.begin(), entry.second.end(), [](const Fragment& a, const Fragment& b) {
            return a.start < b.start;
        });
    }
}

bool PathDirectory::has_path(const std::string& name) const {
    return fragments.count(name) > 0;
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
        int order = compare_intervals(a.interval, b.interval);
        if (order != 0) {
            return order < 0;
        }
        return a.fields < b.fields;
    });
    auto end = std::unique(features.begin(), features.end(), [](const GraphFeature& a, const GraphFeature& b) {
        return compare_intervals(a.interval, b.interval) == 0 && a.fields == b.fields;
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
            const std::vector<handle_t>& walk = features[feature].interval.walk;
            if (step + walk.size() > path.size()) {
                continue;
            }
            bool matches = true;
            for (size_t j = 0; j < walk.size() && matches; j++) {
                handle_t expected = is_reverse ? graph.flip(walk[walk.size() - 1 - j]) : walk[j];
                matches = (GBWTGraph::node_to_handle(path[step + j]) == expected);
            }
            if (matches) {
                visit(feature, step, is_reverse);
            }
        }
    }
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
            size_t high = std::min(node_end, interval.offset + interval.length);
            if (low < high) {
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

std::vector<PathOccurrence> chain_anchors(std::vector<FeatureNodeIndex::Anchor> anchors,
                                          const std::vector<GraphFeature>& features, double min_match) {
    using Anchor = FeatureNodeIndex::Anchor;
    std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
        if (a.feature != b.feature) {
            return a.feature < b.feature;
        }
        if (a.start != b.start) {
            return a.start < b.start;
        }
        return a.from < b.from;
    });
    auto end_of = [](const Anchor& anchor) {
        return anchor.start + (anchor.to - anchor.from);
    };

    std::vector<PathOccurrence> result;
    constexpr size_t NONE = std::numeric_limits<size_t>::max();
    for (size_t begin = 0, end = 0; begin < anchors.size(); begin = end) {
        size_t feature = anchors[begin].feature;
        while (end < anchors.size() && anchors[end].feature == feature) {
            end++;
        }
        size_t length = features[feature].interval.length;

        // Find the chain ending at each anchor that holds the most feature
        // bases. Its previous anchor ends before this one starts on the path,
        // with a gap of at most the feature's length. A path that reads the
        // feature in reverse visits its bases in decreasing order.
        size_t count = end - begin;
        std::vector<size_t> score(count), previous(count, NONE);
        for (size_t i = 0; i < count; i++) {
            const Anchor& anchor = anchors[begin + i];
            score[i] = anchor.to - anchor.from;
            for (size_t j = i; j-- > 0;) {
                const Anchor& before = anchors[begin + j];
                if (before.start + 2 * length < anchor.start) {
                    // This anchor and all earlier ones end more than the
                    // feature's length back.
                    break;
                }
                bool continues = before.is_reverse == anchor.is_reverse && end_of(before) <= anchor.start &&
                                 anchor.start - end_of(before) <= length &&
                                 (anchor.is_reverse ? anchor.to <= before.from : anchor.from >= before.to);
                if (continues && score[j] + (anchor.to - anchor.from) > score[i]) {
                    score[i] = score[j] + (anchor.to - anchor.from);
                    previous[i] = j;
                }
            }
        }

        // Take the chains best first. Each anchor goes in one chain, so a
        // later chain stops where it reaches an anchor already taken.
        std::vector<size_t> order(count);
        for (size_t i = 0; i < count; i++) {
            order[i] = i;
        }
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return score[a] > score[b];
        });
        std::vector<bool> used(count, false);
        for (size_t last : order) {
            std::vector<size_t> chain;
            for (size_t i = last; i != NONE && !used[i]; i = previous[i]) {
                chain.push_back(i);
                used[i] = true;
            }
            if (chain.empty()) {
                continue;
            }
            std::reverse(chain.begin(), chain.end());
            const Anchor& first = anchors[begin + chain.front()];
            PathOccurrence place { first.start, first.start, feature, first.is_reverse, 0, {} };
            const Anchor* before = nullptr;
            for (size_t i : chain) {
                const Anchor& anchor = anchors[begin + i];
                bool adjacent = before != nullptr && anchor.start == place.end &&
                                (anchor.is_reverse ? anchor.to == before->from : anchor.from == before->to);
                if (adjacent) {
                    place.blocks.back().second = end_of(anchor);
                } else {
                    place.blocks.emplace_back(anchor.start, end_of(anchor));
                }
                place.end = end_of(anchor);
                place.matched += anchor.to - anchor.from;
                before = &anchor;
            }
            if (place.matched >= min_match * length) {
                result.push_back(std::move(place));
            }
        }
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
            paths.push_back({ path_id, name, start });
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
            std::vector<uint64_t> coords { path.start };
            for (gbwt::node_type node : steps) {
                coords.push_back(coords.back() + graph.get_length(GBWTGraph::node_to_handle(node)));
            }
            std::vector<PathOccurrence>& places = found[i - batch_start];
            places.clear();
            find(steps, coords, places);
            std::sort(places.begin(), places.end(), [](const PathOccurrence& a, const PathOccurrence& b) {
                if (a.start != b.start) {
                    return a.start < b.start;
                }
                if (a.end != b.end) {
                    return a.end < b.end;
                }
                if (a.feature != b.feature) {
                    return a.feature < b.feature;
                }
                return a.is_reverse < b.is_reverse;
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
            uint64_t walk_start = coords[step], walk_end = coords[step + interval.walk.size()];
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
