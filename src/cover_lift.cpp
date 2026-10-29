/**
 * \file cover_lift.cpp
 * Lifting features from a cover, see cover_lift.hpp.
 */

#include "cover_lift.hpp"

#include <algorithm>
#include <istream>
#include <sstream>
#include <stdexcept>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

std::vector<CoverFeature> read_bed_features(std::istream& in) {
    std::vector<CoverFeature> features;
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
        CoverFeature feature;
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

//------------------------------------------------------------------------------

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

Cover::Cover(const GBWTGraph& graph, const GrefIndex* gref) : graph(graph), gref(gref) {
}

void Cover::add_tigs(std::istream& gaf, const std::function<bool(const std::string&)>& wanted) {
    std::string line;
    std::string name;
    while (std::getline(gaf, line)) {
        if (line.empty() || line[0] == '#') {
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
    if (gref != nullptr) {
        uint32_t path = gref->find_path(name);
        if (path != GrefIndex::NO_PATH) {
            return gref->path_length(path);
        }
    }
    auto found = tigs.find(name);
    if (found != tigs.end()) {
        return found->second.length;
    }
    return std::nullopt;
}

GraphFeature Cover::lift(const CoverFeature& feature) const {
    GraphFeature result;
    if (is_gref_path(feature.sequence)) {
        result.interval = gref->interval(gref->find_path(feature.sequence), feature.start, feature.end);
    } else {
        result.interval = tigs.at(feature.sequence).subinterval(graph, feature.start, feature.end);
    }
    result.fields = feature.fields;
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

void for_each_path_occurrences(const gbwtgraph::GBZ& gbz, const std::vector<GraphFeature>& features,
                               const std::function<bool(const std::string&)>& include,
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

    FeatureMatcher matcher(graph, features);
    // Search a batch of paths in parallel, then report them in order.
    size_t batch_size = 4 * omp_get_max_threads();
    std::vector<std::vector<PathOccurrence>> found(batch_size);
    for (size_t batch_start = 0; batch_start < paths.size(); batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, paths.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = batch_start; i < batch_end; i++) {
            const PathStart& path = paths[i];
            gbwt::vector_type steps = gbz.index.extract(gbwt::Path::encode(path.id, false));
            // Path coordinate of the first base of each step, and of the end of the path.
            std::vector<uint64_t> coords { path.start };
            for (gbwt::node_type node : steps) {
                coords.push_back(coords.back() + graph.get_length(GBWTGraph::node_to_handle(node)));
            }
            std::vector<PathOccurrence>& occurrences = found[i - batch_start];
            occurrences.clear();
            matcher.for_each_match(steps, [&](size_t feature, size_t step, bool is_reverse) {
                const WalkInterval& interval = features[feature].interval;
                uint64_t walk_start = coords[step], walk_end = coords[step + interval.walk.size()];
                uint64_t start = is_reverse ? walk_end - interval.offset - interval.length : walk_start + interval.offset;
                occurrences.push_back({ start, feature, is_reverse });
            });
            std::sort(occurrences.begin(), occurrences.end(), [](const PathOccurrence& a, const PathOccurrence& b) {
                if (a.start != b.start) {
                    return a.start < b.start;
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
