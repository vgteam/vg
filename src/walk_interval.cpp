/**
 * \file walk_interval.cpp
 * Walk intervals and path steps, see walk_interval.hpp.
 */

#include "walk_interval.hpp"

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

std::string WalkInterval::sequence(const GBWTGraph& graph) const {
    std::string result;
    result.reserve(length);
    size_t skip = offset;
    for (const handle_t& handle : walk) {
        std::string_view view = graph.get_sequence_view(handle);
        size_t take = std::min(view.size() - skip, length - result.size());
        result.append(view.substr(skip, take));
        skip = 0;
    }
    return result;
}

WalkInterval WalkInterval::reverse(const GBWTGraph& graph) const {
    WalkInterval result;
    result.walk.reserve(walk.size());
    for (auto it = walk.rbegin(); it != walk.rend(); ++it) {
        result.walk.push_back(graph.flip(*it));
    }
    result.offset = walk_length(graph) - (offset + length);
    result.length = length;
    return result;
}

WalkInterval WalkInterval::subinterval(const GBWTGraph& graph, size_t start, size_t end) const {
    WalkInterval result;
    result.length = end - start;
    // Walk coordinates of the first base and of the base after the last.
    size_t from = offset + start, to = offset + end;
    size_t node_start = 0;
    for (const handle_t& handle : walk) {
        size_t node_end = node_start + graph.get_length(handle);
        if (node_end > from) {
            if (result.walk.empty()) {
                result.offset = from - node_start;
            }
            result.walk.push_back(handle);
        }
        if (node_end >= to) {
            break;
        }
        node_start = node_end;
    }
    return result;
}

size_t WalkInterval::walk_length(const GBWTGraph& graph) const {
    size_t total = 0;
    for (const handle_t& handle : walk) {
        total += graph.get_length(handle);
    }
    return total;
}

void write_gaf_columns(std::ostream& out, const GBWTGraph& graph, const std::string& name,
                       const WalkInterval& interval) {
    out << name << '\t' << interval.length << "\t0\t" << interval.length << "\t+\t";
    for (const handle_t& handle : interval.walk) {
        out << (graph.get_is_reverse(handle) ? '<' : '>') << graph.get_id(handle);
    }
    out << '\t' << interval.walk_length(graph) << '\t' << interval.offset << '\t'
        << (interval.offset + interval.length) << '\t' << interval.length << '\t' << interval.length << "\t255";
}

WalkInterval parse_gaf_columns(const std::string& line, const GBWTGraph& graph, std::string& name) {
    auto fail = [&](const std::string& problem) {
        throw std::runtime_error(problem + " in GAF line: " + line);
    };
    std::vector<std::string> columns;
    std::stringstream stream(line);
    std::string column;
    while (columns.size() < 12 && std::getline(stream, column, '\t')) {
        columns.push_back(column);
    }
    if (columns.size() < 12) {
        fail("fewer than 12 columns");
    }
    auto number = [&](const std::string& text) -> size_t {
        size_t used = 0;
        size_t value = 0;
        try {
            value = std::stoull(text, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used == 0 || used != text.size()) {
            fail("bad number " + text);
        }
        return value;
    };

    name = columns[0];
    WalkInterval interval;
    const std::string& walk = columns[5];
    for (size_t i = 0; i < walk.size();) {
        if (walk[i] != '>' && walk[i] != '<') {
            fail("not a walk of node IDs");
        }
        size_t j = walk.find_first_of("<>", i + 1);
        if (j == std::string::npos) {
            j = walk.size();
        }
        nid_t id = number(walk.substr(i + 1, j - i - 1));
        if (!graph.has_node(id)) {
            fail("node " + std::to_string(id) + " not in the graph");
        }
        interval.walk.push_back(graph.get_handle(id, walk[i] == '<'));
        i = j;
    }
    size_t start = number(columns[7]), end = number(columns[8]);
    if (interval.walk.empty() || start >= end || end > interval.walk_length(graph)) {
        fail("bad walk or path coordinates");
    }
    interval.offset = start;
    interval.length = end - start;
    // Trim the walk to the nodes that hold the interval's bases.
    return interval.subinterval(graph, 0, interval.length);
}

void for_each_path_start(const gbwt::GBWT& index,
                         const std::function<void(gbwt::size_type, const std::string&, uint64_t)>& visit) {
    if (!index.hasMetadata() || !index.metadata.hasPathNames()) {
        return;
    }
    using handlegraph::PathMetadata;
    using handlegraph::PathSense;
    gbwtgraph::sample_name_set reference_samples = gbwtgraph::parse_reference_samples_tag(index);
    for (gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++) {
        PathSense sense = gbwtgraph::get_path_sense(index, path_id, reference_samples);
        std::string sample = gbwtgraph::get_path_sample_name(index, path_id, sense);
        std::string locus = gbwtgraph::get_path_locus_name(index, path_id, sense);
        size_t haplotype = gbwtgraph::get_path_haplotype(index, path_id, sense);
        uint64_t start = 0;
        if (sense == PathSense::HAPLOTYPE) {
            start = gbwtgraph::get_path_phase_block(index, path_id, sense);
            // Compose the name as a reference path name. That is the haplotype
            // path name without the phase block.
            sense = PathSense::REFERENCE;
        } else {
            handlegraph::subrange_t subrange = gbwtgraph::get_path_subrange(index, path_id, sense);
            if (subrange != PathMetadata::NO_SUBRANGE) {
                start = subrange.first;
            }
        }
        std::string name = PathMetadata::create_path_name(sense, sample, locus, haplotype,
                                                          PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE);
        visit(path_id, name, start);
    }
}

//------------------------------------------------------------------------------

PathSteps::PathSteps(const GBWTGraph& graph, gbwt::size_type path_id) {
    for (gbwt::node_type node : graph.index->extract(gbwt::Path::encode(path_id, false))) {
        handle_t handle = GBWTGraph::node_to_handle(node);
        step_list.push_back({ path_length, handle });
        path_length += graph.get_length(handle);
    }
}

namespace {

/// Return the index of the step that holds the given path coordinate.
size_t step_at(const std::vector<PathSteps::Step>& steps, uint64_t coord) {
    auto after = std::upper_bound(steps.begin(), steps.end(), coord, [](uint64_t c, const PathSteps::Step& step) {
        return c < step.start;
    });
    return after - steps.begin() - 1;
}

}

WalkInterval PathSteps::interval(uint64_t start, uint64_t end) const {
    size_t first = step_at(step_list, start), last = step_at(step_list, end - 1);
    WalkInterval result;
    for (size_t step = first; step <= last; step++) {
        result.walk.push_back(step_list[step].handle);
    }
    result.offset = start - step_list[first].start;
    result.length = end - start;
    return result;
}

char PathSteps::base(const GBWTGraph& graph, uint64_t coord) const {
    const Step& step = step_list[step_at(step_list, coord)];
    return graph.get_sequence_view(step.handle)[coord - step.start];
}

}
