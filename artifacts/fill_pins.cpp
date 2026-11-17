#include "fill_pins.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

#include <omp.h>

#include "parallel_sort.hpp"
#include "site_values.hpp"
#include "../src/statistics.hpp"

namespace vg {

using namespace std;

void FillPointTally::add(const FillPointTally& other) {
    unitigs += other.unitigs;
    placed += other.placed;
    dropped_at_end += other.dropped_at_end;
    merged += other.merged;
}

namespace {

/// Return the one handle across the right side of `h`, or the left side with `go_left`, or nothing
/// where that side has no edge or more than one.
std::optional<handle_t> only_neighbour(const HandleGraph& graph, const handle_t& h, bool go_left) {
    std::optional<handle_t> found;
    size_t count = 0;
    graph.follow_edges(h, go_left, [&](const handle_t& next) {
        found = next;
        return ++count < 2;
    });
    if (count != 1) {
        return std::nullopt;
    }
    return found;
}

/// Return the handle that the unitig holding `h` continues into across the right side of `h`, or
/// nothing where the unitig ends there.
std::optional<handle_t> unitig_next(const HandleGraph& graph, const handle_t& h) {
    const std::optional<handle_t> next = only_neighbour(graph, h, false);
    if (!next || graph.get_id(*next) == graph.get_id(h)
        || !only_neighbour(graph, *next, true)) {
        return std::nullopt;
    }
    return next;
}

/// Return the unitig that starts at `start` and runs across its right side, stopping before it
/// would come back to the node of `start`.
vector<handle_t> unitig_from(const HandleGraph& graph, const handle_t& start) {
    vector<handle_t> unitig{start};
    const nid_t first = graph.get_id(start);
    // A join uses one side of each of its two nodes, and a side is in at most one join, so a
    // unitig can come back only to the node it started from. The cap on its length is a guard.
    const size_t most = graph.get_node_count();
    while (const std::optional<handle_t> next = unitig_next(graph, unitig.back())) {
        if (graph.get_id(*next) == first || unitig.size() >= most) {
            break;
        }
        unitig.push_back(*next);
    }
    return unitig;
}

/**
 * Place the points of the unitigs that are cycles, which have no end to start from.
 *
 * Each cycle is read forward from its node with the smallest ID. Runs on one thread, and visits
 * every node whose two sides both join it to other nodes.
 */
void place_cycle_points(const HandleGraph& graph, size_t spacing, vector<FillPoint>& out,
                        FillPointTally& tally) {
    const nid_t low = graph.min_node_id();
    vector<bool> seen((size_t)(graph.max_node_id() - low + 1), false);
    for (nid_t id = low; id <= graph.max_node_id(); ++id) {
        if (!graph.has_node(id) || seen[(size_t)(id - low)]) {
            continue;
        }
        const handle_t forward = graph.get_handle(id, false);
        if (!unitig_next(graph, forward) || !unitig_next(graph, graph.flip(forward))) {
            // The node ends its unitig, which is a path, not a cycle.
            continue;
        }
        // Nodes are visited in ID order, so the first node of a cycle reached is its smallest.
        vector<handle_t> unitig{forward};
        seen[(size_t)(id - low)] = true;
        bool cycle = false;
        while (const std::optional<handle_t> next = unitig_next(graph, unitig.back())) {
            const nid_t next_id = graph.get_id(*next);
            if (next_id == id) {
                cycle = true;
                break;
            }
            if (seen[(size_t)(next_id - low)]) {
                break;
            }
            seen[(size_t)(next_id - low)] = true;
            unitig.push_back(*next);
        }
        if (cycle) {
            place_unitig_points(graph, unitig, spacing, out, tally);
        }
    }
}

}

vector<handle_t> canonical_unitig(const HandleGraph& graph, vector<handle_t> unitig) {
    if (unitig.size() == 1) {
        unitig[0] = graph.forward(unitig[0]);
    } else if (unitig.size() > 1 && graph.get_id(unitig.back()) < graph.get_id(unitig.front())) {
        std::ranges::reverse(unitig);
        for (handle_t& h : unitig) {
            h = graph.flip(h);
        }
    }
    return unitig;
}

void place_unitig_points(const HandleGraph& graph, const vector<handle_t>& unitig, size_t spacing,
                         vector<FillPoint>& out, FillPointTally& tally) {
    if (spacing == 0 || unitig.empty()) {
        return;
    }
    // Where each node starts along the unitig.
    vector<size_t> starts(unitig.size());
    size_t length = 0;
    for (size_t i = 0; i < unitig.size(); ++i) {
        starts[i] = length;
        length += graph.get_length(unitig[i]);
    }
    if (length <= spacing) {
        return;
    }
    ++tally.unitigs;
    // The number of equal parts the points divide the unitig into, one more than the points.
    const size_t parts = (length + spacing - 1) / spacing;
    // The node holding the current point, which only moves forward.
    size_t node = 0;
    // The bases before the last point placed, or 0 before the first.
    size_t last = 0;
    for (size_t i = 1; i < parts; ++i) {
        // i/parts of the length, rounded to the nearest base, half up.
        size_t bases_before = (2 * i * length + parts) / (2 * parts);
        if (bases_before <= last) {
            // The point before it moved forward off a join to here or past it, and this point
            // would move to the same place.
            ++tally.merged;
            continue;
        }
        while (bases_before < length) {
            while (node + 1 < unitig.size() && starts[node + 1] <= bases_before) {
                ++node;
            }
            if (bases_before > starts[node]) {
                break;
            }
            // The point is at the join before the node, so move it forward one base.
            ++bases_before;
        }
        if (bases_before >= length) {
            ++tally.dropped_at_end;
            continue;
        }
        const handle_t& h = unitig[node];
        const size_t into = bases_before - starts[node];
        const size_t node_length = graph.get_length(h);
        out.push_back(FillPoint{graph.get_id(h),
                                (uint32_t)(graph.get_is_reverse(h) ? node_length - into : into)});
        ++tally.placed;
        last = bases_before;
    }
}

vector<FillPoint> place_fill_points(const HandleGraph& graph, size_t spacing,
                                    FillPointTally& tally) {
    tally = FillPointTally();
    vector<FillPoint> points;
    if (spacing == 0 || graph.get_node_count() == 0) {
        return points;
    }
    const nid_t low = graph.min_node_id();
    const nid_t high = graph.max_node_id();
    const auto threads = (size_t)max(1, omp_get_max_threads());
    vector<vector<FillPoint>> found(threads);
    vector<FillPointTally> tallies(threads);
    std::atomic<size_t> covered{0};
    // Each unitig that is a path is placed once, from its end whose node has the smaller ID. Nodes
    // inside a unitig are skipped, so a unitig that is a cycle, which has no end, is not reached.
#pragma omp parallel for schedule(dynamic, 4096)
    for (nid_t id = low; id <= high; ++id) {
        if (!graph.has_node(id)) {
            continue;
        }
        const handle_t forward = graph.get_handle(id, false);
        const bool joins_left = unitig_next(graph, graph.flip(forward)).has_value();
        if (joins_left && unitig_next(graph, forward).has_value()) {
            continue;
        }
        const vector<handle_t> unitig = unitig_from(graph, joins_left ? graph.flip(forward) : forward);
        if (unitig.size() > 1 && graph.get_id(unitig.back()) < id) {
            continue;
        }
        covered += unitig.size();
        const auto thread = (size_t)omp_get_thread_num() % threads;
        place_unitig_points(graph, canonical_unitig(graph, unitig), spacing, found[thread],
                            tallies[thread]);
    }
    for (size_t t = 0; t < threads; ++t) {
        points.insert(points.end(), found[t].begin(), found[t].end());
        tally.add(tallies[t]);
    }
    if (covered.load() < graph.get_node_count()) {
        place_cycle_points(graph, spacing, points, tally);
    }
    parallel_sort(points.begin(), points.end(), std::less<FillPoint>());
    return points;
}

vector<pair<int, uint8_t>> fill_marks_of_call(const vector<int>& genotype, int nested_strand) {
    vector<pair<int, uint8_t>> marks;
    if (genotype.size() == 1) {
        const bool on_one_strand = nested_strand == 0 || nested_strand == 1;
        marks.emplace_back(genotype[0],
                           on_one_strand ? (uint8_t)(1 << nested_strand) : FILL_EVERY_STRAND);
    } else if (genotype.size() == 2 && genotype[0] != genotype[1]) {
        marks.emplace_back(genotype[0], FILL_STRAND_0);
        marks.emplace_back(genotype[1], FILL_STRAND_1);
    } else {
        // A homozygote's one allele is on every strand. So, for want of a slot for each, is each
        // allele of a genotype of more than two strands.
        for (int allele : genotype) {
            if (std::ranges::none_of(marks, [&](const auto& mark) { return mark.first == allele; })) {
                marks.emplace_back(allele, FILL_EVERY_STRAND);
            }
        }
    }
    return marks;
}

FillSlot classify_fill_point(bool held_by_root, uint8_t mark) {
    const bool strand_0 = (mark & FILL_STRAND_0) != 0;
    const bool strand_1 = (mark & FILL_STRAND_1) != 0;
    if (held_by_root || (mark & FILL_EVERY_STRAND) != 0 || (strand_0 && strand_1)) {
        return FillSlot{AnchorKind::fill_shared, 0};
    }
    if (strand_0 || strand_1) {
        return FillSlot{AnchorKind::fill_strand, strand_0 ? 0 : 1};
    }
    return FillSlot{AnchorKind::fill_off_call, 0};
}

float fill_read_score(int32_t mapq, const FillScoreParams& params) {
    double mismap = params.use_mapping_quality ? phred_to_prob((double)mapq) : params.min_mismap;
    mismap = min(max(mismap, params.min_mismap), params.max_mismap);
    return mismap > 0.0 ? (float)(-10.0 * log10(mismap)) : 99.0f;
}

void FillAnchors::configure(size_t point_spacing, const FillScoreParams& read_scores,
                            double min_score) {
    spacing = point_spacing;
    scores = read_scores;
    min_read_score = min_score;
}

void FillAnchors::place(const HandleGraph& point_graph, const SnarlDecomposition& decomposition) {
    graph = &point_graph;
    sites = &decomposition;
    fill_points = place_fill_points(point_graph, spacing, point_tally);
    owners.assign(fill_points.size(), net_handle_t());
    held_by_root.assign(fill_points.size(), 0);
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < fill_points.size(); ++i) {
        const net_handle_t site =
            multipass::placement_of_node(decomposition, point_graph, fill_points[i].node).site;
        owners[i] = decomposition.canonical(site);
        held_by_root[i] = decomposition.is_root(site) ? 1 : 0;
    }
    marks = vector<std::atomic<uint8_t>>(fill_points.size());
}

namespace {

/// Set `visited` to the nodes of `nodes`, which is sorted, that `read` visits, once each and
/// sorted.
void nodes_visited(const SiteRead& read, const vector<nid_t>& nodes, vector<nid_t>& visited) {
    const Path& path = read.aln->path();
    visited.clear();
    if (read.indexed()) {
        // The index holds the read's mappings on the queried nodes, which are `nodes`.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            visited.push_back(path.mapping((int)read.mappings[k]).position().node_id());
        }
    } else {
        for (const Mapping& m : path.mapping()) {
            if (std::ranges::binary_search(nodes, m.position().node_id())) {
                visited.push_back(m.position().node_id());
            }
        }
    }
    std::ranges::sort(visited);
    const auto repeats = std::ranges::unique(visited);
    visited.erase(repeats.begin(), repeats.end());
}

}

void FillAnchors::place_read(const SiteRead& read, const vector<nid_t>& visited, size_t begin,
                             size_t end, vector<FillPlacement>& out) {
    const float score = fill_read_score(read.aln->mapping_quality(), scores);
    // The read's index in read_names(), found at its first placement.
    auto name = numeric_limits<uint32_t>::max();
    const auto first = fill_points.begin() + (ptrdiff_t)begin;
    const auto last = fill_points.begin() + (ptrdiff_t)end;
    for (nid_t id : visited) {
        for (auto point = std::lower_bound(first, last, FillPoint{id, 0});
             point != last && point->node == id; ++point) {
            const AnchorPlacement placed =
                resolve_anchor_pin(read, *graph, id, false, true, read_counters, point->pin);
            if (!placed.placed()) {
                continue;
            }
            ++resolved;
            // A check of the read around the point, such as exact matches to the graph over a
            // few bases on each side, would go here.
            if (score < min_read_score) {
                ++low_score;
                continue;
            }
            if (name == numeric_limits<uint32_t>::max()) {
                name = read_names().intern(read.aln->name());
            }
            out.push_back(FillPlacement{(uint32_t)(point - fill_points.begin()), name,
                                        placed.alignment, (uint32_t)placed.offset, score,
                                        placed.direction});
        }
    }
}

vector<FillPlacement> FillAnchors::resolve_window(const SiteReadSource& reads, size_t begin,
                                                  size_t end) {
    // The points' nodes, which are sorted, and the ranges of node IDs they make.
    vector<nid_t> nodes;
    vector<pair<nid_t, nid_t>> ranges;
    for (size_t i = begin; i < end; ++i) {
        const nid_t id = fill_points[i].node;
        if (!nodes.empty() && nodes.back() == id) {
            continue;
        }
        nodes.push_back(id);
        if (!ranges.empty() && ranges.back().second + 1 == id) {
            ranges.back().second = id;
        } else {
            ranges.emplace_back(id, id);
        }
    }
    vector<FillPlacement> out;
    vector<nid_t> visited;
    reads.for_each_read_uncounted(ranges, [&](const SiteRead& read) {
        nodes_visited(read, nodes, visited);
        place_read(read, visited, begin, end, out);
    });
    return out;
}

void FillAnchors::resolve_reads(const SiteReadSource& reads, size_t window) {
    placed_reads.clear();
    if (fill_points.empty()) {
        return;
    }
    window = max<size_t>(window, 1);
    // The points of each window of node IDs, as [begin, end) in `fill_points`.
    vector<pair<size_t, size_t>> windows;
    for (size_t i = 0; i < fill_points.size();) {
        const nid_t w = fill_points[i].node / (nid_t)window;
        size_t j = i;
        while (j < fill_points.size() && fill_points[j].node / (nid_t)window == w) {
            ++j;
        }
        windows.emplace_back(i, j);
        i = j;
    }
    // The placements of each window, kept apart so that their order does not depend on the
    // threads.
    vector<vector<FillPlacement>> found(windows.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t w = 0; w < windows.size(); ++w) {
        found[w] = resolve_window(reads, windows[w].first, windows[w].second);
    }
    size_t total = 0;
    for (const auto& placements : found) {
        total += placements.size();
    }
    placed_reads.reserve(total);
    for (auto& placements : found) {
        placed_reads.insert(placed_reads.end(), placements.begin(), placements.end());
        vector<FillPlacement>().swap(placements);
    }
    std::ranges::stable_sort(placed_reads, {}, &FillPlacement::point);
}

void FillAnchors::mark_walk(const net_handle_t& site, const Traversal& walk, uint8_t bits) {
    if (fill_points.empty() || sites == nullptr) {
        return;
    }
    const net_handle_t owner = sites->canonical(site);
    for (const handle_t& h : walk) {
        const nid_t id = graph->get_id(h);
        for (auto point = std::ranges::lower_bound(fill_points, FillPoint{id, 0});
             point != fill_points.end() && point->node == id; ++point) {
            const auto i = (size_t)(point - fill_points.begin());
            if (held_by_root[i] == 0 && owners[i] == owner) {
                marks[i].fetch_or(bits, std::memory_order_relaxed);
            }
        }
    }
}

FillSlot FillAnchors::slot_of(size_t point) const {
    return classify_fill_point(held_by_root[point] != 0,
                               marks[point].load(std::memory_order_relaxed));
}

void FillAnchors::add_anchors(AnchorCollector& anchors) {
    kind_counts.fill(0);
    const ReadNameTable& names = read_names();
    for (size_t i = 0; i < placed_reads.size();) {
        const uint32_t p = placed_reads[i].point;
        size_t j = i;
        while (j < placed_reads.size() && placed_reads[j].point == p) {
            ++j;
        }
        const FillSlot fill_slot = slot_of(p);
        ++kind_counts[(size_t)fill_slot.kind - (size_t)AnchorKind::fill_shared];

        AnchorWriter::Anchor anchor;
        anchor.node = fill_points[p].node;
        anchor.pin = fill_points[p].pin;
        anchor.snarl = "fill" + to_string(anchor.node) + "." + to_string(anchor.pin);
        anchor.slot = fill_slot.slot;
        anchor.allele = -1;
        anchor.gqn = numeric_limits<double>::quiet_NaN();
        anchor.explained = numeric_limits<double>::quiet_NaN();
        anchor.reliability = -1.0;
        anchor.source = PinSource::fill;
        anchor.kind = fill_slot.kind;
        anchor.reads.reserve(j - i);
        for (size_t k = i; k < j; ++k) {
            const FillPlacement& placement = placed_reads[k];
            AnchorWriter::ReadRow row;
            row.read = placement.read;
            row.alignment = placement.alignment;
            row.offset = placement.offset;
            row.score = placement.score;
            row.direction = placement.direction;
            anchor.reads.push_back(row);
        }
        // Order the reads by everything the file writes for them, so that the file does not
        // depend on the order the threads found them in.
        std::ranges::sort(anchor.reads, [&names](const AnchorWriter::ReadRow& x,
                                                 const AnchorWriter::ReadRow& y) {
            if (x.read != y.read) {
                return names.name(x.read) < names.name(y.read);
            }
            return tie(x.alignment, x.direction, x.offset) < tie(y.alignment, y.direction, y.offset);
        });
        anchors.add(std::move(anchor));
        i = j;
    }
    vector<FillPlacement>().swap(placed_reads);
}

void FillAnchors::report(ostream& out) const {
    out << "[vg call] fill anchors: " << point_tally.placed << " points every " << spacing
        << " bp in " << point_tally.unitigs << " unbranched node runs longer than that; "
        << point_tally.dropped_at_end << " dropped at a run's end and " << point_tally.merged
        << " at the next point, when moved forward off a join between nodes" << endl;
    const AnchorPinCounters& pins = read_counters.pins;
    out << "[vg call] fill anchors: " << resolved.load() << " read placements resolved, "
        << pins.verified.load() << " checked against the graph's base and "
        << pins.verify_failed.load() << " failing the check; refused: "
        << pins.repeat_visit.load() << " reads visit the point's node more than once, "
        << pins.unaligned_base.load() << " have no read base aligned to the node base before the "
        << "point; " << low_score.load() << " placements below the minimum read score" << endl;
    out << "[vg call] fill anchors: points with reads: " << kind_counts[0] << " fill-shared, "
        << kind_counts[1] << " fill-strand, " << kind_counts[2] << " fill-off-call" << endl;
}

}
