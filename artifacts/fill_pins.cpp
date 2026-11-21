#include "fill_pins.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

#include <omp.h>

#include <handlegraph/util.hpp>

#include "decomposition_sites.hpp"
#include "parallel_sort.hpp"
#include "site_values.hpp"
#include "../src/path.hpp"
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

vector<pair<int, uint8_t>> fill_marks_of_call(const vector<int>& genotype, int nested_strand,
                                              bool nested) {
    vector<pair<int, uint8_t>> marks;
    if (genotype.size() == 1) {
        if (nested_strand == 0 || nested_strand == 1) {
            marks.emplace_back(genotype[0], (uint8_t)(1 << nested_strand));
        } else if (!nested) {
            marks.emplace_back(genotype[0], FILL_EVERY_STRAND);
        }
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

namespace {

/// Order net handles by their bits, which is enough to sort and search them.
bool net_before(const net_handle_t& a, const net_handle_t& b) {
    return handlegraph::as_integer(a) < handlegraph::as_integer(b);
}

}

void FillAnchors::place_points(const HandleGraph& point_graph, size_t window) {
    graph = &point_graph;
    fill_points = place_fill_points(point_graph, spacing, point_tally);
    window_size = max<size_t>(window, 1);
    windows.clear();
    for (size_t i = 0; i < fill_points.size();) {
        const nid_t number = fill_points[i].node / (nid_t)window_size;
        size_t j = i;
        while (j < fill_points.size() && fill_points[j].node / (nid_t)window_size == number) {
            ++j;
        }
        windows.push_back(PointWindow{number, i, j});
        i = j;
    }
    window_claimed = vector<std::atomic<bool>>(windows.size());
    window_placements = vector<vector<FillPlacement>>(windows.size());
}

void FillAnchors::place(const SnarlDecomposition& decomposition) {
    const HandleGraph& point_graph = *graph;
    sites = &decomposition;
    // The sites holding each point's node, innermost first. Most points lie in no site.
    vector<vector<net_handle_t>> held_in(fill_points.size());
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < fill_points.size(); ++i) {
        net_handle_t site =
            multipass::placement_of_node(decomposition, point_graph, fill_points[i].node).site;
        // The depth of a site among a point's sites is kept in a byte (see `marks`).
        while (!decomposition.is_root(site) && held_in[i].size() < 255) {
            held_in[i].push_back(decomposition.canonical(site));
            site = enclosing_site(decomposition, site);
        }
    }
    holding_sites.clear();
    for (const auto& sites_of_point : held_in) {
        holding_sites.insert(holding_sites.end(), sites_of_point.begin(), sites_of_point.end());
    }
    std::ranges::sort(holding_sites, net_before);
    const auto repeats = std::ranges::unique(holding_sites);
    holding_sites.erase(repeats.begin(), repeats.end());
    holders.clear();
    holder_start.assign(1, 0);
    for (auto& sites_of_point : held_in) {
        for (const net_handle_t& site : sites_of_point) {
            holders.push_back((uint32_t)(std::ranges::lower_bound(holding_sites, site, net_before)
                                         - holding_sites.begin()));
        }
        holder_start.push_back((uint32_t)holders.size());
        vector<net_handle_t>().swap(sites_of_point);
    }
    site_marked = vector<std::atomic<uint8_t>>(holding_sites.size());
    marks = vector<std::atomic<uint16_t>>(fill_points.size());
}

void FillAnchors::nodes_visited(const SiteRead& read, const vector<nid_t>& nodes,
                                vector<VisitedNode>& visited) {
    const Path& path = read.aln->path();
    visited.clear();
    auto add = [&](int i) {
        visited.push_back(VisitedNode{path.mapping(i).position().node_id(), i, 1});
    };
    if (read.indexed()) {
        // The index holds the read's mappings on the queried nodes, which are `nodes`.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            add((int)read.mappings[k]);
        }
    } else {
        for (int i = 0; i < path.mapping_size(); ++i) {
            if (std::ranges::binary_search(nodes, path.mapping(i).position().node_id())) {
                add(i);
            }
        }
    }
    std::ranges::sort(visited, {}, [](const VisitedNode& v) { return make_pair(v.node, v.mapping); });
    // Merge the visits to each node into the first.
    size_t kept = 0;
    for (size_t i = 0; i < visited.size(); ++i) {
        if (kept > 0 && visited[kept - 1].node == visited[i].node) {
            ++visited[kept - 1].mappings;
        } else {
            visited[kept++] = visited[i];
        }
    }
    visited.resize(kept);
}

namespace {

/// Return whether mapping `m` aligns its read to node base `base`, which is counted along the
/// forward strand of the node, of `length` bases. A deleted node base counts as aligned.
bool mapping_reaches(const Mapping& m, size_t length, size_t base) {
    const size_t visited_base = m.position().is_reverse() ? length - 1 - base : base;
    const auto first = (size_t)m.position().offset();
    return visited_base >= first && visited_base < first + (size_t)mapping_from_length(m);
}

/// Return whether the read of a placement at an exit pin inside a node has a base after the
/// point, in the direction the read crosses it, outside its soft clips.
bool goes_on_past(const Alignment& aln, const AnchorPlacement& placed) {
    const int64_t after = placed.direction == 0 ? placed.offset + 1 : placed.offset - 1;
    return after >= softclip_start(aln)
        && after < (int64_t)aln.sequence().size() - softclip_end(aln);
}

}

void FillAnchors::place_read(const SiteRead& read, const vector<VisitedNode>& visited,
                             size_t begin, size_t end, vector<FillPlacement>& out) {
    const float score = fill_read_score(read.aln->mapping_quality(), scores);
    const Path& path = read.aln->path();
    // The read's index in read_names(), found at its first placement.
    auto name = numeric_limits<uint32_t>::max();
    const auto first = fill_points.begin() + (ptrdiff_t)begin;
    const auto last = fill_points.begin() + (ptrdiff_t)end;
    for (const VisitedNode& visit : visited) {
        const nid_t id = visit.node;
        const size_t node_length = graph->get_length(graph->get_handle(id));
        for (auto point = std::lower_bound(first, last, FillPoint{id, 0});
             point != last && point->node == id; ++point) {
            if (visit.mappings == 1
                && !mapping_reaches(path.mapping(visit.mapping), node_length, point->pin - 1)) {
                // The read's one mapping on the node ends before the point or starts after it. A
                // read that visits the node more than once is left to `resolve_anchor_pin`.
                ++not_reached;
                continue;
            }
            const AnchorPlacement placed =
                resolve_anchor_pin(read, *graph, id, false, true, read_counters, point->pin);
            if (!placed.placed()) {
                continue;
            }
            // An exit pin needs a read base only before the point, so a read must also be checked
            // to cross it. A check of the read around the point, such as exact matches to the graph
            // over a few bases on each side, would go here too.
            if (!goes_on_past(*read.aln, placed)) {
                ++stopped;
                continue;
            }
            ++resolved;
            if (score < min_read_score) {
                ++low_score;
                continue;
            }
            if (name == numeric_limits<uint32_t>::max()) {
                name = read_names().intern(read.aln->name());
            }
            out.push_back(FillPlacement{(uint32_t)(point - fill_points.begin()),
                                        {name, placed.alignment, placed.offset, score,
                                         placed.direction}});
        }
    }
}

vector<FillPlacement> FillAnchors::resolve_window(const SiteReadSource::WindowReads& reads,
                                                  size_t begin, size_t end) {
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
    vector<VisitedNode> visited;
    reads(ranges, [&](const SiteRead& read) {
        nodes_visited(read, nodes, visited);
        place_read(read, visited, begin, end, out);
    });
    return out;
}

void FillAnchors::resolve_window_reads(nid_t first, nid_t last,
                                       const SiteReadSource::WindowReads& reads) {
    // Only a range that is exactly one window holds all of that window's points and none of
    // another's.
    if (window_size == 0 || first % (nid_t)window_size != 0
        || last - first + 1 != (nid_t)window_size) {
        return;
    }
    const nid_t number = first / (nid_t)window_size;
    const auto found = std::ranges::lower_bound(windows, number, {}, &PointWindow::number);
    if (found == windows.end() || found->number != number) {
        return;
    }
    const auto w = (size_t)(found - windows.begin());
    if (window_claimed[w].exchange(true)) {
        return;
    }
    window_placements[w] = resolve_window(reads, found->begin, found->end);
    ++windows_observed;
}

void FillAnchors::resolve_reads(const SiteReadSource& reads) {
    placed_reads.clear();
    if (fill_points.empty()) {
        return;
    }
    const SiteReadSource::WindowReads fetch = [&reads](
        const vector<pair<nid_t, nid_t>>& ranges, const function<void(const SiteRead&)>& iteratee) {
        reads.for_each_read_uncounted(ranges, iteratee);
    };
    // The placements of each window are kept apart, so that their order depends neither on the
    // threads nor on which windows `resolve_window_reads` took.
    std::atomic<size_t> fetched{0};
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t w = 0; w < windows.size(); ++w) {
        if (!window_claimed[w].exchange(true)) {
            window_placements[w] = resolve_window(fetch, windows[w].begin, windows[w].end);
            ++fetched;
        }
    }
    windows_fetched = fetched.load();
    size_t total = 0;
    for (const auto& placements : window_placements) {
        total += placements.size();
    }
    placed_reads.reserve(total);
    for (auto& placements : window_placements) {
        placed_reads.insert(placed_reads.end(), placements.begin(), placements.end());
        vector<FillPlacement>().swap(placements);
    }
    std::ranges::stable_sort(placed_reads, {}, &FillPlacement::point);
}

namespace {

/// Set `bits` on a point's mark from a site at `depth` among the point's sites, where no site
/// deeper than it has marked the point.
void set_mark(std::atomic<uint16_t>& mark, uint16_t depth, uint8_t bits) {
    uint16_t seen = mark.load(std::memory_order_relaxed);
    while (true) {
        const auto seen_depth = (uint16_t)(seen >> 8);
        if (seen_depth != 0 && seen_depth < depth) {
            return;
        }
        // A site at the same depth adds its bits, and a deeper one replaces the bits of the site
        // around it.
        const auto wanted = (uint16_t)(seen_depth == depth ? (seen | bits) : ((depth << 8) | bits));
        if (wanted == seen
            || mark.compare_exchange_weak(seen, wanted, std::memory_order_relaxed)) {
            return;
        }
    }
}

}

void FillAnchors::mark_walk(const net_handle_t& site, const Traversal& walk, uint8_t bits) {
    if (fill_points.empty() || sites == nullptr) {
        return;
    }
    const net_handle_t canonical_site = sites->canonical(site);
    const auto found = std::ranges::lower_bound(holding_sites, canonical_site, net_before);
    if (found == holding_sites.end() || *found != canonical_site) {
        // The site holds no point.
        return;
    }
    const auto site_index = (uint32_t)(found - holding_sites.begin());
    site_marked[site_index].store(1, std::memory_order_relaxed);
    for (const handle_t& h : walk) {
        const nid_t id = graph->get_id(h);
        for (auto point = std::ranges::lower_bound(fill_points, FillPoint{id, 0});
             point != fill_points.end() && point->node == id; ++point) {
            const auto i = (size_t)(point - fill_points.begin());
            // The walk visits the site's boundary nodes, which the site does not hold, so look
            // for it among the point's sites.
            for (uint32_t k = holder_start[i]; k < holder_start[i + 1]; ++k) {
                if (holders[k] == site_index) {
                    set_mark(marks[i], (uint16_t)(k - holder_start[i] + 1), bits);
                    break;
                }
            }
        }
    }
}

FillSlot FillAnchors::slot_of(size_t point) const {
    const uint32_t first = holder_start[point];
    const uint32_t past = holder_start[point + 1];
    if (first == past) {
        return classify_fill_point(true, 0);
    }
    const uint16_t mark = marks[point].load(std::memory_order_relaxed);
    for (uint32_t k = first; k < past; ++k) {
        if (site_marked[holders[k]].load(std::memory_order_relaxed) != 0) {
            // This is the innermost site with a called genotype. Its walks may not cross the
            // point, and then a site around it that marked the point does not count.
            const bool its_bits = (mark >> 8) == k - first + 1;
            return classify_fill_point(false, its_bits ? (uint8_t)(mark & 0xff) : 0);
        }
    }
    return classify_fill_point(false, 0);
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
            anchor.reads.push_back(placed_reads[k].row);
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
    out << "[vg call] fill anchors: reads placed from " << windows_observed.load() << " windows of "
        << window_size << " node IDs as the direct pass fetched them, and from " << windows_fetched
        << " windows fetched after it" << endl;
    const AnchorPinCounters& pins = read_counters.pins;
    out << "[vg call] fill anchors: " << resolved.load() << " read placements resolved, "
        << pins.verified.load() << " checked against the graph's base and "
        << pins.verify_failed.load() << " failing the check; refused: "
        << not_reached.load() << " reads on the point's node end before the point or start "
        << "after it, " << pins.repeat_visit.load() << " visit the node more than once, "
        << pins.unaligned_base.load() << " have no read base aligned to the node base before the "
        << "point, " << stopped.load() << " stop at the point or go on past it only in a soft "
        << "clip; " << low_score.load() << " placements below the minimum read score" << endl;
    out << "[vg call] fill anchors: points with reads: " << kind_counts[0] << " fill-shared, "
        << kind_counts[1] << " fill-strand, " << kind_counts[2] << " fill-off-call" << endl;
}

}
