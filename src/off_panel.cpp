#include "off_panel.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <map>

#include <omp.h>

#include "alignment.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

double EditCandidate::fraction() const {
    const uint32_t covering = alt_fragments + ref_fragments + other_fragments;
    return covering == 0 ? 0.0 : (double)alt_fragments / (double)covering;
}

namespace {

/// One read base on a node base: `base` along the node's forward strand, or 0 where the read
/// matches the node.
struct Observation {
    nid_t node;
    uint32_t offset;
    char base;
    uint64_t fragment;
    int32_t mapq;

    bool operator<(const Observation& other) const {
        return tie(node, offset, base, fragment) < tie(other.node, other.offset, other.base,
                                                        other.fragment);
    }
};

/// Whether every node the read's path visits is in the graph. A read aligned to a larger graph
/// can leave this one; such a read cannot be walked, so it adds no off-panel evidence.
bool path_in_graph(const HandleGraph& graph, const Alignment& aln) {
    for (const Mapping& m : aln.path().mapping()) {
        if (!graph.has_node(m.position().node_id())) {
            return false;
        }
    }
    return true;
}

/// Walk the read's aligned bases on the nodes `wanted` accepts, each of quality at least `min_q`,
/// giving each one's node, forward offset and forward base (0 for a match).
void for_each_aligned_base(const HandleGraph& graph, const Alignment& aln, int min_q,
                           const function<bool(nid_t)>& wanted,
                           const function<void(nid_t, uint32_t, char)>& visit) {
    const string& qual = aln.quality();
    size_t read_pos = 0;
    for (const Mapping& m : aln.path().mapping()) {
        const nid_t node = m.position().node_id();
        const bool take = wanted(node) && graph.has_node(node);
        const bool reverse = m.position().is_reverse();
        const size_t node_len = take ? graph.get_length(graph.get_handle(node)) : 0;
        size_t node_pos = (size_t)m.position().offset();
        for (const Edit& e : m.edit()) {
            const size_t from = (size_t)e.from_length();
            const size_t to = (size_t)e.to_length();
            if (take && from == to && from > 0) {
                for (size_t i = 0; i < from; ++i) {
                    const size_t r = read_pos + i;
                    if (!qual.empty() && (r >= qual.size() || (int)(uint8_t)qual[r] < min_q)) {
                        continue;
                    }
                    const size_t visited = node_pos + i;
                    if (visited >= node_len) {
                        continue;
                    }
                    const uint32_t forward = (uint32_t)(reverse ? node_len - 1 - visited : visited);
                    char base = 0;
                    if (!e.sequence().empty()) {
                        base = (char)toupper(e.sequence()[i]);
                        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                            continue;
                        }
                        if (reverse) {
                            base = reverse_complement(base);
                        }
                    }
                    visit(node, forward, base);
                }
            }
            node_pos += from;
            read_pos += to;
        }
    }
}

/// Whether a node beside `node`, sharing a neighbour with it, spells the node with the base at
/// `offset` replaced by `alt`.
bool spelled_by_sibling(const HandleGraph& graph, nid_t node, uint32_t offset, char alt) {
    const handle_t handle = graph.get_handle(node);
    string wanted = graph.get_sequence(handle);
    wanted[offset] = alt;
    bool found = false;
    for (bool go_left : {true, false}) {
        graph.follow_edges(handle, go_left, [&](const handle_t& neighbour) {
            graph.follow_edges(neighbour, !go_left, [&](const handle_t& sibling) {
                if (graph.get_id(sibling) != node && graph.get_length(sibling) == wanted.size()
                    && graph.get_sequence(graph.get_is_reverse(sibling) ? graph.flip(sibling)
                                                                        : sibling) == wanted) {
                    found = true;
                }
                return !found;
            });
            return !found;
        });
    }
    return found;
}

}

namespace {

vector<EditCandidate> count_snv_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params) {

    vector<EditCandidate> kept;
    std::hash<string> name_hash;

    // The substitutions, by fragment.
    vector<Observation> alts;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(graph, aln, params.min_base_quality, owns,
                              [&](nid_t node, uint32_t offset, char base) {
                                  if (base != 0) {
                                      alts.push_back(Observation{node, offset, base, fragment,
                                                                 aln.mapping_quality()});
                                  }
                              });
    });
    if (alts.empty()) {
        return kept;
    }
    std::sort(alts.begin(), alts.end());

    // The (node, offset) of each substitution that enough fragments carry, before the fragments
    // covering the base are counted.
    vector<pair<nid_t, uint32_t>> positions;
    unordered_set<nid_t> position_nodes;
    for (size_t i = 0; i < alts.size();) {
        size_t j = i;
        size_t fragments = 0;
        for (; j < alts.size() && alts[j].node == alts[i].node && alts[j].offset == alts[i].offset
               && alts[j].base == alts[i].base;
             ++j) {
            if (j == i || alts[j].fragment != alts[j - 1].fragment) {
                ++fragments;
            }
        }
        if (fragments >= params.min_fragments) {
            positions.emplace_back(alts[i].node, alts[i].offset);
            position_nodes.insert(alts[i].node);
        }
        i = j;
    }
    if (positions.empty()) {
        return kept;
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    // Every base the reads put at those positions, by fragment.
    vector<Observation> covering;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(
            graph, aln, params.min_base_quality,
            [&](nid_t node) { return owns(node) && position_nodes.count(node); },
            [&](nid_t node, uint32_t offset, char base) {
                if (std::binary_search(positions.begin(), positions.end(),
                                       make_pair(node, offset))) {
                    covering.push_back(Observation{node, offset, base, fragment, 0});
                }
            });
    });
    std::sort(covering.begin(), covering.end(),
              [](const Observation& a, const Observation& b) {
                  return tie(a.node, a.offset, a.fragment, a.base)
                         < tie(b.node, b.offset, b.fragment, b.base);
              });

    for (size_t i = 0; i < covering.size();) {
        // One position.
        size_t end = i;
        while (end < covering.size() && covering[end].node == covering[i].node
               && covering[end].offset == covering[i].offset) {
            ++end;
        }
        // Each fragment's one base, or 1 for a fragment whose reads disagree.
        vector<char> fragment_base;
        for (size_t f = i; f < end;) {
            size_t g = f;
            bool agree = true;
            while (g < end && covering[g].fragment == covering[f].fragment) {
                agree = agree && covering[g].base == covering[f].base;
                ++g;
            }
            fragment_base.push_back(agree ? covering[f].base : (char)1);
            f = g;
        }
        const nid_t node = covering[i].node;
        const uint32_t offset = covering[i].offset;
        const char ref = (char)toupper(graph.get_base(graph.get_handle(node), offset));
        for (char alt : {'A', 'C', 'G', 'T'}) {
            if (alt == ref) {
                continue;
            }
            EditCandidate c;
            c.node = node;
            c.offset = offset;
            c.ref = ref;
            c.alt = alt;
            for (char b : fragment_base) {
                if (b == 1) {
                    ++c.discordant_fragments;
                } else if (b == 0) {
                    ++c.ref_fragments;
                } else if (b == alt) {
                    ++c.alt_fragments;
                } else {
                    ++c.other_fragments;
                }
            }
            if (c.alt_fragments < params.min_fragments || c.fraction() < params.min_fraction) {
                continue;
            }
            if (spelled_by_sibling(graph, node, offset, alt)) {
                continue;
            }
            // The mean mapping quality of the reads carrying the ALT.
            Observation probe{node, offset, alt, 0, 0};
            double total = 0.0;
            size_t n = 0;
            for (auto it = std::lower_bound(alts.begin(), alts.end(), probe);
                 it != alts.end() && it->node == node && it->offset == offset && it->base == alt;
                 ++it) {
                total += it->mapq;
                ++n;
            }
            c.mean_alt_mapq = n == 0 ? 0.0f : (float)(total / (double)n);
            kept.push_back(c);
        }
        i = end;
    }
    std::sort(kept.begin(), kept.end());
    return kept;
}

/// A read's walk: its mappings, read so that the first reads its node forward, the walk base each
/// mapping's node starts at, and its alignment along the walk, indels left-aligned.
class ReadWalk {
public:
    explicit ReadWalk(const HandleGraph& graph) : graph(graph) {}

    /// Load the read's walk. Returns false, loading nothing, for a read whose path leaves the
    /// graph.
    bool load(const Alignment& aln) {
        if (!path_in_graph(graph, aln)) {
            return false;
        }
        const Alignment* use = &aln;
        if (aln.path().mapping_size() > 0 && aln.path().mapping(0).position().is_reverse()) {
            flipped = reverse_complement_alignment(aln, [&](nid_t id) {
                return (int64_t)graph.get_length(graph.get_handle(id));
            });
            use = &flipped;
        }
        read = use;
        mappings.clear();
        handles.clear();
        node_start.clear();
        int64_t at = 0;
        for (const Mapping& m : use->path().mapping()) {
            const handle_t h = graph.get_handle(m.position().node_id(), m.position().is_reverse());
            mappings.push_back(&m);
            handles.push_back(h);
            node_start.push_back(at);
            at += (int64_t)graph.get_length(h);
        }
        node_start.push_back(at);
        walked = walk_alignment(mappings, node_start, 0, use->sequence());
        left_align(walked, base_fn);
        return true;
    }

    /// The mapping holding walk base `x`.
    size_t mapping_at(int64_t x) const {
        return (size_t)(std::upper_bound(node_start.begin(), node_start.end(), x)
                        - node_start.begin()) - 1;
    }

    char base(int64_t x) const {
        const size_t i = mapping_at(x);
        return (char)toupper(graph.get_base(handles[i], (size_t)(x - node_start[i])));
    }

    /// A hash of the nodes, with their orientations, of mappings `first` to `last`, or 0 if any
    /// reads its node backward.
    uint64_t walk_hash(size_t first, size_t last) const {
        uint64_t h = 1469598103934665603ULL;
        for (size_t i = first; i <= last; ++i) {
            if (graph.get_is_reverse(handles[i])) {
                return 0;
            }
            h = (h ^ (uint64_t)graph.get_id(handles[i])) * 1099511628211ULL;
        }
        return h == 0 ? 1 : h;
    }

    /// Whether the read base on walk base `x`, if one is, has quality at least `min_q`.
    bool quality_ok(int64_t x, int min_q) const {
        const string& qual = read->quality();
        if (qual.empty()) {
            return true;
        }
        const int64_t r = walked.read_offset_at(x);
        return r < 0 || ((size_t)r < qual.size() && (int)(uint8_t)qual[r] >= min_q);
    }

    const Alignment* read = nullptr;
    vector<const Mapping*> mappings;
    vector<handle_t> handles;
    vector<int64_t> node_start;
    WalkAlignment walked;
    const function<char(int64_t)> base_fn = [this](int64_t x) { return base(x); };

private:
    const HandleGraph& graph;
    Alignment flipped;
};

bool has_indel(const Alignment& aln) {
    for (const Mapping& m : aln.path().mapping()) {
        for (const Edit& e : m.edit()) {
            if (e.from_length() != e.to_length()) {
                return true;
            }
        }
    }
    return false;
}

/// One read's indel, keyed by its anchor, with the walk to the end of its tract where the read
/// covers that far (`span` 0 where it does not).
struct IndelObservationRow {
    nid_t node;
    uint32_t offset;
    char kind;
    uint32_t length;
    uint64_t sequence_hash;
    uint64_t fragment;
    uint32_t span;
    uint64_t walk_hash;

    bool same_key(const IndelObservationRow& o) const {
        return node == o.node && offset == o.offset && kind == o.kind && length == o.length
               && sequence_hash == o.sequence_hash;
    }
    bool operator<(const IndelObservationRow& o) const {
        return tie(node, offset, kind, length, sequence_hash, fragment, span, walk_hash)
               < tie(o.node, o.offset, o.kind, o.length, o.sequence_hash, o.fragment, o.span,
                     o.walk_hash);
    }
};

/// What one read says about one key, by its fragment.
struct Seen {
    size_t key;
    uint64_t fragment;
    uint8_t observed;
    bool operator<(const Seen& o) const {
        return tie(key, fragment, observed) < tie(o.key, o.fragment, o.observed);
    }
};

/// Count the fragments of `seen[i, end)`, sorted and all of one key, into `c`: each fragment once,
/// as the class all its reads agree on, or as discordant.
void count_fragments(const vector<Seen>& seen, size_t i, size_t end, EditCandidate& c) {
    for (size_t f = i; f < end;) {
        size_t g = f;
        bool agree = true;
        while (g < end && seen[g].fragment == seen[f].fragment) {
            agree = agree && seen[g].observed == seen[f].observed;
            ++g;
        }
        if (!agree) {
            ++c.discordant_fragments;
        } else if (seen[f].observed == (uint8_t)IndelObservation::ALT) {
            ++c.alt_fragments;
        } else if (seen[f].observed == (uint8_t)IndelObservation::REF) {
            ++c.ref_fragments;
        } else if (seen[f].observed == (uint8_t)IndelObservation::OTHER) {
            ++c.other_fragments;
        }
        f = g;
    }
}

/// Whether a walk of the graph from `walk`'s first handle to its last other than `walk` spells
/// `alt`.
bool spelled_by_other_walk(const HandleGraph& graph, const vector<handle_t>& walk,
                           const string& alt) {
    vector<pair<handle_t, size_t>> stack{{walk.front(), 0}};
    size_t budget = 4096;
    while (!stack.empty() && budget-- > 0) {
        const pair<handle_t, size_t> at = stack.back();
        stack.pop_back();
        const string seq = graph.get_sequence(at.first);
        if (at.second + seq.size() > alt.size() || alt.compare(at.second, seq.size(), seq) != 0) {
            continue;
        }
        const size_t next = at.second + seq.size();
        if (next == alt.size()) {
            if (at.first == walk.back()) {
                return true;
            }
            continue;
        }
        graph.follow_edges(at.first, false, [&](const handle_t& n) {
            stack.emplace_back(n, next);
        });
    }
    return false;
}

/// Whether another walk of the graph spells the indel candidate `c`'s change, as a panel allele.
bool spelled_by_panel(const HandleGraph& graph, const EditCandidate& c) {
    vector<handle_t> handles;
    string alt;
    for (nid_t id : c.walk) {
        handles.push_back(graph.get_handle(id));
        alt += graph.get_sequence(handles.back());
    }
    if (c.kind == 'D') {
        alt.erase(c.offset + 1, c.sequence.size());
    } else {
        alt.insert(c.offset + 1, c.sequence);
    }
    return spelled_by_other_walk(graph, handles, alt);
}

/// One read's indel of `LARGE_INDEL_MIN_LENGTH` bases or more, anchored on a node the count owns,
/// with the walk from its anchor's node to that of the first base after its tract where the read
/// reaches that far, forward (`walk` empty where it does not).
struct LargeRow {
    char kind;
    nid_t node;
    uint32_t offset;
    uint32_t length;
    uint64_t fragment;
    string inserted;
    vector<nid_t> walk;
    uint32_t span;

    bool operator<(const LargeRow& o) const {
        return tie(kind, node, offset, length, fragment, inserted, span, walk)
               < tie(o.kind, o.node, o.offset, o.length, o.fragment, o.inserted, o.span, o.walk);
    }
};

/// One read's clip of `min_clip` bases or more: 'R' where it follows the read's last aligned base,
/// at `node` and `offset`, and 'L' where it precedes its first; with the `min_clip` clipped bases
/// nearest that base.
struct ClipRow {
    char side;
    nid_t node;
    uint32_t offset;
    string seed;
    uint64_t fragment;

    bool operator<(const ClipRow& o) const {
        return tie(side, node, offset, seed, fragment)
               < tie(o.side, o.node, o.offset, o.seed, o.fragment);
    }
};

string upper(string bases) {
    for (char& b : bases) {
        b = (char)toupper(b);
    }
    return bases;
}

/// File `indel`, one of `walk`'s of `LARGE_INDEL_MIN_LENGTH` bases or more, in `rows` if its anchor
/// is on a node `owns` accepts, read forward.
void add_large_row(const HandleGraph& graph, const ReadWalk& walk, const WalkIndel& indel,
                   uint64_t fragment, const function<bool(nid_t)>& owns, vector<LargeRow>& rows) {
    const int64_t anchor = indel.at - 1;
    const size_t i = walk.mapping_at(anchor);
    const handle_t h = walk.handles[i];
    if (graph.get_is_reverse(h) || !owns(graph.get_id(h))) {
        return;
    }
    LargeRow row{indel.deleted > 0 ? 'D' : 'I', graph.get_id(h),
                 (uint32_t)(anchor - walk.node_start[i]), indel.length(), fragment,
                 indel.inserted, {}, 0};
    const int64_t after = tract_end(indel, walk.walked.end, walk.base_fn);
    if (after < walk.walked.end && walk.walk_hash(i, walk.mapping_at(after)) != 0) {
        for (size_t m = i; m <= walk.mapping_at(after); ++m) {
            row.walk.push_back(graph.get_id(walk.handles[m]));
        }
        row.span = (uint32_t)(after - anchor);
    }
    rows.push_back(std::move(row));
}

/// File `walk`'s clips of `min_clip` bases or more in `rows`, where the aligned base beside each is
/// on a node read forward.
void add_clip_rows(const HandleGraph& graph, const ReadWalk& walk, uint64_t fragment,
                   const LargeIndelParams& params, vector<ClipRow>& rows) {
    const WalkAlignment& w = walk.walked;
    if (w.blocks.empty()) {
        return;
    }
    const string& sequence = walk.read->sequence();
    if (w.clip_end >= params.min_clip) {
        const size_t i = walk.mapping_at(w.end - 1);
        const size_t from = w.blocks.back().read + w.blocks.back().length;
        if (!graph.get_is_reverse(walk.handles[i]) && from + params.min_clip <= sequence.size()) {
            rows.push_back(ClipRow{'R', graph.get_id(walk.handles[i]),
                                   (uint32_t)(w.end - 1 - walk.node_start[i]),
                                   upper(sequence.substr(from, params.min_clip)), fragment});
        }
    }
    if (w.clip_begin >= params.min_clip) {
        const size_t i = walk.mapping_at(w.begin);
        const size_t to = w.blocks.front().read;
        if (!graph.get_is_reverse(walk.handles[i]) && to >= params.min_clip
            && to <= sequence.size()) {
            rows.push_back(ClipRow{'L', graph.get_id(walk.handles[i]),
                                   (uint32_t)(w.begin - walk.node_start[i]),
                                   upper(sequence.substr(to - params.min_clip, params.min_clip)),
                                   fragment});
        }
    }
}

/// The indels among `kept`'s candidates: counted as `count_snv_edits` counts substitutions, with
/// the observation classes of `observe_indel`.
vector<EditCandidate> count_indel_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params,
    vector<LargeRow>* large_rows, vector<ClipRow>* clip_rows) {

    vector<EditCandidate> kept;
    std::hash<string> name_hash;
    std::hash<string> sequence_hasher;
    ReadWalk walk(graph);

    // Every read's indels, anchored on a node this count owns, wherever the aligner put them.
    vector<IndelObservationRow> rows;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq || !has_indel(aln)) {
            return;
        }
        if (!walk.load(aln)) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for (const WalkIndel& indel : walk.walked.indels) {
            const int64_t anchor = indel.at - 1;
            if (large_rows != nullptr && anchor >= walk.walked.begin
                && indel.length() >= LARGE_INDEL_MIN_LENGTH
                && indel.length() <= params.large.max_length) {
                add_large_row(graph, walk, indel, fragment, owns, *large_rows);
            }
            if (anchor < walk.walked.begin || indel.length() > params.max_indel_length) {
                continue;
            }
            const size_t i = walk.mapping_at(anchor);
            const handle_t h = walk.handles[i];
            if (graph.get_is_reverse(h) || !owns(graph.get_id(h))) {
                continue;
            }
            IndelObservationRow row{graph.get_id(h), (uint32_t)(anchor - walk.node_start[i]),
                                    indel.deleted > 0 ? 'D' : 'I', indel.length(),
                                    indel.deleted > 0 ? 0 : sequence_hasher(indel.inserted),
                                    fragment, 0, 0};
            const int64_t after = tract_end(indel, walk.walked.end, walk.base_fn);
            if (after < walk.walked.end) {
                const uint64_t h_walk = walk.walk_hash(i, walk.mapping_at(after));
                if (h_walk != 0) {
                    row.span = (uint32_t)(after - anchor);
                    row.walk_hash = h_walk;
                }
            }
            rows.push_back(row);
        }
        if (clip_rows != nullptr) {
            add_clip_rows(graph, walk, fragment, params.large, *clip_rows);
        }
    });
    if (rows.empty()) {
        return kept;
    }
    std::sort(rows.begin(), rows.end());

    // The keys that enough fragments carry, each with the walk most of its reads take to the end
    // of its tract.
    vector<IndelObservationRow> keys;
    for (size_t i = 0; i < rows.size();) {
        size_t j = i;
        size_t fragments = 0;
        map<pair<uint32_t, uint64_t>, size_t> walks;
        for (; j < rows.size() && rows[j].same_key(rows[i]); ++j) {
            if (j == i || rows[j].fragment != rows[j - 1].fragment) {
                ++fragments;
            }
            if (rows[j].span > 0) {
                ++walks[make_pair(rows[j].span, rows[j].walk_hash)];
            }
        }
        if (fragments >= params.min_fragments && !walks.empty()) {
            auto best = walks.begin();
            for (auto it = walks.begin(); it != walks.end(); ++it) {
                if (it->second > best->second) {
                    best = it;
                }
            }
            IndelObservationRow key = rows[i];
            key.span = best->first.first;
            key.walk_hash = best->first.second;
            keys.push_back(key);
        }
        i = j;
    }
    if (keys.empty()) {
        return kept;
    }
    unordered_set<nid_t> anchor_nodes;
    for (const IndelObservationRow& key : keys) {
        anchor_nodes.insert(key.node);
    }

    // What each read covering a key's anchor through the end of its tract, on the key's walk, says
    // about it, by fragment.
    vector<Seen> seen;
    vector<EditCandidate> found(keys.size());
    vector<double> mapq_total(keys.size(), 0.0);
    vector<size_t> mapq_count(keys.size(), 0);
    vector<const WalkIndel*> touching;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        bool near = false;
        for (const Mapping& m : aln.path().mapping()) {
            if (anchor_nodes.count(m.position().node_id())) {
                near = true;
                break;
            }
        }
        if (!near) {
            return;
        }
        if (!walk.load(aln)) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for (size_t i = 0; i + 1 < walk.node_start.size(); ++i) {
            const handle_t h = walk.handles[i];
            const nid_t id = graph.get_id(h);
            if (graph.get_is_reverse(h) || !anchor_nodes.count(id)) {
                continue;
            }
            IndelObservationRow probe{};
            probe.node = id;
            auto k = std::lower_bound(keys.begin(), keys.end(), probe,
                                      [](const IndelObservationRow& a,
                                         const IndelObservationRow& b) {
                                          return a.node < b.node;
                                      });
            for (; k != keys.end() && k->node == id; ++k) {
                const size_t key = k - keys.begin();
                const int64_t anchor = walk.node_start[i] + k->offset;
                const int64_t after = anchor + k->span;
                if (anchor < walk.walked.begin || after >= walk.walked.end
                    || walk.walk_hash(i, walk.mapping_at(after)) != k->walk_hash
                    || !walk.quality_ok(anchor, params.min_base_quality)
                    || !walk.quality_ok(after, params.min_base_quality)) {
                    continue;
                }
                WalkIndel wanted;
                wanted.at = anchor + 1;
                if (k->kind == 'D') {
                    wanted.deleted = k->length;
                }
                IndelObservation observed = observe_indel(walk.walked, wanted, after, &touching);
                if (observed == IndelObservation::OTHER && touching.size() == 1
                    && k->kind == 'I' && touching[0]->at == wanted.at
                    && touching[0]->deleted == 0 && touching[0]->inserted.size() == k->length
                    && sequence_hasher(touching[0]->inserted) == k->sequence_hash) {
                    observed = IndelObservation::ALT;
                }
                seen.push_back(Seen{key, fragment, (uint8_t)observed});
                if (observed != IndelObservation::ALT) {
                    continue;
                }
                mapq_total[key] += aln.mapping_quality();
                ++mapq_count[key];
                EditCandidate& c = found[key];
                if (c.walk.empty()) {
                    // The key fixes all of these, whichever read gives them.
                    c.node = id;
                    c.offset = k->offset;
                    c.kind = k->kind;
                    c.alt = 0;
                    c.ref = walk.base(anchor);
                    c.span = k->span;
                    if (k->kind == 'D') {
                        for (int64_t x = anchor + 1; x <= anchor + (int64_t)k->length; ++x) {
                            c.sequence.push_back(walk.base(x));
                        }
                    } else {
                        c.sequence = touching[0]->inserted;
                    }
                    for (size_t m = i; m <= walk.mapping_at(after); ++m) {
                        c.walk.push_back(graph.get_id(walk.handles[m]));
                    }
                    const char first = c.sequence[0];
                    for (int64_t x = anchor + 1; x < walk.walked.end && walk.base(x) == first;
                         ++x) {
                        ++c.homopolymer;
                    }
                }
            }
        }
    });
    std::sort(seen.begin(), seen.end());
    for (size_t i = 0; i < seen.size();) {
        const size_t key = seen[i].key;
        size_t end = i;
        while (end < seen.size() && seen[end].key == key) {
            ++end;
        }
        EditCandidate& c = found[key];
        count_fragments(seen, i, end, c);
        i = end;
        if (c.walk.empty() || c.alt_fragments < params.min_fragments
            || c.fraction() < params.min_fraction) {
            continue;
        }
        if (spelled_by_panel(graph, c)) {
            continue;
        }
        c.mean_alt_mapq =
            mapq_count[key] == 0 ? 0.0f : (float)(mapq_total[key] / (double)mapq_count[key]);
        kept.push_back(std::move(c));
    }
    return kept;
}

/// A large indel to count reads for: its anchor, kind, length and ALT, and the walk from its
/// anchor's node to that of the first base after its tract, with where each node starts along it
/// and its bases.
struct LargeKey {
    char kind = 'D';
    nid_t node = 0;
    uint32_t offset = 0;
    uint32_t length = 0;
    string inserted;
    vector<nid_t> walk;
    vector<int64_t> starts;
    string spelling;
    uint32_t span = 0;
    /// The fragments that found it, and whether only clips did.
    size_t support = 0;
    bool clips_only = false;
};

/// The nodes met from `from`, read forward, going forward (or backward) along the only edge out of
/// each node into a node with no other edge in, until they hold `limit` bases; in walk order.
vector<handle_t> unbranched_walk(const HandleGraph& graph, handle_t from, bool forward,
                                 size_t limit) {
    vector<handle_t> out{from};
    size_t total = graph.get_length(from);
    handle_t at = from;
    while (total < limit) {
        size_t n = 0;
        handle_t next = at;
        graph.follow_edges(at, !forward, [&](const handle_t& h) {
            next = h;
            return ++n < 2;
        });
        if (n != 1 || graph.get_is_reverse(next) || graph.get_id(next) == graph.get_id(from)) {
            break;
        }
        size_t in = 0;
        graph.follow_edges(next, forward, [&](const handle_t&) { return ++in < 2; });
        if (in != 1) {
            break;
        }
        out.push_back(next);
        total += graph.get_length(next);
        at = next;
    }
    if (!forward) {
        std::reverse(out.begin(), out.end());
    }
    return out;
}

/// The key for `indel`, given in the coordinates of the walk `nodes` read forward, once
/// left-aligned along it; false if its tract runs off the walk's end.
bool make_large_key(const HandleGraph& graph, const vector<nid_t>& nodes, WalkIndel indel,
                    size_t support, bool clips_only, LargeKey& out) {
    string spelling;
    vector<int64_t> starts;
    for (nid_t id : nodes) {
        starts.push_back((int64_t)spelling.size());
        spelling += graph.get_sequence(graph.get_handle(id));
    }
    spelling = upper(std::move(spelling));
    const function<char(int64_t)> base = [&](int64_t x) { return spelling[x]; };
    left_align(indel, 1, base);
    const int64_t after = tract_end(indel, (int64_t)spelling.size(), base);
    if (after >= (int64_t)spelling.size() || indel.at < 1) {
        return false;
    }
    const int64_t anchor = indel.at - 1;
    auto node_of = [&](int64_t x) {
        return (size_t)(std::upper_bound(starts.begin(), starts.end(), x) - starts.begin()) - 1;
    };
    const size_t first = node_of(anchor);
    const size_t last = node_of(after);
    out = LargeKey();
    out.kind = indel.deleted > 0 ? 'D' : 'I';
    out.node = nodes[first];
    out.offset = (uint32_t)(anchor - starts[first]);
    out.length = indel.length();
    out.inserted = indel.inserted;
    out.walk.assign(nodes.begin() + first, nodes.begin() + last + 1);
    for (size_t i = first; i <= last; ++i) {
        out.starts.push_back(starts[i] - starts[first]);
    }
    const int64_t end = last + 1 < starts.size() ? starts[last + 1] : (int64_t)spelling.size();
    out.spelling = spelling.substr(starts[first], end - starts[first]);
    out.span = (uint32_t)(after - anchor);
    out.support = support;
    out.clips_only = clips_only;
    return true;
}

/// The sorted 8-mer codes of `bases`, for a cheap distance between two reads' insertions.
vector<uint32_t> kmer_profile(const string& bases) {
    vector<uint32_t> out;
    uint32_t code = 0;
    size_t valid = 0;
    for (char b : bases) {
        const char* at = strchr("ACGT", b);
        if (b == 0 || at == nullptr) {
            valid = 0;
            continue;
        }
        code = ((code << 2) | (uint32_t)(at - "ACGT")) & 0xFFFF;
        if (++valid >= 8) {
            out.push_back(code);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

/// The size of the multiset difference of two sorted profiles, both ways.
size_t profile_distance(const vector<uint32_t>& a, const vector<uint32_t>& b) {
    size_t i = 0, j = 0, shared = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            ++shared;
            ++i;
            ++j;
        }
    }
    return a.size() + b.size() - 2 * shared;
}

/// Where `key`'s anchor lies along `other`'s walk, if `other`'s walk holds its anchor's node, or
/// `other`'s anchor along `key`'s, as a difference `key` minus `other`; false if neither does.
bool anchor_difference(const LargeKey& key, const LargeKey& other, int64_t& difference) {
    for (size_t j = 0; j < other.walk.size(); ++j) {
        if (other.walk[j] == key.node) {
            difference = other.starts[j] + key.offset - (int64_t)other.offset;
            return true;
        }
    }
    for (size_t j = 0; j < key.walk.size(); ++j) {
        if (key.walk[j] == other.node) {
            difference = (int64_t)key.offset - (key.starts[j] + other.offset);
            return true;
        }
    }
    return false;
}

/// The large indels among the reads `for_each_read` visits, from the rows the indel count filed:
/// each read's own large indels (`large_rows`), clustered by anchor within `slop` and length within
/// `length_ratio`, the medoid read's taken as the ALT; and deletions found from clips
/// (`clip_rows`), whose clipped bases the graph spells further along. Each is then counted as
/// `observe_large_indel` matches the reads to it.
vector<EditCandidate> count_large_indel_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params,
    vector<LargeRow> rows, vector<ClipRow> clips) {
    const LargeIndelParams& large = params.large;
    vector<LargeKey> keys;
    auto fragments_of = [](vector<uint64_t> f) {
        std::sort(f.begin(), f.end());
        return (size_t)(std::unique(f.begin(), f.end()) - f.begin());
    };

    // The reads' own large indels: on one node, by anchor, then by length.
    std::sort(rows.begin(), rows.end());
    for (size_t i = 0; i < rows.size();) {
        size_t j = i + 1;
        while (j < rows.size() && rows[j].kind == rows[i].kind && rows[j].node == rows[i].node
               && (int64_t)rows[j].offset - (int64_t)rows[j - 1].offset <= large.slop) {
            ++j;
        }
        vector<size_t> by_length;
        for (size_t k = i; k < j; ++k) {
            by_length.push_back(k);
        }
        std::stable_sort(by_length.begin(), by_length.end(), [&](size_t a, size_t b) {
            return rows[a].length < rows[b].length;
        });
        for (size_t a = 0; a < by_length.size();) {
            size_t b = a + 1;
            while (b < by_length.size()
                   && (double)rows[by_length[b]].length
                          <= (double)rows[by_length[b - 1]].length * (1.0 + large.length_ratio)) {
                ++b;
            }
            vector<uint64_t> fragments;
            vector<size_t> walked;
            for (size_t k = a; k < b; ++k) {
                fragments.push_back(rows[by_length[k]].fragment);
                if (!rows[by_length[k]].walk.empty() && walked.size() < 40) {
                    walked.push_back(by_length[k]);
                }
            }
            const size_t support = fragments_of(std::move(fragments));
            if (support >= params.min_fragments && !walked.empty()) {
                // The medoid: the read whose indel is nearest all the others'.
                vector<vector<uint32_t>> profiles;
                if (rows[i].kind == 'I') {
                    for (size_t k : walked) {
                        profiles.push_back(kmer_profile(rows[k].inserted));
                    }
                }
                size_t best = 0;
                size_t best_total = numeric_limits<size_t>::max();
                for (size_t x = 0; x < walked.size(); ++x) {
                    size_t total = 0;
                    for (size_t y = 0; y < walked.size(); ++y) {
                        const LargeRow& p = rows[walked[x]];
                        const LargeRow& q = rows[walked[y]];
                        total += rows[i].kind == 'I'
                                     ? profile_distance(profiles[x], profiles[y])
                                     : (size_t)(std::abs((int64_t)p.offset - (int64_t)q.offset)
                                                + std::abs((int64_t)p.length - (int64_t)q.length));
                    }
                    if (total < best_total) {
                        best_total = total;
                        best = x;
                    }
                }
                const LargeRow& medoid = rows[walked[best]];
                WalkIndel indel;
                indel.at = (int64_t)medoid.offset + 1;
                if (medoid.kind == 'D') {
                    indel.deleted = medoid.length;
                } else {
                    indel.inserted = medoid.inserted;
                }
                LargeKey key;
                if (make_large_key(graph, medoid.walk, indel, support, false, key)) {
                    keys.push_back(std::move(key));
                }
            }
            a = b;
        }
        i = j;
    }

    // Deletions from clips: each clip position many fragments share, with the clipped bases most
    // of them carry, found again within `max_clip_deletion` along the graph past it.
    std::sort(clips.begin(), clips.end());
    for (size_t i = 0; i < clips.size();) {
        size_t j = i;
        while (j < clips.size() && clips[j].side == clips[i].side && clips[j].node == clips[i].node
               && clips[j].offset == clips[i].offset) {
            ++j;
        }
        size_t best_support = 0;
        size_t best = i;
        for (size_t k = i; k < j;) {
            size_t l = k;
            vector<uint64_t> fragments;
            while (l < j && clips[l].seed == clips[k].seed) {
                fragments.push_back(clips[l].fragment);
                ++l;
            }
            const size_t support = fragments_of(std::move(fragments));
            if (support > best_support) {
                best_support = support;
                best = k;
            }
            k = l;
        }
        const ClipRow& clip = clips[best];
        i = j;
        if (best_support < params.min_fragments) {
            continue;
        }
        if (!graph.has_node(clip.node)) {
            continue;
        }
        const handle_t h = graph.get_handle(clip.node);
        const size_t seed = clip.seed.size();
        const size_t reach = (size_t)large.max_clip_deletion + seed + 1;
        vector<handle_t> handles;
        vector<nid_t> nodes;
        string spelling;
        WalkIndel indel;
        if (clip.side == 'R') {
            handles = unbranched_walk(graph, h, true, clip.offset + reach);
            for (const handle_t& n : handles) {
                nodes.push_back(graph.get_id(n));
                spelling += graph.get_sequence(n);
            }
            spelling = upper(std::move(spelling));
            const size_t from = (size_t)clip.offset + 1 + LARGE_INDEL_MIN_LENGTH;
            const size_t q = spelling.find(clip.seed, from);
            if (q == string::npos || q - clip.offset - 1 > large.max_clip_deletion) {
                continue;
            }
            indel.at = (int64_t)clip.offset + 1;
            indel.deleted = (uint32_t)(q - clip.offset - 1);
        } else {
            handles = unbranched_walk(graph, h, false, reach + graph.get_length(h));
            for (const handle_t& n : handles) {
                nodes.push_back(graph.get_id(n));
                spelling += graph.get_sequence(n);
            }
            spelling = upper(std::move(spelling));
            const int64_t q = (int64_t)(spelling.size() - graph.get_length(h)) + clip.offset;
            const int64_t last_start = q - (int64_t)LARGE_INDEL_MIN_LENGTH - (int64_t)seed;
            if (last_start < 1) {
                continue;
            }
            const size_t s = spelling.rfind(clip.seed, (size_t)last_start);
            if (s == string::npos || s < 1
                || q - (int64_t)(s + seed) > (int64_t)large.max_clip_deletion) {
                continue;
            }
            indel.at = (int64_t)(s + seed);
            indel.deleted = (uint32_t)(q - indel.at);
        }
        LargeKey key;
        if (make_large_key(graph, nodes, indel, best_support, true, key)) {
            keys.push_back(std::move(key));
        }
    }

    // One key per change: the best supported of those matching, reads' own indels first.
    vector<size_t> order(keys.size());
    for (size_t k = 0; k < keys.size(); ++k) {
        order[k] = k;
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return tie(keys[a].clips_only, keys[b].support) < tie(keys[b].clips_only, keys[a].support);
    });
    vector<LargeKey> chosen;
    for (size_t k : order) {
        const LargeKey& key = keys[k];
        bool duplicate = !owns(key.node);
        for (size_t c = 0; c < chosen.size() && !duplicate; ++c) {
            int64_t difference = 0;
            const LargeKey& other = chosen[c];
        // The tract lets either slide right this far.
        const int64_t room = std::max<int64_t>(
            0, (int64_t)other.span - 1 - (other.kind == 'D' ? (int64_t)other.length : 0));
        duplicate = other.kind == key.kind && anchor_difference(key, other, difference)
                    && std::abs(difference) <= large.slop + room
                    && std::abs((double)key.length - (double)other.length)
                           <= large.length_ratio * (double)other.length;
        }
        if (!duplicate) {
            chosen.push_back(key);
        }
    }
    keys = std::move(chosen);
    std::sort(keys.begin(), keys.end(), [](const LargeKey& a, const LargeKey& b) {
        return tie(a.node, a.offset, a.kind, a.length, a.inserted)
               < tie(b.node, b.offset, b.kind, b.length, b.inserted);
    });
    vector<EditCandidate> kept;
    if (keys.empty()) {
        return kept;
    }

    // What each read on a key's walk says about it, by fragment.
    unordered_map<nid_t, vector<pair<size_t, size_t>>> on_node;
    for (size_t k = 0; k < keys.size(); ++k) {
        for (size_t j = 0; j < keys[k].walk.size(); ++j) {
            on_node[keys[k].walk[j]].emplace_back(k, j);
        }
    }
    vector<Seen> seen;
    // The fragments clipped at each end of a deletion, and those with it whole.
    vector<vector<uint64_t>> clipped_left(keys.size()), clipped_right(keys.size());
    vector<size_t> whole(keys.size(), 0);
    vector<double> mapq_total(keys.size(), 0.0);
    vector<size_t> mapq_count(keys.size(), 0);
    ReadWalk walk(graph);
    std::hash<string> name_hash;
    vector<bool> done(keys.size(), false);
    vector<size_t> touched;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        bool near = false;
        for (const Mapping& m : aln.path().mapping()) {
            if (on_node.count(m.position().node_id())) {
                near = true;
                break;
            }
        }
        if (!near) {
            return;
        }
        if (!walk.load(aln)) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for (size_t k : touched) {
            done[k] = false;
        }
        touched.clear();
        for (size_t i = 0; i + 1 < walk.node_start.size(); ++i) {
            const handle_t h = walk.handles[i];
            auto found = on_node.find(graph.get_id(h));
            if (found == on_node.end() || graph.get_is_reverse(h)) {
                continue;
            }
            for (const pair<size_t, size_t>& entry : found->second) {
                const size_t k = entry.first;
                if (done[k]) {
                    continue;
                }
                done[k] = true;
                touched.push_back(k);
                const LargeKey& key = keys[k];
                // The read must take the key's walk wherever the two overlap.
                bool agrees = true;
                for (size_t t = 0; t < key.walk.size() && agrees; ++t) {
                    const int64_t m = (int64_t)i + (int64_t)t - (int64_t)entry.second;
                    if (m >= 0 && m + 1 < (int64_t)walk.node_start.size()) {
                        agrees = walk.handles[m] == graph.get_handle(key.walk[t]);
                    }
                }
                if (!agrees) {
                    continue;
                }
                const int64_t anchor = walk.node_start[i] - key.starts[entry.second] + key.offset;
                WalkIndel wanted;
                wanted.at = anchor + 1;
                if (key.kind == 'D') {
                    wanted.deleted = key.length;
                } else {
                    wanted.inserted = key.inserted;
                }
                int clip = 0;
                const IndelObservation observed =
                    observe_large_indel(walk.walked, wanted, anchor + key.span, large,
                                        walk.walked.begin, walk.walked.end, &clip);
                if (observed == IndelObservation::NONE) {
                    continue;
                }
                seen.push_back(Seen{k, fragment, (uint8_t)observed});
                if (observed != IndelObservation::ALT) {
                    continue;
                }
                if (clip == 1) {
                    clipped_left[k].push_back(fragment);
                } else if (clip == 2) {
                    clipped_right[k].push_back(fragment);
                } else {
                    ++whole[k];
                }
                mapq_total[k] += aln.mapping_quality();
                ++mapq_count[k];
            }
        }
    });
    std::sort(seen.begin(), seen.end());
    for (size_t i = 0; i < seen.size();) {
        const size_t k = seen[i].key;
        size_t end = i;
        while (end < seen.size() && seen[end].key == k) {
            ++end;
        }
        const LargeKey& key = keys[k];
        EditCandidate c;
        count_fragments(seen, i, end, c);
        i = end;
        if (c.alt_fragments < params.min_fragments || c.fraction() < params.min_fraction) {
            continue;
        }
        if (whole[k] == 0
            && (fragments_of(clipped_left[k]) < large.min_fragments_per_end
                || fragments_of(clipped_right[k]) < large.min_fragments_per_end)) {
            continue;
        }
        c.node = key.node;
        c.offset = key.offset;
        c.kind = key.kind;
        c.alt = 0;
        c.ref = key.spelling[key.offset];
        c.span = key.span;
        c.sequence = key.kind == 'D' ? key.spelling.substr(key.offset + 1, key.length)
                                     : key.inserted;
        c.walk = key.walk;
        for (size_t x = key.offset + 1; x < key.spelling.size() && key.spelling[x] == c.sequence[0];
             ++x) {
            ++c.homopolymer;
        }
        if (spelled_by_panel(graph, c)) {
            continue;
        }
        c.mean_alt_mapq = mapq_count[k] == 0 ? 0.0f : (float)(mapq_total[k] / (double)mapq_count[k]);
        kept.push_back(std::move(c));
    }
    return kept;
}

/// Merge `kept`'s touching candidates, as `count_edits` says.
void merge_touching(const HandleGraph& graph,
                    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
                    const EditCountParams& params, vector<EditCandidate>& kept) {
    // The changed bases [from, to) along each admitted candidate's node, by node.
    vector<tuple<nid_t, uint32_t, uint32_t, size_t>> changes;
    for (size_t i = 0; i < kept.size(); ++i) {
        const EditCandidate& c = kept[i];
        if (params.bar->admits(c) && !c.is_large()) {
            const uint32_t from = c.kind == 'S' ? c.offset : c.offset + 1;
            changes.emplace_back(c.node, from, c.kind == 'S' ? from + 1 : c.offset + c.span, i);
        }
    }
    std::sort(changes.begin(), changes.end());
    // Runs of them that touch, each with its bases [lo, hi) and the longest of its walks.
    struct Group {
        vector<size_t> members;
        uint32_t lo;
        uint32_t hi;
        vector<nid_t> walk;
    };
    vector<Group> groups;
    for (size_t i = 0, j = 0; i < changes.size(); i = j) {
        Group g{{}, get<1>(changes[i]), get<2>(changes[i]), {get<0>(changes[i])}};
        bool one_anchor = true;
        for (; j < changes.size() && get<0>(changes[j]) == get<0>(changes[i])
               && get<1>(changes[j]) <= g.hi;
             ++j) {
            const EditCandidate& c = kept[get<3>(changes[j])];
            g.hi = std::max(g.hi, get<2>(changes[j]));
            g.members.push_back(get<3>(changes[j]));
            g.walk = c.walk.size() > g.walk.size() ? c.walk : g.walk;
            one_anchor = one_anchor && c.is_indel() && c.offset == kept[g.members[0]].offset;
        }
        if (g.members.size() > 1 && g.hi > g.lo && !one_anchor) {
            groups.push_back(std::move(g));
        }
    }
    if (groups.empty()) {
        return;
    }

    // What each read following a group's walk spells between the bases on either side of its
    // bases, where it aligns both, by fragment.
    vector<tuple<size_t, uint64_t, string>> spelled;
    std::hash<string> name_hash;
    ReadWalk walk(graph);
    const int min_q = params.min_base_quality;
    unordered_set<nid_t> nodes;
    for (const Group& g : groups) {
        nodes.insert(g.walk[0]);
    }
    for_each_read([&](const Alignment& aln) {
        const auto& maps = aln.path().mapping();
        if (aln.mapping_quality() < params.min_mapq
            || std::none_of(maps.begin(), maps.end(),
                            [&](const Mapping& m) { return nodes.count(m.position().node_id()); })
            || !walk.load(aln)) {
            return;
        }
        for (size_t i = 0; i < walk.handles.size(); ++i) {
            const nid_t id = graph.get_id(walk.handles[i]);
            for (size_t k = std::lower_bound(groups.begin(), groups.end(), id,
                                             [](const Group& g, nid_t n) { return g.walk[0] < n; })
                            - groups.begin();
                 k < groups.size() && groups[k].walk[0] == id; ++k) {
                const Group& g = groups[k];
                bool follows = i + g.walk.size() <= walk.handles.size();
                for (size_t s = 0; follows && s < g.walk.size(); ++s) {
                    follows = walk.handles[i + s] == graph.get_handle(g.walk[s]);
                }
                const int64_t before = walk.node_start[i] + g.lo - 1;
                const int64_t after = walk.node_start[i] + g.hi;
                const int64_t r0 = follows ? walk.walked.read_offset_at(before) : -1;
                const int64_t r1 = follows ? walk.walked.read_offset_at(after) : -1;
                bool good = r0 >= 0 && r1 > r0 && walk.quality_ok(before, min_q)
                            && walk.quality_ok(after, min_q);
                for (size_t m : g.members) {
                    good = good && (kept[m].kind != 'S'
                                    || walk.quality_ok(walk.node_start[i] + kept[m].offset, min_q));
                }
                if (good) {
                    spelled.emplace_back(k, name_hash(aln.name()),
                                         upper(walk.read->sequence().substr(r0 + 1, r1 - r0 - 1)));
                }
            }
        }
    });
    std::sort(spelled.begin(), spelled.end());
    spelled.erase(std::unique(spelled.begin(), spelled.end()), spelled.end());

    vector<EditCandidate> merged;
    for (size_t i = 0, j = 0; i < spelled.size(); i = j) {
        const Group& g = groups[get<0>(spelled[i])];
        EditCandidate m;
        // Each fragment counted once, under the spelling all its reads give, or as discordant.
        map<string, uint32_t> fragments;
        uint32_t total = 0;
        for (size_t f = i; f < spelled.size() && get<0>(spelled[f]) == get<0>(spelled[i]); f = j) {
            for (j = f + 1; j < spelled.size() && get<0>(spelled[j]) == get<0>(spelled[f])
                            && get<1>(spelled[j]) == get<1>(spelled[f]);
                 ++j) {
            }
            if (j > f + 1) {
                ++m.discordant_fragments;
            } else {
                ++fragments[get<2>(spelled[f])];
                ++total;
            }
        }
        string bases;
        for (nid_t id : g.walk) {
            bases += graph.get_sequence(graph.get_handle(id));
        }
        m.node = g.walk[0];
        m.offset = g.lo;
        m.kind = 'M';
        m.span = g.hi - g.lo;
        m.walk = g.walk;
        m.sequence = upper(bases.substr(g.lo, m.span));
        m.ref_fragments = fragments.count(m.sequence) ? fragments[m.sequence] : 0;
        for (const auto& spelling : fragments) {
            // Admitted on the counts a candidate of its own would need, though not checked for a
            // homopolymer, and changing the length only between two walk bases.
            EditCandidate alone;
            alone.kind = spelling.first.size() == m.span ? 'S' : 'I';
            alone.alt_fragments = spelling.second;
            alone.ref_fragments = m.ref_fragments;
            alone.other_fragments = total - spelling.second - m.ref_fragments;
            if (spelling.first != m.sequence && params.bar->admits(alone)
                && (alone.kind == 'S' || (g.lo > 0 && g.hi < bases.size()))) {
                m.alts.push_back(spelling.first);
                m.alt_fragments += spelling.second;
            }
        }
        m.other_fragments = total - m.ref_fragments - m.alt_fragments;
        // Each candidate in some ALT: a substitution's base at its place from either end, an
        // indel's length change.
        bool carried = !m.alts.empty();
        for (size_t c : g.members) {
            const EditCandidate& e = kept[c];
            const size_t left = e.offset - g.lo, right = g.hi - e.offset;
            carried = carried && std::any_of(m.alts.begin(), m.alts.end(), [&](const string& a) {
                return e.kind == 'S' ? (left < a.size() && a[left] == e.alt)
                                           || (right <= a.size() && a[a.size() - right] == e.alt)
                                     : (int64_t)a.size() - (int64_t)m.span
                                           == (e.kind == 'I' ? 1 : -1) * (int64_t)e.sequence.size();
            });
            m.mean_alt_mapq = std::max(m.mean_alt_mapq, e.mean_alt_mapq);
        }
        if (carried) {
            for (size_t c : g.members) {
                kept[c].kind = 0;
            }
            merged.push_back(std::move(m));
        }
    }
    // The merged candidates in place of their members, which are now of kind 0.
    kept.erase(std::remove_if(kept.begin(), kept.end(),
                              [](const EditCandidate& c) { return c.kind == 0; }),
               kept.end());
    kept.insert(kept.end(), std::make_move_iterator(merged.begin()),
                std::make_move_iterator(merged.end()));
    std::sort(kept.begin(), kept.end());
}

}

vector<EditCandidate> count_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params) {
    vector<EditCandidate> kept = count_snv_edits(graph, for_each_read, owns, params);
    vector<LargeRow> large_rows;
    vector<ClipRow> clip_rows;
    vector<EditCandidate> indels =
        count_indel_edits(graph, for_each_read, owns, params,
                          params.large_indels ? &large_rows : nullptr,
                          params.large_indels ? &clip_rows : nullptr);
    if (params.large_indels) {
        vector<EditCandidate> large = count_large_indel_edits(graph, for_each_read, owns, params,
                                                              std::move(large_rows),
                                                              std::move(clip_rows));
        indels.insert(indels.end(), std::make_move_iterator(large.begin()),
                      std::make_move_iterator(large.end()));
    }
    if (!indels.empty()) {
        kept.insert(kept.end(), std::make_move_iterator(indels.begin()),
                    std::make_move_iterator(indels.end()));
        std::sort(kept.begin(), kept.end());
    }
    if (params.bar != nullptr) {
        merge_touching(graph, for_each_read, params, kept);
    }
    return kept;
}

int64_t WalkAlignment::read_offset_at(int64_t x) const {
    auto it = std::upper_bound(blocks.begin(), blocks.end(), x,
                               [](int64_t v, const Block& b) { return v < b.walk; });
    if (it == blocks.begin()) {
        return -1;
    }
    --it;
    if (x >= it->walk + (int64_t)it->length) {
        return -1;
    }
    return (int64_t)it->read + (x - it->walk);
}

WalkAlignment walk_alignment(const vector<const Mapping*>& mappings,
                             const vector<int64_t>& node_start, size_t read_offset,
                             const string& sequence) {
    WalkAlignment out;
    size_t read_pos = read_offset;
    bool started = false;
    // Indels wait for an aligned base after them, so that a trailing clip is left out.
    vector<WalkIndel> pending;
    size_t read_after_aligned = read_offset;
    auto aligned = [&](int64_t walk_pos, size_t length) {
        if (!started) {
            out.begin = walk_pos;
            out.clip_begin = read_pos - read_offset;
            started = true;
        }
        read_after_aligned = read_pos + length;
        for (WalkIndel& p : pending) {
            out.indels.push_back(std::move(p));
        }
        pending.clear();
        if (!out.blocks.empty() && out.blocks.back().walk + (int64_t)out.blocks.back().length
                                       == walk_pos
            && out.blocks.back().read + out.blocks.back().length == read_pos) {
            out.blocks.back().length += length;
        } else {
            out.blocks.push_back(WalkAlignment::Block{walk_pos, read_pos, length});
        }
        out.end = walk_pos + (int64_t)length;
    };
    for (size_t i = 0; i < mappings.size(); ++i) {
        const Mapping& m = *mappings[i];
        int64_t walk_pos = node_start[i] + (int64_t)m.position().offset();
        for (const Edit& e : m.edit()) {
            size_t from = (size_t)e.from_length();
            size_t to = (size_t)e.to_length();
            const size_t both = std::min(from, to);
            if (both > 0) {
                aligned(walk_pos, both);
                walk_pos += (int64_t)both;
                read_pos += both;
                from -= both;
                to -= both;
            }
            if (from > 0) {
                if (started) {
                    if (!pending.empty() && pending.back().deleted > 0
                        && pending.back().at + (int64_t)pending.back().deleted == walk_pos) {
                        pending.back().deleted += (uint32_t)from;
                    } else {
                        WalkIndel d;
                        d.at = walk_pos;
                        d.deleted = (uint32_t)from;
                        pending.push_back(std::move(d));
                    }
                }
                walk_pos += (int64_t)from;
            } else if (to > 0) {
                if (started) {
                    string bases = read_pos + to <= sequence.size() ? sequence.substr(read_pos, to)
                                                                     : string();
                    for (char& b : bases) {
                        b = (char)toupper(b);
                    }
                    if (!pending.empty() && pending.back().deleted == 0
                        && pending.back().at == walk_pos) {
                        pending.back().inserted += bases;
                    } else {
                        WalkIndel ins;
                        ins.at = walk_pos;
                        ins.inserted = std::move(bases);
                        pending.push_back(std::move(ins));
                    }
                }
                read_pos += to;
            }
        }
    }
    if (started) {
        out.clip_end = read_pos - read_after_aligned;
    }
    return out;
}

void left_align(WalkIndel& indel, int64_t floor, const function<char(int64_t)>& base) {
    if (indel.deleted > 0) {
        while (indel.at > floor && base(indel.at - 1) == base(indel.at + indel.deleted - 1)) {
            --indel.at;
        }
    } else if (!indel.inserted.empty()) {
        while (indel.at > floor && base(indel.at - 1) == indel.inserted.back()) {
            indel.inserted.pop_back();
            indel.inserted.insert(indel.inserted.begin(), base(indel.at - 1));
            --indel.at;
        }
    }
}

void left_align(WalkAlignment& aln, const function<char(int64_t)>& base) {
    int64_t floor = aln.begin + 1;
    for (WalkIndel& indel : aln.indels) {
        if (indel.at > floor) {
            left_align(indel, floor, base);
        }
        floor = std::max(floor, indel.at + (int64_t)indel.deleted + 1);
    }
}

int64_t tract_end(const WalkIndel& indel, int64_t limit, const function<char(int64_t)>& base) {
    if (indel.deleted > 0) {
        int64_t p = indel.at;
        while (p + (int64_t)indel.deleted < limit && base(p) == base(p + indel.deleted)) {
            ++p;
        }
        return p + (int64_t)indel.deleted;
    }
    string moving = indel.inserted;
    int64_t p = indel.at;
    while (!moving.empty() && p < limit && base(p) == moving[0]) {
        moving.erase(moving.begin());
        moving.push_back(base(p));
        ++p;
    }
    return p;
}

IndelObservation observe_indel(const WalkAlignment& aln, const WalkIndel& indel, int64_t after,
                               vector<const WalkIndel*>* touching) {
    const int64_t anchor = indel.at - 1;
    if (touching != nullptr) {
        touching->clear();
    }
    if (anchor < aln.begin || after >= aln.end) {
        return IndelObservation::NONE;
    }
    // The indels are in walk order and do not overlap, so their ends are too.
    auto it = std::lower_bound(aln.indels.begin(), aln.indels.end(), anchor,
                               [](const WalkIndel& d, int64_t a) {
                                   return d.at + (int64_t)std::max<uint32_t>(d.deleted, 1) <= a;
                               });
    size_t count = 0;
    const WalkIndel* only = nullptr;
    for (; it != aln.indels.end() && it->at <= after; ++it) {
        const bool touches = it->deleted > 0 ? it->at + (int64_t)it->deleted > anchor
                                             : it->at > anchor;
        if (touches) {
            ++count;
            only = &*it;
            if (touching != nullptr) {
                touching->push_back(&*it);
            }
        }
    }
    if (count == 0) {
        return IndelObservation::REF;
    }
    return count == 1 && *only == indel ? IndelObservation::ALT : IndelObservation::OTHER;
}

IndelObservation observe_large_indel(const WalkAlignment& aln, const WalkIndel& indel,
                                     int64_t after, const LargeIndelParams& params,
                                     int64_t reach_begin, int64_t reach_end, int* by_clip) {
    if (by_clip != nullptr) {
        *by_clip = 0;
    }
    const bool deletion = indel.deleted > 0;
    const int64_t length = (int64_t)indel.length();
    // The change's two ends, the bases after the anchor and after the deleted bases, and how far
    // its tract lets it slide right.
    const int64_t left = indel.at;
    const int64_t right = deletion ? indel.at + length : indel.at;
    const int64_t room = std::max<int64_t>(0, after - right);
    const int64_t lo = left - 1 - params.slop;
    const int64_t hi = after + params.slop;
    // The read's bases of the change's kind near it.
    int64_t near = 0;
    for (const WalkIndel& d : aln.indels) {
        if (d.at > hi) {
            break;
        }
        if (deletion && d.deleted > 0 && d.at + (int64_t)d.deleted >= lo) {
            near += d.deleted;
        } else if (!deletion && d.deleted == 0 && d.at >= lo) {
            near += (int64_t)d.inserted.size();
        }
    }
    const double tolerance = params.length_ratio * (double)length;
    if (std::abs((double)(near - length)) <= tolerance) {
        return IndelObservation::ALT;
    }
    const bool clip_before = aln.clip_end > 0 && aln.end >= lo && aln.end <= hi;
    const bool clip_after = aln.clip_begin > 0 && aln.begin >= lo && aln.begin <= hi;
    if (deletion) {
        // Clipped where the read's aligned bases meet one end, the rest of the read being the
        // bases past the other.
        if (aln.clip_end >= params.min_clip && aln.end >= left - params.slop
            && aln.end <= left + room + params.slop) {
            if (by_clip != nullptr) {
                *by_clip = 1;
            }
            return IndelObservation::ALT;
        }
        if (aln.clip_begin >= params.min_clip && aln.begin >= right - params.slop
            && aln.begin <= right + room + params.slop) {
            if (by_clip != nullptr) {
                *by_clip = 2;
            }
            return IndelObservation::ALT;
        }
    }
    auto spans = [&](int64_t from, int64_t to) {
        return reach_begin <= from - params.flank && reach_end >= to + params.flank;
    };
    if (!spans(left, left + room) && !(deletion && spans(right, right + room))) {
        return IndelObservation::NONE;
    }
    if ((double)near < tolerance && !clip_before && !clip_after) {
        return IndelObservation::REF;
    }
    return IndelObservation::OTHER;
}

void EditOutcomes::add(EditOutcome&& outcome) {
    by_thread.at(omp_get_thread_num()).push_back(std::move(outcome));
}

vector<EditOutcome> EditOutcomes::sorted() const {
    vector<EditOutcome> all;
    for (const auto& list : by_thread) {
        all.insert(all.end(), list.begin(), list.end());
    }
    std::sort(all.begin(), all.end(), [](const EditOutcome& a, const EditOutcome& b) {
        return tie(a.node, a.offset, a.kind, a.alt, a.sequence)
               < tie(b.node, b.offset, b.kind, b.alt, b.sequence);
    });
    return all;
}

void for_each_candidate_on(const vector<EditCandidate>& sorted, const vector<nid_t>& nodes,
                           const function<void(const EditCandidate&)>& iteratee) {
    auto it = sorted.begin();
    for (nid_t node : nodes) {
        it = std::lower_bound(it, sorted.end(), node,
                              [](const EditCandidate& c, nid_t n) { return c.node < n; });
        for (; it != sorted.end() && it->node == node; ++it) {
            iteratee(*it);
        }
    }
}

}
