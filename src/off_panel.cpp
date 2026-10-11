#include "off_panel.hpp"

#include <algorithm>
#include <cctype>
#include <tuple>
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

/// Walk the read's aligned bases on the nodes `wanted` accepts, each of quality at least `min_q`,
/// giving each one's node, forward offset and forward base (0 for a match).
void for_each_aligned_base(const HandleGraph& graph, const Alignment& aln, int min_q,
                           const function<bool(nid_t)>& wanted,
                           const function<void(nid_t, uint32_t, char)>& visit) {
    const string& qual = aln.quality();
    size_t read_pos = 0;
    for (const Mapping& m : aln.path().mapping()) {
        const nid_t node = m.position().node_id();
        const bool take = wanted(node);
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

    void load(const Alignment& aln) {
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

/// The indels among `kept`'s candidates: counted as `count_snv_edits` counts substitutions, with
/// the observation classes of `observe_indel`.
vector<EditCandidate> count_indel_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params) {

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
        walk.load(aln);
        const uint64_t fragment = name_hash(aln.name());
        for (const WalkIndel& indel : walk.walked.indels) {
            const int64_t anchor = indel.at - 1;
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
    struct Seen {
        size_t key;
        uint64_t fragment;
        uint8_t observed;
        bool operator<(const Seen& o) const {
            return tie(key, fragment, observed) < tie(o.key, o.fragment, o.observed);
        }
    };
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
        walk.load(aln);
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
        i = end;
        if (c.walk.empty() || c.alt_fragments < params.min_fragments
            || c.fraction() < params.min_fraction) {
            continue;
        }
        // Another walk of the graph spelling the same change is a panel allele.
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
        if (spelled_by_other_walk(graph, handles, alt)) {
            continue;
        }
        c.mean_alt_mapq =
            mapq_count[key] == 0 ? 0.0f : (float)(mapq_total[key] / (double)mapq_count[key]);
        kept.push_back(std::move(c));
    }
    return kept;
}

}

vector<EditCandidate> count_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params) {
    vector<EditCandidate> kept = count_snv_edits(graph, for_each_read, owns, params);
    vector<EditCandidate> indels = count_indel_edits(graph, for_each_read, owns, params);
    if (!indels.empty()) {
        kept.insert(kept.end(), std::make_move_iterator(indels.begin()),
                    std::make_move_iterator(indels.end()));
        std::sort(kept.begin(), kept.end());
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
    auto aligned = [&](int64_t walk_pos, size_t length) {
        if (!started) {
            out.begin = walk_pos;
            started = true;
        }
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
