#include "off_panel.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <map>

#include <omp.h>

#include "../src/alignment.hpp"
#include "../src/utility.hpp"

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

    /// Order observations by node, offset, base and fragment, so that each fragment's reads of one
    /// ALT at one base are adjacent.
    std::strong_ordering operator<=>(const Observation& other) const {
        return tie(node, offset, base, fragment)
               <=> tie(other.node, other.offset, other.base, other.fragment);
    }

    /// Whether two observations are of one base by one fragment, whatever their mapping qualities.
    bool operator==(const Observation& other) const {
        return (*this <=> other) == 0;
    }
};

/// Whether every node the read's path visits is in the graph. A read aligned to a larger graph
/// can leave this one; such a read cannot be walked, so it adds no off-panel evidence.
bool path_in_graph(const HandleGraph& graph, const Alignment& aln) {
    return std::ranges::all_of(aln.path().mapping(), [&](const Mapping& m) {
        return graph.has_node(m.position().node_id());
    });
}

/// Whether read base `r` has quality at least `min_q`. A read with no qualities passes.
bool base_quality_ok(const string& qual, size_t r, int min_q) {
    return qual.empty() || (r < qual.size() && (int)(uint8_t)qual[r] >= min_q);
}

/// Read base `i` of the aligned edit `e` along the node's forward strand into `base`: 0 for a
/// match, or the substituted base. False for a substituted base other than A, C, G or T.
bool forward_read_base(const Edit& e, size_t i, bool reverse, char& base) {
    base = 0;
    if (e.sequence().empty()) {
        return true;
    }
    base = (char)toupper(e.sequence()[i]);
    if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
        return false;
    }
    if (reverse) {
        base = reverse_complement(base);
    }
    return true;
}

/// Give `visit` the node, forward offset and forward base (0 for a match) of each aligned base of
/// mapping `m` whose quality is at least `min_q`. The node has `node_len` bases, and the mapping's
/// first base is at read offset `read_pos`. Return the read offset after the mapping.
///
/// With `node_len` 0, visit nothing.
template<typename Visit>
size_t visit_mapping_bases(const Mapping& m, size_t node_len, size_t read_pos, const string& qual,
                           int min_q, Visit& visit) {
    const nid_t node = m.position().node_id();
    const bool reverse = m.position().is_reverse();
    auto node_pos = (size_t)m.position().offset();
    for (const Edit& e : m.edit()) {
        const auto from = (size_t)e.from_length();
        const auto to = (size_t)e.to_length();
        const bool aligned = node_len > 0 && from == to && from > 0;
        for (size_t i = 0; aligned && i < from; ++i) {
            const size_t visited = node_pos + i;
            char base = 0;
            if (base_quality_ok(qual, read_pos + i, min_q) && visited < node_len
                && forward_read_base(e, i, reverse, base)) {
                visit(node, (uint32_t)(reverse ? node_len - 1 - visited : visited), base);
            }
        }
        node_pos += from;
        read_pos += to;
    }
    return read_pos;
}

/// Walk the read's aligned bases on the nodes `wanted` accepts, each of quality at least `min_q`,
/// giving each one's node, forward offset and forward base (0 for a match) to `visit`.
template<typename Wanted, typename Visit>
void for_each_aligned_base(const HandleGraph& graph, const Alignment& aln, int min_q,
                           const Wanted& wanted, Visit&& visit) {
    size_t read_pos = 0;
    for (const Mapping& m : aln.path().mapping()) {
        const nid_t node = m.position().node_id();
        const bool take = wanted(node) && graph.has_node(node);
        const size_t node_len = take ? graph.get_length(graph.get_handle(node)) : 0;
        read_pos = visit_mapping_bases(m, node_len, read_pos, aln.quality(), min_q, visit);
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

/// The source of the reads to count, which calls its argument on each read.
using ReadSource = function<void(const function<void(const Alignment&)>&)>;

/// Collect the substitutions the reads make on the nodes `owns` accepts, sorted.
vector<Observation> collect_substitutions(const HandleGraph& graph, const ReadSource& for_each_read,
                                          const function<bool(nid_t)>& owns,
                                          const EditCountParams& params) {
    std::hash<string> name_hash;
    vector<Observation> alts;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(graph, aln, params.min_base_quality, owns,
                              [&](nid_t node, uint32_t offset, char base) {
                                  if (base != 0) {
                                      alts.emplace_back(node, offset, base, fragment,
                                                        aln.mapping_quality());
                                  }
                              });
    });
    std::ranges::sort(alts);
    return alts;
}

/// Return the (node, offset) of each substitution in the sorted `alts` that at least
/// `min_fragments` fragments carry, sorted and without duplicates, and put their nodes in `nodes`.
vector<pair<nid_t, uint32_t>> recurring_positions(const vector<Observation>& alts,
                                                  size_t min_fragments,
                                                  unordered_set<nid_t>& nodes) {
    vector<pair<nid_t, uint32_t>> positions;
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
        if (fragments >= min_fragments) {
            positions.emplace_back(alts[i].node, alts[i].offset);
            nodes.insert(alts[i].node);
        }
        i = j;
    }
    std::ranges::sort(positions);
    const auto duplicates = std::ranges::unique(positions);
    positions.erase(duplicates.begin(), duplicates.end());
    return positions;
}

/// Collect every base the reads put at the sorted `positions`, whose nodes are `nodes`, sorted by
/// position and then by fragment.
vector<Observation> collect_covering(const HandleGraph& graph, const ReadSource& for_each_read,
                                     const function<bool(nid_t)>& owns,
                                     const EditCountParams& params,
                                     const vector<pair<nid_t, uint32_t>>& positions,
                                     const unordered_set<nid_t>& nodes) {
    std::hash<string> name_hash;
    vector<Observation> covering;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(
            graph, aln, params.min_base_quality,
            [&](nid_t node) { return owns(node) && nodes.contains(node); },
            [&](nid_t node, uint32_t offset, char base) {
                if (std::ranges::binary_search(positions, make_pair(node, offset))) {
                    covering.emplace_back(node, offset, base, fragment, 0);
                }
            });
    });
    std::ranges::sort(covering, [](const Observation& a, const Observation& b) {
        return tie(a.node, a.offset, a.fragment, a.base)
               < tie(b.node, b.offset, b.fragment, b.base);
    });
    return covering;
}

/// Return each fragment's one base among `covering[i, end)`, which is one position sorted by
/// fragment, or 1 for a fragment whose reads disagree.
vector<char> fragment_bases(const vector<Observation>& covering, size_t i, size_t end) {
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
    return fragment_base;
}

/// Return the mean mapping quality of the reads in the sorted `alts` that read `alt` at `offset`
/// on `node`.
float mean_substitution_mapq(const vector<Observation>& alts, nid_t node, uint32_t offset,
                             char alt) {
    const Observation probe{node, offset, alt, 0, 0};
    double total = 0.0;
    size_t n = 0;
    for (auto it = std::ranges::lower_bound(alts, probe);
         it != alts.end() && it->node == node && it->offset == offset && it->base == alt; ++it) {
        total += it->mapq;
        ++n;
    }
    return n == 0 ? 0.0f : (float)(total / (double)n);
}

/// Return the substitution of `ref` by `alt` at `offset` on `node`, with its fragments counted
/// from each fragment's base there (`fragment_bases`).
EditCandidate tally_substitution(nid_t node, uint32_t offset, char ref, char alt,
                                 const vector<char>& fragment_base) {
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
    return c;
}

/// Count the substitutions, in two passes over the reads. The first finds the positions that
/// enough fragments alter, and the second counts every fragment's base at each of them.
vector<EditCandidate> count_snv_edits(const HandleGraph& graph, const ReadSource& for_each_read,
                                      const function<bool(nid_t)>& owns,
                                      const EditCountParams& params) {
    vector<EditCandidate> kept;
    const vector<Observation> alts = collect_substitutions(graph, for_each_read, owns, params);
    if (alts.empty()) {
        return kept;
    }
    unordered_set<nid_t> position_nodes;
    const vector<pair<nid_t, uint32_t>> positions =
        recurring_positions(alts, params.min_fragments, position_nodes);
    if (positions.empty()) {
        return kept;
    }
    const vector<Observation> covering =
        collect_covering(graph, for_each_read, owns, params, positions, position_nodes);

    for (size_t i = 0; i < covering.size();) {
        size_t end = i;
        while (end < covering.size() && covering[end].node == covering[i].node
               && covering[end].offset == covering[i].offset) {
            ++end;
        }
        const vector<char> fragment_base = fragment_bases(covering, i, end);
        const nid_t node = covering[i].node;
        const uint32_t offset = covering[i].offset;
        const auto ref = (char)toupper(graph.get_base(graph.get_handle(node), offset));
        for (char alt : {'A', 'C', 'G', 'T'}) {
            if (alt == ref) {
                continue;
            }
            EditCandidate c = tally_substitution(node, offset, ref, alt, fragment_base);
            if (c.alt_fragments < params.min_fragments || c.fraction() < params.min_fraction
                || spelled_by_sibling(graph, node, offset, alt)) {
                continue;
            }
            c.mean_alt_mapq = mean_substitution_mapq(alts, node, offset, alt);
            kept.push_back(c);
        }
        i = end;
    }
    std::ranges::sort(kept);
    return kept;
}

/// A cursor over one read's walk through the graph at a time: its mappings, read so that the first
/// reads its node forward, the walk base each mapping's node starts at, and its alignment along the
/// walk, indels left-aligned.
///
/// Walk bases count from the first base of the first mapping's node. Each `load` replaces the walk.
class ReadWalk {
public:
    explicit ReadWalk(const HandleGraph& graph) : graph(graph) {}
    // `base_fn` captures `this`, so a copy would read the original's walk.
    ReadWalk(const ReadWalk&) = delete;
    ReadWalk& operator=(const ReadWalk&) = delete;

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
        loaded = use;
        mappings.clear();
        walk_handles.clear();
        starts.clear();
        int64_t at = 0;
        for (const Mapping& m : use->path().mapping()) {
            const handle_t h = graph.get_handle(m.position().node_id(), m.position().is_reverse());
            mappings.push_back(&m);
            walk_handles.push_back(h);
            starts.push_back(at);
            at += (int64_t)graph.get_length(h);
        }
        starts.push_back(at);
        walked = walk_alignment(mappings, starts, 0, use->sequence());
        left_align(walked, base_fn);
        return true;
    }

    /// The read loaded, as its walk reads it.
    const Alignment& read() const {
        return *loaded;
    }

    /// The handle of each mapping's node, oriented as the walk reads it.
    const vector<handle_t>& handles() const {
        return walk_handles;
    }

    /// The walk base each mapping's node starts at, then the walk's length.
    const vector<int64_t>& node_start() const {
        return starts;
    }

    /// The read's alignment along the walk, indels left-aligned.
    const WalkAlignment& alignment() const {
        return walked;
    }

    /// The mapping holding walk base `x`.
    size_t mapping_at(int64_t x) const {
        return (size_t)(std::ranges::upper_bound(starts, x) - starts.begin()) - 1;
    }

    /// The walk's base at `x`, upper-cased.
    char base(int64_t x) const {
        const size_t i = mapping_at(x);
        return (char)toupper(graph.get_base(walk_handles[i], (size_t)(x - starts[i])));
    }

    /// Return `base` as a `std::function`, the form `left_align` and `tract_end` take.
    const function<char(int64_t)>& base_function() const {
        return base_fn;
    }

    /// A hash of the nodes, with their orientations, of mappings `first` to `last`, or 0 if any
    /// reads its node backward.
    uint64_t walk_hash(size_t first, size_t last) const {
        uint64_t h = 1469598103934665603ULL;
        for (size_t i = first; i <= last; ++i) {
            if (graph.get_is_reverse(walk_handles[i])) {
                return 0;
            }
            h = (h ^ (uint64_t)graph.get_id(walk_handles[i])) * 1099511628211ULL;
        }
        return h == 0 ? 1 : h;
    }

    /// Whether the read base on walk base `x`, if one is, has quality at least `min_q`.
    bool quality_ok(int64_t x, int min_q) const {
        const int64_t r = walked.read_offset_at(x);
        return r < 0 || base_quality_ok(loaded->quality(), (size_t)r, min_q);
    }

private:
    const HandleGraph& graph;
    /// The read reverse-complemented, when its first mapping reads its node backward.
    Alignment flipped;
    /// The read loaded: the caller's, or `flipped`.
    const Alignment* loaded = nullptr;
    vector<const Mapping*> mappings;
    vector<handle_t> walk_handles;
    vector<int64_t> starts;
    WalkAlignment walked;
    /// Calls `base` on this walk.
    const function<char(int64_t)> base_fn = [this](int64_t x) { return base(x); };
};

/// Whether the read's alignment has an insertion or a deletion.
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
    /// The anchor's node, read forward.
    nid_t node;
    /// The anchor's offset on `node`.
    uint32_t offset;
    /// 'D' for a deletion, 'I' for an insertion.
    char kind;
    uint32_t length;
    /// A hash of the inserted bases, or 0 for a deletion.
    uint64_t sequence_hash;
    uint64_t fragment;
    /// The bases from the anchor to the first base after the tract.
    uint32_t span;
    /// The `ReadWalk::walk_hash` of the nodes from the anchor's to that of the first base after the
    /// tract.
    uint64_t walk_hash;

    /// Whether two rows are of the same indel.
    bool same_key(const IndelObservationRow& o) const {
        return node == o.node && offset == o.offset && kind == o.kind && length == o.length
               && sequence_hash == o.sequence_hash;
    }
    /// Order rows by key, then by fragment, then by walk.
    std::strong_ordering operator<=>(const IndelObservationRow& o) const {
        return tie(node, offset, kind, length, sequence_hash, fragment, span, walk_hash)
               <=> tie(o.node, o.offset, o.kind, o.length, o.sequence_hash, o.fragment, o.span,
                       o.walk_hash);
    }
    bool operator==(const IndelObservationRow& o) const = default;
};

/// What one read says about one key, by its fragment.
struct Seen {
    size_t key;
    uint64_t fragment;
    /// The `IndelObservation`.
    uint8_t observed;
    auto operator<=>(const Seen& o) const = default;
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

/// Whether a walk of the graph from the first handle of `walk` to its last spells `alt`.
///
/// The search gives up, answering false, after taking 4096 handles off its stack.
bool spelled_by_other_walk(const HandleGraph& graph, const vector<handle_t>& walk,
                           string_view alt) {
    vector<pair<handle_t, size_t>> stack{{walk.front(), 0}};
    size_t budget = 4096;
    while (!stack.empty() && budget > 0) {
        --budget;
        const auto [handle, spelled] = stack.back();
        stack.pop_back();
        const string seq = graph.get_sequence(handle);
        if (spelled + seq.size() > alt.size() || alt.compare(spelled, seq.size(), seq) != 0) {
            continue;
        }
        const size_t next = spelled + seq.size();
        if (next == alt.size()) {
            if (handle == walk.back()) {
                return true;
            }
            continue;
        }
        graph.follow_edges(handle, false, [&](const handle_t& n) {
            stack.emplace_back(n, next);
        });
    }
    return false;
}

/// Whether another walk of the graph spells the change of the indel candidate `c`, as a panel
/// allele.
///
/// The change alters the length of the candidate's walk, so a walk that spells it is another walk.
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

/// One read's indel of `LARGE_INDEL_MIN_LENGTH` bases or more, anchored on a node the count owns.
struct LargeRow {
    char kind;
    nid_t node;
    uint32_t offset;
    uint32_t length;
    uint64_t fragment;
    string inserted;
    /// The nodes, read forward, from the anchor's to that of the first base after the tract, or
    /// empty where the read does not reach that far.
    vector<nid_t> walk;
    /// The bases from the anchor to the first base after the tract, or 0 where `walk` is empty.
    uint32_t span;

    /// Order rows by kind, anchor and length, so that a cluster of one indel's rows is a run.
    std::strong_ordering operator<=>(const LargeRow& o) const {
        return tie(kind, node, offset, length, fragment, inserted, span, walk)
               <=> tie(o.kind, o.node, o.offset, o.length, o.fragment, o.inserted, o.span, o.walk);
    }
    bool operator==(const LargeRow& o) const = default;
};

/// One read's clip of `min_clip` bases or more.
struct ClipRow {
    /// 'R' for a clip after the read's last aligned base, 'L' for one before its first.
    char side;
    /// The node of the aligned base beside the clip, read forward.
    nid_t node;
    /// The offset of the aligned base beside the clip on `node`.
    uint32_t offset;
    /// The `min_clip` clipped bases nearest the aligned bases, upper-cased.
    string seed;
    uint64_t fragment;

    auto operator<=>(const ClipRow& o) const = default;
};

/// Return `bases` upper-cased.
string upper(string bases) {
    for (char& b : bases) {
        b = (char)toupper(b);
    }
    return bases;
}

/// File the large indel `indel` of the loaded `walk` in `rows`, if its anchor lies on a node that
/// `owns` accepts and the walk reads forward.
template<typename Owns>
void add_large_row(const HandleGraph& graph, const ReadWalk& walk, const WalkIndel& indel,
                   uint64_t fragment, const Owns& owns, vector<LargeRow>& rows) {
    const int64_t anchor = indel.at - 1;
    const size_t i = walk.mapping_at(anchor);
    const handle_t h = walk.handles()[i];
    if (graph.get_is_reverse(h) || !owns(graph.get_id(h))) {
        return;
    }
    LargeRow row{indel.deleted > 0 ? 'D' : 'I', graph.get_id(h),
                 (uint32_t)(anchor - walk.node_start()[i]), indel.length(), fragment,
                 indel.inserted, {}, 0};
    if (const int64_t after = tract_end(indel, walk.alignment().end, walk.base_function());
        after < walk.alignment().end && walk.walk_hash(i, walk.mapping_at(after)) != 0) {
        for (size_t m = i; m <= walk.mapping_at(after); ++m) {
            row.walk.push_back(graph.get_id(walk.handles()[m]));
        }
        row.span = (uint32_t)(after - anchor);
    }
    rows.push_back(std::move(row));
}

/// File the clips of the loaded `walk` of `min_clip` bases or more in `rows`, where the aligned base
/// beside the clip is on a node the walk reads forward.
void add_clip_rows(const HandleGraph& graph, const ReadWalk& walk, uint64_t fragment,
                   const LargeIndelParams& params, vector<ClipRow>& rows) {
    const WalkAlignment& w = walk.alignment();
    if (w.blocks.empty()) {
        return;
    }
    const string& sequence = walk.read().sequence();
    if (w.clip_end >= params.min_clip) {
        const size_t i = walk.mapping_at(w.end - 1);
        const size_t from = w.blocks.back().read + w.blocks.back().length;
        if (!graph.get_is_reverse(walk.handles()[i]) && from + params.min_clip <= sequence.size()) {
            rows.emplace_back('R', graph.get_id(walk.handles()[i]),
                              (uint32_t)(w.end - 1 - walk.node_start()[i]),
                              upper(sequence.substr(from, params.min_clip)), fragment);
        }
    }
    if (w.clip_begin >= params.min_clip) {
        const size_t i = walk.mapping_at(w.begin);
        const size_t to = w.blocks.front().read;
        if (!graph.get_is_reverse(walk.handles()[i]) && to >= params.min_clip
            && to <= sequence.size()) {
            rows.emplace_back('L', graph.get_id(walk.handles()[i]),
                              (uint32_t)(w.begin - walk.node_start()[i]),
                              upper(sequence.substr(to - params.min_clip, params.min_clip)),
                              fragment);
        }
    }
}

/// File the loaded read's indels whose anchors lie on nodes `owns` accepts, read forward: those
/// up to `max_indel_length` bases in `rows`, and, unless `large_rows` is null, those of
/// `LARGE_INDEL_MIN_LENGTH` to `large.max_length` bases in `large_rows`.
template<typename Owns>
void file_read_indels(const HandleGraph& graph, const ReadWalk& walk, uint64_t fragment,
                      const Owns& owns, const EditCountParams& params,
                      vector<IndelObservationRow>& rows, vector<LargeRow>* large_rows) {
    std::hash<string> sequence_hasher;
    const WalkAlignment& walked = walk.alignment();
    for (const WalkIndel& indel : walked.indels) {
        const int64_t anchor = indel.at - 1;
        if (large_rows != nullptr && anchor >= walked.begin
            && indel.length() >= LARGE_INDEL_MIN_LENGTH
            && indel.length() <= params.large.max_length) {
            add_large_row(graph, walk, indel, fragment, owns, *large_rows);
        }
        if (anchor < walked.begin || indel.length() > params.max_indel_length) {
            continue;
        }
        const size_t i = walk.mapping_at(anchor);
        const handle_t h = walk.handles()[i];
        if (graph.get_is_reverse(h) || !owns(graph.get_id(h))) {
            continue;
        }
        IndelObservationRow row{graph.get_id(h), (uint32_t)(anchor - walk.node_start()[i]),
                                indel.deleted > 0 ? 'D' : 'I', indel.length(),
                                indel.deleted > 0 ? 0 : sequence_hasher(indel.inserted),
                                fragment, 0, 0};
        // Record the walk to the end of the tract only where the read covers that far.
        if (const int64_t after = tract_end(indel, walked.end, walk.base_function());
            after < walked.end) {
            if (const uint64_t h_walk = walk.walk_hash(i, walk.mapping_at(after)); h_walk != 0) {
                row.span = (uint32_t)(after - anchor);
                row.walk_hash = h_walk;
            }
        }
        rows.push_back(row);
    }
}

/// Return one row for each key in the sorted `rows` that at least `min_fragments` fragments carry,
/// with the walk most of the key's reads take to the end of its tract. Keys no read covers to the
/// end of the tract are left out. The keys come out sorted by node.
vector<IndelObservationRow> recurring_indel_keys(const vector<IndelObservationRow>& rows,
                                                 size_t min_fragments) {
    vector<IndelObservationRow> keys;
    for (size_t i = 0; i < rows.size();) {
        size_t j = i;
        size_t fragments = 0;
        // Reads by the walk, as span and hash, that they take to the end of the tract.
        map<pair<uint32_t, uint64_t>, size_t> walks;
        for (; j < rows.size() && rows[j].same_key(rows[i]); ++j) {
            if (j == i || rows[j].fragment != rows[j - 1].fragment) {
                ++fragments;
            }
            if (rows[j].span > 0) {
                ++walks[make_pair(rows[j].span, rows[j].walk_hash)];
            }
        }
        if (fragments >= min_fragments && !walks.empty()) {
            // Ties go to the first walk in map order.
            const auto best = std::ranges::max_element(
                walks, {}, [](const auto& entry) { return entry.second; });
            IndelObservationRow key = rows[i];
            key.span = best->first.first;
            key.walk_hash = best->first.second;
            keys.push_back(key);
        }
        i = j;
    }
    return keys;
}

/// A tally of what the reads say about each recurring short indel, by fragment.
///
/// The candidate for each indel is filled in from the first read that carries it, since the key
/// fixes everything but the counts.
class IndelTally {
public:
    /// Tally reads against `keys`, sorted by node, each with the walk its reads take.
    explicit IndelTally(vector<IndelObservationRow> keys)
        : keys(std::move(keys)), by_key(this->keys.size()) {
        for (const IndelObservationRow& key : this->keys) {
            anchor_nodes.insert(key.node);
        }
    }

    /// Whether the read visits the anchor node of some indel.
    bool near(const Alignment& aln) const {
        return std::ranges::any_of(aln.path().mapping(), [&](const Mapping& m) {
            return anchor_nodes.contains(m.position().node_id());
        });
    }

    /// Record what the loaded `walk`, from `fragment` at mapping quality `mapq`, says about each
    /// indel anchored on a node it reads forward. Only reads whose bases at the anchor and at the
    /// end of the tract have quality at least `min_q` are counted.
    void observe(const HandleGraph& graph, const ReadWalk& walk, uint64_t fragment, int32_t mapq,
                 int min_q) {
        for (size_t i = 0; i + 1 < walk.node_start().size(); ++i) {
            const handle_t h = walk.handles()[i];
            const nid_t id = graph.get_id(h);
            if (graph.get_is_reverse(h) || !anchor_nodes.contains(id)) {
                continue;
            }
            for (auto k = std::ranges::lower_bound(keys, id, {}, &IndelObservationRow::node);
                 k != keys.end() && k->node == id; ++k) {
                observe_key(graph, walk, i, (size_t)(k - keys.begin()), fragment, mapq, min_q);
            }
        }
    }

    /// Return the candidates that pass the floors of `params` and that no other walk of the graph
    /// spells.
    vector<EditCandidate> candidates(const HandleGraph& graph, const EditCountParams& params) {
        vector<EditCandidate> kept;
        std::ranges::sort(seen);
        for (size_t i = 0; i < seen.size();) {
            const size_t key = seen[i].key;
            size_t end = i;
            while (end < seen.size() && seen[end].key == key) {
                ++end;
            }
            KeyTally& tally = by_key[key];
            EditCandidate& c = tally.candidate;
            count_fragments(seen, i, end, c);
            i = end;
            if (c.walk.empty() || c.alt_fragments < params.min_fragments
                || c.fraction() < params.min_fraction || spelled_by_panel(graph, c)) {
                continue;
            }
            c.mean_alt_mapq = tally.mapq_count == 0
                                  ? 0.0f
                                  : (float)(tally.mapq_total / (double)tally.mapq_count);
            kept.push_back(std::move(c));
        }
        return kept;
    }

private:
    /// What the reads carrying one indel have given so far.
    struct KeyTally {
        /// The indel's candidate, its `walk` empty until a read carries it.
        EditCandidate candidate;
        double mapq_total = 0.0;
        size_t mapq_count = 0;
    };

    /// Record what `walk` says about key number `key`, anchored on the node of mapping `i`.
    void observe_key(const HandleGraph& graph, const ReadWalk& walk, size_t i, size_t key,
                     uint64_t fragment, int32_t mapq, int min_q) {
        const IndelObservationRow& k = keys[key];
        const int64_t anchor = walk.node_start()[i] + k.offset;
        const int64_t after = anchor + k.span;
        if (anchor < walk.alignment().begin || after >= walk.alignment().end
            || walk.walk_hash(i, walk.mapping_at(after)) != k.walk_hash
            || !walk.quality_ok(anchor, min_q) || !walk.quality_ok(after, min_q)) {
            return;
        }
        WalkIndel wanted;
        wanted.at = anchor + 1;
        if (k.kind == 'D') {
            wanted.deleted = k.length;
        }
        IndelObservation observed = observe_indel(walk.alignment(), wanted, after, &touching);
        // `wanted` has no inserted bases, so an insertion of the key's length and bases at its
        // place reads as another indel. It is this one.
        if (observed == IndelObservation::OTHER && touching.size() == 1 && k.kind == 'I'
            && touching[0]->at == wanted.at && touching[0]->deleted == 0
            && touching[0]->inserted.size() == k.length
            && sequence_hasher(touching[0]->inserted) == k.sequence_hash) {
            observed = IndelObservation::ALT;
        }
        seen.emplace_back(key, fragment, (uint8_t)observed);
        if (observed != IndelObservation::ALT) {
            return;
        }
        KeyTally& tally = by_key[key];
        tally.mapq_total += mapq;
        ++tally.mapq_count;
        if (tally.candidate.walk.empty()) {
            fill_candidate(tally.candidate, graph, walk, i, k, anchor);
        }
    }

    /// Fill in `c` for the indel `k` from `walk`, which carries it anchored at walk base `anchor`
    /// on the node of mapping `i`.
    void fill_candidate(EditCandidate& c, const HandleGraph& graph, const ReadWalk& walk, size_t i,
                        const IndelObservationRow& k, int64_t anchor) const {
        const int64_t after = anchor + k.span;
        c.node = k.node;
        c.offset = k.offset;
        c.kind = k.kind;
        c.alt = 0;
        c.ref = walk.base(anchor);
        c.span = k.span;
        if (k.kind == 'D') {
            for (int64_t x = anchor + 1; x <= anchor + (int64_t)k.length; ++x) {
                c.sequence.push_back(walk.base(x));
            }
        } else {
            c.sequence = touching[0]->inserted;
        }
        for (size_t m = i; m <= walk.mapping_at(after); ++m) {
            c.walk.push_back(graph.get_id(walk.handles()[m]));
        }
        const char first = c.sequence[0];
        for (int64_t x = anchor + 1; x < walk.alignment().end && walk.base(x) == first; ++x) {
            ++c.homopolymer;
        }
    }

    /// The indels, sorted by node.
    vector<IndelObservationRow> keys;
    /// The tally for each of `keys`.
    vector<KeyTally> by_key;
    /// The nodes of the indels' anchors.
    unordered_set<nid_t> anchor_nodes;
    /// What each read said about each indel.
    vector<Seen> seen;
    /// The read's indels touching the indel being observed, reused between reads.
    vector<const WalkIndel*> touching;
    std::hash<string> sequence_hasher;
};

/// Count the short indels: those of at most `params.max_indel_length` bases, in two passes over
/// the reads. The first finds the indels that at least `params.min_fragments` fragments carry, and
/// the second classes what each read says about each of them with `observe_indel`.
///
/// Unless they are null, the first pass also files the reads' large indels in `large_rows` and
/// their long clips in `clip_rows`, for `count_large_indel_edits`.
vector<EditCandidate> count_indel_edits(const HandleGraph& graph, const ReadSource& for_each_read,
                                        const function<bool(nid_t)>& owns,
                                        const EditCountParams& params,
                                        vector<LargeRow>* large_rows, vector<ClipRow>* clip_rows) {
    std::hash<string> name_hash;
    ReadWalk walk(graph);

    // Every read's indels, wherever the aligner put them.
    vector<IndelObservationRow> rows;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq || !has_indel(aln) || !walk.load(aln)) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        file_read_indels(graph, walk, fragment, owns, params, rows, large_rows);
        if (clip_rows != nullptr) {
            add_clip_rows(graph, walk, fragment, params.large, *clip_rows);
        }
    });
    if (rows.empty()) {
        return {};
    }
    std::ranges::sort(rows);
    vector<IndelObservationRow> keys = recurring_indel_keys(rows, params.min_fragments);
    if (keys.empty()) {
        return {};
    }

    // What each read covering a key's anchor through the end of its tract, on the key's walk, says
    // about it.
    IndelTally tally(std::move(keys));
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq || !tally.near(aln) || !walk.load(aln)) {
            return;
        }
        tally.observe(graph, walk, name_hash(aln.name()), aln.mapping_quality(),
                      params.min_base_quality);
    });
    return tally.candidates(graph, params);
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
    /// Where each node of `walk` starts along it.
    vector<int64_t> starts;
    /// The bases of `walk`, upper-cased.
    string spelling;
    /// The bases from the anchor to the first base after the tract.
    uint32_t span = 0;
    /// The number of fragments that found it.
    size_t support = 0;
    /// Whether only clips found it.
    bool clips_only = false;
};

/// Find in `next` the node, read forward, that an unbranched walk going forward (or backward) from
/// `at` goes on to: the end of the only edge out of `at` that way, when that node has no other edge
/// in and is not the node `start`. False if there is none.
bool sole_successor(const HandleGraph& graph, handle_t at, bool forward, nid_t start,
                    handle_t& next) {
    size_t n = 0;
    next = at;
    graph.follow_edges(at, !forward, [&](const handle_t& h) {
        next = h;
        return ++n < 2;
    });
    if (n != 1 || graph.get_is_reverse(next) || graph.get_id(next) == start) {
        return false;
    }
    size_t in = 0;
    graph.follow_edges(next, forward, [&](const handle_t&) { return ++in < 2; });
    return in == 1;
}

/// Return the nodes met from `from`, read forward, going forward (or backward) along the only edge
/// out of each node into a node with no other edge in, until they hold `limit` bases. Order them
/// as a forward walk reads them.
vector<handle_t> unbranched_walk(const HandleGraph& graph, handle_t from, bool forward,
                                 size_t limit) {
    vector<handle_t> out{from};
    size_t total = graph.get_length(from);
    handle_t at = from;
    handle_t next = from;
    while (total < limit && sole_successor(graph, at, forward, graph.get_id(from), next)) {
        out.push_back(next);
        total += graph.get_length(next);
        at = next;
    }
    if (!forward) {
        std::ranges::reverse(out);
    }
    return out;
}

/// Make in `out` the key for `indel`, given in the coordinates of the walk `nodes` read forward,
/// once left-aligned along it. False if its tract runs off the end of the walk.
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
        return (size_t)(std::ranges::upper_bound(starts, x) - starts.begin()) - 1;
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

/// Return the sorted 8-mer codes of `bases`, for a cheap distance between two reads' insertions.
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
    std::ranges::sort(out);
    return out;
}

/// Return the size of the symmetric multiset difference of two sorted profiles.
size_t profile_distance(const vector<uint32_t>& a, const vector<uint32_t>& b) {
    size_t i = 0;
    size_t j = 0;
    size_t shared = 0;
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

/// Put in `difference` how far the anchor of `key` lies right of the anchor of `other`.
///
/// It is measured along the walk of `other` if that holds the anchor node of `key`, or else along
/// the walk of `key` if that holds the anchor node of `other`. False if neither does.
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

/// Return the number of distinct fragments in `f`.
size_t distinct_fragments(vector<uint64_t> f) {
    std::ranges::sort(f);
    return (size_t)(std::ranges::unique(f).begin() - f.begin());
}

/// Return the index into `walked` of the medoid of the rows of `rows` it lists, all of kind
/// `kind`: the row whose indel is nearest all the others'.
///
/// Insertions are compared by the k-mer profiles of their bases, and deletions by anchor and
/// length. Ties go to the first.
size_t medoid_of(const vector<LargeRow>& rows, const vector<size_t>& walked, char kind) {
    vector<vector<uint32_t>> profiles;
    if (kind == 'I') {
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
            total += kind == 'I' ? profile_distance(profiles[x], profiles[y])
                                 : (size_t)(std::abs((int64_t)p.offset - (int64_t)q.offset)
                                            + std::abs((int64_t)p.length - (int64_t)q.length));
        }
        if (total < best_total) {
            best_total = total;
            best = x;
        }
    }
    return best;
}

/// Add to `keys` the key for one cluster of read indels, the rows `members` of `rows`, if at least
/// `min_fragments` fragments carry it. The key is the medoid's indel, chosen among the first 40
/// members whose reads reach the end of the tract.
void add_cluster_key(const HandleGraph& graph, const vector<LargeRow>& rows,
                     const vector<size_t>& members, size_t min_fragments, vector<LargeKey>& keys) {
    vector<uint64_t> fragments;
    vector<size_t> walked;
    for (size_t k : members) {
        fragments.push_back(rows[k].fragment);
        if (!rows[k].walk.empty() && walked.size() < 40) {
            walked.push_back(k);
        }
    }
    const size_t support = distinct_fragments(std::move(fragments));
    if (support < min_fragments || walked.empty()) {
        return;
    }
    const LargeRow& medoid = rows[walked[medoid_of(rows, walked, rows[members[0]].kind)]];
    WalkIndel indel;
    indel.at = (int64_t)medoid.offset + 1;
    if (medoid.kind == 'D') {
        indel.deleted = medoid.length;
    } else {
        indel.inserted = medoid.inserted;
    }
    if (LargeKey key; make_large_key(graph, medoid.walk, indel, support, false, key)) {
        keys.push_back(std::move(key));
    }
}

/// Cluster the rows `rows[i, j)`, one run of anchors, by length, and add a key for each cluster
/// to `keys`. In a cluster, each indel's length is within `length_ratio` of the next shorter one's.
void add_length_cluster_keys(const HandleGraph& graph, const vector<LargeRow>& rows, size_t i,
                             size_t j, const EditCountParams& params, vector<LargeKey>& keys) {
    vector<size_t> by_length;
    for (size_t k = i; k < j; ++k) {
        by_length.push_back(k);
    }
    std::ranges::stable_sort(by_length,
                             [&](size_t a, size_t b) { return rows[a].length < rows[b].length; });
    const double ratio = params.large.length_ratio;
    for (size_t a = 0; a < by_length.size();) {
        size_t b = a + 1;
        while (b < by_length.size()
               && (double)rows[by_length[b]].length
                      <= (double)rows[by_length[b - 1]].length * (1.0 + ratio)) {
            ++b;
        }
        const vector<size_t> members(by_length.begin() + (int64_t)a,
                                     by_length.begin() + (int64_t)b);
        add_cluster_key(graph, rows, members, params.min_fragments, keys);
        a = b;
    }
}

/// Return keys for the reads' own large indels in `rows`, which this sorts. The rows are clustered
/// on one node by anchor, each within `slop` of the next, then by length.
vector<LargeKey> read_indel_keys(const HandleGraph& graph, vector<LargeRow>& rows,
                                 const EditCountParams& params) {
    vector<LargeKey> keys;
    std::ranges::sort(rows);
    for (size_t i = 0; i < rows.size();) {
        size_t j = i + 1;
        while (j < rows.size() && rows[j].kind == rows[i].kind && rows[j].node == rows[i].node
               && (int64_t)rows[j].offset - (int64_t)rows[j - 1].offset <= params.large.slop) {
            ++j;
        }
        add_length_cluster_keys(graph, rows, i, j, params, keys);
        i = j;
    }
    return keys;
}

/// Return the index of the clip in `clips[i, j)`, one clip position sorted by seed, whose clipped
/// bases the most fragments carry, and put that number of fragments in `support`. Ties go to the
/// first.
size_t commonest_seed(const vector<ClipRow>& clips, size_t i, size_t j, size_t& support) {
    support = 0;
    size_t best = i;
    for (size_t k = i; k < j;) {
        size_t l = k;
        vector<uint64_t> fragments;
        while (l < j && clips[l].seed == clips[k].seed) {
            fragments.push_back(clips[l].fragment);
            ++l;
        }
        if (const size_t found = distinct_fragments(std::move(fragments)); found > support) {
            support = found;
            best = k;
        }
        k = l;
    }
    return best;
}

/// Return the upper-cased bases of `handles`, and put their node IDs in `nodes`.
string spell_handles(const HandleGraph& graph, const vector<handle_t>& handles,
                     vector<nid_t>& nodes) {
    string spelling;
    for (const handle_t& n : handles) {
        nodes.push_back(graph.get_id(n));
        spelling += graph.get_sequence(n);
    }
    return upper(std::move(spelling));
}

/// Find the deletion that `clip` shows, its clipped bases found again along the unbranched walk
/// from the clip's node, at least `LARGE_INDEL_MIN_LENGTH` and at most `max_clip_deletion` bases
/// away. Put the walk's nodes in `nodes` and the deletion, in the walk's coordinates, in `indel`.
/// False if there is none.
bool deletion_from_clip(const HandleGraph& graph, const ClipRow& clip,
                        const LargeIndelParams& large, vector<nid_t>& nodes, WalkIndel& indel) {
    const handle_t h = graph.get_handle(clip.node);
    const size_t seed = clip.seed.size();
    const size_t reach = (size_t)large.max_clip_deletion + seed + 1;
    if (clip.side == 'R') {
        // The deleted bases follow the clip's base, and the seed comes after them.
        const string spelling =
            spell_handles(graph, unbranched_walk(graph, h, true, clip.offset + reach), nodes);
        const size_t from = (size_t)clip.offset + 1 + LARGE_INDEL_MIN_LENGTH;
        const size_t q = spelling.find(clip.seed, from);
        if (q == string::npos || q - clip.offset - 1 > large.max_clip_deletion) {
            return false;
        }
        indel.at = (int64_t)clip.offset + 1;
        indel.deleted = (uint32_t)(q - clip.offset - 1);
        return true;
    }
    // The deleted bases precede the clip's base, and the seed comes before them.
    const string spelling = spell_handles(
        graph, unbranched_walk(graph, h, false, reach + graph.get_length(h)), nodes);
    const int64_t q = (int64_t)(spelling.size() - graph.get_length(h)) + clip.offset;
    const int64_t last_start = q - (int64_t)LARGE_INDEL_MIN_LENGTH - (int64_t)seed;
    if (last_start < 1) {
        return false;
    }
    const size_t s = spelling.rfind(clip.seed, (size_t)last_start);
    if (s == string::npos || s < 1 || q - (int64_t)(s + seed) > (int64_t)large.max_clip_deletion) {
        return false;
    }
    indel.at = (int64_t)(s + seed);
    indel.deleted = (uint32_t)(q - indel.at);
    return true;
}

/// Return keys for the deletions that the clips in `clips`, which this sorts, show. Each clip
/// position that enough fragments share gives the deletion its commonest clipped bases show.
vector<LargeKey> clip_deletion_keys(const HandleGraph& graph, vector<ClipRow>& clips,
                                    const EditCountParams& params) {
    vector<LargeKey> keys;
    std::ranges::sort(clips);
    for (size_t i = 0; i < clips.size();) {
        size_t j = i;
        while (j < clips.size() && clips[j].side == clips[i].side && clips[j].node == clips[i].node
               && clips[j].offset == clips[i].offset) {
            ++j;
        }
        size_t support = 0;
        const ClipRow& clip = clips[commonest_seed(clips, i, j, support)];
        i = j;
        if (support < params.min_fragments || !graph.has_node(clip.node)) {
            continue;
        }
        vector<nid_t> nodes;
        WalkIndel indel;
        if (!deletion_from_clip(graph, clip, params.large, nodes, indel)) {
            continue;
        }
        if (LargeKey key; make_large_key(graph, nodes, indel, support, true, key)) {
            keys.push_back(std::move(key));
        }
    }
    return keys;
}

/// Whether `key` makes the same change as `other`: an indel of the same kind, its anchor within
/// `slop` of where `other`'s tract lets `other` slide, and its length within `length_ratio` of
/// `other`'s.
bool same_large_change(const LargeKey& key, const LargeKey& other, const LargeIndelParams& large) {
    const int64_t room = std::max<int64_t>(
        0, (int64_t)other.span - 1 - (other.kind == 'D' ? (int64_t)other.length : 0));
    int64_t difference = 0;
    return other.kind == key.kind && anchor_difference(key, other, difference)
           && std::abs(difference) <= large.slop + room
           && std::abs((double)key.length - (double)other.length)
                  <= large.length_ratio * (double)other.length;
}

/// Return one key per change among `keys`, leaving out those anchored on nodes `owns` refuses,
/// sorted by anchor. Of keys making the same change, the reads' own indels win over clips, and
/// then the best supported wins.
template<typename Owns>
vector<LargeKey> distinct_large_keys(const vector<LargeKey>& keys, const Owns& owns,
                                     const LargeIndelParams& large) {
    vector<size_t> order(keys.size());
    for (size_t k = 0; k < keys.size(); ++k) {
        order[k] = k;
    }
    std::ranges::stable_sort(order, [&](size_t a, size_t b) {
        return tie(keys[a].clips_only, keys[b].support) < tie(keys[b].clips_only, keys[a].support);
    });
    vector<LargeKey> chosen;
    for (size_t k : order) {
        const LargeKey& key = keys[k];
        const bool duplicate = !owns(key.node)
                               || std::ranges::any_of(chosen, [&](const LargeKey& other) {
                                      return same_large_change(key, other, large);
                                  });
        if (!duplicate) {
            chosen.push_back(key);
        }
    }
    std::ranges::sort(chosen, [](const LargeKey& a, const LargeKey& b) {
        return tie(a.node, a.offset, a.kind, a.length, a.inserted)
               < tie(b.node, b.offset, b.kind, b.length, b.inserted);
    });
    return chosen;
}

/// A tally of what the reads say about each large indel, by fragment, as `observe_large_indel`
/// matches them to it.
///
/// A deletion that no read carries whole must have reads clipped at each of its ends.
class LargeTally {
public:
    /// Tally reads against `keys`, sorted by anchor, matched as `large` says.
    LargeTally(vector<LargeKey> keys, const LargeIndelParams& large)
        : keys(std::move(keys)), by_key(this->keys.size()), large(large) {
        for (size_t k = 0; k < this->keys.size(); ++k) {
            for (size_t j = 0; j < this->keys[k].walk.size(); ++j) {
                on_node[this->keys[k].walk[j]].emplace_back(k, j);
            }
        }
    }

    /// Whether the read visits a node of some key's walk.
    bool near(const Alignment& aln) const {
        return std::ranges::any_of(aln.path().mapping(), [&](const Mapping& m) {
            return on_node.contains(m.position().node_id());
        });
    }

    /// Record what the loaded `walk`, from `fragment` at mapping quality `mapq`, says about each
    /// key whose walk it visits forward, once per key.
    void observe(const HandleGraph& graph, const ReadWalk& walk, uint64_t fragment, int32_t mapq) {
        for (size_t k : touched) {
            by_key[k].done = false;
        }
        touched.clear();
        for (size_t i = 0; i + 1 < walk.node_start().size(); ++i) {
            const handle_t h = walk.handles()[i];
            auto found = on_node.find(graph.get_id(h));
            if (found == on_node.end() || graph.get_is_reverse(h)) {
                continue;
            }
            for (const auto& [k, j] : found->second) {
                if (by_key[k].done) {
                    continue;
                }
                by_key[k].done = true;
                touched.push_back(k);
                observe_key(graph, walk, i, k, j, fragment, mapq);
            }
        }
    }

    /// Return the candidates that pass the floors of `params` and that no other walk of the graph
    /// spells.
    vector<EditCandidate> candidates(const HandleGraph& graph, const EditCountParams& params) {
        vector<EditCandidate> kept;
        std::ranges::sort(seen);
        for (size_t i = 0; i < seen.size();) {
            const size_t k = seen[i].key;
            size_t end = i;
            while (end < seen.size() && seen[end].key == k) {
                ++end;
            }
            EditCandidate c;
            count_fragments(seen, i, end, c);
            i = end;
            const KeyTally& tally = by_key[k];
            if (c.alt_fragments < params.min_fragments || c.fraction() < params.min_fraction) {
                continue;
            }
            if (tally.whole == 0
                && (distinct_fragments(tally.clipped_left) < large.min_fragments_per_end
                    || distinct_fragments(tally.clipped_right) < large.min_fragments_per_end)) {
                continue;
            }
            fill_candidate(c, keys[k]);
            if (spelled_by_panel(graph, c)) {
                continue;
            }
            c.mean_alt_mapq = tally.mapq_count == 0
                                  ? 0.0f
                                  : (float)(tally.mapq_total / (double)tally.mapq_count);
            kept.push_back(std::move(c));
        }
        return kept;
    }

private:
    /// What the reads have said about one key.
    struct KeyTally {
        /// The fragments with an ALT read clipped at the deletion's left end.
        vector<uint64_t> clipped_left;
        /// The fragments with an ALT read clipped at the deletion's right end.
        vector<uint64_t> clipped_right;
        /// The ALT reads that carry the indel whole.
        size_t whole = 0;
        double mapq_total = 0.0;
        size_t mapq_count = 0;
        /// Whether the current read has been matched to the key.
        bool done = false;
    };

    /// Whether `walk` takes the walk of `key` wherever the two overlap, where its mapping `i`
    /// is on node `j` of the key's walk.
    static bool follows_key_walk(const HandleGraph& graph, const ReadWalk& walk, size_t i,
                                 const LargeKey& key, size_t j) {
        bool agrees = true;
        for (size_t t = 0; t < key.walk.size() && agrees; ++t) {
            const int64_t m = (int64_t)i + (int64_t)t - (int64_t)j;
            if (m >= 0 && m + 1 < (int64_t)walk.node_start().size()) {
                agrees = walk.handles()[m] == graph.get_handle(key.walk[t]);
            }
        }
        return agrees;
    }

    /// Record what `walk` says about key number `k`, where its mapping `i` is on node `j` of the
    /// key's walk.
    void observe_key(const HandleGraph& graph, const ReadWalk& walk, size_t i, size_t k, size_t j,
                     uint64_t fragment, int32_t mapq) {
        const LargeKey& key = keys[k];
        if (!follows_key_walk(graph, walk, i, key, j)) {
            return;
        }
        const int64_t anchor = walk.node_start()[i] - key.starts[j] + key.offset;
        WalkIndel wanted;
        wanted.at = anchor + 1;
        if (key.kind == 'D') {
            wanted.deleted = key.length;
        } else {
            wanted.inserted = key.inserted;
        }
        int clip = 0;
        const IndelObservation observed =
            observe_large_indel(walk.alignment(), wanted, anchor + key.span, large,
                                walk.alignment().begin, walk.alignment().end, &clip);
        if (observed == IndelObservation::NONE) {
            return;
        }
        seen.emplace_back(k, fragment, (uint8_t)observed);
        if (observed != IndelObservation::ALT) {
            return;
        }
        KeyTally& tally = by_key[k];
        if (clip == 1) {
            tally.clipped_left.push_back(fragment);
        } else if (clip == 2) {
            tally.clipped_right.push_back(fragment);
        } else {
            ++tally.whole;
        }
        tally.mapq_total += mapq;
        ++tally.mapq_count;
    }

    /// Fill in the change `c` makes, from its key.
    static void fill_candidate(EditCandidate& c, const LargeKey& key) {
        c.node = key.node;
        c.offset = key.offset;
        c.kind = key.kind;
        c.alt = 0;
        c.ref = key.spelling[key.offset];
        c.span = key.span;
        c.sequence = key.kind == 'D' ? key.spelling.substr(key.offset + 1, key.length)
                                     : key.inserted;
        c.walk = key.walk;
        for (size_t x = key.offset + 1;
             x < key.spelling.size() && key.spelling[x] == c.sequence[0]; ++x) {
            ++c.homopolymer;
        }
    }

    /// The keys, sorted by anchor.
    vector<LargeKey> keys;
    /// The tally for each of `keys`.
    vector<KeyTally> by_key;
    const LargeIndelParams& large;
    /// For each node, each key whose walk visits it, and where along the walk.
    unordered_map<nid_t, vector<pair<size_t, size_t>>> on_node;
    /// What each read said about each key.
    vector<Seen> seen;
    /// The keys the current read has been matched to.
    vector<size_t> touched;
};

/// Count the large indels among the reads, from the rows the short indel count filed.
///
/// The keys come from each read's own large indels (`rows`), clustered by anchor and length with
/// the medoid read's indel taken as the ALT, and from deletions that clips show (`clips`). Each is
/// then counted as `observe_large_indel` matches the reads to it.
vector<EditCandidate> count_large_indel_edits(const HandleGraph& graph,
                                              const ReadSource& for_each_read,
                                              const function<bool(nid_t)>& owns,
                                              const EditCountParams& params, vector<LargeRow> rows,
                                              vector<ClipRow> clips) {
    vector<LargeKey> keys = read_indel_keys(graph, rows, params);
    vector<LargeKey> from_clips = clip_deletion_keys(graph, clips, params);
    keys.insert(keys.end(), std::make_move_iterator(from_clips.begin()),
                std::make_move_iterator(from_clips.end()));
    keys = distinct_large_keys(keys, owns, params.large);
    if (keys.empty()) {
        return {};
    }

    LargeTally tally(std::move(keys), params.large);
    ReadWalk walk(graph);
    std::hash<string> name_hash;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq || !tally.near(aln) || !walk.load(aln)) {
            return;
        }
        tally.observe(graph, walk, name_hash(aln.name()), aln.mapping_quality());
    });
    return tally.candidates(graph, params);
}

/// A run of admitted candidates on one node whose changed bases touch or overlap, to be merged.
struct TouchGroup {
    /// The indices of its candidates in the candidate list.
    vector<size_t> members;
    /// The changed bases [lo, hi) along the node.
    uint32_t lo;
    uint32_t hi;
    /// The longest walk of its candidates, which starts at the node.
    vector<nid_t> walk;
};

/// Return the groups of the candidates in `kept` that `bar` admits, large indels aside, that are
/// worth merging, sorted by node.
///
/// A group is worth merging when it has more than one candidate and they are not all indels at
/// one anchor.
vector<TouchGroup> touching_groups(const vector<EditCandidate>& kept, const EditCallBar& bar) {
    // The changed bases [from, to) along each admitted candidate's node, by node.
    vector<tuple<nid_t, uint32_t, uint32_t, size_t>> changes;
    for (size_t i = 0; i < kept.size(); ++i) {
        const EditCandidate& c = kept[i];
        if (bar.admits(c) && !c.is_large()) {
            const uint32_t from = c.kind == 'S' ? c.offset : c.offset + 1;
            changes.emplace_back(c.node, from, c.kind == 'S' ? from + 1 : c.offset + c.span, i);
        }
    }
    std::ranges::sort(changes);
    vector<TouchGroup> groups;
    for (size_t i = 0, j = 0; i < changes.size(); i = j) {
        TouchGroup g{{}, get<1>(changes[i]), get<2>(changes[i]), {get<0>(changes[i])}};
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
    return groups;
}

/// Put in `spelling` what the loaded `walk` spells for group `g`, whose walk it follows from its
/// mapping `i`: the read's bases between the bases on either side of the group's bases. False if
/// the read does not follow the group's walk, does not align both of those bases, or has a base
/// below quality `min_q` there or at a substitution of `kept`.
bool group_spelling(const HandleGraph& graph, const ReadWalk& walk, size_t i, const TouchGroup& g,
                    const vector<EditCandidate>& kept, int min_q, string& spelling) {
    bool follows = i + g.walk.size() <= walk.handles().size();
    for (size_t s = 0; follows && s < g.walk.size(); ++s) {
        follows = walk.handles()[i + s] == graph.get_handle(g.walk[s]);
    }
    const int64_t before = walk.node_start()[i] + g.lo - 1;
    const int64_t after = walk.node_start()[i] + g.hi;
    const int64_t r0 = follows ? walk.alignment().read_offset_at(before) : -1;
    const int64_t r1 = follows ? walk.alignment().read_offset_at(after) : -1;
    bool good = r0 >= 0 && r1 > r0 && walk.quality_ok(before, min_q)
                && walk.quality_ok(after, min_q);
    for (size_t m : g.members) {
        good = good && (kept[m].kind != 'S'
                        || walk.quality_ok(walk.node_start()[i] + kept[m].offset, min_q));
    }
    if (good) {
        spelling = upper(walk.read().sequence().substr(r0 + 1, r1 - r0 - 1));
    }
    return good;
}

/// What a read spells for a group: the group's index, the read's fragment, and the bases.
using GroupSpelling = tuple<size_t, uint64_t, string>;

/// Return what each read following a group's walk spells for it, sorted and without duplicates.
vector<GroupSpelling> collect_group_spellings(const HandleGraph& graph,
                                              const ReadSource& for_each_read,
                                              const EditCountParams& params,
                                              const vector<EditCandidate>& kept,
                                              const vector<TouchGroup>& groups) {
    vector<GroupSpelling> spelled;
    std::hash<string> name_hash;
    ReadWalk walk(graph);
    unordered_set<nid_t> nodes;
    for (const TouchGroup& g : groups) {
        nodes.insert(g.walk[0]);
    }
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq
            || std::ranges::none_of(aln.path().mapping(),
                                    [&](const Mapping& m) {
                                        return nodes.contains(m.position().node_id());
                                    })
            || !walk.load(aln)) {
            return;
        }
        for (size_t i = 0; i < walk.handles().size(); ++i) {
            const nid_t id = graph.get_id(walk.handles()[i]);
            for (auto k = std::ranges::lower_bound(
                     groups, id, {}, [](const TouchGroup& g) { return g.walk[0]; });
                 k != groups.end() && k->walk[0] == id; ++k) {
                if (string spelling;
                    group_spelling(graph, walk, i, *k, kept, params.min_base_quality, spelling)) {
                    spelled.emplace_back((size_t)(k - groups.begin()), name_hash(aln.name()),
                                         std::move(spelling));
                }
            }
        }
    });
    std::ranges::sort(spelled);
    const auto duplicates = std::ranges::unique(spelled);
    spelled.erase(duplicates.begin(), duplicates.end());
    return spelled;
}

/// The fragments' spellings of one group.
struct SpellingCount {
    /// The fragments giving each spelling, each counted under the one spelling all its reads give.
    map<string, uint32_t> fragments;
    /// The fragments counted in `fragments`.
    uint32_t total = 0;
    /// The fragments whose reads give more than one spelling.
    uint32_t discordant = 0;
};

/// Count the fragments' spellings in the run of `spelled` from `i` that is for one group, into
/// `count`. Return the end of the run.
size_t count_spellings(const vector<GroupSpelling>& spelled, size_t i, SpellingCount& count) {
    size_t j = i;
    for (size_t f = i; f < spelled.size() && get<0>(spelled[f]) == get<0>(spelled[i]); f = j) {
        j = f + 1;
        while (j < spelled.size() && get<0>(spelled[j]) == get<0>(spelled[f])
               && get<1>(spelled[j]) == get<1>(spelled[f])) {
            ++j;
        }
        if (j > f + 1) {
            ++count.discordant;
        } else {
            ++count.fragments[get<2>(spelled[f])];
            ++count.total;
        }
    }
    return j;
}

/// Return the merged candidate for group `g`, whose walk spells `bases`, from its fragments'
/// spellings. Its ALTs are the spellings `bar` admits as a candidate of their own, though not
/// checked for a homopolymer, and changing the length only between two walk bases.
EditCandidate merged_candidate(const TouchGroup& g, const string& bases,
                               const SpellingCount& count, const EditCallBar& bar) {
    EditCandidate m;
    m.discordant_fragments = count.discordant;
    m.node = g.walk[0];
    m.offset = g.lo;
    m.kind = 'M';
    m.span = g.hi - g.lo;
    m.walk = g.walk;
    m.sequence = upper(bases.substr(g.lo, m.span));
    const auto ref = count.fragments.find(m.sequence);
    m.ref_fragments = ref != count.fragments.end() ? ref->second : 0;
    for (const auto& [spelling, fragments] : count.fragments) {
        EditCandidate alone;
        alone.kind = spelling.size() == m.span ? 'S' : 'I';
        alone.alt_fragments = fragments;
        alone.ref_fragments = m.ref_fragments;
        alone.other_fragments = count.total - fragments - m.ref_fragments;
        if (spelling != m.sequence && bar.admits(alone)
            && (alone.kind == 'S' || (g.lo > 0 && g.hi < bases.size()))) {
            m.alts.push_back(spelling);
            m.alt_fragments += fragments;
        }
    }
    m.other_fragments = count.total - m.ref_fragments - m.alt_fragments;
    return m;
}

/// Whether the ALT `a` of a merged change of `span` bases from `lo` to `hi` carries the member
/// candidate `e`. It carries a substitution if it has its base at its place counted from either
/// end, and an indel if it changes the length by the indel's.
bool alt_carries(const string& a, const EditCandidate& e, uint32_t lo, uint32_t hi, uint32_t span) {
    if (e.kind == 'S') {
        const size_t left = e.offset - lo;
        const size_t right = hi - e.offset;
        return (left < a.size() && a[left] == e.alt)
               || (right <= a.size() && a[a.size() - right] == e.alt);
    }
    const int64_t change = (e.kind == 'I' ? 1 : -1) * (int64_t)e.sequence.size();
    return (int64_t)a.size() - (int64_t)span == change;
}

/// Merge `kept`'s touching candidates, as `count_edits` says, keeping `kept` sorted.
void merge_touching(const HandleGraph& graph, const ReadSource& for_each_read,
                    const EditCountParams& params, vector<EditCandidate>& kept) {
    const vector<TouchGroup> groups = touching_groups(kept, *params.bar);
    if (groups.empty()) {
        return;
    }
    const vector<GroupSpelling> spelled =
        collect_group_spellings(graph, for_each_read, params, kept, groups);

    vector<EditCandidate> merged;
    for (size_t i = 0; i < spelled.size();) {
        const TouchGroup& g = groups[get<0>(spelled[i])];
        SpellingCount count;
        i = count_spellings(spelled, i, count);
        string bases;
        for (nid_t id : g.walk) {
            bases += graph.get_sequence(graph.get_handle(id));
        }
        EditCandidate m = merged_candidate(g, bases, count, *params.bar);
        // Merge only where every member is in some ALT.
        bool carried = !m.alts.empty();
        for (size_t c : g.members) {
            const EditCandidate& e = kept[c];
            carried = carried && std::ranges::any_of(m.alts, [&](const string& a) {
                return alt_carries(a, e, g.lo, g.hi, m.span);
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
    // Put the merged candidates in place of their members, which are now of kind 0.
    std::erase_if(kept, [](const EditCandidate& c) { return c.kind == 0; });
    kept.insert(kept.end(), std::make_move_iterator(merged.begin()),
                std::make_move_iterator(merged.end()));
    std::ranges::sort(kept);
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
        std::ranges::sort(kept);
    }
    if (params.bar != nullptr) {
        merge_touching(graph, for_each_read, params, kept);
    }
    return kept;
}

int64_t WalkAlignment::read_offset_at(int64_t x) const {
    auto it = std::ranges::upper_bound(blocks, x, {}, &Block::walk);
    if (it == blocks.begin()) {
        return -1;
    }
    --it;
    if (x >= it->walk + (int64_t)it->length) {
        return -1;
    }
    return (int64_t)it->read + (x - it->walk);
}

namespace {

/// A builder of a read's `WalkAlignment`, fed the read's edits in walk order.
///
/// Indels wait for an aligned base after them, so that the insertion of a trailing clip is left
/// out. Indels before the first aligned base are left out too, as part of the leading clip.
class WalkAlignmentBuilder {
public:
    /// Start a read whose first edit's first base is at read offset `read_offset` in `sequence`.
    WalkAlignmentBuilder(size_t read_offset, string_view sequence)
        : read_offset(read_offset), sequence(sequence), read_pos(read_offset),
          read_after_aligned(read_offset) {}

    /// Add `length` aligned bases from walk base `walk_pos`.
    void aligned(int64_t walk_pos, size_t length) {
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
        if (!out.blocks.empty()
            && out.blocks.back().walk + (int64_t)out.blocks.back().length == walk_pos
            && out.blocks.back().read + out.blocks.back().length == read_pos) {
            out.blocks.back().length += length;
        } else {
            out.blocks.emplace_back(walk_pos, read_pos, length);
        }
        out.end = walk_pos + (int64_t)length;
        read_pos += length;
    }

    /// Add a deletion of the `length` walk bases from walk base `walk_pos`, stitched onto a
    /// pending deletion that ends there.
    void deleted(int64_t walk_pos, size_t length) {
        if (!started) {
            return;
        }
        if (!pending.empty() && pending.back().deleted > 0
            && pending.back().at + (int64_t)pending.back().deleted == walk_pos) {
            pending.back().deleted += (uint32_t)length;
        } else {
            WalkIndel d;
            d.at = walk_pos;
            d.deleted = (uint32_t)length;
            pending.push_back(std::move(d));
        }
    }

    /// Add an insertion of the next `length` read bases before walk base `walk_pos`, joined onto a
    /// pending insertion there.
    void inserted(int64_t walk_pos, size_t length) {
        if (started) {
            // Bases past the end of the sequence are left empty.
            string bases = read_pos + length <= sequence.size()
                               ? string(sequence.substr(read_pos, length))
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
        read_pos += length;
    }

    /// Return the alignment, with its trailing clip measured.
    WalkAlignment finish() {
        if (started) {
            out.clip_end = read_pos - read_after_aligned;
        }
        return std::move(out);
    }

private:
    const size_t read_offset;
    const string_view sequence;
    WalkAlignment out;
    /// The read offset of the next edit's first base.
    size_t read_pos;
    /// Whether an aligned base has been added.
    bool started = false;
    /// The indels since the last aligned base.
    vector<WalkIndel> pending;
    /// The read offset after the last aligned base.
    size_t read_after_aligned;
};

}

WalkAlignment walk_alignment(const vector<const Mapping*>& mappings,
                             const vector<int64_t>& node_start, size_t read_offset,
                             string_view sequence) {
    WalkAlignmentBuilder builder(read_offset, sequence);
    for (size_t i = 0; i < mappings.size(); ++i) {
        const Mapping& m = *mappings[i];
        int64_t walk_pos = node_start[i] + (int64_t)m.position().offset();
        for (const Edit& e : m.edit()) {
            auto from = (size_t)e.from_length();
            auto to = (size_t)e.to_length();
            if (const size_t both = std::min(from, to); both > 0) {
                builder.aligned(walk_pos, both);
                walk_pos += (int64_t)both;
                from -= both;
                to -= both;
            }
            if (from > 0) {
                builder.deleted(walk_pos, from);
                walk_pos += (int64_t)from;
            } else if (to > 0) {
                builder.inserted(walk_pos, to);
            }
        }
    }
    return builder.finish();
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
    // Find the first indel ending after the anchor. The indels are in walk order and do not
    // overlap, so their ends are in order too.
    auto it = std::ranges::lower_bound(aln.indels, anchor + 1, {}, [](const WalkIndel& d) {
        return d.at + (int64_t)std::max<uint32_t>(d.deleted, 1);
    });
    size_t count = 0;
    const WalkIndel* only = nullptr;
    for (; it != aln.indels.end() && it->at <= after; ++it) {
        const bool touches = it->deleted > 0 ? it->at + (int64_t)it->deleted > anchor
                                             : it->at > anchor;
        if (touches) {
            ++count;
            only = std::to_address(it);
            if (touching != nullptr) {
                touching->push_back(only);
            }
        }
    }
    if (count == 0) {
        return IndelObservation::REF;
    }
    return count == 1 && *only == indel ? IndelObservation::ALT : IndelObservation::OTHER;
}

namespace {

/// Return the number of bases of the indels of `aln` that are of the kind `deletion` says and lie
/// near the walk bases [lo, hi]. A deletion is near if it ends at or after `lo`, and an insertion
/// if it is at or after `lo`, and either if it starts no later than `hi`.
int64_t near_indel_bases(const WalkAlignment& aln, bool deletion, int64_t lo, int64_t hi) {
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
    return near;
}

/// Return which end of a deletion the read `aln` is clipped at: 1 if its aligned bases end at the
/// deletion's left end `left`, 2 if they start at its right end `right`, and 0 if neither. The
/// clip must be `min_clip` bases or more, and the aligned bases must meet the end within `slop`,
/// or `room` bases further right, as far as the deletion's tract lets it slide.
int deletion_clip_end(const WalkAlignment& aln, int64_t left, int64_t right, int64_t room,
                      const LargeIndelParams& params) {
    if (aln.clip_end >= params.min_clip && aln.end >= left - params.slop
        && aln.end <= left + room + params.slop) {
        return 1;
    }
    if (aln.clip_begin >= params.min_clip && aln.begin >= right - params.slop
        && aln.begin <= right + room + params.slop) {
        return 2;
    }
    return 0;
}

}

IndelObservation observe_large_indel(const WalkAlignment& aln, const WalkIndel& indel,
                                     int64_t after, const LargeIndelParams& params,
                                     int64_t reach_begin, int64_t reach_end, int* by_clip) {
    if (by_clip != nullptr) {
        *by_clip = 0;
    }
    const bool deletion = indel.deleted > 0;
    const auto length = (int64_t)indel.length();
    // The change's two ends, the bases after the anchor and after the deleted bases, and how far
    // its tract lets it slide right.
    const int64_t left = indel.at;
    const int64_t right = deletion ? indel.at + length : indel.at;
    const int64_t room = std::max<int64_t>(0, after - right);
    const int64_t lo = left - 1 - params.slop;
    const int64_t hi = after + params.slop;
    const int64_t near = near_indel_bases(aln, deletion, lo, hi);
    const double tolerance = params.length_ratio * (double)length;
    if (std::abs((double)(near - length)) <= tolerance) {
        return IndelObservation::ALT;
    }
    // A deletion's read can be clipped where its aligned bases meet one end, the rest of the read
    // being the bases past the other end.
    if (deletion) {
        if (const int end = deletion_clip_end(aln, left, right, room, params); end != 0) {
            if (by_clip != nullptr) {
                *by_clip = end;
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
    const bool clip_before = aln.clip_end > 0 && aln.end >= lo && aln.end <= hi;
    const bool clip_after = aln.clip_begin > 0 && aln.begin >= lo && aln.begin <= hi;
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
    std::ranges::sort(all, [](const EditOutcome& a, const EditOutcome& b) {
        return tie(a.node, a.offset, a.kind, a.alt, a.sequence)
               < tie(b.node, b.offset, b.kind, b.alt, b.sequence);
    });
    return all;
}

void for_each_candidate_on(const vector<EditCandidate>& sorted, const vector<nid_t>& nodes,
                           const function<void(const EditCandidate&)>& iteratee) {
    auto it = sorted.begin();
    for (nid_t node : nodes) {
        it = std::ranges::lower_bound(it, sorted.end(), node, {}, &EditCandidate::node);
        for (; it != sorted.end() && it->node == node; ++it) {
            iteratee(*it);
        }
    }
}

}
