#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <tuple>
#include <unordered_set>

#include <omp.h>

#include "mosaic_writer.hpp"
#include "allele_splicer.hpp"
#include "../src/sorted_vcf_writer.hpp"
#include "../src/graph_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "../src/algorithms/expand_context.hpp"
#include "../src/annotation.hpp"
#include "../src/gref.hpp"
#include "../src/traversal_clusters.hpp"
#include "../src/utility.hpp"

namespace vg::multipass {

namespace {

using PhaseCall = LinkageCollector::PhaseCall;

/// How far a walk may run before it is abandoned. Walks go only in a direction already
/// established, so this limits a long run rather than a wrong-way search.
const size_t WALK_LIMIT = 1u << 17;

/// Finds and follows panel haplotypes in the GBWT, for one contig of a mosaic. Each contig is
/// written by one thread with its own walker, so its caches need no locking.
class PanelWalker {
public:
    explicit PanelWalker(const PanelLookup& panel)
        : index(panel.gbwt()), haplotype_of_sequence(panel.sequence_to_haplotype()) {}

    /// The GBWT, or null when there is no panel.
    const gbwt::GBWT* gbwt() const { return index; }

    /// Where `hap` sits at one oriented node, or invalid.
    gbwt::edge_type position_at(gbwt::node_type node, size_t hap) {
        if (index == nullptr || haplotype_of_sequence == nullptr) {
            return gbwt::invalid_edge();
        }
        // Finding a position costs a `locate` for each sequence in the node's range, far more
        // than an LF step, and the same (node, haplotype) is asked for repeatedly, so positions
        // are cached. A node's range is scanned once in all: a search goes on from where the last
        // one at the node stopped, and caches every haplotype it passes at that haplotype's first
        // position, which is the one a search for it would have returned.
        if (hap >= POSITION_KEY_HAPLOTYPES) {
            return search_uncached(node, hap);
        }
        const uint64_t key = position_key(node, hap);
        if (auto hit = position_cache.find(key); hit != position_cache.end()) {
            return hit->second;
        }
        Scan& scan = scans[node];
        if (!scan.started) {
            gbwt::SearchState state = index->find(node);
            scan.started = true;
            scan.next = state.range.first;
            scan.end = state.empty() ? state.range.first : state.range.second + 1;
        }
        while (scan.next < scan.end) {
            const gbwt::edge_type here(node, scan.next++);
            const gbwt::size_type seq = sequence_at(here);
            if (seq < haplotype_of_sequence->size()) {
                const size_t found = (*haplotype_of_sequence)[seq];
                if (found < POSITION_KEY_HAPLOTYPES) {
                    position_cache.emplace(position_key(node, found), here);
                }
                if (found == hap) {
                    return here;
                }
            }
        }
        position_cache[key] = gbwt::invalid_edge();
        return gbwt::invalid_edge();
    }

    /// The GBWT sequence at a position, as `locate` gives it, remembered: positions are compared
    /// by sequence again and again, and most of them were located while a node was scanned.
    gbwt::size_type sequence_at(gbwt::edge_type position) {
        const bool cacheable = position.second < (1u << 24);
        const uint64_t key = ((uint64_t)position.first << 24) | (uint64_t)position.second;
        if (cacheable) {
            auto hit = sequence_cache.find(key);
            if (hit != sequence_cache.end()) {
                return hit->second;
            }
        }
        const gbwt::size_type seq = index->locate(position);
        if (cacheable) {
            sequence_cache.emplace(key, seq);
        }
        return seq;
    }

    /// Follow the haplotype at `start` to `to_node`, and report where it arrives.
    ///
    /// The caller gives the direction. A GBWT stores each path in both orientations, so where a
    /// haplotype visits a node has two answers, while where it gets to along a known walk has one.
    /// Given the oriented node, `position_at` finds the position, and `LF` continues in the same
    /// direction. A local guess at the direction, such as the node's forward orientation, assumes
    /// the walk advances in reference order, which fails where the sample's walk does not follow
    /// the reference, as at large balanced structural variants.
    bool follow(gbwt::edge_type start, int64_t to_node, gbwt::node_type* out_end) const {
        if (index == nullptr || start == gbwt::invalid_edge()) {
            return false;
        }
        // Finding a position is far more costly than an LF step, so the caller passes the
        // position in, and the walk itself is cheap.
        if ((int64_t)gbwt::Node::id(start.first) == to_node) {
            if (out_end != nullptr) *out_end = start.first;
            return true;
        }
        gbwt::edge_type at = start;
        for (size_t step = 0; step < WALK_LIMIT; ++step) {
            at = index->LF(at);
            if (at == gbwt::invalid_edge() || at.first == gbwt::ENDMARKER) {
                return false;
            }
            if ((int64_t)gbwt::Node::id(at.first) == to_node) {
                if (out_end != nullptr) *out_end = at.first;
                return true;
            }
        }
        return false;
    }

    /// Append the oriented nodes the haplotype at `start` visits, from there to where `follow`
    /// stops for the node of `end`, to `out`. Return false, appending nothing, where `follow`
    /// would fail.
    bool nodes_to(gbwt::edge_type start, gbwt::node_type end,
                  vector<gbwt::node_type>& out) const {
        if (index == nullptr || start == gbwt::invalid_edge()) {
            return false;
        }
        const size_t before = out.size();
        out.push_back(start.first);
        gbwt::edge_type at = start;
        for (size_t step = 0; gbwt::Node::id(at.first) != gbwt::Node::id(end); ++step) {
            at = index->LF(at);
            if (step >= WALK_LIMIT || at == gbwt::invalid_edge()
                || at.first == gbwt::ENDMARKER) {
                out.resize(before);
                return false;
            }
            out.push_back(at.first);
        }
        return true;
    }

    /// Where `hap` sits at a node in either orientation, forward first, or invalid.
    ///
    /// Forward first, since snarl boundaries are stored oriented along the reference, so the
    /// reverse orientation is the exception. The answer is ambiguous, since a GBWT stores each
    /// path in both orientations. It serves only callers with no direction to work from, which ask
    /// whether the haplotype is at the node at all; a position on a walk comes from `follow`.
    gbwt::edge_type gbwt_position(int64_t node_id, size_t hap) {
        for (int orientation = 0; orientation < 2; ++orientation) {
            const gbwt::edge_type at =
                position_at(gbwt::Node::encode(node_id, orientation == 1), hap);
            if (at != gbwt::invalid_edge()) {
                return at;
            }
        }
        return gbwt::invalid_edge();
    }

private:
    /// The number of haplotypes `position_key` can tell apart. Positions of later haplotypes are
    /// not cached.
    static constexpr size_t POSITION_KEY_HAPLOTYPES = (size_t)1 << 20;

    /// Return the key of a (node, haplotype) pair in `position_cache`.
    ///
    /// @pre `hap` is less than POSITION_KEY_HAPLOTYPES.
    static uint64_t position_key(gbwt::node_type node, size_t hap) {
        return ((uint64_t)node << 20) | (uint64_t)hap;
    }

    /// Return the first position of `hap` at an oriented node, or invalid, scanning the whole
    /// range of the node.
    ///
    /// A haplotype too large for `position_key` needs this. The incremental scan in
    /// `position_at` would already have passed it without caching it.
    gbwt::edge_type search_uncached(gbwt::node_type node, size_t hap) {
        gbwt::SearchState state = index->find(node);
        if (state.empty()) {
            return gbwt::invalid_edge();
        }
        for (gbwt::size_type offset = state.range.first; offset <= state.range.second; ++offset) {
            const gbwt::edge_type here(node, offset);
            const gbwt::size_type seq = sequence_at(here);
            if (seq < haplotype_of_sequence->size() && (*haplotype_of_sequence)[seq] == hap) {
                return here;
            }
        }
        return gbwt::invalid_edge();
    }

    /// How far a node's range has been scanned: offsets [next, end) are still to be located.
    struct Scan {
        bool started = false;
        gbwt::size_type next = 0;
        gbwt::size_type end = 0;
    };

    const gbwt::GBWT* index;
    const vector<size_t>* haplotype_of_sequence;
    /// (oriented node, haplotype) -> GBWT position.
    unordered_map<uint64_t, gbwt::edge_type> position_cache;
    /// Oriented node -> how far its range has been scanned.
    unordered_map<gbwt::node_type, Scan> scans;
    /// GBWT position -> its sequence.
    unordered_map<uint64_t, gbwt::size_type> sequence_cache;
};

/// What a strand holds at a site. A run is cut where the kind changes, so each row has one kind.
enum class StrandKind {
    /// The strand's allele, on a panel haplotype the linkage model named.
    Carried,
    /// Nothing. The strand is the other strand of a nested ploidy-1 site, whose parent allele
    /// bypasses the site.
    Empty,
    /// Sequence the panel cannot attribute to a haplotype.
    Unexplained
};

/// Return the panel haplotype that strand 0 or 1 copies at a site, or NO_HAPLOTYPE.
size_t strand_haplotype(const PhaseCall& call, int strand) {
    return strand == 0 ? call.hap_first : call.hap_second;
}

/// Return what strand 0 or 1 holds at a site.
StrandKind strand_kind(const PhaseCall& call, int strand) {
    if (strand_haplotype(call, strand) != LinkageModel::NO_HAPLOTYPE) {
        return StrandKind::Carried;
    }
    if (call.nested_strand >= 0 && (int)call.nested_strand != strand) {
        return StrandKind::Empty;
    }
    return StrandKind::Unexplained;
}

/// Return which of a site's two called alleles strand 0 or 1 carries, 0 for the allele of
/// `trav_first` and 1 for that of `trav_second`, or -1 if the site is not on the strand.
///
/// A ploidy-1 site has one allele, on the strand `nested_strand` names, or on strand 0 where it
/// names none, which is the strand the `H` rows put such a site on.
int called_slot(const PhaseCall& call, int strand) {
    if (call.ploidy == 1) {
        return strand == (call.nested_strand >= 0 ? (int)call.nested_strand : 0) ? 0 : -1;
    }
    return strand;
}

/// What an exact mosaic's `A` rows are written from, for one contig.
struct ExactInputs {
    /// The phase call of every genotyped site of the contig, with or without a VCF line, in
    /// reference order.
    std::span<const PhaseCall> genotyped;
    /// The called alleles of those sites.
    const CalledAlleles* alleles = nullptr;
    /// The graph, for the sequences of nodes.
    const HandleGraph* graph = nullptr;
};

/**
 * One strand's walk through the graph, as its `H` rows spell it and its `A` rows change it.
 *
 * This is the walk a consumer of an exact mosaic holds: the `H` rows of each fragment joined on
 * their shared nodes, then changed by each `A` row in turn. The writer keeps it so that it can
 * tell where the walk does not take a called allele, and what an `A` row there has to say.
 *
 * The walk is a set of slots, one per visit of a node, each linked to the next, so that a
 * stretch can be replaced without moving the rest. A replaced slot is dead.
 */
class StrandWalk {
public:
    static constexpr uint32_t NO_SLOT = std::numeric_limits<uint32_t>::max();

    /// One crossing of a site: the slots of its two bound nodes, in walk order.
    using Crossing = pair<uint32_t, uint32_t>;

    /// Forget the walk.
    void clear() {
        node.clear();
        next.clear();
        owner.clear();
        live.clear();
        row_index.clear();
        added_index.clear();
        tail = NO_SLOT;
    }

    /// Add the walk of an `H` row of fragment `fragment`, and return the slot of its last node.
    ///
    /// The row joins the row before it where the two are in one fragment and share the node
    /// they meet at, which is then one slot.
    uint32_t append_row(size_t fragment, const vector<gbwt::node_type>& row_nodes) {
        const bool joins = tail != NO_SLOT && tail_fragment == fragment
                           && node[tail] == row_nodes.front();
        if (!joins) {
            tail = NO_SLOT;
        }
        for (size_t i = joins ? 1 : 0; i < row_nodes.size(); ++i) {
            const uint32_t slot = add_slot(row_nodes[i], -1);
            if (tail != NO_SLOT) {
                next[tail] = slot;
            }
            tail = slot;
        }
        tail_fragment = fragment;
        return tail;
    }

    /// Index the slots of the `H` rows by node. Call it once, after the last `append_row`.
    void index_rows() {
        row_index.reserve(node.size());
        for (uint32_t slot = 0; slot < node.size(); ++slot) {
            row_index.emplace_back(node[slot], slot);
        }
        std::ranges::sort(row_index);
    }

    /// Return the `A` row whose walk slot `slot` is from, or -1 if it is from an `H` row.
    int32_t owner_of(uint32_t slot) const { return owner[slot]; }

    /// List the stretches of the walk that run from oriented node `from` to the next visit of
    /// oriented node `to`, without coming back to `from`'s node on the way.
    vector<Crossing> crossings(gbwt::node_type from, gbwt::node_type to) const {
        vector<Crossing> found;
        auto try_from = [&](uint32_t slot) {
            if (!live[slot]) {
                return;
            }
            uint32_t at = slot;
            for (size_t step = 0; node[at] != to; ++step) {
                at = next[at];
                if (at == NO_SLOT || step >= WALK_LIMIT
                    || gbwt::Node::id(node[at]) == gbwt::Node::id(from)) {
                    return;
                }
            }
            found.emplace_back(slot, at);
        };
        const auto rows = std::ranges::equal_range(row_index, from, {},
                                                   &pair<gbwt::node_type, uint32_t>::first);
        for (const auto& [visited, slot] : rows) {
            try_from(slot);
        }
        const auto added = added_index.equal_range(from);
        const size_t from_rows = found.size();
        for (auto it = added.first; it != added.second; ++it) {
            try_from(it->second);
        }
        // An unordered map lists equal keys in no fixed order.
        std::sort(found.begin() + from_rows, found.end());
        return found;
    }

    /// Return whether the walk takes exactly `walk` over a crossing.
    bool spells(const Crossing& crossing, std::span<const gbwt::node_type> walk) const {
        uint32_t at = crossing.first;
        for (size_t i = 0; i < walk.size(); ++i) {
            if (node[at] != walk[i]) {
                return false;
            }
            if (at == crossing.second) {
                return i + 1 == walk.size();
            }
            at = next[at];
        }
        return false;
    }

    /// Make the walk take `walk` over a crossing, whose bound nodes `walk` starts and ends
    /// with. The new slots are from `A` row `row`.
    void replace(const Crossing& crossing, std::span<const gbwt::node_type> walk, int32_t row) {
        if (crossing.first == crossing.second) {
            return;
        }
        for (uint32_t at = next[crossing.first]; at != crossing.second; at = next[at]) {
            live[at] = 0;
        }
        uint32_t before = crossing.first;
        for (size_t i = 1; i + 1 < walk.size(); ++i) {
            const uint32_t slot = add_slot(walk[i], row);
            added_index.emplace(walk[i], slot);
            next[before] = slot;
            before = slot;
        }
        next[before] = crossing.second;
    }

    /// Add `walk` as a piece of the walk that joins nothing, from `A` row `row`.
    void add_piece(std::span<const gbwt::node_type> walk, int32_t row) {
        uint32_t before = NO_SLOT;
        for (const gbwt::node_type& visited : walk) {
            const uint32_t slot = add_slot(visited, row);
            added_index.emplace(visited, slot);
            if (before != NO_SLOT) {
                next[before] = slot;
            }
            before = slot;
        }
    }

private:
    /// Add a live slot for a visit of `visited`, linked to nothing, and return it.
    uint32_t add_slot(gbwt::node_type visited, int32_t row) {
        node.push_back(visited);
        next.push_back(NO_SLOT);
        owner.push_back(row);
        live.push_back(1);
        return (uint32_t)(node.size() - 1);
    }

    /// The oriented node of each slot.
    vector<gbwt::node_type> node;
    /// The slot after each slot, or NO_SLOT at the end of a piece of the walk.
    vector<uint32_t> next;
    /// The `A` row each slot is from, or -1 for a slot from an `H` row.
    vector<int32_t> owner;
    /// Whether each slot is still on the walk.
    vector<uint8_t> live;
    /// The slots from `H` rows, as (node, slot), sorted. They are numbered in walk order.
    vector<pair<gbwt::node_type, uint32_t>> row_index;
    /// The slots from `A` rows, by node.
    unordered_multimap<gbwt::node_type, uint32_t> added_index;
    /// The last slot of the last `H` row added, and that row's fragment.
    uint32_t tail = NO_SLOT;
    size_t tail_fragment = 0;
};

/**
 * A writer of one contig's mosaic rows.
 *
 * It cuts each strand's sites into runs on one panel haplotype, and writes each run as one row
 * per GBWT fragment, since a row's one GBWT position must reach the whole row. It carries the
 * walk from one row of a strand to the next, so that consecutive rows join on the same oriented
 * node and the strand reads as one walk.
 *
 * A contig's rows depend only on its own sites, so each contig can get its own writer on its
 * own thread.
 */
class ContigRowWriter {
public:
    /// Prepare to write the rows of `sites`, one contig's phase calls in reference order, to
    /// `out`. `reference_hap` is the reference's index in the panel, or NO_HAPLOTYPE.
    ///
    /// With `exact`, each strand's `A` rows are written among its `H` rows.
    ContigRowWriter(std::span<const PhaseCall> sites, const PanelLookup& panel,
                    size_t reference_hap, const MosaicParams& params, MosaicCounters& counters,
                    ostream& contig_out, const ExactInputs* exact = nullptr)
        : sites(sites), walk(panel), reference_hap(reference_hap), params(params),
          counters(counters), contig_out(contig_out), exact(exact) {}

    /// Write the rows of every strand of the contig.
    void write();

    /// The number of rows written.
    size_t rows_written() const { return total_rows; }

    /// The number of rows written that name no haplotype.
    size_t unexplained_rows() const { return unexplained; }

private:
    /// How a row ends against the next run, and the reference fill written after it.
    struct RowEnd {
        /// The node the row's walk ends at.
        int64_t to_node = 0;
        /// The node a reference fill after the row ends at, or -1 for no fill.
        int64_t fill_to = -1;
        /// The reference positions written for the fill.
        size_t fill_from_position = 0;
        size_t fill_to_position = 0;
        /// True when nothing crosses from this row to the next, so the next row starts a new
        /// fragment.
        bool gap_open = false;
    };

    /// The walk found for a row.
    struct RowWalk {
        /// The row's GBWT position, or invalid when the row cannot be walked.
        gbwt::edge_type position = gbwt::invalid_edge();
        /// The oriented node where the walk ends.
        gbwt::node_type end = gbwt::ENDMARKER;
        /// True when the row is walked on the reference in place of its own haplotype.
        bool as_ref = false;
        bool walkable = false;
        /// True when the row is walked against the direction the previous row arrived in.
        bool direction_broken = false;
    };

    /// Return the phase call of site `pos` of the current strand's sites.
    const PhaseCall& site(size_t pos) const { return sites[strand_sites[pos]]; }

    /// Return the name of panel haplotype `hap`, or "?" for an index with no name.
    string haplotype_name(size_t hap) const {
        return hap < params.haplotype_names.size() ? params.haplotype_names[hap] : string("?");
    }

    /// Return the number of strands to write: 1 on a haploid contig, where a second would claim
    /// a copy the sample lacks, and 2 otherwise.
    int strand_count() const;

    /// Fill `strand_sites` with the sites the current strand's rows cover.
    void collect_strand_sites();

    /// Write the current strand's rows, one run of one haplotype and kind at a time.
    void write_runs();

    /// Write the run over sites [from, to] of the strand on haplotype `hap`, as one row per GBWT
    /// fragment.
    void write_span(size_t from, size_t to, size_t hap);

    /// Write a run on haplotype `hap` whose first site the haplotype does not reach. The sites
    /// up to the first one it reaches become a row with no position, and the rest is written as
    /// its own run.
    void write_clipped_span(size_t from, size_t to, size_t hap);

    /// Return the last site in [from, to) on the GBWT sequence `sequence` of haplotype `hap`,
    /// given that `from` is on it and `to` is not.
    size_t last_site_on_sequence(size_t from, size_t to, size_t hap, gbwt::size_type sequence);

    /// Write sites [a_idx, b_idx] of the strand as one row on haplotype `hap`, starting from GBWT
    /// position `position`, followed by any reference fill to the next row.
    void write_row(size_t a_idx, size_t b_idx, size_t hap, gbwt::edge_type position);

    /// Decide where the row ending at site `b_idx` ends, given the next site of the strand.
    ///
    /// A row normally ends where the next one begins, so that the stretch between their sites is
    /// covered. Where the row's haplotype cannot be followed there, the next row's start may move
    /// back instead (`pending_from_node`), or a reference fill closes the gap, or the gap is left.
    RowEnd plan_row_end(size_t b_idx, size_t hap, gbwt::edge_type position);

    /// Return whether haplotype `hap`, at GBWT position `position`, reaches the start of site
    /// `next` on the same GBWT sequence, in the direction the next row's haplotype passes it.
    bool extends_right(const PhaseCall& next, size_t hap, gbwt::edge_type position);

    /// Try to start the next row, on haplotype `next_hap`, back at the end node of site `b`,
    /// which this row's haplotype `hap` passes in the same direction. Return true on success.
    bool extend_left(const PhaseCall& b, int64_t next_start, size_t hap, size_t next_hap);

    /// Plan a reference fill from the end of site `b` to the start of site `next`, in `end`, if
    /// the reference crosses that stretch on one GBWT sequence.
    void plan_reference_fill(const PhaseCall& b, const PhaseCall& next, RowEnd& end);

    /// Find the walk of a row from `from_node` to `to_node` on haplotype `hap`.
    ///
    /// The row is walked on its haplotype in the carried direction if it can be. Failing that,
    /// with gap patching on, it is walked on the reference. As a last resort, each is tried
    /// against the carried direction, which happens at an inversion.
    RowWalk find_row_walk(size_t hap, int64_t from_node, int64_t to_node);

    /// Try to walk haplotype `h` from `from_node` to `to_node`, setting `row`'s position and end
    /// on success.
    ///
    /// Where the previous row ended at `from_node`, the walk starts in the direction it arrived
    /// in, unless `ignore_carry` is set. Otherwise both orientations of `from_node` are tried,
    /// forward first.
    bool start_and_walk(size_t h, int64_t from_node, int64_t to_node, bool ignore_carry,
                        RowWalk& row);

    /// Mark `row` as broken if it starts against the carried direction, which it can only do
    /// where `carry_applies`.
    void check_direction(RowWalk& row, bool carry_applies);

    /// Print the row from site `a` to site `b`, starting at `from_node`, which covers
    /// `site_count` sites.
    void print_row(const PhaseCall& a, const PhaseCall& b, int64_t from_node, size_t hap,
                   const RowWalk& row, size_t site_count);

    /// Write the reference fill `end` plans after the row from site `a` to site `b`, or count a
    /// gap if the reference cannot be walked from where the row ended.
    void write_reference_fill(const PhaseCall& a, const PhaseCall& b, const RowEnd& end);

    /// An `H` row of the current strand, as the `A` rows need it.
    struct RowNote {
        /// Where the row's line ends in `out`.
        size_t text_end = 0;
        /// The row's fragment.
        size_t fragment = 0;
        /// The slot of the row's last node in `strand_walk`, if the row could be walked.
        uint32_t last_slot = StrandWalk::NO_SLOT;
    };

    /// An `A` row of the current strand.
    struct AlleleRow {
        /// The index in `row_notes` of the `H` row this row is written after.
        size_t host = 0;
        /// The index of the outermost `A` row whose walk this row applies to, or its own.
        size_t root = 0;
        /// The called allele's walk, read the way the strand crosses the site.
        vector<gbwt::node_type> walk;
        /// Whether the walk does not spell the allele, so that the row gives its sequence.
        bool with_sequence = false;
    };

    /// A change to the sequence of a node, in the node's forward orientation: `removed` bases
    /// from offset `at` on become `added`.
    struct NodeEdit {
        size_t at = 0;
        size_t removed = 0;
        string added;
        /// The record key of the site whose allele makes the change.
        size_t key = 0;
    };

    /// Note the `H` row just printed, which starts at `position` and ends at oriented node
    /// `end`, and add its walk to `strand_walk`.
    void note_row(gbwt::edge_type position, gbwt::node_type end);

    /// Write the current strand's `H` rows from `out`, with its `A` rows among them.
    ///
    /// Sites are taken enclosing site first. Each strand's called allele is compared with what
    /// `strand_walk` takes between the site's bounds, and where they differ an `A` row is made
    /// and `strand_walk` is changed as the row says, so that a site inside is compared with
    /// what its enclosing site's row left.
    void write_allele_rows();

    /// Compare one called allele with `strand_walk` at site `call`, and add any `A` rows.
    void place_allele(const PhaseCall& call, const CalledAlleles::Allele& allele,
                      vector<AlleleRow>& rows);

    /// Add the `A` row that puts `target` on `found`, the crossings of one site in one
    /// direction, if any of them takes another walk or if `with_sequence`.
    void add_allele_row(std::span<const gbwt::node_type> target,
                        const vector<StrandWalk::Crossing>& found, bool with_sequence,
                        vector<AlleleRow>& rows);

    /// Return the index in `row_notes` of the first `H` row whose walk has slot `slot`, which
    /// is a slot from an `H` row.
    size_t row_of_slot(uint32_t slot) const;

    /// Record the changes of a called allele of the site with record key `key` by node, in
    /// `node_edits`.
    void note_edits(const CalledAlleles::Allele& allele, size_t key);

    /// Return the bases of oriented node `visited` on the strand, with the changes of
    /// `node_edits` made. The bases a change puts before the first base the visit reads are
    /// left out unless `with_before`, and those it puts after the last unless `with_after`.
    string strand_sequence(gbwt::node_type visited, bool with_before, bool with_after) const;

    /// Print `A` row `row`, as a row of fragment `row_fragment`.
    void print_allele_row(const AlleleRow& row, size_t row_fragment);

    /// The contig's phase calls, in reference order.
    std::span<const PhaseCall> sites;
    /// Finds and follows haplotypes in the panel's GBWT, with caches for this contig only.
    PanelWalker walk;
    /// The reference's index in the panel, or NO_HAPLOTYPE.
    size_t reference_hap;
    const MosaicParams& params;
    MosaicCounters& counters;
    /// Where the contig's rows go, a strand at a time.
    ostream& contig_out;
    /// What the `A` rows are written from, or null when the mosaic is not exact.
    const ExactInputs* exact;
    /// The current strand's `H` rows. They are held until the strand is done, since its `A`
    /// rows are known only then and go among them.
    ostringstream out;
    /// The current strand's `H` rows, in the order written. Kept only for an exact mosaic, as
    /// are the four members after it.
    vector<RowNote> row_notes;
    /// The rows of `row_notes` that could be walked, as (last slot, index in `row_notes`).
    vector<pair<uint32_t, size_t>> walked_rows;
    /// The current strand's walk.
    StrandWalk strand_walk;
    /// The record key of each site that an `H` row with no position covers, with that row's
    /// index in `row_notes`.
    unordered_map<size_t, size_t> unwalked_sites;
    /// The changes the strand's called alleles make to nodes, by node ID.
    unordered_map<nid_t, vector<NodeEdit>> node_edits;
    size_t total_rows = 0;
    size_t unexplained = 0;

    /// The strand being written, 0 or 1.
    int strand = 0;
    /// The sites the strand's rows cover, as indices into `sites`, in reference order.
    vector<size_t> strand_sites;
    /// The strand's current fragment. (contig, strand, fragment) names one contiguous walk, which
    /// a loader makes into one path.
    size_t fragment = 0;
    /// The node the next row starts at, moved back by a left extension, or -1 if not moved.
    int64_t pending_from_node = -1;
    /// The GBWT position of the next row's moved start.
    gbwt::edge_type pending_from_pos = gbwt::invalid_edge();
    /// The oriented node this strand's walk has reached, or ENDMARKER when the walk is broken. It
    /// carries the walk's direction from one row to the next.
    gbwt::node_type carry = gbwt::ENDMARKER;
};

void ContigRowWriter::write() {
    const int strands = strand_count();
    for (strand = 0; strand < strands; ++strand) {
        pending_from_node = -1;
        carry = gbwt::ENDMARKER;
        fragment = 0;
        out.str(string());
        row_notes.clear();
        walked_rows.clear();
        strand_walk.clear();
        unwalked_sites.clear();
        node_edits.clear();
        collect_strand_sites();
        write_runs();
        if (exact != nullptr) {
            write_allele_rows();
        } else {
            contig_out << out.str();
        }
    }
}

void ContigRowWriter::note_row(gbwt::edge_type position, gbwt::node_type end) {
    RowNote note;
    note.text_end = (size_t)out.tellp();
    note.fragment = fragment;
    if (vector<gbwt::node_type> row_nodes; walk.nodes_to(position, end, row_nodes)) {
        note.last_slot = strand_walk.append_row(fragment, row_nodes);
        walked_rows.emplace_back(note.last_slot, row_notes.size());
    }
    row_notes.push_back(note);
}

void ContigRowWriter::write_allele_rows() {
    strand_walk.index_rows();
    const std::span<const PhaseCall> genotyped = exact->genotyped;
    // An enclosing site has a lower level than the sites inside it.
    vector<size_t> order(genotyped.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::stable_sort(order, {}, [&](size_t g) { return genotyped[g].level; });
    unordered_set<size_t> written;
    written.reserve(sites.size() * 2);
    for (const PhaseCall& call : sites) {
        written.insert(call.record_key);
    }

    vector<AlleleRow> rows;
    for (size_t g : order) {
        const PhaseCall& call = genotyped[g];
        const int slot = called_slot(call, strand);
        if (slot < 0) {
            continue;
        }
        CalledAlleles::Allele allele;
        if (!exact->alleles->find(call.record_key, slot, allele) || allele.walk.empty()
            || allele.traversal != (slot == 0 ? call.trav_first : call.trav_second)) {
            ++counters.strands_without_allele;
            continue;
        }
        if (allele.no_change) {
            continue;
        }
        const size_t rows_before = rows.size();
        place_allele(call, allele, rows);
        if (written.count(call.record_key) == 0) {
            counters.allele_rows_unwritten += rows.size() - rows_before;
        }
    }

    // Each `H` row is followed by the `A` rows it hosts: each outermost row, then the rows that
    // apply to its walk, in the order made, which is enclosing site first.
    vector<size_t> row_order(rows.size());
    std::iota(row_order.begin(), row_order.end(), 0);
    std::ranges::sort(row_order, {}, [&](size_t r) {
        return std::make_tuple(rows[r].host, rows[r].root, r);
    });
    const string h_rows = out.str();
    size_t text_from = 0;
    auto next_row = row_order.begin();
    for (size_t h = 0; h < row_notes.size(); ++h) {
        contig_out.write(h_rows.data() + text_from, row_notes[h].text_end - text_from);
        text_from = row_notes[h].text_end;
        for (; next_row != row_order.end() && rows[*next_row].host == h; ++next_row) {
            print_allele_row(rows[*next_row], row_notes[h].fragment);
            counters.allele_rows_nested += (size_t)(rows[*next_row].root != *next_row);
        }
    }
}

void ContigRowWriter::place_allele(const PhaseCall& call, const CalledAlleles::Allele& allele,
                                   vector<AlleleRow>& rows) {
    ++counters.strands_compared;
    // A strand can cross a site from its end to its start, as inside an inverted allele of an
    // enclosing site, and then it takes the allele's walk read backwards.
    vector<gbwt::node_type> backwards(allele.walk.rbegin(), allele.walk.rend());
    for (gbwt::node_type& visited : backwards) {
        visited = gbwt::Node::reverse(visited);
    }
    bool crossed = false;
    for (bool reversed : {false, true}) {
        const std::span<const gbwt::node_type> target =
            reversed ? std::span<const gbwt::node_type>(backwards) : allele.walk;
        if (target.size() > 1 && target.front() == target.back()) {
            // A walk that comes back to the node it starts on has no stretch to stand for.
            continue;
        }
        const vector<StrandWalk::Crossing> found =
            strand_walk.crossings(target.front(), target.back());
        if (found.empty()) {
            continue;
        }
        if (!crossed && allele.edits != nullptr) {
            note_edits(allele, call.record_key);
        }
        crossed = true;
        const size_t rows_before = rows.size();
        add_allele_row(target, found, allele.edits != nullptr, rows);
        counters.allele_rows_reversed += (size_t)(reversed && rows.size() > rows_before);
    }
    if (crossed) {
        return;
    }
    const auto unwalked = unwalked_sites.find(call.record_key);
    if (unwalked == unwalked_sites.end()) {
        ++counters.strands_not_crossed;
        return;
    }
    // The site's `H` row has no walk, so the allele's walk stands alone, and sites inside it
    // are compared with it.
    if (allele.edits != nullptr) {
        note_edits(allele, call.record_key);
    }
    AlleleRow row;
    row.host = unwalked->second;
    row.root = rows.size();
    row.walk.assign(allele.walk.begin(), allele.walk.end());
    row.with_sequence = allele.edits != nullptr;
    strand_walk.add_piece(allele.walk, (int32_t)rows.size());
    rows.push_back(std::move(row));
    ++counters.allele_rows_standalone;
}

void ContigRowWriter::add_allele_row(std::span<const gbwt::node_type> target,
                                     const vector<StrandWalk::Crossing>& found,
                                     bool with_sequence, vector<AlleleRow>& rows) {
    const size_t self = rows.size();
    AlleleRow row;
    row.root = self;
    row.with_sequence = with_sequence;
    bool differs = false;
    bool placed = false;
    for (const StrandWalk::Crossing& crossing : found) {
        differs = differs || !strand_walk.spells(crossing, target);
        // A crossing on the walk of an earlier `A` row needs that row first, so this row goes
        // with it. A crossing on the `H` rows needs the row its end bound is in. The row goes
        // where the last of its crossings needs it.
        const int32_t owner = std::max(strand_walk.owner_of(crossing.first),
                                       strand_walk.owner_of(crossing.second));
        const size_t host = owner >= 0 ? rows[owner].host : row_of_slot(crossing.second);
        const size_t root = owner >= 0 ? rows[owner].root : self;
        if (!placed || std::tie(host, root) > std::tie(row.host, row.root)) {
            row.host = host;
            row.root = root;
            placed = true;
        }
    }
    if (!differs && !with_sequence) {
        return;
    }
    for (const StrandWalk::Crossing& crossing : found) {
        if (!strand_walk.spells(crossing, target)) {
            strand_walk.replace(crossing, target, (int32_t)self);
        }
    }
    row.walk.assign(target.begin(), target.end());
    rows.push_back(std::move(row));
}

size_t ContigRowWriter::row_of_slot(uint32_t slot) const {
    // Slots are numbered along the rows, so the walked rows' last slots are in order.
    const auto found = std::ranges::lower_bound(walked_rows, slot, {},
                                                &pair<uint32_t, size_t>::first);
    return found != walked_rows.end() ? found->second : walked_rows.back().second;
}

void ContigRowWriter::note_edits(const CalledAlleles::Allele& allele, size_t key) {
    for (const AlleleEdit& edit : *allele.edits) {
        if (edit.step >= allele.walk.size()) {
            continue;
        }
        const gbwt::node_type visited = allele.walk[edit.step];
        const nid_t id = gbwt::Node::id(visited);
        const size_t length = exact->graph->get_length(exact->graph->get_handle(id));
        if (edit.at + edit.removed > length) {
            continue;
        }
        // The allele's walk may read the node backwards, and so may another allele that changes
        // the same node, so the changes are kept in the node's forward orientation.
        NodeEdit change;
        change.key = key;
        change.removed = edit.removed;
        if (gbwt::Node::is_reverse(visited)) {
            change.at = length - edit.at - edit.removed;
            change.added = reverse_complement(edit.added);
        } else {
            change.at = edit.at;
            change.added = edit.added;
        }
        node_edits[id].push_back(std::move(change));
    }
}

string ContigRowWriter::strand_sequence(gbwt::node_type visited, bool with_before,
                                        bool with_after) const {
    const nid_t id = gbwt::Node::id(visited);
    const bool reversed = gbwt::Node::is_reverse(visited);
    string own = exact->graph->get_sequence(exact->graph->get_handle(id));
    if (reversed) {
        own = reverse_complement(own);
    }
    const auto found = node_edits.find(id);
    if (found == node_edits.end()) {
        return own;
    }
    // The changes as the visit reads the node, made as a record makes them where it spells a
    // strand through the node (`AlleleSplicer::apply_step_edits`), so that two changes that
    // take one base leave the mosaic and the records with the same one of them.
    vector<AlleleSplicer::StepEdit> edits;
    for (const NodeEdit& change : found->second) {
        const size_t at = reversed ? own.size() - change.at - change.removed : change.at;
        if (change.removed == 0
            && ((at == 0 && !with_before) || (at == own.size() && !with_after))) {
            continue;
        }
        edits.push_back(AlleleSplicer::StepEdit{
            0, at, change.removed, reversed ? reverse_complement(change.added) : change.added,
            change.key});
    }
    StepSequences steps{std::move(own)};
    AlleleSplicer::apply_step_edits(std::move(edits), steps);
    return steps[0];
}

void ContigRowWriter::print_allele_row(const AlleleRow& row, size_t row_fragment) {
    contig_out << "A\t" << sites.front().contig << "\t" << strand << "\t" << row_fragment << "\t"
               << row.walk.front() << "\t" << row.walk.back() << "\t";
    for (const gbwt::node_type& visited : row.walk) {
        contig_out << (gbwt::Node::is_reverse(visited) ? '<' : '>') << gbwt::Node::id(visited);
    }
    contig_out << "\t";
    if (row.with_sequence) {
        // Bases put before the first base of the walk, or after its last, are not over it.
        for (size_t i = 0; i < row.walk.size(); ++i) {
            contig_out << strand_sequence(row.walk[i], i > 0, i + 1 < row.walk.size());
        }
        ++counters.allele_rows_with_sequence;
    } else {
        contig_out << ".";
    }
    contig_out << "\n";
    ++counters.allele_rows;
}

int ContigRowWriter::strand_count() const {
    // Taken over the whole contig rather than from its first site, since a diploid contig can
    // begin with a nested ploidy-1 site.
    return std::ranges::any_of(sites, [](const PhaseCall& call) { return call.ploidy != 1; })
               ? 2
               : 1;
}

void ContigRowWriter::collect_strand_sites() {
    strand_sites.clear();
    for (size_t t = 0; t < sites.size(); ++t) {
        const PhaseCall& call = sites[t];
        // A nested ploidy-1 chain is on one of its parent's strands. The other strand takes the
        // parent's other allele, which bypasses the chain, so the site is not on its walk.
        if (strand_kind(call, strand) == StrandKind::Empty) {
            continue;
        }
        // A switch of haplotype inside a nested chain is real, but it follows the parent's
        // route and a consumer counting recombinations may not want it. Dropping the nested
        // sites merges the runs across them, so the walk follows the parent's haplotype
        // through the child snarl. The header records which was done.
        if (!params.keep_nested && call.level > 0) {
            continue;
        }
        // Dropping a site with no haplotype carries the flanking haplotype through it. That
        // keeps the strand one path, but writes the haplotype's sequence across the site rather
        // than the called allele.
        if (params.connect_unexplained && strand_kind(call, strand) == StrandKind::Unexplained) {
            continue;
        }
        strand_sites.push_back(t);
    }
}

void ContigRowWriter::write_runs() {
    size_t run_start = 0;
    for (size_t t = 0; t < strand_sites.size(); ++t) {
        const size_t hap = strand_haplotype(site(t), strand);
        const StrandKind kind = strand_kind(site(t), strand);
        const bool last = t + 1 == strand_sites.size();
        const bool changes = !last
                             && (strand_haplotype(site(t + 1), strand) != hap
                                 || strand_kind(site(t + 1), strand) != kind);
        if (last || changes) {
            write_span(run_start, t, hap);
            run_start = t + 1;
        }
    }
}

void ContigRowWriter::write_span(size_t from, size_t to, size_t hap) {
    const gbwt::edge_type pos = (hap == LinkageModel::NO_HAPLOTYPE)
                                    ? gbwt::invalid_edge()
                                    : walk.gbwt_position(site(from).start_node, hap);
    if (pos == gbwt::invalid_edge() && hap != LinkageModel::NO_HAPLOTYPE && from != to) {
        write_clipped_span(from, to, hap);
        return;
    }
    if (pos == gbwt::invalid_edge() || from == to) {
        if (pos == gbwt::invalid_edge() && hap != LinkageModel::NO_HAPLOTYPE) {
            ++counters.unwalkable;
        }
        write_row(from, to, hap, pos);
        return;
    }
    // We locate the run's two ends, and search for a fragment boundary only when they are on
    // different GBWT sequences. This misses a haplotype that leaves a fragment and comes back to
    // it within one run, which fragments of one path cannot do in reference order.
    if (const gbwt::edge_type end_pos = walk.gbwt_position(site(to).start_node, hap);
        end_pos == gbwt::invalid_edge() || walk.sequence_at(pos) == walk.sequence_at(end_pos)) {
        // Same fragment at both ends, or no way to tell. One row.
        write_row(from, to, hap, pos);
        return;
    }
    const size_t last = last_site_on_sequence(from, to, hap, walk.sequence_at(pos));
    write_row(from, last, hap, pos);
    write_span(last + 1, to, hap);
}

void ContigRowWriter::write_clipped_span(size_t from, size_t to, size_t hap) {
    // We write the head the haplotype does not reach as small as it really is, rather than
    // giving up on the whole run.
    size_t first_ok = from;
    while (first_ok <= to
           && walk.gbwt_position(site(first_ok).start_node, hap) == gbwt::invalid_edge()) {
        ++first_ok;
    }
    if (first_ok > to) {
        // Clipped across the whole run.
        write_row(from, to, hap, gbwt::invalid_edge());
        ++counters.unwalkable;
        return;
    }
    if (first_ok > from) {
        write_row(from, first_ok - 1, hap, gbwt::invalid_edge());
        ++counters.unwalkable;
        ++counters.head_clipped;
    }
    write_span(first_ok, to, hap);
}

size_t ContigRowWriter::last_site_on_sequence(size_t from, size_t to, size_t hap,
                                              gbwt::size_type sequence) {
    // A site the haplotype does not reach counts as past the boundary, which keeps the search
    // monotone.
    size_t lo = from;
    size_t hi = to;
    while (hi - lo > 1) {
        const size_t mid = std::midpoint(lo, hi);
        const gbwt::edge_type p = walk.gbwt_position(site(mid).start_node, hap);
        if (p != gbwt::invalid_edge() && walk.sequence_at(p) == sequence) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

void ContigRowWriter::write_row(size_t a_idx, size_t b_idx, size_t hap,
                                gbwt::edge_type position) {
    const PhaseCall& a = site(a_idx);
    const PhaseCall& b = site(b_idx);
    // A left extension planned by the previous row moves this row's start back, and its
    // position with it.
    int64_t from_node = a.start_node;
    if (pending_from_node >= 0) {
        from_node = pending_from_node;
        position = pending_from_pos;
    }
    pending_from_node = -1;

    const RowEnd end = plan_row_end(b_idx, hap, position);
    const RowWalk row = find_row_walk(hap, from_node, end.to_node);
    carry = row.walkable ? row.end : gbwt::ENDMARKER;
    // A row whose direction was broken cannot join the row before it, so it starts a new
    // fragment. The row after it can join it as usual.
    if (hap == LinkageModel::NO_HAPLOTYPE || !row.walkable || row.direction_broken) {
        ++fragment;
    }
    print_row(a, b, from_node, hap, row, b_idx - a_idx + 1);
    if (exact != nullptr && row.position == gbwt::invalid_edge()) {
        for (size_t t = a_idx; t <= b_idx; ++t) {
            unwalked_sites.emplace(site(t).record_key, row_notes.size() - 1);
        }
    }
    if (end.fill_to >= 0) {
        write_reference_fill(a, b, end);
    }
    if (end.gap_open) {
        ++fragment;
    }
    // No consumer could walk across a row with no position, so it also ends the fragment.
    if (hap == LinkageModel::NO_HAPLOTYPE || !row.walkable) {
        ++fragment;
    }
}

ContigRowWriter::RowEnd ContigRowWriter::plan_row_end(size_t b_idx, size_t hap,
                                                      gbwt::edge_type position) {
    const PhaseCall& b = site(b_idx);
    RowEnd end;
    end.to_node = b.end_node;
    if (b_idx + 1 >= strand_sites.size()) {
        return end;
    }
    const PhaseCall& next = site(b_idx + 1);
    const int64_t next_start = next.start_node;
    const size_t next_hap = strand_haplotype(next, strand);
    // A change of haplotype between a parent snarl Ps..Pe and a child snarl Cs..Ce has two
    // boundaries, both the child's, since the parent's walk is Ps -> ... -> Cs -> [child] -> Ce
    // -> ... -> Pe. Every traversal of the child visits Cs and Ce, so the two haplotypes meet
    // there.
    //
    // Sites arrive in reference order, and a parent is always recorded, so a site deeper than
    // the one before it is inside that one. Node IDs would not tell, since they need not follow
    // the walk, as in an inversion.
    const bool entering = next.level > b.level;
    const bool leaving = next.level < b.level;
    const bool right = position != gbwt::invalid_edge() && extends_right(next, hap, position);
    if (entering) {
        // The walk reaches the parent's end only after the child, so the row ends where the
        // child begins.
        end.to_node = next_start;
        ++counters.nested_enter;
    } else if (leaving) {
        // The next row starts at the child's end, so the parent's haplotype, whose called
        // allele governs the stretch from there to the parent's end, covers it.
        end.to_node = b.end_node;
        pending_from_node = b.end_node;
        pending_from_pos = next_hap == LinkageModel::NO_HAPLOTYPE
                               ? gbwt::invalid_edge()
                               : walk.gbwt_position(b.end_node, next_hap);
        ++counters.nested_leave;
    } else if (right) {
        // Nothing shows which of the two haplotypes covers the stretch between the runs, so
        // extending this one is a convention, which the header states.
        end.to_node = next_start;
        ++counters.extended;
    } else if (!extend_left(b, next_start, hap, next_hap)) {
        if (params.patch_gaps && reference_hap != LinkageModel::NO_HAPLOTYPE) {
            plan_reference_fill(b, next, end);
        }
        if (end.fill_to < 0) {
            ++counters.gap_left;
            end.gap_open = true;
        }
    }
    return end;
}

bool ContigRowWriter::extends_right(const PhaseCall& next, size_t hap,
                                    gbwt::edge_type position) {
    const gbwt::edge_type np = walk.gbwt_position(next.start_node, hap);
    // The next row's haplotype must pass the junction in the same direction as this one, not
    // only through the same node.
    const size_t next_hap = strand_haplotype(next, strand);
    const gbwt::edge_type entry = next_hap == LinkageModel::NO_HAPLOTYPE
                                      ? gbwt::invalid_edge()
                                      : walk.gbwt_position(next.start_node, next_hap);
    return np != gbwt::invalid_edge()
           && walk.sequence_at(np) == walk.sequence_at(position)
           && (entry == gbwt::invalid_edge() || entry.first == np.first);
}

bool ContigRowWriter::extend_left(const PhaseCall& b, int64_t next_start, size_t hap,
                                  size_t next_hap) {
    if (next_hap == LinkageModel::NO_HAPLOTYPE) {
        return false;
    }
    const gbwt::edge_type here = walk.gbwt_position(b.end_node, next_hap);
    const gbwt::edge_type there = walk.gbwt_position(next_start, next_hap);
    const gbwt::edge_type mine = walk.gbwt_position(b.end_node, hap);
    if (here != gbwt::invalid_edge() && there != gbwt::invalid_edge()
        && walk.sequence_at(here) == walk.sequence_at(there)
        && (mine == gbwt::invalid_edge() || mine.first == here.first)) {
        pending_from_node = b.end_node;
        pending_from_pos = here;
        ++counters.extended_left;
        return true;
    }
    return false;
}

void ContigRowWriter::plan_reference_fill(const PhaseCall& b, const PhaseCall& next,
                                          RowEnd& end) {
    const gbwt::edge_type left = walk.gbwt_position(b.end_node, reference_hap);
    const gbwt::edge_type right = walk.gbwt_position(next.start_node, reference_hap);
    if (left != gbwt::invalid_edge() && right != gbwt::invalid_edge()
        && walk.sequence_at(left) == walk.sequence_at(right)) {
        end.fill_to = next.start_node;
        end.fill_from_position = b.position;
        end.fill_to_position = next.position;
        ++counters.patched;
    }
}

ContigRowWriter::RowWalk ContigRowWriter::find_row_walk(size_t hap, int64_t from_node,
                                                        int64_t to_node) {
    RowWalk row;
    row.end = gbwt::Node::encode(to_node, false);
    // The carried direction constrains the row only where the previous row ended at its start,
    // so not on a strand's first row, nor where an extension moved the start.
    const bool carry_applies = carry != gbwt::ENDMARKER
                               && (int64_t)gbwt::Node::id(carry) == from_node;
    const bool can_use_reference = params.patch_gaps
                                   && reference_hap != LinkageModel::NO_HAPLOTYPE;
    if (hap == LinkageModel::NO_HAPLOTYPE) {
        return row;
    }
    // A row whose own haplotype does not cross it becomes a reference substitution, marked `ref`
    // like a fill but keeping its site count, since it covers called sites.
    row.walkable = start_and_walk(hap, from_node, to_node, false, row);
    if (!row.walkable && can_use_reference
        && start_and_walk(reference_hap, from_node, to_node, false, row)) {
        row.as_ref = true;
        row.walkable = true;
        ++counters.row_to_ref;
    }
    // A row about to have no position has no contiguity left to break, so a walk against the
    // carried direction is better than none.
    if (!row.walkable && start_and_walk(hap, from_node, to_node, true, row)) {
        row.walkable = true;
        check_direction(row, carry_applies);
    }
    // A panel haplotype can be clipped across a site, which the linkage model allows, so the
    // reference gets the same last resort.
    if (!row.walkable && can_use_reference
        && start_and_walk(reference_hap, from_node, to_node, true, row)) {
        row.as_ref = true;
        row.walkable = true;
        ++counters.row_to_ref;
        check_direction(row, carry_applies);
    }
    if (!row.walkable) {
        row.position = gbwt::invalid_edge();
        row.end = gbwt::Node::encode(to_node, false);
    }
    return row;
}

bool ContigRowWriter::start_and_walk(size_t h, int64_t from_node, int64_t to_node,
                                     bool ignore_carry, RowWalk& row) {
    if (!ignore_carry && carry != gbwt::ENDMARKER
        && (int64_t)gbwt::Node::id(carry) == from_node) {
        // The carried direction decides. If the haplotype cannot be followed from it, the row
        // does not continue the walk, rather than taking the other orientation and breaking
        // contiguity.
        const gbwt::edge_type at = walk.position_at(carry, h);
        if (at == gbwt::invalid_edge() || !walk.follow(at, to_node, &row.end)) {
            return false;
        }
        row.position = at;
        return true;
    }
    for (int o = 0; o < 2; ++o) {
        const gbwt::edge_type at = walk.position_at(gbwt::Node::encode(from_node, o == 1), h);
        if (at != gbwt::invalid_edge() && walk.follow(at, to_node, &row.end)) {
            row.position = at;
            return true;
        }
    }
    return false;
}

void ContigRowWriter::check_direction(RowWalk& row, bool carry_applies) {
    row.direction_broken = carry_applies && row.position.first != carry;
    if (row.direction_broken) {
        ++counters.direction_broken;
    }
}

void ContigRowWriter::print_row(const PhaseCall& a, const PhaseCall& b, int64_t from_node,
                                size_t hap, const RowWalk& row, size_t site_count) {
    // Nodes are oriented, `id * 2 + is_reverse`, since two rows can share a node and pass it in
    // opposite directions. A row with no position takes the reference orientation.
    const gbwt::node_type row_start = row.position != gbwt::invalid_edge()
                                          ? row.position.first
                                          : gbwt::Node::encode(from_node, false);
    out << "H\t" << a.contig << "\t" << strand << "\t" << fragment << "\t"
        << a.position << "\t" << b.position << "\t"
        << row_start << "\t" << row.end << "\t";
    if (row.as_ref) {
        out << "ref\t" << haplotype_name(reference_hap);
    } else if (hap == LinkageModel::NO_HAPLOTYPE) {
        out << "*\t*";
        ++unexplained;
    } else {
        out << hap << "\t" << haplotype_name(hap);
    }
    out << "\t" << site_count << "\t";
    if (row.position == gbwt::invalid_edge()) {
        // "." rather than 0, which would be a valid offset.
        out << ".";
    } else {
        out << row.position.second;
    }
    out << "\n";
    ++total_rows;
    if (exact != nullptr) {
        note_row(row.position, row.end);
    }
}

void ContigRowWriter::write_reference_fill(const PhaseCall& a, const PhaseCall& b,
                                           const RowEnd& end) {
    // The fill is a walk too, starting where the row ended, so the carried direction runs
    // through it.
    gbwt::node_type fill_end = gbwt::Node::encode(end.fill_to, false);
    const gbwt::node_type start = carry != gbwt::ENDMARKER
                                      ? carry
                                      : gbwt::Node::encode(b.end_node, false);
    const gbwt::edge_type fill_pos = walk.position_at(start, reference_hap);
    if (fill_pos != gbwt::invalid_edge() && walk.follow(fill_pos, end.fill_to, &fill_end)) {
        // The fill covers no called site, so its site count is ".".
        out << "H\t" << a.contig << "\t" << strand << "\t" << fragment << "\t"
            << end.fill_from_position << "\t" << end.fill_to_position << "\t"
            << fill_pos.first << "\t" << fill_end << "\t"
            << "ref\t" << haplotype_name(reference_hap)
            << "\t.\t" << fill_pos.second << "\n";
        ++total_rows;
        if (exact != nullptr) {
            note_row(fill_pos, fill_end);
        }
        carry = fill_end;
    } else {
        ++counters.gap_left;
        --counters.patched;
        carry = gbwt::ENDMARKER;
        ++fragment;
    }
}

/// Write the mosaic's header lines, for `sample_name`, to `out`.
///
/// The header states two things the rows do not: which reference the positions are in, since a
/// graph can hold several references with the same contig names, and what each hap_index means,
/// since the index is internal to the run.
void write_header(ostream& out, const MosaicParams& params, const string& sample_name) {
    out << "#mosaic-version\t" << (params.exact ? 6 : 5) << "\n";
    out << "#graph\t" << params.graph_name << "\n";
    out << "#sample\t" << sample_name << "\n";
    // gRef fragments are counted, not listed: a cover can name thousands of contigs, and each row
    // names its own.
    size_t gref_fragments = 0;
    for (const string& ref : params.reference_paths) {
        if (GrefCover::is_gref_name(ref)) {
            ++gref_fragments;
            continue;
        }
        out << "#reference\t" << ref << "\n";
    }
    if (gref_fragments > 0) {
        out << "#gref-fragments\t" << gref_fragments << "\n";
    }
    out << "#decoding\tconstrained-viterbi\n";
    out << "#patch\t" << (params.patch_gaps ? "reference" : "none") << "\n";
    out << "#nested\t" << (params.keep_nested ? "kept" : "merged") << "\n";
    out << "#unexplained\t" << (params.connect_unexplained ? "connected" : "broken") << "\n";
    out << "#note\tref_start/ref_end are advisory, in the #reference coordinate system; "
        << "start_node/end_node are the authoritative anchors and are intrinsic to the graph.\n";
    out << "#note\tsegments are maximal runs on one panel haplotype; walk the haplotype from "
        << "start_node to end_node to reconstruct it. * means the strand traverses this stretch "
        << "and the panel cannot name a haplotype for it. A site a strand does not traverse is "
        << "simply absent from that strand's rows.\n";
    out << "#note\thap_index ref marks a stretch FILLED WITH THE REFERENCE because no panel "
        << "haplotype could be carried across it. It is walkable like any other row. With sites=. "
        << "it filled a boundary between two segments and covers no called site; with a site count "
        << "it replaced a haplotype the graph does not carry across that segment. Either way the "
        << "reference is a poor proxy for the sample, so the fill is marked rather than blended in; "
        << "its span is start_node..end_node. --no-mosaic-patch-gaps leaves the gap instead.\n";
    out << "#note\tend_node is the NEXT segment's start_node, so consecutive segments of one "
        << "strand meet at a shared node and the strand is ONE WALK: concatenate them, counting "
        << "each junction once. (contig, strand, fragment) is the identity of that walk.\n";
    out << "#note\tWhich haplotype covers the stretch between two segments' sites is ARBITRARY: no "
        << "called site lies in it, so nothing distinguishes the earlier haplotype from the later "
        << "one and a recombination anywhere inside is equally consistent. Extending the earlier "
        << "one is a convention. The crossover is BRACKETED by that stretch, not located within "
        << "it, and a consumer reading the boundary as the crossover point will over-trust it.\n";
    out << "#note\thap_index is internal to this run; haplotype (sample#phase) is the portable "
        << "identifier and names a haplotype, not a single GBWT path.\n";
    for (size_t h = 0; h < params.haplotype_names.size(); ++h) {
        out << "#haplotype\t" << h << "\t" << params.haplotype_names[h] << "\n";
    }
    out << "#note\tstart_node and end_node are ORIENTED node ids, id * 2 + is_reverse. Node "
        << "identity alone does not make a walk: two segments can share a node and traverse it in "
        << "opposite directions. (start_node, gbwt_offset) is therefore the GBWT position "
        << "outright -- extract() it and follow LF() to end_node, with no locate and no "
        << "r-index.\n";
    out << "#note\ta segment never spans a GBWT fragment boundary, so one position walks the "
        << "whole of it; a haplotype in several fragments yields several segments.\n";
    out << "#note\ta fragment is a PATH: its rows join end to start, on the same oriented node, "
        << "and expand to an exact walk in the graph. A row that cannot be walked in the direction "
        << "the fragment has reached starts a new fragment instead -- an inversion boundary is the "
        << "usual reason, and X+ followed by X- is not a walk. A row with no position at all "
        << "(gbwt_offset .) is alone in its fragment, so it never sits inside a walk.\n";
    if (params.exact) {
        out << "#note\tan A row gives the allele a strand was called with at one site, where the "
            << "walk the rows before it spell does not take that allele. start_node and end_node "
            << "are the site's two bound nodes, oriented as in H rows and in the order the strand "
            << "passes them. walk is the allele's node walk from start_node to end_node, both "
            << "included, written >id for a node read forwards and <id for one read backwards.\n";
        out << "#note\tto apply an A row, find every stretch of the strand's walk that runs from "
            << "start_node to its next visit of end_node, and make it walk. If the strand's walk "
            << "has no such stretch, walk is a piece of the strand that joins no row: the site's "
            << "H row has no position. Apply the A rows in file order. Each follows the H rows "
            << "and the A rows it needs: a site's row comes before the rows of the sites inside "
            << "it. The strand then takes the called allele at every genotyped site its walk "
            << "crosses, including sites with no VCF line.\n";
        out << "#note\tsequence is . where walk spells the allele. Otherwise the allele changes "
            << "bases inside nodes, as an allele no panel haplotype carries does, and sequence is "
            << "the strand's bases over walk, from the first base of start_node to the last base "
            << "of end_node, read along walk. It has every such change the strand's alleles make "
            << "to those nodes. Bases that an allele puts before the first base of start_node or "
            << "after the last base of end_node are not over walk, and are in the row of the walk "
            << "they lie inside.\n";
        out << "#note\tan A row's fragment is that of the H row it follows.\n";
    }
    out << "#H\tcontig\tstrand\tfragment\tref_start\tref_end\tstart_node\tend_node"
        << "\thap_index\thaplotype\tsites\tgbwt_offset\n";
    if (params.exact) {
        out << "#A\tcontig\tstrand\tfragment\tstart_node\tend_node\twalk\tsequence\n";
    }
}

/// Return the reference's index in the panel, for filling a gap no haplotype can cross, or
/// NO_HAPLOTYPE if the reference is not a panel haplotype.
///
/// The index follows GBWT metadata order, so it is found by name. The panel names haplotypes as
/// sample#phase (CHM13#0), so the contig is dropped from each reference path name
/// (CHM13#0#chr20) before matching.
size_t reference_haplotype(const MosaicParams& params) {
    for (const string& full : params.reference_paths) {
        // A gRef cover is stitched together from many donors, and it is not in the panel.
        if (GrefCover::is_gref_derived(full)) {
            continue;
        }
        const size_t h1 = full.find('#');
        const size_t h2 = h1 == string::npos ? string::npos : full.find('#', h1 + 1);
        const string base = h2 == string::npos ? full : full.substr(0, h2);
        for (size_t k = 0; k < params.haplotype_names.size(); ++k) {
            if (params.haplotype_names[k] == base) {
                return k;
            }
        }
    }
    return LinkageModel::NO_HAPLOTYPE;
}

/// Return the range of `phasing` that each contig's sites occupy, as [first, past-the-end)
/// pairs. The phasing is grouped by contig.
vector<pair<size_t, size_t>> contig_ranges(const vector<PhaseCall>& phasing) {
    vector<pair<size_t, size_t>> contigs;
    for (size_t i = 0; i < phasing.size();) {
        size_t j = i;
        while (j < phasing.size() && phasing[j].contig == phasing[i].contig) {
            ++j;
        }
        contigs.emplace_back(i, j);
        i = j;
    }
    return contigs;
}

/// Print the counters of the boundaries and rows the mosaic writer handled to stderr.
void report_counters(const MosaicCounters& counters) {
    cerr << "[vg call] mosaic: " << counters.extended.load()
         << " segment boundaries closed by extending right, " << counters.extended_left.load()
         << " by extending left instead, " << counters.patched.load()
         << " filled with the reference because neither haplotype could be carried across, "
         << counters.gap_left.load() << " left as a gap" << endl;
    cerr << "[vg call] mosaic: " << counters.nested_enter.load()
         << " boundaries where the walk enters a child snarl, " << counters.nested_leave.load()
         << " where it leaves one -- both stated at the CHILD's boundary node" << endl;
    cerr << "[vg call] mosaic: " << counters.direction_broken.load()
         << " rows walked against the carried direction, each standing alone (inversions)" << endl;
    cerr << "[vg call] mosaic: " << counters.unwalkable.load()
         << " segments name a haplotype the graph does not carry across them, of which "
         << counters.head_clipped.load() << " are a clipped head whose remainder is walkable; "
         << counters.row_to_ref.load() << " rewritten as a reference substitution" << endl;
}

/// Print the counters of an exact mosaic's `A` rows to stderr.
void report_allele_rows(const MosaicCounters& counters) {
    cerr << "[vg call] mosaic: " << counters.allele_rows.load() << " allele rows for "
         << counters.strands_compared.load() << " called alleles compared with the walk, of which "
         << counters.allele_rows_with_sequence.load() << " give a sequence, "
         << counters.allele_rows_nested.load() << " apply to an earlier allele row's walk, "
         << counters.allele_rows_unwritten.load() << " are at sites with no VCF line, "
         << counters.allele_rows_reversed.load() << " cross their site backwards and "
         << counters.allele_rows_standalone.load() << " join no row" << endl;
    cerr << "[vg call] mosaic: " << counters.strands_not_crossed.load()
         << " called alleles at sites the strand's walk does not cross, and "
         << counters.strands_without_allele.load()
         << " with no walk kept, so no allele row could be written for them" << endl;
}

}

void CalledAlleles::add(const HandleGraph& graph, const PhaseCall& phase,
                        const CandidateAlleles& alleles, bool edit_site) {
    const auto [found, fresh] = entries.try_emplace(phase.record_key);
    if (!fresh) {
        return;
    }
    Entry& entry = found->second;
    entry.edit_site = edit_site;
    const int traversals[2] = {phase.trav_first, phase.ploidy == 2 ? phase.trav_second : -1};
    for (int slot = 0; slot < 2; ++slot) {
        const int traversal = traversals[slot];
        if (traversal < 0 || (size_t)traversal >= alleles.walks.size()) {
            continue;
        }
        entry.traversal[slot] = traversal;
        if (slot == 1 && entry.traversal[0] == traversal) {
            // A homozygous site keeps one copy.
            entry.offset[1] = entry.offset[0];
            entry.length[1] = entry.length[0];
            entry.edits[1] = entry.edits[0];
            continue;
        }
        entry.offset[slot] = nodes.size();
        entry.length[slot] = (uint32_t)alleles.walks[traversal].size();
        for (const handle_t& step : alleles.walks[traversal]) {
            nodes.push_back(gbwt::Node::encode(graph.get_id(step), graph.get_is_reverse(step)));
        }
        if (!alleles.spelled_by_walk(traversal)) {
            entry.edits[slot] = (int32_t)edit_sets.size();
            edit_sets.push_back(alleles.edits[traversal]);
        }
    }
}

bool CalledAlleles::find(size_t record_key, int slot, Allele& out) const {
    const auto found = entries.find(record_key);
    if (found == entries.end() || slot < 0 || slot > 1 || found->second.traversal[slot] < 0) {
        return false;
    }
    const Entry& entry = found->second;
    out.traversal = entry.traversal[slot];
    out.walk = std::span<const gbwt::node_type>(nodes).subspan(entry.offset[slot],
                                                               entry.length[slot]);
    out.edits = entry.edits[slot] >= 0 ? &edit_sets[entry.edits[slot]] : nullptr;
    out.no_change = entry.edit_site && out.edits == nullptr;
    return true;
}

void MosaicWriter::write(const vector<PhaseCall>& phasing, const PanelLookup& panel,
                         const string& sample_name, const vector<PhaseCall>* genotyped,
                         const HandleGraph* graph) {
    ofstream out(params.path);
    if (!out) {
        cerr << "error [vg call]: could not open " << params.path << " for the mosaic output"
             << endl;
        return;
    }
    write_header(out, params, sample_name);

    const size_t reference_hap = reference_haplotype(params);
    cerr << "[vg call] mosaic: reference ";
    if (reference_hap == LinkageModel::NO_HAPLOTYPE) {
        cerr << "is NOT a panel haplotype, so gaps cannot be patched with it";
    } else {
        cerr << "is panel haplotype " << reference_hap;
    }
    cerr << endl;

    // Each contig is written by one thread into a buffer of its own, and the buffers are written
    // out in the phasing's order, so the file is as one pass writes it. Contigs are started
    // largest first, so that the longest ones do not start last.
    const vector<pair<size_t, size_t>> contigs = contig_ranges(phasing);
    vector<size_t> order(contigs.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::stable_sort(order, [&](size_t a, size_t b) {
        return contigs[a].second - contigs[a].first > contigs[b].second - contigs[b].first;
    });
    // An exact mosaic also takes each contig's genotyped sites, which are grouped as the
    // phasing is.
    if (params.exact && (genotyped == nullptr || graph == nullptr)) {
        cerr << "error [vg call]: an exact mosaic needs the genotyped sites and the graph" << endl;
        return;
    }
    unordered_map<string, pair<size_t, size_t>> genotyped_of_contig;
    if (params.exact) {
        for (const auto& [first, past] : contig_ranges(*genotyped)) {
            genotyped_of_contig.emplace((*genotyped)[first].contig, make_pair(first, past));
        }
    }
    vector<string> rows(contigs.size());
    vector<size_t> segments(contigs.size(), 0);
    vector<size_t> unexplained(contigs.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t k = 0; k < order.size(); ++k) {
        const size_t c = order[k];
        ostringstream contig_rows;
        ExactInputs exact;
        if (params.exact) {
            const auto [first, past] = genotyped_of_contig.at(phasing[contigs[c].first].contig);
            exact.genotyped = std::span<const PhaseCall>(*genotyped).subspan(first, past - first);
            exact.alleles = &alleles;
            exact.graph = graph;
        }
        ContigRowWriter writer(
            std::span<const PhaseCall>(phasing).subspan(contigs[c].first,
                                                        contigs[c].second - contigs[c].first),
            panel, reference_hap, params, counters, contig_rows,
            params.exact ? &exact : nullptr);
        writer.write();
        rows[c] = contig_rows.str();
        segments[c] = writer.rows_written();
        unexplained[c] = writer.unexplained_rows();
    }
    size_t total_segments = 0;
    size_t unexplained_segments = 0;
    for (size_t c = 0; c < contigs.size(); ++c) {
        out << rows[c];
        total_segments += segments[c];
        unexplained_segments += unexplained[c];
    }
    cerr << "[vg call] mosaic: " << total_segments << " segments over " << phasing.size()
         << " sites, written to " << params.path << endl;
    cerr << "[vg call] mosaic: " << unexplained_segments
         << " segments the panel cannot name a haplotype for" << endl;
    report_counters(counters);
    if (params.exact) {
        report_allele_rows(counters);
    }
}

}
