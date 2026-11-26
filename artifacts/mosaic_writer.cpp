#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>
#include <map>
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

/// What an exact mosaic's `A` rows are written from, for one contig.
struct ExactInputs {
    /// The phase call of every genotyped site of the contig, with or without a VCF line, in
    /// reference order.
    std::span<const PhaseCall> genotyped;
    /// What each strand carries at each top-level site.
    const StrandAlleles* alleles = nullptr;
    /// The graph, for the sequences of nodes.
    const HandleGraph* graph = nullptr;
};

/**
 * One strand's walk through the graph as its `H` rows spell it: the visit of each node, in walk
 * order.
 *
 * The rows of one fragment that join on a shared node are one piece of the walk. A stretch of the
 * walk lies in one piece.
 */
class RowVisits {
public:
    static constexpr uint32_t NO_VISIT = std::numeric_limits<uint32_t>::max();

    /// One crossing of a site: the visits of its two bound nodes, in walk order.
    using Crossing = pair<uint32_t, uint32_t>;

    /// Forget the walk.
    void clear() {
        node.clear();
        piece.clear();
        piece_fragment.clear();
        index.clear();
    }

    /// Add the walk of an `H` row of fragment `fragment`.
    ///
    /// The row joins the row before it where the two are in one fragment and share the node
    /// they meet at, which is then one visit.
    void append_row(size_t fragment, const vector<gbwt::node_type>& row_nodes) {
        const bool joins = !node.empty() && piece_fragment.back() == fragment
                           && node.back() == row_nodes.front();
        if (!joins) {
            piece_fragment.push_back(fragment);
        }
        node.insert(node.end(), row_nodes.begin() + (joins ? 1 : 0), row_nodes.end());
        piece.resize(node.size(), (uint32_t)(piece_fragment.size() - 1));
    }

    /// Index the visits by node. Call it once, after the last `append_row`.
    void index_visits() {
        index.reserve(node.size());
        for (uint32_t visit = 0; visit < node.size(); ++visit) {
            index.emplace_back(node[visit], visit);
        }
        std::ranges::sort(index);
    }

    /// Count the visits.
    size_t size() const { return node.size(); }

    /// Return the fragment of the row that visit `visit` is from.
    size_t fragment_of(uint32_t visit) const { return piece_fragment[piece[visit]]; }

    /// List the stretches of the walk that run from oriented node `from` to the next visit of
    /// oriented node `to`, without coming back to `from`'s node on the way, in walk order. With
    /// `to` equal to `from`, they are the visits of that node.
    vector<Crossing> crossings(gbwt::node_type from, gbwt::node_type to) const {
        vector<Crossing> found;
        for (uint32_t start = next_visit(from, 0); start != NO_VISIT;
             start = next_visit(from, start + 1)) {
            const uint32_t end = from == to ? start : next_visit(to, start + 1);
            if (end == NO_VISIT || piece[end] != piece[start]) {
                continue;
            }
            if (from == to
                || std::min(next_visit(from, start + 1),
                            next_visit(gbwt::Node::reverse(from), start + 1)) > end) {
                found.emplace_back(start, end);
            }
        }
        return found;
    }

    /// Return whether the walk takes exactly `walk` over a crossing.
    bool spells(const Crossing& crossing, std::span<const gbwt::node_type> walk) const {
        return crossing.second - crossing.first + 1 == walk.size()
               && std::equal(walk.begin(), walk.end(), node.begin() + crossing.first);
    }

private:
    /// Return the first visit of oriented node `visited` at or after visit `from`, or NO_VISIT.
    uint32_t next_visit(gbwt::node_type visited, uint32_t from) const {
        const auto found = std::ranges::lower_bound(index, make_pair(visited, from));
        return found != index.end() && found->first == visited ? found->second : NO_VISIT;
    }

    /// The oriented node of each visit.
    vector<gbwt::node_type> node;
    /// The piece of the walk each visit is in.
    vector<uint32_t> piece;
    /// The fragment of each piece.
    vector<size_t> piece_fragment;
    /// The visits, as (node, visit), sorted.
    vector<pair<gbwt::node_type, uint32_t>> index;
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
    /// With `exact`, each strand's `A` rows are written after its `H` rows.
    ContigRowWriter(std::span<const PhaseCall> sites, const PanelLookup& panel,
                    size_t reference_hap, const MosaicParams& params, MosaicCounters& counters,
                    ostream& out, const ExactInputs* exact = nullptr)
        : sites(sites), walk(panel), reference_hap(reference_hap), params(params),
          counters(counters), out(out), exact(exact) {}

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

    /// Add the walk of the `H` row just printed, which starts at `position` and ends at oriented
    /// node `end`, to `visits`.
    void note_row(gbwt::edge_type position, gbwt::node_type end);

    /// Write the current strand's `A` rows.
    ///
    /// A top-level snarl gets a row where the strand's walk through it, or the bases of a step
    /// of that walk, are not what `visits` has between the snarl's bounds. A node of the
    /// top-level chain gets a row where the strand's edit sites change its bases.
    void write_allele_rows();

    /// An `A` row to print, with the visit its walk starts at.
    using AlleleRow = pair<uint32_t, string>;

    /// Add an `A` row that puts `walk`, with the bases `changed` gives, on the crossings
    /// `found`, for each fragment that has a crossing the row changes.
    ///
    /// @return The number of rows added.
    size_t add_allele_rows(const vector<RowVisits::Crossing>& found,
                           std::span<const gbwt::node_type> walk,
                           std::span<const StrandAlleles::StepBases> changed,
                           vector<AlleleRow>& rows);

    /// The contig's phase calls, in reference order.
    std::span<const PhaseCall> sites;
    /// Finds and follows haplotypes in the panel's GBWT, with caches for this contig only.
    PanelWalker walk;
    /// The reference's index in the panel, or NO_HAPLOTYPE.
    size_t reference_hap;
    const MosaicParams& params;
    MosaicCounters& counters;
    /// Where the contig's rows go.
    ostream& out;
    /// What the `A` rows are written from, or null when the mosaic is not exact.
    const ExactInputs* exact;
    /// The current strand's walk, as its `H` rows spell it. Kept only for an exact mosaic.
    RowVisits visits;
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
        visits.clear();
        collect_strand_sites();
        write_runs();
        if (exact != nullptr) {
            write_allele_rows();
        }
    }
}

void ContigRowWriter::note_row(gbwt::edge_type position, gbwt::node_type end) {
    if (vector<gbwt::node_type> row_nodes; walk.nodes_to(position, end, row_nodes)) {
        visits.append_row(fragment, row_nodes);
    }
}

void ContigRowWriter::write_allele_rows() {
    visits.index_visits();
    const HandleGraph& graph = *exact->graph;
    unordered_set<size_t> written;
    written.reserve(sites.size() * 2);
    for (const PhaseCall& call : sites) {
        written.insert(call.record_key);
    }

    vector<AlleleRow> rows;
    // Whether each visit lies inside a top-level snarl, between its two bounds.
    vector<uint8_t> inside(visits.size(), 0);
    // The changes the strand's edit sites make, by node.
    map<nid_t, vector<AlleleSplicer::StepEdit>> node_edits;
    for (const PhaseCall& call : exact->genotyped) {
        for (const StrandAlleles::NodeChange& change :
             exact->alleles->changes_of(call.record_key, strand)) {
            node_edits[change.node].push_back(AlleleSplicer::StepEdit{
                0, change.at, change.removed, change.added, call.record_key});
        }
        StrandAlleles::Allele allele;
        if (!exact->alleles->find(call.record_key, strand, allele)) {
            continue;
        }
        ++counters.strands_compared;
        const gbwt::node_type start = allele.walk.front();
        const gbwt::node_type end = allele.walk.back();
        const vector<RowVisits::Crossing> found =
            start == end ? vector<RowVisits::Crossing>() : visits.crossings(start, end);
        if (found.empty()) {
            const bool backward =
                start != end && !visits.crossings(gbwt::Node::reverse(end),
                                                  gbwt::Node::reverse(start)).empty();
            ++(backward ? counters.strands_backward : counters.strands_not_crossed);
            continue;
        }
        for (const RowVisits::Crossing& crossing : found) {
            std::fill(inside.begin() + crossing.first + 1, inside.begin() + crossing.second, 1);
        }
        const size_t added = add_allele_rows(found, allele.walk, allele.changed, rows);
        if (written.count(call.record_key) == 0) {
            counters.allele_rows_unwritten += added;
        }
    }

    for (auto& [id, edits] : node_edits) {
        // The changes are made as a record makes them where it spells a strand through a node,
        // so that of two changes that take one base, the same one is made here and there.
        StepSequences bases{graph.get_sequence(graph.get_handle(id))};
        if (!AlleleSplicer::apply_step_edits(std::move(edits), bases)) {
            continue;
        }
        size_t added = 0;
        for (bool reversed : {false, true}) {
            const gbwt::node_type visited = gbwt::Node::encode(id, reversed);
            vector<RowVisits::Crossing> found = visits.crossings(visited, visited);
            std::erase_if(found, [&](const RowVisits::Crossing& at) { return inside[at.first]; });
            const StrandAlleles::StepBases changed{
                0, reversed ? reverse_complement(bases[0]) : bases[0]};
            added += add_allele_rows(found, std::span(&visited, 1), std::span(&changed, 1), rows);
        }
        counters.allele_rows_node += added;
        counters.node_changes_unwritten += (size_t)(added == 0);
    }

    std::ranges::sort(rows);
    for (const AlleleRow& row : rows) {
        out << row.second;
    }
    counters.allele_rows += rows.size();
}

size_t ContigRowWriter::add_allele_rows(const vector<RowVisits::Crossing>& found,
                                        std::span<const gbwt::node_type> walk,
                                        std::span<const StrandAlleles::StepBases> changed,
                                        vector<AlleleRow>& rows) {
    size_t added = 0;
    // A row puts its walk on every crossing of its fragment, so a fragment gets one row.
    size_t last_fragment = std::numeric_limits<size_t>::max();
    for (const RowVisits::Crossing& crossing : found) {
        const size_t row_fragment = visits.fragment_of(crossing.first);
        if (row_fragment == last_fragment || (changed.empty() && visits.spells(crossing, walk))) {
            continue;
        }
        last_fragment = row_fragment;
        ostringstream row;
        row << "A\t" << sites.front().contig << "\t" << strand << "\t" << row_fragment << "\t"
            << walk.front() << "\t" << walk.back() << "\t";
        for (const gbwt::node_type& visited : walk) {
            row << (gbwt::Node::is_reverse(visited) ? '<' : '>') << gbwt::Node::id(visited);
        }
        row << "\t";
        if (changed.empty()) {
            row << ".";
        } else {
            auto next = changed.begin();
            for (size_t step = 0; step < walk.size(); ++step) {
                row << (step > 0 ? "," : "");
                if (next != changed.end() && next->step == step) {
                    row << (next++)->entry;
                } else {
                    row << "*";
                }
            }
            ++counters.allele_rows_with_bases;
        }
        row << "\n";
        rows.emplace_back(crossing.first, row.str());
        ++added;
    }
    return added;
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
        out << "#note\tan A row gives a stretch of a strand where the strand is not what its H "
            << "rows spell. start_node and end_node are oriented as in H rows. walk is the "
            << "strand's node walk from start_node to end_node, both included, written >id for a "
            << "node read forwards and <id for one read backwards. A row is for a top-level site, "
            << "with its two bound nodes as start_node and end_node, or for one node, which is "
            << "both. The strand it gives is the one the VCF's records spell there, with the "
            << "called alleles of the sites nested in the site.\n";
        out << "#note\tto apply an A row, take the walk its fragment's H rows spell, and at each "
            << "visit of start_node find the next visit of end_node. Replace the visits between "
            << "the two with the steps of walk between its first and last. Then give each step of "
            << "walk its entry of bases. No two A rows of a strand change the same visit, so "
            << "they apply in any order.\n";
        out << "#note\tbases is . where walk is a walk of the graph and every step of it has "
            << "the bases of its node. Otherwise it has one entry per step of walk, separated by "
            << "commas: the bases the strand has at that step, read along walk. Within an entry "
            << "* stands for the bases the visit already has, so * alone leaves a step as it is. "
            << "An entry can end with + and bases, or with + alone: those bases, or none, lie "
            << "after the step on no node, and the next step follows them without an edge of the "
            << "graph. A site's row leaves the bases of its two bound nodes to the rows of those "
            << "nodes.\n";
    }
    out << "#H\tcontig\tstrand\tfragment\tref_start\tref_end\tstart_node\tend_node"
        << "\thap_index\thaplotype\tsites\tgbwt_offset\n";
    if (params.exact) {
        out << "#A\tcontig\tstrand\tfragment\tstart_node\tend_node\twalk\tbases\n";
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
void report_allele_rows(const MosaicCounters& counters, const StrandAlleles& alleles) {
    cerr << "[vg call] mosaic: " << counters.allele_rows.load() << " allele rows, of which "
         << counters.allele_rows_node.load() << " are of one node, "
         << counters.allele_rows_with_bases.load() << " give bases and "
         << counters.allele_rows_unwritten.load() << " are at sites with no VCF line; "
         << counters.strands_compared.load()
         << " strands of top-level sites compared with the walk" << endl;
    cerr << "[vg call] mosaic: no allele row could be written for "
         << counters.strands_not_crossed.load()
         << " strands of sites the strand's walk does not cross, "
         << counters.strands_backward.load() << " of sites it crosses backwards, "
         << alleles.bound_changed << " that change a bound node of their site, and "
         << counters.node_changes_unwritten.load()
         << " nodes an edit site changes inside a site or off the walk; "
         << alleles.unphased_parents << " sites with sites under them have no phase" << endl;
}

}

bool StrandAlleles::add_site(const HandleGraph& graph, size_t record_key, int strand,
                             const Traversal& walk, const StepSequences& steps) {
    vector<StepBases> changed;
    for (size_t i = 0; i < steps.size(); ++i) {
        const string own = graph.get_sequence(walk[i]);
        const bool last = i + 1 == steps.size();
        // Whether the strand goes on to the next step without an edge of the graph, as where
        // an allele built from the reads joins the two bound nodes of its site.
        const bool joined = !last && !graph.has_edge(walk[i], walk[i + 1]);
        const bool kept = steps[i].starts_with(own);
        if (steps[i] == own && !joined) {
            continue;
        }
        // A bound node keeps its bases. Only the start bound can have bases put after it.
        if ((i == 0 || last) && (last || !kept)) {
            ++bound_changed;
            return false;
        }
        string entry = kept ? "*" : steps[i];
        if (joined || (kept && steps[i].size() > own.size())) {
            entry += "+" + (kept ? steps[i].substr(own.size()) : string());
        }
        changed.push_back(StepBases{(uint32_t)i, std::move(entry)});
    }
    SiteEntry& entry = sites[record_key];
    entry.walk[strand] = Slice{nodes.size(), walk.size()};
    for (const handle_t& step : walk) {
        nodes.push_back(gbwt::Node::encode(graph.get_id(step), graph.get_is_reverse(step)));
    }
    entry.changed[strand] = Slice{step_bases.size(), changed.size()};
    std::ranges::move(changed, std::back_inserter(step_bases));
    return true;
}

void StrandAlleles::add_edit_site(const HandleGraph& graph, size_t record_key, int strand,
                                  const Traversal& walk, const vector<AlleleEdit>& edits) {
    Slice& kept = edit_sites[record_key][strand];
    kept.at = node_changes.size();
    for (const AlleleEdit& edit : edits) {
        if (edit.step >= walk.size()) {
            continue;
        }
        const handle_t visited = walk[edit.step];
        const size_t length = graph.get_length(visited);
        if (edit.at + edit.removed > length) {
            continue;
        }
        // The walk may read the node backwards, and so may another site that changes the same
        // node, so the change is kept along the node's forward strand.
        const bool reversed = graph.get_is_reverse(visited);
        node_changes.push_back(NodeChange{
            graph.get_id(visited), reversed ? length - edit.at - edit.removed : edit.at,
            edit.removed, reversed ? reverse_complement(edit.added) : edit.added});
    }
    kept.size = node_changes.size() - kept.at;
}

bool StrandAlleles::find(size_t record_key, int strand, Allele& out) const {
    const auto found = sites.find(record_key);
    if (found == sites.end() || found->second.walk[strand].size == 0) {
        return false;
    }
    const SiteEntry& entry = found->second;
    out.walk = std::span(nodes).subspan(entry.walk[strand].at, entry.walk[strand].size);
    out.changed = std::span(step_bases).subspan(entry.changed[strand].at,
                                                entry.changed[strand].size);
    return true;
}

std::span<const StrandAlleles::NodeChange> StrandAlleles::changes_of(size_t record_key,
                                                                     int strand) const {
    const auto found = edit_sites.find(record_key);
    if (found == edit_sites.end()) {
        return {};
    }
    return std::span(node_changes).subspan(found->second[strand].at, found->second[strand].size);
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
        report_allele_rows(counters, alleles);
    }
}

}
