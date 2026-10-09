#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "mosaic_writer.hpp"
#include "vcf_output_caller.hpp"
#include "graph_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

namespace {

/// How far a walk may run before it is abandoned. Walks go only in a direction already
/// established, so this limits a long run rather than a wrong-way search.
const size_t WALK_LIMIT = 1u << 17;

/// Finds and follows panel haplotypes in the GBWT, for one mosaic. The mosaic is written
/// serially, so one position cache with no locking is enough, and it lives only for that mosaic.
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
        // are cached.
        const uint64_t key = ((uint64_t)node << 20) | (uint64_t)(hap & 0xFFFFF);
        auto hit = position_cache.find(key);
        if (hit != position_cache.end()) {
            return hit->second;
        }
        gbwt::SearchState state = index->find(node);
        if (!state.empty()) {
            for (gbwt::size_type i = state.range.first; i <= state.range.second; ++i) {
                gbwt::size_type seq = index->locate(node, i);
                if (seq < haplotype_of_sequence->size() && (*haplotype_of_sequence)[seq] == hap) {
                    position_cache[key] = gbwt::edge_type(node, i);
                    return gbwt::edge_type(node, i);
                }
            }
        }
        position_cache[key] = gbwt::invalid_edge();
        return gbwt::invalid_edge();
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
    const gbwt::GBWT* index;
    const vector<size_t>* haplotype_of_sequence;
    /// (oriented node, haplotype) -> GBWT position.
    unordered_map<uint64_t, gbwt::edge_type> position_cache;
};

}

void MosaicWriter::write(const vector<LinkageCollector::PhaseCall>& phasing,
                         const PanelLookup& panel, const string& sample_name) {
    PanelWalker walk(panel);
    ofstream out(params.path);
    if (!out) {
        cerr << "error [vg call]: could not open " << params.path << " for the mosaic output"
             << endl;
        return;
    }

    // A segment is a maximal run over which one strand stays on one panel haplotype. A consumer
    // rebuilds a haplotype by walking it from start_node to end_node, so segments are located by
    // node ID rather than by reference position.
    //
    // The header states two things the rows do not: which reference the positions are in, since
    // a graph can hold several references with the same contig names, and what each hap_index
    // means, since the index is internal to the run. A haplotype's name is its (sample, phase)
    // pair, which with the row's contig is enough to find its paths.
    out << "#mosaic-version\t5\n";
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
    // The node IDs define the segment; the positions are derived from them.
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
    out << "#H\tcontig\tstrand\tfragment\tref_start\tref_end\tstart_node\tend_node"
        << "\thap_index\thaplotype\tsites\tgbwt_offset\n";

    // The phasing is grouped by contig and in reference order. Each strand is written separately,
    // and a switch on one strand ends only that strand's segment.
    size_t i = 0;
    size_t total_segments = 0;

    // Emit sites [from, to] on one strand, all on haplotype `hap`, as one row per GBWT fragment.
    //
    // A row carries one GBWT position, from which the whole segment must be walkable, so a run is
    // cut wherever the fragment under it changes. We resolve the positions at the run's two ends;
    // only when they are on different fragments do we binary-search the sites for the boundary.
    // This misses a haplotype that leaves a fragment and comes back to it within one run, which
    // fragments of one path cannot do in reference order.
    //
    // What a strand holds at a site: its allele on a known haplotype (Carried); nothing, because it
    // is the other strand of a nested ploidy-1 site (Empty); or sequence the panel cannot attribute
    // to a haplotype (Unexplained). A run is cut where the kind changes, so each row has one.
    enum class StrandKind { Carried, Empty, Unexplained };
    auto strand_kind = [&](size_t t, int strand) -> StrandKind {
        const size_t hap = strand == 0 ? phasing[t].hap_first : phasing[t].hap_second;
        if (hap != LinkageModel::WILDCARD) {
            return StrandKind::Carried;
        }
        if (phasing[t].nested_strand >= 0 && (int)phasing[t].nested_strand != strand) {
            return StrandKind::Empty;
        }
        return StrandKind::Unexplained;
    };
    size_t unexplained_segments = 0;
    // The sites this strand passes through, as indices into `phasing`, in reference order. A nested
    // ploidy-1 chain is on one of its parent's strands; the other strand takes the parent's other
    // allele, which bypasses the chain, so the site is not on that strand's walk and has no row
    // there. `emit_span` and `emit_row` index into this list; `site()` maps back to `phasing`.
    //
    // The reference's index in the panel, for filling a gap no haplotype can cross. Looked up by
    // name, since the index follows GBWT metadata order. `params.reference_paths` holds full path
    // names (CHM13#0#chr20) and the panel names haplotypes as sample#phase (CHM13#0), so the
    // contig is dropped before matching.
    size_t reference_hap = LinkageModel::WILDCARD;
    for (const string& full : params.reference_paths) {
        // Never a gRef path: a gRef cover is stitched together from many donors, and it is not in the
        // panel.
        if (GrefCover::is_gref_derived(full)) {
            continue;
        }
        size_t h1 = full.find('#');
        size_t h2 = h1 == string::npos ? string::npos : full.find('#', h1 + 1);
        const string base = h2 == string::npos ? full : full.substr(0, h2);
        for (size_t k = 0; k < params.haplotype_names.size(); ++k) {
            if (params.haplotype_names[k] == base) {
                reference_hap = k;
                break;
            }
        }
        if (reference_hap != LinkageModel::WILDCARD) {
            break;
        }
    }
    cerr << "[vg call] mosaic: reference "
         << (reference_hap == LinkageModel::WILDCARD
                 ? string("is NOT a panel haplotype, so gaps cannot be patched with it")
                 : "is panel haplotype " + std::to_string(reference_hap))
         << endl;

    const bool patch_gaps = params.patch_gaps;
    const bool keep_nested = params.keep_nested;
    const bool connect_unexplained = params.connect_unexplained;
    // Which contiguous walk a row belongs to. (contig, strand, fragment) identifies a path: a loader
    // makes one path per triple from its rows. Incremented only where a gap is left unfilled, so
    // with gap patching each strand is one path.
    size_t fragment = 0;
    vector<size_t> strand_sites;
    auto site = [&](size_t pos) -> const LinkageCollector::PhaseCall& {
        return phasing[strand_sites[pos]];
    };
    // A left extension made by the previous row: this segment begins at the previous segment's last
    // node rather than at its own first site, and its GBWT position moves with it. Reset for each
    // strand.
    int64_t pending_from_node = -1;
    gbwt::edge_type pending_from_pos = gbwt::invalid_edge();
    // The oriented node this strand's walk has reached, which carries the walk's direction forward.
    // It is found once per strand and then followed, so consecutive rows join by construction.
    gbwt::node_type carry = gbwt::ENDMARKER;

    std::function<void(size_t, size_t, int, size_t, StrandKind)> emit_span =
        [&](size_t from, size_t to, int strand, size_t hap, StrandKind kind) {
        gbwt::edge_type pos = (hap == LinkageModel::WILDCARD)
                                  ? gbwt::invalid_edge()
                                  : walk.gbwt_position(site(from).start_node, hap);

        auto emit_row = [&](size_t a_idx, size_t b_idx, gbwt::edge_type p) {
            const LinkageCollector::PhaseCall& a = site(a_idx);
            const LinkageCollector::PhaseCall& b = site(b_idx);
            // Right extension: a segment ends where the next one begins, rather than at its own last
            // site's end, so that the stretch between two segments is covered. Nothing there shows
            // which of the two haplotypes covers it, so extending rightward is a convention; the
            // header says so. The extension is made only if the haplotype reaches the next
            // segment's first node on the same GBWT fragment; otherwise the row ends at its own last
            // site and the gap is counted, for the caller to patch or break.
            //
            // A left extension made by the previous row moves this segment's start back, and its
            // position with it.
            int64_t from_node = a.start_node;
            if (pending_from_node >= 0) {
                from_node = pending_from_node;
                p = pending_from_pos;
            }
            pending_from_node = -1;

            int64_t to_node = b.end_node;
            // A gap this row could not close, and the reference stretch that fills it, written just
            // after this row.
            int64_t patch_to = -1;
            gbwt::edge_type patch_pos = gbwt::invalid_edge();
            size_t patch_from_pos = 0, patch_to_pos = 0;
            bool boundary_open = false;
            if (b_idx + 1 < strand_sites.size()) {
                const LinkageCollector::PhaseCall& nx = site(b_idx + 1);
                const int64_t next_start = nx.start_node;
                // A nested snarl is contained in its parent, so the parent's walk is
                // Ps -> ... -> Cs -> [child] -> Ce -> ... -> Pe, and a change of haplotype between a
                // parent and a child has two boundaries, both the child's: Cs, where the walk enters
                // the child, and Ce, where it leaves. Every traversal of the child passes through Cs
                // and Ce, so the two haplotypes meet there.
                //
                // Depth alone decides entering and leaving. Sites arrive in reference order, and a
                // parent is always recorded, so a site deeper than the one before it is inside that
                // one. Comparing node IDs would fail wherever IDs do not follow the walk, as in an
                // inversion.
                const bool entering = nx.level > b.level;
                const bool leaving = nx.level < b.level;
                // A right extension needs this segment to be walkable.
                bool right = false;
                if (p != gbwt::invalid_edge()) {
                    const gbwt::edge_type np = walk.gbwt_position(next_start, hap);
                    // The next segment's haplotype must pass the junction in the same direction as
                    // this one, not only through the same node.
                    const size_t nh2 = strand == 0 ? site(b_idx + 1).hap_first
                                                   : site(b_idx + 1).hap_second;
                    const gbwt::edge_type entry = nh2 == LinkageModel::WILDCARD
                                                      ? gbwt::invalid_edge()
                                                      : walk.gbwt_position(next_start, nh2);
                    right = np != gbwt::invalid_edge()
                            && walk.gbwt()->locate(np) == walk.gbwt()->locate(p)
                            && (entry == gbwt::invalid_edge() || entry.first == np.first);
                }
                if (entering) {
                    // The row ends where the child's snarl begins, since the walk reaches the
                    // parent's end only after the child.
                    to_node = next_start;
                    ++counters.nested_enter;
                } else if (leaving) {
                    // The row ends at the child's own end, and the next row starts there, so the
                    // stretch from Ce to the parent's end is covered by the parent's haplotype, whose
                    // called allele governs it.
                    to_node = b.end_node;
                    pending_from_node = b.end_node;
                    const size_t nh = strand == 0 ? nx.hap_first : nx.hap_second;
                    pending_from_pos = nh == LinkageModel::WILDCARD
                                           ? gbwt::invalid_edge()
                                           : walk.gbwt_position(b.end_node, nh);
                    ++counters.nested_leave;
                } else if (right) {
                    to_node = next_start;
                    ++counters.extended;
                } else {
                    // Left extension: this segment's haplotype cannot be carried forward, so try
                    // carrying the next segment's haplotype back to this segment's last node, which
                    // closes the gap with a panel haplotype rather than the reference.
                    const size_t nh = strand == 0 ? site(b_idx + 1).hap_first
                                                  : site(b_idx + 1).hap_second;
                    bool closed = false;
                    if (nh != LinkageModel::WILDCARD) {
                        const gbwt::edge_type here = walk.gbwt_position(b.end_node, nh);
                        const gbwt::edge_type there = walk.gbwt_position(next_start, nh);
                        const gbwt::edge_type mine = walk.gbwt_position(b.end_node, hap);
                        if (here != gbwt::invalid_edge() && there != gbwt::invalid_edge()
                            && walk.gbwt()->locate(here) == walk.gbwt()->locate(there)
                            && (mine == gbwt::invalid_edge() || mine.first == here.first)) {
                            pending_from_node = b.end_node;
                            pending_from_pos = here;
                            ++counters.extended_left;
                            closed = true;
                        }
                    }
                    if (!closed) {
                        // Neither haplotype crosses the gap, so fill it with the reference, if it
                        // crosses. The fill is contiguous but says little about the sample, so the
                        // row records what it filled.
                        if (patch_gaps && reference_hap != LinkageModel::WILDCARD) {
                            const gbwt::edge_type rl =
                                walk.gbwt_position(b.end_node, reference_hap);
                            const gbwt::edge_type rr =
                                walk.gbwt_position(next_start, reference_hap);
                            if (rl != gbwt::invalid_edge() && rr != gbwt::invalid_edge()
                                && walk.gbwt()->locate(rl) == walk.gbwt()->locate(rr)) {
                                patch_to = next_start;
                                patch_pos = rl;
                                patch_from_pos = b.position;
                                patch_to_pos = site(b_idx + 1).position;
                                ++counters.patched;
                            }
                        }
                        if (patch_to < 0) {
                            ++counters.gap_left;
                            boundary_open = true;
                        }
                    }
                }
            }
            // A row names its haplotype only if that haplotype crosses the row's whole node range on
            // one GBWT fragment. Otherwise the row cannot be walked, so it becomes a reference
            // substitution, marked `ref` like a gap fill but keeping its site count, since it covers
            // called sites.
            //
            // The walk's start is oriented. The carried direction applies when the previous row
            // ended at this node; otherwise, for a strand's first row or a moved start, both
            // orientations are tried and the one that reaches the row's far end is kept. Each
            // position found is reused for the walk.
            const auto start_and_walk = [&](size_t h, gbwt::edge_type* pos,
                                            gbwt::node_type* end,
                                            bool ignore_carry = false) -> bool {
                // Where the carried direction applies, it decides: if the haplotype cannot be
                // followed from it, the row does not continue the walk, rather than taking the
                // other orientation and breaking contiguity.
                if (!ignore_carry && carry != gbwt::ENDMARKER
                    && (int64_t)gbwt::Node::id(carry) == from_node) {
                    const gbwt::edge_type at = walk.position_at(carry, h);
                    if (at == gbwt::invalid_edge() || !walk.follow(at, to_node, end)) {
                        return false;
                    }
                    *pos = at;
                    return true;
                }
                // No direction yet, at the strand's first row: try both.
                for (int o = 0; o < 2; ++o) {
                    const gbwt::edge_type at =
                        walk.position_at(gbwt::Node::encode(from_node, o == 1), h);
                    if (at != gbwt::invalid_edge() && walk.follow(at, to_node, end)) {
                        *pos = at;
                        return true;
                    }
                }
                return false;
            };
            gbwt::edge_type row_pos = gbwt::invalid_edge();
            gbwt::node_type row_end = gbwt::Node::encode(to_node, false);
            bool as_ref = false;
            bool walkable = false;
            // Whether the carried direction constrains this row. It does not for a strand's first row,
            // or where an extension moved the start.
            const bool carry_applies = carry != gbwt::ENDMARKER
                                       && (int64_t)gbwt::Node::id(carry) == from_node;
            if (hap != LinkageModel::WILDCARD) {
                walkable = start_and_walk(hap, &row_pos, &row_end);
                if (!walkable && patch_gaps && reference_hap != LinkageModel::WILDCARD
                    && start_and_walk(reference_hap, &row_pos, &row_end)) {
                    as_ref = true;
                    walkable = true;
                    ++counters.row_to_ref;
                }
            }
            // Last resort: try the row without the carried direction. A row about to have no
            // position has no contiguity left to break, and a walk in the other direction is better
            // than none. This happens at an inversion, whose ends the haplotype traverses in
            // reverse; such a row cannot join the row before it, so the fragment breaks there.
            bool direction_broken = false;
            if (!walkable && hap != LinkageModel::WILDCARD
                && start_and_walk(hap, &row_pos, &row_end, true)) {
                walkable = true;
                direction_broken = carry_applies && row_pos.first != carry;
                if (direction_broken) {
                    ++counters.direction_broken;
                }
            }
            // The same last resort for the reference substitution. A panel haplotype can be clipped
            // across a site, which the linkage model allows, so the row falls back to the
            // reference, and the carried direction, inherited from an inverted row, may not match
            // the reference's.
            if (!walkable && patch_gaps && hap != LinkageModel::WILDCARD
                && reference_hap != LinkageModel::WILDCARD
                && start_and_walk(reference_hap, &row_pos, &row_end, true)) {
                as_ref = true;
                walkable = true;
                ++counters.row_to_ref;
                direction_broken = carry_applies && row_pos.first != carry;
                if (direction_broken) {
                    ++counters.direction_broken;
                }
            }
            if (!walkable) {
                row_pos = gbwt::invalid_edge();
                row_end = gbwt::Node::encode(to_node, false);
            }
            carry = walkable ? row_end : gbwt::ENDMARKER;
            // A row whose direction was broken starts a new fragment, since it cannot join the row
            // before it; the row after it can join it as usual. A row with no position stands
            // alone, breaking on both sides.
            if (hap == LinkageModel::WILDCARD || !walkable || direction_broken) {
                ++fragment;
            }
            // Oriented node IDs, `id * 2 + is_reverse`, as vg encodes them, since two segments can
            // share a node and pass it in opposite directions. The orientation comes from the
            // resolved position; without one, as for an unexplained row, the reference orientation
            // is used.
            const gbwt::node_type row_start = row_pos != gbwt::invalid_edge()
                                                  ? row_pos.first
                                                  : gbwt::Node::encode(from_node, false);
            out << "H\t" << a.contig << "\t" << strand << "\t" << fragment << "\t"
                << a.position << "\t" << b.position << "\t"
                << row_start << "\t" << row_end << "\t";
            if (as_ref) {
                out << "ref\t"
                    << (reference_hap < params.haplotype_names.size()
                            ? params.haplotype_names[reference_hap] : string("?"));
            } else if (hap == LinkageModel::WILDCARD) {
                // The strand passes through here, and the panel cannot name a haplotype for it.
                out << "*\t*";
                ++unexplained_segments;
            } else {
                out << hap << "\t"
                    << (hap < params.haplotype_names.size()
                            ? params.haplotype_names[hap] : string("?"));
            }
            out << "\t" << (b_idx - a_idx + 1) << "\t";
            if (row_pos == gbwt::invalid_edge()) {
                // No position: the strand is the wildcard, or the haplotype does not cross this run
                // in the graph. "." rather than 0, which would be a valid offset. Panel haplotypes
                // are often clipped, so the second case is common.
                out << ".";
            } else {
                out << row_pos.second;
            }
            out << "\n";
            ++total_segments;

            // The reference fill, between the two segments it joins. Its haplotype is `ref`, so a
            // consumer can tell it from a panel haplotype, and its site columns are ".", since it
            // covers no called site.
            if (patch_to >= 0) {
                // The fill is a walk too, starting where this row ended, so the carried direction
                // runs through it.
                gbwt::edge_type ps = gbwt::invalid_edge();
                gbwt::node_type pe = gbwt::Node::encode(patch_to, false);
                const gbwt::node_type s = carry != gbwt::ENDMARKER
                                              ? carry
                                              : gbwt::Node::encode(b.end_node, false);
                ps = walk.position_at(s, reference_hap);
                if (ps != gbwt::invalid_edge() && walk.follow(ps, patch_to, &pe)) {
                    out << "H\t" << a.contig << "\t" << strand << "\t" << fragment << "\t"
                        << patch_from_pos << "\t" << patch_to_pos << "\t"
                        << ps.first << "\t" << pe << "\t"
                        << "ref\t"
                        << (reference_hap < params.haplotype_names.size()
                                ? params.haplotype_names[reference_hap] : string("?"))
                        << "\t.\t" << ps.second << "\n";
                    ++total_segments;
                    carry = pe;
                } else {
                    ++counters.gap_left;
                    --counters.patched;
                    carry = gbwt::ENDMARKER;
                    ++fragment;
                }
            }
            // Only a gap left open ends the fragment. A left extension closes the gap by moving the
            // next row's start, so this row's end is unchanged and the gap is not open.
            if (boundary_open) {
                ++fragment;
            }
            // A row with no position ends the fragment, since no consumer could walk across it.
            if (hap == LinkageModel::WILDCARD || !walkable) {
                ++fragment;
            }
        };

        if (pos == gbwt::invalid_edge() && hap != LinkageModel::WILDCARD && from != to) {
            // The run's first site is not in the graph for this haplotype, but a later one may be, so
            // find the first site that resolves and write the walkable rest separately, rather than
            // giving up on the run or patching sites whose haplotype is known.
            size_t first_ok = from;
            while (first_ok <= to
                   && walk.gbwt_position(site(first_ok).start_node, hap)
                          == gbwt::invalid_edge()) {
                ++first_ok;
            }
            if (first_ok > to) {
                emit_row(from, to, pos);        // clipped across the whole run
                ++counters.unwalkable;
                return;
            }
            // The unresolvable head, as small as it really is, then the walkable remainder.
            if (first_ok > from) {
                emit_row(from, first_ok - 1, gbwt::invalid_edge());
                ++counters.unwalkable;
                ++counters.head_clipped;
            }
            emit_span(first_ok, to, strand, hap, kind);
            return;
        }
        if (pos == gbwt::invalid_edge() || from == to) {
            if (pos == gbwt::invalid_edge() && hap != LinkageModel::WILDCARD) {
                ++counters.unwalkable;
            }
            emit_row(from, to, pos);
            return;
        }
        gbwt::edge_type end_pos = walk.gbwt_position(site(to).start_node, hap);
        if (end_pos == gbwt::invalid_edge()
            || walk.gbwt()->locate(pos) == walk.gbwt()->locate(end_pos)) {
            // Same fragment at both ends, or no way to tell. One row.
            emit_row(from, to, pos);
            return;
        }
        // The fragment changes somewhere in (from, to]. Binary search for the last site still on
        // the starting fragment; a site the haplotype does not reach is treated as past the
        // boundary, which keeps the search monotone.
        gbwt::size_type seq = walk.gbwt()->locate(pos);
        size_t lo = from, hi = to;
        while (hi - lo > 1) {
            size_t mid = lo + (hi - lo) / 2;
            gbwt::edge_type p = walk.gbwt_position(site(mid).start_node, hap);
            if (p != gbwt::invalid_edge() && walk.gbwt()->locate(p) == seq) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        emit_row(from, lo, pos);
        emit_span(hi, to, strand, hap, kind);
    };

    while (i < phasing.size()) {
        size_t j = i;
        while (j < phasing.size() && phasing[j].contig == phasing[i].contig) {
            ++j;
        }
        // One strand on a haploid chain, since a second would claim a copy the sample lacks. Taken
        // over the whole run rather than from its first site, since a diploid contig can begin with
        // a nested ploidy-1 site.
        int strands = 1;
        for (size_t t = i; t < j && strands == 1; ++t) {
            if (phasing[t].ploidy != 1) {
                strands = 2;
            }
        }
        for (int strand = 0; strand < strands; ++strand) {
            pending_from_node = -1;
            carry = gbwt::ENDMARKER;
            fragment = 0;
            // A site this strand does not traverse is not on its walk, so it is not in its list.
            strand_sites.clear();
            for (size_t t = i; t < j; ++t) {
                if (strand_kind(t, strand) == StrandKind::Empty) {
                    continue;
                }
                // A switch of haplotype inside a nested chain is real, but it follows the parent's
                // route and a consumer counting recombinations may not want it. Dropping the nested
                // sites merges the runs across them, so the walk follows the parent's haplotype
                // through the child snarl. The header records which was done.
                if (!keep_nested && phasing[t].level > 0) {
                    continue;
                }
                // A stretch the panel cannot explain with few switches: the wildcard. Its alleles may
                // all be carried by panel haplotypes; what is missing is a panel walk through the
                // stretch. By default the flanking haplotype is carried through, keeping the strand
                // one path, at the cost of writing that haplotype's sequence across those sites
                // rather than the called alleles. With `connect_unexplained` off, the hole is left.
                if (connect_unexplained
                    && strand_kind(t, strand) == StrandKind::Unexplained) {
                    continue;
                }
                strand_sites.push_back(t);
            }
            size_t seg_start = 0;
            for (size_t t = 0; t < strand_sites.size(); ++t) {
                size_t hap = strand == 0 ? site(t).hap_first : site(t).hap_second;
                StrandKind kind = strand_kind(strand_sites[t], strand);
                bool last = (t + 1 == strand_sites.size());
                // Cut where the kind changes too, so that each run has one.
                bool changes = !last
                               && ((strand == 0 ? site(t + 1).hap_first
                                                : site(t + 1).hap_second) != hap
                                   || strand_kind(strand_sites[t + 1], strand) != kind);
                if (last || changes) {
                    // One run of one haplotype, possibly over several GBWT fragments; emit_span
                    // writes one row per fragment, so each row can be walked from its position.
                    emit_span(seg_start, t, strand, hap, kind);
                    seg_start = t + 1;
                }
            }
        }
        i = j;
    }
    cerr << "[vg call] mosaic: " << total_segments << " segments over " << phasing.size()
         << " sites, written to " << params.path << endl;
    // Empty and unexplained segments are reported separately. The unexplained count compares with
    // the phasing report's.
    cerr << "[vg call] mosaic: " << unexplained_segments
         << " segments the panel cannot name a haplotype for" << endl;
    // Segments naming a haplotype the graph does not carry across them, which a consumer has to
    // patch or break at.
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

}
