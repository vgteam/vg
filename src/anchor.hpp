#ifndef VG_ANCHOR_HPP_INCLUDED
#define VG_ANCHOR_HPP_INCLUDED

/** \file anchor.hpp
 *
 * Anchors for pangenome-guided assembly.
 *
 * An anchor is a **zero-length pin**: a point *between* two adjacent positions, carrying the set of
 * reads that cross it and, for each, where it sits in that read. There is no anchor sequence, so
 * there is no question of whether a read matches one, no minimum length, and no way for the members
 * to disagree about what the anchor spells.
 *
 * Two pins per site, both canonical -- no data-dependent choice and no sweep over the members:
 *
 *   - the **S** pin is at the junction between the snarl's start boundary node and the site
 *     interior: it sits immediately after that node's last base, reading in the site's direction;
 *   - the **E** pin is at the junction between the site interior and the end boundary node: it sits
 *     immediately before that node's first base, reading in the site's direction.
 *
 * ## A pin has a side, and that is what makes it unique
 *
 * A pin is not "at node n" but at one *end* of n. Two snarls sharing a boundary node in a chain
 * therefore pin at opposite ends of it and never coincide, whatever the node's length. Measured on
 * chr20 of a 34-haplotype HPRC graph: all 529,304 pins land on 529,304 distinct (node, side) sites,
 * not one claimed by two snarls. Drop the side and pin at the node instead and 37.3% of pins
 * collide -- 98,715 boundary nodes are shared by two snarls. 12.9% of boundary nodes are exactly
 * 1 bp, which is why this is stated outright rather than left implicit: a 1 bp node has no interior
 * position to pin at, only its two ends, and those still belong to two different snarls.
 *
 * ## One position, one anchor
 *
 * A read position is a junction between two consecutive nodes of the read's walk. Such a junction
 * can host at most two pins: the S pin of a snarl starting at the left node with its interior facing
 * right, and the E pin of a snarl ending at the right node with its interior facing left. Each
 * (node, side) has exactly one owning snarl, so if a junction hosts both, they belong to the *same*
 * snarl. Cross-snarl collision is therefore impossible, and one case remains: a read whose walk goes
 * directly from a snarl's start boundary to its end boundary carries an allele that deletes the
 * whole site, so it has no bases inside it and that snarl's own two pins land on the same read
 * offset. 19.8% of chr20's snarls have such an edge. Such a read is kept in the S anchor and dropped
 * from the E anchor; see `AnchorCounters::coincident`.
 */

#include <atomic>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "handle.hpp"
#include "site_read_source.hpp"
#include "vg/vg.pb.h"

namespace vg {

using namespace std;

/**
 * Where one read crosses one pin.
 *
 * `offset` is the 0-based index, **in the read as sequenced**, of the last base before the pin,
 * reading in the site's direction. So for `strand == 0` the pin lies immediately *after*
 * `read[offset]`, and for `strand == 1` immediately *before* it.
 *
 * The read-as-sequenced frame is deliberate. vg expresses reverse-strand alignment through
 * reverse-oriented node visits rather than by reverse-complementing the stored sequence (libvgio
 * writes `gaf.strand = '+'`, "always positive relative to the path"), so this is the frame the read
 * file is in and the frame every offset in vg already uses. No conversion happens inside vg, which
 * is where a strand bug would otherwise live. A consumer wanting the oriented read does
 * `len - 1 - offset` itself.
 */
struct AnchorPlacement {
    /// Negative when the read does not cross this pin.
    int64_t offset = -1;
    /// 0 when the read crosses the pin in the site's direction, 1 when against it.
    uint8_t strand = 0;
    bool placed() const { return offset >= 0; }
};

/// Why a read could not be pinned, and whether the invariant held. Shared across threads.
struct AnchorCounters {
    /// Pins resolved and checked against the graph's own base.
    atomic<size_t> verified{0};
    /// **Must stay zero.** A coordinate or strand bug; see `resolve_anchor_pin`.
    atomic<size_t> verify_failed{0};
    /// The read's alignment never visits the pin's node.
    atomic<size_t> no_visit{0};
    /// It visits it more than once, so which visit the pin belongs to is ambiguous.
    atomic<size_t> repeat_visit{0};
    /// The node base the pin is defined against has no read base aligned to it -- a deletion covers
    /// it, or the alignment stops short of the node's end. Walking back to an earlier base is
    /// exactly what must not happen: on a 1 bp boundary node it would walk straight past the node
    /// into the neighbouring snarl's pin position, which is how a single deleted base would
    /// manufacture a shared pin.
    atomic<size_t> unaligned_base{0};
    /// E pin only: the read has no base immediately upstream on its own walk -- it begins here, or
    /// the node one step upstream is deleted outright in this read, so there is no base adjacent to
    /// the pin. Refused rather than resolved to whatever lies beyond that node, which would be the
    /// neighbouring snarl's pin position.
    atomic<size_t> no_neighbour{0};
    /// A read whose two pins at one snarl resolved to the same position; kept at S, dropped at E.
    atomic<size_t> coincident{0};
    /// Reads excluded because their best-fitting allele was not a called one.
    atomic<size_t> off_call{0};
    /// Snarls whose two boundary nodes are the same node, where only the S pin is emitted.
    atomic<size_t> degenerate_site{0};
    /// End anchors suppressed because every read on them is also on their slot's start anchor.
    atomic<size_t> single_pin{0};
    /// Reads sharing a name with another read at the same site. Paired-end mates share a name, so
    /// this is a property of the read file rather than a fault -- but it means the name does not
    /// identify an alignment, and a consumer keying on it will merge two different reads. Reported
    /// because it is invisible otherwise, and because long reads, the target here, are unpaired.
    atomic<size_t> shared_name{0};

    void report(ostream& out) const;
};

/**
 * Resolve where one read crosses one pin, and check the result against the graph.
 *
 * `site_backward` is the orientation the site's alleles visit `node_id` in, which is the snarl's own
 * start or end visit -- the anchor's direction is the site's direction, not the node's arbitrary
 * one. `exit_pin` selects the S pin (the junction leaving the node, so the reference base is the
 * node's last base in the site's direction) from the E pin (the junction entering it, reference base
 * the node's first base).
 *
 * The invariant is checked here, where the read sequence is live, because the read is gone by the
 * time anchors are written: the read base the pin is resolved against, complemented when the read
 * visits the node in reverse, must equal the graph's base. One comparison, and it catches every
 * off-by-one and every strand inversion. `counters.verify_failed` must read zero.
 */
AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters);

/// The run's counters. Process-wide, like the mosaic's, because the sweep is parallel over node-ID
/// windows and nothing about a counter is per site.
AnchorCounters& anchor_counters();

/// One read's contribution to a site's anchors.
struct AnchorRead {
    string name;
    /// The MAPQ-derived mismapping probability, clamped as the genotype model clamps it.
    float mismap = 0.0f;
    AnchorPlacement start_pin;
    AnchorPlacement end_pin;
};

/**
 * Everything a site's anchors need, retained from the sweep to the render.
 *
 * No read is retained: the alignment does not survive the callback it arrives in, and the pin is
 * resolved to `(strand, offset)` while it is live. Deferring that resolution to render is the one
 * change that would force the alignments to stay around, so it is not made.
 *
 * Held on the site's `ReadLikelihoodCallInfo`, which the caller already retains for every record, so
 * this needs no container of its own.
 */
struct AnchorSiteEvidence {
    vector<AnchorRead> reads;
    /// rel(r, a), row major, reads x alleles, row-normalised into [0,1].
    vector<float> rel;
    size_t n_alleles = 0;
    /// The alleles' spelled lengths and the site's mean read length, for the mixture weights the
    /// per-read score uses. Whole traversal length rather than the unique content the genotype
    /// model's weights use: unique content is pairwise, so it cannot be reduced to one number per
    /// allele before the genotype is known.
    vector<uint32_t> allele_length;
    float mean_read_length = 0.0f;
    bool length_weighted = true;
    nid_t start_node = 0;
    nid_t end_node = 0;

    float rel_at(size_t read, size_t allele) const {
        return rel[read * n_alleles + allele];
    }
    /// Retained bytes, for the progress line. Reported rather than estimated.
    size_t bytes() const;
};

/// Thresholds and switches, from the --anchors-* options.
struct AnchorParams {
    bool enabled = false;
    /// Restrict to heterozygous sites, dropping the ones where the pair degenerates to a single
    /// anchor holding every read.
    ///
    /// Off by default, so homozygous and haploid sites DO get an anchor. They carry no haplotype
    /// information -- there is nothing to partition -- but they carry connectivity, and an anchor
    /// graph is built out of contiguity as much as out of phasing. Excluding them leaves the
    /// consumer with anchors only where the sample happens to be heterozygous, which is a sparse and
    /// uneven scaffold.
    bool het_only = false;

    /// Emit a slot's end anchor only where it holds at least this many reads that slot's start
    /// anchor does not. 0 emits it always, which is the default.
    ///
    /// The two pins of a site carry the **same** partition -- membership and slot are decided per
    /// site, not per pin -- so the second pin is never extra haplotype information. What it is, is a
    /// second vertex. And if every read at the end pin is also at the start pin, the two are joined
    /// by *all* the same reads: every linkage the end pin could offer, the start pin already offers,
    /// and an anchor graph sees a trivial chain. Dropping it then cannot remove a read-linkage.
    ///
    /// That is the condition itself rather than a proxy for it. An earlier version thresholded on
    /// the widest called allele's interior -- how far apart the pins *can* get -- which is a graph
    /// property standing in for "do different reads reach them". It measures the same thing at one
    /// remove and mispredicts wherever coverage or read length disagrees with allele length.
    ///
    /// The cost of the proxy is that it is stable and this is not: which reads reach a boundary
    /// depends on depth and read length, so the same graph and sample at different coverage give
    /// different anchor sets. For feeding one run's reads to an assembler that is exactly right;
    /// for comparing anchor sets between runs it is not, and the length rule would have been.
    size_t end_pin_min_new = 0;
    /// Only leaf snarls, which is closest to the published construction.
    bool leaf_only = false;
    size_t min_reads = 2;
    double min_gqn = 0.0;
    double min_read_score = 0.0;
    /// Keep reads whose best-fitting allele over ALL scored traversals is not a called one. A read
    /// that fits neither called allele should not be asserted onto one, so this is off by default:
    /// for assembly anchors purity beats yield.
    bool keep_off_call = false;
};

/**
 * Accumulates anchors during the render pass and writes the TSV.
 *
 * Rows are ordered by node ID, with each anchor's read rows immediately following it. Read rows do
 * not repeat the anchor key -- the single biggest saving available, since the snarl ID is around 24
 * characters and would otherwise repeat per read -- so a read row means nothing without the anchor
 * row above it. The header says so. Nothing needs to sort the file, because it is written sorted.
 */
class AnchorWriter {
public:
    struct ReadRow {
        string name;
        int64_t offset = -1;
        float score = 0.0f;
        uint8_t strand = 0;
    };
    struct Anchor {
        nid_t node = 0;
        string snarl;
        /// Which haplotype of the settled genotype this anchor partitions to, and which candidate
        /// traversal that haplotype is.
        ///
        /// `slot` is an index into the settled *phased* pair, so slot i is field i of the record's
        /// `GT` at the same snarl id -- slot 0 the left allele, slot 1 the right. A homozygote
        /// collapses to one slot holding every read, and then there is no haplotype information to
        /// join on.
        ///
        /// A nested HAPLOID site collapses to one slot too, but for the opposite reason, and its
        /// slot is *not* always 0: the chain sits on one strand of a diploid locus because the
        /// parent's other allele deletes it, so its single slot IS a haplotype -- the one
        /// `nested_strand` names, which the VCF writes as `a|.` or `.|a`. Slot 0 for a `.|a` site
        /// would name the wrong haplotype, and nothing in the VCF could see it. That was v4's
        /// remaining half of the bug the v3 -> v4 bump was made for.
        ///
        /// `allele` is the index into the site's candidate traversal set, which is NOT the
        /// VCF's ALT numbering: the ALT list is chosen later, in emit_variant, and anchors are built
        /// before it. Two slots carrying the same `allele` is what a homozygote looks like before
        /// collapsing; it cannot otherwise happen.
        int slot = 0;
        int allele = -1;
        double gqn = -1.0;
        double explained = 1.0;
        /// The site's mean per-read `score`, over the reads it actually emitted, each read once.
        ///
        /// Site-level like `gqn` and `explained`, so it repeats across the site's anchors. It is
        /// the same quantity `--phase-min-q` thresholds -- `PhaseSite::reliability` in
        /// read_phasing.hpp -- and it is low exactly where the reads cannot tell the site's alleles
        /// apart, which on ONT means a 1 bp indel. Derivable from the R rows, and written anyway
        /// for two reasons: deriving it costs a pass over every read row plus knowing to count a
        /// read once across both pins and both slots -- and the R rows are written to one decimal,
        /// so a derived value is only accurate to about 0.05 where this one is exact.
        ///
        /// -1 where the site emitted no read, which cannot happen for an anchor that is written.
        double reliability = -1.0;
        vector<ReadRow> reads;
    };

    explicit AnchorWriter(size_t threads);

    /// Add an anchor. Indexes its queue by `omp_get_thread_num()` rather than by any index the
    /// caller supplies: the render pass is an OpenMP loop, and a caller-supplied index that did not
    /// happen to equal the thread number could put two threads on one queue.
    void add(Anchor&& anchor);

    size_t anchor_count() const;
    size_t read_row_count() const;

    /// Sort by (node, snarl, slot) and write. Returns false if the file could not be opened.
    bool write(const string& path, const string& graph_name, const string& sample,
               const string& reads_source, double mismap_min, const AnchorParams& params);

private:
    vector<vector<Anchor>> queues;
};

/**
 * Turn one site's retained evidence and its settled genotype into anchors.
 *
 * The read is assigned to the called allele with the larger responsibility
 *
 *     resp_i = (1 - e_r) * w_i * rel(r, a_i),   resp_x = e_r
 *
 * where `resp_x` is "this read is from somewhere else", and the score is the phred of the
 * complement of the winner's share. It is bounded above by the mismap floor -- at `--mismap-min
 * 0.02` a perfectly discriminating read on a balanced het scores about 14, not 60 -- which is honest
 * and makes that flag directly meaningful here: read it as "P(this read's evidence at this site is
 * unreliable)".
 */
/// The mixture weights the per-read responsibilities use, for a settled slot-to-allele mapping.
///
/// Shared with read-backed phasing rather than reimplemented there: the weights decide how a read is
/// split between the two haplotypes, so two copies that drifted apart would make the anchor file and
/// the phase disagree about the same read.
vector<double> site_slot_weights(const vector<uint32_t>& allele_length, size_t n_alleles,
                                 float mean_read_length, bool length_weighted,
                                 const vector<int>& slot_allele);

/// `haploid_slot` is which strand a one-allele `genotype` sits on, 0 or 1, and is ignored for any
/// other genotype size. Supplied by the caller because `nested_strand` lives on the phasing record,
/// which this layer has no access to -- the same reason the pair arrives already phase-ordered.
void build_site_anchors(const AnchorSiteEvidence& evidence, const vector<int>& genotype,
                        const string& snarl_id, double gqn, double explained, int haploid_slot,
                        const AnchorParams& params, AnchorCounters& counters,
                        vector<AnchorWriter::Anchor>& out);

}

#endif
