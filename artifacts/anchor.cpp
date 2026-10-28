#include "anchor.hpp"
#include "read_phasing.hpp"
#include "read_strand_table.hpp"

#include "../src/version.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <functional>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

#include "parallel_sort.hpp"
#include "../src/path.hpp"
#include "../src/utility.hpp"

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

void AnchorPinCounters::report(ostream& out) const {
    out << "[vg call] anchors: " << verified.load() << " pins verified against the graph, "
        << verify_failed.load() << " failed" << endl;
    out << "[vg call] anchor pins refused: " << no_visit.load()
        << " read does not visit the boundary node, " << repeat_visit.load()
        << " visits it more than once, " << unaligned_base.load()
        << " no read base aligned to the pin's node base, " << no_neighbour.load()
        << " nothing upstream of an entry pin" << endl;
}

void AnchorPhaseCheckCounters::report(ostream& out) const {
    if (checked.load() == 0) {
        return;
    }
    const size_t n = checked.load();
    const size_t ok = agree.load();
    const size_t cn = confident.load();
    const size_t ck = confident_agree.load();
    out << "[vg call] anchors: cross-site phase reproduces the allele partition at het "
           "sites (held out, the split's own accuracy -- reported only when het placement "
           "does NOT use the strand, or it would measure itself): "
        << ok << "/" << n << " agree (" << (100.0 * (double)ok / (double)n) << "%)";
    if (cn > 0) {
        out << ", confident " << ck << "/" << cn << " ("
            << (100.0 * (double)ck / (double)cn) << "%)";
    }
    out << ", " << no_opinion.load() << " reads with no cross-site opinion" << endl;
}

void AnchorPhaseCheckCounters::count(double lo, int slot, double phase_min) {
    if (std::isnan(lo) || lo == 0.0) {
        no_opinion.fetch_add(1);
        return;
    }
    const bool agrees = (lo > 0.0 ? 0 : 1) == slot;
    checked.fetch_add(1);
    if (agrees) {
        agree.fetch_add(1);
    }
    if (std::abs(lo) >= phase_min) {
        confident.fetch_add(1);
        if (agrees) {
            confident_agree.fetch_add(1);
        }
    }
}

void AnchorCounters::report(ostream& out) const {
    pins.report(out);
    out << "[vg call] anchor reads: " << coincident.load()
        << " with both pins at one position (kept at S, dropped at E), " << off_call.load()
        << " fitting no called allele";
    if (shared_name.load() > 0) {
        out << ", " << shared_name.load()
            << " sharing a read name with another read at the same site (paired mates do; the name "
               "then does not identify an alignment)";
    }
    if (degenerate_site.load() > 0) {
        out << ", " << degenerate_site.load() << " sites whose two pins are one point";
    }
    if (single_pin.load() > 0) {
        out << ", " << single_pin.load()
            << " sites given only their start pin because the end pin held no reads of its own";
    }
    out << endl;
    if (hom_split.load() > 0 || hom_unsplit.load() > 0) {
        out << "[vg call] anchors: " << hom_split.load() << " homozygous sites split by read phase, "
            << hom_unsplit.load()
            << " left collapsed for want of a confident partition on both strands";
        if (hom_split_no_opinion.load() > 0) {
            out << "; " << hom_split_no_opinion.load()
                << " read placements at split sites had no cross-site opinion; "
                << hom_split_coin.load() << " were assigned by the per-read coin and "
                << (hom_split_no_opinion.load() - hom_split_coin.load())
                << " dropped for having a strand only in another phase set";
        }
        out << endl;
    }
    if (het_phase_tilted.load() > 0) {
        out << "[vg call] anchors: " << het_phase_tilted.load()
            << " read placements at heterozygous sites used the read's cross-site strand";
        if (het_strict_moved.load() > 0) {
            out << " (--anchors-strict-hets, which moved " << het_strict_moved.load()
                << " of them off the slot the allele match alone would have chosen)";
        }
        out << endl;
    }
    phase_check.report(out);
}

////////////////////////////////////////////////////////////////////////////////
// Pin resolution
////////////////////////////////////////////////////////////////////////////////

/**
 * The first and last read base of one mapping that consume a node base.
 *
 * Kept for each mapping rather than as one run over the read, because the entry pin must step
 * exactly one node back along the read's walk. A merged run cannot tell the last base of the
 * previous node from the last base before several deleted nodes.
 */
struct MappingExtent {
    /// The first read base aligned to a node base, or -1 when the mapping consumes no node base,
    /// as for a node the read deletes.
    int64_t first_consuming = -1;
    /// The last read base aligned to a node base, or -1 when the mapping consumes no node base.
    int64_t last_consuming = -1;
};

/// Return the extent of a mapping whose first read base is `read_start`.
static MappingExtent extent_of(const Mapping& m, size_t read_start) {
    MappingExtent out;
    size_t read_pos = read_start;
    for (int j = 0; j < m.edit_size(); ++j) {
        const Edit& e = m.edit(j);
        if (e.from_length() > 0 && e.to_length() > 0 && e.from_length() == e.to_length()) {
            if (out.first_consuming < 0) {
                out.first_consuming = (int64_t)read_pos;
            }
            out.last_consuming = (int64_t)(read_pos + (size_t)e.to_length() - 1);
        }
        read_pos += (size_t)e.to_length();
    }
    return out;
}

/// Return a hash of where an alignment starts, which tells apart alignments of one read name.
static uint32_t alignment_tag(const Path& path) {
    if (path.mapping_size() == 0) {
        return 0;
    }
    const Position& first = path.mapping(0).position();
    return (uint32_t)(((uint64_t)first.node_id() * 0x9E3779B97F4A7C15ULL) >> 32)
           ^ (uint32_t)(first.offset() * 2 + (first.is_reverse() ? 1 : 0));
}

/// A read's visits to the node a pin lies on.
struct PinNodeVisit {
    /// How many mappings of the read visit the node.
    size_t visits = 0;
    /// The index of the last mapping that visits the node, or -1 for none.
    int mapping = -1;
    /// The read base where that mapping starts.
    size_t read_start = 0;
    /// The read base where the mapping before it starts, or 0 when it is the first mapping.
    size_t before_read_start = 0;
};

/**
 * Find a read's visits to a node from the read's index of mappings near the site.
 *
 * The pin node is a boundary of the site, so every visit to it is in the index.
 */
static PinNodeVisit find_indexed_visit(const SiteRead& read, nid_t node_id) {
    const Path& path = read.aln->path();
    PinNodeVisit found;
    for (size_t k = 0; k < read.mapping_count; ++k) {
        const auto i = (int)read.mappings[k];
        if (path.mapping(i).position().node_id() != node_id) {
            continue;
        }
        ++found.visits;
        found.mapping = i;
        found.read_start = read.read_offsets[i];
    }
    if (found.mapping > 0) {
        found.before_read_start = read.read_offsets[found.mapping - 1];
    }
    return found;
}

/// Find a read's visits to a node by walking its whole alignment.
static PinNodeVisit find_walked_visit(const Path& path, nid_t node_id) {
    PinNodeVisit found;
    size_t read_pos = 0;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& m = path.mapping(i);
        if (m.position().node_id() == node_id) {
            ++found.visits;
            found.mapping = i;
            found.read_start = read_pos;
        }
        read_pos += (size_t)mapping_to_length(m);
    }
    if (found.mapping > 0) {
        found.before_read_start =
            found.read_start - (size_t)mapping_to_length(path.mapping(found.mapping - 1));
    }
    return found;
}

/**
 * Return the node base a pin is defined against, in node-forward coordinates.
 *
 * An exit (S) pin follows the last base of the node in the site's direction. An entry (E) pin
 * precedes the first base of the node in the site's direction. A pin inside the node uses the
 * node base before it in the site's direction for an exit pin, and the node base after it for an
 * entry pin.
 *
 * @returns -1 when a pin inside the node has no such base.
 */
static int64_t pin_reference_base(size_t node_len, bool site_backward, bool exit_pin, int64_t pin) {
    if (pin >= 0) {
        const int64_t inside = (exit_pin != site_backward) ? pin - 1 : pin;
        if (inside < 0 || inside >= (int64_t)node_len) {
            return -1;
        }
        return inside;
    }
    // An exit pin with the site reading the node forward, or an entry pin with the site reading it
    // backward, is defined against the last base of the node's forward strand.
    if (exit_pin != site_backward) {
        return (int64_t)node_len - 1;
    }
    return 0;
}

/// The read base aligned to one node base, as found in the mapping's edits.
struct AlignedReadBase {
    /// The read base, or -1 when no read base is aligned to the node base.
    int64_t index = -1;
    /// Whether the edit holding the base is a match.
    bool exact = false;
    /// The base the read has in a mismatch edit, or 0 when the edit records none.
    char substituted = 0;
};

/// Return the read base aligned to the node base `offset` bases into an edit of equal from and to
/// lengths, which starts at read base `read_start`.
static AlignedReadBase read_base_in_edit(const Edit& e, size_t offset, size_t read_start) {
    AlignedReadBase found;
    found.index = (int64_t)(read_start + offset);
    if (e.sequence().empty()) {
        found.exact = true;
    } else if (e.sequence().size() == (size_t)e.to_length()) {
        found.substituted = e.sequence()[offset];
    }
    return found;
}

/**
 * Find the read base aligned to node base `q_vis` of a mapping, by walking its edits.
 *
 * @param q_vis The node base, in the orientation the read visits the node in, in which the
 * mapping's offsets are measured.
 * @param read_start The read base where the mapping starts.
 */
static AlignedReadBase aligned_read_base(const Mapping& m, size_t q_vis, size_t read_start) {
    AlignedReadBase found;
    auto node_pos = (size_t)m.position().offset();
    size_t rp = read_start;
    for (int j = 0; j < m.edit_size(); ++j) {
        const Edit& e = m.edit(j);
        const auto from = (size_t)e.from_length();
        const auto to = (size_t)e.to_length();
        if (from > 0 && q_vis >= node_pos && q_vis < node_pos + from) {
            // A deleted node base (to == 0), or one inside a replacement of unequal length, has no
            // corresponding read base.
            if (from == to) {
                found = read_base_in_edit(e, q_vis - node_pos, rp);
            }
            break;
        }
        node_pos += from;
        rp += to;
    }
    return found;
}

/**
 * Check a read base against the graph, and count the outcome.
 *
 * On a match edit the read base must be the graph's base, complemented when the read visits the
 * node in reverse, which catches offset and strand errors. On a mismatch edit the read base
 * differs from the graph's, so it is compared with the edit's recorded sequence, which still
 * catches offset errors. A mismatch edit with no recorded sequence cannot be checked.
 *
 * @returns false if the check failed.
 */
static bool verify_read_base(char read_base, const AlignedReadBase& aligned, const HandleGraph& graph,
                             handle_t handle, size_t q_fwd, bool is_rev, AnchorPinCounters& counters) {
    char want = 0;
    if (aligned.exact) {
        const char node_base = graph.get_base(handle, q_fwd);
        want = (char)toupper(is_rev ? reverse_complement(node_base) : node_base);
    } else if (aligned.substituted != 0) {
        want = (char)toupper(aligned.substituted);
    } else {
        return true;
    }
    if (read_base != want) {
        ++counters.verify_failed;
        return false;
    }
    ++counters.verified;
    return true;
}

/**
 * Return the read base before an entry pin at a node end, or -1 if the read has none.
 *
 * The reference base of an entry pin is downstream of the pin, so the offset is the read base
 * just upstream. That is the adjacent base on the previous node of the read's walk. Stepping back
 * by read index instead would pass over a node the read deletes, which has no read bases, and land
 * on the pin of the neighbouring snarl. So we step one mapping along the walk, skipping insertions
 * and soft clips, which consume no node base. If that mapping consumes no node base at all, the
 * node upstream is deleted in this read and the read has no placement.
 */
static int64_t entry_pin_neighbour(const Alignment& aln, const PinNodeVisit& visit,
                                   uint8_t direction) {
    const Path& path = aln.path();
    const int neighbour = (direction == 0) ? visit.mapping - 1 : visit.mapping + 1;
    if (neighbour < 0 || neighbour >= path.mapping_size()) {
        // The read begins (or ends) here, so it does not cross the pin.
        return -1;
    }
    size_t neighbour_read_start = visit.before_read_start;
    if (direction != 0) {
        neighbour_read_start = visit.read_start + (size_t)mapping_to_length(path.mapping(visit.mapping));
    }
    const MappingExtent adjacent = extent_of(path.mapping(neighbour), neighbour_read_start);
    const int64_t candidate = (direction == 0) ? adjacent.last_consuming : adjacent.first_consuming;
    if (candidate < 0 || candidate >= (int64_t)aln.sequence().size()) {
        return -1;
    }
    return candidate;
}

AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters, int64_t pin) {
    AnchorPlacement out;
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();
    out.alignment = alignment_tag(path);

    const PinNodeVisit visit = read.indexed() ? find_indexed_visit(read, node_id)
                                              : find_walked_visit(path, node_id);
    if (visit.visits == 0) {
        // A read that starts inside the site does not visit the start boundary, and one that ends
        // inside it does not visit the end boundary.
        ++counters.pins.no_visit;
        return out;
    }
    if (visit.visits > 1) {
        ++counters.pins.repeat_visit;
        return out;
    }

    const Mapping& m = path.mapping(visit.mapping);
    const bool is_rev = m.position().is_reverse();
    handle_t handle = graph.get_handle(node_id, false);
    const size_t node_len = graph.get_length(handle);
    if (node_len == 0) {
        ++counters.pins.unaligned_base;
        return out;
    }
    const int64_t reference_base = pin_reference_base(node_len, site_backward, exit_pin, pin);
    if (reference_base < 0) {
        ++counters.pins.unaligned_base;
        return out;
    }
    const auto q_fwd = (size_t)reference_base;
    const size_t q_vis = is_rev ? (node_len - 1 - q_fwd) : q_fwd;

    const AlignedReadBase aligned = aligned_read_base(m, q_vis, visit.read_start);
    if (aligned.index < 0 || (size_t)aligned.index >= aln.sequence().size()) {
        ++counters.pins.unaligned_base;
        return out;
    }
    if (const auto read_base = (char)toupper(aln.sequence()[aligned.index]);
        !verify_read_base(read_base, aligned, graph, handle, q_fwd, is_rev, counters.pins)) {
        return out;
    }

    // `is_rev` is relative to the node, and `site_backward` is the site's orientation on it, so
    // the read crosses the pin in the site's direction when they agree.
    out.direction = (is_rev == site_backward) ? 0 : 1;

    if (exit_pin) {
        // The reference base is just upstream of the pin, which is what `offset` reports.
        out.offset = aligned.index;
        return out;
    }
    if (pin >= 0) {
        // An entry pin inside the node takes the read base just upstream of the reference base. A
        // read whose alignment begins (or ends) at the reference base, soft clipped or not, does
        // not cross the pin.
        if (const int64_t upstream = out.direction == 0 ? aligned.index - 1 : aligned.index + 1;
            upstream < softclip_start(path)
            || upstream >= (int64_t)aln.sequence().size() - softclip_end(path)) {
            ++counters.pins.no_neighbour;
        } else {
            out.offset = upstream;
        }
        return out;
    }
    const int64_t neighbour = entry_pin_neighbour(aln, visit, out.direction);
    if (neighbour < 0) {
        ++counters.pins.no_neighbour;
        return out;
    }
    out.offset = neighbour;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Read names
////////////////////////////////////////////////////////////////////////////////

uint32_t ReadNameTable::intern(string_view name) {
    size_t s = std::hash<string_view>{}(name) % SHARDS;
    Shard& shard = shards[s];
    std::scoped_lock guard(shard.mutex);
    if (auto found = shard.index.find(name); found != shard.index.end()) {
        return found->second;
    }
    size_t position = shard.count;
    if (position / CHUNK >= MAX_CHUNKS || position * SHARDS + s > UINT32_MAX) {
        throw runtime_error("too many distinct read names for a 32-bit index");
    }
    string* chunk = shard.chunks[position / CHUNK].load(std::memory_order_acquire);
    if (chunk == nullptr) {
        chunk = new string[CHUNK];
        shard.chunks[position / CHUNK].store(chunk, std::memory_order_release);
    }
    chunk[position % CHUNK] = name;
    ++shard.count;
    auto index = (uint32_t)(position * SHARDS + s);
    shard.index.try_emplace(string_view(chunk[position % CHUNK]), index);
    return index;
}

string_view ReadNameTable::name(uint32_t index) const {
    const Shard& shard = shards[index % SHARDS];
    size_t position = index / SHARDS;
    const string* chunk = shard.chunks[position / CHUNK].load(std::memory_order_acquire);
    return chunk[position % CHUNK];
}

ReadNameTable::Shard::~Shard() {
    for (const auto& chunk : chunks) {
        delete[] chunk.load();
    }
}

ReadNameTable& read_names() {
    // Leak the table, so that the process does not free millions of names one by one at exit.
    static auto* table = new ReadNameTable();
    return *table;
}

size_t AnchorSiteEvidence::bytes() const {
    // Names are held once, in read_names(), not per site.
    size_t total = sizeof(AnchorSiteEvidence);
    total += reads.capacity() * sizeof(AnchorRead);
    total += rel.capacity() * sizeof(float);
    total += allele_length.capacity() * sizeof(uint32_t);
    return total;
}

////////////////////////////////////////////////////////////////////////////////
// Building anchors from a chosen genotype
////////////////////////////////////////////////////////////////////////////////

/// Return whether every allele of a genotype is one of the site's candidate traversals.
///
/// A star or missing allele has no traversal to assign reads to.
static bool genotype_has_traversals(const vector<int>& genotype, size_t n_alleles) {
    return std::ranges::all_of(genotype, [n_alleles](int a) {
        return a >= 0 && (size_t)a < n_alleles;
    });
}

/// Return whether every allele of a genotype is the same one, as in a homozygote or a haploid.
static bool all_alleles_equal(const vector<int>& genotype) {
    return std::ranges::all_of(genotype, [&genotype](int a) { return a == genotype[0]; });
}

/**
 * Return whether a diploid homozygous site's reads divide confidently between the two strands.
 *
 * A read counts when it crosses a pin and its strand log-odds are at least `params.phase_min` in
 * size. Each strand needs at least `params.phase_min_side` such reads.
 */
static bool reads_divide_by_strand(const AnchorSiteEvidence& evidence,
                                   const vector<double>& read_strand, const AnchorParams& params) {
    size_t side0 = 0;
    size_t side1 = 0;
    for (size_t r = 0; r < evidence.reads.size(); ++r) {
        if (!evidence.reads[r].start_pin.placed() && !evidence.reads[r].end_pin.placed()) {
            continue;
        }
        const double lo = read_strand[r];
        // NaN is a strand from another phase set, and 0 is no strand log-odds at all. Neither is
        // evidence for strand 1, though the `lo > 0.0` test below would count a 0 there.
        if (std::isnan(lo) || lo == 0.0 || std::abs(lo) < params.phase_min) {
            continue;
        }
        if (lo > 0.0) {
            ++side0;
        } else {
            ++side1;
        }
    }
    return side0 >= params.phase_min_side && side1 >= params.phase_min_side;
}

/**
 * Return whether to split a homozygous site into two slots, one for each strand, and count the
 * decision.
 *
 * Only a diploid site with the strand log-odds of all its reads is considered, and only under
 * `params.hom_split`.
 */
static bool split_homozygote(const AnchorSiteEvidence& evidence, const vector<int>& genotype,
                             const vector<double>* read_strand, const AnchorParams& params,
                             AnchorCounters& counters) {
    if (!params.hom_split || genotype.size() != 2 || read_strand == nullptr
        || read_strand->size() != evidence.reads.size()) {
        return false;
    }
    if (reads_divide_by_strand(evidence, *read_strand, params)) {
        ++counters.hom_split;
        return true;
    }
    ++counters.hom_unsplit;
    return false;
}

/// Count the reads of a site that share a name with another of its reads, which the file cannot
/// tell apart.
static void count_shared_names(const vector<AnchorRead>& reads, AnchorCounters& counters) {
    unordered_map<uint32_t, size_t> name_count;
    name_count.reserve(reads.size() * 2);
    for (const AnchorRead& read : reads) {
        ++name_count[read.read];
    }
    for (const auto& [name, count] : name_count) {
        if (count > 1) {
            counters.shared_name += count;
        }
    }
}

/// The two anchors of one slot of a site, filled as reads are placed.
struct SlotAnchors {
    /// The anchor at the start pin.
    AnchorWriter::Anchor start;
    /// The anchor at the end pin.
    AnchorWriter::Anchor end;
    /// The number of reads placed in the slot that cross the end pin but not the start pin.
    ///
    /// Counted per slot rather than per site, since anchors are written and filtered by
    /// `min_reads` per slot. See AnchorParams::end_pin_min_new.
    size_t new_at_end = 0;
};

/// The slot a read joins, with the responsibilities its confidence comes from.
struct SlotChoice {
    /// The slot.
    size_t slot = 0;
    /// The untilted responsibility of the slot's allele for the read.
    double resp = 0.0;
    /// The sum of the read's untilted responsibilities, the mismapping outcome included.
    double total = 0.0;
    /// The read's phred-scaled confidence in the slot, from its share `resp / total`.
    double score = 0.0;
};

/**
 * Choose the slot of a read at a homozygous site split by strand, or return nothing to drop the
 * read.
 *
 * Both slots spell the same allele, so the responsibilities cannot choose between them, and the
 * read's strand log-odds `lo` decide instead. The read's share is computed over the site's one
 * distinct allele rather than over the two slots. So a read's confidence, and with it
 * `min_read_score` and the reliability column, mean the same at split sites as elsewhere.
 */
static std::optional<SlotChoice> choose_split_slot(const AnchorSiteEvidence& evidence, size_t r,
                                                   int allele, double lo, AnchorCounters& counters) {
    const AnchorRead& read = evidence.reads[r];
    SlotChoice choice;
    if (std::isnan(lo)) {
        // The read has a strand only in another phase set, or in more than one, and phase sets
        // number their strands independently. Where the read also crosses a heterozygous site of
        // this phase set, a coin would contradict the slot that site's allele gives it half the
        // time. Elsewhere a coin would join two phase sets that the phasing left apart.
        ++counters.hom_split_no_opinion;
        return std::nullopt;
    }
    if (lo == 0.0) {
        // The read reached no other phased site, or none whose alleles it can tell apart. Either
        // slot spells its allele, so place it by a hash of its name. The read then gets the same
        // slot at every site and both pins, so it stays on one strand.
        ++counters.hom_split_no_opinion;
        ++counters.hom_split_coin;
        choice.slot = (size_t)(std::hash<string_view>{}(read_names().name(read.read)) & 1ULL);
    } else {
        choice.slot = lo > 0.0 ? 0 : 1;
    }
    const auto mismap = (double)read.mismap;
    choice.resp = (1.0 - mismap) * (double)evidence.rel_at(r, (size_t)allele);
    choice.total = mismap + choice.resp;
    return choice;
}

/**
 * Choose the slot of a read from its allele match, tilted by its strand log-odds as `params` say.
 *
 * Under `params.phase_hets` the disfavoured slot's weight is multiplied by exp(-|strand|). This
 * gives the weights `phase_aware_correction` uses in re-genotyping, without overflow. Under
 * `params.strict_hets` the sign of `strand` alone chooses. Whichever rule picks the slot, the
 * read's confidence comes from the untilted responsibilities. The confidence feeds the
 * reliability, which measures whether the site's own reads tell its alleles apart, so the strand
 * log-odds must not enter it.
 *
 * @param strand The read's strand log-odds, or 0 where they cannot be used here.
 */
static SlotChoice choose_allele_slot(const AnchorSiteEvidence& evidence, size_t r,
                                     const vector<int>& slot_allele, const vector<double>& weight,
                                     double strand, const AnchorParams& params,
                                     AnchorCounters& counters) {
    const auto mismap = (double)evidence.reads[r].mismap;
    const bool has_opinion = strand != 0.0;
    // The slot whose weight is tilted down, or none, and the factor applied to it.
    size_t disfavoured = slot_allele.size();
    double tilt = 1.0;
    if (params.phase_hets && !params.strict_hets && has_opinion) {
        tilt = std::exp(-std::abs(strand));
        disfavoured = strand > 0.0 ? 1 : 0;
        ++counters.het_phase_tilted;
    }
    SlotChoice choice;
    choice.total = mismap;
    double best_resp = -1.0;
    for (size_t i = 0; i < slot_allele.size(); ++i) {
        const auto rel = (double)evidence.rel_at(r, (size_t)slot_allele[i]);
        const double plain = (1.0 - mismap) * weight[i] * rel;
        const double tilted_weight = i == disfavoured ? weight[i] * tilt : weight[i];
        const double resp = (1.0 - mismap) * tilted_weight * rel;
        choice.total += plain;
        // Break ties on the allele, not on the slot. `slot_allele` is in phase order, which read
        // phasing can change, and the division of reads must not depend on the phase. A tie is
        // the same `rel` value for both alleles, so it is exact.
        if (resp > best_resp
            || (resp == best_resp && slot_allele[i] < slot_allele[choice.slot])) {
            best_resp = resp;
            choice.resp = plain;
            choice.slot = i;
        }
    }
    if (params.strict_hets && has_opinion) {
        const size_t want = strand > 0.0 ? 0 : 1;
        if (want != choice.slot) {
            ++counters.het_strict_moved;
        }
        choice.slot = want;
        choice.resp = (1.0 - mismap) * weight[want]
                      * (double)evidence.rel_at(r, (size_t)slot_allele[want]);
        ++counters.het_phase_tilted;
    }
    return choice;
}

/// Return whether a read fits a called allele at least as well as it fits every candidate
/// allele.
static bool called_allele_fits_best(const AnchorSiteEvidence& evidence, size_t r,
                                    const vector<int>& slot_allele) {
    double row_best = 0.0;
    for (size_t a = 0; a < evidence.n_alleles; ++a) {
        row_best = max(row_best, (double)evidence.rel_at(r, a));
    }
    bool called_is_best = false;
    for (int allele : slot_allele) {
        if ((double)evidence.rel_at(r, (size_t)allele) >= row_best) {
            called_is_best = true;
        }
    }
    return called_is_best;
}

/// Return the phred-scaled confidence of a read whose slot has `share` of its responsibility, or
/// 99 when the share is whole.
static double share_confidence(double share) {
    double score = 99.0;
    if (share < 1.0) {
        score = -10.0 * log10(1.0 - share);
    }
    return score;
}

/**
 * A judge that decides which slot of a site each read joins, and how confident the read is there.
 *
 * At a homozygous site split by strand, the strand log-odds of the read decide. Elsewhere the
 * allele match decides, tilted by the strand log-odds as the parameters say. The judge turns away
 * a read that fits no called allele, or that fits an uncalled allele better unless
 * `AnchorParams::keep_off_call`, or whose confidence is below `AnchorParams::min_read_score`.
 */
class SlotJudge {
public:
    /**
     * Make a judge for one site.
     *
     * @param slot_allele The allele each slot holds.
     * @param weight The allele-length weight of each slot.
     * @param read_strand Each read's strand log-odds, or null. Must be present when
     * `phase_split` is true.
     * @param phase_split Whether the site is a homozygote split by strand.
     */
    SlotJudge(const AnchorSiteEvidence& evidence, const vector<int>& slot_allele,
              const vector<double>& weight, const vector<double>* read_strand, bool phase_split,
              const AnchorParams& params, AnchorCounters& counters)
        : evidence(evidence), slot_allele(slot_allele), weight(weight), read_strand(read_strand),
          phase_split(phase_split),
          tiltable(read_strand != nullptr && slot_allele.size() == 2
                   && read_strand->size() == evidence.reads.size()
                   && slot_allele[0] != slot_allele[1]),
          params(params), counters(counters) {}

    /// Return the slot read `r` joins and its confidence, or nothing if the read is turned away.
    std::optional<SlotChoice> judge(size_t r) const;

private:
    const AnchorSiteEvidence& evidence;
    const vector<int>& slot_allele;
    const vector<double>& weight;
    const vector<double>* read_strand;
    bool phase_split;
    /// Whether the strand log-odds can tilt a read's choice between two slots with different
    /// alleles.
    bool tiltable;
    const AnchorParams& params;
    /// Counts the reads turned away and how slots were chosen.
    AnchorCounters& counters;
};

/// Return a read's strand log-odds if they name a strand, or 0 if they are 0 or NaN.
static double strand_opinion(double lo) {
    return std::isnan(lo) ? 0.0 : lo;
}

std::optional<SlotChoice> SlotJudge::judge(size_t r) const {
    std::optional<SlotChoice> choice;
    if (phase_split) {
        choice = choose_split_slot(evidence, r, slot_allele[0], (*read_strand)[r], counters);
        if (!choice) {
            return std::nullopt;
        }
    } else {
        const double strand = tiltable ? strand_opinion((*read_strand)[r]) : 0.0;
        choice = choose_allele_slot(evidence, r, slot_allele, weight, strand, params, counters);
    }
    if (!(choice->resp > 0.0) || !(choice->total > 0.0)) {
        // The read fits no called allele at all.
        ++counters.off_call;
        return std::nullopt;
    }
    if (!called_allele_fits_best(evidence, r, slot_allele)) {
        ++counters.off_call;
        if (!params.keep_off_call) {
            return std::nullopt;
        }
    }
    choice->score = share_confidence(choice->resp / choice->total);
    if (choice->score < params.min_read_score) {
        return std::nullopt;
    }
    return choice;
}

/**
 * Add a read to the anchors of its slot, at each pin it crosses.
 *
 * @param degenerate Whether the site's two pins are one point, so that only the start anchor is
 * used.
 */
static void add_read_rows(const AnchorRead& read, float score, bool degenerate, SlotAnchors& slot,
                          AnchorCounters& counters) {
    AnchorWriter::ReadRow row;
    row.read = read.read;
    row.alignment = read.start_pin.alignment;
    row.score = score;
    if (read.start_pin.placed()) {
        row.offset = read.start_pin.offset;
        row.direction = read.start_pin.direction;
        slot.start.reads.push_back(row);
    }
    if (degenerate || !read.end_pin.placed()) {
        return;
    }
    // A read that crosses straight from the start boundary to the end boundary carries an allele
    // deleting the site's interior, so both its pins are at one read position. Keep it at the start
    // pin only, so that a read position is in at most one anchor.
    if (read.start_pin.placed() && read.end_pin.offset == read.start_pin.offset
        && read.end_pin.direction == read.start_pin.direction) {
        ++counters.coincident;
        return;
    }
    row.offset = read.end_pin.offset;
    row.direction = read.end_pin.direction;
    slot.end.reads.push_back(row);
}

/// Append the anchors of each slot that hold enough reads to `out`.
///
/// A slot's end anchor is left out where the site is degenerate, or where it has too few reads of
/// its own; see AnchorParams::end_pin_min_new.
static void emit_slot_anchors(vector<SlotAnchors>& slots, bool degenerate,
                              const AnchorParams& params, AnchorCounters& counters,
                              vector<AnchorWriter::Anchor>& out) {
    for (SlotAnchors& slot : slots) {
        if (slot.start.reads.size() >= params.min_reads) {
            out.push_back(std::move(slot.start));
        }
        bool drop_end = degenerate;
        if (!drop_end && params.end_pin_min_new > 0 && slot.new_at_end < params.end_pin_min_new) {
            drop_end = true;
            ++counters.single_pin;
        }
        if (!drop_end && slot.end.reads.size() >= params.min_reads) {
            out.push_back(std::move(slot.end));
        }
    }
}

/**
 * Set the reliability of a site's anchors, `out[begin]` onward, to the mean confidence of the
 * reads they hold.
 *
 * Reads are keyed by name, which counts a read once across both pins and slots. Paired mates also
 * count once, at the higher of their scores, as `merge_mates` keeps the more confident mate.
 */
static void set_site_reliability(vector<AnchorWriter::Anchor>& out, size_t begin) {
    // Key by the name itself, not by its index in read_names(). The scores are summed in the map's
    // iteration order, which depends on the hashes of the keys, and the last digits of a
    // floating-point sum depend on that order. Keying by name keeps the file as it was.
    unordered_map<string_view, float> per_read;
    for (size_t i = begin; i < out.size(); ++i) {
        for (const AnchorWriter::ReadRow& row : out[i].reads) {
            auto [entry, added] = per_read.try_emplace(read_names().name(row.read), row.score);
            if (!added && row.score > entry->second) {
                entry->second = row.score;
            }
        }
    }
    double sum = 0.0;
    for (const auto& [name, score] : per_read) {
        sum += (double)score;
    }
    const double site_reliability = per_read.empty() ? -1.0 : sum / (double)per_read.size();
    for (size_t i = begin; i < out.size(); ++i) {
        out[i].reliability = site_reliability;
    }
}

/// Return whether a site's two pins are one point, as for a snarl whose two boundaries are the
/// same node, or an insertion with no repeat tract.
static bool pins_coincide(const AnchorSiteEvidence& evidence) {
    return evidence.start_node == evidence.end_node
           && (!evidence.inside_node || evidence.start_pin == evidence.end_pin);
}

/**
 * Make the empty anchors of each slot of a site.
 *
 * @param first_slot The slot number of the first slot.
 * @param site_anchor An anchor with the fields every anchor of the site shares.
 */
static vector<SlotAnchors> make_slot_anchors(const AnchorSiteEvidence& evidence,
                                             const vector<int>& slot_allele, int first_slot,
                                             const AnchorWriter::Anchor& site_anchor) {
    vector<SlotAnchors> slots(slot_allele.size());
    for (size_t i = 0; i < slots.size(); ++i) {
        AnchorWriter::Anchor& start = slots[i].start;
        start = site_anchor;
        start.node = evidence.start_node;
        start.pin = evidence.start_pin;
        start.slot = (int)i + first_slot;
        start.allele = slot_allele[i];
        slots[i].end = start;
        slots[i].end.node = evidence.end_node;
        slots[i].end.pin = evidence.end_pin;
    }
    return slots;
}

void build_site_anchors(const AnchorSiteEvidence& evidence, const AnchorCall& call,
                        string_view snarl_id, const AnchorParams& params,
                        vector<AnchorWriter::Anchor>& out,
                        const vector<double>* read_strand) {
    const vector<int>& genotype = call.genotype;
    const double gqn = call.gqn;
    AnchorCounters& counters = *params.counters;

    if (genotype.empty() || evidence.reads.empty() || evidence.n_alleles == 0
        || !genotype_has_traversals(genotype, evidence.n_alleles)) {
        return;
    }
    // A haploid genotype counts as homozygous here, so that `het_only` drops it too.
    const bool hom = all_alleles_equal(genotype);
    if (hom && params.het_only) {
        return;
    }
    // NaN means the site had no gap to normalise, which is different from zero, and is filtered
    // out. A negative GQN, on a record whose genotype the linkage model changed, is a real value
    // and is compared like any other.
    if (params.min_gqn > 0.0 && (std::isnan(gqn) || gqn < params.min_gqn)) {
        return;
    }

    // Make one slot per strand of the genotype, or one slot holding every read for a homozygote
    // that is not split.
    const bool phase_split = hom && split_homozygote(evidence, genotype, read_strand, params, counters);
    vector<int> slot_allele = genotype;
    if (hom) {
        slot_allele.assign((size_t)(phase_split ? 2 : 1), genotype[0]);
    }

    // A homozygote's single slot is slot 0 and names no strand. A nested haploid site's single
    // slot is the strand the phasing gave it, which is 1 for a `.|a` site. Only 0 or 1 is
    // written, since a slot indexes a GT field.
    const int base_slot = (genotype.size() == 1 && call.haploid_slot == 1) ? 1 : 0;

    const vector<double> weight = allele_length_weights(
        evidence.allele_length, evidence.n_alleles, evidence.mean_read_length,
        evidence.length_weighted, slot_allele);

    // Skip the end pin where the two pins are one point.
    const bool degenerate = pins_coincide(evidence);
    if (degenerate) {
        ++counters.degenerate_site;
    }
    AnchorWriter::Anchor site_anchor;
    site_anchor.snarl = snarl_id;
    site_anchor.gqn = gqn;
    site_anchor.explained = call.explained;
    site_anchor.edit = evidence.inside_node;
    vector<SlotAnchors> slots = make_slot_anchors(evidence, slot_allele, base_slot, site_anchor);

    count_shared_names(evidence.reads, counters);

    const SlotJudge judge(evidence, slot_allele, weight, read_strand, phase_split, params, counters);
    for (size_t r = 0; r < evidence.reads.size(); ++r) {
        const AnchorRead& read = evidence.reads[r];
        if (!read.start_pin.placed() && !read.end_pin.placed()) {
            continue;
        }
        const std::optional<SlotChoice> choice = judge.judge(r);
        if (!choice) {
            continue;
        }
        if (read.end_pin.placed() && !read.start_pin.placed()) {
            ++slots[choice->slot].new_at_end;
        }
        add_read_rows(read, (float)choice->score, degenerate, slots[choice->slot], counters);
    }

    const size_t out_begin = out.size();
    emit_slot_anchors(slots, degenerate, params, counters, out);
    if (out.size() > out_begin) {
        // Average only the reads a consumer can see, after the filters removed slots and end
        // anchors.
        set_site_reliability(out, out_begin);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Writer
////////////////////////////////////////////////////////////////////////////////

AnchorWriter::AnchorWriter(size_t threads) : queues(max((size_t)1, threads)) {
}

void AnchorWriter::add(Anchor&& anchor) {
    auto thread = (size_t)max(0, omp_get_thread_num());
    queues[thread % queues.size()].push_back(std::move(anchor));
}

size_t AnchorWriter::anchor_count() const {
    size_t total = 0;
    for (const auto& q : queues) {
        total += q.size();
    }
    return total;
}

/// Return -1, 0 or 1 as `a` is below, equal to or above `b`, with NaN below every number and
/// equal to itself.
static int compare_values(double a, double b) {
    const bool a_nan = std::isnan(a);
    const bool b_nan = std::isnan(b);
    if (a_nan && b_nan) {
        return 0;
    }
    if (a_nan) {
        return -1;
    }
    if (b_nan) {
        return 1;
    }
    if (a < b) {
        return -1;
    }
    if (b < a) {
        return 1;
    }
    return 0;
}

/**
 * Return whether anchor `a` comes before `b` on what the file writes for them besides node,
 * snarl and slot.
 *
 * That is the allele, gqn, explained, reliability, then the read rows in order (read name,
 * direction, offset, score). Two anchors equal on all of it are written identically, so their
 * order cannot show.
 */
static bool written_before(const AnchorWriter::Anchor& a, const AnchorWriter::Anchor& b,
                           const ReadNameTable& table) {
    if (a.allele != b.allele) {
        return a.allele < b.allele;
    }
    for (int c : {compare_values(a.gqn, b.gqn), compare_values(a.explained, b.explained),
                  compare_values(a.reliability, b.reliability)}) {
        if (c != 0) {
            return c < 0;
        }
    }
    if (a.reads.size() != b.reads.size()) {
        return a.reads.size() < b.reads.size();
    }
    for (size_t i = 0; i < a.reads.size(); ++i) {
        const AnchorWriter::ReadRow& x = a.reads[i];
        const AnchorWriter::ReadRow& y = b.reads[i];
        if (x.read != y.read && table.name(x.read) != table.name(y.read)) {
            return table.name(x.read) < table.name(y.read);
        }
        if (x.direction != y.direction) {
            return x.direction < y.direction;
        }
        if (x.offset != y.offset) {
            return x.offset < y.offset;
        }
        if (const int c = compare_values(x.score, y.score); c != 0) {
            return c < 0;
        }
    }
    return false;
}

vector<AnchorWriter::Anchor> AnchorWriter::take_sorted() {
    vector<Anchor> all;
    all.reserve(anchor_count());
    for (auto& q : queues) {
        std::ranges::move(q, back_inserter(all));
        q.clear();
    }
    // Sort the reads within each anchor by read name, so that the file does not depend on thread
    // scheduling. Do it before the anchors are sorted, so that their sort can compare whole
    // anchors.
    const ReadNameTable& table = read_names();
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < all.size(); ++i) {
        std::ranges::sort(all[i].reads, [&table](const ReadRow& x, const ReadRow& y) {
            if (x.read != y.read) {
                return table.name(x.read) < table.name(y.read);
            }
            return x.offset < y.offset;
        });
    }
    // Order the anchors by node ID, then by snarl and slot, since two snarls can share a boundary
    // node. Order anchors that agree on node, snarl and slot by the rest of what is written for
    // them, so that the order cannot depend on how the parallel sort shares out the work, nor on
    // which thread collected which anchor.
    parallel_sort(all.begin(), all.end(), [&table](const Anchor& a, const Anchor& b) {
        if (a.node != b.node) {
            return a.node < b.node;
        }
        if (a.snarl != b.snarl) {
            return a.snarl < b.snarl;
        }
        if (a.slot != b.slot) {
            return a.slot < b.slot;
        }
        if (a.pin != b.pin) {
            return a.pin < b.pin;
        }
        return written_before(a, b, table);
    });
    return all;
}

/**
 * A point of a read alignment: the alignment, as its read index and alignment hash, and the read
 * base the point precedes.
 */
using ReadPoint = pair<uint64_t, int64_t>;

/// Return the point of a read alignment at which a read row crosses its pin.
static ReadPoint point_of(const AnchorWriter::ReadRow& r) {
    return ReadPoint(((uint64_t)r.read << 32) | r.alignment, r.direction == 0 ? r.offset + 1 : r.offset);
}

/// A hash function for read points.
struct ReadPointHash {
    size_t operator()(const ReadPoint& p) const {
        return std::hash<uint64_t>{}(p.first ^ ((uint64_t)p.second * 0x9E3779B97F4A7C15ULL));
    }
};

/**
 * Remove read rows from off-panel anchors so that no point of a read alignment is in two
 * anchors, then remove the off-panel anchors left with too few reads.
 *
 * In file order, an off-panel anchor leaves out a read whose point a snarl's anchor holds, or an
 * earlier off-panel anchor kept. Snarl anchors are not checked against one another. Two snarls'
 * pins face into their own snarls, so no read crosses one point at both.
 *
 * @param tally Where the read rows and anchors left out are counted.
 */
static void leave_out_held_points(vector<AnchorWriter::Anchor>& all, size_t min_reads,
                                  AnchorWriter::Tally& tally) {
    // Who holds each point of an off-panel anchor: 0 for no anchor yet, 1 for a snarl's anchor,
    // and 2 for an off-panel anchor that was kept.
    unordered_map<ReadPoint, uint8_t, ReadPointHash> held(1024);
    for (const AnchorWriter::Anchor& a : all) {
        if (!a.edit) {
            continue;
        }
        for (const AnchorWriter::ReadRow& r : a.reads) {
            held.try_emplace(point_of(r), 0);
        }
    }
    for (const AnchorWriter::Anchor& a : all) {
        if (a.edit || held.empty()) {
            continue;
        }
        for (const AnchorWriter::ReadRow& r : a.reads) {
            if (auto found = held.find(point_of(r)); found != held.end()) {
                found->second = 1;
            }
        }
    }
    tally.placements_left_out = 0;
    for (AnchorWriter::Anchor& a : all) {
        if (!a.edit) {
            continue;
        }
        tally.placements_left_out +=
            std::erase_if(a.reads, [&held](const AnchorWriter::ReadRow& r) { return held[point_of(r)] != 0; });
        if (a.reads.size() >= min_reads) {
            for (const AnchorWriter::ReadRow& r : a.reads) {
                held[point_of(r)] = 2;
            }
        }
    }
    tally.anchors_left_out = std::erase_if(all, [min_reads](const AnchorWriter::Anchor& a) {
        return a.edit && a.reads.size() < min_reads;
    });
}

/**
 * The anchor file's table of read names: the distinct reads of the anchors written, sorted by
 * name.
 *
 * Each name is written once, in the table, and read rows refer to it by its row, since a read is
 * usually placed at several anchors. Rows follow sorted-name order, so that the file does not
 * depend on thread scheduling and the table can be binary-searched.
 */
struct FileReadTable {
    /// The names, in table order.
    vector<string_view> names;
    /// The row of each read in the table, by its index in read_names(). Meaningful only for
    /// reads in the table.
    vector<uint32_t> row_of;
};

/**
 * Make the table of the reads that the anchors hold.
 *
 * There are a few million names but hundreds of millions of read rows. A read's index in
 * read_names() is about as large as the number of names, so the reads that appear are marked in
 * a flat array, on several threads, and `row_of` is a flat array too.
 */
static FileReadTable build_file_read_table(const vector<AnchorWriter::Anchor>& all) {
    uint32_t max_read = 0;
#pragma omp parallel for schedule(dynamic, 1024) reduction(max : max_read)
    for (size_t i = 0; i < all.size(); ++i) {
        for (const AnchorWriter::ReadRow& r : all[i].reads) {
            max_read = max(max_read, r.read);
        }
    }
    const size_t read_slots = all.empty() ? 0 : (size_t)max_read + 1;
    vector<uint32_t> reads;
    {
        vector<std::atomic<uint8_t>> seen(read_slots);
#pragma omp parallel for schedule(dynamic, 1024)
        for (size_t i = 0; i < all.size(); ++i) {
            for (const AnchorWriter::ReadRow& r : all[i].reads) {
                seen[r.read].store(1, std::memory_order_relaxed);
            }
        }
        for (size_t read = 0; read < read_slots; ++read) {
            if (seen[read].load(std::memory_order_relaxed) != 0) {
                reads.push_back((uint32_t)read);
            }
        }
    }
    // Two reads never share a name, so the sorted order does not depend on how the parallel sort
    // shares out the work.
    const ReadNameTable& names = read_names();
    parallel_sort(reads.begin(), reads.end(), [&names](uint32_t a, uint32_t b) {
        return names.name(a) < names.name(b);
    });
    FileReadTable table;
    table.names.reserve(reads.size());
    table.row_of.resize(read_slots);
    for (size_t i = 0; i < reads.size(); ++i) {
        table.names.push_back(names.name(reads[i]));
        table.row_of[reads[i]] = (uint32_t)i;
    }
    return table;
}

/// Write the anchor file's header, which describes the run and the columns, and ends with the
/// table of read names.
static void write_header(ostream& out, const string& graph_name, const string& sample,
                         const string& reads_source, double mismap_min, const AnchorParams& params,
                         const FileReadTable& table) {
    // The format version, raised whenever the columns or their meaning change.
    out << "#anchors-version\t9\n";
    // Which vg wrote the file, since values can change between versions of vg while the format
    // version does not. Consumers skip unknown '#' lines.
    out << "#vg-version\t" << Version::get_version() << "\n";
    out << "#graph\t" << graph_name << "\n";
    out << "#sample\t" << sample << "\n";
    out << "#reads\t" << reads_source << "\n";
    out << "#mismap-min\t" << mismap_min << "\n";
    out << "#sites\t" << (params.het_only ? "het" : "het+hom")
        << (params.leaf_only ? ",leaf" : "") << "\n";
    // Name the rule that placed each read in its slot at a heterozygous site, so that files built
    // with different rules can be told apart.
    //   allele  -- the site's own allele match alone (`phase_hets` off)
    //   tilt    -- the allele match, weighted by the read's strand log-odds
    //   strand  -- the sign of the strand log-odds alone (`strict_hets`)
    const char* het_placement = "allele";
    if (params.strict_hets) {
        het_placement = "strand";
    } else if (params.phase_hets) {
        het_placement = "tilt";
    }
    // What was filtered out, so a consumer seeing no anchor at a site it expected one at can tell a
    // threshold from an absence. Written whether or not the defaults were changed.
    out << "#filters\tmin-reads=" << params.min_reads << " min-gqn=" << params.min_gqn
        << " min-read-score=" << params.min_read_score
        << " off-call=" << (params.keep_off_call ? "kept" : "dropped")
        << " end-pin-min-new=" << params.end_pin_min_new
        << " hom-split=" << (params.hom_split ? "on" : "off")
        << " het-placement=" << het_placement << "\n";
    if (params.hom_split) {
        // The two slots of a split homozygous site look like a heterozygote's, except that their
        // `allele` column is equal, so the header says that splitting was on.
        out << "#note\thom-split is ON: a homozygous site whose reads partition confidently by "
               "cross-site phase is written as TWO slots carrying the SAME allele. Slot is then a "
               "haplotype claim inferred from other sites, not read off this one -- held out on "
               "chr20 ONT it agrees with the allele partition 94.7% of the time, so roughly one "
               "read in twenty is on the wrong strand. Equal alleles across two slots is the tell.\n";
    }
    out << "#note\tan anchor is a zero-length pin at a point between two adjacent graph bases. A "
           "snarl pins at the junctions between its boundary nodes and its interior. An off-panel "
           "site (snarl column op...) pins around the bases its alleles change: after the last "
           "unchanged base before them and before the first unchanged base after them, an indel's "
           "whole repeat tract counting as changed. An insertion with no repeat tract changes no "
           "base and has one pin. There is no anchor sequence and no length.\n";
    out << "#note\tnode is the graph's own node ID. snarl is the VCF ID, which -N or -O writes in "
           "translated segment names, so under them node is not one of the IDs in snarl.\n";
    out << "#note\tR rows name their read by an integer id into the #read table above them, which "
           "is sorted by name. The table is the file's only copy of each name.\n";
    out << "#note\toffset is the 0-based index, in the read AS SEQUENCED, of the last base before "
           "the pin read in the site's direction. strand 0: the pin follows that base. strand 1: "
           "it precedes it. For the oriented read, use (read length - 1 - offset).\n";
    out << "#note\tR rows belong to the A row above them and do not repeat its key. The file is "
           "already ordered by node; do not sort it.\n";
    out << "#note\tno read point is in two anchors. For an R row the point is gap = offset+1 on "
           "strand 0 and offset on strand 1, the point between read[gap-1] and read[gap]. Where an "
           "off-panel site and a snarl would hold one alignment's point, only the snarl's anchor "
           "does; where two off-panel sites would, only the anchor written first does. vg reports "
           "on stderr how many placements it left out. Two alignments of one read name (paired-end "
           "mates, or a long read's supplementary alignments) share a read_id, so one read_id can "
           "name one gap twice, once for each alignment.\n";
    out << "#note\tslot indexes the settled PHASED pair, so slot i is field i of that snarl's GT in "
        << "the VCF -- slot 0 the left allele, slot 1 the right. A homozygous site collapses to one "
        << "slot holding every read, and carries no haplotype information to join on\n";
    out << "#note\ta nested HAPLOID site also holds one slot, but it is a haplotype: the chain is "
        << "on one strand of a diploid locus because the parent's other allele deletes it. Its slot "
        << "is the strand, matching the non-'.' field of the site's a|. or .|a GT. From v5 -- v4 "
        << "wrote 0 for both strands\n";
    out << "#note\tallele is the index into the site's candidate traversal set, which is NOT the "
        << "VCF ALT number: the ALT list is chosen after anchors are built. Use slot to join to GT\n";
    out << "#note\tgqn is SIGNED in [-1,1] and matches the VCF's GQN, including on the records the "
           "linkage layer moved: negative means linkage called against the reads, which is the "
           "population that carries a 37.8% false-positive rate against 8.6% overall. `.` means the "
           "site offered no gap to normalise, which is not the same as a negative value.\n";
    out << "#note\tgqn and explained are the site's, so they repeat across its slots. Anything else "
           "about the site is in the VCF, joinable on the snarl column, which is its ID.\n";
    out << "#note\tscore is the read's phred-scaled complement of the winning slot's share of its "
        << "own responsibility, the mismapping probability included in the denominator. So its "
        << "CEILING DEPENDS ON THE SITE. A site with one distinct allele -- a homozygote, split "
        << "or not, and a haploid -- gives the winner the whole (1-e) weight and caps at "
        << "phred(--mismap-min) = "
        << std::fixed << std::setprecision(2)
        << (mismap_min > 0.0 ? -10.0 * log10(mismap_min) : 99.0)
        << ". A HETEROZYGOTE splits that weight between two alleles and caps at "
        << "phred(e/(e+(1-e)/2)) = "
        << (mismap_min > 0.0
                ? -10.0 * log10(mismap_min / (mismap_min + (1.0 - mismap_min) / 2.0))
                : 99.0)
        << std::defaultfloat
        << ", measured p99 on chr20 ONT. Neither is 60 or 99. A threshold set from the "
        << "homozygous cap therefore discards EVERY heterozygous site, which is the opposite of "
        << "what a consumer filtering for phase information wants. 99 means the winner took the "
        << "whole share and no cap applied\n";
    out << "#note\treliability is the site's mean score over the reads it emitted, each read once "
        << "across both pins and both slots -- the same quantity --phase-min-q thresholds, default "
        << "9.5. It is low exactly where the reads cannot tell the site's alleles apart, which on "
        << "ONT means a 1 bp indel. Site-level, so it repeats across the site's rows. It is the "
        << "UNROUNDED mean, while the R rows it averages are written to one decimal, so "
        << "re-deriving it from them lands within about 0.05 rather than exactly. An off-panel "
        << "site's reliability counts the reads it placed before the rule above removed any. From "
        << "v6\n";
    out << "#reads-interned\t" << table.names.size() << "\n";
    out << "#note\tpin is where the anchor lies in its node: the node's bases before it, along the "
           "node's forward strand. A snarl pins at a node end; an off-panel site can pin inside a "
           "node. From v8. Off-panel changes on one node that the call bar admits, large indels "
           "aside, whose changed bases touch or overlap are written as one site where the reads "
           "spell alleles that carry each of them and they are not all indels at one anchor, with ID "
           "op<node>.<offset>.M<length>: offset is its first changed base along the node's forward "
           "strand, length the bases it replaces, and its ALTs are what the reads spell there. "
           "From v9\n";
    out << "#H\tA\tnode\tsnarl\tslot\tallele\tgqn\texplained\treliability\tpin\n";
    out << "#H\tR\tread_id\tstrand\toffset\tscore\n";
    for (size_t i = 0; i < table.names.size(); ++i) {
        out << "#read\t" << i << "\t" << table.names[i] << "\n";
    }
}

/// Write an anchor's A row and its R rows, to a stream set to fixed-point notation.
static void format_anchor(ostream& os, const AnchorWriter::Anchor& a, const vector<uint32_t>& row_of) {
    os << "A\t" << a.node << "\t" << a.snarl << "\t" << a.slot << "\t" << a.allele << "\t";
    if (std::isnan(a.gqn)) {
        os << ".";
    } else {
        os << std::setprecision(3) << a.gqn;
    }
    os << "\t" << std::setprecision(3) << a.explained << "\t";
    if (a.reliability < 0.0) {
        os << ".";
    } else {
        os << std::setprecision(2) << a.reliability;
    }
    os << "\t" << a.pin << "\n";
    for (const AnchorWriter::ReadRow& r : a.reads) {
        os << "R\t" << row_of[r.read] << "\t" << (int)r.direction << "\t" << r.offset << "\t"
           << std::setprecision(1) << r.score << "\n";
    }
}

/**
 * Write the anchors' rows, formatting them on several threads.
 *
 * Anchors are taken a batch at a time, and each batch is cut into pieces that threads format into
 * strings of their own. So the text held at once is one batch's rather than the whole file's,
 * since this runs when the run already uses the most memory.
 */
static void write_anchor_rows(ostream& out, const vector<AnchorWriter::Anchor>& all,
                              const vector<uint32_t>& row_of) {
    const size_t piece_anchors = 2048;
    const size_t batch_anchors = piece_anchors * (size_t)max(1, omp_get_max_threads());
    vector<string> pieces;
    for (size_t batch_start = 0; batch_start < all.size(); batch_start += batch_anchors) {
        size_t batch_end = min(all.size(), batch_start + batch_anchors);
        size_t piece_count = (batch_end - batch_start + piece_anchors - 1) / piece_anchors;
        pieces.assign(piece_count, string());
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t p = 0; p < piece_count; ++p) {
            ostringstream os;
            os << std::fixed;
            size_t from = batch_start + p * piece_anchors;
            size_t to = min(batch_end, from + piece_anchors);
            for (size_t i = from; i < to; ++i) {
                format_anchor(os, all[i], row_of);
            }
            pieces[p] = os.str();
        }
        for (const string& piece : pieces) {
            out.write(piece.data(), piece.size());
        }
    }
}

bool AnchorWriter::write(const string& path, const string& graph_name, const string& sample,
                         const string& reads_source, double mismap_min,
                         const AnchorParams& params) {
    vector<Anchor> all = take_sorted();
    leave_out_held_points(all, params.min_reads, written);
    written.anchors_written = all.size();
    written.placements_written = 0;
    for (const Anchor& a : all) {
        written.placements_written += a.reads.size();
    }
    const FileReadTable table = build_file_read_table(all);

    ofstream out(path);
    if (!out) {
        cerr << "error [vg call]: could not open " << path << " for the anchor output" << endl;
        return false;
    }
    write_header(out, graph_name, sample, reads_source, mismap_min, params, table);
    write_anchor_rows(out, all, table.row_of);
    return (bool)out;
}

void AnchorCollector::configure(string_view anchor_path, const AnchorParams& anchor_params,
                                string_view graph, string_view reads, double mismap_floor) {
    path = anchor_path;
    params = anchor_params;
    if (params.counters == nullptr) {
        params.counters = &owned_counters;
    }
    graph_name = graph;
    reads_source = reads;
    mismap_min = mismap_floor;
    if (!path.empty()) {
        // Make one queue per OpenMP thread, since the passes that fill it are parallel.
        writer = make_unique<AnchorWriter>((size_t)max(1, get_thread_count()));
    }
}

/**
 * Count, at one heterozygous site, how often each read's strand log-odds favour the slot its
 * allele put it in.
 *
 * This checks `AnchorParams::hom_split`, since at a heterozygous site the alleles show which
 * strand each read is on. Each read is counted once, though both pins of its slot hold it.
 *
 * @param strand_log_odds Called with a read name to get the read's strand log-odds, leaving the
 * site out.
 */
/// Return whether slots 0 and 1 of a site's anchors both exist and hold different alleles.
static bool slots_hold_two_alleles(const vector<AnchorWriter::Anchor>& anchors) {
    std::array<int, 2> allele_of_slot{-1, -1};
    for (const AnchorWriter::Anchor& anchor : anchors) {
        if (anchor.slot >= 0 && anchor.slot < 2) {
            allele_of_slot[anchor.slot] = anchor.allele;
        }
    }
    return allele_of_slot[0] >= 0 && allele_of_slot[1] >= 0 && allele_of_slot[0] != allele_of_slot[1];
}

template<typename StrandLogOdds>
static void check_phase_against_alleles(const vector<AnchorWriter::Anchor>& anchors,
                                        double phase_min, StrandLogOdds& strand_log_odds,
                                        AnchorPhaseCheckCounters& counters) {
    if (!slots_hold_two_alleles(anchors)) {
        return;
    }
    unordered_set<uint32_t> counted;
    for (const AnchorWriter::Anchor& anchor : anchors) {
        if (anchor.slot < 0 || anchor.slot > 1) {
            continue;
        }
        for (const AnchorWriter::ReadRow& row : anchor.reads) {
            if (counted.insert(row.read).second) {
                counters.count(strand_log_odds(read_names().name(row.read)), anchor.slot, phase_min);
            }
        }
    }
}

void AnchorCollector::collect(const AnchorSiteEvidence& evidence, const AnchorCall& call,
                              const string& site_name, bool is_leaf,
                              const multipass::ReadStrandTable& strands, size_t record_key) {
    if (path.empty() || writer == nullptr) {
        return;
    }
    if (params.leaf_only && !is_leaf) {
        return;
    }
    // The site's own log-odds for each read, built on the first lookup and shared by the rest, so
    // that looking up every read at a site costs one pass over its reads rather than one per
    // read.
    unordered_map<uint64_t, double> site_own;
    // 0 when `site_own` is not built yet, 1 when it is, and 2 when the site has nothing to leave
    // out.
    int site_own_state = 0;
    auto strand_log_odds = [&](std::string_view read_name) {
        if (site_own_state == 0) {
            site_own_state = strands.site_own_strand_log_odds(record_key, site_own) ? 1 : 2;
        }
        return strands.read_strand_log_odds(record_key, read_name,
                                            site_own_state == 1 ? &site_own : nullptr);
    };
    const vector<int>& genotype = call.genotype;
    // Each read's strand log-odds, leaving out this site.
    vector<double> read_strand;
    // Fill it only where `build_site_anchors` reads it: at a diploid homozygote that may be split,
    // or at a heterozygous site with `phase_hets` or `strict_hets` on. `build_site_anchors` checks
    // the length of the vector, so an empty one is ignored.
    const bool splittable_hom = genotype.size() == 2 && genotype[0] == genotype[1];
    if (const bool tiltable_het = (params.phase_hets || params.strict_hets)
                                  && genotype.size() == 2 && genotype[0] != genotype[1];
        (params.hom_split && splittable_hom) || tiltable_het) {
        read_strand.reserve(evidence.reads.size());
        for (const AnchorRead& read : evidence.reads) {
            read_strand.push_back(strand_log_odds(read_names().name(read.read)));
        }
    }
    vector<AnchorWriter::Anchor> anchors;
    build_site_anchors(evidence, call, site_name, params, anchors,
                       (params.hom_split || params.phase_hets || params.strict_hets)
                           ? &read_strand : nullptr);
    // Check the strand log-odds only when splitting is on, and not when the heterozygous
    // placement itself uses them, since the check would then compare the strand with itself.
    if (params.hom_split && !params.phase_hets && !params.strict_hets && anchors.size() >= 2) {
        check_phase_against_alleles(anchors, params.phase_min, strand_log_odds,
                                    params.counters->phase_check);
    }
    for (AnchorWriter::Anchor& anchor : anchors) {
        writer->add(std::move(anchor));
    }
}

void AnchorCollector::write(const string& sample_name) {
    if (path.empty() || writer == nullptr) {
        return;
    }
    if (!writer->write(path, graph_name, sample_name, reads_source, mismap_min, params)) {
        return;
    }
    cerr << "[vg call] anchors: " << writer->tally().anchors_written << " written over "
         << writer->tally().placements_written << " read placements to " << path << "; "
         << writer->tally().placements_left_out << " off-panel read placements left out because a "
         << "snarl's anchor or an earlier off-panel anchor holds the same read point, and then "
         << writer->tally().anchors_left_out << " off-panel anchors below the minimum reads" << endl;
    if (params.counters != nullptr) {
        params.counters->report(cerr);
    }
}

}
