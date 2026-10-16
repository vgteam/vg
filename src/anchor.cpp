#include "anchor.hpp"
#include "read_phasing.hpp"
#include "read_strand_table.hpp"

#include "version.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <functional>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

#include "parallel_sort.hpp"
#include "path.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

void AnchorCounters::report(ostream& out) const {
    out << "[vg call] anchors: " << verified.load() << " pins verified against the graph, "
        << verify_failed.load() << " failed" << endl;
    out << "[vg call] anchor pins refused: " << no_visit.load()
        << " read does not visit the boundary node, " << repeat_visit.load()
        << " visits it more than once, " << unaligned_base.load()
        << " no read base aligned to the pin's node base, " << no_neighbour.load()
        << " nothing upstream of an entry pin" << endl;
    out << "[vg call] anchor reads: " << coincident.load()
        << " with both pins at one position (kept at S, dropped at E), " << off_call.load()
        << " fitting no called allele";
    if (shared_name.load() > 0) {
        out << ", " << shared_name.load()
            << " sharing a read name with another read at the same site (paired mates do; the name "
               "then does not identify an alignment)";
    }
    if (degenerate_site.load() > 0) {
        out << ", " << degenerate_site.load() << " sites whose two boundaries are one node";
    }
    if (single_pin.load() > 0) {
        out << ", " << single_pin.load()
            << " sites given only their start pin because the end pin held no reads of its own";
    }
    out << endl;
    if (hom_split.load() > 0 || hom_unsplit.load() > 0) {
        const size_t s = hom_split.load(), u = hom_unsplit.load();
        out << "[vg call] anchors: " << s << " homozygous sites split by read phase, " << u
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
    if (phase_checked.load() > 0) {
        const size_t n = phase_checked.load(), ok = phase_agree.load();
        const size_t cn = phase_confident.load(), ck = phase_confident_agree.load();
        out << "[vg call] anchors: cross-site phase reproduces the allele partition at het "
               "sites (held out, the split's own accuracy -- reported only when het placement "
               "does NOT use the strand, or it would measure itself): "
            << ok << "/" << n << " agree (" << (100.0 * (double)ok / (double)n) << "%)";
        if (cn > 0) {
            out << ", confident " << ck << "/" << cn << " ("
                << (100.0 * (double)ck / (double)cn) << "%)";
        }
        out << ", " << phase_no_opinion.load() << " reads with no cross-site opinion" << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Pin resolution
////////////////////////////////////////////////////////////////////////////////

/// The first and last read base of one mapping that consumes a node base, or -1 for a mapping that
/// consumes none, as for a node the read deletes.
///
/// Kept per mapping rather than as one run over the read, because the entry pin must step exactly
/// one node back along the read's walk, and a merged run cannot tell the previous node's last base
/// from the last base before several deleted nodes.
struct MappingExtent {
    int64_t first_consuming = -1;
    int64_t last_consuming = -1;
};

static MappingExtent extent_of(const Mapping& m, size_t read_start) {
    MappingExtent out;
    size_t read_pos = read_start;
    for (int64_t j = 0; j < m.edit_size(); ++j) {
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

AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters, int64_t pin) {
    AnchorPlacement out;
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();

    // Find the pin node's mapping and where in the read it starts. The pin node is a boundary of
    // the site, so every visit to it is in the site's index when there is one. A read that visits
    // it more than once is refused.
    int64_t hit = -1;
    size_t hit_read_start = 0;
    size_t visits = 0;
    // Read start of the mapping before the pin's, needed only by the entry pin, and noted here
    // because the unindexed walk passes it only once.
    size_t before_read_start = 0;
    if (read.indexed()) {
        for (size_t k = 0; k < read.mapping_count; ++k) {
            int64_t i = (int64_t)read.mappings[k];
            if (path.mapping(i).position().node_id() != node_id) {
                continue;
            }
            ++visits;
            hit = i;
            hit_read_start = read.read_offsets[i];
        }
        if (hit > 0) {
            before_read_start = read.read_offsets[hit - 1];
        }
    } else {
        size_t read_pos = 0;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            const Mapping& m = path.mapping(i);
            if (m.position().node_id() == node_id) {
                ++visits;
                hit = i;
                hit_read_start = read_pos;
            }
            read_pos += (size_t)mapping_to_length(m);
        }
        if (hit > 0) {
            before_read_start = hit_read_start
                                - (size_t)mapping_to_length(path.mapping(hit - 1));
        }
    }

    if (visits == 0) {
    // A read that starts inside the site does not visit the start boundary, and one that ends
    // inside it does not visit the end boundary.
        ++counters.no_visit;
        return out;
    }
    if (visits > 1) {
    // Which visit the pin belongs to is ambiguous, so the read is refused.
        ++counters.repeat_visit;
        return out;
    }

    const Mapping& m = path.mapping(hit);
    const bool is_rev = m.position().is_reverse();
    handle_t handle = graph.get_handle(node_id, false);
    const size_t node_len = graph.get_length(handle);
    if (node_len == 0) {
        ++counters.unaligned_base;
        return out;
    }

    // The node base the pin is defined against, in node-forward coordinates.
    //   exit  (S) pin: the node's last base in the site's direction; the pin follows it.
    //   entry (E) pin: the node's first base in the site's direction; the pin precedes it.
    if (pin >= 0) {
        // A pin inside the node is resolved as an exit pin, against the base before it.
        exit_pin = true;
        const int64_t before = site_backward ? pin : pin - 1;
        if (before < 0 || before >= (int64_t)node_len) {
            ++counters.unaligned_base;
            return out;
        }
    }
    const size_t q_fwd = pin >= 0 ? (size_t)(site_backward ? pin : pin - 1)
                                  : exit_pin ? (site_backward ? 0 : node_len - 1)
                                             : (site_backward ? node_len - 1 : 0);
    // The same base in the orientation the read visits the node in, in which the mapping's offsets
    // are measured.
    const size_t q_vis = is_rev ? (node_len - 1 - q_fwd) : q_fwd;

    // Walk the mapping's edits for the read base aligned to it.
    size_t node_pos = (size_t)m.position().offset();
    size_t rp = hit_read_start;
    int64_t base_index = -1;
    bool exact_edit = false;
    char edit_base = 0;
    for (int64_t j = 0; j < m.edit_size(); ++j) {
        const Edit& e = m.edit(j);
        const size_t from = (size_t)e.from_length();
        const size_t to = (size_t)e.to_length();
        if (from > 0 && q_vis >= node_pos && q_vis < node_pos + from) {
            if (from == to) {
                base_index = (int64_t)(rp + (q_vis - node_pos));
                if (e.sequence().empty()) {
                    exact_edit = true;
                } else if (e.sequence().size() == to) {
                    edit_base = e.sequence()[q_vis - node_pos];
                }
            }
            // Otherwise the node base is deleted (to == 0) or sits inside a replacement of unequal
            // length, and no read base corresponds to it.
            break;
        }
        node_pos += from;
        rp += to;
    }
    if (base_index < 0 || (size_t)base_index >= aln.sequence().size()) {
        ++counters.unaligned_base;
        return out;
    }

    // The check, made here because the read is not kept. On a match edit the read base must be the
    // graph's base, complemented when the read visits the node in reverse, which catches offset and
    // strand errors. On a mismatch edit the base differs from the graph's, so the edit's recorded
    // sequence is compared instead, which still catches offset errors.
    const char read_base = (char)toupper(aln.sequence()[base_index]);
    if (exact_edit) {
        const char node_base = graph.get_base(handle, q_fwd);
        const char want = (char)toupper(is_rev ? reverse_complement(node_base) : node_base);
        if (read_base != want) {
            ++counters.verify_failed;
            return out;
        }
        ++counters.verified;
    } else if (edit_base != 0) {
        if (read_base != (char)toupper(edit_base)) {
            ++counters.verify_failed;
            return out;
        }
        ++counters.verified;
    }

    // Does the read cross the pin in the site's direction? `is_rev` is relative to the node, the
    // site's own orientation is `site_backward`, and the strand is the comparison of the two.
    out.direction = (is_rev == site_backward) ? 0 : 1;

    if (exit_pin) {
        // The reference base is immediately upstream of the pin, which is what `offset` reports.
        out.offset = base_index;
        return out;
    }

    // Entry pin: the reference base is downstream of the pin, so the offset is the read base just
    // upstream, which is the adjacent base on the previous node of the read's walk. Stepping back by
    // read index instead would pass over a node the read deletes, which has no read bases, and land
    // on the neighbouring snarl's pin. So we step one mapping along the walk, skipping insertions
    // and soft clips, which consume no node base; if that mapping consumes no node base at all, the
    // node upstream is deleted in this read and the read is refused.
    const int64_t neighbour = (out.direction == 0) ? hit - 1 : hit + 1;
    if (neighbour < 0 || neighbour >= (int64_t)path.mapping_size()) {
        // The read begins (or ends) here, so it does not cross the pin.
        ++counters.no_neighbour;
        return out;
    }
    const size_t neighbour_read_start =
        (out.direction == 0) ? before_read_start
                          : hit_read_start + (size_t)mapping_to_length(path.mapping(hit));
    const MappingExtent adjacent = extent_of(path.mapping(neighbour), neighbour_read_start);
    const int64_t candidate =
        (out.direction == 0) ? adjacent.last_consuming : adjacent.first_consuming;
    if (candidate < 0 || candidate >= (int64_t)aln.sequence().size()) {
        ++counters.no_neighbour;
        return out;
    }
    out.offset = candidate;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Retention
////////////////////////////////////////////////////////////////////////////////

uint32_t ReadNameTable::intern(const string& name) {
    size_t s = std::hash<string_view>{}(name) % SHARDS;
    Shard& shard = shards[s];
    lock_guard<std::mutex> guard(shard.mutex);
    auto found = shard.index.find(name);
    if (found != shard.index.end()) {
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
    uint32_t index = (uint32_t)(position * SHARDS + s);
    shard.index.emplace(string_view(chunk[position % CHUNK]), index);
    return index;
}

string_view ReadNameTable::name(uint32_t index) const {
    const Shard& shard = shards[index % SHARDS];
    size_t position = index / SHARDS;
    const string* chunk = shard.chunks[position / CHUNK].load(std::memory_order_acquire);
    return chunk[position % CHUNK];
}

ReadNameTable::Shard::~Shard() {
    for (auto& chunk : chunks) {
        delete[] chunk.load();
    }
}

ReadNameTable& read_names() {
    // Never destroyed: the table lives until the process ends, and freeing a few million names one
    // by one at exit only delays it.
    static ReadNameTable* table = new ReadNameTable();
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

void build_site_anchors(const AnchorSiteEvidence& evidence, const vector<int>& genotype,
                        const string& snarl_id, double gqn, double explained, int haploid_slot,
                        const AnchorParams& params, AnchorCounters& counters,
                        vector<AnchorWriter::Anchor>& out,
                        const vector<double>* read_strand) {

    if (genotype.empty() || evidence.reads.empty() || evidence.n_alleles == 0) {
        return;
    }
    for (int a : genotype) {
        if (a < 0 || (size_t)a >= evidence.n_alleles) {
            // A star or missing allele: there is no traversal to assign reads to.
            return;
        }
    }
    // `hom` is also true for a haploid call, which `het_only` should drop too. `haploid` is kept
    // separately because its one slot names a strand.
    const bool haploid = genotype.size() == 1;
    bool hom = true;
    for (size_t i = 1; i < genotype.size(); ++i) {
        if (genotype[i] != genotype[0]) {
            hom = false;
        }
    }
    if (hom && params.het_only) {
        return;
    }
    // NaN means the site had no gap to normalise, which is different from zero, and is filtered
    // out. A negative GQN, on a record whose genotype the linkage model changed, is a real value
    // and is compared like any other.
    if (params.min_gqn > 0.0 && (std::isnan(gqn) || gqn < params.min_gqn)) {
        return;
    }

    // One slot per strand of the genotype, or one slot holding every read for a homozygote.
    vector<int> slot_allele;
    bool phase_split = false;
    if (hom) {
        // Split only when enough confidently placed reads fall on each side.
        if (params.hom_split && !haploid && genotype.size() == 2 && read_strand != nullptr
            && read_strand->size() == evidence.reads.size()) {
            size_t side0 = 0, side1 = 0;
            for (size_t r = 0; r < evidence.reads.size(); ++r) {
                // Only reads that will be written: the placement loop below drops a read pinned at
                // neither end.
                if (!evidence.reads[r].start_pin.placed()
                    && !evidence.reads[r].end_pin.placed()) {
                    continue;
                }
                const double lo = (*read_strand)[r];
                if (std::isnan(lo) || lo == 0.0) {
                    // NaN is a strand from another phase set, and exactly zero is no strand
                    // log-odds at all: no table, or no phased site reached. Neither is evidence for
                    // strand 1, though the strict `> 0.0` test below would count a zero there.
                    continue;
                }
                if (std::abs(lo) < params.phase_min) {
                    continue;
                }
                if (lo > 0.0) {
                    ++side0;
                } else {
                    ++side1;
                }
            }
            if (side0 >= params.phase_min_side && side1 >= params.phase_min_side) {
                phase_split = true;
                ++counters.hom_split;
            } else {
                ++counters.hom_unsplit;
            }
        }
        slot_allele.push_back(genotype[0]);
        if (phase_split) {
            slot_allele.push_back(genotype[0]);   // same allele on both strands, by definition
        }
    } else {
        slot_allele.assign(genotype.begin(), genotype.end());
    }
    const size_t n_slots = slot_allele.size();

    // Where a single slot goes. A homozygote's is slot 0 and names no strand; a nested haploid
    // site's is the strand `nested_strand` gave, which is 1 for a `.|a` site. Only 0 or 1 is
    // written, since slot indexes a GT field.
    const int base_slot = (haploid && haploid_slot == 1) ? 1 : 0;

    const vector<double> weight = allele_length_weights(
        evidence.allele_length, evidence.n_alleles, evidence.mean_read_length,
        evidence.length_weighted, slot_allele);

    // Two anchors per slot: one at each pin. The end pin is skipped where the site's two boundaries
    // are the same node, since the anchors would then be indistinguishable; a site inside one node
    // instead has its pins on either side of its base.
    bool degenerate = evidence.start_node == evidence.end_node && !evidence.inside_node;
    if (degenerate) {
        ++counters.degenerate_site;
    }
    vector<AnchorWriter::Anchor> start_anchors(n_slots), end_anchors(n_slots);
    for (size_t i = 0; i < n_slots; ++i) {
        start_anchors[i].node = evidence.start_node;
        start_anchors[i].pin = evidence.start_pin;
        start_anchors[i].snarl = snarl_id;
        start_anchors[i].slot = (int)i + base_slot;
        start_anchors[i].allele = slot_allele[i];
        start_anchors[i].gqn = gqn;
        start_anchors[i].explained = explained;
        end_anchors[i] = start_anchors[i];
        end_anchors[i].node = evidence.end_node;
        end_anchors[i].pin = evidence.end_pin;
    }

    // Count reads sharing a name at this site, which cannot be told apart once the file is
    // written.
    {
        unordered_map<uint32_t, size_t> name_count;
        name_count.reserve(evidence.reads.size() * 2);
        for (const AnchorRead& read : evidence.reads) {
            ++name_count[read.read];
        }
        for (const auto& entry : name_count) {
            if (entry.second > 1) {
                counters.shared_name += entry.second;
            }
        }
    }

    // For each slot, the reads its end anchor holds that its start anchor does not; see
    // AnchorParams::end_pin_min_new. Counted per slot rather than per site, since anchors are
    // written and filtered by `min_reads` per slot.
    vector<size_t> new_at_end(n_slots, 0);

    for (size_t r = 0; r < evidence.reads.size(); ++r) {
        const AnchorRead& read = evidence.reads[r];
        if (!read.start_pin.placed() && !read.end_pin.placed()) {
            continue;
        }
        const double mismap = (double)read.mismap;

        // Responsibilities over the called alleles, plus the mismapping outcome, which bounds the
        // confidence at the heterozygous score ceiling.
        double best_resp = -1.0;
        double best_plain = 0.0;
        size_t best_slot = 0;
        double total = mismap;
        if (phase_split) {
            // Both slots spell the same allele, so the responsibilities cannot choose between them,
            // and the read's strand log-odds decide instead.
            //
            // The share is computed over the site's one distinct allele rather than over the two
            // slots, so that a read's confidence, and so `min_read_score` and the reliability
            // column, mean the same at split sites as elsewhere.
            const double lo = (*read_strand)[r];
            int coin = -1;
            if (std::isnan(lo)) {
                // The read has a strand only in another phase set, or in more than one, and phase
                // sets number their strands independently. A coin would put it on a strand of this
                // phase set with nothing behind it: where the read also crosses a heterozygous site
                // of this phase set, whose allele places it, the coin contradicts that slot half
                // the time, and otherwise it joins two phase sets that the phasing left apart. It
                // is dropped.
                ++counters.hom_split_no_opinion;
                continue;
            }
            if (lo == 0.0) {
                // No strand log-odds: the read reached no other phased site, or none whose alleles
                // it can tell apart. The site is homozygous, so either slot spells the read's
                // allele, and the read is placed by a coin flip derived from a hash of its name.
                // The same read then gets the same slot at every site and both pins, so it stays on
                // one strand, and it is not missing from some anchors while present at others.
                ++counters.hom_split_no_opinion;
                ++counters.hom_split_coin;
                coin = (int)(std::hash<string_view>{}(read_names().name(read.read)) & 1ull);
            }
            best_resp = (1.0 - mismap) * (double)evidence.rel_at(r, (size_t)slot_allele[0]);
            total = mismap + best_resp;
            best_slot = coin >= 0 ? (size_t)coin : (lo > 0.0 ? 0 : 1);
        } else {
        // Which slot this read joins, under one of three rules: the allele match alone, the
        // allele match tilted towards the strand the read's log-odds favour, or the strand alone.
        // Whichever rule picks the slot, the read's confidence, and so the reliability column,
        // comes from the untilted share `plain`. Reliability measures whether the site's own
        // reads tell its alleles apart, so the strand log-odds must not enter it.
        const bool tiltable = read_strand != nullptr && n_slots == 2
                              && read_strand->size() == evidence.reads.size()
                              && slot_allele[0] != slot_allele[1];
        const double lo = tiltable ? (*read_strand)[r] : 0.0;
        const bool has_opinion = tiltable && !std::isnan(lo) && lo != 0.0;

        // The default weighting: multiply the disfavoured slot's weight by exp(-|lo|), which gives the
        // same weights `phase_aware_correction` uses in re-genotyping, computed without overflow.
        double w0 = weight[0];
        double w1 = n_slots > 1 ? weight[1] : 0.0;
        if (params.phase_hets && !params.strict_hets && has_opinion) {
            const double f = std::exp(-std::abs(lo));
            if (lo > 0.0) {
                w1 *= f;
            } else {
                w0 *= f;
            }
            ++counters.het_phase_tilted;
        }
        for (size_t i = 0; i < n_slots; ++i) {
            const double plain = (1.0 - mismap) * weight[i]
                                 * (double)evidence.rel_at(r, (size_t)slot_allele[i]);
            const double resp = (1.0 - mismap) * (i == 0 ? w0 : (i == 1 ? w1 : weight[i]))
                                * (double)evidence.rel_at(r, (size_t)slot_allele[i]);
            total += plain;
            // Ties break on the allele, not on the slot position. `slot_allele` is in phase order,
            // which read phasing can change, so a tie-break by position would make the division of
            // reads depend on the phase. Exact equality is intended: a tie is the same `rel` value
            // for both alleles.
            if (resp > best_resp
                || (resp == best_resp && slot_allele[i] < slot_allele[best_slot])) {
                best_resp = resp;
                best_plain = plain;
                best_slot = i;
            }
        }
        // With `strict_hets` on, use only the sign of the strand log-odds, ignoring the allele
        // match. A read with no strand log-odds keeps its allele-match slot, so both rules place the
        // same reads.
        if (params.strict_hets && has_opinion) {
            const size_t want = lo > 0.0 ? 0 : 1;
            if (want != best_slot) {
                ++counters.het_strict_moved;
            }
            best_slot = want;
            best_plain = (1.0 - mismap) * weight[want]
                         * (double)evidence.rel_at(r, (size_t)slot_allele[want]);
            ++counters.het_phase_tilted;
        }
        best_resp = best_plain;   // from here on it is the untilted share's numerator
        }
        if (!(best_resp > 0.0) || !(total > 0.0)) {
            // The read fits no called allele at all.
            ++counters.off_call;
            continue;
        }

        // Is the read's best-fitting allele, over every scored traversal, a called one? A read that
        // prefers an uncalled allele is not placed on a called one.
        double row_best = 0.0;
        for (size_t a = 0; a < evidence.n_alleles; ++a) {
            row_best = max(row_best, (double)evidence.rel_at(r, a));
        }
        bool called_is_best = false;
        for (size_t i = 0; i < n_slots; ++i) {
            if ((double)evidence.rel_at(r, (size_t)slot_allele[i]) >= row_best) {
                called_is_best = true;
            }
        }
        if (!called_is_best) {
            ++counters.off_call;
            if (!params.keep_off_call) {
                continue;
            }
        }

        const double share = best_resp / total;
        double score = 99.0;
        if (share < 1.0) {
            score = -10.0 * log10(1.0 - share);
        }
        if (score < params.min_read_score) {
            continue;
        }

        AnchorWriter::ReadRow row;
        row.read = read.read;
        row.score = (float)score;

        if (read.end_pin.placed() && !read.start_pin.placed()) {
            ++new_at_end[best_slot];
        }

        bool wrote_start = false;
        if (read.start_pin.placed()) {
            row.offset = read.start_pin.offset;
            row.direction = read.start_pin.direction;
            start_anchors[best_slot].reads.push_back(row);
            wrote_start = true;
        }
        if (!degenerate && read.end_pin.placed()) {
            // A read that crosses straight from the start boundary to the end boundary carries an
            // allele deleting the site's interior, so both its pins are at one read position. It is
            // kept at the start pin and dropped here, so that a read position is in at most one anchor.
            if (wrote_start && read.end_pin.offset == read.start_pin.offset
                && read.end_pin.direction == read.start_pin.direction) {
                ++counters.coincident;
            } else {
                row.offset = read.end_pin.offset;
                row.direction = read.end_pin.direction;
                end_anchors[best_slot].reads.push_back(row);
            }
        }
    }

    const size_t out_begin = out.size();
    for (size_t i = 0; i < n_slots; ++i) {
        if (start_anchors[i].reads.size() >= params.min_reads) {
            out.push_back(std::move(start_anchors[i]));
        }
        // Drop this slot's end anchor where it has too few reads of its own, decided from the reads
        // actually placed.
        bool drop_end = degenerate;
        if (!drop_end && params.end_pin_min_new > 0
            && new_at_end[i] < params.end_pin_min_new) {
            drop_end = true;
            ++counters.single_pin;
        }
        if (!drop_end && end_anchors[i].reads.size() >= params.min_reads) {
            out.push_back(std::move(end_anchors[i]));
        }
    }

    // The site's reliability: the mean confidence of the reads it wrote, each read once. Computed
    // from `out`, after the filters above have removed slots and end anchors, so that it averages
    // the reads a consumer can see. Reads are keyed by name, which counts a read once across both
    // pins and slots, and paired mates once, at the higher of their scores, as `merge_mates` keeps
    // the more confident mate.
    if (out.size() > out_begin) {
        // Keyed by the name itself, not its index in read_names(): the scores are summed in the
        // map's iteration order, which depends on the keys' hashes, and a floating-point sum's last
        // digits depend on that order. Keying by name keeps the sum, and the file, as before.
        unordered_map<string_view, float> per_read;
        for (size_t i = out_begin; i < out.size(); ++i) {
            for (const AnchorWriter::ReadRow& row : out[i].reads) {
                auto placed = per_read.emplace(read_names().name(row.read), row.score);
                if (!placed.second && row.score > placed.first->second) {
                    placed.first->second = row.score;
                }
            }
        }
        double sum = 0.0;
        for (const auto& entry : per_read) {
            sum += (double)entry.second;
        }
        const double site_reliability =
            per_read.empty() ? -1.0 : sum / (double)per_read.size();
        for (size_t i = out_begin; i < out.size(); ++i) {
            out[i].reliability = site_reliability;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Writer
////////////////////////////////////////////////////////////////////////////////

AnchorWriter::AnchorWriter(size_t threads) : queues(max((size_t)1, threads)) {
}

void AnchorWriter::add(Anchor&& anchor) {
    size_t thread = (size_t)max(0, omp_get_thread_num());
    queues[thread % queues.size()].push_back(std::move(anchor));
}

size_t AnchorWriter::anchor_count() const {
    size_t total = 0;
    for (const auto& q : queues) {
        total += q.size();
    }
    return total;
}

size_t AnchorWriter::read_row_count() const {
    size_t total = 0;
    for (const auto& q : queues) {
        for (const Anchor& a : q) {
            total += a.reads.size();
        }
    }
    return total;
}

/// -1, 0 or 1 as `a` is below, equal to or above `b`, with NaN below every number and equal to
/// itself, so that it can order values that may be missing.
static int compare_values(double a, double b) {
    const bool a_nan = std::isnan(a), b_nan = std::isnan(b);
    if (a_nan || b_nan) {
        return a_nan == b_nan ? 0 : (a_nan ? -1 : 1);
    }
    return a < b ? -1 : (b < a ? 1 : 0);
}

/// Whether anchor `a` comes before `b` on what the file writes for them besides node, snarl and
/// slot: allele, gqn, explained, reliability, then the read rows in order (read name, direction,
/// offset, score). Two anchors equal on all of it are written identically, so their order cannot
/// show.
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
        if (x.read != y.read) {
            const string_view xn = table.name(x.read), yn = table.name(y.read);
            if (xn != yn) {
                return xn < yn;
            }
        }
        if (x.direction != y.direction) {
            return x.direction < y.direction;
        }
        if (x.offset != y.offset) {
            return x.offset < y.offset;
        }
        const int c = compare_values(x.score, y.score);
        if (c != 0) {
            return c < 0;
        }
    }
    return false;
}

bool AnchorWriter::write(const string& path, const string& graph_name, const string& sample,
                         const string& reads_source, double mismap_min,
                         const AnchorParams& params) {
    vector<Anchor> all;
    all.reserve(anchor_count());
    for (auto& q : queues) {
        std::move(q.begin(), q.end(), back_inserter(all));
        q.clear();
    }
    // Within an anchor, by read name, so that the file does not depend on thread scheduling. Each
    // anchor's reads are sorted on their own, so the anchors are shared out among threads. Done
    // before the anchors are sorted, so that the sort below can compare whole anchors.
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < all.size(); ++i) {
        const ReadNameTable& table = read_names();
        sort(all[i].reads.begin(), all[i].reads.end(), [&](const ReadRow& x, const ReadRow& y) {
            if (x.read != y.read) {
                return table.name(x.read) < table.name(y.read);
            }
            return x.offset < y.offset;
        });
    }
    // Ordered by node ID, then by snarl and slot, since two snarls can share a boundary node. The
    // sort runs on several threads. Anchors that agree on node, snarl and slot are put in order by
    // the rest of what is written for them, so that the order is fully determined: it cannot depend
    // on how the sort shares out the work, nor on which thread collected which anchor.
    const ReadNameTable& names_table = read_names();
    parallel_sort(all.begin(), all.end(), [&](const Anchor& a, const Anchor& b) {
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
        return written_before(a, b, names_table);
    });

    // Intern the read names: each name is written once, in a table, and read rows refer to it by
    // index, since a read is usually placed at several anchors. Indices follow sorted-name order,
    // so that the file does not depend on thread scheduling and the table can be binary-searched.
    //
    // The table is the distinct names, sorted. There are a few million names but hundreds of
    // millions of rows, and a read's index in read_names() is about as large as the number of
    // names, so the reads that appear are marked in a flat array, on several threads, and a row
    // finds its place in the table through another flat array.
    uint32_t max_read = 0;
#pragma omp parallel for schedule(dynamic, 1024) reduction(max : max_read)
    for (size_t i = 0; i < all.size(); ++i) {
        for (const ReadRow& r : all[i].reads) {
            max_read = max(max_read, r.read);
        }
    }
    const size_t read_slots = all.empty() ? 0 : (size_t)max_read + 1;
    vector<uint32_t> reads;
    {
        vector<std::atomic<uint8_t>> seen(read_slots);
#pragma omp parallel for schedule(dynamic, 1024)
        for (size_t i = 0; i < all.size(); ++i) {
            for (const ReadRow& r : all[i].reads) {
                seen[r.read].store(1, std::memory_order_relaxed);
            }
        }
        for (size_t read = 0; read < read_slots; ++read) {
            if (seen[read].load(std::memory_order_relaxed) != 0) {
                reads.push_back((uint32_t)read);
            }
        }
    }
    {
        // Two reads never share a name, so the sorted order does not depend on how the parallel
        // sort shares out the work.
        const ReadNameTable& table = read_names();
        parallel_sort(reads.begin(), reads.end(), [&](uint32_t a, uint32_t b) {
            return table.name(a) < table.name(b);
        });
    }
    vector<string_view> names;
    names.reserve(reads.size());
    // Each read's place in the table, by its index in read_names().
    vector<uint32_t> table_index(read_slots);
    for (size_t i = 0; i < reads.size(); ++i) {
        names.push_back(read_names().name(reads[i]));
        table_index[reads[i]] = (uint32_t)i;
    }

    ofstream out(path);
    if (!out) {
        cerr << "error [vg call]: could not open " << path << " for the anchor output" << endl;
        return false;
    }
    // The format version, raised whenever the columns or their meaning change.
    out << "#anchors-version\t8\n";
    // Which vg wrote the file, since values can change between versions of vg while the format
    // version does not. Consumers skip unknown '#' lines.
    out << "#vg-version\t" << Version::get_version() << "\n";
    out << "#graph\t" << graph_name << "\n";
    out << "#sample\t" << sample << "\n";
    out << "#reads\t" << reads_source << "\n";
    out << "#mismap-min\t" << mismap_min << "\n";
    out << "#sites\t" << (params.het_only ? "het" : "het+hom")
        << (params.leaf_only ? ",leaf" : "") << "\n";
    // What was filtered out, so a consumer seeing no anchor at a site it expected one at can tell a
    // threshold from an absence. Written whether or not the defaults were changed.
    out << "#filters\tmin-reads=" << params.min_reads << " min-gqn=" << params.min_gqn
        << " min-read-score=" << params.min_read_score
        << " off-call=" << (params.keep_off_call ? "kept" : "dropped")
        << " end-pin-min-new=" << params.end_pin_min_new
        << " hom-split=" << (params.hom_split ? "on" : "off")
        // Which rule placed each read in its slot at a heterozygous site, so that files built with
        // different rules can be told apart:
        //   allele  -- the site's own allele match alone (`phase_hets` off)
        //   tilt    -- the allele match, weighted by the read's strand log-odds
        //   strand  -- the sign of the strand log-odds alone (`strict_hets`)
        << " het-placement="
        << (params.strict_hets ? "strand" : (params.phase_hets ? "tilt" : "allele")) << "\n";
    if (params.hom_split) {
        // The two slots of a split homozygous site look like a heterozygote's, except that their
        // `allele` column is equal, so the header says that splitting was on.
        out << "#note\thom-split is ON: a homozygous site whose reads partition confidently by "
               "cross-site phase is written as TWO slots carrying the SAME allele. Slot is then a "
               "haplotype claim inferred from other sites, not read off this one -- held out on "
               "chr20 ONT it agrees with the allele partition 94.7% of the time, so roughly one "
               "read in twenty is on the wrong strand. Equal alleles across two slots is the tell.\n";
    }
    out << "#note\tan anchor is a zero-length pin at the junction between a snarl boundary node "
           "and the site interior. There is no anchor sequence and no length.\n";
    out << "#note\tnode is the graph's own node ID. snarl is the VCF ID, which -N or -O writes in "
           "translated segment names, so under them node is not one of the IDs in snarl.\n";
    out << "#note\tR rows name their read by an integer id into the #read table above them, which "
           "is sorted by name. The table is the file's only copy of each name.\n";
    out << "#note\toffset is the 0-based index, in the read AS SEQUENCED, of the last base before "
           "the pin read in the site's direction. strand 0: the pin follows that base. strand 1: "
           "it precedes it. For the oriented read, use (read length - 1 - offset).\n";
    out << "#note\tR rows belong to the A row above them and do not repeat its key. The file is "
           "already ordered by node; do not sort it.\n";
    out << "#note\tno (read, strand, offset) appears in more than one anchor -- for distinct "
           "alignments. A read NAME shared by two alignments (paired-end mates share one) can "
           "repeat a position; the shared_name counter reports how many such reads there were.\n";
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
        << "re-deriving it from them lands within about 0.05 rather than exactly. From v6\n";
    out << "#reads-interned\t" << reads.size() << "\n";
    out << "#note\tpin is where the anchor lies in its node: the node's bases before it, along the "
           "node's forward strand. A site between two nodes pins at a node end; a site inside one "
           "node, as an off-panel SNV is, pins on either side of its base. From v8\n";
    out << "#H\tA\tnode\tsnarl\tslot\tallele\tgqn\texplained\treliability\tpin\n";
    out << "#H\tR\tread_id\tstrand\toffset\tscore\n";
    for (size_t i = 0; i < names.size(); ++i) {
        out << "#read\t" << i << "\t" << names[i] << "\n";
    }
    // An anchor's A row and its R rows, formatted as the file has always written them.
    auto format_anchor = [&](ostream& os, const Anchor& a) {
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
        for (const ReadRow& r : a.reads) {
            os << "R\t" << table_index[r.read] << "\t" << (int)r.direction << "\t" << r.offset << "\t"
               << std::setprecision(1) << r.score << "\n";
        }
    };

    // The rows are formatted on several threads and written in order. Anchors are taken a batch at
    // a time, each batch cut into pieces that threads format into strings of their own, so that the
    // text held at once is one batch's rather than the whole file's: this runs when the run already
    // uses the most memory.
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
                format_anchor(os, all[i]);
            }
            pieces[p] = os.str();
        }
        for (const string& piece : pieces) {
            out.write(piece.data(), piece.size());
        }
    }
    return (bool)out;
}

void AnchorCollector::configure(const string& path, const AnchorParams& params,
                                const string& graph_name, const string& reads_source,
                                double mismap_min) {
    this->path = path;
    this->params = params;
    if (this->params.counters == nullptr) {
        // A caller that did not supply counters, such as a unit test, counts into this
        // instance's own, so that every use below can dereference without a check.
        this->params.counters = &this->owned_counters;
    }
    this->graph_name = graph_name;
    this->reads_source = reads_source;
    this->mismap_min = mismap_min;
    if (!path.empty()) {
        // One queue per OpenMP thread, since the passes that fill it are parallel.
        this->writer = make_unique<AnchorWriter>((size_t)max(1, get_thread_count()));
    }
}

void AnchorCollector::collect(const AnchorSiteEvidence& evidence, double explained_share,
                              const vector<int>& genotype, int haploid_slot,
                              const string& site_name, bool is_leaf, double gqn,
                              const multipass::ReadStrandTable& strands, size_t record_key) {
    if (path.empty() || writer == nullptr) {
        return;
    }
    if (params.leaf_only && !is_leaf) {
        return;
    }
    vector<AnchorWriter::Anchor> anchors;
    // The site's own log-odds per read, built on the first lookup and shared by the rest, so that
    // looking up every read at a site costs one pass over its reads rather than one per read.
    unordered_map<uint64_t, double> site_own;
    int site_own_state = 0; // 0: not built yet; 1: built; 2: nothing to leave out
    auto strand_log_odds = [&](std::string_view read_name) {
        if (site_own_state == 0) {
            site_own_state = strands.site_own_strand_log_odds(record_key, site_own) ? 1 : 2;
        }
        return strands.read_strand_log_odds(record_key, read_name,
                                            site_own_state == 1 ? &site_own : nullptr);
    };
    // Each read's strand log-odds, leaving out this site.
    vector<double> read_strand;
    // Built only where `build_site_anchors` reads it: at a diploid homozygote that may be split, or
    // at a heterozygous site with `phase_hets` or `strict_hets` on. The test must match its gate,
    // which also checks the vector's length against `evidence.reads`.
    const bool splittable_hom = genotype.size() == 2 && genotype[0] == genotype[1];
    const bool tiltable_het = (params.phase_hets || params.strict_hets)
                              && genotype.size() == 2
                              && genotype[0] != genotype[1];
    if ((params.hom_split && splittable_hom) || tiltable_het) {
        read_strand.reserve(evidence.reads.size());
        for (const AnchorRead& read : evidence.reads) {
            read_strand.push_back(strand_log_odds(read_names().name(read.read)));
        }
    }
    build_site_anchors(evidence, genotype, site_name, gqn, explained_share, haploid_slot, params,
                       *params.counters, anchors,
                       (params.hom_split || params.phase_hets || params.strict_hets)
                           ? &read_strand : nullptr);
    // A check for `hom_split`, reported per run: at heterozygous sites, whose alleles show which
    // strand each read is on, how often the read's strand log-odds agree. The log-odds leave the
    // site out. Computed only when splitting is on, and not when the heterozygous placement itself
    // uses the strand log-odds, since the check would then compare the strand with itself.
    if (params.hom_split && !params.phase_hets && !params.strict_hets
        && anchors.size() >= 2) {
        int slot_of_allele[2] = {-1, -1};
        int allele_of_slot[2] = {-1, -1};
        for (const AnchorWriter::Anchor& anchor : anchors) {
            if (anchor.slot >= 0 && anchor.slot < 2) {
                allele_of_slot[anchor.slot] = anchor.allele;
            }
        }
        if (allele_of_slot[0] >= 0 && allele_of_slot[1] >= 0
            && allele_of_slot[0] != allele_of_slot[1]) {
            (void)slot_of_allele;
            unordered_set<uint32_t> counted;
            for (const AnchorWriter::Anchor& anchor : anchors) {
                if (anchor.slot < 0 || anchor.slot > 1) {
                    continue;
                }
                for (const AnchorWriter::ReadRow& row : anchor.reads) {
                    if (!counted.insert(row.read).second) {
                        continue;   // both pins carry the same partition; count each read once
                    }
                    const double lo = strand_log_odds(read_names().name(row.read));
                    if (std::isnan(lo) || lo == 0.0) {
                        params.counters->phase_no_opinion.fetch_add(1);
                        continue;
                    }
                    const int phase_slot = lo > 0.0 ? 0 : 1;
                    const bool agree = phase_slot == anchor.slot;
                    params.counters->phase_checked.fetch_add(1);
                    if (agree) {
                        params.counters->phase_agree.fetch_add(1);
                    }
                    // The same threshold the split uses, `phase_min`, in natural-log units.
                    if (std::abs(lo) >= params.phase_min) {
                        params.counters->phase_confident.fetch_add(1);
                        if (agree) {
                            params.counters->phase_confident_agree.fetch_add(1);
                        }
                    }
                }
            }
        }
    }
    for (AnchorWriter::Anchor& anchor : anchors) {
        writer->add(std::move(anchor));
    }
}

void AnchorCollector::write(const string& sample_name) {
    if (path.empty() || writer == nullptr) {
        return;
    }
    size_t anchors = writer->anchor_count();
    size_t rows = writer->read_row_count();
    if (!writer->write(path, graph_name, sample_name, reads_source, mismap_min, params)) {
        return;
    }
    cerr << "[vg call] anchors: " << anchors << " written over " << rows
         << " read placements to " << path << endl;
    if (params.counters != nullptr) {
        params.counters->report(cerr);
    }
}

}
