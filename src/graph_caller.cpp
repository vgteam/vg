#include <atomic>
#include <chrono>

#include <omp.h>

#include "graph_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"

//#define debug

namespace vg {

/// Stage 0 instrumentation for post-linkage nested descent (eval task #52). Two things the current
/// design cannot answer about itself:
///
/// - How deep symbolic descent actually goes. A design that resolves linkage between levels needs
///   one barrier per level, and the level count has only ever been assumed.
/// - How many children descent skips because the *pre-linkage* parent genotype crosses them zero
///   times. Those are dropped at `copies <= 0` and never reconsidered, so a parent linkage moves
///   onto an allele that does cross the chain leaves a call nobody makes. The post-linkage half of
///   that count is in `LinkageCollector::resolve`, which is where the final genotype lives.
///
/// Namespace-scope, so they are zero-initialised before any dynamic initialisation and no
/// constructor has to know about them. Reported under --progress and otherwise inert.
static std::atomic<size_t> g_descent_depth_hist[16];
static std::atomic<size_t> g_descent_skipped_no_ref(0);
static std::atomic<size_t> g_descent_skipped_no_copy(0);
/// Children a called traversal enters more than once. Visits after the first are masked: one copy for
/// ploidy, and the first crossing for distance. A chain crossed twice by ONE traversal is one
/// haplotype carrying two copies, not two haplotypes carrying one each, so it must not become ploidy
/// 2 -- and representing the second copy at all is stage 17's question, deliberately deferred.
static std::atomic<size_t> g_child_multi_crossing(0);

// Stage 2 of planning/symbolic-diff-decomposition.md: measure, change nothing. Every counter here
// answers a question the plan currently answers with an offline Python proxy over INFO/AT, and the
// proxy cannot see what the caller sees -- notably that projection is inert for a flipped snarl.
static std::atomic<size_t> g_atomize_records(0);          // records with a symbolic site and a called ALT
static std::atomic<size_t> g_atomize_site_unresolvable(0); // flip_snarl left projection with no symbols
static std::atomic<size_t> g_atomize_site_reversed(0);     // resolved only via the reversed pairing
static std::atomic<size_t> g_atomize_blocks_hist[16];      // blocks per called ALT, index 15 = 15+
static std::atomic<size_t> g_atomize_blocks_total(0);
static std::atomic<size_t> g_atomize_multi_block(0);       // called ALTs a diff splits in two or more
static std::atomic<size_t> g_atomize_degraded(0);          // DP refused, fell back to one block
static std::atomic<size_t> g_atomize_case_c(0);            // chain in both alleles, unmatched: the only true double-report
static std::atomic<size_t> g_atomize_offref_alts(0);       // called ALTs carrying a chain the reference does not cross
static std::atomic<size_t> g_atomize_offref_bases(0);      // bases those chains account for
static std::atomic<size_t> g_atomize_alt_bases(0);         // total bases of the ALTs in the line above
static std::atomic<size_t> g_atomize_small_block_on_sv(0); // a >=50 bp ALT whose diff has a <50 bp block
// Stage 5: what the emitter actually did, as opposed to what stage 2 says it could do.
static std::atomic<size_t> g_atomize_split_sites(0);
static std::atomic<size_t> g_atomize_split_lines(0);
// Stage 6: chains whose own record is suppressed because a block ALT already spells them.
static std::atomic<size_t> g_atomize_child_inlined(0);
static thread_local int g_descent_depth = 0;

void GraphCaller::report_descent_instrumentation() const {
    size_t total = 0;
    for (int d = 0; d < 16; ++d) {
        total += g_descent_depth_hist[d].load();
    }
    if (total == 0) {
        return;   // no symbolic descent in this run
    }
    cerr << "[vg call] descent depth:";
    for (int d = 1; d < 16; ++d) {
        size_t n = g_descent_depth_hist[d].load();
        if (n > 0) {
            cerr << " " << d << "=" << n;
        }
    }
    cerr << " (" << total << " child calls)" << endl;
    if (g_child_multi_crossing.load() > 0) {
        cerr << "[vg call] descent: " << g_child_multi_crossing.load()
             << " children a called traversal enters more than once; visits after the first are"
             << " masked, so each contributes one copy and its first crossing's distance" << endl;
    }
    cerr << "[vg call] descent skipped: " << g_descent_skipped_no_copy.load()
         << " children no called allele reaches, " << g_descent_skipped_no_ref.load()
         << " with no reference path through them" << endl;
}

void report_atomize_instrumentation() {
    size_t records = g_atomize_records.load();
    size_t unresolvable = g_atomize_site_unresolvable.load();
    if (records == 0 && unresolvable == 0) {
        return;
    }

    cerr << "[vg call] atomize: " << records << " called ALTs projected, "
         << g_atomize_blocks_total.load() << " difference blocks, "
         << g_atomize_multi_block.load() << " ALTs a diff would split (>=2 blocks)" << endl;
    cerr << "[vg call] atomize blocks/ALT:";
    for (int b = 0; b < 16; ++b) {
        size_t n = g_atomize_blocks_hist[b].load();
        if (n > 0) {
            cerr << " " << b << (b == 15 ? "+" : "") << "=" << n;
        }
    }
    cerr << endl;

    // `unresolvable` is expected to be zero for the ordinary path, where every site is a managed
    // snarl and a reversed one now resolves through the reversed boundary pairing. It is NOT
    // expected to be zero under -I/--chains, which builds a fake snarl spanning a whole chain
    // (see the chain-piece construction in this file): that is not a managed snarl, so resolving it
    // to null is correct rather than a failure of the fix. The second number says how often the
    // reversed branch ran, which is what makes its coverage a measurement rather than an assumption.
    cerr << "[vg call] atomize: " << unresolvable
         << " sites where projection is inert because the snarl does not resolve, "
         << g_atomize_site_reversed.load()
         << " resolved as the reversal flip_snarl produces" << endl;

    cerr << "[vg call] atomize: " << g_atomize_case_c.load()
         << " ALTs carrying a chain the reference also crosses but the alignment did not match"
         << " -- the only population that would be reported twice" << endl;

    size_t offref = g_atomize_offref_alts.load();
    size_t offref_bases = g_atomize_offref_bases.load();
    size_t alt_bases = g_atomize_alt_bases.load();
    cerr << "[vg call] atomize: " << offref
         << " ALTs carry a chain the reference does not cross, " << offref_bases << " of "
         << alt_bases << " bases";
    if (alt_bases > 0) {
        cerr << " (" << (100.0 * (double)offref_bases / (double)alt_bases) << "%)";
    }
    cerr << " -- variation with no record of its own, reachable only inside this allele" << endl;

    if (g_atomize_child_inlined.load() > 0) {
        cerr << "[vg call] atomize: " << g_atomize_child_inlined.load()
             << " child chains not descended into because a block ALT already spells them" << endl;
    }
    if (g_atomize_split_sites.load() > 0) {
        cerr << "[vg call] atomize: " << g_atomize_split_sites.load()
             << " sites emitted as blocks instead of one record, " << g_atomize_split_lines.load()
             << " lines" << endl;
    }
    if (g_atomize_small_block_on_sv.load() > 0 || g_atomize_degraded.load() > 0) {
        cerr << "[vg call] atomize: " << g_atomize_small_block_on_sv.load()
             << " ALTs >=50 bp whose diff contains a block under 50 bp (truvari size-filter"
             << " exposure), " << g_atomize_degraded.load() << " alignments refused as too large"
             << endl;
    }
}


GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager), show_progress(false) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::set_show_progress(bool show_progress) {
    this->show_progress = show_progress;
}

void GraphCaller::set_node_id_ordering(bool ordered, size_t window_size) {
    node_id_ordering = ordered;
    node_id_window = max<size_t>(1, window_size);
}

/// Key a snarl by the lower of its two boundary node IDs. Sorting on this groups
/// snarls that a node-ID-range read source would fetch together.
static nid_t snarl_node_key(const Snarl* snarl) {
    return min(snarl->start().node_id(), snarl->end().node_id());
}

/// Sort snarls by node ID, so a read source fetching by node-ID range sees each
/// window once instead of re-querying per site.
static void sort_snarls_by_node_id(vector<const Snarl*>& snarls) {
    std::sort(snarls.begin(), snarls.end(), [](const Snarl* a, const Snarl* b) {
        return snarl_node_key(a) < snarl_node_key(b);
    });
}

void GraphCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {

    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<const Snarl*>> snarl_queue(thread_count);

    std::atomic<std::int64_t> top_snarl_count(0);
    std::atomic<std::int64_t> nested_snarl_count(0);
    bool top_level = true;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {

        if (!snarl_manager.is_trivial(snarl, graph)) {

#ifdef debug
            cerr << "GraphCaller running call_snarl on " << pb2json(*snarl) << endl;
#endif

            bool was_called = call_snarl(*snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
                vector<const Snarl*>& thread_queue = snarl_queue[omp_get_thread_num()];
                thread_queue.insert(thread_queue.end(), children.begin(), children.end());
            }
            
            if (show_progress) {
                if (top_level) {
                    ++top_snarl_count;
                    if (top_snarl_count % 100000 == 0) {
#pragma omp critical (cerr)
                        cerr << "[vg call]: Processed " << top_snarl_count << " top-level snarls" << endl;
                    }
                } else {
                    ++nested_snarl_count;
                    if (nested_snarl_count % 100000 == 0) {
#pragma omp critical (cerr)                    
                        cerr << "[vg call]: Processed " << top_snarl_count << " nested snarls" << endl;
                    }
                }
            }
        }
    };

    // Start with the top level snarls
    if (node_id_ordering) {
        // Visit in node-ID order, grouped into windows, so a read source that fetches
        // by node-ID range touches each window once and can release it. roots is
        // already a materialised vector, so this is a sort rather than a traversal.
        vector<const Snarl*> roots;
        snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            roots.push_back(snarl);
        });
        sort_snarls_by_node_id(roots);

        // Partition into contiguous windows. Parallelism becomes one task per window
        // rather than per snarl, which is what keeps a window's reads useful for the
        // whole time they are resident.
        vector<pair<size_t, size_t>> windows;
        size_t begin = 0;
        while (begin < roots.size()) {
            size_t window = (size_t)(snarl_node_key(roots[begin]) / (nid_t)node_id_window);
            size_t end = begin + 1;
            while (end < roots.size() &&
                   (size_t)(snarl_node_key(roots[end]) / (nid_t)node_id_window) == window) {
                ++end;
            }
            windows.emplace_back(begin, end);
            begin = end;
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int w = 0; w < (int)windows.size(); ++w) {
            for (size_t i = windows[w].first; i < windows[w].second; ++i) {
                process_snarl(roots[i]);
            }
        }
    } else {
        snarl_manager.for_each_top_level_snarl_parallel(process_snarl);
    }
    if (show_progress) cerr << "[vg call]: Finished processing " << top_snarl_count << " top-level snarls" << endl;

    top_level = false;

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(snarl_queue.begin(), snarl_queue.end(),
                        [](const vector<const Snarl*>& snarl_vec) {return snarl_vec.empty();})) {
        vector<const Snarl*> cur_queue;
        for (vector<const Snarl*>& thread_queue : snarl_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

        if (node_id_ordering) {
            // Keep queued children window-ordered as well, or the recursion rounds
            // undo the ordering the top-level pass established.
            sort_snarls_by_node_id(cur_queue);
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
    if (show_progress && nested_snarl_count > 0) cerr << "[vg call]: Finished processing " << nested_snarl_count << " nested snarls" << endl;
    if (show_progress) {
        report_descent_instrumentation();
    }
  
}

static void flip_snarl(Snarl& snarl) {
    Visit v = snarl.start();
    *snarl.mutable_start() = reverse(snarl.end());
    *snarl.mutable_end() = reverse(v);
}

void GraphCaller::call_top_level_chains(const HandleGraph& graph, size_t max_edges, size_t max_trivial, RecurseType recurse_type) {
    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<Chain>> chain_queue(thread_count);

    // Run the snarl caller on a chain. queue up the children if it fails
    auto process_chain = [&](const Chain* chain) {

#ifdef debug
        cerr << "calling top level chain ";
        for (const auto& i : *chain) {
            cerr << pb2json(*i.first) << "," << i.second << ",";
        }
        cerr << endl;
#endif
        // Break up the chain
        vector<Chain> chain_pieces = break_chain(graph, *chain, max_edges, max_trivial);

        for (Chain& chain_piece : chain_pieces) {
            // Make a fake snarl spanning the chain
            // It is important to remember that along with not actually being a snarl,
            // it's not managed by the snarl manager so functions looking into its nesting
            // structure will not work
            Snarl fake_snarl;
            *fake_snarl.mutable_start() = chain_piece.front().second == true ? reverse(chain_piece.front().first->end()) :
                chain_piece.front().first->start();
            *fake_snarl.mutable_end() = chain_piece.back().second == true ? reverse(chain_piece.back().first->start()) :
                chain_piece.back().first->end();

#ifdef debug
            cerr << "calling fake snarl " << pb2json(fake_snarl) << endl;
#endif
            
            bool was_called = call_snarl(fake_snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                vector<Chain>& thread_queue = chain_queue[omp_get_thread_num()];                
                for (pair<const Snarl*, bool> chain_link : chain_piece) {
                    const deque<Chain>& child_chains = snarl_manager.chains_of(chain_link.first);
                    thread_queue.insert(thread_queue.end(), child_chains.begin(), child_chains.end());
                }
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_chain_parallel(process_chain);

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(chain_queue.begin(), chain_queue.end(),
                        [](const vector<Chain>& chain_vec) {return chain_vec.empty();})) {
        vector<Chain> cur_queue;
        for (vector<Chain>& thread_queue : chain_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_chain(&cur_queue[i]);
        }
    
    }
}

vector<Chain> GraphCaller::break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial) {
    
    vector<Chain> chain_frags;

    // keep track of the current fragment and add it to chain_frags as soon as it gets too big
    Chain frag;
    size_t frag_edge_count = 0;
    size_t frag_triv_count = 0;
    
    for (const pair<const Snarl*, bool>& link : chain) {
        // todo: we're getting the contents here as well as within the caller.
        auto contents = snarl_manager.deep_contents(link.first, graph, false);

        // todo: use annotation from snarl itself?
        bool trivial = contents.second.empty();

        if ((trivial && frag_triv_count > max_trivial) ||
            (contents.second.size() + frag_edge_count > max_edges)) {
            // adding anything more to the chain would make it too long, so we
            // add it to the output and clear the current fragment
            if (!frag.empty() && frag_triv_count < frag.size()) {
                chain_frags.push_back(frag);
            }
            frag.clear();
            frag_edge_count = 0;
            frag_triv_count = 0;
        }

        if (!trivial || (frag_triv_count < max_trivial)) {
            // we start a new fragment or add to an existing fragment
            frag.push_back(link);
            frag_edge_count += contents.second.size();
            if (trivial) {
                ++frag_triv_count;
            }
        }
    }

    // and the last one
    if (!frag.empty()) {
        chain_frags.push_back(frag);
    }

    return chain_frags;
}
    
VCFOutputCaller::VCFOutputCaller(const string& sample_name) : sample_name(sample_name), translation(nullptr), include_nested(false)
{
    output_variants.resize(get_thread_count());
}

VCFOutputCaller::~VCFOutputCaller() {
}

string VCFOutputCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;    
    for (int i = 0; i < contigs.size(); ++i) {
        const string& contig = contigs[i];
        size_t length;
        if (i < contig_length_overrides.size()) {
            // length override provided
            length = contig_length_overrides[i];
        } else {
            length = 0;
            for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
                length += graph.get_length(handle);
            }
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    if (include_nested) {
        ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree counting only ancestors whose record is on this record's own reference contig (0=top level for this contig)\">" << endl;
        ss << "##INFO=<ID=CH,Number=1,Type=Integer,Description=\"Number of ancestors in the VCF whose record is on a different reference contig than the one below it, ie nesting steps between VCF reference contigs\">" << endl;
        ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
        ss << "##INFO=<ID=RC,Number=1,Type=String,Description=\"Reference contig of top-level containing site\">" << endl;
        ss << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Reference start position of top-level containing site\">" << endl;
        ss << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"Reference end position of top-level containing site\">" << endl;
    }
    if (emit_phasing) {
        // FORMAT/PS is the VCF standard phase set and is what phasing tools look for. Note the
        // deliberate name clash with INFO/PS above, which is vg's parent-snarl pointer under -A:
        // different namespaces, so both are legal in one file, but a reader skimming for "PS"
        // will find two unrelated things and the descriptions have to say which is which.
        ss << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: the phase of a "
           << "genotype is comparable only with others carrying the same PS. One phase set per "
           << "chain, so blocks are chromosome-scale -- much longer than a read-based phaser "
           << "gives, because the phase comes from the haplotype panel rather than from reads "
           << "spanning consecutive sites. Not the INFO/PS emitted under -A, which is a parent "
           << "snarl pointer\">" << endl;
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    if (atomize_blocks) {
        ss << "##INFO=<ID=SB,Number=2,Type=Integer,Description=\"Index and count of this "
           << "difference block within its snarl: this record is one of several the same snarl "
           << "emitted, because the reference and the called haplotypes differ from each other in "
           << "more than one place inside it. All records of one snarl share an ID. "
           << "DOUBLE COUNTING: the per-sample evidence is the SNARL's, repeated on every block, "
           << "not apportioned between them -- AD, GL, GQ, GQI, GP and QUAL are identical across "
           << "the set, because the genotype likelihood was computed over whole-snarl traversals "
           << "and has no per-block decomposition. DP, DR and BL are per-site read counts and are "
           << "site-level by definition. So any consumer that sums, averages or otherwise "
           << "aggregates evidence across records must group by ID first and count each snarl "
           << "once. Records without SB are unaffected: they are the only record their snarl "
           << "emitted.\">" << endl;
    }
    if (allele_merge_threshold < 1.0) {
        ss << "##INFO=<ID=MAT,Number=.,Type=String,Description=\"Merged Allele Traversal: "
           << "ALT alleles merged after genotyping by -L/--cluster, as OLD>NEW:SIMILARITY using "
           << "pre-merge allele numbers. In a nested run this record gives the collapsed view of "
           << "the site and its child records the precise one, so they disagree by design. "
           << "pre-merge allele numbers. AD and GL are folded onto the surviving allele and MAD is "
           << "recomputed; DP, QUAL, GQ, GP and FILTER are as computed over the pre-merge allele set.\">"
           << endl;
    }
    return ss.str();
}

void VCFOutputCaller::set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                                  const vector<size_t>* sequence_to_haplotype) {
    this->linkage_collector = collector;
    this->linkage_gbwt = gbwt;
    this->linkage_sequence_to_haplotype = sequence_to_haplotype;
    this->linkage_gbwt_cache.clear();
    this->linkage_gbwt_cache_anchor.clear();
    if (gbwt != nullptr) {
        // One per thread, built here so the parallel region never allocates one.
        this->linkage_gbwt_cache.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            this->linkage_gbwt_cache.emplace_back(*gbwt);
        }
        this->linkage_gbwt_cache_anchor.assign(omp_get_max_threads(), 0);
    }
}

vector<int> VCFOutputCaller::panel_alleles(const HandleGraph& graph,
                                          const vector<SnarlTraversal>& travs) const {
    vector<int> out;
    if (linkage_gbwt == nullptr || linkage_sequence_to_haplotype == nullptr) {
        return out;
    }
    // -1 means "this haplotype carries no allele here", which is not the same as carrying the
    // reference: a haplotype whose path ends inside the site genuinely has nothing to say, and
    // recording it as reference would invent evidence.
    out.assign(linkage_sequence_to_haplotype->size(), -1);

    // The cache, not the index: same results, but records stay decompressed between sites.
    // Falls back to the index itself if set_linkage was never given one to size the vector.
    int thread = omp_get_thread_num();
    const bool cached = (size_t)thread < linkage_gbwt_cache.size();

    // CachedGBWT only grows -- the gbwt header recommends short-lived instances -- and with
    // node-ID-ordered windows a thread never revisits a retired window, so an uncleared cache
    // accumulates every decompressed record the contig ever touched. Clearing when the site's
    // neighbourhood moves past the fetch-window width keeps the hit rate the cache exists for
    // (adjacent snarls share records) while bounding residency to about one window.
    //
    // Measured on chr20, isolated runs of the same binary: 4.17 GB peak without this, 3.83 GB with
    // it, and the two outputs agree as multisets -- the cache is an accelerator, so anything else
    // would be a bug. Runtime is unchanged (179.3 s against 175.8 s, inside run-to-run noise).
    if (cached && (size_t)thread < linkage_gbwt_cache_anchor.size() && !travs.empty()) {
        static const nid_t CACHE_ANCHOR_SPAN = 4096;
        nid_t lead = 0;
        for (int64_t i = 0; i < travs[0].visit_size() && lead == 0; ++i) {
            lead = travs[0].visit(i).node_id();
        }
        if (lead != 0) {
            nid_t& anchor = linkage_gbwt_cache_anchor[thread];
            if (anchor == 0 || lead > anchor + CACHE_ANCHOR_SPAN
                || lead + CACHE_ANCHOR_SPAN < anchor) {
                linkage_gbwt_cache[thread].clearCache();
                anchor = lead;
            }
        }
    }

    for (size_t a = 0; a < travs.size(); ++a) {
        const SnarlTraversal& trav = travs[a];
        if (trav.visit_size() < 1) {
            continue;
        }
        gbwt::SearchState state;
        bool ok = true;
        for (int64_t i = 0; i < trav.visit_size(); ++i) {
            const Visit& visit = trav.visit(i);
            if (visit.node_id() == 0) {
                // A visit to a child snarl rather than a node: the traversal is not expanded, so
                // it cannot be looked up in the GBWT.
                ok = false;
                break;
            }
            gbwt::node_type node = gbwt::Node::encode(visit.node_id(), visit.backward());
            if (cached) {
                const gbwt::CachedGBWT& c = linkage_gbwt_cache[thread];
                state = (i == 0) ? c.find(node) : c.extend(state, node);
            } else {
                state = (i == 0) ? linkage_gbwt->find(node) : linkage_gbwt->extend(state, node);
            }
            if (state.empty()) {
                ok = false;
                break;
            }
        }
        if (!ok || state.empty()) {
            continue;
        }
        vector<gbwt::size_type> seqs = cached ? linkage_gbwt_cache[thread].locate(state)
                                              : linkage_gbwt->locate(state);
        for (gbwt::size_type seq : seqs) {
            if (seq < linkage_sequence_to_haplotype->size()) {
                size_t hap = (*linkage_sequence_to_haplotype)[seq];
                if (hap < out.size()) {
                    // A haplotype in several fragments can reach the same site twice; the allele
                    // is the same either way, so first writer wins rather than being an error.
                    out[hap] = (int)a;
                }
            }
        }
    }
    return out;
}

void VCFOutputCaller::set_ploidy_regions(const string& bed_path) {
    ifstream in(bed_path);
    if (!in) {
        cerr << "error [vg call]: could not open --ploidy-bed file " << bed_path << endl;
        exit(1);
    }
    string line;
    size_t line_number = 0;
    while (getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0
            || line.compare(0, 7, "browser") == 0) {
            continue;
        }
        istringstream ss(line);
        string chrom;
        long long start = -1, end = -1;
        int ploidy = -1;
        if (!(ss >> chrom >> start >> end >> ploidy)) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " is not CHROM START END PLOIDY: " << line << endl;
            exit(1);
        }
        if (start < 0 || end < start) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " has a negative or reversed interval: " << line << endl;
            exit(1);
        }
        // The callers implement ploidy 1 and 2 only. Caught here rather than left to reach a
        // caller as an unsupported ploidy, which aborts in a much less informative place.
        if (ploidy != 1 && ploidy != 2) {
            cerr << "error [vg call]: --ploidy-bed " << bed_path << " line " << line_number
                 << " has ploidy " << ploidy << ", which must be 1 or 2" << endl;
            exit(1);
        }
        if (start == end) {
            // Covers nothing. Keeping it would leave a region no lookup can ever hit.
            continue;
        }
        ploidy_regions[chrom].push_back({(size_t)start, (size_t)end, ploidy});
    }

    // Sorted so lookups can binary-search, and checked for overlap while they are in order.
    for (auto& entry : ploidy_regions) {
        auto& regions = entry.second;
        sort(regions.begin(), regions.end(),
             [](const PloidyRegion& a, const PloidyRegion& b) { return a.start < b.start; });
        for (size_t i = 1; i < regions.size(); ++i) {
            if (regions[i].start < regions[i - 1].end) {
                cerr << "error [vg call]: --ploidy-bed " << bed_path << " has overlapping "
                     << "intervals on " << entry.first << ": [" << regions[i - 1].start << ","
                     << regions[i - 1].end << ") and [" << regions[i].start << ","
                     << regions[i].end << "). Two ploidies for one base has no correct reading, "
                     << "so this is not resolved by precedence." << endl;
                exit(1);
            }
        }
    }
}

int VCFOutputCaller::region_ploidy(const string& ref_path_name, size_t position,
                                   int fallback) const {
    if (ploidy_regions.empty()) {
        return fallback;
    }
    // Match on the contig as the VCF spells it, so a BED written against the output works.
    // Same reduction emit_variant applies when it sets sequenceName.
    string contig = Paths::strip_subrange(ref_path_name);
    string locus = PathMetadata::parse_locus_name(contig);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        contig = locus;
    }
    auto found = ploidy_regions.find(contig);
    if (found == ploidy_regions.end()) {
        return fallback;
    }
    const vector<PloidyRegion>& regions = found->second;
    // First region starting after the position; its predecessor is the only one that can cover,
    // since the regions are non-overlapping.
    auto it = upper_bound(regions.begin(), regions.end(), position,
                          [](size_t p, const PloidyRegion& r) { return p < r.start; });
    if (it == regions.begin()) {
        return fallback;
    }
    --it;
    return (position >= it->start && position < it->end) ? it->ploidy : fallback;
}

int VCFOutputCaller::ploidy_at(const string& ref_path_name, int64_t interval_start,
                               int64_t ref_offset, int fallback) const {
    if (ploidy_regions.empty()) {
        return fallback;
    }
    // Same arithmetic emit_variant uses for POS, minus the +1 that makes VCF 1-based: the BED is
    // 0-based, so the two agree on which base an interval boundary falls on.
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    int64_t position = interval_start + ref_offset + basepath_offset;
    if (position < 0) {
        return fallback;
    }
    return region_ploidy(ref_path_name, (size_t)position, fallback);
}

size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout) {
    if (layout == GLLayout::Colexicographic) {
        // The VCF spec's order: genotypes sorted by their larger allele, then their smaller. Does
        // not depend on the allele count at all, which is why it is the one a reader can trust.
        return j * (j + 1) / 2 + i;
    }
    // i-major: all genotypes with smaller allele 0, then all with 1, and so on.
    return i * n_alleles - (i * (i - 1)) / 2 + (j - i);
}

vector<double> fold_genotype_likelihoods(const vector<double>& old_gl,
                                         const vector<int>& new_index,
                                         size_t n_new, GLLayout layout) {
    const size_t n_old = new_index.size();
    vector<double> folded(n_new * (n_new + 1) / 2, -std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < n_old; ++i) {
        for (size_t j = i; j < n_old; ++j) {
            const size_t from = gl_genotype_index(i, j, n_old, layout);
            if (from >= old_gl.size()) {
                continue;
            }
            size_t ni = (size_t)new_index[i], nj = (size_t)new_index[j];
            if (ni > nj) {
                std::swap(ni, nj);
            }
            double& slot = folded[gl_genotype_index(ni, nj, n_new, layout)];
            slot = std::max(slot, old_gl[from]);
        }
    }
    return folded;
}

bool buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b) {
    if (a.contig != b.contig) {
        return a.contig < b.contig;
    }
    if (a.position != b.position) {
        return a.position < b.position;
    }
    if (a.id != b.id) {
        return a.id < b.id;
    }
    return a.block < b.block;
}

bool VCFOutputCaller::add_variant(vcflib::Variant& var, size_t block) const {
    var.setVariantCallFile(output_vcf);
    stringstream ss;
    ss << var;
    string dest;
    if (ss.str().length() > VCFOutputCaller::max_vcf_line_length) {
        return false;
    }         
    int ret = zstdutil::CompressString(ss.str(), dest);
    assert(ret == 0);
    // the Variant object is too big to keep in memory when there are many genotypes, so we
    // store it in a zstd-compressed string
    output_variants[omp_get_thread_num()].push_back(
        make_pair(BufferedRecordKey{var.sequenceName, (size_t)var.position, var.id, block}, dest));
    return true;
}

void VCFOutputCaller::resolve_linkage() {
    if (linkage_resolved) {
        return;
    }
    if (linkage_collector == nullptr) {
        resolve_linkage_generation(0, true);
        return;
    }
    // Every generation, not just the first. Chain construction skips entries whose generation is
    // above the one being resolved, so resolving 0 alone dropped every nested site from linkage,
    // from phasing and from the mosaic. That is latent rather than live -- `call_main` arms deferred
    // descent for every nested run, and `run_deferred_descent` does loop the generations, leaving
    // this a no-op via `linkage_resolved` -- but the invariant the header states is that a nested
    // site is always settled, and it should not depend on which caller got there first.
    //
    // `max_generation()` is re-read each pass for the same reason the barrier re-reads it: a pass
    // can gain a chain at a deeper generation than anything recorded before it, and a bound
    // snapshotted once would leave that chain emitted but never settled.
    for (size_t gen = 0;; ++gen) {
        const size_t deepest = linkage_collector->max_generation();
        resolve_linkage_generation(gen, gen >= deepest);
        if (gen >= deepest) {
            break;
        }
    }
}

size_t VCFOutputCaller::record_key_of(const Snarl& snarl) const {
    return std::hash<string>{}(print_snarl(snarl, false));
}

void VCFOutputCaller::build_render_phases() {
    // Built between the barrier and the render, from the phasing the barrier accumulated.
    //
    // No `emitted` filter here, unlike the mosaic and unlike the patch index this replaces. That
    // filter existed to keep sites with no VCF line out of a structure keyed to lines, and it is what
    // forced the whole bookkeeping to run *after* the render -- because "does this site have a line?"
    // is false for everything while the genotypes are still being decided. A render-time lookup needs
    // no such filter: a site with no line simply never looks itself up.
    render_phases.clear();
    if (!emit_phasing) {
        return;
    }
    render_phases.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        // Last writer wins, which is the later generation: a site revised at the barrier gets a
        // second PhaseCall and the later one describes the genotype it ends up with.
        render_phases[pc.record_key] = pc;
    }
}

void VCFOutputCaller::finalise_linkage_outputs() {
    // Built here, after every record has been rendered, and not while the genotypes are being
    // resolved.
    //
    // Both of these read `PhaseCall::emitted` to tell a record from a site that wrote no line,
    // and once the record is built AFTER the decision that is not known during resolution: every
    // entry still says unemitted there. Building the patch index then produced an empty map and
    // unphased the whole output; building the mosaic then would have filled it with the 100k
    // sites that never become records.
    if (linkage_collector == nullptr) {
        return;
    }
    // Live, not the snapshot each PhaseCall carries: `linkage_phased` holds copies taken during
    // resolution, before any line existed, so every copy's `emitted` is false.
    const std::unordered_set<size_t> emitted_records = linkage_collector->emitted_records();
    size_t unexplained = 0;
    size_t order_arbitrary = 0;
    // The patch index is gone: phasing is applied while each record is rendered, so there is nothing
    // here to key to a line. What survives is the mosaic and the counters, both of which genuinely
    // need to know which sites became records -- read live from `emitted_records`, because a
    // PhaseCall's own `emitted` is a snapshot taken before any record existed.
    size_t phased_unwritten = 0;
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        if (emitted_records.count(pc.record_key) == 0) {
            // Phased, and deliberately so -- its children read their strand off it -- but it is not
            // a record, so it must not enter the mosaic or any count of records. Counted on its own.
            ++phased_unwritten;
            continue;
        }
        // Count only the strands a chain actually has. On a haploid chain hap_second is the
        // wildcard by construction, so counting it reported every site as unexplained while
        // the mosaic was naming real haplotypes throughout.
        // A haploid site has one strand and one wildcard, and *which* side holds the
        // wildcard is not fixed: a haploid contig fills the first, while a nested site hanging
        // off the parent's second strand fills the second. Testing hap_first alone therefore
        // reported every nested site on the parent's second strand as unexplained when the
        // panel named its haplotype perfectly well -- 612 of chr20's 2020, which is why that
        // count fell to 1408 the moment the strands were derived correctly rather than because
        // any site became better explained.
        unexplained += (pc.ploidy == 1)
                       ? (pc.hap_first == LinkageModel::WILDCARD
                          && pc.hap_second == LinkageModel::WILDCARD)
                       : (pc.hap_first == LinkageModel::WILDCARD
                          || pc.hap_second == LinkageModel::WILDCARD);
        order_arbitrary += pc.order_arbitrary;
    }
    cerr << "[vg call] linkage: " << linkage_collector->num_sites() << " sites, "
         << (linkage_collector->bytes() / (1024.0 * 1024.0)) << " MB retained, "
         << linkage_changed << " genotypes moved by linkage, " << linkage_seconds << " s" << endl;
    if (emit_phasing) {
        // The wildcard count is the honest caveat on a chromosome-length phase block: at
        // those sites the panel does not name a strand, so the phase either side of them
        // rests on the transition model alone.
        cerr << "[vg call] phasing: " << (linkage_phased.size() - phased_unwritten)
             << " sites phased, " << unexplained
             << " with a strand the panel does not explain" << endl;
        if (phased_unwritten > 0) {
            // Sites that wrote no VCF line and are phased anyway. A parent whose alleles differ
            // only inside its children collapses to the reference and emits nothing, and its
            // children still need to know which of its two haplotypes carries the chain. This is
            // the population that used to be absent from the layer altogether.
            cerr << "[vg call] phasing: " << phased_unwritten
                 << " collapsed sites phased with no line of their own, so their children can"
                 << " inherit a strand" << endl;
        }
        if (order_arbitrary > 0) {
            // A heterozygous site where no panel haplotype on either strand spells either called
            // allele. The record still comes out phased and in the block, because it has a
            // position in it, but which allele went on which strand was decided by sorting the
            // pair -- so the orientation there is a placeholder, not a call. Reported because a
            // reader cannot tell these from the rest.
            cerr << "[vg call] phasing: " << order_arbitrary
                 << " heterozygous sites carry an allele order the panel does not determine"
                 << endl;
        }
    }
    if (!mosaic_path.empty()) {
        // Records only. The mosaic's segments are runs over *sites in the call set*, and its site
        // counts are index arithmetic over the vector it is handed, so a collapsed site with no
        // line would inflate every run it fell inside and break the invariant that the mosaic
        // accounts for exactly the emitted records. Tracing a path through those sites is what
        // stage 5 of the traversal-space plan is for, and it needs more than a row count.
        vector<LinkageCollector::PhaseCall> written;
        written.reserve(linkage_phased.size());
        for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
            if (emitted_records.count(pc.record_key) != 0) {
                written.push_back(pc);
            }
        }
        write_mosaic(written);
    }
}

void VCFOutputCaller::resolve_linkage_generation(size_t generation, bool last) {
    linkage_resolved = true;
    if (linkage_collector == nullptr) {
        return;
    }
    // Reported rather than estimated. The retained-bytes figure in the LinkageCollector
    // header comment was arithmetic -- sites times a per-site size -- and `bytes()` exists so
    // that it can be an observation instead; it had never been called. The elapsed time
    // answers the other question the design asserted without checking: this pass is serial,
    // between calling and writing, in a caller that is otherwise parallel over snarls.
    auto start = std::chrono::steady_clock::now();
    // `linkage_phased` accumulates across generations rather than being replaced. The model needs
    // the earlier generations back: a nested site's strand is read off its parent's PhaseCall, and
    // a clamped site's phase is pinned to the pair already emitted for it.
    const size_t moved =
        linkage_collector->resolve_generation(generation, last,
                                              emit_phasing ? &linkage_phased : nullptr);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    linkage_seconds += seconds;
    // How many sites the model moved off the genotype the reads alone chose. This used to be the
    // number of *patches* produced, which since the record is built from the settled genotype would
    // now read zero at every generation -- a counter saying linkage changed nothing, about a pass
    // that decides every genotype in the output.
    linkage_changed += moved;
    if (!last) {
        // One line per intermediate generation, so a deferred-descent run shows its own shape:
        // how many sites each barrier settled and what it cost.
        cerr << "[vg call] linkage generation " << generation << ": "
             << linkage_collector->num_sites_at(generation) << " sites, "
             << moved << " genotypes moved by linkage, " << seconds << " s" << endl;
        return;
    }

}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (include_nested) {
        update_nesting_info_tags(snarl_manager);
    }
    vector<pair<BufferedRecordKey, string>> all_variants;
    // Reserve once: doing it inside the loop below reallocates per thread buffer.
    size_t total_variants = 0;
    for (const auto& buf : output_variants) {
        total_variants += buf.size();
    }
    all_variants.reserve(total_variants);
    // `buf` must not be const: std::move() over const_iterators silently degrades to a copy,
    // which duplicated every compressed record at the one point where the whole VCF is in
    // memory at once.  Free the buffer as we go for the same reason.
    //
    // This makes write_variants() single-use, which it already effectively was -- a real move
    // leaves the buffers empty either way.  All three callers (deconstructor.cpp,
    // call_main.cpp, mcmc_main.cpp) call it exactly once.
    for (auto& buf : output_variants) {
        std::move(buf.begin(), buf.end(), std::back_inserter(all_variants));
        buf.clear();
        buf.shrink_to_fit();
    }
    std::sort(all_variants.begin(), all_variants.end(),
              [](const pair<BufferedRecordKey, string>& v1,
                 const pair<BufferedRecordKey, string>& v2) {
                  return buffered_record_key_less(v1.first, v2.first);
              });
    // Phase two of the linkage pass. The records are already all here, compressed, with
    // (contig, position) uncompressed as the sort key -- so a change can be matched without ever
    // having kept the record itself, and only the records that actually move are re-parsed.
    resolve_linkage();
    finalise_linkage_outputs();


    for (const auto& v : all_variants) {
        if (v.second.empty()) {
            // A tombstone: the barrier retracted or replaced this line in place, by the handle it
            // captured at emit time. Identifying stale lines here by re-hashing the ID column and
            // counting GT separators mistook a phased haploid replacement ("1|.") for a diploid
            // line, and was gated on sets that could be empty while replacements existed.
            continue;
        }
        string dest;
        int ret = zstdutil::DecompressString(v.second, dest);
        assert(ret == 0);
        // The record key is the hash of the ID column, which is how the linkage layer keyed the
        // site -- see `record_key_of`, which every producer goes through -- so the identity is
        // recoverable from the line itself and nothing extra has to be carried through the
        // compressed buffer. This is the producer that cannot be changed, so it is the one that
        // fixes the form for the other six. Computed once, lazily: several records can share
        // a (contig, position), and every patch below must land on its own record, not the first
        // line at the position.
        size_t line_key = 0;
        bool have_line_key = false;
        auto id_key = [&]() -> size_t {
            if (!have_line_key) {
                size_t a = dest.find('\t');
                size_t b = a == string::npos ? string::npos : dest.find('\t', a + 1);
                size_t c = b == string::npos ? string::npos : dest.find('\t', b + 1);
                if (c != string::npos) {
                    line_key = std::hash<string>{}(dest.substr(b + 1, c - b - 1));
                }
                have_line_key = true;
            }
            return line_key;
        };
        if (linkage_collector != nullptr) {
            // Quality before phasing, and no genotype patch between them any more: the genotype the
            // line carries IS the settled one, because the line was built from it.
            const auto& quality = linkage_collector->moved_quality();
            if (!quality.empty()) {
                auto found = quality.find(id_key());
                if (found != quality.end()) {
                    if (!apply_linkage_quality(dest, found->second.first, found->second.second)) {
                        ++quality_declined;
                    }
                }
            }
        }
        out_stream << dest << endl;
    }
    if (phase_declined.load() > 0 || quality_declined.load() > 0) {
        cerr << "[vg call] linkage: " << phase_declined.load()
             << " phases refused by the record they were rendered onto, and "
             << quality_declined.load() << " quality rewrites refused" << endl;
    }
    // Here rather than with the descent report: rendering happens after the calling sweep, so a
    // report printed there sees only the records that emitted inline.
    report_atomize_instrumentation();
}


gbwt::edge_type VCFOutputCaller::mosaic_gbwt_position(int64_t node_id, size_t hap) const {
    if (linkage_gbwt == nullptr || linkage_sequence_to_haplotype == nullptr) {
        return gbwt::invalid_edge();
    }
    // Forward first: snarl boundaries are stored oriented along the reference, so the reverse
    // orientation is the exception rather than a coin flip.
    for (int orientation = 0; orientation < 2; ++orientation) {
        gbwt::node_type node = gbwt::Node::encode(node_id, orientation == 1);
        gbwt::SearchState state = linkage_gbwt->find(node);
        if (state.empty()) {
            continue;
        }
        for (gbwt::size_type i = state.range.first; i <= state.range.second; ++i) {
            gbwt::size_type seq = linkage_gbwt->locate(node, i);
            if (seq < linkage_sequence_to_haplotype->size()
                && (*linkage_sequence_to_haplotype)[seq] == hap) {
                return gbwt::edge_type(node, i);
            }
        }
    }
    return gbwt::invalid_edge();
}

void VCFOutputCaller::write_mosaic(const vector<LinkageCollector::PhaseCall>& phasing) const {
    ofstream out(mosaic_path);
    if (!out) {
        cerr << "error [vg call]: could not open " << mosaic_path << " for the mosaic output"
             << endl;
        return;
    }

    // A segment is a maximal run over which one strand stays on one panel haplotype. The consumer
    // reconstructs a haplotype by walking it from start_node to end_node, so the anchors are node
    // IDs rather than reference positions: a node ID is intrinsic to the graph, while a position
    // is a statement about one reference path.
    //
    // Two things the header has to state outright, because a consumer cannot recover either from
    // the segment rows alone.
    //
    // **Which reference the positions are in.** The segment rows carry the contig as the VCF
    // spells it -- the locus part of the path name -- and a graph can hold several references
    // under the same locus. The HPRC graphs do: their GBZ tags name CHM13 and GRCh38 both as
    // reference samples and both appear in the panel, so `chr20` alone does not say whether
    // position 24 is CHM13 chr20 or GRCh38 chr20. Those are different coordinate systems and a
    // consumer guessing wrong would see nothing amiss.
    //
    // **What a hap_index means.** It is assigned in GBWT metadata order by the run that produced
    // the file and is meaningless outside it, so the whole mapping is written out. The name is
    // the portable identifier: a (sample, phase) pair, which is the unit the linkage model works
    // in -- a haplotype present in several GBWT fragments is one haplotype, so no single GBWT path
    // name would do. Name plus the segment's contig is enough to find the paths again.
    out << "#mosaic-version\t3\n";
    out << "#graph\t" << mosaic_graph_name << "\n";
    out << "#sample\t" << sample_name << "\n";
    for (const string& ref : mosaic_reference_paths) {
        out << "#reference\t" << ref << "\n";
    }
    out << "#decoding\tconstrained-viterbi\n";
    // Positions are derived; the node IDs are not. Said plainly so a consumer knows which to
    // trust when a file is read against a graph whose reference paths have moved.
    out << "#note\tref_start/ref_end are advisory, in the #reference coordinate system; "
        << "start_node/end_node are the authoritative anchors and are intrinsic to the graph.\n";
    out << "#note\tsegments are maximal runs on one panel haplotype; walk the haplotype "
        << "from start_node to end_node to reconstruct it. * means the panel does not "
        << "explain that strand there; . means that strand carries no sequence there at all, "
        << "which is what a nested haploid site's other strand looks like. Version 2 spelled "
        << "both with *.\n";
    out << "#note\thap_index is internal to this run; haplotype (sample#phase) is the portable "
        << "identifier and names a haplotype, not a single GBWT path.\n";
    for (size_t h = 0; h < mosaic_haplotype_names.size(); ++h) {
        out << "#haplotype\t" << h << "\t" << mosaic_haplotype_names[h] << "\n";
    }
    out << "#note\tgbwt_node/gbwt_offset is the GBWT position of the haplotype at start_node: "
        << "extract({gbwt_node, gbwt_offset}) and follow LF() to end_node, with no locate and no "
        << "r-index. gbwt_node carries the orientation, start_node does not.\n";
    out << "#note\ta segment never spans a GBWT fragment boundary, so one position walks the "
        << "whole of it; a haplotype in several fragments yields several segments.\n";
    out << "#H\tcontig\tstrand\tref_start\tref_end\tstart_node\tend_node"
        << "\thap_index\thaplotype\tsites\tgbwt_node\tgbwt_offset\n";

    // The phasing arrives grouped by contig and in reference order, which is how resolve() builds
    // it. Both strands are emitted, and a switch on either one closes only its own segment.
    size_t i = 0;
    size_t total_segments = 0;

    // Emit sites [from, to] on one strand, all on haplotype `hap`, as one row per GBWT fragment.
    //
    // A row carries one GBWT position, which is a claim that the position walks the whole segment.
    // That fails across a fragment boundary, so a run is cut wherever the fragment under it changes.
    // Finding those cuts is the entire cost of the feature, and how it is found matters a great deal
    // -- measured on chr20 against a 146 s baseline with positions written but no splitting:
    //
    //   * Resolving a position at every site: 332 s. Seven million resolves.
    //   * Following the haplotype with LF from site to site sounds far cheaper and is not: 172 s,
    //     210 million LF steps. Sites sit about twenty nodes apart in reference terms, but where a
    //     haplotype runs in the reverse orientation, walking forward moves *away* from the next site
    //     and burns the whole step budget before giving up. 6,640 of 210,000 transitions did.
    //   * Resolving only the two *ends* of a run and comparing them: 150 s. Further work is needed
    //     only where they differ, and then a binary search over the sites finds the boundary in log
    //     time. On chr20 that is 7,344 resolves and two searches.
    //
    // The third is what this does. The first two are recorded because both look obviously cheaper
    // than they are.
    //
    // The limitation this accepts: comparing endpoints detects a fragment that *changes* across a
    // run, not one that leaves and returns within it. A haplotype re-entering its starting fragment
    // before the run ends would be emitted as one row, and the position would not walk the middle of
    // it. Fragments partition a path, so this needs the run to leave and re-enter the same fragment
    // in reference order; it does not arise on this panel, and detecting it would cost the 2.3x
    // above.
    // What a strand actually holds at a site. The file used to spell two different facts with one
    // character: a nested haploid site's *other* strand has no sequence there at all, while an
    // unexplained strand has sequence the panel cannot attribute to a haplotype. Conflating them is
    // why the wildcard count is not a usable metric -- raw wildcard segments rose 437 -> 616 across
    // the traversal-space work while the count that means only the second thing fell 463 -> 239.
    //
    // The discriminator is already on the record: a nested haploid site names the strand its allele
    // sits on, so the *other* strand is the empty one. It has to reach segmentation and not only the
    // writer, because a run is cut where the haplotype changes -- and without this a run could mix
    // both kinds and then have no single character to print.
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
    size_t empty_segments = 0, unexplained_segments = 0;

    std::function<void(size_t, size_t, int, size_t, StrandKind)> emit_span =
        [&](size_t from, size_t to, int strand, size_t hap, StrandKind kind) {
        gbwt::edge_type pos = (hap == LinkageModel::WILDCARD)
                                  ? gbwt::invalid_edge()
                                  : mosaic_gbwt_position(phasing[from].start_node, hap);

        auto emit_row = [&](size_t a_idx, size_t b_idx, gbwt::edge_type p) {
            const LinkageCollector::PhaseCall& a = phasing[a_idx];
            const LinkageCollector::PhaseCall& b = phasing[b_idx];
            out << "H\t" << a.contig << "\t" << strand << "\t"
                << a.position << "\t" << b.position << "\t"
                << a.start_node << "\t" << b.end_node << "\t";
            if (hap == LinkageModel::WILDCARD) {
                // '.' for a strand with no sequence here, matching what the VCF already means by a
                // '.' in a phased GT; '*' stays "the panel cannot name a haplotype".
                if (kind == StrandKind::Empty) {
                    out << ".\t.";
                    ++empty_segments;
                } else {
                    out << "*\t*";
                    ++unexplained_segments;
                }
            } else {
                out << hap << "\t"
                    << (hap < mosaic_haplotype_names.size()
                            ? mosaic_haplotype_names[hap] : string("?"));
            }
            out << "\t" << (b_idx - a_idx + 1) << "\t";
            if (p == gbwt::invalid_edge()) {
                // No position to give: either the strand is the wildcard, or the panel does not
                // place this haplotype on this node. '.' rather than a zero, which would read as a
                // valid offset. On chr20 ten segments are in the latter case, all of them in the
                // centromere where node IDs stop following reference order.
                out << ".\t.";
            } else {
                out << p.first << "\t" << p.second;
            }
            out << "\n";
            ++total_segments;
        };

        if (pos == gbwt::invalid_edge() || from == to) {
            emit_row(from, to, pos);
            return;
        }
        gbwt::edge_type end_pos = mosaic_gbwt_position(phasing[to].start_node, hap);
        if (end_pos == gbwt::invalid_edge()
            || linkage_gbwt->locate(pos) == linkage_gbwt->locate(end_pos)) {
            // Same fragment at both ends, or no way to tell. One row.
            emit_row(from, to, pos);
            return;
        }
        // The fragment changes somewhere in (from, to]. Binary search for the last site still on
        // the starting fragment; a site the haplotype does not reach is treated as past the
        // boundary, which keeps the search monotone.
        gbwt::size_type seq = linkage_gbwt->locate(pos);
        size_t lo = from, hi = to;
        while (hi - lo > 1) {
            size_t mid = lo + (hi - lo) / 2;
            gbwt::edge_type p = mosaic_gbwt_position(phasing[mid].start_node, hap);
            if (p != gbwt::invalid_edge() && linkage_gbwt->locate(p) == seq) {
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
        // One strand on a haploid chain. Emitting a second would assert a homozygote where the
        // sample has one copy, which is the opposite of what the file is for.
        //
        // Over the whole run, not from its first site. A diploid contig can open with a nested
        // haploid site now that the phasing is sorted into reference order, and reading the ploidy
        // off phasing[i] would then drop the contig's entire second strand -- silently, since a
        // one-strand mosaic is exactly what a haploid contig is meant to look like. chr20 happens
        // to open at position 24 with a diploid record, which is the only reason this was not
        // already visible.
        int strands = 1;
        for (size_t t = i; t < j && strands == 1; ++t) {
            if (phasing[t].ploidy != 1) {
                strands = 2;
            }
        }
        for (int strand = 0; strand < strands; ++strand) {
            size_t seg_start = i;
            for (size_t t = i; t < j; ++t) {
                size_t hap = strand == 0 ? phasing[t].hap_first : phasing[t].hap_second;
                StrandKind kind = strand_kind(t, strand);
                bool last = (t + 1 == j);
                // Cut on the kind too: two adjacent wildcard sites can mean different things, and a
                // run holding both has no single character to print for its haplotype column.
                bool changes = !last
                               && ((strand == 0 ? phasing[t + 1].hap_first
                                                : phasing[t + 1].hap_second) != hap
                                   || strand_kind(t + 1, strand) != kind);
                if (last || changes) {
                    // One run of a single haplotype, but possibly several GBWT fragments of it.
                    // A segment must be walkable from one position, so the run is cut wherever the
                    // fragment under it changes: emit_span walks the sites, re-resolving whenever
                    // the position stops belonging to the same fragment, and emits one row per
                    // fragment. Without this a consumer following LF() from the segment's position
                    // would hit the endmarker partway and have no way to pick up the rest.
                    emit_span(seg_start, t, strand, hap, kind);
                    seg_start = t + 1;
                }
            }
        }
        i = j;
    }
    cerr << "[vg call] mosaic: " << total_segments << " segments over " << phasing.size()
         << " sites, written to " << mosaic_path << endl;
    // Reported apart because they are different facts and the sum is the figure the old file
    // reported as one. The unexplained half is the one comparable to the phasing report's count.
    cerr << "[vg call] mosaic: " << empty_segments << " segments on a strand with no sequence there, "
         << unexplained_segments << " the panel cannot name a haplotype for ("
         << (empty_segments + unexplained_segments) << " were one figure before)" << endl;
}

bool VCFOutputCaller::apply_linkage_quality(string& line, double posterior,
                                            double explained_share) const {
    // The quality half of what `apply_linkage_change` used to do, kept when the genotype half
    // became unnecessary. It is not decoration, and losing it is invisible to F1: GQ is not a
    // filter here, so a run whose quality silently reverted to the per-site value scores
    // identically and is differently calibrated. That is exactly what happened when the genotype
    // patch stopped being produced -- 9,980 records carried posterior-derived quality at stage 9,
    // 3,524 at stage 10, and none at all once the record was built from the settled genotype.
    vector<string> fields;
    {
        size_t start = 0;
        while (true) {
            size_t tab = line.find('\t', start);
            fields.push_back(line.substr(start, tab == string::npos ? string::npos : tab - start));
            if (tab == string::npos) {
                break;
            }
            start = tab + 1;
        }
    }
    if (fields.size() < 10) {
        return false;
    }
    vector<string> keys, values;
    {
        size_t start = 0;
        while (true) {
            size_t colon = fields[8].find(':', start);
            keys.push_back(fields[8].substr(start,
                                            colon == string::npos ? string::npos : colon - start));
            if (colon == string::npos) {
                break;
            }
            start = colon + 1;
        }
        start = 0;
        while (true) {
            size_t colon = fields[9].find(':', start);
            values.push_back(fields[9].substr(start,
                                              colon == string::npos ? string::npos : colon - start));
            if (colon == string::npos) {
                break;
            }
            start = colon + 1;
        }
    }
    if (keys.size() != values.size()) {
        return false;
    }
    size_t gq_field = keys.size(), gqi_field = keys.size(), gqn_field = keys.size();
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "GQ") {
            gq_field = i;
        } else if (keys[i] == "GQI") {
            gqi_field = i;
        } else if (keys[i] == "GQN") {
            gqn_field = i;
        }
    }
    if (gq_field != keys.size()) {
        // GQ becomes the phred-scaled complement of the posterior, then discounted by the
        // explained-read share exactly as the per-site GQ was: reads whose best allele lies outside
        // the call enter every genotype's likelihood and cancel, so neither the likelihood ratio nor
        // the posterior built from that emission can see them.
        //
        // The cap at GQI is the larger effect and is not a consistency nicety. The posterior is
        // computed under a strong prior (the panel frequency exponent defaults to 5), so
        // `1 - posterior` understates uncertainty, and it does so exactly where the per-site
        // evidence was weakest -- which is where the linkage layer acts at all. Capping re-anchors
        // the reported quality to the read evidence: linkage may lower confidence and may not raise
        // it above what the reads alone supported. Worth about +0.003 AUC and 1-2% fewer surviving
        // false calls, against +0.0001 to +0.0009 for the share discount alone. It is also what
        // makes `GQ <= GQI` hold on every record (test/t/18_vg_call.t); the discount alone does not,
        // because the posterior quality has no arithmetic relation to GQI to be bounded by.
        double q = posterior >= 1.0 ? 256.0 : -10.0 * log10(max(1.0 - posterior, 1e-26));
        q *= explained_share;
        if (gqi_field != keys.size()) {
            try {
                q = min(q, stod(values[gqi_field]));
            } catch (const std::exception&) {
                // GQI absent or unparsable: leave the discounted posterior uncapped rather than
                // dropping the quality entirely.
            }
        }
        // Truncated, not rounded, because that is what the per-site emission does -- two records
        // with the same underlying quality must not print different integers depending on whether
        // linkage touched one of them.
        values[gq_field] = std::to_string((int)min(256.0, max(0.0, q)));
    }
    if (gqn_field != keys.size()) {
        // GQN is the per-site likelihood-ratio gap as a fraction of the site's maximum, and it
        // described the genotype linkage just moved away from. It is not derivable from the
        // posterior, so it is blanked -- "." means no measurement, which is now the truth -- rather
        // than left describing a call the record no longer makes.
        values[gqn_field] = ".";
    }
    // `lowconf` was set from the pre-linkage GQN against --min-confidence, so it labels the
    // abandoned call. The record's confidence is the GQ written above; the stale flag is cleared
    // rather than left to gate downstream filtering on the wrong genotype's confidence.
    if (fields[6] == "lowconf") {
        fields[6] = "PASS";
    }

    string format;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) {
            format += ":";
        }
        format += values[i];
    }
    fields[9] = format;
    line.clear();
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) {
            line += "\t";
        }
        line += fields[i];
    }
    return true;
}


static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                if (g->first == alleleIndex) {
                    alts += g->second;
                }
            }
        }
    }
    return alts;
}

static int countAlleles(vcflib::Variant& var) {
    int alleles = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
		    alleles += g->second;
		}
            }
        }
    }
    return alleles;
}

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in the file again and again
static int countSamplesWithData(vcflib::Variant& var) {
    int samples_with_data = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        bool has_data = false;
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
                    has_data = true;
                    break;
		}
            }
        }
        if (has_data) {
            ++samples_with_data;
        }
    }
    return samples_with_data;
}

void VCFOutputCaller::vcf_fixup(vcflib::Variant& var) const {
    // copied from https://github.com/vgteam/vcflib/blob/master/src/vcffixup.cpp
    
    stringstream ns;
    ns << countSamplesWithData(var);
    var.info["NS"].clear();
    var.info["NS"].push_back(ns.str());

    var.info["AC"].clear();
    var.info["AF"].clear();
    var.info["AN"].clear();

    int allelecount = countAlleles(var);
    stringstream an;
    an << allelecount;
    var.info["AN"].push_back(an.str());

    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& allele = *a;
        int altcount = countAlts(var, var.getAltAlleleIndex(allele) + 1);
        stringstream ac;
        ac << altcount;
        var.info["AC"].push_back(ac.str());
        stringstream af;
        double faf = (double) altcount / (double) allelecount;
        if(faf != faf) faf = 0;
        af << faf;
        var.info["AF"].push_back(af.str());
    }
}

void VCFOutputCaller::set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation) {
    this->translation = translation;
}

void VCFOutputCaller::set_nested(bool nested) {
    include_nested = nested;
}

void VCFOutputCaller::set_allele_merge(double threshold, int64_t min_len) {
    allele_merge_threshold = threshold;
    allele_merge_min_len = min_len;
}

bool VCFOutputCaller::snarl_traversal_to_handles(const HandleGraph& graph, const SnarlTraversal& trav,
                                                 Traversal& out_trav) {
    // cluster_traversals asserts size() >= 2, and a Visit carrying a child Snarl has no single
    // handle.  Both are real inputs here (the "*" placeholder, and NestedFlowCaller traversals via
    // SnarlGraph::embed_snarl), so refuse rather than fabricate something.
    if (trav.visit_size() < 2) {
        return false;
    }
    out_trav.clear();
    out_trav.reserve(trav.visit_size());
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() <= 0) {
            return false;
        }
        out_trav.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return true;
}

namespace {
/// Parse a VCF FORMAT value without throwing and without exiting.  vg::parse<double> exits on
/// failure and the 2-argument vg::parse can throw; merge_similar_alleles runs inside an OpenMP
/// region, where an escaping exception is std::terminate rather than something a caller can handle,
/// and exit() would abandon whatever the other threads had already buffered.  Missing values
/// ("." and "") are ordinary input here.
bool parse_vcf_double(const string& field, double& value) {
    try {
        size_t after;
        value = std::stod(field, &after);
        return after == field.size();
    } catch (const std::exception&) {
        return false;
    }
}
}

int64_t VCFOutputCaller::allele_core_length(const vector<string>& alleles) {
    vector<const string*> seqs;
    for (const string& a : alleles) {
        if (a != "*") {
            seqs.push_back(&a);
        }
    }
    if (seqs.empty()) {
        return 0;
    }
    size_t min_len = seqs[0]->length();
    size_t max_len = 0;
    for (const string* s : seqs) {
        min_len = std::min(min_len, s->length());
        max_len = std::max(max_len, s->length());
    }
    // The prefix and the suffix may not overlap, exactly as in flatten_common_allele_ends: a shared
    // region can only be counted once, or {"AAAA","AAAAA"} would come out at -3 instead of 1.
    // Case-insensitive to match flatten's own toupper and deconstruct's toUppercase.
    auto shared = [&](size_t skip, bool from_back) {
        auto at = [&](const string* s, size_t i) {
            return std::toupper((*s)[from_back ? s->length() - 1 - i : i]);
        };
        size_t n = 0;
        while (skip + n < min_len) {
            int ch = at(seqs[0], n);
            bool match = true;
            for (size_t j = 1; j < seqs.size() && match; ++j) {
                match = at(seqs[j], n) == ch;
            }
            if (!match) {
                break;
            }
            ++n;
        }
        return n;
    };
    size_t prefix = shared(0, false);
    size_t suffix = shared(prefix, true);
    // non-negative structurally, not by clamping: the loop caps give prefix + suffix <= min_len <= max_len
    return (int64_t)(max_len - prefix - suffix);
}

bool VCFOutputCaller::merge_similar_alleles(const PathPositionHandleGraph& graph,
                                            const vector<SnarlTraversal>& site_traversals,
                                            vector<int>& site_genotype,
                                            const string& sample_name,
                                            vcflib::Variant& out_variant,
                                            GLLayout gl_layout) const {
    if (!(allele_merge_threshold < 1.0)) {
        return false;
    }
    // we only collapse a genotype that actually calls two distinct ALTs.  This also keeps the -a
    // padding block above (which adds uncalled alleles) out of scope: those alleles are advertised,
    // not called, and rewriting them without a genotype change would be a silent surprise.
    set<int> called_alts;
    for (int g : site_genotype) {
        if (g > 0) {
            called_alts.insert(g);
        }
    }
    if (called_alts.size() < 2) {
        return false;
    }
    // Per-site gate, decided over the alleles this record actually emits.  NOT over the traversal
    // finder's candidate list: that is up to max_yens_traversals (50) speculative paths, most of
    // which never become an allele, so gating on them lets an invisible branch with no reads and no
    // AT entry decide whether merging happens.  deconstruct's equivalent gate is decided over the
    // set that becomes ITS alleles -- the reference plus everything get_traversal_order kept.  Both
    // tools gate on what they emit, but those sets differ (we see only the called genotype's
    // alleles, deconstruct sees every haplotype), so the two can disagree at a site whose uncalled
    // haplotypes are much larger than its called ones.
    // The quantity is CORE LENGTH (see allele_core_length): the longest allele once the prefix and
    // suffix shared by every allele are stripped.  Raw string length would answer differently from
    // vg deconstruct on the same variant, because this record has been flattened down to an anchor
    // base and deconstruct's has not.
    if (allele_merge_min_len > 0 &&
        allele_core_length(out_variant.alleles) < allele_merge_min_len) {
        return false;
    }

    // ALT-vs-ALT only.  Absorbing an ALT into allele 0 would empty out_variant.alt and the record
    // would then be dropped entirely by the caller, turning a het call into no call at all.
    // (vg deconstruct does fold near-reference alleles into the reference cluster and drop the
    //  record; that is deliberate there and deliberately not copied here.)
    vector<Traversal> alt_travs;
    vector<int> alt_to_allele;
    for (size_t i = 1; i < site_traversals.size(); ++i) {
        if (!called_alts.count((int)i)) {
            continue;
        }
        Traversal trav;
        if (!snarl_traversal_to_handles(graph, site_traversals[i], trav)) {
            // star placeholder or a child-snarl visit: leave this allele alone
            continue;
        }
        alt_travs.push_back(std::move(trav));
        alt_to_allele.push_back((int)i);
    }
    if (alt_travs.size() < 2) {
        return false;
    }

    // same clustering call deconstruct makes, so the metric, the endpoint pruning and the
    // >= comparison are inherited rather than reimplemented.
    //
    // Cluster in descending allele-depth order, so each cluster's head -- the allele that survives,
    // and the one MAT's similarity is measured against -- is its best-supported member.  Identity order
    // would instead inherit the traversal finder's ranking, which FlowCaller::call_snarl_internal (and NestedFlowCaller's copy of it) switches to
    // length-weighted average flow once a snarl's interior passes the average-support threshold.
    // That ranking can put a short, lightly-supported allele ahead of a long, heavily-supported one,
    // and merging into it emits the minority sequence as a homozygous call carrying the pooled depth.
    vector<int> order(alt_travs.size());
    std::iota(order.begin(), order.end(), 0);
    {
        auto& sample_fields = out_variant.samples[sample_name];
        auto ad_it = sample_fields.find("AD");
        if (ad_it != sample_fields.end() && ad_it->second.size() == out_variant.alleles.size()) {
            vector<double> ad(alt_travs.size(), 0);
            bool usable = true;
            for (size_t k = 0; k < alt_travs.size() && usable; ++k) {
                usable = parse_vcf_double(ad_it->second.at(alt_to_allele[k]), ad[k]);
            }
            if (usable) {
                // stable, so equal depths keep the finder's own ranking
                std::stable_sort(order.begin(), order.end(),
                                 [&](int a, int b) { return ad[a] > ad[b]; });
            }
        }
    }
    // The VCF reference allele sets the scale a pure deletion is measured against.  It is
    // site_traversals[0] and is deliberately not clustered (the loop above starts at 1).
    // The nullptr fallback is currently unreachable: the only producer of a Visit without a node is
    // NestedFlowCaller, and --bottom-up is rejected with -L.  It is kept because the consequence of
    // being wrong about that is a crash, and because falling back to pairwise scoring can only
    // merge less, never more.
    Traversal ref_trav;
    const Traversal* site_ref_trav = nullptr;
    if (!site_traversals.empty() && snarl_traversal_to_handles(graph, site_traversals[0], ref_trav)) {
        site_ref_trav = &ref_trav;
    }
    vector<pair<double, int64_t>> cluster_info;
    vector<int> unused_child_snarl_mapping;
    vector<vector<int>> clusters = cluster_traversals(&graph, alt_travs, order,
                                                      vector<pair<handle_t, handle_t>>(),
                                                      allele_merge_threshold,
                                                      cluster_info, unused_child_snarl_mapping,
                                                      site_ref_trav);

    // merge_to[a] == a for a surviving allele, else the allele it collapses into
    vector<int> merge_to(out_variant.alleles.size());
    std::iota(merge_to.begin(), merge_to.end(), 0);
    vector<string> mat_entries;
    bool merged_any = false;
    for (const vector<int>& cluster : clusters) {
        if (cluster.size() < 2) {
            continue;
        }
        int survivor = alt_to_allele[cluster.front()];
        for (size_t j = 1; j < cluster.size(); ++j) {
            int absorbed = alt_to_allele[cluster[j]];
            merge_to[absorbed] = survivor;
            merged_any = true;
            stringstream ss;
            ss.precision(3);
            ss << absorbed << ">" << survivor << ":" << cluster_info[cluster[j]].first;
            mat_entries.push_back(ss.str());
        }
    }
    if (!merged_any) {
        return false;
    }

    // dense renumbering of the survivors, preserving order
    vector<int> new_index(merge_to.size(), -1);
    int next = 0;
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_index[a] = next++;
        }
    }
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] != (int)a) {
            new_index[a] = new_index[merge_to[a]];
        }
    }
    int n_new = next;

    // alleles / alt
    vector<string> new_alleles(n_new);
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_alleles[new_index[a]] = out_variant.alleles[a];
        }
    }
    out_variant.alleles = new_alleles;
    out_variant.alt.assign(new_alleles.begin() + 1, new_alleles.end());

    // AT is Number=R, so it is indexed by allele just like alleles
    auto at_it = out_variant.info.find("AT");
    if (at_it != out_variant.info.end() && at_it->second.size() == merge_to.size()) {
        vector<string> new_at(n_new);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            if (merge_to[a] == (int)a) {
                new_at[new_index[a]] = at_it->second[a];
            }
        }
        at_it->second = new_at;
    }

    auto& sample = out_variant.samples[sample_name];

    // AD is a per-allele count, so the absorbed allele's reads move onto the survivor.  sum(AD) is
    // therefore unchanged by the merge; it can slightly over-count when the merged alleles share
    // interior nodes, whose depth was proportionally split between them.
    auto ad_it = sample.find("AD");
    if (ad_it != sample.end() && ad_it->second.size() == merge_to.size()) {
        vector<double> summed(n_new, 0);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            double v = 0;
            // treat an unparseable entry as 0 rather than bailing: the merge is already committed,
            // and dropping AD would leave the record with no per-allele depth at all
            parse_vcf_double(ad_it->second[a], v);
            summed[new_index[a]] += v;
        }
        vector<string> new_ad(n_new);
        for (int a = 0; a < n_new; ++a) {
            new_ad[a] = std::to_string((int64_t)std::llround(summed[a]));
        }
        ad_it->second = new_ad;
        // MAD is the min allele depth over the called alleles; recompute so it agrees with the AD
        // and GT printed beside it.  FILTER is left as computed pre-merge, and since the new MAD is
        // >= the old one, a depth filter can only over-filter, never under-filter.
        auto mad_it = sample.find("MAD");
        if (mad_it != sample.end() && mad_it->second.size() == 1) {
            double min_ad = -1;
            for (int g : site_genotype) {
                if (g >= 0 && g < (int)merge_to.size()) {
                    double v = summed[new_index[g]];
                    if (min_ad < 0 || v < min_ad) {
                        min_ad = v;
                    }
                }
            }
            if (min_ad >= 0) {
                mad_it->second[0] = std::to_string((int64_t)std::llround(min_ad));
            }
        }
    }

    // GL is Number=G, and which order it is in depends on which caller wrote it. The comment that
    // stood here said vg emits i-major and named PoissonSupportSnarlCaller as the only writer of GL.
    // That is not true: ReadLikelihoodSnarlCaller::update_vcf_info writes it in
    // AlleleReadLikelihoods::enumerate_genotypes order, which is the VCF spec's colexicographic one,
    // and at three alleles the two disagree at indices 2 and 3 -- (1,1) against (0,2). So a merged
    // record under --read-likelihood had two of its six likelihoods transposed, which silently
    // max-marginalises the het-with-allele-2 class into the hom-allele-1 class and back.
    //
    // The layout comes from the caller rather than being assumed, because the Poisson path really is
    // i-major and still live: indexing everything through enumerate_genotypes would fix one caller by
    // corrupting the other.
    auto gl_it = sample.find("GL");
    if (gl_it != sample.end()) {
        size_t n_old = merge_to.size();
        // The diploid layout is the only one that can occur: merging needs at least two distinct
        // called ALTs, site_genotype has one entry per ploidy, and both GL writers assert ploidy is
        // 1 or 2. So n_old is always 3 here.
        assert(gl_it->second.size() == n_old * (n_old + 1) / 2);
        bool gl_usable = true;
        vector<double> old_gl(gl_it->second.size(), 0.0);
        for (size_t g = 0; g < gl_it->second.size() && gl_usable; ++g) {
            gl_usable = parse_vcf_double(gl_it->second[g], old_gl[g]);
        }
        vector<double> folded;
        if (gl_usable) {
            folded = fold_genotype_likelihoods(old_gl, new_index, (size_t)n_new, gl_layout);
        }
        if (gl_usable) {
            vector<string> new_gl(folded.size());
            for (size_t i = 0; i < folded.size(); ++i) {
                new_gl[i] = std::to_string(folded[i]);
            }
            gl_it->second = new_gl;
        } else {
            // A value we cannot parse.  Leaving GL alone would emit a Number=G field whose length
            // disagrees with the new allele count, so drop it rather than lie.
            sample.erase(gl_it);
            auto& fmt = out_variant.format;
            fmt.erase(std::remove(fmt.begin(), fmt.end(), string("GL")), fmt.end());
        }
    }
    // GQ and GP are deliberately untouched: they come from the caller's own CallInfo, computed over
    // its candidate set rather than from the emitted GL, so recomputing them here would silently
    // swap one statistic for another.

    // GT, from the renumbered genotype
    for (int& g : site_genotype) {
        if (g >= 0 && g < (int)merge_to.size()) {
            g = new_index[g];
        }
    }
    stringstream vcf_gt;
    for (size_t i = 0; i < site_genotype.size(); ++i) {
        if (site_genotype[i] == MISSING_ALLELE_MARKER) {
            vcf_gt << ".";
        } else {
            vcf_gt << site_genotype[i];
        }
        if (i != site_genotype.size() - 1) {
            vcf_gt << "/";
        }
    }
    sample["GT"] = {vcf_gt.str()};

    // record what happened: without this a merged 1/1 is indistinguishable from a real hom-alt
    out_variant.info["MAT"] = mat_entries;

    out_variant.updateAlleleIndexes();
    return true;
}

unordered_set<string> VCFOutputCaller::get_output_contigs() const {
    unordered_set<string> contigs;
    // The sort key is (sequenceName, position) (see add_variant), so the contig is right
    // there and nothing has to be decompressed.
    for (const auto& thread_buf : output_variants) {
        for (const auto& output_variant_record : thread_buf) {
            contigs.insert(output_variant_record.first.contig);
        }
    }
    return contigs;
}

string VCFOutputCaller::prune_header_contigs(const string& header,
                                             const unordered_set<string>& keep) const {
    static const string contig_prefix = "##contig=<ID=";
    stringstream pruned;
    vector<string> lines = split_delims(header, "\n");
    for (const string& line : lines) {
        if (line.compare(0, contig_prefix.size(), contig_prefix) == 0) {
            // Parse the ID back out the same way it was written, rather than scanning for a
            // delimiter: contig names are path names and nothing stops one containing ',' or
            // '>'.  Both producers emit exactly ##contig=<ID=NAME,length=N> -- vcf_header()
            // above and Deconstructor::add_contigs_to_vcf_header().
            static const string contig_suffix = ",length=";
            size_t id_start = contig_prefix.size();
            size_t id_end = line.rfind(contig_suffix);
            if (id_end == string::npos || id_end < id_start) {
                // not a shape we wrote; leave it alone rather than guess
                pruned << line << "\n";
                continue;
            }
            string id = line.substr(id_start, id_end - id_start);
            if (!keep.count(id)) {
                continue;
            }
        }
        pruned << line << "\n";
    }
    string result = pruned.str();
    if (!header.empty() && header.back() != '\n' && !result.empty()) {
        // input had no trailing newline, so don't invent one
        result.pop_back();
    }
    return result;
}

void VCFOutputCaller::add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele, const Traversal& trav,
                                              bool reversed, bool one_based) const {
    SnarlTraversal proto_trav;
    for (const handle_t& handle : trav) {
        Visit* visit = proto_trav.add_visit();
        visit->set_node_id(graph->get_id(handle));
        visit->set_backward(graph->get_is_reverse(handle));
    }
    this->add_allele_path_to_info(v, allele, proto_trav, reversed, one_based);
}

void VCFOutputCaller::add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav,
                                              bool reversed, bool one_based) const {
    auto& trav_info = v.info["AT"];
    assert(allele < trav_info.size());

    vector<int> nodes;
    nodes.reserve(trav.visit_size());
    const Visit* prev_visit = nullptr;
    unordered_map<nid_t, pair<string, size_t>>::const_iterator prev_trans;
    
    for (size_t i = 0; i < trav.visit_size(); ++i) {
        size_t j = !reversed ? i : trav.visit_size() - 1 - i;
        const Visit& visit = trav.visit(j);
        nid_t node_id = visit.node_id();
        string node_name = std::to_string(node_id);
        bool skip = false;
        // todo: check one_based? (we kind of ignore that when writing the snarl name, so maybe not pertienent)
        if (translation) {
            auto i = translation->find(node_id);
            if (i == translation->end()) {
                throw runtime_error("Error [vg deconstruct]: Unable to find node " + node_name + " in translation file");
            }
            if (prev_visit) {
                nid_t prev_node_id = prev_visit->node_id();
                if (prev_trans->second.first == i->second.first && node_id != prev_node_id) {
                    // here is a case where we have two consecutive nodes that map back to
                    // the same source node.
                    // todo: check if translation node properly covered
                    skip = true;
                }
            }
            node_name = i->second.first;
            prev_trans = i;
        }

        if (!skip) {
            bool vrev = visit.backward() != reversed;
            trav_info[allele] += (vrev ? "<" : ">");
            trav_info[allele] += node_name;
        }
        prev_visit = &visit;
    }
    if (trav_info[allele].empty()) {
        // note: * alleles get empty traversals
        trav_info[allele] = ".";
    }
}

string VCFOutputCaller::trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const {
    string seq;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() > 0) {
            seq += graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        } else {
            seq += print_snarl(visit.snarl(), true);
        }
    }
    return seq;    
}

thread_local VCFOutputCaller::NestedContext VCFOutputCaller::nested_context;
thread_local size_t VCFOutputCaller::current_generation = 0;
thread_local VCFOutputCaller::EmittedAlleles VCFOutputCaller::last_emitted;

bool VCFOutputCaller::chain_reported_inline(const Snarl& snarl,
                                            const vector<SnarlTraversal>& travs,
                                            const vector<int>& genotype, int ref_trav_idx,
                                            const Snarl& child) const {
    // Inert unless block emission is on. With one record per snarl the parent's ALT spans the whole
    // snarl anyway, so there is no sense in which a chain is "inside a block".
    if (!atomize_blocks || symbolic_manager == nullptr) {
        return false;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size() || genotype.empty()) {
        return false;
    }
    // A snarl whose projection has no symbols at all cannot answer this question: every child would
    // read as "not matched" and the whole subtree would be dropped rather than delegated. That is
    // the flipped-snarl case, and getting it wrong turns double reporting into non-reporting.
    if (!symbolic_site_resolvable(snarl, *symbolic_manager)) {
        return false;
    }

    const Snarl* managed_child = symbolic_manager->into_which_snarl(child.start().node_id(),
                                                                   child.start().backward());
    if (managed_child == nullptr) {
        return false;
    }
    pair<nid_t, nid_t> bounds = chain_bounds_of(managed_child, *symbolic_manager);

    SymbolicAllele sref = symbolic_allele(travs[ref_trav_idx], snarl, *symbolic_manager);
    // Where the chain sits in the reference projection. If it is nowhere, the reference does not
    // cross it and the caller's own reference gate has already dealt with that.
    vector<size_t> chain_steps;
    for (size_t i = 0; i < sref.size(); ++i) {
        if (sref[i].is_chain() && sref[i].id == bounds.first && sref[i].end_id == bounds.second) {
            chain_steps.push_back(i);
        }
    }
    if (chain_steps.empty()) {
        return false;
    }

    bool any_crossing = false;
    for (int allele : genotype) {
        if (allele < 0 || (size_t)allele >= travs.size()) {
            continue;
        }
        if (allele == ref_trav_idx) {
            // This haplotype IS the reference here, so every reference step matches, the chain
            // among them. Delegated, and nothing further to check.
            return false;
        }
        SymbolicAllele salt = symbolic_allele(travs[allele], snarl, *symbolic_manager);
        bool degraded = false;
        vector<DiffBlock> blocks = symbolic_diff(sref, salt, &degraded);
        if (degraded) {
            return false;   // no reliable block structure, so no reliable conclusion
        }
        // This haplotype's OWN crossings, found in its own projection.
        //
        // Testing the reference's chain step instead is wrong, and wrong in a way that fires 60x
        // too often: a haplotype that DELETES the chain puts the reference's step inside a block
        // while crossing the chain zero times, so it would read as "already reported" when it
        // contributes no copy at all. Those chains are exactly the ones the caller retains for
        // possible revision when no called allele reaches them yet, and suppressing them here
        // deleted 399 records that linkage was still entitled to move.
        for (size_t j = 0; j < salt.size(); ++j) {
            if (!salt[j].is_chain() || salt[j].id != bounds.first ||
                salt[j].end_id != bounds.second) {
                continue;
            }
            any_crossing = true;
            bool inside = false;
            for (const DiffBlock& b : blocks) {
                if ((size_t)b.alt_begin <= j && j < (size_t)b.alt_end) {
                    inside = true;
                    break;
                }
            }
            if (!inside) {
                return false;   // matched on this haplotype: its own record reports it
            }
        }
    }

    if (!any_crossing) {
        // No called haplotype crosses this chain, so there is no copy for a block ALT to have
        // spelled and nothing here to suppress. Whether it is visited anyway is the caller's
        // decision, not this predicate's.
        return false;
    }

    // Every crossing by every called haplotype falls inside a difference block, so the block's ALT
    // already spells the route through it. Reporting it again would be the same variation twice.
    ++g_atomize_child_inlined;
    return true;
}

bool VCFOutputCaller::is_symbolically_reference(const vector<SnarlTraversal>& called_traversals,
                                                int trav_idx, int ref_trav_idx,
                                                const Snarl& snarl) const {
    // Off unless a snarl hierarchy was supplied, so the default path compares alleles by sequence
    // exactly as it did before.
    if (symbolic_manager == nullptr || ref_trav_idx < 0 || trav_idx < 0 ||
        ref_trav_idx >= (int)called_traversals.size() ||
        trav_idx >= (int)called_traversals.size()) {
        return false;
    }
    return symbolically_equal(called_traversals[trav_idx], called_traversals[ref_trav_idx],
                              snarl, *symbolic_manager);
}


/// Stage-2 instrumentation for planning/symbolic-diff-decomposition.md. Projects the reference and
/// each distinct called ALT traversal, aligns them, and tallies. Writes only atomics and reads
/// nothing it can alter, so a run with this compiled in is byte-identical to one without.
static void tally_atomize(const PathPositionHandleGraph& graph, const SnarlManager* mgr,
                          const Snarl& snarl, const vector<SnarlTraversal>& travs,
                          const vector<int>& genotype, int ref_trav_idx) {
    if (mgr == nullptr || ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size()) {
        return;
    }
    bool site_reversed = false;
    if (!symbolic_site_resolvable(snarl, *mgr, &site_reversed)) {
        // Projection would report a bare node list here, so a block count from it would measure
        // node-level shredding rather than chain structure. Counted and skipped, not folded in.
        ++g_atomize_site_unresolvable;
        return;
    }
    if (site_reversed) {
        // Counted HERE, once per record, so it is commensurable with the counter above -- the one
        // that read 9,279 before the reversed pairing was accepted. Counting inside the resolver
        // instead would count calls: projection runs per traversal, so a single site would bump it
        // once per allele per haplotype and the number would look like an over-fire.
        ++g_atomize_site_reversed;
    }

    vector<pair<int, int>> ref_ranges;
    SymbolicAllele sref = symbolic_allele(travs[ref_trav_idx], snarl, *mgr, &ref_ranges);

    // Which chains the reference itself crosses, so an unmatched chain the reference also visits
    // (double-reportable) can be told from one it never visits (reportable only inside this allele).
    set<pair<nid_t, nid_t>> ref_chains;
    for (const SymbolicStep& s : sref) {
        if (s.is_chain()) {
            ref_chains.emplace(s.id, s.end_id);
        }
    }

    auto span_bases = [&](const SnarlTraversal& t, const vector<pair<int, int>>& ranges,
                          int step_begin, int step_end) -> size_t {
        size_t n = 0;
        for (int k = step_begin; k < step_end && k < (int)ranges.size(); ++k) {
            for (int v = ranges[k].first; v < ranges[k].second && v < t.visit_size(); ++v) {
                const Visit& vis = t.visit(v);
                if (vis.node_id() > 0) {
                    n += graph.get_length(graph.get_handle(vis.node_id()));
                }
            }
        }
        return n;
    };

    set<int> seen;
    for (int gt : genotype) {
        if (gt < 0 || gt == ref_trav_idx || (size_t)gt >= travs.size()) {
            continue;
        }
        if (!seen.insert(gt).second) {
            continue;
        }

        vector<pair<int, int>> alt_ranges;
        SymbolicAllele salt = symbolic_allele(travs[gt], snarl, *mgr, &alt_ranges);

        bool degraded = false;
        vector<DiffBlock> blocks = symbolic_diff(sref, salt, &degraded);
        if (degraded) {
            ++g_atomize_degraded;
        }

        ++g_atomize_records;
        g_atomize_blocks_total += blocks.size();
        ++g_atomize_blocks_hist[std::min<size_t>(blocks.size(), 15)];
        if (blocks.size() >= 2) {
            ++g_atomize_multi_block;
        }

        size_t alt_total = span_bases(travs[gt], alt_ranges, 0, (int)salt.size());
        bool any_offref = false;
        bool any_case_c = false;
        size_t offref_bases = 0;
        size_t smallest_block = std::numeric_limits<size_t>::max();

        for (const DiffBlock& b : blocks) {
            size_t rb = span_bases(travs[ref_trav_idx], ref_ranges, b.ref_begin, b.ref_end);
            size_t ab = span_bases(travs[gt], alt_ranges, b.alt_begin, b.alt_end);
            smallest_block = std::min(smallest_block, std::max(rb, ab));
            for (int k = b.alt_begin; k < b.alt_end && k < (int)salt.size(); ++k) {
                if (!salt[k].is_chain()) {
                    continue;
                }
                if (ref_chains.count(make_pair(salt[k].id, salt[k].end_id))) {
                    any_case_c = true;
                } else {
                    any_offref = true;
                    offref_bases += span_bases(travs[gt], alt_ranges, k, k + 1);
                }
            }
        }
        if (any_case_c) {
            ++g_atomize_case_c;
        }
        if (any_offref) {
            ++g_atomize_offref_alts;
            g_atomize_offref_bases += offref_bases;
            g_atomize_alt_bases += alt_total;
        }
        // Only a genuine split can move a variant across truvari's 50 bp filter. A single block
        // spans everything but the two boundary steps, so comparing it against the whole
        // traversal's length measures how long the boundary nodes are and nothing else -- which is
        // what the first version of this counter did, reporting 74,885 of 115,996 ALTs as
        // "exposed" when the real exposure cannot exceed the number of ALTs that split at all.
        if (blocks.size() >= 2 && alt_total >= 50 &&
            smallest_block != std::numeric_limits<size_t>::max() && smallest_block < 50) {
            ++g_atomize_small_block_on_sv;
        }
    }
}

/// A nested haploid genotype: one allele on a named strand, with "." on the other.
///
/// "." rather than "0" because the other haplotype does not carry the reference sequence here -- it
/// carries nothing, the parent's other allele having deleted the chain. One place rather than two,
/// so a site record and the difference blocks it decomposes into cannot drift on how one call is
/// spelled.
static string nested_strand_genotype(int allele, int strand) {
    const string a = std::to_string(allele);
    return strand == 0 ? a + "|." : "." + ("|" + a);
}

int VCFOutputCaller::emit_block_records(const PathPositionHandleGraph& graph, const Snarl& snarl,
                                        const vector<SnarlTraversal>& called_traversals,
                                        const vector<int>& genotype, int ref_trav_idx,
                                        const string& sample_name, const vcflib::Variant& site,
                                        const map<int, int>& trav_to_allele,
                                        int64_t site_position, GLLayout gl_layout,
                                        bool genotype_snarls) const {
    // Every refusal below returns -1, meaning "the site record stands". That is the safe direction:
    // -1 reproduces today's output exactly, so a case this does not understand degrades to the
    // behaviour it was going to have anyway rather than to a wrong record.
    if (!atomize_blocks || symbolic_manager == nullptr || genotype_snarls || genotype.empty()) {
        return -1;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= called_traversals.size()) {
        return -1;
    }
    if (!symbolic_site_resolvable(snarl, *symbolic_manager)) {
        // Projection would see no child chains here, so a diff over it measures node-level
        // shredding rather than route structure. Counted in the stage-2 instrumentation.
        return -1;
    }

    const SnarlTraversal& ref_trav = called_traversals[ref_trav_idx];
    vector<pair<int, int>> ref_ranges;
    SymbolicAllele sref = symbolic_allele(ref_trav, snarl, *symbolic_manager, &ref_ranges);
    const size_t m = sref.size();
    if (m == 0 || ref_ranges.size() != m) {
        return -1;
    }

    // Base offset of every visit boundary of the reference traversal, measured from the snarl's
    // first base. The reference traversal is consecutive reference-path steps in path order (the
    // caller asserts that where it builds it), so a cumulative sum of node lengths IS the offset.
    vector<size_t> ref_visit_off(ref_trav.visit_size() + 1, 0);
    for (int v = 0; v < ref_trav.visit_size(); ++v) {
        size_t len = 0;
        if (ref_trav.visit(v).node_id() > 0) {
            len = graph.get_length(graph.get_handle(ref_trav.visit(v).node_id()));
        }
        ref_visit_off[v + 1] = ref_visit_off[v] + len;
    }

    auto visit_of_step = [](const vector<pair<int, int>>& ranges, size_t step,
                            const SnarlTraversal& t) -> int {
        // The ranges partition the visits contiguously, so the visit index at step boundary k is
        // ranges[k].first, and one past the end is the traversal length.
        return step < ranges.size() ? ranges[step].first : t.visit_size();
    };
    // `max(vb, 0)`, not `vb`: the anchor path below asks for [vb - 1, vb), so a caller that let
    // vb == 0 through would read visit(-1). The `vb <= 0` refusal is what stops that today, and it
    // stays -- attributing the decline to the right reason -- but a helper should not depend on one
    // caller's bounds check to be memory-safe.
    auto seq_of = [&](const SnarlTraversal& t, int vb, int ve) -> string {
        string s;
        for (int v = std::max(vb, 0); v < ve && v < t.visit_size(); ++v) {
            const Visit& vis = t.visit(v);
            if (vis.node_id() > 0) {
                s += graph.get_sequence(graph.get_handle(vis.node_id(), vis.backward()));
            }
        }
        return s;
    };

    struct HapAlign {
        int trav = -1;
        SymbolicAllele sym;
        vector<pair<int, int>> ranges;
        vector<DiffBlock> blocks;
        vector<int> alt_before_ref;
    };
    vector<HapAlign> haps(genotype.size());
    for (size_t s = 0; s < genotype.size(); ++s) {
        haps[s].trav = genotype[s];
        if (genotype[s] < 0 || (size_t)genotype[s] >= called_traversals.size()) {
            continue;   // a star or missing allele: nothing to align
        }
        if (genotype[s] == ref_trav_idx) {
            continue;   // the reference itself: every step matches, so no blocks
        }
        haps[s].sym = symbolic_allele(called_traversals[genotype[s]], snarl, *symbolic_manager,
                                      &haps[s].ranges);
        bool degraded = false;
        haps[s].blocks = symbolic_diff(sref, haps[s].sym, &degraded, &haps[s].alt_before_ref);
        if (degraded || haps[s].alt_before_ref.size() != m + 1) {
            return -1;
        }
    }

    // Cluster every haplotype's blocks by overlap of reference step range, treating touching ranges
    // as overlapping. That is what merges a deletion on one haplotype abutting an insertion on the
    // other: their ranges share only an endpoint, but they cannot be reported as separate records
    // because the two alleles would have to disagree about the same reference span.
    vector<pair<int, int>> ivs;
    for (const HapAlign& h : haps) {
        for (const DiffBlock& b : h.blocks) {
            ivs.emplace_back(b.ref_begin, b.ref_end);
        }
    }
    if (ivs.empty()) {
        return -1;
    }
    sort(ivs.begin(), ivs.end());
    vector<pair<int, int>> clusters;
    clusters.push_back(ivs[0]);
    for (size_t k = 1; k < ivs.size(); ++k) {
        if (ivs[k].first <= clusters.back().second) {
            clusters.back().second = std::max(clusters.back().second, ivs[k].second);
        } else {
            clusters.push_back(ivs[k]);
        }
    }

    // Build every record before writing any, so the decision to split can be taken on the finished
    // set. A half-written split is the one outcome with no safe fallback.
    struct Built {
        vcflib::Variant var;
        bool wanted = false;
    };
    vector<Built> built;
    built.reserve(clusters.size());

    for (const pair<int, int>& cluster : clusters) {
        const size_t rb = (size_t)cluster.first;
        const size_t re = (size_t)cluster.second;
        int vb = visit_of_step(ref_ranges, rb, ref_trav);
        int ve = visit_of_step(ref_ranges, re, ref_trav);
        string ref_str = seq_of(ref_trav, vb, ve);

        // Each haplotype's allele over this same reference span. A haplotype with no block here is
        // matched throughout, so its span is the aligned one and its string equals the reference's.
        vector<string> slot_str(genotype.size());
        vector<bool> slot_marker(genotype.size(), false);
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (haps[s].trav < 0) {
                slot_marker[s] = true;
                continue;
            }
            if (haps[s].trav == ref_trav_idx) {
                slot_str[s] = ref_str;
                continue;
            }
            int ab = haps[s].alt_before_ref[rb];
            int ae = haps[s].alt_before_ref[re];
            for (const DiffBlock& b : haps[s].blocks) {
                if ((size_t)b.ref_begin <= re && rb <= (size_t)b.ref_end) {
                    ab = std::min(ab, b.alt_begin);
                    ae = std::max(ae, b.alt_end);
                }
            }
            const SnarlTraversal& t = called_traversals[haps[s].trav];
            slot_str[s] = seq_of(t, visit_of_step(haps[s].ranges, (size_t)ab, t),
                                 visit_of_step(haps[s].ranges, (size_t)ae, t));
        }

        // VCF has no representation for an empty allele, so an indel borrows the base before it.
        // One base, which is what flatten_common_allele_ends leaves on every indel it touches, so
        // the two agree about how an indel is spelled.
        bool needs_anchor = ref_str.empty();
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (!slot_marker[s] && slot_str[s].empty()) {
                needs_anchor = true;
            }
        }
        int64_t pos = site_position + (int64_t)ref_visit_off[vb];
        if (needs_anchor) {
            if (vb <= 0) {
                return -1;   // no base to the left inside this snarl; the site record stands
            }
            string left = seq_of(ref_trav, vb - 1, vb);
            if (left.empty()) {
                return -1;
            }
            string base(1, left.back());
            ref_str = base + ref_str;
            for (size_t s = 0; s < genotype.size(); ++s) {
                if (!slot_marker[s]) {
                    slot_str[s] = base + slot_str[s];
                }
            }
            pos -= 1;
        }

        // Dedup by string, exactly as the site record does, so two haplotypes taking different
        // routes to the same sequence come out homozygous rather than as two identical ALTs.
        map<string, int> allele_to_gt;
        vector<string> alleles;
        vector<int> site_of_block;
        allele_to_gt[ref_str] = 0;
        alleles.push_back(ref_str);
        site_of_block.push_back(0);
        vector<int> block_gt(genotype.size(), MISSING_ALLELE_MARKER);
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (slot_marker[s]) {
                continue;
            }
            auto found = allele_to_gt.find(slot_str[s]);
            if (found != allele_to_gt.end()) {
                block_gt[s] = found->second;
                continue;
            }
            int a = (int)alleles.size();
            allele_to_gt[slot_str[s]] = a;
            alleles.push_back(slot_str[s]);
            // Which site allele this block allele inherits its evidence from: the one the same
            // haplotype carries at the site. That is the whole of the replication decision.
            auto sa = trav_to_allele.find(haps[s].trav);
            site_of_block.push_back(sa != trav_to_allele.end() ? sa->second : 0);
            block_gt[s] = a;
        }
        if (alleles.size() < 2) {
            continue;   // every haplotype is the reference here; nothing to report
        }

        Built b;
        b.var = site;   // inherit INFO, FORMAT and every site-level field
        b.var.position = pos;
        b.var.ref = alleles[0];
        b.var.alt.assign(alleles.begin() + 1, alleles.end());
        b.var.alleles = alleles;
        b.var.info["AT"].clear();
        b.var.info["AT"].resize(alleles.size());
        b.var.info["SB"] = {std::to_string(built.size()), std::to_string(clusters.size())};

        // AT per block allele, over the visit range this record actually spells.
        {
            SnarlTraversal ref_span;
            for (int v = vb; v < ve && v < ref_trav.visit_size(); ++v) {
                *ref_span.add_visit() = ref_trav.visit(v);
            }
            add_allele_path_to_info(b.var, 0, ref_span, false, false);
            for (size_t a = 1; a < alleles.size(); ++a) {
                SnarlTraversal span;
                for (size_t s = 0; s < genotype.size(); ++s) {
                    if (slot_marker[s] || block_gt[s] != (int)a) {
                        continue;
                    }
                    const SnarlTraversal& t = called_traversals[haps[s].trav];
                    int ab = haps[s].trav == ref_trav_idx
                                 ? vb : visit_of_step(haps[s].ranges,
                                                      (size_t)haps[s].alt_before_ref[rb], t);
                    int ae = haps[s].trav == ref_trav_idx
                                 ? ve : visit_of_step(haps[s].ranges,
                                                      (size_t)haps[s].alt_before_ref[re], t);
                    for (const DiffBlock& blk : haps[s].blocks) {
                        if ((size_t)blk.ref_begin <= re && rb <= (size_t)blk.ref_end) {
                            ab = std::min(ab, visit_of_step(haps[s].ranges,
                                                            (size_t)blk.alt_begin, t));
                            ae = std::max(ae, visit_of_step(haps[s].ranges,
                                                            (size_t)blk.alt_end, t));
                        }
                    }
                    for (int v = ab; v < ae && v < t.visit_size(); ++v) {
                        *span.add_visit() = t.visit(v);
                    }
                    break;
                }
                add_allele_path_to_info(b.var, a, span, false, false);
            }
        }

        // GT, inheriting the site's phase rather than dropping it.
        //
        // The site's GT is in phased order when it carries a separator, and that order is in SITE
        // allele space; this record's slots are in genotyper order. So the two are compared by
        // mapping each slot to the site allele its haplotype carries. Dropping the phase instead
        // would leave a record carrying PS beside an unphased GT -- a contradiction on its face --
        // and would score these records against a phased baseline on unequal terms, which is a
        // confound rather than a conservative choice.
        //
        // THREE shapes carry a phase set, and two of them are haploid, which is what this used to
        // get wrong by testing for a diploid pair before looking at anything else:
        //
        //   - "a|b" is a phased diploid pair, and the orientation transfers by slot.
        //   - "a|." or ".|a" is a nested chain called at ploidy 1: one strand of a diploid locus,
        //     where the parent's other allele deletes the chain. Which strand is the only thing
        //     this form carries that a bare "a" does not, and it is the reason the strand reaches
        //     the VCF at all instead of living only in the mosaic. A block of such a call is a
        //     piece of the same allele on the same strand, so the strand transfers unchanged.
        //   - "a" is a genuinely haploid locus (chrY, chrX outside the PARs). No orientation to
        //     inherit, but PS is a block label and still applies.
        //
        // Only a slash-separated GT is unphased, and only there does PS mean nothing.
        bool keep_phase_set = false;
        int nested_strand = -1;   // haploid record: which side of "a|." this allele sits on
        vector<int> slot_order(block_gt.size());
        for (size_t s = 0; s < block_gt.size(); ++s) {
            slot_order[s] = (int)s;
        }
        const string* site_gt_text = nullptr;
        {
            auto site_gt = site.samples.find(sample_name);
            if (site_gt != site.samples.end()) {
                auto gt_field = site_gt->second.find("GT");
                if (gt_field != site_gt->second.end() && !gt_field->second.empty()) {
                    site_gt_text = &gt_field->second[0];
                }
            }
        }
        if (site_gt_text != nullptr && site_gt_text->find('/') == string::npos) {
            const size_t bar = site_gt_text->find('|');
            if (bar == string::npos) {
                keep_phase_set = block_gt.size() == 1 && block_gt[0] >= 0;
            } else {
                const string left = site_gt_text->substr(0, bar);
                const string right = site_gt_text->substr(bar + 1);
                if (block_gt.size() == 1 && block_gt[0] >= 0 && (left == "." || right == ".")) {
                    keep_phase_set = true;
                    nested_strand = left == "." ? 1 : 0;
                } else if (block_gt.size() == 2 && block_gt[0] >= 0 && block_gt[1] >= 0) {
                    int a = -1;
                    int b2 = -1;
                    try {
                        a = std::stoi(left);
                        b2 = std::stoi(right);
                    } catch (const std::exception&) {
                        a = -1;
                    }
                    auto site_allele_of_slot = [&](size_t sl) {
                        auto found = trav_to_allele.find(haps[sl].trav);
                        return found != trav_to_allele.end() ? found->second : -1;
                    };
                    int s0 = site_allele_of_slot(0);
                    int s1 = site_allele_of_slot(1);
                    if (a >= 0 && s0 >= 0 && s1 >= 0) {
                        if (s0 == a && s1 == b2) {
                            keep_phase_set = true;
                        } else if (s0 == b2 && s1 == a) {
                            keep_phase_set = true;
                            slot_order[0] = 1;
                            slot_order[1] = 0;
                        }
                    }
                }
            }
        }
        {
            string gt;
            if (nested_strand >= 0) {
                gt = nested_strand_genotype(block_gt[0], nested_strand);
            } else {
                for (size_t s = 0; s < block_gt.size(); ++s) {
                    const int value = block_gt[slot_order[s]];
                    gt += value == MISSING_ALLELE_MARKER ? string(".") : std::to_string(value);
                    if (s + 1 != block_gt.size()) {
                        gt += keep_phase_set ? '|' : '/';
                    }
                }
            }
            b.var.samples[sample_name]["GT"] = {gt};
        }
        if (!keep_phase_set) {
            // PS labels a phase block, so it means nothing on an unphased genotype.
            b.var.samples[sample_name].erase("PS");
            b.var.format.erase(std::remove(b.var.format.begin(), b.var.format.end(), "PS"),
                               b.var.format.end());
        }

        // AD and GL are looked up through site_of_block rather than recomputed. Every block of a
        // snarl therefore reports the same evidence, which is honest only about arity: see D2 in
        // planning/symbolic-diff-decomposition.md. INFO/SB is what makes the replicated set
        // recoverable by a consumer that must not double-count it.
        auto& fmt = b.var.samples[sample_name];
        auto site_fmt = site.samples.find(sample_name);
        if (site_fmt != site.samples.end()) {
            auto site_ad = site_fmt->second.find("AD");
            if (site_ad != site_fmt->second.end()) {
                vector<string> ad;
                for (size_t a = 0; a < alleles.size(); ++a) {
                    size_t si = (size_t)site_of_block[a];
                    ad.push_back(si < site_ad->second.size() ? site_ad->second[si] : "0");
                }
                fmt["AD"] = ad;
            }
            auto site_gl = site_fmt->second.find("GL");
            bool has_marker = std::any_of(block_gt.begin(), block_gt.end(),
                                          [](int a) { return a < 0; });
            if (site_gl != site_fmt->second.end() && !has_marker) {
                const size_t k = alleles.size();
                const size_t ks = site.alleles.size();
                vector<string> gl;
                bool ok = true;
                if (block_gt.size() == 1) {
                    gl.resize(k);
                    for (size_t a = 0; a < k && ok; ++a) {
                        size_t si = (size_t)site_of_block[a];
                        ok = si < site_gl->second.size();
                        if (ok) {
                            gl[a] = site_gl->second[si];
                        }
                    }
                } else if (block_gt.size() == 2) {
                    gl.resize(k * (k + 1) / 2);
                    for (size_t j = 0; j < k && ok; ++j) {
                        for (size_t i = 0; i <= j && ok; ++i) {
                            size_t si = (size_t)site_of_block[i];
                            size_t sj = (size_t)site_of_block[j];
                            size_t src = gl_genotype_index(std::min(si, sj), std::max(si, sj), ks,
                                                           gl_layout);
                            size_t dst = gl_genotype_index(i, j, k, gl_layout);
                            ok = src < site_gl->second.size() && dst < gl.size();
                            if (ok) {
                                gl[dst] = site_gl->second[src];
                            }
                        }
                    }
                } else {
                    ok = false;
                }
                if (ok) {
                    fmt["GL"] = gl;
                } else {
                    fmt.erase("GL");
                    b.var.format.erase(std::remove(b.var.format.begin(), b.var.format.end(), "GL"),
                                       b.var.format.end());
                }
            }
        }

        b.var.updateAlleleIndexes();
        flatten_common_allele_ends(b.var, true, 0);
        flatten_common_allele_ends(b.var, false, 0);
        b.wanted = !b.var.alt.empty();
        built.push_back(std::move(b));
    }

    size_t wanted = 0;
    for (const Built& b : built) {
        if (b.wanted) {
            ++wanted;
        }
    }
    if (wanted == 0) {
        return -1;
    }
    // Only take over from the site record where doing so changes the answer: more than one record,
    // or one record whose allele set is smaller than the site's -- which is the two-haplotypes-
    // one-route case coming out homozygous instead of as two near-identical ALTs. Everything else
    // falls through so its bytes are unchanged.
    bool collapses = wanted == 1 && built.size() == 1 &&
                     built[0].var.alleles.size() < site.alleles.size();
    if (wanted < 2 && !collapses) {
        return -1;
    }

    int added = 0;
    size_t block_index = 0;
    for (Built& b : built) {
        if (!b.wanted) {
            continue;
        }
        if (add_variant(b.var, block_index)) {
            ++added;
        }
        ++block_index;
    }
    return added;
}

bool VCFOutputCaller::emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                                   const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                                   const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                   const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                                   function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string) {
    
#ifdef debug
    cerr << "emitting variant for " << pb2json(snarl) << endl;
    for (int i = 0; i < called_traversals.size(); ++i) {
        if (i == ref_trav_idx) {
            cerr << "*";
        }
        cerr << "ct[" << i << "]=" << pb2json(called_traversals[i]) << endl;
    }
    for (int i = 0; i < genotype.size(); ++i) {
        cerr << "gt[" << i << "]=" << genotype[i] << endl;
    }
#endif

    // Stale from the previous emit on this thread until this one fills it in. Cleared rather than
    // left, so a descent after an emit that wrote nothing cannot read the last snarl's mapping.
    last_emitted.valid = false;
    last_emitted.num_alleles = 0;
    last_emitted.trav_to_allele.clear();
    last_emitted.contig.clear();
    last_emitted.position = 0;
    last_emitted.buffer_thread = -1;
    last_emitted.buffer_index = 0;

    if (trav_to_string == nullptr) {
        trav_to_string = [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {
            return trav_string(graph, travs[trav_allele]);    
        };
    }

    vcflib::Variant out_variant;

    vector<SnarlTraversal> site_traversals = {called_traversals[ref_trav_idx]};
    vector<int> site_genotype;
    auto ref_gt_it = std::find(genotype.begin(), genotype.end(), ref_trav_idx);
    out_variant.ref = trav_to_string(called_traversals, genotype, ref_trav_idx,
                                     ref_gt_it != genotype.end() ? ref_gt_it - genotype.begin() : 0,
                                     ref_trav_idx);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;
    // Which VCF allele each *called traversal* became. The two numberings differ -- alleles are
    // deduplicated by string, and only called traversals reach the record at all -- and the
    // linkage pass needs the VCF numbering, because that is what GL and GT are written in.
    map<int, int> trav_to_allele;
    allele_to_gt[out_variant.ref] = 0;
    trav_to_allele[ref_trav_idx] = 0;
    int star_allele_idx = -1;  // index for star allele in allele_to_gt, if needed
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == STAR_ALLELE_MARKER) {
            // Star allele: haplotype spans this site but has no defined traversal here
            if (star_allele_idx < 0) {
                // Add star allele to allele list
                star_allele_idx = allele_to_gt.size();
                allele_to_gt["*"] = star_allele_idx;
                // Add empty traversal as placeholder (won't be used for AT info)
                site_traversals.push_back(SnarlTraversal());
            }
            site_genotype.push_back(star_allele_idx);
        } else if (genotype[i] == MISSING_ALLELE_MARKER) {
            // Missing allele: parent doesn't traverse this child, output as '.' in VCF
            site_genotype.push_back(MISSING_ALLELE_MARKER);
        } else if (genotype[i] == ref_trav_idx || is_symbolically_reference(called_traversals,
                                                                            genotype[i], ref_trav_idx,
                                                                            snarl)) {
            // Either literally the reference traversal, or one that takes the same route through
            // this snarl and differs only inside child chains. The second kind is the reference
            // allele *here*; its differences belong to the nested sites that contain them, and
            // emitting it as a long ALT is what buries them.
            site_genotype.push_back(0);
            trav_to_allele[genotype[i]] = 0;
        } else {
            string allele_string = trav_to_string(called_traversals, genotype, genotype[i], i, ref_trav_idx);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_traversals.push_back(called_traversals[genotype[i]]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
            trav_to_allele[genotype[i]] = site_genotype.back();
        }
    }

    tally_atomize(graph, symbolic_manager, snarl, called_traversals, genotype, ref_trav_idx);

    // add on fixed number of uncalled traversals if we're making a ref-call
    // with genotype_snarls set to true
    if (genotype_snarls && site_traversals.size() <= 1) {
        // note: we're adding all the strings here and sorting to make this deterministic
        // at the cost of speed
        map<string, const SnarlTraversal*> allele_map;
        for (int i = 0; i < called_traversals.size(); ++i) {
            // todo: verify index below.  it's for uncalled traversals so not important tho
            string allele_string = trav_to_string(called_traversals, genotype, i, max(0, (int)genotype.size() - 1), ref_trav_idx);
            if (!allele_map.count(allele_string)) {
                allele_map[allele_string] = &called_traversals[i];
            }
        }
        // pick out the first "max_uncalled_alleles" traversals to add
        int i = 0;
        for (auto ai = allele_map.begin(); i < max_uncalled_alleles && ai != allele_map.end(); ++i, ++ai) {
            if (!allele_to_gt.count(ai->first)) {
                allele_to_gt[ai->first] = allele_to_gt.size();
                site_traversals.push_back(*ai->second);
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    
    // init the traversal info
    out_variant.info["AT"].resize(allele_to_gt.size());

    for (auto& allele_gt : allele_to_gt) {
#ifdef debug
        cerr << "allele " << allele_gt.first << " -> gt " << allele_gt.second << endl;
#endif
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;

        // update the traversal info
        add_allele_path_to_info(out_variant, allele_gt.second, site_traversals.at(allele_gt.second), false, false); 
    }

    // resolve subpath naming
    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    // in VCF we usually just want a contig
    string contig_name = PathMetadata::parse_locus_name(basepath_name);
    if (contig_name != PathMetadata::NO_LOCUS_NAME) {
        basepath_name = contig_name;
    }
    // fill out the rest of the variant    
    out_variant.sequenceName = basepath_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_offset + 1 + basepath_offset;
    // Kept before flattening moves it. This is the position of the reference traversal's FIRST
    // base, which is what per-block offsets are measured from; the flattened value is not, because
    // the shared prefix it removes depends on the emitted allele set.
    const int64_t site_position_unflattened = out_variant.position;
    out_variant.id = print_snarl(snarl, false);
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    
    stringstream vcf_gt;
    if (!genotype.empty()) {
        for (int i = 0; i < site_genotype.size(); ++i) {
            if (site_genotype[i] == MISSING_ALLELE_MARKER) {
                vcf_gt << ".";
            } else {
                vcf_gt << site_genotype[i];
            }
            if (i != site_genotype.size() - 1) {
                vcf_gt << "/";
            }
        }
    } else {
        for (int i = 0; i < ploidy; ++i) {
            vcf_gt << ".";
            if (i != ploidy - 1) {
                vcf_gt << "/";
            }
        }
    }
                    
    genotype_vector.push_back(vcf_gt.str());

    int64_t phase_set_to_write = -1;

    // Phased here, while the traversal-to-allele mapping this record was just built with is still in
    // hand, rather than by rewriting the line afterwards.
    //
    // That mapping is the whole reason phasing used to be a patch. A PhaseCall names a *traversal*
    // pair; the VCF needs allele numbers; and the map between them did not exist until the record was
    // built. So the phase was resolved against whatever numbering was available at resolution time --
    // which under decide-then-render is none at all, making `render_phase_pair` fall back on every one
    // of chr20's 192,045 top-level sites -- and then patched in against the line. Here the map is
    // `trav_to_allele`, complete and correct, a few lines above.
    if (!genotype.empty() && emit_phasing && !render_phases.empty()) {
        auto found = render_phases.find(record_key_of(snarl));
        if (found != render_phases.end()) {
            const LinkageCollector::PhaseCall& phase = found->second;
            // `find`, not `operator[]` and not a bounds check against `size()`. This is a
            // std::map keyed BY traversal, so its size is the number of alleles the record carries,
            // not a bound on traversal indices -- comparing against it refused 101,947 phases, and
            // `operator[]` on a miss would have inserted a default 0, writing allele 0 for a
            // traversal the record does not carry into the very map that `set_allele_map` is then
            // handed.
            const auto found_a = trav_to_allele.find(phase.trav_first);
            const auto found_b = trav_to_allele.find(phase.trav_second);
            const int a = (phase.trav_first >= 0 && found_a != trav_to_allele.end())
                              ? found_a->second : -1;
            const int b = (phase.trav_second >= 0 && found_b != trav_to_allele.end())
                              ? found_b->second : -1;
            // The phased genotype must be a permutation of the one this record carries. Checked, not
            // assumed: phasing that silently substituted a genotype would be invisible in the output,
            // and two records can share a position.
            bool same = false;
            if (phase.ploidy == 1 && site_genotype.size() == 1) {
                same = (a >= 0 && a == site_genotype[0]);
            } else if (phase.ploidy == 2 && site_genotype.size() == 2) {
                same = (a >= 0 && b >= 0)
                       && ((a == site_genotype[0] && b == site_genotype[1])
                           || (a == site_genotype[1] && b == site_genotype[0]));
            }

            if (same) {
                if (phase.ploidy == 1 && phase.nested_strand >= 0) {
                    // A nested haploid site is one strand of a diploid locus, not a haploid locus: it
                    // is called at ploidy 1 because the parent's *other* allele deletes the chain, so
                    // there is no sequence on that strand to genotype. Written as a phased pair with
                    // the empty strand as ".", the only place the VCF can carry which strand the
                    // allele is on -- a bare "a" names none, which is why the strand lived only in the
                    // mosaic and no phasing tool could read it.
                    genotype_vector[0] = nested_strand_genotype(a, phase.nested_strand);
                } else if (phase.ploidy == 1) {
                    // A genuinely haploid locus. One allele, no phase; only PS is meaningful, and
                    // only as a block label. "a|a" would claim a homozygous diploid call.
                    genotype_vector[0] = std::to_string(a);
                } else {
                    genotype_vector[0] = std::to_string(a) + "|" + std::to_string(b);
                }
                // PS is added after update_vcf_info below, not here: adding it now would put it
                // second in FORMAT, ahead of DP/AD/GL, where the patch that used to append it left
                // it last. Same fields and same values either way, but every line differs, which
                // costs the byte comparison against the previous arm for no gain.
                phase_set_to_write = (int64_t)phase.phase_set;
                ++phased_records;
            } else {
                ++phase_declined;
            }
        }
    }

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name, out_variant);

    // PS last in FORMAT, which is where appending it as a patch used to leave it.
    if (phase_set_to_write >= 0) {
        out_variant.format.push_back("PS");
        out_variant.samples[sample_name]["PS"].push_back(std::to_string(phase_set_to_write));
    }

    // if genotype_snarls, then we only flatten up to the snarl endpoints
    // (this is when we are in genotyping mode and want consistent calls regardless of the sample)
    int64_t flatten_len_s = 0;
    int64_t flatten_len_e = 0;
    if (genotype_snarls) {
        flatten_len_s = graph.get_length(graph.get_handle(snarl.start().node_id()));
        assert(flatten_len_s >= 0);
        flatten_len_e = graph.get_length(graph.get_handle(snarl.end().node_id()));
    }
    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true, flatten_len_e);
    flatten_common_allele_ends(out_variant, false, flatten_len_s);

    // Merge near-identical called ALT alleles (vg call -L), turning 1/2 into 1/1.  Placed here on
    // purpose: after update_vcf_info so the genotyper saw every candidate, and after flattening so
    // the surviving allele's string, POS and REF are byte-identical to a run without -L (a shorter
    // allele list can share a longer prefix and flatten further).  The missing-allele fixup below
    // still runs after it, and is unaffected: merging is ALT-vs-ALT so it never empties alt.
    // Which order this record's GL is in follows from which caller wrote it, and nothing in the
    // record says so -- hence passing it rather than sniffing it.
    const GLLayout gl_layout =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get())
            != nullptr
            ? GLLayout::Colexicographic
            : GLLayout::IMajor;
    merge_similar_alleles(graph, site_traversals, site_genotype, sample_name, out_variant,
                          gl_layout);
#ifdef debug
    for (int i = 0; i < site_traversals.size(); ++i) {
        cerr << " site trav[" << i << "]=" << pb2json(site_traversals[i]) << endl;
    }
    for (int i = 0; i < site_genotype.size(); ++i) {
        cerr << " site geno[" << i << "]=" << site_genotype[i] << endl;
    }
#endif

    // If genotype contains missing allele but no ALT, add * as ALT to emit valid VCF
    // This happens when one parent haplotype doesn't traverse a nested child snarl
    bool has_missing = std::find(site_genotype.begin(), site_genotype.end(), MISSING_ALLELE_MARKER) != site_genotype.end();
    if (has_missing && out_variant.alt.empty()) {
        out_variant.alt.push_back("*");
        out_variant.alleles.push_back("*");
        out_variant.info["AT"].push_back(".");
        // AD is Number=R and update_vcf_info wrote it above, over the allele set as it stood --
        // one entry, for the reference. Appending an allele without appending its AD entry leaves
        // the record contradicting its own header, which is what it did before this line existed.
        // GL needs no matching fixup: it is already omitted whenever the genotype carries a
        // marker, and MISSING_ALLELE_MARKER is one (read_likelihood_caller.cpp, genotype_has_marker).
        auto ad_it = out_variant.samples[sample_name].find("AD");
        if (ad_it != out_variant.samples[sample_name].end()) {
            ad_it->second.push_back("0");
        }
    }

    // Hand the traversal-to-allele mapping to symbolic descent, which runs next on this thread and
    // needs a child's crossing pattern in allele space. Filled whether or not the record survives
    // add_variant: a parent that collapses to the reference emits nothing and still has children to
    // descend into, which is the case nested calling exists for.
    if (symbolic_manager != nullptr) {
        last_emitted.valid = true;
        last_emitted.num_alleles = out_variant.alleles.size();
        last_emitted.trav_to_allele.clear();
        for (const auto& kv : trav_to_allele) {
            last_emitted.trav_to_allele.emplace(kv.first, kv.second);
        }
        // Post-flatten, which is the point: this is the (contig, POS) add_variant keys the line
        // under, so a linkage entry recorded from these is one write_variants can actually find.
        // get_ref_position reproduces the pre-flatten value and must not be used for this.
        last_emitted.contig = out_variant.sequenceName;
        last_emitted.position = out_variant.position;
    }

    // One record per difference block, where that changes the answer. Placed here on purpose: the
    // site record above is finished, so every field a block does not redefine is inherited from a
    // record that already passed update_vcf_info, flattening and merging. A return of -1 means this
    // site was declined, and then the site record below is written exactly as it always was.
    const int block_lines = emit_block_records(graph, snarl, called_traversals, genotype,
                                               ref_trav_idx, sample_name, out_variant,
                                               trav_to_allele, site_position_unflattened,
                                               gl_layout, genotype_snarls);
    if (block_lines >= 0) {
        ++g_atomize_split_sites;
        g_atomize_split_lines += (size_t)block_lines;
        // last_emitted above still describes this snarl's traversals, which is what descent reads;
        // the buffer handle is deliberately left unset, because there is no single line to retract.
        return block_lines > 0;
    }

    // Whether this site wants a line at all. A traversal pair differing from the reference only
    // inside child chains collapses to allele 0 here and leaves `alt` empty -- that is what
    // symbolic collapsing means -- and such a site is exactly the one whose children need it most.
    const bool wants_line = genotype_snarls || !out_variant.alt.empty();
    bool added = false;
    if (wants_line) {
        added = add_variant(out_variant);
        if (added) {
            // Exactly which buffered line this record became, so the barrier can retract or
            // replace it without re-deriving its identity from the text.
            last_emitted.buffer_thread = omp_get_thread_num();
            last_emitted.buffer_index = output_variants[omp_get_thread_num()].size() - 1;
        }
    }
    // The linkage layer is fed whether or not a line exists, which is the whole of stage 2.
    //
    // A parent that collapses to the reference still HAS two alleles; they differ only inside its
    // children, which is precisely the information those children need to know which of its
    // haplotypes carries the chain. Entering it only when it wrote a line meant it was absent from
    // the model entirely, so 289 of chr20's 292 strandless haploid records had no phased parent to
    // inherit from -- and nothing distinguished that from a strand that was genuinely undecidable.
    //
    // Recording it under the *emitted* allele numbering would not have worked and is worth saying
    // so: collapsing maps every one of its called traversals to allele 0, so it would enter as
    // 0/0, and a homozygous-reference site has no strand distinction to inherit. It is only in
    // traversal space that it is a real heterozygous site.
    if (linkage_collector != nullptr && !suppress_linkage_record) {
        // The site was recorded at the genotyping site, before any of this ran. All that is left is
        // the traversal-to-VCF-allele map, which is a function of the allele list chosen just above
        // and so cannot exist before the record is built, and whether a line was actually written.
        //
        // Both exist only to serve the patch path: the map is what renders a settled compact allele
        // into a number a GT can name, and the flag is what says there is a line to patch at all.
        // Stage 11 deletes patching and both go with it.
        vector<int> trav_to_allele_vec(called_traversals.size(), -1);
        for (const auto& kv : trav_to_allele) {
            if (kv.first >= 0 && (size_t)kv.first < trav_to_allele_vec.size()) {
                trav_to_allele_vec[kv.first] = kv.second;
            }
        }
        linkage_collector->set_allele_map(record_key_of(snarl), trav_to_allele_vec,
                                          added);
    }
    if (wants_line && !added) {
        stringstream ss;
        ss << out_variant;
        cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":" << out_variant.position
             << " with ID=" << out_variant.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
             << VCFOutputCaller::max_vcf_line_length << endl;
    }
    // Unchanged contract, and the barrier depends on the distinction: true when the record
    // flattened to nothing, false only when add_variant refused a line it wanted to write.
    return wants_line ? added : true;
}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> VCFOutputCaller::get_ref_interval(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) const {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    map<size_t, step_handle_t> start_steps;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                start_steps[graph.get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    map<size_t, step_handle_t> end_steps;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                end_steps[graph.get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    // just because we found a pair of steps on our path that correspond to the snarl ends, doesn't
    // mean the path threads the snarl.  verify that we can actaully walk, either forwards or backwards
    // along the path from the start node and hit then end node in the right orientation. 
    bool start_rev = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
    bool end_rev = graph.get_is_reverse(graph.get_handle_of_step(end_step)) != snarl.end().backward();
    bool found_end = start_rev == end_rev && start_rev == start_steps.begin()->first > end_steps.begin()->first;
        
    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        found_end = false;
        // try each start step
        for (auto i = start_steps.begin(); i != start_steps.end() && !found_end; ++i) {
            start_step = i->second;
            bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
            if (scan_backward) {
                // if we're going backward, we expect to reach the end backward
                end_handle = graph.get_handle(snarl.end().node_id(), !snarl.end().backward());
            }            
            if (scan_backward) {
                for (step_handle_t cur_step = start_step; graph.has_previous_step(cur_step) && !found_end;
                     cur_step = graph.get_previous_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            } else {
                for (step_handle_t cur_step = start_step; graph.has_next_step(cur_step) && !found_end;
                     cur_step = graph.get_next_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            }
        }
    }
    int64_t start_position = start_steps.begin()->first;
    step_handle_t out_start_step = start_step;
    int64_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : graph.get_position_of_step(end_step);
    step_handle_t out_end_step = end_step == end_steps.begin()->second ? end_steps.begin()->second : end_step;
    bool backward = end_position < start_position;
    

    if (!found_end) {
        // oops, once of the above checks failed.  we tell caller we coudlnt find by hacking in a -1 coordinate.
        start_position = -1;
        end_position = -1;
    }

    if (backward) {
        return make_tuple(end_position, start_position, backward, out_end_step, out_start_step);
    } else {
        return make_tuple(start_position, end_position, backward, out_start_step, out_end_step);
    }
}

pair<string, int64_t> VCFOutputCaller::get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                                        int64_t ref_path_offset) const {

    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    // +1 to convert to 1-based VCF
    int64_t position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_path_offset + 1 + basepath_offset;
    return make_pair(basepath_name, position);
}

void VCFOutputCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const {
    if (variant.alt.size() == 0) {
        return;
    }

    // find the minimum allele length to make sure we don't delete an entire allele
    size_t min_allele_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_allele_len = std::min(min_allele_len, variant.alleles[i].length());
    }

    // An empty allele makes the arithmetic below underflow. max_flatten_len is size_t, and
    // min_allele_len == 0 turns the "leave one base in the reference" decrement into SIZE_MAX,
    // after which the backward pass evaluates alleles[j][length() - 1 - i] on an empty string --
    // an out-of-bounds read, not a wrong answer. The early return above does not save us: a pure
    // deletion has one ALT, which is empty. There is nothing to zip up in this case anyway.
    if (min_allele_len == 0) {
        return;
    }

    // the maximum number of bases we want ot zip up, applying override if provided
    size_t max_flatten_len = len_override > 0 ? len_override : min_allele_len;
    
    // want to leave at least one in the reference position
    if (max_flatten_len == min_allele_len) {
        --max_flatten_len;
    }
    
    bool match = true;
    int shared_prefix_len = 0;
    for (int i = 0; i < max_flatten_len && match; ++i) {
        char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
        for (int j = 1; j < variant.alleles.size() && match; ++j) {
            char c2 = std::toupper(variant.alleles[j][!backward ? i : variant.alleles[j].length() - 1 - i]);
            match = c1 == c2;
        }
        if (match) {
            ++shared_prefix_len;
        }
    }

    if (!backward) {
        variant.position += shared_prefix_len;
    }
    for (int i = 0; i < variant.alleles.size(); ++i) {
        if (!backward) {
            variant.alleles[i] = variant.alleles[i].substr(shared_prefix_len);
        } else {
            variant.alleles[i] = variant.alleles[i].substr(0, variant.alleles[i].length() - shared_prefix_len);
        }
        if (i == 0) {
            variant.ref = variant.alleles[i];
        } else {
            variant.alt[i - 1] = variant.alleles[i];
        }
    }
}

string VCFOutputCaller::print_snarl(const HandleGraph* graph, const handle_t& snarl_start,
                                    const handle_t& snarl_end, bool in_brackets) const {
    Snarl snarl;
    Visit* start = snarl.mutable_start();
    start->set_node_id(graph->get_id(snarl_start));
    start->set_backward(graph->get_is_reverse(snarl_start));
    Visit* end = snarl.mutable_end();
    end->set_node_id(graph->get_id(snarl_end));
    end->set_backward(graph->get_is_reverse(snarl_end));
    return this->print_snarl(snarl, in_brackets);
}

string VCFOutputCaller::print_snarl(const Snarl& snarl, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    nid_t start_node_id = snarl.start().node_id();
    nid_t end_node_id = snarl.end().node_id();
    string start_node = std::to_string(start_node_id);
    string end_node = std::to_string(end_node_id);
    if (translation) {
        auto i = translation->find(start_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + start_node + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + end_node + " in translation file");
        }
        end_node = i->second.first;
    }
    stringstream ss;
    if (in_brackets) {
        ss << "(";
    }
    ss << (snarl.start().backward() ? "<" : ">") << start_node << (snarl.end().backward() ? "<" : ">") << end_node;
    if (in_brackets) {
        ss << ")";
    }
    return ss.str();
}
string VCFOutputCaller::print_flipped_snarl(const Snarl& snarl, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    Snarl flipped_snarl;
    flipped_snarl.mutable_start()->set_node_id(snarl.end().node_id());
    flipped_snarl.mutable_start()->set_backward(!snarl.end().backward());
    flipped_snarl.mutable_end()->set_node_id(snarl.start().node_id());
    flipped_snarl.mutable_end()->set_backward(!snarl.start().backward());
    return print_snarl(flipped_snarl, in_brackets);
}

void VCFOutputCaller::scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const {
    int left = -1;
    int last = 0;
    Snarl snarl;
    string frag;
    for (int i = 0; i < allele_string.length(); ++i) {
        if (allele_string[i] == '(') {
            assert(left == -1);
            if (last < i) {
                frag = allele_string.substr(last, i-last);
                callback(frag, snarl);
            }
            left = i;
        } else if (allele_string[i] == ')') {
            assert(left >= 0 && i > left + 3);
            frag = allele_string.substr(left + 1, i - left - 1);
            auto toks = split_delims(frag, "><");
            assert(toks.size() == 2);
            assert(frag[0] == '<' || frag[0] == '>');
            int64_t start = std::stoi(toks[0]);
            snarl.mutable_start()->set_node_id(start);
            snarl.mutable_start()->set_backward(frag[0] == '<');
            assert(frag[toks[0].size() + 1] == '<' || frag[toks[0].size() + 1] == '>');
            int64_t end = std::stoi(toks[1]);
            snarl.mutable_end()->set_node_id(abs(end));
            snarl.mutable_end()->set_backward(frag[toks[0].size() + 1] == '<');
            callback("", snarl);
            left = -1;
            last = i + 1;
        }
    }
    if (last == 0) {
        callback(allele_string, snarl);
    } else {
        frag = allele_string.substr(last);
        callback(frag, snarl);
    }
}

GAFOutputCaller::GAFOutputCaller(AlignmentEmitter* emitter, const string& sample_name, const vector<string>& ref_paths,
                                 size_t trav_padding) :
    emitter(emitter),
    gaf_sample_name(sample_name),
    ref_paths(ref_paths.begin(), ref_paths.end()),
    trav_padding(trav_padding) {
    
}

GAFOutputCaller::~GAFOutputCaller() {
}

void GAFOutputCaller::emit_gaf_traversals(const PathHandleGraph& graph, const string& snarl_name,
                                          const vector<SnarlTraversal>& travs,
                                          int64_t ref_trav_idx,
                                          const string& ref_path_name, int64_t ref_path_position,
                                          const TraversalSupportFinder* support_finder) {
    assert(emitter != nullptr);
    vector<Alignment> aln_batch;
    aln_batch.reserve(travs.size());

    stringstream ss;
    if (!ref_path_name.empty()) {
        ss << ref_path_name << "#" << ref_path_position << "#";
    }
    ss << snarl_name << "#" << gaf_sample_name;
    string variant_id = ss.str();

    // create allele ordering where reference is 0
    vector<int> alleles;
    if (ref_trav_idx >= 0 && ref_trav_idx < (int64_t)travs.size()) {
        alleles.push_back(ref_trav_idx);
    }
    for (int i = 0; i < travs.size(); ++i) {
        if (i != ref_trav_idx) {
            alleles.push_back(i);
        }
    }
    // make an alignment for each traversal
    for (int i = 0; i < alleles.size(); ++i) {
        const SnarlTraversal& trav = travs[alleles[i]];
        Alignment trav_aln;
        if (trav_padding > 0) {
            trav_aln = to_alignment(pad_traversal(graph, trav), graph);
        } else {
            trav_aln = to_alignment(trav, graph);
        }
        trav_aln.set_name(variant_id + "#" + std::to_string(i));
        if (support_finder) {
            int64_t support = support_finder->support_val(support_finder->get_traversal_support(trav));
            set_annotation(trav_aln, "support", std::to_string(support));
        }        
        aln_batch.push_back(trav_aln);
    }
    emitter->emit_singles(std::move(aln_batch)); 
}

void GAFOutputCaller::emit_gaf_variant(const PathHandleGraph& graph, const string& snarl_name,
                                       const vector<SnarlTraversal>& travs,
                                       const vector<int>& genotype,
                                       int64_t ref_trav_idx,
                                       const string& ref_path_name, int64_t ref_path_position,
                                       const TraversalSupportFinder* support_finder) {
    assert(emitter != nullptr);

    // pretty bare bones for now, just output the genotype as a pair of traversals
    // todo: we could embed some basic information (likelihood, ploidy, sample etc) in the gaf
    // gt_travs is a NEW vector holding one entry per called allele, so ref_trav_idx -- an index into
    // travs -- does not address it.  Passing it through unremapped made emit_gaf_traversals index a
    // 2-element vector with an index into the full traversal list.  The genotype can also carry the
    // negative star/missing markers, which address nothing at all.
    vector<SnarlTraversal> gt_travs;
    int64_t gt_ref_trav_idx = -1;
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        if (allele == ref_trav_idx && gt_ref_trav_idx < 0) {
            gt_ref_trav_idx = (int64_t)gt_travs.size();
        }
        gt_travs.push_back(travs[allele]);
    }
    if (gt_travs.empty()) {
        return;
    }
    emit_gaf_traversals(graph, snarl_name, gt_travs, gt_ref_trav_idx, ref_path_name, ref_path_position, support_finder);
}

SnarlTraversal GAFOutputCaller::pad_traversal(const PathHandleGraph& graph, const SnarlTraversal& trav) const {

    assert(trav.visit_size() >= 2);

    SnarlTraversal out_trav;

    // traversal endpoints
    handle_t start_handle = graph.get_handle(trav.visit(0).node_id(), trav.visit(0).backward());
    handle_t end_handle = graph.get_handle(trav.visit(trav.visit_size() - 1).node_id(), trav.visit(trav.visit_size() - 1).backward());

    // find a reference path that touches the start node
    // todo: we could be more clever by finding the longest one or something
    path_handle_t reference_path;
    step_handle_t reference_step;
    bool found = false;
    size_t padding = 0;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            reference_path = graph.get_path_handle_of_step(step_handle);
            string name = graph.get_path_name(reference_path);
            if (!Paths::is_alt(name) && (ref_paths.empty() || ref_paths.count(name))) {
                reference_step = step_handle;
                found = true;
            }
            return !found;
        });

    // add left padding
    if (found) {
        deque<Visit> left_padding;

        if (graph.get_is_reverse(start_handle) == graph.get_is_reverse(graph.get_handle_of_step(reference_step))) {
            // path and handle oriented the same, we can just backtrack along the path to get previous stuff
            for (step_handle_t step = graph.get_previous_step(reference_step);
                 step != graph.path_front_end(reference_path) && padding < trav_padding;
                 step = graph.get_previous_step(step)) {
                left_padding.push_front(to_visit(graph, graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        } else {
            // path and handle oriented differently, we go forward in the path, flipping each step
            for (step_handle_t step = graph.get_next_step(reference_step);
                 step != graph.path_end(reference_path) && padding < trav_padding;
                 step = graph.get_next_step(step)) {
                left_padding.push_front(to_visit(graph, graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        }

        for (const Visit& visit : left_padding) {
            *out_trav.add_visit() = visit;
        }
    }

    // copy over center
    for (int i = 0; i < trav.visit_size(); ++i) {
        *out_trav.add_visit() = trav.visit(i);
    }

    // go through the whole thing again with the end
    found = false;
    padding = 0;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
            reference_path = graph.get_path_handle_of_step(step_handle);
            string name = graph.get_path_name(reference_path);
            if (!Paths::is_alt(name) && (ref_paths.empty() || ref_paths.count(name))) {
                reference_step = step_handle;
                found = true;
            }
            return !found;
        });

    // add right padding
    if (found) {
        if (graph.get_is_reverse(end_handle) == graph.get_is_reverse(graph.get_handle_of_step(reference_step))) {
            // path and handle oriented the same, we can just continue along the path to get next stuff
            for (step_handle_t step = graph.get_next_step(reference_step);
                 step != graph.path_end(reference_path) && padding < trav_padding;
                 step = graph.get_next_step(step)) {
                Visit* visit = out_trav.add_visit();
                *visit = to_visit(graph, graph.get_handle_of_step(step));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        } else {
            // path and handle oriented differently, we go backward in the path, flipping each step
            for (step_handle_t step = graph.get_previous_step(reference_step);
                 step != graph.path_front_end(reference_path) && padding < trav_padding;
                 step = graph.get_previous_step(step)) {
                Visit* visit = out_trav.add_visit();
                *visit = to_visit(graph, graph.flip(graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        }
    }
        
    return out_trav;
}

void VCFOutputCaller::update_nesting_info_tags(const SnarlManager* snarl_manager) {

    // index the snarl tree by name
    unordered_map<string, const Snarl*> name_to_snarl;
    Snarl flipped_snarl;
    snarl_manager->for_each_snarl_preorder([&](const Snarl* snarl) {
            name_to_snarl[print_snarl(*snarl)] = snarl;
            // also add a map from the flipped snarl (as call sometimes messes with orientation)
            flipped_snarl.mutable_start()->set_node_id(snarl->end().node_id());
            flipped_snarl.mutable_start()->set_backward(!snarl->end().backward());
            flipped_snarl.mutable_end()->set_node_id(snarl->start().node_id());
            flipped_snarl.mutable_end()->set_backward(!snarl->start().backward());
            name_to_snarl[print_snarl(flipped_snarl)] = snarl;
        });

    // pass 1) index sites in vcf
    // (todo: this could be done more quickly upstream)
    //
    // One index, not two: presence in chrom_of_name IS "this snarl name is in the VCF", and
    // the value is which reference contig its record landed on.  Keeping a separate
    // names_in_vcf set alongside would store all 400k snarl-ID strings twice, which measured
    // as +70 MB of peak RSS on chr22 -- the keys, not the values, are what costs.
    //
    // Contig names are interned rather than stored per record for the same reason: there are
    // at most a few thousand distinct ones, and all we ever ask is whether two are the same.
    unordered_map<string, uint32_t> chrom_index;
    auto intern_chrom = [&](const string& chrom) -> uint32_t {
        return chrom_index.emplace(chrom, (uint32_t)chrom_index.size()).first->second;
    };
    // One entry per snarl name.  A snarl ID is not unique -- a cyclic reference path that
    // traverses the same snarl twice emits two records with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa) -- but both occurrences are traversals of one
    // path, so they are on the same contig and it does not matter which one wins here.
    unordered_map<string, uint32_t> chrom_of_name;
    for (auto& thread_buf : output_variants) {
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 4);
            chrom_of_name.emplace(toks[2], intern_chrom(toks[0]));
        }
    }

    // pass 2) identify top-level snarls (those with no ancestors in VCF)
    // and store reference info only for them
    struct RefInfo {
        string chrom;
        size_t pos;
        size_t ref_len;
    };
    // Keyed by snarl name, then by (chrom, pos), because a snarl ID can carry more than one
    // record: a cyclic reference emits two, both with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa, which gives two <5<1 records at POS 20 and 44).
    // A plain name -> RefInfo map was last-write-wins, so both records were handed the
    // surviving one's interval and the record at POS 20 reported RS=44.  The inner map is
    // ordered so that picking begin() is deterministic regardless of thread scheduling.
    unordered_map<string, map<pair<string, size_t>, size_t>> top_level_ref_info;

    // Helper to check if a snarl is top-level (no ancestors in VCF)
    auto is_top_level = [&](const string& name) -> bool {
        auto it = name_to_snarl.find(name);
        if (it == name_to_snarl.end()) return true; // not found, treat as top-level
        const Snarl* snarl = it->second;
        while ((snarl = snarl_manager->parent_of(snarl))) {
            string cur_name = print_snarl(*snarl);
            string flipped_name = print_flipped_snarl(*snarl);
            if (chrom_of_name.count(cur_name) || chrom_of_name.count(flipped_name)) {
                return false; // has ancestor in VCF
            }
        }
        return true; // no ancestors in VCF
    };

    // Second pass through variants to extract ref info only for top-level snarls
    for (auto& thread_buf : output_variants) {
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            const string& name = toks[2];
            if (is_top_level(name)) {
                top_level_ref_info[name][make_pair(toks[0], static_cast<size_t>(stoul(toks[1])))] =
                    toks[3].length();
            }
        }
    }

    // determine the tags from the index
    //
    // There are exactly two ways a snarl can nest inside its parent's record, and they need
    // to be told apart.  A site inside a *deletion* is covered by its parent contig's own
    // reference allele, so it has coordinates on that contig and its record's CHROM is the
    // same.  A site inside an *insertion* has no path of the parent's contig through it at
    // all, so it is only callable once some other reference (a gref fragment) covers the
    // inserted allele -- and its record's CHROM is therefore different.  So:
    //
    //   contig_level    ancestors whose record is on this record's own CHROM, i.e. how deep
    //                   the site is in its own coordinate system
    //   contig_hops     steps in the chain where CHROM changed, i.e. how many insertions deep
    //                   the site is
    //
    // Returns: (contig_level, contig_hops, parent_name, top_level_name)
    function<tuple<size_t, size_t, string, string>(const string&, const string&)> get_nesting_tags =
        [&](const string& name, const string& my_chrom) {
        string parent_name;
        string top_level_name = name;  // default to self (for the top-level case)
        size_t contig_level = 0;
        size_t contig_hops = 0;
        // Our own contig, and the contig of the previously visited link in the chain.
        // chrom_index is complete after pass 1, so this lookup always hits.
        uint32_t my_chrom_id = chrom_index.at(my_chrom);
        uint32_t prev_chrom_id = my_chrom_id;
        const Snarl* snarl = name_to_snarl.at(name);

        assert(snarl != nullptr);
        // walk up the snarl tree
        while ((snarl = snarl_manager->parent_of(snarl))) {
            string cur_name = print_snarl(*snarl);

            // Since it is possible that the snarl is actually flipped in the vcf, check for the flipped version too
            string flipped_name = print_flipped_snarl(*snarl);
            const string* hit = nullptr;
            if (chrom_of_name.count(cur_name)) {
                // only count snarls that are in the vcf
                hit = &cur_name;
            } else if (chrom_of_name.count(flipped_name)) {
                // snarl is in vcf under flipped orientation
                hit = &flipped_name;
            }
            if (hit == nullptr) {
                continue;
            }

            auto chrom_it = chrom_of_name.find(*hit);
            uint32_t anc_chrom_id = chrom_it == chrom_of_name.end() ? my_chrom_id
                                                                    : chrom_it->second;
            if (anc_chrom_id == my_chrom_id) {
                ++contig_level;
            }
            if (anc_chrom_id != prev_chrom_id) {
                ++contig_hops;
            }
            prev_chrom_id = anc_chrom_id;

            if (parent_name.empty()) {
                // remember the first parent
                parent_name = *hit;
            }
            // keep updating top_level to find the topmost ancestor in VCF
            top_level_name = *hit;
        }
        return make_tuple(contig_level, contig_hops, parent_name, top_level_name);
    };

    // pass 3) add the LV, PS, RC, RS, RD tags
#pragma omp parallel for
    for (uint64_t i = 0; i < output_variants.size(); ++i) {
        auto& thread_buf = output_variants[i];
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            //string& output_variant_string = output_variant_record.second;
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            const string& name = toks[2];

            auto [contig_level, contig_hops, parent_name, top_level_name] =
                get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig.  It used to count
            // ancestors across every contig: for a VCF with a single reference contig the two
            // are identical, but once gref fragments give the insides of insertions their own
            // contigs, the whole-file count is not what a level filter wants.  No gref-contig
            // record was ever at LV=0 under the old definition, so `vcfbub -l 0` deleted every
            // one of them.
            string nesting_tags = ";LV=" + std::to_string(contig_level);
            nesting_tags += ";CH=" + std::to_string(contig_hops);
            if (!parent_name.empty()) {
                // Not "if (lv != 0)": those were equivalent only while LV was the absolute
                // count.  A record can now legitimately be at LV=0 and still have a parent on
                // another contig, and it must keep PS -- vcfbub's rescue of the children of
                // popped bubbles is keyed on it.
                nesting_tags += ";PS=" + parent_name;
            }

            // Add RC, RS, RD tags (reference info from top-level snarl)
            RefInfo top_ref;
            if (top_level_name == name) {
                // We are our own top-level site, so the answer is our own interval.  Looking
                // it up by name would risk picking a different record that shares our ID.
                top_ref = {toks[0], static_cast<size_t>(stoul(toks[1])), toks[3].length()};
            } else {
                const auto& candidates = top_level_ref_info.at(top_level_name);
                // If the ancestor produced several records, prefer one on our own contig;
                // failing that take the smallest (chrom, pos).  Which one is "right" is
                // genuinely ambiguous, so pick deterministically rather than by chance.
                auto chosen = candidates.begin();
                for (auto it = candidates.begin(); it != candidates.end(); ++it) {
                    if (it->first.first == toks[0]) {
                        chosen = it;
                        break;
                    }
                }
                top_ref = {chosen->first.first, chosen->first.second, chosen->second};
            }
            nesting_tags += ";RC=" + top_ref.chrom;
            nesting_tags += ";RS=" + std::to_string(top_ref.pos);
            nesting_tags += ";RD=" + std::to_string(top_ref.pos + top_ref.ref_len);

            // rewrite the output string using the updated info toks
            output_variant_string.clear();
            for (size_t i = 0; i < toks.size(); ++i) {
                output_variant_string += toks[i];
                if (i == 7) {
                    output_variant_string += nesting_tags;
                }
                if (i != toks.size() - 1) {
                    output_variant_string += "\t";
                }
            }
            output_variant_record.second.clear();
            ret = zstdutil::CompressString(output_variant_string, output_variant_record.second);
            assert(ret == 0);
        }
    }
}

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<int>& ref_path_ploidies,
                           FastaReference* ref_fasta,
                           FastaReference* ins_fasta,
                           AlignmentEmitter* aln_emitter,
                           bool traversals_only,
                           bool gaf_output,
                           size_t trav_padding) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    input_vcf(variant_file),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths, ref_fasta, ins_fasta, snarl_caller.get_skip_allele_fn()),
    traversals_only(traversals_only),
    gaf_output(gaf_output) {

    scan_contig_lengths();

    assert(ref_paths.size() == ref_path_ploidies.size());
    for (int i = 0; i < ref_paths.size(); ++i) {
        path_to_ploidy[ref_paths[i]] = ref_path_ploidies[i];
    }
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // could be that our graph is a subgraph of the graph the snarls were computed from
    // so bypass snarls we can't process
    if (!graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        return false;
    }

    // get our traversals out of the finder
    vector<pair<SnarlTraversal, vector<int>>> alleles;
    vector<vcflib::Variant*> variants;
    std::tie(alleles, variants) = traversal_finder.find_allele_traversals(snarl);

    if (!alleles.empty()) {

        // hmm, maybe find a way not to copy?
        vector<SnarlTraversal> travs;
        travs.reserve(alleles.size());
        for (const auto& ta : alleles) {
            travs.push_back(ta.first);
        }

        // find the reference traversal
        // todo: is it the reference always first?
        int ref_trav_idx = -1;
        for (int i = 0; i < alleles.size() && ref_trav_idx < 0; ++i) {
            if (std::all_of(alleles[i].second.begin(), alleles[i].second.end(), [](int x) {return x == 0;})) {
                ref_trav_idx = i;
            }
        }

        // find a path range corresponding to our snarl by way of the VCF variants.
        tuple<string, size_t, size_t> ref_positions = get_ref_positions(variants);

        // just print the traversals if requested
        if (traversals_only) {
            assert(gaf_output);
            // todo: can't get ref position here without pathposition graph
            emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, "", -1);
            return true;
        }
        
        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype;
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, path_to_ploidy[get<0>(ref_positions)],
                                                                        get<0>(ref_positions),  make_pair(get<1>(ref_positions), get<2>(ref_positions)));

        assert(trav_genotype.size() <= 2);

        if (gaf_output) {
            // todo: can't get ref position here without pathposition graph
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, "", -1);
            return true;
        }

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            vector<int> vcf_alleles;
            set<int> used_vcf_alleles;
            string vcf_genotype;
            vector<SnarlTraversal> vcf_traversals(variants[i]->alleles.size());            
            if (trav_genotype.empty()) {
                vcf_genotype = "./.";
            } else {
                // map our traversal genotype to a vcf variant genotype
                // using the information out of the traversal finder
                for (int j = 0; j < trav_genotype.size(); ++j) {
                    int trav_allele = trav_genotype[j];
                    int vcf_allele = alleles[trav_allele].second[i];
                    vcf_genotype += std::to_string(vcf_allele);
                    if (j < trav_genotype.size() - 1) {
                        vcf_genotype += "/";
                    }
                    if (!used_vcf_alleles.count(vcf_allele)) {                    
                        vcf_alleles.push_back(vcf_allele);
                        used_vcf_alleles.insert(vcf_allele);
                        vcf_traversals[vcf_allele] = travs[trav_allele];
                    }
                }
                // add traversals that correspond to vcf genotypes that are not
                // present in the traversal_genotypes
                for (int j = 0; j < travs.size(); ++j) {
                    int vcf_allele = alleles[j].second[i];
                    if (!used_vcf_alleles.count(vcf_allele)) {
                        vcf_traversals[vcf_allele] = travs[j];
                        used_vcf_alleles.insert(vcf_allele);
                    }
                }
            }
            // create an output variant from the input one
            vcflib::Variant out_variant;
            out_variant.sequenceName = variants[i]->sequenceName;
            out_variant.position = variants[i]->position;
            out_variant.id = variants[i]->id;
            out_variant.ref = variants[i]->ref;
            out_variant.alt = variants[i]->alt;
            out_variant.alleles = variants[i]->alleles;
            out_variant.filter = "PASS";
            out_variant.updateAlleleIndexes();

            // add the genotype
            out_variant.format.push_back("GT");
            auto& genotype_vector = out_variant.samples[sample_name]["GT"];
            genotype_vector.push_back(vcf_genotype);

            // add some info
            snarl_caller.update_vcf_info(snarl, vcf_traversals, vcf_alleles, trav_call_info, sample_name, out_variant);

            // print the variant
            add_variant(out_variant);
        }
        return true;
    }
    
    return false;

}

string VCFGenotyper::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    assert(contig_length_overrides.empty()); // using this override makes no sense

    // get the contig length overrides from the VCF
    vector<size_t> vcf_contig_lengths;
    auto length_map = scan_contig_lengths();
    for (int i = 0; i < ref_paths.size(); ++i) {
        vcf_contig_lengths.push_back(length_map[ref_paths[i]]);
    }
    
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, vcf_contig_lengths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

tuple<string, size_t, size_t> VCFGenotyper::get_ref_positions(const vector<vcflib::Variant*>& variants) const {
    // if there is more than one path in our snarl (unlikely for most graphs we'll vcf-genoetype)
    // then we return the one with the biggest interval
    map<string, pair<size_t, size_t>> path_offsets;
    for (const vcflib::Variant* var : variants) {
        if (path_offsets.count(var->sequenceName)) {
            pair<size_t, size_t>& record = path_offsets[var->sequenceName];
            record.first = std::min((size_t)var->position, record.first);
            record.second = std::max((size_t)var->position + var->ref.length(), record.second);
        } else {
            path_offsets[var->sequenceName] = make_pair(var->position, var->position + var->ref.length());
        }
    }

    string ref_path;
    size_t ref_range_size = 0;
    pair<size_t, size_t> ref_range;
    for (auto& path_offset : path_offsets) {
        size_t len = path_offset.second.second - path_offset.second.first;
        if (len > ref_range_size) {
            ref_range_size = len;
            ref_path = path_offset.first;
            ref_range = path_offset.second;
        }
    }

    return make_tuple(ref_path, ref_range.first, ref_range.second);
}

unordered_map<string, size_t> VCFGenotyper::scan_contig_lengths() const {

    unordered_map<string, size_t> ref_lengths;
    
    // copied from dumpContigsFromHeader.cpp in vcflib
    vector<string> headerLines = split(input_vcf.header, "\n");
    for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++) {
        if((*it).substr(0,8) == "##contig"){
            string contigInfo = (*it).substr(10, (*it).length() -11);
            vector<string> info = split(contigInfo, ",");
            string id;
            int64_t length = -1;
            for(vector<string>::iterator sub = info.begin(); sub != info.end(); sub++) {
                vector<string> subfield = split((*sub), "=");
                if(subfield[0] == "ID"){
                    id = subfield[1];
                }
                if(subfield[0] == "length"){
                    length = parse<int>(subfield[1]);
                }
            }
            if (!id.empty() && length >= 0) {
                ref_lengths[id] = length;
            }
        }
    }

    return ref_lengths;
}


LegacyCaller::LegacyCaller(const PathPositionHandleGraph& graph,
                           SupportBasedSnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets,
                           const vector<int>& ref_path_ploidies) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    graph(graph),
    ref_paths(ref_paths) {

    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    
    is_vg = dynamic_cast<const VG*>(&graph) != nullptr;
    if (is_vg) {
        // our graph is in vg format.  we index the paths and make a traversal finder just
        // like in the old call code
        for (auto ref_path : ref_paths) {
            path_indexes.push_back(new PathIndex(graph, ref_path));
        }
        // map snarl to the first reference path that spans it
        function<PathIndex*(const Snarl&)> get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, path_indexes).second;
        };
        // initialize our traversal finder
        traversal_finder = new RepresentativeTraversalFinder(graph, snarl_manager,
                                                             max_search_depth,
                                                             max_search_width,
                                                             max_bubble_paths,
                                                             0,
                                                             0,
                                                             get_path_index,
                                                             [&](id_t id) { return snarl_caller.get_support_finder().get_min_node_support(id);},
                                                             [&](edge_t edge) { return snarl_caller.get_support_finder().get_edge_support(edge);});

    } else {
        // our graph is not in vg format.  we will make graphs for each site as needed and work with those
        traversal_finder = nullptr;
    }
}

LegacyCaller::~LegacyCaller() {
    delete traversal_finder;
    for (PathIndex* path_index : path_indexes) {
        delete path_index;
    }
}

bool LegacyCaller::call_snarl(const Snarl& snarl) {

    // if we can't handle the snarl, then the GraphCaller framework will recurse on its children
    if (!is_traversable(snarl)) {
        return false;
    }
           
    RepresentativeTraversalFinder* rep_trav_finder;
    vector<PathIndex*> site_path_indexes;
    function<PathIndex*(const Snarl&)> get_path_index;
    VG vg_graph;
    SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
    bool was_called = false;
    
    if (is_vg) {
        // our graph is in VG format, so we've sorted this out in the constructor
        rep_trav_finder = traversal_finder;
        get_path_index = [&](const Snarl& site) {
            return find_index(site, path_indexes).second;
        };
        
    } else {
        // our graph isn't in VG format.  we are using a (hopefully temporary) workaround
        // of converting the subgraph into VG.
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&snarl, graph, true);
        size_t total_snarl_length = 0;
        for (auto node_id : contents.first) {
            handle_t new_handle = vg_graph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            if (node_id != snarl.start().node_id() && node_id != snarl.end().node_id()) {
                total_snarl_length += vg_graph.get_length(new_handle);
            }
        }
        for (auto edge : contents.second) {
            vg_graph.create_edge(vg_graph.get_handle(graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first)),
                                 vg_graph.get_handle(graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second)));
            total_snarl_length += 1;
        }
        // add the paths to the subgraph
        algorithms::expand_context_with_paths(&graph, &vg_graph, 1);
        // and index them
        for (auto& ref_path : ref_paths) {
            if (vg_graph.has_path(ref_path)) {
                site_path_indexes.push_back(new PathIndex(vg_graph, ref_path));
            } else {
                site_path_indexes.push_back(nullptr);
            }
        }
        get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, site_path_indexes).second;
        };
        // determine the support threshold for the traversal finder.  if we're using average
        // support, then we don't use any (set to 0), other wise, use the minimum support for a call
        SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
        size_t threshold = support_caller.get_support_finder().get_average_traversal_support_switch_threshold();
        double support_cutoff = total_snarl_length <= threshold ? support_caller.get_min_total_support_for_call() : 0;
        rep_trav_finder = new RepresentativeTraversalFinder(vg_graph, snarl_manager,
                                                            max_search_depth,
                                                            max_search_width,
                                                            max_bubble_paths,
                                                            support_cutoff,
                                                            support_cutoff,
                                                            get_path_index,
                                                            [&](id_t id) { return support_caller.get_support_finder().get_min_node_support(id);},
                                                            // note: because our traversal finder and support caller have
                                                            // different graphs, they can't share edge handles
                                                            [&](edge_t edge) { return support_caller.get_support_finder().get_edge_support(
                                                                    vg_graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first),
                                                                    vg_graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second));});
                                                            
    }

    PathIndex* path_index = get_path_index(snarl);
    if (path_index != nullptr) {
        string path_name = find_index(snarl, is_vg ? path_indexes : site_path_indexes).first;

        // orient the snarl along the reference path
        tuple<size_t, size_t, bool, step_handle_t, step_handle_t> ref_interval = get_ref_interval(graph, snarl, path_name);
        if (get<2>(ref_interval) == true) {
            snarl_manager.flip(&snarl);
        }

        // recursively genotype the site beginning here at the top level snarl
        vector<SnarlTraversal> called_traversals;
        // these integers map the called traversals to their positions in the list of all traversals
        // of the top level snarl.  
        vector<int> genotype;
        int ploidy = ploidy_at(path_name, get<0>(ref_interval),
                               ref_offsets.count(path_name) ? ref_offsets.at(path_name) : 0,
                               ref_ploidies[path_name]);
        std::tie(called_traversals, genotype) = top_down_genotype(snarl, *rep_trav_finder, ploidy,
                                                                  path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)));
    
        if (!called_traversals.empty()) {
            // regenotype our top-level traversals now that we know they aren't nested, and we have a
            // good idea of all the sizes
            unique_ptr<SnarlCaller::CallInfo> call_info;
            std::tie(called_traversals, genotype, call_info) = re_genotype(snarl, *rep_trav_finder, called_traversals, genotype, ploidy,
                                                                           path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)));

            // emit our vcf variant
            was_called = emit_variant(graph, snarl_caller, snarl, called_traversals, genotype, 0, call_info, path_name, ref_offsets.find(path_name)->second, false,
                         ploidy);

        }
    }        
    if (!is_vg) {
        // delete the temporary vg subgraph and traversal finder we created for this snarl
        delete rep_trav_finder;
        for (PathIndex* path_index : site_path_indexes) {
            delete path_index;
        }
    }

    return was_called;
}

string LegacyCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

pair<vector<SnarlTraversal>, vector<int>> LegacyCaller::top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy,
                                                                          const string& ref_path_name, pair<size_t, size_t> ref_interval) const {

    // get the traversals through the site
    vector<SnarlTraversal> traversals = trav_finder.find_traversals(snarl);

    // use our support caller to choose our genotype
    vector<int> trav_genotype;
    unique_ptr<SnarlCaller::CallInfo> trav_call_info;
    std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, traversals, 0, ploidy, ref_path_name, ref_interval);
    if (trav_genotype.empty()) {
        return make_pair(vector<SnarlTraversal>(), vector<int>());
    }

    assert(trav_genotype.size() == ploidy);

    vector<SnarlTraversal> called_travs(ploidy);

    // do we have two paths going through a given traversal?  This is handled
    // as a special case below
    bool hom = trav_genotype.size() == 2 && trav_genotype[0] == trav_genotype[1];
    
    for (int i = 0; i < trav_genotype.size() && (!hom || i < 1); ++i) {
        int allele = trav_genotype[i];
        const SnarlTraversal& traversal = traversals[allele];
        Visit prev_end;
        for (int j = 0; j < traversal.visit_size(); ++j) {
            if (traversal.visit(j).node_id() > 0) {
                *called_travs[i].add_visit() = traversal.visit(j);
                if (hom && i == 0) {
                    *called_travs[1].add_visit() = traversal.visit(j);
                }
            } else {
                // recursively determine the traversal
                const Snarl* into_snarl = snarl_manager.into_which_snarl(traversal.visit(j));
                bool flipped = traversal.visit(j).backward();
                if (flipped) {
                    // we're always processing our snarl from start to end, so make sure
                    // it lines up with the parent (note that we've oriented the root along the ref path)
                    snarl_manager.flip(into_snarl);
                }
                vector<SnarlTraversal> child_genotype = top_down_genotype(*into_snarl,
                                                                          trav_finder, hom ? 2: 1, ref_path_name, ref_interval).first;                
                if (child_genotype.empty()) {
                    return make_pair(vector<SnarlTraversal>(), vector<int>());
                }
                bool back_to_back = j > 0 && traversal.visit(j - 1).node_id() == 0 && prev_end == into_snarl->start();

                for (int k = back_to_back ? 1 : 0; k < child_genotype[0].visit_size(); ++k) {
                    *called_travs[i].add_visit() = child_genotype[0].visit(k);
                }
                if (hom) {
                    assert(child_genotype.size() == 2 && i == 0);
                    for (int k = back_to_back ? 1 : 0; k < child_genotype[1].visit_size(); ++k) {
                        *called_travs[1].add_visit() = child_genotype[1].visit(k);
                    }
                }
                prev_end = into_snarl->end();
                if (flipped) {
                    // leave our snarl like we found it
                    snarl_manager.flip(into_snarl);
                }
            }
        }
    }

    return make_pair(called_travs, trav_genotype);
}

SnarlTraversal LegacyCaller::get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const {

    // get the ref traversal through the site
    // todo: don't avoid so many traversal recomputations
    SnarlTraversal traversal = trav_finder.find_traversals(snarl)[0];
    SnarlTraversal out_traversal;

    Visit prev_end;
    for (int i = 0; i < traversal.visit_size(); ++i) {
        const Visit& visit = traversal.visit(i);
        if (visit.node_id() != 0) {
            *out_traversal.add_visit() = visit;
        } else {
            const Snarl* into_snarl = snarl_manager.into_which_snarl(visit);
            if (visit.backward()) {
                snarl_manager.flip(into_snarl);
            }
            bool back_to_back = i > 0 && traversal.visit(i - 1).node_id() == 0 && prev_end == into_snarl->start();

            SnarlTraversal child_ref = get_reference_traversal(*into_snarl, trav_finder);
            for (int j = back_to_back ? 1 : 0; j < child_ref.visit_size(); ++j) {
                *out_traversal.add_visit() = child_ref.visit(j);
            }
            prev_end = into_snarl->end();
            if (visit.backward()) {
                // leave our snarl like we found it
                snarl_manager.flip(into_snarl);
            }
        }
    }
    return out_traversal;    
}

tuple<vector<SnarlTraversal>, vector<int>, unique_ptr<SnarlCaller::CallInfo>>
LegacyCaller::re_genotype(const Snarl& snarl, TraversalFinder& trav_finder,
                          const vector<SnarlTraversal>& in_traversals,
                          const vector<int>& in_genotype,
                          int ploidy,
                          const string& ref_path_name,
                          pair<size_t, size_t> ref_interval) const {
    
    assert(in_traversals.size() == in_genotype.size());
    
    // create a set of unique traversal candidates that must include the reference first
    vector<SnarlTraversal> rg_traversals;
    // add our reference traversal to the front
    for (int i = 0; i < in_traversals.size() && !rg_traversals.empty(); ++i) {
        if (in_genotype[i] == 0) {
            rg_traversals.push_back(in_traversals[i]);
        }
    }
    if (rg_traversals.empty()) {
        rg_traversals.push_back(get_reference_traversal(snarl, trav_finder));
    }
    set<int> gt_set = {0};
    for (int i = 0; i < in_traversals.size(); ++i) {
        if (!gt_set.count(in_genotype[i])) {
            rg_traversals.push_back(in_traversals[i]);
            gt_set.insert(in_genotype[i]);
        }
    }
    
    // re-genotype the candidates
    vector<int> rg_genotype;
    unique_ptr<SnarlCaller::CallInfo> rg_call_info;
    std::tie(rg_genotype, rg_call_info) = snarl_caller.genotype(snarl, rg_traversals, 0, ploidy, ref_path_name, ref_interval);

    return make_tuple(rg_traversals, rg_genotype, std::move(rg_call_info));
}

bool LegacyCaller::is_traversable(const Snarl& snarl) {
    // we need this to be true all the way down to use the RepresentativeTraversalFinder on our snarl.
    bool ret = snarl.start_end_reachable() && snarl.directed_acyclic_net_graph() &&
       graph.has_node(snarl.start().node_id()) && graph.has_node(snarl.end().node_id());
    if (ret == true) {
        const vector<const Snarl*>& children = snarl_manager.children_of(&snarl);
        for (int i = 0; i < children.size() && ret; ++i) {
            ret = is_traversable(*children[i]);
        }
    }
    return ret;
}

pair<string, PathIndex*> LegacyCaller::find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const {
    assert(path_indexes.size() == ref_paths.size());
    for (int i = 0; i < path_indexes.size(); ++i) {
        PathIndex* path_index = path_indexes[i];
        if (path_index != nullptr &&
            path_index->by_id.count(snarl.start().node_id()) &&
            path_index->by_id.count(snarl.end().node_id())) {
            // This path threads through this site
            return make_pair(ref_paths[i], path_index);
        }
    }
    return make_pair("", nullptr);
}

FlowCaller::FlowCaller(const PathPositionHandleGraph& graph,
                       SupportBasedSnarlCaller& snarl_caller,
                       SnarlManager& snarl_manager,
                       const string& sample_name,
                       TraversalFinder& traversal_finder,
                       const vector<string>& ref_paths,
                       const vector<size_t>& ref_path_offsets,
                       const vector<int>& ref_path_ploidies,
                       AlignmentEmitter* aln_emitter,
                       bool traversals_only,
                       bool gaf_output,
                       size_t trav_padding,
                       bool genotype_snarls,
                       const pair<size_t, size_t>& allele_length_range) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }

}
   
FlowCaller::FlowCaller(const PathPositionHandleGraph& graph,
                       SupportBasedSnarlCaller& snarl_caller,
                       SnarlManager& snarl_manager,
                       const string& sample_name,
                       TraversalFinder& traversal_finder,
                       const vector<string>& ref_paths,
                       const vector<size_t>& ref_path_offsets,
                       const vector<int>& ref_path_ploidies,
                       AlignmentEmitter* aln_emitter,
                       bool traversals_only,
                       bool gaf_output,
                       size_t trav_padding,
                       bool genotype_snarls,
                       const pair<size_t, size_t>& allele_length_range,
                       bool nested,
                       bool star_allele) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    allele_length_range(allele_length_range),
    nested(nested),
    star_allele(star_allele)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
}

FlowCaller::~FlowCaller() {

}

bool FlowCaller::call_snarl(const Snarl& managed_snarl) {
    // Entry point: call with no parent context
    return call_snarl_internal(managed_snarl, "", make_pair(0, 0), nullptr);
}

TraversalSet FlowCaller::find_child_traversal_set(const SnarlTraversal& parent_trav,
                                                   const Snarl& child) const {
    TraversalSet result;

    // First, check if the parent traversal goes through this child snarl
    // by finding the child's start and end nodes in the parent
    nid_t child_start_id = child.start().node_id();
    nid_t child_end_id = child.end().node_id();
    bool found_start = false, found_end = false;

    for (int i = 0; i < parent_trav.visit_size(); ++i) {
        nid_t visit_id = parent_trav.visit(i).node_id();
        if (visit_id == child_start_id) found_start = true;
        if (visit_id == child_end_id) found_end = true;
    }

    // If parent doesn't traverse the child, return empty set (star allele case)
    if (!found_start || !found_end) {
        return result;
    }

    // Use the traversal finder to enumerate all traversals through the child
    FlowTraversalFinder* flow_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_finder != nullptr) {
        auto weighted_travs = flow_finder->find_weighted_traversals(child, false);
        result = std::move(weighted_travs.first);
    } else {
        result = traversal_finder.find_traversals(child);
    }

    return result;
}

int FlowCaller::crossings_of_child(const SnarlTraversal& trav, const Snarl& child) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    // Count crossings: an entry at one boundary followed by the other. Order matters -- testing
    // for the two boundaries independently would count a traversal that touches both on
    // unrelated excursions, which is the bug in find_child_traversal_set.
    int crossings = 0;
    nid_t open = 0;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        nid_t node = trav.visit(i).node_id();
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
        } else if (open != 0 && node == open) {
            ++crossings;
            open = 0;
        }
    }
    return crossings;
}

/// The bp offset of a child's start along a parent traversal, and the child's bp span along it.
///
/// Stage 14 of planning/decide-then-render.md: the linkage model measures the distance between
/// adjacent sites along the REFERENCE, but a haplotype that carries an insertion has travelled
/// further between the same two sites than the reference has, and a haplotype carrying a deletion
/// has travelled less. This is the distance in the haplotype's own frame.
///
/// Returns false where the child's boundaries are not both crossed by this traversal, which is not
/// an error: a traversal that does not enter the child has no offset along it to report.
static bool traversal_offset_span(const PathHandleGraph& graph, const SnarlTraversal& trav,
                                  const Snarl& child, int64_t* offset, int64_t* span,
                                  int64_t* total = nullptr, bool* entered_at_start = nullptr) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    int64_t walked = 0;
    int64_t entered_at = -1;
    bool at_start = true;
    nid_t open = 0;
    bool found = false;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        const nid_t node = trav.visit(i).node_id();
        if (!found) {
            if (open == 0 && (node == start || node == end)) {
                open = (node == start) ? end : start;
                entered_at = walked;
                // WHICH boundary was entered. Without this a grandchild under a crossing the parent
                // makes end-to-start is measured from the wrong end, so its whole subtree is
                // mirrored and its sibling order reverses. Masked today only because v1 descends
                // solely where the reference also goes, so every child is flipped onto its reference
                // path -- and stage 16 removes exactly that gate.
                at_start = (node == start);
            } else if (open != 0 && node == open) {
                // Span measured boundary-to-boundary inclusive of the closing node, matching how the
                // reference span of a snarl is taken. NOTE this is the child's extent along the
                // PARENT, not the child's own traversal length.
                *offset = entered_at;
                *span = walked + (int64_t)graph.get_length(graph.get_handle(node)) - entered_at;
                if (entered_at_start != nullptr) {
                    *entered_at_start = at_start;
                }
                found = true;
                // Visits after the first are masked, by decision: one copy for ploidy and the first
                // crossing for distance. Representing a second copy is stage 17's question. So the
                // walk continues only to total the traversal, not to look for another crossing.
                if (total == nullptr) {
                    return true;
                }
            }
        }
        walked += (int64_t)graph.get_length(graph.get_handle(node));
    }
    if (total != nullptr) {
        // The traversal's whole length, so the distance from this child to one under a LATER parent
        // can be formed as (tail of this traversal) + (anchor gap) + (offset in the next), which is a
        // pairwise quantity and therefore cannot reorder anything.
        *total = walked;
    }
    return found;
}

/// Stage 14 instrumentation, inert: does the haplotype frame reorder or re-space adjacent sites?
///
/// Aggregated rather than stored. The stage's off-ramp is a decision about whether stage 15 should
/// carry per-haplotype distances at all, and that decision needs distributions, not per-site fields
/// -- so nothing is added to `Entry` or to `record()` until the answer says it is worth it.
static std::atomic<size_t> g_frame_pairs{0};          // adjacent sibling pairs measured
static std::atomic<size_t> g_frame_reordered{0};      // ... that swap order in the haplotype frame
static std::atomic<size_t> g_frame_gap_within_5pct{0};
static std::atomic<size_t> g_frame_gap_measured{0};
static std::atomic<size_t> g_frame_no_offset{0};      // children with no offset on a called traversal
static std::atomic<size_t> g_frame_disagree{0};       // diploid parents whose traversals disagree
static std::atomic<size_t> g_frame_offset_delta_sum{0};
static std::atomic<size_t> g_frame_offset_delta_max{0};

uint64_t FlowCaller::child_crossing_mask(const vector<SnarlTraversal>& travs,
                                         const Snarl& child, bool* known) {
    if (known != nullptr) {
        *known = true;
    }
    // One bit per *candidate traversal*, not per VCF allele. The linkage layer settles a site on a
    // traversal pair, so that is the space a crossing question has to be asked in; asking it in
    // emitted-allele space meant the answer had to be mapped, and the mapping is what two of this
    // caller's worst bugs lived in. There is nothing to map now -- bit i is "travs[i] crosses child".
    if (travs.size() > 64) {
        // Unknown rather than none: the mask cannot index this site's candidates.
        if (known != nullptr) {
            *known = false;
        }
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < travs.size(); ++i) {
        if (crossings_of_child(travs[i], child) > 0) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

int FlowCaller::child_ploidy(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                             const Snarl& child, int cap) const {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;   // star or missing: that haplotype contributes no copy here
        }
        int crossings = crossings_of_child(travs[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;   // a cycle or tandem duplication; see the header comment
        }
        copies += crossings;
    }
    if (capped) {
        // Counted, not printed per occurrence. Masking visits after the first is a decision, not an
        // accident (see planning/decide-then-render.md 15'), and the size of what it masks is the
        // size of the deferred copy-number question -- so it needs a number reported once a run, not
        // a line per site gated on --progress that has to be grepped out of 24 logs. Measured that
        // way: 0 on chr20, 242 on chrX.
        ++g_child_multi_crossing;
    }
    return min(copies, cap);
}

void FlowCaller::set_defer_nested_descent(bool defer) {
    this->defer_nested_descent = defer;
    if (defer) {
        // Sized once, here, rather than lazily inside the parallel region that writes it.
        size_t threads = max((size_t)get_thread_count(), (size_t)omp_get_max_threads());
        // resize, not assign: PendingRecord owns a unique_ptr and so is move-only, and
        // assign(n, {}) would need to copy its prototype into each slot.
        pending_records.clear();
        pending_records.resize(max(threads, (size_t)1));
        render_records.clear();
        render_records.resize(max(threads, (size_t)1));
    }
}

size_t FlowCaller::pending_record_count() const {
    size_t n = 0;
    for (const auto& queue : pending_records) {
        n += queue.size();
    }
    return n;
}

size_t FlowCaller::render_record_count() const {
    size_t n = 0;
    for (const auto& queue : render_records) {
        n += queue.size();
    }
    return n;
}

/// Stage the inputs a record could be rendered from, for a snarl the barrier will not revise.
///
/// Everything `emit_variant` needs that is not recoverable from the graph: the traversals, the
/// genotype over them, and the CallInfo -- which is not optional. `update_vcf_info` maps emitted
/// alleles back to matrix columns by structural comparison of SnarlTraversal objects, indexes GL by
/// the sorted matrix-column multiset, and derives QUAL by renormalising the all-reference posterior
/// over that same map. Retaining allele strings alone cannot rebuild any of it.
unique_ptr<FlowCaller::PendingRecord> FlowCaller::stage_render_record(
        const Snarl& snarl, const vector<int>& trav_genotype, int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info,
        const string& ref_path_name, int ref_offset, int ploidy, bool emitted) {
    if (render_records.empty()) {
        return nullptr;
    }
    unique_ptr<PendingRecord> rec(new PendingRecord());
    rec->snarl = snarl;
    rec->ref_path_name = ref_path_name;
    rec->ref_offset = ref_offset;
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->record_key = record_key_of(snarl);
    rec->generation = 0;
    rec->emitted = emitted;
    rec->call_info = std::move(call_info);
    // `travs` is NOT taken here. Descent runs after every emit branch, top-level included, and reads
    // `travs` to work out which children the called alleles reach -- so moving it out at emit time
    // empties it before that loop and every child comes back with a copy number of zero. That is the
    // failure the nested branch's staging discipline exists to avoid (2,494 chr20 records, commit
    // 906812957); doing it in the top-level branch cost 12,302. Staged here, completed after descent.
    return rec;
}

pair<string, size_t> FlowCaller::site_ref_key(const Snarl& snarl, const string& ref_path_name,
                                             int ref_offset) const {
    // Pre-flatten, and locus-spelled. The flattened POS depends on which alleles the line carries,
    // so it does not exist before the record does; what the model uses position for is ordering and
    // the transition gaps, which a pre-flatten position serves exactly as well. The locus reduction
    // matters because `get_ref_position` answers with the base path name, "CHM13#0#chr20", where the
    // rest of the layer spells it "chr20".
    pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset);
    const string locus = PathMetadata::parse_locus_name(pos_info.first);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        pos_info.first = locus;
    }
    return make_pair(pos_info.first, (size_t)max((int64_t)0, pos_info.second));
}

void FlowCaller::record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                            const vector<int>& trav_genotype,
                            const unique_ptr<SnarlCaller::CallInfo>& call_info,
                            const string& ref_path_name, int ref_offset) {
    if (linkage_collector == nullptr || suppress_linkage_record) {
        return;
    }
    const auto* rl_info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
    if (rl_info == nullptr) {
        return;
    }
    // The same admission test the emitter used, in traversal space rather than emitted-allele space:
    // a genotype of one or two alleles, none of them a missing or star marker. Haploid chains are
    // included -- dropping them once cost chrY and non-pseudoautosomal chrX the linkage layer and
    // the mosaic entirely, about 5% of a genome and the part where the mosaic is the whole answer.
    const size_t site_ploidy = trav_genotype.size();
    if (site_ploidy != 1 && site_ploidy != 2) {
        return;
    }
    for (int allele : trav_genotype) {
        if (allele < 0) {
            return;
        }
    }
    // Position is pre-flatten now, and that is fine only because the patch index no longer keys on
    // it (6d8fef2c3). What the model uses position for is ordering and the transition gaps.
    pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset);
    // The contig must be spelled the way the VCF spells it, because the patch index keys on it.
    // `get_ref_position` returns the base path name -- "CHM13#0#chr20" -- while `emit_variant`
    // additionally reduces it to the locus, "chr20". Recording the unreduced form makes every
    // lookup miss: chr20 came out with every record unphased and no PS at all, which is what this
    // motion's byte-identity gate is for.
    {
        const string locus = PathMetadata::parse_locus_name(pos_info.first);
        if (locus != PathMetadata::NO_LOCUS_NAME) {
            pos_info.first = locus;
        }
    }
    const int called_i = trav_genotype[0];
    const int called_j = site_ploidy > 1 ? trav_genotype[1] : called_i;
    // No allele map: the emitted allele list does not exist yet, and it is chosen while the record is
    // built. `set_allele_map` supplies it afterwards, and stage 11 removes the need for it.
    static const vector<int> no_allele_map;
    linkage_collector->record(
        pos_info.first, (size_t)max((int64_t)0, pos_info.second),
        rl_info->genotype_lls,
        panel_alleles(graph, travs),
        called_i, called_j, no_allele_map,
        record_key_of(snarl),
        rl_info->explained_share, site_ploidy,
        (int64_t)snarl.start().node_id(), (int64_t)snarl.end().node_id(),
        nested_context.active, nested_context.parent_record_key,
        nested_context.parent_trav, nested_context.parent_crossing,
        current_generation, /*emitted*/ false);
}

void FlowCaller::render_retained_records() {
    // The phase, before any record is built: every generation has settled by now, so the phasing is
    // complete, and each record is phased as it is rendered rather than patched afterwards.
    build_render_phases();
    if (render_records.empty()) {
        return;
    }
    // Thread-locals are the hazard here, not the records. `emit_variant` reads `nested_context` and
    // `current_generation` when it records the site, and in a batch pass those hold whatever the last
    // snarl this thread happened to genotype left behind. A stale `nested_context.active` would file a
    // top-level site as nested -- which does not fail loudly: it would put the site into the nested
    // strand population instead of the diploid chain, silently. Every record here is one the barrier
    // does not revise, so the context is reset per record rather than trusted.
    const size_t n_threads = render_records.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < n_threads; ++t) {
        NestedContext saved_ctx = nested_context;
        size_t saved_gen = current_generation;
        nested_context = NestedContext();
        current_generation = 0;
        for (PendingRecord& rec : render_records[t]) {
            // The settled pair, not the one the reads alone picked. This is the whole point of the
            // phase: the ALT list, the symbolic-reference test that decides whether a line exists at
            // all, QUAL, and the arity of AD/GL/GQI are all built by iterating the genotype handed in
            // -- so handing in the settled one makes every one of them agree with the call instead of
            // being patched towards it afterwards. A settled traversal is renderable by construction,
            // because the allele list is chosen from it.
            vector<int> genotype = rec.genotype;
            int settled_a = -1, settled_b = -1;
            size_t settled_ploidy = 0;
            if (linkage_collector != nullptr
                && linkage_collector->settled_traversals(rec.record_key, &settled_a, &settled_b,
                                                         &settled_ploidy)
                && settled_ploidy == genotype.size()) {
                genotype.assign(1, settled_a);
                if (settled_ploidy > 1) {
                    genotype.push_back(settled_b);
                }
            }
            emit_variant(graph, snarl_caller, rec.snarl, rec.travs, genotype, rec.ref_trav_idx,
                         rec.call_info, rec.ref_path_name, rec.ref_offset, genotype_snarls,
                         rec.ploidy);
        }
        nested_context = saved_ctx;
        current_generation = saved_gen;
    }
    if (show_progress) {
        cerr << "[vg call] rendered " << render_record_count()
             << " retained records after the sweep" << endl;
    }
}

void FlowCaller::run_deferred_descent() {
    if (!defer_nested_descent) {
        return;
    }
    // Descent already happened, inline, during the one sweep the reads were resident for. What is
    // left is to settle the chains in the order their ploidies depend on: a generation's parents
    // before its children. Nothing here touches the reads.
    size_t generations = 0;
    if (linkage_collector != nullptr) {
        generations = linkage_collector->max_generation();
    }
    // Merged once: the sweep filled these per thread, and the barrier walks them in generation
    // order.
    vector<PendingRecord> pending;
    pending.reserve(pending_record_count());
    for (auto& queue : pending_records) {
        std::move(queue.begin(), queue.end(), std::back_inserter(pending));
        queue.clear();
    }

    // Blank the exact line the sweep buffered for a chain, so a retraction leaves nothing and a
    // replacement leaves exactly one copy. The handle was captured at emit time; identifying the
    // line later by re-hashing its ID column and counting GT separators mistook a phased haploid
    // replacement ("1|.") for a diploid line, and its gate could be closed while replacements
    // existed -- both lines then reached the output.
    auto blank_buffered_line = [&](PendingRecord& pr) {
        if (pr.buffer_thread >= 0 && (size_t)pr.buffer_thread < output_variants.size()
            && pr.buffer_index < output_variants[pr.buffer_thread].size()) {
            output_variants[pr.buffer_thread][pr.buffer_index].second.clear();
        }
        pr.buffer_thread = -1;
    };

    // parent record key -> indices of its pending children, so that dropping a chain can drop
    // everything under it. Built once: `pending` does not grow during the barrier.
    unordered_map<size_t, vector<size_t>> children_of;
    children_of.reserve(pending.size() * 2);
    for (size_t i = 0; i < pending.size(); ++i) {
        children_of[pending[i].parent_record_key].push_back(i);
    }

    // record key -> the record itself, over BOTH containers. `children_of` is parent-to-children and
    // cannot answer the question a frame needs, which is "give me my parent's traversals". And the
    // parents of generation 1 -- the largest slice of nested sites by far -- are top-level records in
    // `render_records`, which the barrier otherwise never indexes at all, so a `pending`-only map
    // would leave exactly that slice with nothing to measure along.
    // Frames not written, by reason. Counted rather than defaulted: an unset frame must be visible,
    // because a site that silently keeps a default sorts to the head of its group and hands its
    // neighbour a nonsense gap.
    size_t frame_written = 0, frame_no_entry = 0, frame_not_crossed = 0;
    size_t frame_no_parent = 0, frame_no_single_trav = 0;
    unordered_map<size_t, PendingRecord*> record_by_key;
    record_by_key.reserve((pending.size() + render_record_count()) * 2);
    for (PendingRecord& pr : pending) {
        record_by_key[pr.record_key] = &pr;
    }
    for (auto& queue : render_records) {
        for (PendingRecord& pr : queue) {
            record_by_key[pr.record_key] = &pr;
        }
    }
    // Drop a chain and its whole subtree: the settled parent does not carry the chain, so the
    // sample has no copy of it, and nothing nested inside a sequence the sample lacks exists
    // either. Returns how many entries were actually retracted, for the report.
    //
    // Iterative rather than recursive because the depth is data, not a constant, and breadth-first
    // over an explicit stack cannot blow the C++ stack on a pathological hierarchy.
    std::function<size_t(size_t)> drop_subtree = [&](size_t root) -> size_t {
        size_t dropped_here = 0;
        vector<size_t> stack{root};
        while (!stack.empty()) {
            size_t idx = stack.back();
            stack.pop_back();
            PendingRecord& victim = pending[idx];
            if (victim.dropped) {
                continue;
            }
            victim.dropped = true;
            if (linkage_collector != nullptr && linkage_collector->retract(victim.record_key)) {
                ++dropped_here;
            }
            if (victim.emitted) {
                blank_buffered_line(victim);
                victim.emitted = false;
            }
            auto kids = children_of.find(victim.record_key);
            if (kids != children_of.end()) {
                for (size_t k : kids->second) {
                    if (k != idx) {
                        stack.push_back(k);
                    }
                }
            }
        }
        return dropped_here;
    };

    size_t revised = 0, retracted = 0, gained = 0, crossing_unknown = 0, stale_respecify = 0;
    // `generations` is re-read at the end of each pass rather than snapshotted once: gaining a
    // chain records a linkage entry at a generation the collector may never have held before, and
    // a fixed bound would leave that entry -- and every pending chain below it -- outside every
    // resolve pass: emitted but never settled, never phased, absent from the mosaic.
    for (size_t gen = 0; gen <= generations; ++gen) {
        // The final pass carries last=true: it builds the phasing map and the mosaic from the full
        // accumulated set. If this pass then gains a deeper chain, the bound grows and a later
        // iteration re-runs that bookkeeping over the fuller set; the phasing map is rebuilt from
        // scratch there, so nothing is double-counted. The old shape -- a loop of last=false passes
        // and one extra last=true call -- resolved the final generation twice and appended every
        // final-generation site's PhaseCall to the mosaic input a second time.
        resolve_linkage_generation(gen, gen == generations);

        // This generation's parents are settled, so every chain hanging off one can be put on the
        // ploidy its parent's *final* genotype implies, before its own generation resolves. That is
        // the whole coherence guarantee, and it is a revision rather than a re-call because the
        // genotyping kept both ploidies' answers when the reads were resident.
        unordered_map<size_t, const LinkageCollector::PhaseCall*> settled;
        settled.reserve(linkage_phased.size() * 2);
        for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
            settled[pc.record_key] = &pc;
        }
        for (size_t i = 0; i < pending.size(); ++i) {
            PendingRecord& pr = pending[i];
            if (pr.generation != gen + 1 || pr.dropped) {
                continue;
            }
            if (!pr.crossing_known) {
                // The sweep could not compute this chain's crossing mask: its parent emitted
                // nothing (a retained chain -- re-emitting the parent below recomputes these), or
                // the parent carries more alleles than a 64-bit mask can index. Left exactly as
                // the sweep left it, and counted, because reading an unknown mask as "no allele
                // crosses" silently exempted these chains from revision.
                ++crossing_unknown;
                continue;
            }
            if (pr.parent_crossing == 0) {
                continue;   // no emitted parent allele crosses: no settled genotype can reach it
            }
            auto found = settled.find(pr.parent_record_key);
            if (found == settled.end()) {
                // The parent emitted no record, so linkage never touched its genotype: the ploidy
                // this chain was called at came from a genotype that was already final.
                continue;
            }
            const LinkageCollector::PhaseCall& parent = *found->second;
            // The settled pair as traversals, because that is what the mask indexes. Testing the
            // compact allele index here retracted 3,615 chains against a true 190 on chr20: the two
            // agree only when every allele at the parent is panel-carried.
            bool first = parent.trav_first >= 0 && parent.trav_first < 64
                         && ((pr.parent_crossing >> parent.trav_first) & 1);
            bool second = parent.ploidy == 2 && parent.trav_second >= 0 && parent.trav_second < 64
                          && ((pr.parent_crossing >> parent.trav_second) & 1);
            int copies = (int)first + (int)second;

            // The one derivation. Which of the parent's settled traversals carries this chain is the
            // same fact as how many copies of it the sample has, so both come from `first`/`second`
            // and nothing downstream re-decides either. Recorded even when the ploidy needs no
            // change: the descent-time value was computed against the parent's *pre-linkage*
            // genotype, and leaving it stale is what let the phasing pass and the barrier disagree
            // about a chain neither had any reason to doubt.
            // -1 means no settled traversal carries the chain; -2 means both do. Both are recorded
            // here rather than recomputed when the child is phased, which is the point: the count
            // and the identity come from the same two booleans, so they cannot disagree.
            pr.parent_trav = copies == 1 ? (first ? parent.trav_first : parent.trav_second)
                                         : (copies == 2 ? -2 : -1);
            linkage_collector->set_parent_trav(pr.record_key, pr.parent_trav);

            // Stage 15': where this chain sits along the traversals its parent SETTLED on,
            // measured now because that is when the parent's genotype is known. One slot per settled
            // traversal, indexed by the parent's traversal ORDER rather than by strand -- a haploid
            // parent has trav_first == trav_second, so both slots are written and neither can be read
            // unset, which is the failure a strand-keyed pair invites.
            //
            // Not from the CALLED traversals at descent: `carrying_trav` is set only under
            // `copies == 1`, so the chains a called allele never reached -- the ones the barrier may
            // now have moved the parent onto -- would have had nothing to measure along there.
            {
                auto parent_rec = record_by_key.find(pr.parent_record_key);
                if (parent_rec == record_by_key.end()) {
                    ++frame_no_parent;
                } else {
                    const vector<SnarlTraversal>& ptravs = parent_rec->second->travs;
                    const int settled[2] = {parent.trav_first, parent.trav_second};
                    bool any = false;
                    for (int slot = 0; slot < 2; ++slot) {
                        const int t = settled[slot];
                        if (t < 0 || (size_t)t >= ptravs.size()) {
                            continue;
                        }
                        int64_t off = 0, span = 0, total = 0;
                        bool at_start = true;
                        if (!traversal_offset_span(graph, ptravs[t], pr.snarl, &off, &span, &total,
                                                   &at_start)) {
                            continue;   // this settled traversal does not enter the chain
                        }
                        if (linkage_collector->set_frame(pr.record_key, slot, (int)off,
                                                        (int)(off + span), (int)total, !at_start)) {
                            any = true;
                        } else {
                            ++frame_no_entry;
                        }
                    }
                    if (any) {
                        ++frame_written;
                    } else {
                        ++frame_not_crossed;
                    }
                }
            }

            if (copies == 0) {
                // A call on a haplotype the sample turns out not to have -- and everything inside
                // it is on that same absent haplotype, so the whole subtree goes with it.
                //
                // Two things here were wrong. The retraction was conditional on a line existing,
                // which since collapsed sites started being recorded left line-less entries in the
                // layer at a ploidy their parent contradicts. And it never reached descendants: a
                // grandchild kept its own line and pointed at an entry that no longer existed,
                // which is precisely what the phasing pass reports as "no phased parent" -- and it
                // explains why that count was zero at generation 1, where the parent is top-level
                // and so was never a pending record that could be retracted.
                retracted += drop_subtree(i);
                continue;
            }
            if (copies == pr.ploidy && linkage_collector != nullptr
                && linkage_collector->has_entry(pr.record_key)) {
                // The ploidy the chain was called at is the one its parent's settled genotype
                // implies, so its entry and its staged record are both already right. Nothing to
                // revise; the render writes it from the settled pair like every other record.
                //
                // `has_entry` is not redundant. A chain no called parent allele reached at sweep
                // time is staged but deliberately NOT recorded (2,713 of them on chr20), so for
                // those the ploidy matching proves nothing: falling through is what files them in
                // the layer, via the record() fallback below. Skipping on ploidy alone would render
                // a nested record that had never entered the linkage layer at all.
                continue;
            }

            // emit_variant indexes its traversal vector directly -- `called_traversals[ref_trav_idx]`
            // with no bounds check -- and a chain held back during the sweep has never been through
            // it, so a missing reference traversal or an empty candidate list reaches it here for the
            // first time. That segfaulted. Anything unrenderable is left alone rather than guessed at.
            if (pr.travs.empty() || pr.ref_trav_idx < 0
                || (size_t)pr.ref_trav_idx >= pr.travs.size()) {
                continue;
            }
            bool genotype_in_range = !pr.genotype.empty();
            for (int allele : pr.genotype) {
                if (allele >= 0 && (size_t)allele >= pr.travs.size()) {
                    genotype_in_range = false;
                }
            }
            if (!genotype_in_range) {
                continue;
            }

            // Build the record at the ploidy the settled parent implies, from the genotyping kept
            // when the reads were resident. `alt_ploidy_info` holds the other ploidy's whole answer,
            // so this needs no re-reading and no re-scoring.
            ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* rl =
                dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(pr.call_info.get());
            unique_ptr<SnarlCaller::CallInfo> use_info;
            vector<int> use_genotype;
            if (copies == pr.ploidy) {
                use_genotype = pr.genotype;
            } else if (rl != nullptr && rl->alt_ploidy_info != nullptr
                       && (int)rl->alt_ploidy_info->ploidy == copies
                       && (int)rl->alt_ploidy_best.size() == copies) {
                bool ok = true;
                for (int allele : rl->alt_ploidy_best) {
                    if (allele < 0 || (size_t)allele >= pr.travs.size()) {
                        ok = false;
                    }
                }
                if (!ok) {
                    continue;
                }
                use_info.reset(rl->alt_ploidy_info.release());
                use_genotype = rl->alt_ploidy_best;
            } else {
                // No answer at the ploidy wanted -- too few traversals for a second genotype, or the
                // alternate was never computed. Left as it stands rather than invented.
                continue;
            }
            // "Gained" now means "was not in the linkage layer", not "had no line": no nested
            // chain has a line at this point, so the old test was true for every record and reported
            // 0 revised against 2,950 gained. A chain no called parent allele reached at sweep time
            // is the one that is genuinely new to the layer here.
            const bool was_gained = linkage_collector == nullptr
                                    || !linkage_collector->has_entry(pr.record_key);
            // Revised, not re-emitted. The chain's ploidy has changed, so the record it will be
            // rendered from has to change with it -- but the record does not exist yet, and that is
            // the point: there is no line to blank, no replacement to register, and no patch to
            // apply. The staged inputs are updated and the render pass builds the line once, at the
            // end, from whatever the layer finally settles on.
            pr.genotype = use_genotype;
            pr.ploidy = copies;
            if (use_info != nullptr) {
                pr.call_info = std::move(use_info);
            }
            const unique_ptr<SnarlCaller::CallInfo>& info = pr.call_info;

            const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* used =
                dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(info.get());
            if (used != nullptr) {
                // Traversal space here as at the inline site, so the barrier and the sweep describe
                // a site identically. The remap this block used to do -- GLs, panel and genotype all
                // pushed into the emitted numbering -- is gone: the collector builds its own compact
                // space, and pr.travs is the candidate list those indices already refer to.
                // No allele map, for the same reason `record_site` passes none: the emitted allele
                // list is chosen while the record is built, which has not happened. `set_allele_map`
                // supplies it at render time.
                static const vector<int> no_allele_map;
                const vector<int>& trav_to_allele_vec = no_allele_map;
                const pair<string, size_t> key =
                    site_ref_key(pr.snarl, pr.ref_path_name, pr.ref_offset);
                int called_i = use_genotype.empty() ? -1 : use_genotype[0];
                int called_j = use_genotype.size() > 1 ? use_genotype[1] : called_i;
                vector<int> panel = panel_alleles(graph, pr.travs);
                if (!linkage_collector->respecify(pr.record_key,
                                                  key.first, key.second,
                                                  used->genotype_lls, panel,
                                                  called_i, called_j, trav_to_allele_vec,
                                                  (size_t)copies, copies == 1, pr.parent_trav,
                                                  pr.parent_crossing, /*emitted*/ false)) {
                    // respecify refuses a site whose compact space it cannot build (no called
                    // traversal, no likelihoods, or more than 127 reachable alleles). If such a site
                    // is already in the layer, its entry describes the line this barrier pass has
                    // just tombstoned, and its allele numbering belongs to that dead line -- so
                    // leaving it would patch the *replacement* with the old numbering, which is how
                    // a GT naming an allele the record has no ALT for gets written. Dropped instead:
                    // an unpatched record is a per-site call, which is a correct answer, where a
                    // mis-numbered one is not a VCF.
                    if (linkage_collector->has_entry(pr.record_key)) {
                        linkage_collector->retract(pr.record_key);
                        ++stale_respecify;
                    } else {
                        // Not previously recorded -- a chain nothing was written for during the
                        // sweep -- so it joins the layer now rather than being revised. The contig
                        // and POS are the ones the record was just written with, post-flatten: the
                        // same key add_variant filed the line under, so write_variants' lookup can
                        // find the entry. get_ref_position reproduces the pre-flatten POS and left
                        // every gained record unphased.
                        linkage_collector->record(
                            key.first, key.second,
                            used->genotype_lls, panel, called_i, called_j, trav_to_allele_vec,
                            pr.record_key, 1.0,
                            (size_t)copies, pr.snarl.start().node_id(), pr.snarl.end().node_id(),
                            copies == 1, pr.parent_record_key,
                            pr.parent_trav, pr.parent_crossing, pr.generation, /*emitted*/ false);
                    }
                }
            }
            if (was_gained) {
                ++gained;
            } else {
                ++revised;
            }

            // The record this chain's children were masked against has just changed -- or, for a
            // gained chain, exists for the first time -- so their crossing masks are recomputed
            // from the mapping the emit above produced. The sweep-time masks were built from
            // whatever this thread had last emitted, which for children of a retained chain was a
            // foreign snarl's mapping.
            for (PendingRecord& child : pending) {
                if (child.parent_record_key == pr.record_key) {
                    bool known = true;
                    child.parent_crossing = child_crossing_mask(pr.travs, child.snarl, &known);
                    child.crossing_known = known;
                }
            }
        }
        if (linkage_collector != nullptr) {
            generations = max(generations, linkage_collector->max_generation());
        }
    }
    if (show_progress) {
        // Exact, not estimated, and deterministic -- which matters because peak RSS cannot resolve a
        // delta this size: six runs of one binary on chr20 spread 3.39 to 4.42 GB, wider than the
        // retention itself. Two independent sizing estimates disagreed by 1.9x and both were
        // guesses; this walks the objects.
        size_t retained_bytes = 0, retained_visits = 0, retained_gls = 0;
        auto measure = [&](const PendingRecord& rec) {
            retained_bytes += sizeof(PendingRecord) + rec.ref_path_name.capacity()
                              + rec.genotype.capacity() * sizeof(int);
            retained_bytes += rec.travs.capacity() * sizeof(SnarlTraversal);
            for (const SnarlTraversal& t : rec.travs) {
                retained_visits += (size_t)t.visit_size();
                retained_bytes += (size_t)t.visit_size() * sizeof(Visit);
            }
            const auto* rl = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (rl != nullptr) {
                for (const auto& kv : rl->genotype_lls) {
                    ++retained_gls;
                    retained_bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                }
                if (rl->alt_ploidy_info != nullptr) {
                    for (const auto& kv : rl->alt_ploidy_info->genotype_lls) {
                        ++retained_gls;
                        retained_bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                }
            }
        };
        for (const auto& queue : render_records) {
            for (const auto& rec : queue) {
                measure(rec);
            }
        }
        for (const auto& rec : pending) {
            measure(rec);
        }
        // "Not revised by the barrier" rather than "top-level": recurse-on-fail reaches children with
        // no ploidy override, so they take the same path. On chr20 that is 165,408 top-level snarls
        // plus 26,799 such children.
        // Stage 14's off-ramp, printed so the decision rests on numbers rather than on whether
        // the frame sounds like it should matter. The criteria were fixed before the measurement:
        // if under 1% of adjacent pairs reorder AND 99% of gap ratios sit inside 1.05, the
        // haplotype frame buys nothing measurable and stage 15 drops its distance half.
        if (g_frame_pairs.load() > 0) {
            const double reorder_pct = 100.0 * (double)g_frame_reordered.load()
                                       / (double)g_frame_pairs.load();
            const size_t gm = g_frame_gap_measured.load();
            const double within_pct = gm ? 100.0 * (double)g_frame_gap_within_5pct.load() / (double)gm
                                         : 0.0;
            cerr << "[vg call] haplotype frame: " << g_frame_pairs.load()
                 << " adjacent sibling pairs, " << g_frame_reordered.load()
                 << " reorder (" << reorder_pct << "%); " << gm
                 << " gaps measured, " << within_pct << "% within 1.05 of the reference gap"
                 << endl;
            cerr << "[vg call] haplotype frame: " << g_frame_disagree.load()
                 << " children whose two parent traversals disagree on offset (mean |delta| "
                 << (g_frame_disagree.load()
                         ? (double)g_frame_offset_delta_sum.load() / (double)g_frame_disagree.load()
                         : 0.0)
                 << " bp, max " << g_frame_offset_delta_max.load() << "); "
                 << g_frame_no_offset.load() << " with no offset on any called traversal" << endl;
        }
        cerr << "[vg call] retained for rendering: " << render_record_count()
             << " snarls the barrier will not revise, plus " << pending.size()
             << " nested chains; " << (retained_bytes / (1024.0 * 1024.0)) << " MB over "
             << retained_visits << " traversal visits and " << retained_gls
             << " genotype likelihoods" << endl;
        cerr << "[vg call] frames: " << frame_written << " nested chains measured along their"
             << " settled parent traversal; not written: " << frame_no_single_trav
             << " with no single carrying traversal, " << frame_no_parent
             << " whose parent record was not found, " << frame_not_crossed
             << " the settled traversal does not cross, " << frame_no_entry
             << " with no layer entry" << endl;
        cerr << "[vg call] single sweep: " << pending.size() << " nested chains retained over "
             << (generations + 1) << " generations; " << revised << " revised, " << gained
             << " reachable only under the settled parent, " << retracted << " retracted";
        if (crossing_unknown > 0) {
            cerr << ", " << crossing_unknown << " with a crossing mask the sweep could not compute";
        }
        if (stale_respecify > 0) {
            cerr << ", " << stale_respecify << " dropped from the layer because the replacement "
                 << "record could not be respecified";
        }
        cerr << endl;
    }

    // Hand every surviving chain to the render pass. This is what makes the two populations one:
    // top-level records were already staged and rendered from the settled genotype, and now nested
    // ones are too, so there is a single place a line is written and a single genotype it is written
    // from. A dropped chain is not handed over -- its parent's settled genotype does not carry it,
    // so the sample has no copy of it and it is not a record.
    //
    // Appended round-robin rather than all onto one queue: the render pass is parallel over queues,
    // and 30,416 chains on one thread would serialise it.
    size_t next_queue = 0;
    for (PendingRecord& pr : pending) {
        if (pr.dropped) {
            continue;
        }
        if (render_records.empty()) {
            render_records.resize(max((size_t)1, (size_t)get_thread_count()));
        }
        render_records[next_queue % render_records.size()].push_back(std::move(pr));
        ++next_queue;
    }
    pending.clear();
}

bool FlowCaller::call_snarl_internal(const Snarl& managed_snarl,
                                      const string& parent_ref_path_name,
                                      pair<size_t, size_t> parent_ref_interval,
                                      const ChildTraversalSets* parent_child_trav_sets,
                                    int ploidy_override) {


    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash.
    Snarl snarl = managed_snarl;

    // Staged in the nested branch below and completed after the descent loop, because the loop reads
    // `travs` and this record is what finally takes ownership of it.
    unique_ptr<PendingRecord> pending_this;
    // The same staging, for a snarl the barrier will not revise. Completed at the same point and for
    // the same reason: `travs` cannot be moved until descent has finished reading it.
    unique_ptr<PendingRecord> render_this;
    // Whether THIS invocation ran emit_variant, so the descent block below knows the thread-local
    // last_emitted describes this snarl. A retained chain (and a GAF run) skips the emit, and the
    // stale mapping -- some other snarl's -- must not be used for its children's crossing masks.
    bool emitted_this_call = false;

#ifdef debug
    cerr << "call_snarl_internal on " << pb2json(snarl) << " with parent_ref_path=" << parent_ref_path_name
         << " parent_child_trav_sets=" << (parent_child_trav_sets ? "provided" : "null") << endl;
#endif

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }

    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.deep_contents(&snarl, graph, false);
        if (snarl_contents.second.size() > max_snarl_edges) {
            // size cap needed as non-nested FlowCaller doesn't handle large snarls
            return false;
        }        
        size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
        size_t length = 0;
        for (auto i = snarl_contents.first.begin(); i != snarl_contents.first.end() && length < len_threshold; ++i) {
            length += graph.get_length(graph.get_handle(*i));
        }
        greedy_avg_flow = length > len_threshold;
    }
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    if (common_names.empty()) {
        // No reference path through snarl
        // If we have parent context, we can still process using parent's ref path
        if (parent_child_trav_sets == nullptr || parent_ref_path_name.empty()) {
#ifdef debug
            cerr << "  -> returning false: no common ref path and no parent context" << endl;
#endif
            return false;
        }
#ifdef debug
        cerr << "  -> using parent ref path: " << parent_ref_path_name << endl;
#endif
    }

    // Use parent's ref path if no direct path, otherwise prefer base reference over gref paths
    string ref_path_name;
    if (common_names.empty()) {
        ref_path_name = parent_ref_path_name;
    } else {
        // Prefer base reference paths over derived gref paths.  Test the whole gref
        // namespace, not just the fragment suffix: a gref copy of the reference sorts
        // before the path it was copied from (gref_x < x).
        // common_names is sorted, so we iterate to find first non-gref path
        ref_path_name = common_names.front();  // default to first (lexicographically smallest)
        for (const string& name : common_names) {
            if (!GrefCover::is_gref_derived(name)) {
                ref_path_name = name;
                break;
            }
        }
    }

    // find the reference traversal and coordinates using the path position graph interface
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
    bool use_parent_interval = false;

    if (common_names.empty() && parent_child_trav_sets != nullptr) {
        // No direct reference path - use parent's interval and traversals directly
        ref_interval = make_tuple(parent_ref_interval.first, parent_ref_interval.second, false, step_handle_t(), step_handle_t());
        use_parent_interval = true;
    } else {
        ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        if (get<0>(ref_interval) == -1) {
            // could not find reference path interval consistent with snarl due to orientation conflict
            return false;
        }
        if (get<2>(ref_interval) == true) {
            // calling code assumes snarl forward on reference
            flip_snarl(snarl);
            ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        }
    }

    SnarlTraversal ref_trav;

    if (!use_parent_interval) {
        // Build reference traversal from path steps
        step_handle_t cur_step = get<3>(ref_interval);
        step_handle_t last_step = get<4>(ref_interval);
        if (get<2>(ref_interval)) {
            std::swap(cur_step, last_step);
        }
        bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));

        while (true) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(cur_handle));
            visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
            if (graph.get_id(cur_handle) == snarl.end().node_id()) {
                break;
            } else if (get<2>(ref_interval) == true) {
                if (!graph.has_previous_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_previous_step(cur_step);
            } else {
                if (!graph.has_next_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_next_step(cur_step);
            }
            // todo: we can compute flow at the same time
        }
        assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());
    }
    // If use_parent_interval, ref_trav stays empty - we'll use first parent traversal as pseudo-reference

    vector<SnarlTraversal> travs;
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow);
        travs = std::move(weighted_travs.first);
    } else {
        // find the traversals using the generic interface
        travs = traversal_finder.find_traversals(snarl);
    }

    if (travs.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any traversals through snarl " << pb2json(managed_snarl) << endl;
        return false;
    }
#ifdef debug
    cerr << "  found " << travs.size() << " traversals, use_parent_interval=" << use_parent_interval << endl;
#endif

    // optional traversal length clamp can, ex, avoid trying to resolve a giant snarl    
    if (allele_length_range.first > 0 || allele_length_range.second < numeric_limits<size_t>::max()) {
        size_t max_trav_len = 0;
        for (const SnarlTraversal & trav : travs) {
            size_t trav_len = 0;
            for (size_t i = 1; i < trav.visit_size() - 1; ++i) {
                trav_len += graph.get_length(graph.get_handle(trav.visit(i).node_id()));
            }
            max_trav_len = max(max_trav_len, trav_len);
            if (max_trav_len > allele_length_range.second) {
                return false;
            }
        }
        if (max_trav_len < allele_length_range.first) {
            return false;
        }
    }

    // find the reference traversal in the list of results from the traversal finder
    int ref_trav_idx = -1;

    if (use_parent_interval) {
        // No direct reference path - use first traversal from first non-empty set as pseudo-reference
        if (parent_child_trav_sets != nullptr) {
            for (const auto& tset : *parent_child_trav_sets) {
                if (!tset.empty()) {
                    const SnarlTraversal& first_trav = tset[0];
                    for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
                        if (travs[i] == first_trav) {
                            ref_trav_idx = i;
                        }
                    }
                    if (ref_trav_idx < 0 && first_trav.visit_size() > 0) {
                        ref_trav_idx = travs.size();
                        travs.push_back(first_trav);
                    }
                    break;
                }
            }
        }
        // If still -1, use first traversal from finder (or 0 if none)
        if (ref_trav_idx < 0) {
            ref_trav_idx = travs.empty() ? -1 : 0;
        }
    } else {
        for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
            // todo: is there a way to speed this up?
            if (travs[i] == ref_trav) {
                ref_trav_idx = i;
            }
        }

        if (ref_trav_idx == -1) {
            ref_trav_idx = travs.size();
            // we didn't get the reference traversal from the finder, so we add it here
            travs.push_back(ref_trav);
        }
    }

    bool ret_val = true;
    vector<int> trav_genotype;  // Declared outside block so we can pass to children
    // A propagated ploidy wins over the contig's or the region BED's: it says how many called
    // parent alleles actually reach this child, which is the number of copies present here.
    int ploidy = ploidy_override >= 0
                 ? ploidy_override
                 : ploidy_at(ref_path_name, get<0>(ref_interval),
                             ref_offsets.count(ref_path_name) ? ref_offsets.at(ref_path_name) : 0,
                             ref_ploidies[ref_path_name]);

    // Constants for bounded traversal set handling
    const int MAX_TRAVS_PER_SET = 10;

    if (traversals_only) {
        assert(gaf_output);
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else if (parent_child_trav_sets != nullptr && !parent_child_trav_sets->empty()) {
        // Genotype using bounded search over traversal sets from parent
        // Each set contains traversals consistent with one parent allele
        ploidy = parent_child_trav_sets->size();

        // Track which set each traversal index belongs to (for phase consistency)
        // set_membership[i] = which parent allele set traversal i came from, or -1 if from finder
        vector<int> set_membership(travs.size(), -1);

        // Merge traversals from sets into travs, tracking membership
        // Also rank by support and keep top MAX_TRAVS_PER_SET per set
        vector<vector<int>> set_to_trav_indices(ploidy);  // indices into travs for each set

        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            const TraversalSet& tset = (*parent_child_trav_sets)[set_idx];

            if (tset.empty()) {
                // Empty set means parent allele doesn't traverse this child (star allele)
                continue;
            }

            // Add traversals from this set to travs (avoiding duplicates)
            // Keep track of indices for this set
            vector<pair<double, int>> support_and_idx;  // (support, index in travs)

            for (const SnarlTraversal& trav : tset) {
                // Check if this traversal already exists in travs
                int match_idx = -1;
                for (int i = 0; i < travs.size() && match_idx < 0; ++i) {
                    if (travs[i] == trav) {
                        match_idx = i;
                    }
                }

                if (match_idx < 0) {
                    // New traversal - add it
                    match_idx = travs.size();
                    travs.push_back(trav);
                    set_membership.push_back(set_idx);
                } else if (set_membership[match_idx] < 0) {
                    // Traversal was from finder, now claim it for this set
                    set_membership[match_idx] = set_idx;
                }
                // Note: if already claimed by another set, that's fine (shared region)

                // Get support for ranking
                double support = TraversalSupportFinder::support_val(
                    support_finder.get_traversal_support(trav));
                support_and_idx.push_back({support, match_idx});
            }

            // Sort by support (descending) and keep top MAX_TRAVS_PER_SET
            std::sort(support_and_idx.begin(), support_and_idx.end(),
                      [](const auto& a, const auto& b) { return a.first > b.first; });

            for (int i = 0; i < std::min((int)support_and_idx.size(), MAX_TRAVS_PER_SET); ++i) {
                set_to_trav_indices[set_idx].push_back(support_and_idx[i].second);
            }
        }

        // Which parent haplotypes actually traverse this child? A parent allele with
        // an empty traversal set skips the child entirely, and gets a star or
        // missing allele rather than a genotype.
        vector<int> traversing_sets;
        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            if (!(*parent_child_trav_sets)[set_idx].empty()) {
                traversing_sets.push_back(set_idx);
            }
        }

        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        int marker = star_allele ? STAR_ALLELE_MARKER : MISSING_ALLELE_MARKER;

        if (traversing_sets.empty()) {
            // No parent allele traverses this child at all.
            trav_genotype.assign(ploidy, marker);
        } else {
            // Genotype at the ploidy that actually traverses the site, not at the
            // parent's ploidy.
            //
            // This used to ask for a full-ploidy genotype and then overwrite the
            // positions belonging to empty sets, which was wrong twice over. A site
            // only one haplotype reaches is not diploid, and asking a genotyper for
            // a diploid call there lets a spurious heterozygote absorb noise on a
            // second allele for free, biasing which allele gets picked before any
            // marker is applied. And because genotype() returns a sorted allele
            // multiset with no haplotype identity, overwriting position i had no
            // relationship to which parent haplotype was actually empty -- so which
            // allele got discarded was effectively arbitrary.
            int effective_ploidy = (int)traversing_sets.size();
            vector<int> called_alleles;
            std::tie(called_alleles, trav_call_info) = snarl_caller.genotype(
                snarl, travs, ref_trav_idx, effective_ploidy, ref_path_name,
                make_pair(get<0>(ref_interval), get<1>(ref_interval)));

            // Scatter the called alleles back onto the traversing haplotypes,
            // leaving the others as star/missing.
            trav_genotype.assign(ploidy, marker);
            for (size_t j = 0; j < traversing_sets.size() && j < called_alleles.size(); ++j) {
                trav_genotype[traversing_sets[j]] = called_alleles[j];
            }
        }

        // Emit variant with selected genotype
        bool added = true;

        // Only emit VCF if snarl is on reference path
        if (use_parent_interval) {
            added = true;
        } else if (!gaf_output) {
            // Staged, not emitted: `render_retained_records` writes it after the sweep. `added` is
            // the value emit_variant would have returned, and the only caller that reads it here is
            // the ret_val below, which gates recursion -- so it must not become "the line was
            // written", which is not known yet. A staged record is a record that will be written.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_path_name,
                        ref_offsets[ref_path_name]);
            render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                              ref_path_name, ref_offsets[ref_path_name], ploidy,
                                              true);
            added = render_this != nullptr;
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offsets[ref_path_name],
                                     genotype_snarls, ploidy);
            }
            emitted_this_call = true;
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    } else if (ploidy_override >= 0) {
        // A nested chain: reached by descent, with the ploidy its parent implied.
        //
        // Only a nested chain can have its ploidy revised at the barrier, so only a nested chain
        // needs the other ploidy's answer computed and kept.
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(true);
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(
            snarl, travs, ref_trav_idx, ploidy, ref_path_name,
            make_pair(get<0>(ref_interval), get<1>(ref_interval)));
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(false);

        const bool retain_only = nested_context.retain_only;

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
        bool added = true;
        if (retain_only) {
            // No called parent allele reaches this chain, so nothing about it may reach the VCF --
            // yet. It is genotyped and kept because linkage can still move the parent onto an allele
            // that does reach it, and going back to the reads to find that out is the cost this
            // design exists to remove.
            added = true;
        } else if (!gaf_output) {
            // Recorded here rather than inside emit_variant, and deliberately NOT on the retain_only
            // path above: a retained chain is recorded today only if the barrier later gives it a
            // line, via the fallback record() there. Recording it now would make respecify succeed
            // where it currently falls through, which moves the `gained` count -- a real change, and
            // one for stage 10 to make on purpose rather than for this motion to make by accident.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_path_name,
                        ref_offsets[ref_path_name]);
            // Staged, not emitted -- the same discipline the top-level branch already follows, and
            // the reason both of stage 10's residual counts were non-zero. Emitting here writes the
            // line from a genotype the barrier has not settled yet, and once a line exists the only
            // way to change it is a patch: a patch cannot add an ALT (so a settled genotype naming a
            // traversal the line has no ALT for is dropped and counted -- the 496 `unrenderable`
            // events) and it cannot withdraw a line (so a site that settles on the reference keeps a
            // record with GT 0/0 -- the 1,383). Both counters sit behind `if (!e.emitted) continue`
            // in the resolver, so neither can fire at all once nothing is written before the barrier.
            //
            // `added` is what emit_variant would have returned, and the only thing that reads it is
            // the ret_val below, which gates recursion -- so a staged record counts as one that will
            // be written, exactly as at top level.
            added = defer_nested_descent && !pending_records.empty();
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offsets[ref_path_name],
                                     genotype_snarls, ploidy);
                emitted_this_call = true;
            }
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name,
                                                              ref_offsets[ref_path_name]);
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx,
                             pos_info.first, pos_info.second, &support_finder);
        }

        // Kept for the barrier, but staged rather than stored: `travs` must not be moved out here,
        // because the descent loop further down still reads it to work out which children the called
        // alleles reach. Moving it here emptied it before that loop ran, so every child came back with
        // a copy number of zero and whole subtrees were held back -- 2,494 records off chr20, every one
        // of them nested, top-level untouched. The move happens once descent is done with it.
        if (defer_nested_descent && !pending_records.empty()) {
            pending_this.reset(new PendingRecord());
            pending_this->snarl = snarl;
            pending_this->ref_path_name = ref_path_name;
            pending_this->ref_offset = ref_offsets[ref_path_name];
            pending_this->ref_trav_idx = ref_trav_idx;
            pending_this->genotype = trav_genotype;
            pending_this->ploidy = ploidy;
            pending_this->record_key = record_key_of(snarl);
            pending_this->parent_record_key = nested_context.parent_record_key;
            pending_this->parent_crossing = nested_context.parent_crossing;
            pending_this->crossing_known = nested_context.crossing_known;
            pending_this->parent_trav = nested_context.parent_trav;
            pending_this->generation = (uint8_t)min(current_generation, (size_t)255);
            // Nothing is written for a nested chain during the sweep any more, so there is no
            // line to record the whereabouts of and nothing for the barrier to retract or replace.
            // The record is written once, by the render pass, from the settled genotype.
            pending_this->emitted = false;
            pending_this->call_info = std::move(trav_call_info);
        }
        ret_val = trav_genotype.size() == ploidy && added;
    } else {
        // Top-level snarl or no parent context - genotype from scratch using support
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, ref_path_name,
                                                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)));

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        bool added = true;
        if (!gaf_output) {
            // Staged, not emitted: `render_retained_records` writes it after the sweep. `added` is
            // the value emit_variant would have returned, and the only caller that reads it here is
            // the ret_val below, which gates recursion -- so it must not become "the line was
            // written", which is not known yet. A staged record is a record that will be written.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_path_name,
                        ref_offsets[ref_path_name]);
            render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                              ref_path_name, ref_offsets[ref_path_name], ploidy,
                                              true);
            added = render_this != nullptr;
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offsets[ref_path_name],
                                     genotype_snarls, ploidy);
            }
            emitted_this_call = true;
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    }

    // Symbolic nested calling: descend into each child the called alleles actually reach, at the
    // ploidy they reach it with.
    //
    // This is deliberately not gated on ret_val. Recursion today is a side effect of emission --
    // emit_variant returns true even when it wrote nothing -- so a snarl genotyped hom-ref reports
    // success and its children are never queued. That is what buries nested variants under a parent
    // that looked resolved. Here descent is a decision about ploidy, not about whether a line was
    // written, and a symbolically-reference parent is precisely the case where descending matters
    // most.
    //
    // Children are genotyped independently, with no parent traversal sets: constraining a child to
    // the parent's called alleles is what --top-down does, and it measured worse than the default
    // on every axis including recall, because a child's true allele is then unreachable whenever
    // the parent's call is imperfect.
    // Only where the call succeeded. The driver still descends into the children of a *failed*
    // snarl under RecurseOnFail, so firing here as well would genotype those children twice; and a
    // failed snarl has no genotype to derive a child ploidy from anyway. The case this exists for --
    // a parent called hom-ref, which reports success and today ends the descent -- is covered,
    // because that is a success.
    if (ret_val && symbolic_manager != nullptr && !trav_genotype.empty() &&
        parent_child_trav_sets == nullptr) {
        const Snarl* managed_ptr = snarl_manager.into_which_snarl(snarl.start().node_id(),
                                                                  snarl.start().backward());
        if (managed_ptr != nullptr) {
            // Snapshotted before the loop: each child call runs emit_variant of its own and
            // overwrites the thread's copy, so reading it inside the loop would describe the
            // previous child rather than this parent.
            EmittedAlleles parent_alleles = last_emitted;
            if (!emitted_this_call) {
                // This snarl never went through emit_variant (a retained chain, or a GAF run), so
                // last_emitted still holds some other snarl's mapping. Applying a foreign map to
                // this snarl's traversals produced semantically garbage crossing masks that the
                // barrier then used to gain, drop, and re-ploidy this snarl's descendants. The
                // masks are recomputed at the barrier if this chain is ever actually emitted.
                parent_alleles.valid = false;
            }
            // Deferral turns on exactly one question: can linkage still move this snarl's genotype?
            // Only a snarl that reached the linkage layer can be rewritten, so one with no entry has
            // a final genotype already and its children are visited now, as they always were. That
            // keeps 69% of the children called at ploidy 2 out of the barrier at no cost in
            // coherence, and it is what makes the deferred population the 44% that can actually move.
            // `emit_phasing` belongs in the test even though the option layer refuses the one
            // configuration that would fail it: deferral reads the settled allele pair out of the
            // phasing, so without phasing a deferred child is one whose parent cannot be found, and
            // `descend_pending` would drop it. Kept here so the invariant is a property of this code
            // rather than of a check somewhere else.
            // Stage 14, inert: measure whether the haplotype frame would reorder or re-space
            // these children relative to the reference frame the model currently uses. Done once
            // per parent, before the descent loop, so it costs one pass over the children and no
            // extra graph queries beyond node lengths.
            if (!travs.empty() && ref_trav_idx >= 0 && (size_t)ref_trav_idx < travs.size()) {
                // (offset along the reference, offset along a called traversal) per child that both
                // frames can place.
                vector<pair<int64_t, int64_t>> both;
                for (const Snarl* child : snarl_manager.children_of(managed_ptr)) {
                    if (child == nullptr || snarl_manager.is_trivial(child, graph)) {
                        continue;
                    }
                    int64_t ref_off = 0, ref_span = 0;
                    if (!traversal_offset_span(graph, travs[ref_trav_idx], *child, &ref_off,
                                               &ref_span)) {
                        continue;
                    }
                    // Along each called traversal. A diploid parent gives two, and where they
                    // disagree the frame is genuinely ambiguous -- which is the case stage 15's
                    // consensus rule has to settle.
                    int64_t first = -1, second = -1;
                    for (size_t g = 0; g < trav_genotype.size(); ++g) {
                        const int t = trav_genotype[g];
                        if (t < 0 || (size_t)t >= travs.size()) {
                            continue;
                        }
                        int64_t off = 0, span = 0;
                        if (!traversal_offset_span(graph, travs[t], *child, &off, &span)) {
                            continue;
                        }
                        if (first < 0) {
                            first = off;
                        } else {
                            second = off;
                        }
                    }
                    if (first < 0) {
                        ++g_frame_no_offset;
                        continue;
                    }
                    if (second >= 0 && second != first) {
                        ++g_frame_disagree;
                        const size_t d = (size_t)llabs(second - first);
                        g_frame_offset_delta_sum += d;
                        size_t prev = g_frame_offset_delta_max.load();
                        while (d > prev && !g_frame_offset_delta_max.compare_exchange_weak(prev, d)) {
                        }
                    }
                    both.emplace_back(ref_off, first);
                }
                // Sort by the reference frame, then ask whether the haplotype frame agrees about
                // the order of each adjacent pair and about how far apart they are.
                sort(both.begin(), both.end());
                for (size_t i = 1; i < both.size(); ++i) {
                    ++g_frame_pairs;
                    if (both[i].second < both[i - 1].second) {
                        ++g_frame_reordered;
                    }
                    const int64_t ref_gap = both[i].first - both[i - 1].first;
                    const int64_t hap_gap = llabs(both[i].second - both[i - 1].second);
                    if (ref_gap > 0 && hap_gap > 0) {
                        ++g_frame_gap_measured;
                        const double ratio = (double)hap_gap / (double)ref_gap;
                        if (ratio >= 1.0 / 1.05 && ratio <= 1.05) {
                            ++g_frame_gap_within_5pct;
                        }
                    }
                }
            }

            for (const Snarl* child : snarl_manager.children_of(managed_ptr)) {
                if (child == nullptr || snarl_manager.is_trivial(child, graph)) {
                    continue;
                }
                // v1 descends only where the reference also goes. A chain crossed only by a
                // non-reference allele has no reference path through it, so REF and POS for its
                // record are undefined; --nested-pseudo-ref is where that will be handled.
                if (ref_trav_idx >= 0 && ref_trav_idx < (int)travs.size()) {
                    vector<int> ref_only(1, ref_trav_idx);
                    if (child_ploidy(travs, ref_only, *child, 1) == 0) {
                        ++g_descent_skipped_no_ref;
                        continue;
                    }
                }

                // The exactly-once rule. Under block emission a chain crossed by every called
                // haplotype only inside a difference block has already been spelled out by that
                // block's ALT, so its own record would report the same variation a second time.
                // Inert with the flag off, and deliberately inert for a snarl whose projection has
                // no symbols, where the question cannot be answered.
                if (chain_reported_inline(snarl, travs, trav_genotype, ref_trav_idx, *child)) {
                    continue;
                }

                int copies = child_ploidy(travs, trav_genotype, *child, ploidy);
                bool retain_only = nested_context.retain_only;
                if (copies <= 0) {
                    // No called allele reaches it *yet*. Visited anyway, while the reads for this
                    // window are still resident, because the parent's genotype is not settled and
                    // linkage may move it onto an allele that does reach this chain -- 296 of them on
                    // chr20. Going back to the reads at the barrier to find out is what cost five
                    // sweeps of the contig. Nothing about it is emitted unless the barrier says so.
                    ++g_descent_skipped_no_copy;
                    if (!defer_nested_descent) {
                        continue;   // without retention there is nothing to come back to
                    }
                    retain_only = true;
                }

                // A haploid child hangs off exactly one of the parent's called traversals, and
                // which one decides the strand its allele sits on once the parent is phased. The
                // traversal is recorded, not its index in the genotype: `record` sorts the pair and
                // the Viterbi then orients it against the panel, so an index recorded here means
                // nothing by the time the child is placed, while the traversal still names the same
                // path through the parent.
                int carrying_trav = -1;
                if (copies == 1) {
                    for (size_t g = 0; g < trav_genotype.size(); ++g) {
                        vector<int> one(1, trav_genotype[g]);
                        if (child_ploidy(travs, one, *child, 1) == 1) {
                            carrying_trav = trav_genotype[g];
                            break;
                        }
                    }
                }
                // Saved and restored rather than assigned: a child may descend further, and its own
                // children must see *it* as their parent, not this snarl.
                NestedContext saved = nested_context;
                nested_context.active = (copies == 1);
                nested_context.parent_record_key = record_key_of(snarl);
                nested_context.parent_trav = carrying_trav;
                nested_context.retain_only = retain_only;
                bool crossing_known = parent_alleles.valid;
                // No dependence on the parent having emitted anything: the mask is over this
                // snarl's own candidate traversals, which exist whether or not a line was written.
                // That is what made the old mask unavailable for a collapsed parent.
                nested_context.parent_crossing =
                    child_crossing_mask(travs, *child, &crossing_known);
                nested_context.crossing_known = crossing_known;
                size_t saved_generation = current_generation;
                current_generation = saved_generation + 1;
                ++g_descent_depth;
                if (g_descent_depth < 16) {
                    ++g_descent_depth_hist[g_descent_depth];
                }
                // Falls back to the PARENT's ploidy, not to a literal 2. `copies` is zero here
                // only for a chain no called parent allele reaches -- a retained chain, genotyped
                // now because linkage may still move the parent onto an allele that does reach it --
                // and it has to be genotyped at some ploidy to have an answer at all. Two is the
                // parent's ploidy only on a diploid contig. Under `-d 1` (chrY, or chrX outside the
                // pseudoautosomal regions) it is a ploidy the contig does not have, and a child
                // cannot carry more copies than its parent: `child_ploidy` caps at exactly this
                // value, so the fallback should agree with the cap rather than exceed it.
                //
                // Either ploidy still yields both answers -- `set_want_alt_ploidy` computes the
                // other one and the barrier keeps it -- so this changes which of the two is the
                // primary, not whether the barrier can re-ploidy the chain later.
                call_snarl_internal(*child, ref_path_name,
                                    make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                    nullptr, copies >= 1 ? copies : ploidy);
                --g_descent_depth;
                current_generation = saved_generation;
                nested_context = saved;
            }
        }
    }


    // In nested mode, recursively call child snarls
    if (nested && !trav_genotype.empty()) {
        // Find the managed snarl pointer so we can get its children
        const Snarl* managed_ptr = snarl_manager.into_which_snarl(snarl.start().node_id(), snarl.start().backward());
        if (managed_ptr) {
            const vector<const Snarl*>& children = snarl_manager.children_of(managed_ptr);
            for (const Snarl* child : children) {
                if (child && !snarl_manager.is_trivial(child, graph)) {
                    // Build ChildTraversalSets: one set per parent allele
                    // Each set contains all traversals through child consistent with that parent allele
                    ChildTraversalSets child_trav_sets;
                    bool any_real_traversals = false;

                    for (int allele_idx : trav_genotype) {
                        if (allele_idx >= 0 && allele_idx < travs.size()) {
                            // Find all traversals through child consistent with this parent traversal
                            TraversalSet tset = find_child_traversal_set(travs[allele_idx], *child);
                            if (!tset.empty()) {
                                any_real_traversals = true;
                            }
                            child_trav_sets.push_back(std::move(tset));
                        } else {
                            // Star/missing allele - pass empty set
                            child_trav_sets.push_back(TraversalSet());
                        }
                    }

                    // If no genotyped alleles traverse the child, skip it
                    if (!any_real_traversals) {
                        continue;
                    }

                    // Recursively call child with traversal sets
                    call_snarl_internal(*child, ref_path_name,
                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                        &child_trav_sets);
                }
            }
        }
    }

    // Both traversal readers are behind us now -- symbolic descent above, and the `-A` recursion
    // just above, which builds each child's ChildTraversalSets out of `travs[allele_idx]`. Only here
    // can the staged record take them. Completing it one block earlier broke four -A and --top-down
    // tests; completing it at emit time cost 12,302 chr20 records. At most one of these is set: a
    // snarl is either revised by the barrier or it is not.
    if (pending_this != nullptr) {
        pending_this->travs = std::move(travs);
        pending_records[omp_get_thread_num()].push_back(std::move(*pending_this));
        pending_this.reset();
    } else if (render_this != nullptr) {
        render_this->travs = std::move(travs);
        render_records[omp_get_thread_num()].push_back(std::move(*render_this));
        render_this.reset();
    }


    return ret_val;
}

string FlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

NestedFlowCaller::NestedFlowCaller(const PathPositionHandleGraph& graph,
                                   SupportBasedSnarlCaller& snarl_caller,
                                   SnarlManager& snarl_manager,
                                   const string& sample_name,
                                   TraversalFinder& traversal_finder,
                                   const vector<string>& ref_paths,
                                   const vector<size_t>& ref_path_offsets,
                                   const vector<int>& ref_path_ploidies,
                                   AlignmentEmitter* aln_emitter,
                                   bool traversals_only,
                                   bool gaf_output,
                                   size_t trav_padding,
                                   bool genotype_snarls) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    nested_support_finder(dynamic_cast<NestedCachedPackedTraversalSupportFinder&>(snarl_caller.get_support_finder())){

    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }

}
   
NestedFlowCaller::~NestedFlowCaller() {

}

bool NestedFlowCaller::call_snarl(const Snarl& managed_snarl) {
    
    // remember the calls for each child snarl in this table
    CallTable call_table;

    bool called = call_snarl_recursive(managed_snarl, -1, "", make_pair(0, 0), call_table);

    if (called) { 
        emit_snarl_recursive(managed_snarl, -1, call_table);
    }

    return called;
}

bool NestedFlowCaller::call_snarl_recursive(const Snarl& managed_snarl, int max_ploidy,
                                            const string& parent_ref_path_name, pair<size_t, size_t> parent_ref_interval,
                                            CallTable& call_table) {

    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash.
    Snarl snarl = managed_snarl;

    // hook into our table entry
    CallRecord& record = call_table[managed_snarl];
    
    // get some reference information if possible
    // todo: make a function
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    string ref_path_name;
    SnarlTraversal ref_trav;
    int ref_trav_idx = -1;
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
    string gt_ref_path_name;
    pair<size_t, size_t> gt_ref_interval;
    
    if (!common_names.empty()) {
        // Prefer base reference paths over derived gref paths.  Test the whole gref
        // namespace, not just the fragment suffix: a gref copy of the reference sorts
        // before the path it was copied from (gref_x < x).
        ref_path_name = common_names.front();  // default to first (lexicographically smallest)
        for (const string& name : common_names) {
            if (!GrefCover::is_gref_derived(name)) {
                ref_path_name = name;
                break;
            }
        }

        // find the reference traversal and coordinates using the path position graph interface
        ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        if (get<0>(ref_interval) == -1) {
            // no reference path found due to orientation conflict
            return false;
        }
        if (get<2>(ref_interval) == true) {
            // calling code assumes snarl forward on reference
            flip_snarl(snarl);
            ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        }

        step_handle_t cur_step = get<3>(ref_interval);
        step_handle_t last_step = get<4>(ref_interval);
        if (get<2>(ref_interval)) {
            std::swap(cur_step, last_step);
        }
        bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));
    
        while (true) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(cur_handle));
            visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
            if (graph.get_id(cur_handle) == snarl.end().node_id()) {
                break;
            } else if (get<2>(ref_interval) == true) {
                if (!graph.has_previous_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(snarl) << endl;
                    return false;
                }
                cur_step = graph.get_previous_step(cur_step);
            } else {
                if (!graph.has_next_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(snarl) << endl;
                    return false;
                }
                cur_step = graph.get_next_step(cur_step);
            }
            // todo: we can compute flow at the same time
        }
        assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());

        gt_ref_path_name = ref_path_name;
        gt_ref_interval = make_pair(get<0>(ref_interval), get<1>(ref_interval));
        if (max_ploidy == -1) {
            max_ploidy = ploidy_at(ref_path_name, get<0>(ref_interval),
                                   ref_offsets.count(ref_path_name) ? ref_offsets.at(ref_path_name) : 0,
                                   ref_ploidies[ref_path_name]);
        }        
    } else {
        // if we have no reference infromation, try to get it from the parent snarl
        gt_ref_path_name = parent_ref_path_name;
        gt_ref_interval = parent_ref_interval;
        if (gt_ref_path_name.empty()) {
            // there's just no reference path through this snarl
            return false;
        }
        assert(max_ploidy >= 0);
    }
    
    // recurse on the children
    // todo: do we need to make this iterative for deep snarl trees?
    const vector<const Snarl*>& children = snarl_manager.children_of(&managed_snarl);

    for (const Snarl* child : children) {
        if (!snarl_manager.is_trivial(child, graph)) {
            bool called = call_snarl_recursive(*child, max_ploidy, gt_ref_path_name, gt_ref_interval, call_table);
            if (!called) {
                return false;
            }
        }
    }

#ifdef debug
    cerr << "recursively calling " << pb2json(managed_snarl) << " with " << children.size() << " children"
         << " and ref_path " << gt_ref_path_name << " and parent ref_path " << parent_ref_path_name << endl << endl;
#endif

    // abstract away the child snarls in the graph.  traversals will bypass them via
    // "virtual" edges
    SnarlGraph snarl_graph(&graph, snarl_manager, children);

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }
    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.shallow_contents(&snarl, graph, false);
        if (max(snarl_contents.first.size(), snarl_contents.second.size()) > max_snarl_shallow_size) {
            return false;
        }
        size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
        size_t length = 0;
        for (auto i = snarl_contents.first.begin(); i != snarl_contents.first.end() && length < len_threshold; ++i) {
            length += graph.get_length(graph.get_handle(*i));
        }
        greedy_avg_flow = length > len_threshold;
    }

    vector<SnarlTraversal> travs;
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle and overlay
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow, &snarl_graph);
        travs = std::move(weighted_travs.first);
           
    } else {
        // find the traversals using the generic interface
        assert(false);
        travs = traversal_finder.find_traversals(snarl);
    }

    // todo: we need to make reference traversal nesting aware
#ifdef debug
    for (int i = 0; i < travs.size(); ++i) {
        cerr << "[" << i << "]: " << pb2json(travs[i]) << endl;
    }
#endif
    
    // find the reference traversal in the list of results from the traversal finder
    if (!ref_path_name.empty()) {
        for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
            // todo: is there a way to speed this up?
            if (travs[i] == ref_trav) {
                ref_trav_idx = i;
            }
        }

        if (ref_trav_idx == -1) {
            ref_trav_idx = travs.size();
            // we didn't get the reference traversal from the finder, so we add it here
            travs.push_back(ref_trav);
#ifdef debug
            cerr << "[ref]: " << pb2json(ref_trav) << endl;
#endif
        }
    }
    // store the reference traversal information, which could be empty
    record.ref_path_name = ref_path_name;
    // The interval the ploidy above was decided over. emit_snarl_recursive re-derives the emit-time
    // ploidy from this via ploidy_at, and genotype_by_ploidy is indexed by that ploidy -- the field
    // was declared but never assigned, so every lookup ran at contig position 0 and a --ploidy-bed
    // whose region at base 0 differed could index genotype_by_ploidy out of bounds.
    record.ref_path_interval = make_pair((int64_t)gt_ref_interval.first,
                                         (int64_t)gt_ref_interval.second);
    record.ref_trav_idx = ref_trav_idx;

    // in the snarl graph, snarls a represented by a snarl end point and that's it.  here we fix up the traversals
    // to actually embed the snarls
    // todo: should be able to avoid copy here!
    vector<SnarlTraversal> embedded_travs = travs;
    for (int i = 0; i < embedded_travs.size(); ++i) {
        SnarlTraversal& traversal = embedded_travs[i];
        if (i != ref_trav_idx) {
            snarl_graph.embed_snarls(traversal);
        } else {
            snarl_graph.embed_ref_path_snarls(traversal);
        }
    }

    bool ret_val = true;

    if (traversals_only) {
        assert(gaf_output);
        for (SnarlTraversal& traversal : travs) {
            snarl_graph.embed_snarls(traversal);
        }
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, 0);
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else {
        // use our support caller to choose our genotype
        for (int ploidy = 1; ploidy <= max_ploidy; ++ploidy) {
            vector<int> trav_genotype;
            unique_ptr<SnarlCaller::CallInfo> trav_call_info;
            std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, gt_ref_path_name,  gt_ref_interval);
            assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

            // update the traversal finder with summary support statistics from this call
            // todo: be smarted about ploidy here
            NestedCachedPackedTraversalSupportFinder::SupportMap& child_support_map = nested_support_finder.child_support_map;
            // todo: re-use information that was produced in genotype!!
            int max_trav_size = 0;
            vector<Support> genotype_supports = nested_support_finder.get_traversal_genotype_support(embedded_travs, trav_genotype, {}, ref_trav_idx, &max_trav_size);
            Support total_site_support = std::accumulate(genotype_supports.begin(), genotype_supports.end(), Support());
            // todo: do we want to use max_trav_size, or something derived from the genotype? 
            child_support_map[snarl] = make_tuple(total_site_support, total_site_support, max_trav_size);

            // and now we need to update our own table with the genotype            
            if (record.genotype_by_ploidy.size() < ploidy) {
                record.genotype_by_ploidy.resize(ploidy);
            }
            record.genotype_by_ploidy[ploidy-1].first = trav_genotype;
            record.genotype_by_ploidy[ploidy-1].second.reset(trav_call_info.release());
            record.travs = embedded_travs;
        
            ret_val = trav_genotype.size() == ploidy;
        }
    }

    return ret_val;
}

bool NestedFlowCaller::emit_snarl_recursive(const Snarl& managed_snarl, int ploidy, CallTable& call_table) {
    // fetch the current snarl from the table
    CallRecord& record = call_table[managed_snarl];

    // only emit snarl with reference backbone:
    // todo: emit when no call (at least optionally)
    if (record.ref_trav_idx >= 0 && !record.genotype_by_ploidy.empty() && ploidy != 0) {

        if (ploidy < 0) {
            // Must agree with the ploidy the genotype was decided at, since genotype_by_ploidy is
            // indexed by it -- hence the record's own interval rather than the contig default.
            ploidy = ploidy_at(record.ref_path_name, record.ref_path_interval.first,
                               ref_offsets.count(record.ref_path_name)
                                   ? ref_offsets.at(record.ref_path_name) : 0,
                               ref_ploidies[record.ref_path_name]);
        }
        
        pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>>& genotype = record.genotype_by_ploidy[ploidy - 1];

        // compute count how many times a nested snarl appears in the genotype.  this will be the ploidy
        // it gets emitted with
        // todo: feed into flatten_alt_allele!
        map<Snarl, int, NestedCachedPackedTraversalSupportFinder::snarl_less> nested_ploidy;
        for (int allele : genotype.first) {
            const SnarlTraversal& allele_trav = record.travs[allele];
            for (size_t i = 0; i < allele_trav.visit_size(); ++i) {
                const Visit& visit = allele_trav.visit(i);
                if (visit.node_id() == 0) {
                    ++nested_ploidy[visit.snarl()];
                }
            }
        }

        // recurse on the children
        // todo: do we need to make this iterative for deep snarl trees? 
        const vector<const Snarl*>& children = snarl_manager.children_of(&managed_snarl);
        
        for (const Snarl* child : children) {
            if (!snarl_manager.is_trivial(child, graph)) {
                emit_snarl_recursive(*child, nested_ploidy[*child], call_table);
            }
        }

#ifdef debug
        cerr << "Recursively emitting " << pb2json(managed_snarl) << "with ploidy " << ploidy << endl;
#endif
        function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_flat_string =
            [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {

            string allele_string = trav_string(graph, travs[trav_allele]);
            if (trav_allele == ref_trav_idx) {
                return flatten_reference_allele(allele_string, call_table);
            } else {
                int allele_ploidy = std::max((int)std::count(travs_genotype.begin(), travs_genotype.end(), trav_allele), 1);
                return flatten_alt_allele(allele_string, std::min(allele_ploidy-1, genotype_allele), allele_ploidy, call_table);
            }
        };

        if (!gaf_output) {
            bool added = emit_variant(graph, snarl_caller, managed_snarl, record.travs, genotype.first, record.ref_trav_idx, genotype.second, record.ref_path_name,
                                 ref_offsets[record.ref_path_name], genotype_snarls, ploidy, trav_to_flat_string);
            if (!added) {
                return false;
            }
        } else {
            // todo:
            //    emit_gaf_variant(graph, snarl, travs, trav_genotype);
        }
    }

    return true;
}

string NestedFlowCaller::flatten_reference_allele(const string& nested_allele, const CallTable& call_table) const {

    string flat_allele;

    scan_snarl(nested_allele, [&](const string& fragment, Snarl& snarl) {
            if (!fragment.empty()) {
                flat_allele += fragment;
            } else {
                const CallRecord& record = call_table.at(snarl);
                assert(record.ref_trav_idx >= 0);
                if (record.travs.empty()) {
                    flat_allele += "<***>";
                    assert(false);
                } else{
                    const SnarlTraversal& traversal = record.travs[record.ref_trav_idx];
                    string nested_snarl_allele = trav_string(graph, traversal);
                    flat_allele += flatten_reference_allele(nested_snarl_allele, call_table);
                }
            }       
        });
    
    return flat_allele;
}

string NestedFlowCaller::flatten_alt_allele(const string& nested_allele, int allele, int ploidy, const CallTable& call_table) const {

    string flat_allele;
#ifdef debug
    cerr << "Flattening " << nested_allele << " at allele " << allele << endl;
#endif
    scan_snarl(nested_allele, [&](const string& fragment, Snarl& snarl) {
            if (!fragment.empty()) {
                flat_allele += fragment;
            } else {
                const CallRecord& record = call_table.at(snarl);
#ifdef debug
                cerr << "got record with " << record.travs.size() << " travs and " << record.genotype_by_ploidy.size() << " gts" << endl;
#endif
                int fallback_allele = -1;
                if (record.genotype_by_ploidy[ploidy-1].first.empty()) {
                    // there's no call here. but we really want to emit something, so try picking the reference
                    // or first allele
                    if (record.ref_trav_idx >= 0) {
                        fallback_allele = record.ref_trav_idx;
                    } else if (!record.travs.empty()) {
                        fallback_allele = 0;
                    }
                }
                if (fallback_allele >= (int)record.travs.size()) {
                    flat_allele += "<...>";
                } else {
                    // todo: passing in a single ploidy simplisitic, would need to derive from the calls when
                    // reucrising
                    // in practice, the results will nearly the same but still needs fixing
                    // we try to get the allele from the genotype if possible, but fallback on the fallback_allele
                    int trav_allele = fallback_allele >= 0 ? fallback_allele : record.genotype_by_ploidy[ploidy-1].first[allele];
                    const SnarlTraversal& traversal = record.travs[trav_allele];
                    string nested_snarl_allele = trav_string(graph, traversal);
                    flat_allele += flatten_alt_allele(nested_snarl_allele, allele, ploidy, call_table);
                }                
            }  
        });
    
    return flat_allele;
}



string NestedFlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}


SnarlGraph::SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls) :
    backing_graph(backing_graph),
    snarl_manager(snarl_manager) {
    for (const Snarl* snarl : snarls) {
        if (!snarl_manager.is_trivial(snarl, *backing_graph)) {
            this->snarls[backing_graph->get_handle(snarl->start().node_id(), snarl->start().backward())] =
                make_pair(backing_graph->get_handle(snarl->end().node_id(), snarl->end().backward()), true);
            this->snarls[backing_graph->get_handle(snarl->end().node_id(), !snarl->end().backward())] =
                make_pair(backing_graph->get_handle(snarl->start().node_id(), !snarl->start().backward()), false);
        }
    }
}

pair<bool, handle_t> SnarlGraph::node_to_snarl(handle_t handle) const {
    auto i = snarls.find(handle);
    if (i != snarls.end()) {
        return make_pair(true, i->second.first);
    } else {
        return make_pair(false, handle);
    }
}

tuple<bool, handle_t, edge_t> SnarlGraph::edge_to_snarl_edge(edge_t edge) const {
    auto i = snarls.find(edge.first);
    edge_t out_edge;
    handle_t out_node;
    bool out_found = false;
    if (i != snarls.end()) {
        // edge is from snarl start to after snarl end
        out_edge.first = i->second.first;
        out_edge.second = edge.second;
        out_node = edge.first;
        out_found = true;
    } else {
        // reverse of above
        i = snarls.find(backing_graph->flip(edge.second));
        if (i != snarls.end()) {
            out_edge.first = edge.first;
            out_edge.second = backing_graph->flip(i->second.first);
            out_node = edge.second;
            out_found = true;
        }
    }
    // note that we only have those two cases since our internal map contains
    // both orientations of the snarl.

    return make_tuple(out_found, out_node, out_edge);
}

void SnarlGraph::embed_snarl(Visit& visit) {
    handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
    auto it = snarls.find(handle);
    if (it != snarls.end()) {
        // edit the Visit in place to replace id, with the full snarl
        Snarl* snarl = visit.mutable_snarl();
        snarl->mutable_start()->set_node_id(visit.node_id());
        snarl->mutable_start()->set_backward(visit.backward());
        handle_t other = it->second.first;
        snarl->mutable_end()->set_node_id(backing_graph->get_id(other));
        snarl->mutable_end()->set_backward(backing_graph->get_is_reverse(other));
        if (it->second.second == false) {
            // put the snarl in an orientation consisten with other indexes
            swap(*snarl->mutable_start(), *snarl->mutable_end());
            snarl->mutable_start()->set_backward(!snarl->start().backward());
            snarl->mutable_end()->set_backward(!snarl->end().backward());
        }
        visit.set_node_id(0);
    }    
}

void SnarlGraph::embed_snarls(SnarlTraversal& traversal) {
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        if (visit.node_id() > 0) { 
            embed_snarl(visit);
        }
    }
}

void SnarlGraph::embed_ref_path_snarls(SnarlTraversal& traversal) {
    vector<Visit> out_trav;
    size_t snarl_count = 0;
    bool in_snarl = false;
    handle_t snarl_end;
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
        if (in_snarl) {
            // nothing to do if we're in a snarl except check for the end and come out
            if (handle == snarl_end) {
                in_snarl = false;
            } 
        } else {
            // if we're not in a snarl, check for a new one
            auto it = snarls.find(handle);
            if (it != snarls.end()) {
                embed_snarl(visit);
                snarl_end = it->second.first;
                in_snarl = true;
                ++snarl_count;
            }
            out_trav.push_back(visit);
        }
    }

    // switch in the updated traversal
    if (snarl_count > 0) {
        traversal.clear_visit();
        for (Visit& visit : out_trav) {
            *traversal.add_visit() = visit;
        }
    }
}

bool SnarlGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    if (!go_left) {        
        auto i = snarls.find(handle);
        if (i == snarls.end()) {
            return backing_graph->follow_edges(handle, go_left, iteratee);
        } else {
            return backing_graph->follow_edges(i->second.first, go_left, iteratee);
        }
    } else {
        return this->follow_edges_impl(backing_graph->flip(handle), !go_left, iteratee);
    }
}

// a lot of these don't strictly make sense.  ex, we would want has_node to
// hide stuff inside snarls.  but... we don't want to pay the cost of maintining
// structures for functions that aren't used..
bool SnarlGraph::has_node(nid_t node_id) const {
    return backing_graph->has_node(node_id);
}
handle_t SnarlGraph::get_handle(const nid_t& node_id, bool is_reverse) const {
    return backing_graph->get_handle(node_id, is_reverse);
}
nid_t SnarlGraph::get_id(const handle_t& handle) const {
    return backing_graph->get_id(handle);
}
bool SnarlGraph::get_is_reverse(const handle_t& handle) const {
    return backing_graph->get_is_reverse(handle);
}
handle_t SnarlGraph::flip(const handle_t& handle) const {
    return backing_graph->flip(handle);
}
size_t SnarlGraph::get_length(const handle_t& handle) const {
    return backing_graph->get_length(handle);
}
std::string SnarlGraph::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(handle);
}
size_t SnarlGraph::get_node_count() const {
    return backing_graph->get_node_count();
}
nid_t SnarlGraph::min_node_id() const {
    return backing_graph->min_node_id();
}
nid_t SnarlGraph::max_node_id() const {
    return backing_graph->max_node_id();
}
bool SnarlGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    return backing_graph->for_each_handle(iteratee, parallel);
}

}

