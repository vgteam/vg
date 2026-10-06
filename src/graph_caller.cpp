#include <atomic>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

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

// The names of the AtomizeCounters::refuse reasons, in index order. The initializer sets the
// size, so that the check below fails when a name is missing as well as when one is extra.
static const char* const g_atomize_refuse_name[] = {
    "the genotyper returned no genotype: ploidy 0, or no read the matrix could place",
    "no reference traversal",
    "the snarl does not resolve",
    "the reference projection is empty",
    "the alignment's reference-to-alt step map has the wrong length",
    "no difference blocks: every called haplotype takes the reference route here",
    "an indel needs an anchor base and the snarl has none to its left",
    "the visit left of the anchor has no sequence",
    "every block spelled the reference's own bases: a route difference with no sequence difference",
    "one block, saying what the site record already says",
    "-L merged the called alleles, which blocks would spell apart again",
    "one block that two strands' routes spell differently within one site allele, so the site cannot give it a genotype",
    "one block, where a chain crossed more than once may leave its record short of what the site record says",
};
static_assert(sizeof(g_atomize_refuse_name) / sizeof(g_atomize_refuse_name[0])
                  == sizeof(AtomizeCounters::refuse) / sizeof(AtomizeCounters::refuse[0]),
              "each AtomizeCounters::refuse reason must have a name in g_atomize_refuse_name, "
              "and each name must have a reason");
static thread_local int g_descent_depth = 0;

void FlowCaller::report_descent_instrumentation() const {
    size_t total = 0;
    for (int d = 0; d < 16; ++d) {
        total += descent_counters.depth_hist[d].load();
    }
    if (total == 0) {
        return;   // no symbolic descent in this run
    }
    cerr << "[vg call] descent depth:";
    for (int d = 1; d < 16; ++d) {
        size_t n = descent_counters.depth_hist[d].load();
        if (n > 0) {
            cerr << " " << d << "=" << n;
        }
    }
    cerr << " (" << total << " child calls)" << endl;
    if (descent_counters.child_multi_crossing.load() > 0) {
        cerr << "[vg call] descent: " << descent_counters.child_multi_crossing.load()
             << " children a called traversal enters more than once; visits after the first are"
             << " masked, so each contributes one copy and its first crossing's distance" << endl;
    }
    cerr << "[vg call] descent skipped: " << descent_counters.skipped_no_copy.load()
         << " children no called allele reaches, " << descent_counters.skipped_no_ref.load()
         << " with no reference path through them" << endl;
    if (descent_counters.off_reference.load() > 0 || descent_counters.no_ref_recorded.load() > 0) {
        cerr << "[vg call] off-reference nested: " << descent_counters.off_reference.load()
             << " chains the reference does not cross were descended into, "
             << descent_counters.no_ref_recorded.load() << " recorded into the linkage layer with no line;"
             << " copies 0/1/2 = " << descent_counters.no_ref_copies[0].load() << "/"
             << descent_counters.no_ref_copies[1].load() << "/" << descent_counters.no_ref_copies[2].load() << endl;
    }
}

void VCFOutputCaller::report_atomize_instrumentation() const {
    size_t unresolvable = atomize_counters.site_unresolvable.load();
    // Keyed on whether block emission ran at all, not on any refusal counter, so that the line
    // below is written whenever block emission ran.
    if (atomize_counters.sites.load() == 0) {
        return;
    }


    // `unresolvable` should be zero on the ordinary path, where every site is a managed snarl and a
    // reversed one resolves through its reversed boundaries. Under -I/--chains it need not be: a
    // chain piece is a constructed snarl the manager does not know. The second number counts sites
    // that resolved only through their reversed boundaries.
    cerr << "[vg call] atomize: " << unresolvable
         << " sites where projection is inert because the snarl does not resolve, "
         << atomize_counters.site_reversed.load()
         << " resolved as the reversal flip_snarl produces" << endl;


    if (atomize_counters.child_inlined.load() > 0) {
        cerr << "[vg call] atomize: " << atomize_counters.child_inlined.load()
             << " child chains left without a line because a block ALT already spells them" << endl;
    }
    {
        // One line, listing only the reasons that occurred.
        const size_t reasons = sizeof(AtomizeCounters::refuse) / sizeof(AtomizeCounters::refuse[0]);
        size_t total = 0;
        for (size_t i = 0; i < reasons; ++i) {
            total += atomize_counters.refuse[i].load();
        }
        if (total > 0) {
            cerr << "[vg call] atomize: " << total << " sites declined block emission, so the site"
                 << " record stands:";
            bool first = true;
            for (size_t i = 0; i < reasons; ++i) {
                size_t n = atomize_counters.refuse[i].load();
                if (n > 0) {
                    cerr << (first ? " " : "; ") << n << " " << g_atomize_refuse_name[i];
                    first = false;
                }
            }
            cerr << endl;
        }
    }
    if (atomize_counters.split_sites.load() > 0) {
        cerr << "[vg call] atomize: " << atomize_counters.split_sites.load()
             << " sites written as their difference blocks rather than their site record, "
             << atomize_counters.split_lines.load()
             << " lines" << endl;
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

void GraphCaller::set_snarl_batching(size_t window_size) {
    snarl_batch_window = window_size;
}

/// Key a snarl by the lower of its two boundary node IDs, so that sorting on it puts snarls with
/// nearby node IDs next to each other.
static nid_t snarl_node_key(const Snarl* snarl) {
    return min(snarl->start().node_id(), snarl->end().node_id());
}

/// Sort snarls by `snarl_node_key`, so that snarls with nearby node IDs are called one after
/// another.
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
    if (snarl_batch_window > 0) {
        // Call in node-ID order, batched by window, so that snarls with nearby node IDs are
        // called together.
        vector<const Snarl*> roots;
        snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            roots.push_back(snarl);
        });
        sort_snarls_by_node_id(roots);

        // Split the sorted snarls into runs that share a window, with one parallel job per run
        // rather than per snarl.
        vector<pair<size_t, size_t>> windows;
        size_t begin = 0;
        while (begin < roots.size()) {
            size_t window = (size_t)(snarl_node_key(roots[begin]) / (nid_t)snarl_batch_window);
            size_t end = begin + 1;
            while (end < roots.size() &&
                   (size_t)(snarl_node_key(roots[end]) / (nid_t)snarl_batch_window) == window) {
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

        if (snarl_batch_window > 0) {
            // Keep queued children in node-ID order too.
            sort_snarls_by_node_id(cur_queue);
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
    if (show_progress && nested_snarl_count > 0) cerr << "[vg call]: Finished processing " << nested_snarl_count << " nested snarls" << endl;
  
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
    suppressed_ref_info.resize(get_thread_count());
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
        ss << nesting_info_headers();
    }
    if (emit_phasing) {
        // FORMAT/PS is the VCF phase set, which phasing tools read. It is unrelated to INFO/PS
        // above, vg's parent-snarl field; the two are in different namespaces, so both are legal,
        // and their descriptions say which is which.
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
           << "difference block within its snarl. A snarl is written as one record per difference "
           << "block where the reference and the called haplotypes differ from each other in more "
           << "than one place inside it, or where its own record would repeat a child snarl's, so "
           << "the count can be 1. A block record's ID is the snarl's ID with _ and the index "
           << "appended. DOUBLE COUNTING: the per-sample evidence is the SNARL's, repeated on every "
           << "block, not apportioned between them -- AD, GL, GQ, GQI, GP and QUAL are identical "
           << "across the set, because the genotype likelihood was computed over whole-snarl "
           << "traversals and has no per-block decomposition. DP, DR and BL are per-site read "
           << "counts and are site-level by definition. So any consumer that sums, averages or "
           << "otherwise aggregates evidence across records must group by the snarl's ID first and "
           << "count each snarl once. Records without SB are unaffected: they are the only record "
           << "their snarl emitted.\">" << endl;
    }
    if (allele_merge_threshold < 1.0) {
        ss << "##INFO=<ID=MAT,Number=.,Type=String,Description=\"Merged Allele Traversal: "
           << "ALT alleles merged after genotyping by -L/--cluster, as OLD>NEW:SIMILARITY using "
           << "pre-merge allele numbers. AD and GL are folded onto the surviving allele and MAD is "
           << "recomputed; DP, QUAL, GQ, GP and FILTER are as computed over the pre-merge allele set. "
           << "In a nested run this record gives the collapsed view of the site and its child "
           << "records the precise one, so they disagree by design.\">"
           << endl;
    }
    return ss.str();
}

void VCFOutputCaller::set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                                  const vector<size_t>* sequence_to_haplotype) {
    this->linkage_collector = collector;
    this->linkage_gbwt = gbwt;
    this->linkage_sequence_to_haplotype = sequence_to_haplotype;
    this->linkage_panel_size = collector != nullptr ? collector->panel_size() : 0;
    this->linkage_gbwt_cache.clear();
    this->linkage_gbwt_cache_origin.clear();
    if (gbwt != nullptr) {
        // One per thread, built here so the parallel region never allocates one.
        this->linkage_gbwt_cache.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            this->linkage_gbwt_cache.emplace_back(*gbwt);
        }
        this->linkage_gbwt_cache_origin.assign(omp_get_max_threads(), 0);
    }
}

vector<int> VCFOutputCaller::panel_alleles(const HandleGraph& graph,
                                          const vector<SnarlTraversal>& travs) const {
    vector<int> out;
    if (linkage_gbwt == nullptr || linkage_sequence_to_haplotype == nullptr) {
        return out;
    }
    // -1 means the haplotype carries no allele here, which is different from carrying the
    // reference: a haplotype whose path ends inside the site has nothing to say. Sized by the
    // panel, since the row is indexed by haplotype.
    const size_t row = linkage_panel_size > 0 ? linkage_panel_size
                                              : linkage_sequence_to_haplotype->size();
    out.assign(row, -1);

    // The cache, not the index: same results, but records stay decompressed between sites.
    // Falls back to the index itself if set_linkage was never given one to size the vector.
    int thread = omp_get_thread_num();
    const bool cached = (size_t)thread < linkage_gbwt_cache.size();

    // CachedGBWT only grows, and with node-ID-ordered windows a thread does not come back to an
    // earlier window, so the cache is cleared when the site moves more than a fetch window past
    // where it was filled. Adjacent snarls still share records, and the cache stays to about one
    // window.
    if (cached && (size_t)thread < linkage_gbwt_cache_origin.size() && !travs.empty()) {
        static const nid_t CACHE_ANCHOR_SPAN = 4096;
        nid_t lead = 0;
        for (int64_t i = 0; i < travs[0].visit_size() && lead == 0; ++i) {
            lead = travs[0].visit(i).node_id();
        }
        if (lead != 0) {
            nid_t& anchor = linkage_gbwt_cache_origin[thread];
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
                    // A haplotype stored as several fragments could reach one site twice, with
                    // two traversals; the last one written wins.
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
        // The callers support ploidy 1 and 2 only, so reject anything else here.
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

bool VCFOutputCaller::buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b) {
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
        resolve_linkage_level(0, true);
        return;
    }
    // Resolve every level, since chain construction skips entries of later levels
    // than the one being resolved. `max_level()` is read again on each pass, since a pass can
    // add a chain at a deeper level.
    for (size_t gen = 0;; ++gen) {
        const size_t deepest = linkage_collector->max_level();
        resolve_linkage_level(gen, gen >= deepest);
        if (gen >= deepest) {
            break;
        }
    }
}

/// The ID of the site a record belongs to: a block record's ID without the "_<index>" that
/// tells the site's block records apart, and any other record's ID unchanged.
static string block_site_name(const string& id) {
    size_t underscore = id.rfind('_');
    return underscore == string::npos ? id : id.substr(0, underscore);
}

size_t VCFOutputCaller::record_key_of(const Snarl& snarl) const {
    return std::hash<string>{}(print_snarl(snarl, false));
}

// Each read's strand log-odds for the render, used by the anchors. Built here rather than taken
// from re-genotyping, which may not have run and whose table is built before `phase_sites` is
// final.
size_t VCFOutputCaller::phase_set_id(const string& contig, size_t phase_set) {
    return phase_set_ids.emplace(make_pair(contig, phase_set), phase_set_ids.size()).first->second;
}

void VCFOutputCaller::build_render_lambda() {
    render_lambda.clear();
    render_lambda_site.clear();
    render_lambda_phase_set.clear();
    render_lambda_temper = 0.0;
    render_lambda_ceiling = 1.0;
    if (phase_sites.empty()) {
        return;
    }
    RegenotypeCounters scratch;
    accumulate_lambda(phase_sites, phase_flips, render_lambda, scratch);
    for (const PhaseSite& site : phase_sites) {
        render_lambda_site[site.record_key] = &site;
    }
    // The last PhaseCall written winning, as in `build_render_phases`.
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        render_lambda_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
    }
    // The summed strand log-odds overstate how sure the strand is, so they are tempered. Use the
    // temper re-genotyping fitted, where it ran; otherwise fit one here.
    if (regenotype_counters.fitted_temper > 0.0) {
        render_lambda_temper = regenotype_counters.fitted_temper;
        render_lambda_ceiling = regenotype_counters.fitted_ceiling;
    } else {
        double temper = -1.0;
        double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
        RegenotypeCounters fit_scratch;
        fit_calibration(phase_sites, phase_flips, render_lambda, regenotype_params, temper, ceiling,
                        fit_scratch);
        if (fit_scratch.fitted_temper > 0.0) {
            render_lambda_temper = fit_scratch.fitted_temper;
            render_lambda_ceiling = fit_scratch.fitted_ceiling;
        }
    }
}

double VCFOutputCaller::read_strand_log_odds(size_t record_key, std::string_view read_name) const {
    if (render_lambda.empty() || render_lambda_temper <= 0.0) {
        return 0.0;
    }
    const uint64_t key = (uint64_t)std::hash<std::string_view>{}(read_name);
    const auto found = render_lambda.find(key);
    const auto ps = render_lambda_phase_set.find(record_key);
    const size_t phase_set = ps != render_lambda_phase_set.end() ? ps->second : NO_PHASE_SET;
    if (found == render_lambda.end()) {
        // The read reached no phased site.
        return 0.0;
    }
    if (!read_strand_usable(found->second, phase_set)) {
        // The read has a strand, but in another phase set, whose strands do not correspond to
        // this site's. NaN rather than 0, because a split homozygous site drops such a read but
        // places one with no strand by a coin (see `build_site_anchors`).
        return std::numeric_limits<double>::quiet_NaN();
    }
    double value = found->second.lambda;
    size_t sites = found->second.sites;
    // Subtract this record's own contribution, so that a site is not judged by its own evidence;
    // if it was the only one, there is nothing left.
    const auto site = render_lambda_site.find(record_key);
    if (site != render_lambda_site.end()) {
        unordered_map<uint64_t, double> own;
        site_own_log_odds(*site->second, phase_flips.count(record_key) != 0, own);
        const auto mine = own.find(key);
        if (mine != own.end()) {
            value -= mine->second;
            if (sites > 0) {
                --sites;
            }
        }
    }
    if (sites == 0) {
        return 0.0;
    }
    return calibrated_log_odds(value, render_lambda_temper, render_lambda_ceiling);
}

void VCFOutputCaller::build_render_phases() {
    // Built from the phasing the linkage pass accumulated, as read phasing left it. Sites with no line
    // are included; they are simply never looked up.
    render_phases.clear();
    if (!emit_phasing) {
        return;
    }
    render_phases.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        // Where a site has more than one PhaseCall, the last one written wins.
        render_phases[pc.record_key] = pc;
    }
}

void VCFOutputCaller::finalise_linkage_outputs() {
    // Built after every record has been rendered, since the mosaic needs to know which sites have
    // a line, which is not known while genotypes are being resolved.
    if (linkage_collector == nullptr) {
        return;
    }
    // Read from the collector, since each PhaseCall's `emitted` was copied before any line was
    // written.
    const std::unordered_set<size_t> emitted_records = linkage_collector->emitted_records();
    size_t unexplained = 0;
    size_t order_arbitrary = 0;
    // Count the phased sites, separating those that became records from those that did not.
    size_t phased_unwritten = 0;
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        if (emitted_records.count(pc.record_key) == 0) {
            // Phased, since its children take their strand from it, but not a record, so it is kept
            // out of the mosaic and the record counts.
            ++phased_unwritten;
            continue;
        }
        // Count only the strands a site has. A haploid site has one strand and a wildcard, and the
        // wildcard can be in either slot: a haploid contig fills the first slot, while a nested
        // site on its parent's second strand fills the second.
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
    if (linkage_collector->num_duplicate_live_keys() > 0) {
        // Duplicate keys need not change the output, but `retract` cannot handle those sites, since
        // it retracts only the first live entry.
        cerr << "[vg call] linkage: " << linkage_collector->num_duplicate_live_keys()
             << " sites recorded onto a key that already had a live entry; the retract path cannot"
             << " address these" << endl;
    }
    if (linkage_collector->model_params().hp_prior > 0.0) {
        cerr << "[vg call] linkage: " << linkage_collector->num_site_prior_entries()
             << " live entries decoded at a run-length site's own frequency exponent (--hp-prior)"
             << endl;
    }
    if (emit_phasing) {
        // At sites where a strand is on the wildcard, no panel haplotype names it, so the phase
        // across them rests on the transitions alone.
        cerr << "[vg call] phasing: " << (linkage_phased.size() - phased_unwritten)
             << " sites phased, " << unexplained
             << " with a strand the panel does not explain" << endl;
        if (phased_unwritten > 0) {
            // Sites that wrote no VCF line but are phased. A parent whose alleles differ only inside
            // its children is written as the reference and has no line, and its children still need
            // to know which of its strands carries the chain.
            cerr << "[vg call] phasing: " << phased_unwritten
                 << " collapsed sites phased with no line of their own, so their children can"
                 << " inherit a strand" << endl;
        }
        if (order_arbitrary > 0) {
            // Heterozygous sites where no panel haplotype on either strand carries either called
            // allele. The record is still phased and in the phase set, but its order came from
            // sorting the pair, so it is arbitrary.
            cerr << "[vg call] phasing: " << order_arbitrary
                 << " heterozygous sites carry an allele order the panel does not determine"
                 << endl;
        }
    }
    if (!mosaic_path.empty()) {
        // Records only: the mosaic's segments are runs over sites of the call set, and it accounts
        // for exactly the written records.
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

void VCFOutputCaller::resolve_linkage_level(size_t level, bool last) {
    linkage_resolved = true;
    if (linkage_collector == nullptr) {
        return;
    }
    // Time the pass and report the collector's size.
    auto start = std::chrono::steady_clock::now();
    // `linkage_phased` accumulates across levels, since the model needs the earlier ones: a
    // nested site's strand is read from its parent's PhaseCall, and a clamped site's phase is
    // pinned to its chosen pair.
    const size_t moved =
        linkage_collector->resolve_level(level, last,
                                              emit_phasing ? &linkage_phased : nullptr);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    linkage_seconds += seconds;
    // How many sites the model moved off the genotype the reads alone chose.
    linkage_changed += moved;
    if (!last) {
        // One line per level except the last: its site count, how many of its genotypes the
        // linkage model moved, and the seconds it took.
        cerr << "[vg call] linkage level " << level << ": "
             << linkage_collector->num_sites_at(level) << " sites, "
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
    // `buf` must not be const, since std::move() over const iterators copies, and the whole VCF is
    // in memory here. Each buffer is freed as it is moved. This makes write_variants() usable only
    // once.
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
    // Resolve the linkage model, if it has not been resolved, before the records are written.
    resolve_linkage();
    finalise_linkage_outputs();


    for (const auto& v : all_variants) {
        string dest;
        int ret = zstdutil::DecompressString(v.second, dest);
        assert(ret == 0);
        // The record key is the hash of the site's ID, as `record_key_of` computes it, so the line
        // itself gives the site's identity; a block record's ID carries it before its suffix.
        // Computed once, when first needed; several records can share a (contig, position), and
        // each must get its own site's values.
        size_t line_key = 0;
        bool have_line_key = false;
        auto id_key = [&]() -> size_t {
            if (!have_line_key) {
                size_t a = dest.find('\t');
                size_t b = a == string::npos ? string::npos : dest.find('\t', a + 1);
                size_t c = b == string::npos ? string::npos : dest.find('\t', b + 1);
                if (c != string::npos) {
                    line_key = std::hash<string>{}(block_site_name(dest.substr(b + 1, c - b - 1)));
                }
                have_line_key = true;
            }
            return line_key;
        };
        if (linkage_collector != nullptr) {
            // Quality first, then phasing. The line already carries the chosen genotype, since it
            // was built from it.
            const auto& quality = linkage_collector->moved_quality();
            if (!quality.empty()) {
                auto found = quality.find(id_key());
                if (found != quality.end()) {
                    if (!ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype(
                            dest, found->second, linkage_min_confidence)) {
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
    // Reported after the records are rendered, since block emission happens as they are.
    report_atomize_instrumentation();
}


gbwt::edge_type VCFOutputCaller::mosaic_position_at(gbwt::node_type node, size_t hap) const {
    if (linkage_gbwt == nullptr || linkage_sequence_to_haplotype == nullptr) {
        return gbwt::invalid_edge();
    }
    // Finding a position costs a `locate` for each sequence in the node's range, far more than an
    // LF step, and the same (node, haplotype) is asked for repeatedly, so positions are cached.
    const uint64_t key = ((uint64_t)node << 20) | (uint64_t)(hap & 0xFFFFF);
    auto hit = mosaic_position_cache.find(key);
    if (hit != mosaic_position_cache.end()) {
        return hit->second;
    }
    gbwt::SearchState state = linkage_gbwt->find(node);
    if (!state.empty()) {
        for (gbwt::size_type i = state.range.first; i <= state.range.second; ++i) {
            gbwt::size_type seq = linkage_gbwt->locate(node, i);
            if (seq < linkage_sequence_to_haplotype->size()
                && (*linkage_sequence_to_haplotype)[seq] == hap) {
                mosaic_position_cache[key] = gbwt::edge_type(node, i);
                return gbwt::edge_type(node, i);
            }
        }
    }
    mosaic_position_cache[key] = gbwt::invalid_edge();
    return gbwt::invalid_edge();
}

/// Follow a walk whose direction is already known, rather than guessing the direction. Given the
/// oriented node, `mosaic_position_at` finds the position, and `LF` continues in the same
/// direction. A local guess at the direction, such as the node's forward orientation, assumes the
/// walk advances in reference order, which fails where the sample's walk does not follow the
/// reference, as at large balanced structural variants.
bool VCFOutputCaller::mosaic_follow(gbwt::edge_type start, int64_t to_node,
                                    gbwt::node_type* out_end) const {
    if (linkage_gbwt == nullptr || start == gbwt::invalid_edge()) {
        return false;
    }
    // Finding a position is far more costly than an LF step, so the caller passes the position in,
    // and the walk itself is cheap.
    if ((int64_t)gbwt::Node::id(start.first) == to_node) {
        if (out_end != nullptr) *out_end = start.first;
        return true;
    }
    gbwt::edge_type at = start;
    for (size_t step = 0; step < MOSAIC_WALK_LIMIT; ++step) {
        at = linkage_gbwt->LF(at);
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

gbwt::edge_type VCFOutputCaller::mosaic_gbwt_position(int64_t node_id, size_t hap) const {
    // Forward first: snarl boundaries are stored oriented along the reference, so the reverse
    // orientation is the exception.
    //
    // The answer is ambiguous, since a GBWT stores each path in both orientations. It serves only
    // callers with no direction to work from, which ask whether the haplotype is at the node at
    // all; a position on a walk comes from `mosaic_follow`. Uses `mosaic_position_at`, and so its
    // cache.
    for (int orientation = 0; orientation < 2; ++orientation) {
        const gbwt::edge_type at =
            mosaic_position_at(gbwt::Node::encode(node_id, orientation == 1), hap);
        if (at != gbwt::invalid_edge()) {
            return at;
        }
    }
    return gbwt::invalid_edge();
}

vector<int> VCFOutputCaller::phase_ordered_genotype(size_t record_key,
                                                    const vector<int>& genotype) const {
    vector<int> ordered = genotype;
    if (!emit_phasing || ordered.size() != 2) {
        return ordered;
    }
    const auto found = render_phases.find(record_key);
    // Only on an exact reversal. A PhaseCall that is not a permutation of the chosen pair is left
    // alone, as `emit_variant` refuses to apply one. For a homozygote the swap changes nothing.
    if (found != render_phases.end() && found->second.ploidy == 2
        && found->second.trav_first == ordered[1]
        && found->second.trav_second == ordered[0]) {
        std::swap(ordered[0], ordered[1]);
    }
    return ordered;
}

int VCFOutputCaller::phase_haploid_slot(size_t record_key, const vector<int>& genotype) const {
    if (!emit_phasing || genotype.size() != 1) {
        return 0;
    }
    const auto found = render_phases.find(record_key);
    if (found == render_phases.end() || found->second.ploidy != 1
        || found->second.nested_strand < 0) {
        // No nested strand means a haploid locus, such as chrY or a haploid --ploidy-bed region,
        // where slot 1 means nothing.
        return 0;
    }
    // Only where the phase names the allele chosen for this site, as `phase_ordered_genotype` and
    // `emit_variant` require.
    if (found->second.trav_first != genotype[0]) {
        return 0;
    }
    return (int)found->second.nested_strand;
}

// The anchor gqn column for a staged site.
//
// `gq_fraction` was computed in the direct pass for the reads' best genotype, so on a record whose
// genotype the linkage model changed, it describes the abandoned genotype. Such a record gets the
// signed margin of its chosen genotype instead, as the VCF's GQN does (see
// ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
//
// Returns the direct pass's value when the model did not change the call; the recomputed signed margin
// when it did; and NaN, written as ".", when it did but the margin cannot be recomputed.
double FlowCaller::anchor_gqn_for(const PendingRecord& rec,
                                  const vector<int>& chosen) const {
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    // The direct pass's value, with its "no gap to normalise" value (-1) turned into NaN, written as ".",
    // so that it stays distinct from the signed range [-1, 1].
    const double direct_value = (info == nullptr || info->gq_fraction < 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : info->gq_fraction;
    const double blank = std::numeric_limits<double>::quiet_NaN();
    if (linkage_collector == nullptr) {
        return direct_value;
    }
    const auto& moved = linkage_collector->moved_quality();
    const auto found = moved.find(rec.record_key);
    if (found == moved.end()) {
        return direct_value;   // linkage left the call alone, so the direct pass's value still holds
    }
    // The model changed the call, so the direct pass's value describes the wrong genotype, and any
    // failure below gives NaN rather than falling back to it.
    if (info == nullptr || info->genotype_lls.empty()) {
        return blank;
    }
    // The divisor and share the VCF's GQN uses (see
    // ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype), so that the two agree.
    const LinkageCollector::DirectQuality& direct = found->second.direct;
    if (!(direct.achievable_gap > 0.0)) {
        return blank;   // no scale, and no honest pre-linkage value to fall back on
    }
    const double achievable_phred = 10.0 * direct.achievable_gap / log(10.0);

    // The chosen genotype, not rec.genotype, which is the direct pass's call before the linkage model:
    // the reads prefer that call, so its margin would have the wrong sign.
    vector<int> called = chosen;
    sort(called.begin(), called.end());
    const auto mine = info->genotype_lls.find(called);
    if (mine == info->genotype_lls.end()) {
        return blank;
    }
    // Only genotypes over the written alleles, as the VCF's GL has: the reference traversal and
    // the ones the chosen genotype names. A traversal with no ALT could otherwise beat the call.
    set<int> emitted(called.begin(), called.end());
    if (rec.ref_trav_idx >= 0) {
        emitted.insert(rec.ref_trav_idx);
    }
    double best_other = -numeric_limits<double>::infinity();
    for (const auto& entry : info->genotype_lls) {
        if (entry.first == called) {
            continue;
        }
        bool all_emitted = true;
        for (int a : entry.first) {
            if (emitted.count(a) == 0) {
                all_emitted = false;
                break;
            }
        }
        if (all_emitted) {
            best_other = max(best_other, entry.second);
        }
    }
    if (!std::isfinite(best_other)) {
        return blank;
    }
    // Nats to phred, matching the VCF's GL, which is log10.
    const double margin_phred = 10.0 * (mine->second - best_other) / log(10.0);
    return min(1.0, max(-1.0, margin_phred / achievable_phred * direct.explained_share));
}

/// The genotype the linkage model chose, or the direct pass's own where it chose none. Used by both
/// anchor-collection paths, the render and `hand_off_deferred_records`, so that records with no
/// VCF line (`reported_inline` and `no_reference`) also get anchors for their chosen genotype.
vector<int> FlowCaller::chosen_genotype_for(const PendingRecord& rec) const {
    vector<int> genotype = rec.genotype;
    int chosen_a = -1, chosen_b = -1;
    size_t chosen_ploidy = 0;
    if (linkage_collector != nullptr
        && linkage_collector->chosen_traversals(rec.record_key, &chosen_a, &chosen_b,
                                                 &chosen_ploidy)
        && chosen_ploidy == genotype.size()) {
        genotype.assign(1, chosen_a);
        if (chosen_ploidy > 1) {
            genotype.push_back(chosen_b);
        }
    }
    return genotype;
}

void FlowCaller::collect_anchors_for_record(const PendingRecord& rec,
                                            const vector<int>& genotype) {
    collect_anchors_for(rec.snarl, phase_ordered_genotype(rec.record_key, genotype),
                        phase_haploid_slot(rec.record_key, genotype), rec.call_info,
                        anchors_want_leaf_test() ? snarl_is_leaf(rec.snarl) : true,
                        anchor_gqn_for(rec, genotype), rec.record_key);
}

void VCFOutputCaller::collect_anchors_for(const Snarl& snarl, const vector<int>& genotype,
                                          int haploid_slot,
                                          const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                          bool is_leaf, double gqn, size_t record_key) {
    if (anchor_path.empty() || anchor_writer == nullptr || call_info == nullptr) {
        return;
    }
    if (anchor_params.leaf_only && !is_leaf) {
        return;
    }
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
    if (info == nullptr || info->anchor_evidence == nullptr) {
        // A genotype derived from a parent rather than scored here, or a run whose caller is not the
        // read-likelihood one. There are no per-read responsibilities to partition on.
        return;
    }
    vector<AnchorWriter::Anchor> anchors;
    // Each read's strand log-odds, leaving out this record.
    vector<double> read_strand;
    // Built only where `build_site_anchors` reads it: at a diploid homozygote that may be split, or
    // at a heterozygous site under --anchors-phase-hets or --anchors-strict-hets. The test must
    // match its gate, which also checks the vector's length against `evidence.reads`.
    const bool splittable_hom = genotype.size() == 2 && genotype[0] == genotype[1];
    const bool tiltable_het = (anchor_params.phase_hets || anchor_params.strict_hets)
                              && genotype.size() == 2
                              && genotype[0] != genotype[1];
    if ((anchor_params.hom_split && splittable_hom) || tiltable_het) {
        read_strand.reserve(info->anchor_evidence->reads.size());
        for (const AnchorRead& read : info->anchor_evidence->reads) {
            read_strand.push_back(read_strand_log_odds(record_key, read_names().name(read.read)));
        }
    }
    build_site_anchors(*info->anchor_evidence, genotype, print_snarl(snarl),
                       gqn,
                       info->explained_share, haploid_slot, anchor_params, *anchor_params.counters,
                       anchors,
                       (anchor_params.hom_split || anchor_params.phase_hets
                        || anchor_params.strict_hets) ? &read_strand
                                                                             : nullptr);
    // A check for --anchors-hom-split, reported per run: at heterozygous sites, whose alleles show
    // which strand each read is on, how often the read's strand log-odds agree. The log-odds leave
    // the site out. Computed only when splitting is on, and not when the heterozygous placement
    // itself uses the strand log-odds, since the check would then compare the strand with
    // itself.
    if (anchor_params.hom_split && !anchor_params.phase_hets && !anchor_params.strict_hets
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
                    const double lo = read_strand_log_odds(record_key, read_names().name(row.read));
                    if (std::isnan(lo) || lo == 0.0) {
                        anchor_params.counters->phase_no_opinion.fetch_add(1);
                        continue;
                    }
                    const int phase_slot = lo > 0.0 ? 0 : 1;
                    const bool agree = phase_slot == anchor.slot;
                    anchor_params.counters->phase_checked.fetch_add(1);
                    if (agree) {
                        anchor_params.counters->phase_agree.fetch_add(1);
                    }
                    // The same threshold the split uses, --split-min-q, in natural-log units.
                    if (std::abs(lo) >= anchor_params.phase_min) {
                        anchor_params.counters->phase_confident.fetch_add(1);
                        if (agree) {
                            anchor_params.counters->phase_confident_agree.fetch_add(1);
                        }
                    }
                }
            }
        }
    }
    for (AnchorWriter::Anchor& anchor : anchors) {
        anchor_writer->add(std::move(anchor));
    }
}

void VCFOutputCaller::write_anchors() {
    if (anchor_path.empty() || anchor_writer == nullptr) {
        return;
    }
    size_t anchors = anchor_writer->anchor_count();
    size_t rows = anchor_writer->read_row_count();
    if (!anchor_writer->write(anchor_path, anchor_graph_name, sample_name, anchor_reads_source,
                              anchor_mismap_min, anchor_params)) {
        return;
    }
    cerr << "[vg call] anchors: " << anchors << " written over " << rows
         << " read placements to " << anchor_path << endl;
    if (anchor_params.counters != nullptr) {
        anchor_params.counters->report(cerr);
    }
}

void VCFOutputCaller::write_mosaic(const vector<LinkageCollector::PhaseCall>& phasing) const {
    ofstream out(mosaic_path);
    if (!out) {
        cerr << "error [vg call]: could not open " << mosaic_path << " for the mosaic output"
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
    out << "#graph\t" << mosaic_graph_name << "\n";
    out << "#sample\t" << sample_name << "\n";
    // gRef fragments are counted, not listed: a cover can name thousands of contigs, and each row
    // names its own.
    size_t gref_fragments = 0;
    for (const string& ref : mosaic_reference_paths) {
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
    out << "#patch\t" << (mosaic_patch_gaps ? "reference" : "none") << "\n";
    out << "#nested\t" << (mosaic_keep_nested ? "kept" : "merged") << "\n";
    out << "#unexplained\t" << (mosaic_connect_unexplained ? "connected" : "broken") << "\n";
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
    for (size_t h = 0; h < mosaic_haplotype_names.size(); ++h) {
        out << "#haplotype\t" << h << "\t" << mosaic_haplotype_names[h] << "\n";
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
    // name, since the index follows GBWT metadata order. `mosaic_reference_paths` holds full path
    // names (CHM13#0#chr20) and the panel names haplotypes as sample#phase (CHM13#0), so the
    // contig is dropped before matching.
    size_t reference_hap = LinkageModel::WILDCARD;
    for (const string& full : mosaic_reference_paths) {
        // Never a gRef path: a gRef cover is stitched together from many donors, and it is not in the
        // panel.
        if (GrefCover::is_gref_derived(full)) {
            continue;
        }
        size_t h1 = full.find('#');
        size_t h2 = h1 == string::npos ? string::npos : full.find('#', h1 + 1);
        const string base = h2 == string::npos ? full : full.substr(0, h2);
        for (size_t k = 0; k < mosaic_haplotype_names.size(); ++k) {
            if (mosaic_haplotype_names[k] == base) {
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

    const bool patch_gaps = mosaic_patch_gaps;
    const bool keep_nested = mosaic_keep_nested;
    const bool connect_unexplained = mosaic_connect_unexplained;
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
                                  : mosaic_gbwt_position(site(from).start_node, hap);

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
                    const gbwt::edge_type np = mosaic_gbwt_position(next_start, hap);
                    // The next segment's haplotype must pass the junction in the same direction as
                    // this one, not only through the same node.
                    const size_t nh2 = strand == 0 ? site(b_idx + 1).hap_first
                                                   : site(b_idx + 1).hap_second;
                    const gbwt::edge_type entry = nh2 == LinkageModel::WILDCARD
                                                      ? gbwt::invalid_edge()
                                                      : mosaic_gbwt_position(next_start, nh2);
                    right = np != gbwt::invalid_edge()
                            && linkage_gbwt->locate(np) == linkage_gbwt->locate(p)
                            && (entry == gbwt::invalid_edge() || entry.first == np.first);
                }
                if (entering) {
                    // The row ends where the child's snarl begins, since the walk reaches the
                    // parent's end only after the child.
                    to_node = next_start;
                    ++mosaic_counters.nested_enter;
                } else if (leaving) {
                    // The row ends at the child's own end, and the next row starts there, so the
                    // stretch from Ce to the parent's end is covered by the parent's haplotype, whose
                    // called allele governs it.
                    to_node = b.end_node;
                    pending_from_node = b.end_node;
                    const size_t nh = strand == 0 ? nx.hap_first : nx.hap_second;
                    pending_from_pos = nh == LinkageModel::WILDCARD
                                           ? gbwt::invalid_edge()
                                           : mosaic_gbwt_position(b.end_node, nh);
                    ++mosaic_counters.nested_leave;
                } else if (right) {
                    to_node = next_start;
                    ++mosaic_counters.extended;
                } else {
                    // Left extension: this segment's haplotype cannot be carried forward, so try
                    // carrying the next segment's haplotype back to this segment's last node, which
                    // closes the gap with a panel haplotype rather than the reference.
                    const size_t nh = strand == 0 ? site(b_idx + 1).hap_first
                                                  : site(b_idx + 1).hap_second;
                    bool closed = false;
                    if (nh != LinkageModel::WILDCARD) {
                        const gbwt::edge_type here = mosaic_gbwt_position(b.end_node, nh);
                        const gbwt::edge_type there = mosaic_gbwt_position(next_start, nh);
                        const gbwt::edge_type mine = mosaic_gbwt_position(b.end_node, hap);
                        if (here != gbwt::invalid_edge() && there != gbwt::invalid_edge()
                            && linkage_gbwt->locate(here) == linkage_gbwt->locate(there)
                            && (mine == gbwt::invalid_edge() || mine.first == here.first)) {
                            pending_from_node = b.end_node;
                            pending_from_pos = here;
                            ++mosaic_counters.extended_left;
                            closed = true;
                        }
                    }
                    if (!closed) {
                        // Neither haplotype crosses the gap, so fill it with the reference, if it
                        // crosses. The fill is contiguous but says little about the sample, so the
                        // row records what it filled.
                        if (patch_gaps && reference_hap != LinkageModel::WILDCARD) {
                            const gbwt::edge_type rl =
                                mosaic_gbwt_position(b.end_node, reference_hap);
                            const gbwt::edge_type rr =
                                mosaic_gbwt_position(next_start, reference_hap);
                            if (rl != gbwt::invalid_edge() && rr != gbwt::invalid_edge()
                                && linkage_gbwt->locate(rl) == linkage_gbwt->locate(rr)) {
                                patch_to = next_start;
                                patch_pos = rl;
                                patch_from_pos = b.position;
                                patch_to_pos = site(b_idx + 1).position;
                                ++mosaic_counters.patched;
                            }
                        }
                        if (patch_to < 0) {
                            ++mosaic_counters.gap_left;
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
                    const gbwt::edge_type at = mosaic_position_at(carry, h);
                    if (at == gbwt::invalid_edge() || !mosaic_follow(at, to_node, end)) {
                        return false;
                    }
                    *pos = at;
                    return true;
                }
                // No direction yet, at the strand's first row: try both.
                for (int o = 0; o < 2; ++o) {
                    const gbwt::edge_type at =
                        mosaic_position_at(gbwt::Node::encode(from_node, o == 1), h);
                    if (at != gbwt::invalid_edge() && mosaic_follow(at, to_node, end)) {
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
                    ++mosaic_counters.row_to_ref;
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
                    ++mosaic_counters.direction_broken;
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
                ++mosaic_counters.row_to_ref;
                direction_broken = carry_applies && row_pos.first != carry;
                if (direction_broken) {
                    ++mosaic_counters.direction_broken;
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
                    << (reference_hap < mosaic_haplotype_names.size()
                            ? mosaic_haplotype_names[reference_hap] : string("?"));
            } else if (hap == LinkageModel::WILDCARD) {
                // The strand passes through here, and the panel cannot name a haplotype for it.
                out << "*\t*";
                ++unexplained_segments;
            } else {
                out << hap << "\t"
                    << (hap < mosaic_haplotype_names.size()
                            ? mosaic_haplotype_names[hap] : string("?"));
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
                ps = mosaic_position_at(s, reference_hap);
                if (ps != gbwt::invalid_edge() && mosaic_follow(ps, patch_to, &pe)) {
                    out << "H\t" << a.contig << "\t" << strand << "\t" << fragment << "\t"
                        << patch_from_pos << "\t" << patch_to_pos << "\t"
                        << ps.first << "\t" << pe << "\t"
                        << "ref\t"
                        << (reference_hap < mosaic_haplotype_names.size()
                                ? mosaic_haplotype_names[reference_hap] : string("?"))
                        << "\t.\t" << ps.second << "\n";
                    ++total_segments;
                    carry = pe;
                } else {
                    ++mosaic_counters.gap_left;
                    --mosaic_counters.patched;
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
                   && mosaic_gbwt_position(site(first_ok).start_node, hap)
                          == gbwt::invalid_edge()) {
                ++first_ok;
            }
            if (first_ok > to) {
                emit_row(from, to, pos);        // clipped across the whole run
                ++mosaic_counters.unwalkable;
                return;
            }
            // The unresolvable head, as small as it really is, then the walkable remainder.
            if (first_ok > from) {
                emit_row(from, first_ok - 1, gbwt::invalid_edge());
                ++mosaic_counters.unwalkable;
                ++mosaic_counters.head_clipped;
            }
            emit_span(first_ok, to, strand, hap, kind);
            return;
        }
        if (pos == gbwt::invalid_edge() || from == to) {
            if (pos == gbwt::invalid_edge() && hap != LinkageModel::WILDCARD) {
                ++mosaic_counters.unwalkable;
            }
            emit_row(from, to, pos);
            return;
        }
        gbwt::edge_type end_pos = mosaic_gbwt_position(site(to).start_node, hap);
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
            gbwt::edge_type p = mosaic_gbwt_position(site(mid).start_node, hap);
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
                // rather than the called alleles. --mosaic-break-unexplained leaves the hole.
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
         << " sites, written to " << mosaic_path << endl;
    // Empty and unexplained segments are reported separately. The unexplained count compares with
    // the phasing report's.
    cerr << "[vg call] mosaic: " << unexplained_segments
         << " segments the panel cannot name a haplotype for" << endl;
    // Segments naming a haplotype the graph does not carry across them, which a consumer has to
    // patch or break at.
    cerr << "[vg call] mosaic: " << mosaic_counters.extended.load()
         << " segment boundaries closed by extending right, " << mosaic_counters.extended_left.load()
         << " by extending left instead, " << mosaic_counters.patched.load()
         << " filled with the reference because neither haplotype could be carried across, "
         << mosaic_counters.gap_left.load() << " left as a gap" << endl;
    cerr << "[vg call] mosaic: " << mosaic_counters.nested_enter.load()
         << " boundaries where the walk enters a child snarl, " << mosaic_counters.nested_leave.load()
         << " where it leaves one -- both stated at the CHILD's boundary node" << endl;
    cerr << "[vg call] mosaic: " << mosaic_counters.direction_broken.load()
         << " rows walked against the carried direction, each standing alone (inversions)" << endl;
    cerr << "[vg call] mosaic: " << mosaic_counters.unwalkable.load()
         << " segments name a haplotype the graph does not carry across them, of which "
         << mosaic_counters.head_clipped.load() << " are a clipped head whose remainder is walkable; "
         << mosaic_counters.row_to_ref.load() << " rewritten as a reference substitution" << endl;
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

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in
// the file again and again
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

void VCFOutputCaller::set_gref_levels(map<string, int> levels) {
    this->gref_levels = std::move(levels);
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
    // non-negative structurally, not by clamping: the loop caps give prefix + suffix <= min_len <=
    // max_len
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
    // and the one MAT's similarity is measured against -- is its best-supported member.  Identity
    // order would instead inherit the traversal finder's ranking, which
    // FlowCaller::call_snarl_internal (and NestedFlowCaller's copy of it) switches to
    // length-weighted average flow once a snarl's interior passes the average-support threshold.
    // That ranking can put a short, lightly-supported allele ahead of a long, heavily-supported
    // one, and merging into it emits the minority sequence as a homozygous call carrying the pooled
    // depth.
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

    // GL is Number=G, and its order depends on the caller that wrote it (see GLLayout), so the
    // layout is passed in rather than assumed.
    auto gl_it = sample.find("GL");
    if (gl_it != sample.end()) {
        size_t n_old = merge_to.size();
        // Only the diploid layout can occur: a merge needs two different called ALTs, so the
        // ploidy is 2.
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
        // todo: check one_based? (we kind of ignore that when writing the snarl name, so maybe not
        // pertienent)
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
thread_local size_t VCFOutputCaller::current_level = 0;
thread_local bool VCFOutputCaller::last_emit_valid = false;

VCFOutputCaller::ChainInlineContext VCFOutputCaller::build_chain_inline_context(
    const Snarl& snarl, const vector<SnarlTraversal>& travs,
    const vector<int>& genotype, int ref_trav_idx) const {
    ChainInlineContext ctx;
    // Only under block emission: with one record per snarl, no chain is inside a block.
    if (!atomize_blocks || symbolic_manager == nullptr) {
        return ctx;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size() || genotype.empty()) {
        return ctx;
    }
    // A snarl whose projection has no symbols cannot answer: every child would read as not
    // reported and be dropped.
    if (!symbolic_site_resolvable(snarl, *symbolic_manager)) {
        return ctx;
    }
    // A genotype with the reference allele matches every reference step, including the chain, so
    // the answer is false for every child.
    for (int allele : genotype) {
        if (allele == ref_trav_idx) {
            return ctx;
        }
    }

    ctx.sref = symbolic_allele(travs[ref_trav_idx], snarl, *symbolic_manager);

    for (int allele : genotype) {
        if (allele < 0 || (size_t)allele >= travs.size()) {
            continue;
        }
        ChainInlineContext::Alt alt;
        alt.salt = symbolic_allele(travs[allele], snarl, *symbolic_manager);
        alt.blocks = symbolic_diff(ctx.sref, alt.salt);
        ctx.alts.push_back(std::move(alt));
    }
    ctx.usable = true;
    return ctx;
}

bool VCFOutputCaller::chain_reported_inline(const ChainInlineContext& ctx,
                                            const Snarl& child) const {
    if (!ctx.usable) {
        return false;
    }
    const Snarl* managed_child = symbolic_manager->into_which_snarl(child.start().node_id(),
                                                                   child.start().backward());
    if (managed_child == nullptr) {
        return false;
    }
    pair<nid_t, nid_t> bounds = chain_bounds_of(managed_child, *symbolic_manager);

    // Where the chain sits in the reference projection. If it is not there, the reference does not
    // cross it, which the caller handles.
    bool in_reference = false;
    for (size_t i = 0; i < ctx.sref.size(); ++i) {
        if (ctx.sref[i].is_chain() && ctx.sref[i].id == bounds.first &&
            ctx.sref[i].end_id == bounds.second) {
            in_reference = true;
            break;
        }
    }
    if (!in_reference) {
        return false;
    }

    bool any_crossing = false;
    for (const ChainInlineContext::Alt& alt : ctx.alts) {
        // This strand's own crossings of the chain, found in its own projection. A strand that
        // deletes the chain has none, so it has no copy for a block to report, even when the
        // reference's step for the chain lies inside a block.
        for (size_t j = 0; j < alt.salt.size(); ++j) {
            if (!alt.salt[j].is_chain() || alt.salt[j].id != bounds.first ||
                alt.salt[j].end_id != bounds.second) {
                continue;
            }
            any_crossing = true;
            bool inside = false;
            for (const DiffBlock& b : alt.blocks) {
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
        // No called strand crosses this chain, so no block ALT spells it.
        return false;
    }

    // Every crossing by every called strand falls inside a difference block, whose ALT spells the
    // route through the chain.
    ++atomize_counters.child_inlined;
    return true;
}

bool VCFOutputCaller::chain_reported_inline(const Snarl& snarl,
                                            const vector<SnarlTraversal>& travs,
                                            const vector<int>& genotype, int ref_trav_idx,
                                            const Snarl& child) const {
    return chain_reported_inline(build_chain_inline_context(snarl, travs, genotype, ref_trav_idx),
                                 child);
}

bool VCFOutputCaller::is_symbolically_reference(const vector<SnarlTraversal>& called_traversals,
                                                int trav_idx, int ref_trav_idx,
                                                const Snarl& snarl) const {
    // Only when symbolic collapsing is on.
    if (symbolic_manager == nullptr || ref_trav_idx < 0 || trav_idx < 0 ||
        ref_trav_idx >= (int)called_traversals.size() ||
        trav_idx >= (int)called_traversals.size()) {
        return false;
    }
    return symbolically_equal(called_traversals[trav_idx], called_traversals[ref_trav_idx],
                              snarl, *symbolic_manager);
}


/// Counters for block emission: project the reference and each distinct called ALT traversal,
/// align them, and count. It changes no output.
static void tally_atomize(const PathPositionHandleGraph& graph, const SnarlManager* mgr,
                          const Snarl& snarl, const vector<SnarlTraversal>& travs,
                          const vector<int>& genotype, int ref_trav_idx,
                          AtomizeCounters& atomize_counters) {
    if (mgr == nullptr || ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size()) {
        return;
    }
    ++atomize_counters.sites;
    bool site_reversed = false;
    if (!symbolic_site_resolvable(snarl, *mgr, &site_reversed)) {
        // The projection would be a bare node list here, so the site is counted and skipped.
        ++atomize_counters.site_unresolvable;
        return;
    }
    if (site_reversed) {
        // Counted here, once per record, rather than in the resolver, which runs once per
        // projection.
        ++atomize_counters.site_reversed;
    }

}

/// A nested ploidy-1 genotype: one allele on a named strand, with "." on the other, since the other
/// strand carries nothing here, its parent allele having deleted the chain. Shared by a site record
/// and its block records.
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
                                        bool genotype_snarls, bool alleles_merged) const {
    // Every refusal below returns -1, meaning the site record is written as it is. Block emission
    // being off is not a refusal, so it is not counted.
    if (!atomize_blocks || symbolic_manager == nullptr || genotype_snarls) {
        return -1;
    }
    if (genotype.empty()) {
        ++atomize_counters.refuse[0];
        return -1;
    }
    if (alleles_merged) {
        ++atomize_counters.refuse[10];
        return -1;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= called_traversals.size()) {
        ++atomize_counters.refuse[1];
        return -1;
    }
    if (!symbolic_site_resolvable(snarl, *symbolic_manager)) {
        // The projection would see no child chains here.
        ++atomize_counters.refuse[2];
        return -1;
    }

    const SnarlTraversal& ref_trav = called_traversals[ref_trav_idx];
    vector<pair<int, int>> ref_ranges;
    SymbolicAllele sref = symbolic_allele(ref_trav, snarl, *symbolic_manager, &ref_ranges);
    const size_t m = sref.size();
    if (m == 0 || ref_ranges.size() != m) {
        ++atomize_counters.refuse[3];
        return -1;
    }

    // Base offset of every visit boundary of the reference traversal from the snarl's first base.
    // The reference traversal is consecutive reference-path steps, so the running sum of node
    // lengths is the offset.
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
    // `max(vb, 0)`, so that the helper never reads visit(-1), even though callers already refuse
    // vb <= 0.
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
        haps[s].blocks = symbolic_diff(sref, haps[s].sym, &haps[s].alt_before_ref);
        if (haps[s].alt_before_ref.size() != m + 1) {
            ++atomize_counters.refuse[4];
            return -1;
        }
    }

    // Cluster all strands' blocks by overlap of their reference step ranges, counting touching
    // ranges as overlapping, so that a deletion on one strand next to an insertion on the other is
    // one record: as two records, the alleles would disagree about the same reference span.
    vector<pair<int, int>> ivs;
    for (const HapAlign& h : haps) {
        for (const DiffBlock& b : h.blocks) {
            ivs.emplace_back(b.ref_begin, b.ref_end);
        }
    }
    if (ivs.empty()) {
        ++atomize_counters.refuse[5];
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

    // Build every record before writing any, so that the decision to split is made on the finished
    // set.
    vector<vcflib::Variant> built;
    built.reserve(clusters.size());
    // For each record, whether no two of its alleles stand for one site allele.
    vector<bool> built_one_to_one;

    for (const pair<int, int>& cluster : clusters) {
        const size_t rb = (size_t)cluster.first;
        const size_t re = (size_t)cluster.second;
        int vb = visit_of_step(ref_ranges, rb, ref_trav);
        int ve = visit_of_step(ref_ranges, re, ref_trav);
        string ref_str = seq_of(ref_trav, vb, ve);

        // Each haplotype's allele over this same reference span, as the visits that spell it: its
        // own visits inside its difference blocks, and the reference's over the steps it matches.
        // A matched chain step is only the same chain, which the haplotype may cross by another
        // route, and that route is the chain's own record to report.
        vector<SnarlTraversal> slot_span(genotype.size());
        vector<string> slot_str(genotype.size());
        vector<bool> slot_marker(genotype.size(), false);
        auto append_visits = [](SnarlTraversal& span, const SnarlTraversal& t, int from, int to) {
            for (int v = std::max(from, 0); v < to && v < t.visit_size(); ++v) {
                *span.add_visit() = t.visit(v);
            }
        };
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (haps[s].trav < 0) {
                slot_marker[s] = true;
                continue;
            }
            SnarlTraversal& span = slot_span[s];
            if (haps[s].trav == ref_trav_idx) {
                append_visits(span, ref_trav, vb, ve);
                slot_str[s] = ref_str;
                continue;
            }
            const SnarlTraversal& t = called_traversals[haps[s].trav];
            // Every block that overlaps or touches the cluster lies inside it, since the clusters
            // are the unions of all blocks, so `next` walks the cluster's reference steps in order.
            size_t next = rb;
            for (const DiffBlock& b : haps[s].blocks) {
                if ((size_t)b.ref_begin <= re && rb <= (size_t)b.ref_end) {
                    append_visits(span, ref_trav, visit_of_step(ref_ranges, next, ref_trav),
                                  visit_of_step(ref_ranges, (size_t)b.ref_begin, ref_trav));
                    append_visits(span, t, visit_of_step(haps[s].ranges, (size_t)b.alt_begin, t),
                                  visit_of_step(haps[s].ranges, (size_t)b.alt_end, t));
                    next = (size_t)b.ref_end;
                }
            }
            append_visits(span, ref_trav, visit_of_step(ref_ranges, next, ref_trav), ve);
            slot_str[s] = seq_of(span, 0, span.visit_size());
        }

        // VCF has no empty allele, so an indel takes the base before it, as flatten_common_allele_ends
        // leaves it.
        bool needs_anchor = ref_str.empty();
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (!slot_marker[s] && slot_str[s].empty()) {
                needs_anchor = true;
            }
        }
        int64_t pos = site_position + (int64_t)ref_visit_off[vb];
        if (needs_anchor) {
            if (vb <= 0) {
                ++atomize_counters.refuse[6];
                // Also stops seq_of(ref_trav, -1, 0) below from reading ref_trav.visit(-1), which
                // can happen when the snarl's start node appears twice in the reference traversal.
                return -1;
            }
            string left = seq_of(ref_trav, vb - 1, vb);
            if (left.empty()) {
                ++atomize_counters.refuse[7];
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

        // Merge alleles with the same sequence, as the site record does, so two strands taking
        // different routes to the same sequence are homozygous.
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
            // The site allele this block allele takes its evidence from: the one the same strand
            // carries at the site.
            auto sa = trav_to_allele.find(haps[s].trav);
            site_of_block.push_back(sa != trav_to_allele.end() ? sa->second : 0);
            block_gt[s] = a;
        }
        if (alleles.size() < 2) {
            continue;   // every haplotype is the reference here; nothing to report
        }

        vcflib::Variant b_var;
        b_var = site;   // inherit INFO, FORMAT and every site-level field
        b_var.position = pos;
        b_var.ref = alleles[0];
        b_var.alt.assign(alleles.begin() + 1, alleles.end());
        b_var.alleles = alleles;
        b_var.info["AT"].clear();
        b_var.info["AT"].resize(alleles.size());
        // The second field, the number of records, is known once every cluster is built.
        b_var.info["SB"] = {std::to_string(built.size()), ""};

        // AT per block allele, over the visit range this record actually spells.
        {
            SnarlTraversal ref_span;
            for (int v = vb; v < ve && v < ref_trav.visit_size(); ++v) {
                *ref_span.add_visit() = ref_trav.visit(v);
            }
            add_allele_path_to_info(b_var, 0, ref_span, false, false);
            for (size_t a = 1; a < alleles.size(); ++a) {
                SnarlTraversal span;
                for (size_t s = 0; s < genotype.size(); ++s) {
                    if (!slot_marker[s] && block_gt[s] == (int)a) {
                        span = slot_span[s];
                        break;
                    }
                }
                add_allele_path_to_info(b_var, a, span, false, false);
            }
        }

        // GT, keeping the site's phase. The site's GT is in site allele space and this record's
        // slots are in genotyper order, so each slot is mapped to the site allele its strand
        // carries. Three forms carry a phase set:
        //
        //   - "a|b", a phased diploid pair, whose order carries over by slot.
        //   - "a|." or ".|a", a nested chain at ploidy 1, on one strand of its parent. A block of
        //     it is part of the same allele on the same strand, so the strand carries over.
        //   - "a", a haploid locus (chrY, or chrX outside the pseudoautosomal regions), with no
        //     order, but PS still labels its phase set.
        //
        // Only a slash-separated GT is unphased, and only there is PS removed.
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
            b_var.samples[sample_name]["GT"] = {gt};
        }
        if (!keep_phase_set) {
            // PS labels a phase block, so it means nothing on an unphased genotype.
            b_var.samples[sample_name].erase("PS");
            b_var.format.erase(std::remove(b_var.format.begin(), b_var.format.end(), "PS"),
                               b_var.format.end());
        }

        // AD and GL are taken from the site, so every block of a snarl reports the same evidence;
        // INFO/SB lets a consumer avoid counting it more than once.
        auto& fmt = b_var.samples[sample_name];
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
                    b_var.format.erase(std::remove(b_var.format.begin(), b_var.format.end(), "GL"),
                                       b_var.format.end());
                }
            }
        }

        b_var.updateAlleleIndexes();
        flatten_common_allele_ends(b_var, true, 0);
        flatten_common_allele_ends(b_var, false, 0);
        built.push_back(std::move(b_var));
        built_one_to_one.push_back(set<int>(site_of_block.begin(), site_of_block.end()).size()
                                   == site_of_block.size());
    }

    if (built.empty()) {
        ++atomize_counters.refuse[8];
        return -1;
    }
    // One block replaces the site record only where the site record says more than the block. A
    // block spells a strand's own visits inside its difference blocks and the reference's over the
    // steps the strand matches, so where a strand crosses a matched child chain, the block spells the
    // reference's route through it, and the chain's own record reports the strand's route. The site
    // record spells each strand's site allele in full, so it says more exactly where some strand's
    // allele takes a route through a matched chain that spells other bases, and there it would repeat
    // the chain's record.
    if (built.size() == 1) {
        // A chain is genotyped from each allele's first crossing of it, so where the reference or a
        // strand crosses a chain more than once, the chain's record need not report the crossing the
        // block leaves to it, and the site record stands.
        auto crosses_a_chain_twice = [](const SymbolicAllele& sym) {
            set<pair<nid_t, nid_t>> chains;
            for (const SymbolicStep& step : sym) {
                if (step.is_chain() && !chains.insert(make_pair(step.id, step.end_id)).second) {
                    return true;
                }
            }
            return false;
        };
        bool chain_crossed_twice = crosses_a_chain_twice(sref);
        bool site_says_more = false;
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (haps[s].trav < 0 || haps[s].trav == ref_trav_idx) {
                continue;
            }
            chain_crossed_twice = chain_crossed_twice || crosses_a_chain_twice(haps[s].sym);
            const SnarlTraversal& t = called_traversals[haps[s].trav];
            string as_blocks;
            size_t next = 0;
            for (const DiffBlock& b : haps[s].blocks) {
                as_blocks += seq_of(ref_trav, visit_of_step(ref_ranges, next, ref_trav),
                                    visit_of_step(ref_ranges, (size_t)b.ref_begin, ref_trav));
                as_blocks += seq_of(t, visit_of_step(haps[s].ranges, (size_t)b.alt_begin, t),
                                    visit_of_step(haps[s].ranges, (size_t)b.alt_end, t));
                next = (size_t)b.ref_end;
            }
            as_blocks += seq_of(ref_trav, visit_of_step(ref_ranges, next, ref_trav), ref_trav.visit_size());
            // The site record spells the strand's site allele, which is the reference for a route that
            // differs from it only inside child chains.
            auto allele = trav_to_allele.find(haps[s].trav);
            const string site_allele =
                allele != trav_to_allele.end() && allele->second == 0
                    ? seq_of(ref_trav, visit_of_step(ref_ranges, 0, ref_trav), ref_trav.visit_size())
                    : seq_of(t, visit_of_step(haps[s].ranges, 0, t), t.visit_size());
            site_says_more = site_says_more || as_blocks != site_allele;
        }
        if (!site_says_more) {
            ++atomize_counters.refuse[9];
            return -1;
        }
        if (chain_crossed_twice) {
            ++atomize_counters.refuse[12];
            return -1;
        }
        // The block takes its genotype, likelihoods and phase from the site, through the site allele
        // each of its alleles stands for, so it cannot be written where two of its alleles stand for
        // one. Two routes can spell one site allele and still differ here, where a difference outside
        // a child chain is cancelled by one inside it.
        if (!built_one_to_one[0]) {
            ++atomize_counters.refuse[11];
            return -1;
        }
    }

    int added = 0;
    size_t block_index = 0;
    for (vcflib::Variant& b_var : built) {
        // The number of records the site writes, which leaves out a cluster that every strand
        // spells as the reference.
        b_var.info["SB"][1] = std::to_string(built.size());
        // VCF allows an ID on one record only, so each block record takes the site's ID with its
        // index appended. A snarl's name has no '_', so the site's ID can be read back from it
        // (see block_site_name).
        b_var.id += "_" + b_var.info["SB"][0];
        if (add_variant(b_var, block_index)) {
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

    // Cleared until this emit fills it in, so that descent after an emit that wrote nothing does
    // not read the previous snarl's state.
    last_emit_valid = false;

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
    // Which VCF allele each called traversal became. Alleles with the same sequence are merged, and
    // only called traversals are written.
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
            // The reference traversal, or one that takes the same route through this snarl and
            // differs only inside child chains, whose own records report those differences.
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

    tally_atomize(graph, symbolic_manager, snarl, called_traversals, genotype, ref_trav_idx,
                  atomize_counters);

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
    // Kept before flattening moves it: the position of the reference traversal's first base, from
    // which block offsets are measured.
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

    // Phase the genotype here, where `trav_to_allele` maps the PhaseCall's traversal pair to this
    // record's allele numbers.
    if (!genotype.empty() && emit_phasing && !render_phases.empty()) {
        auto found = render_phases.find(record_key_of(snarl));
        if (found != render_phases.end()) {
            const LinkageCollector::PhaseCall& phase = found->second;
            // `find`, since `operator[]` would insert a default 0 on a miss, and the map's size is not
            // a bound on traversal indices.
            const auto found_a = trav_to_allele.find(phase.trav_first);
            const auto found_b = trav_to_allele.find(phase.trav_second);
            const int a = (phase.trav_first >= 0 && found_a != trav_to_allele.end())
                              ? found_a->second : -1;
            const int b = (phase.trav_second >= 0 && found_b != trav_to_allele.end())
                              ? found_b->second : -1;
            // The phased genotype must be a permutation of the one this record carries, so that
            // phasing cannot change a genotype.
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
                    // A nested ploidy-1 site is one strand of a diploid locus, since the parent's
                    // other allele deletes the chain. Written as a phased pair with "." on the other
                    // strand, which is how the VCF records which strand carries the allele.
                    genotype_vector[0] = nested_strand_genotype(a, phase.nested_strand);
                } else if (phase.ploidy == 1) {
                    // A haploid locus: one allele and no order; PS labels its phase set. "a|a"
                    // would claim a homozygous diploid call.
                    genotype_vector[0] = std::to_string(a);
                } else {
                    genotype_vector[0] = std::to_string(a) + "|" + std::to_string(b);
                }
                // PS is added after update_vcf_info below, so that it comes last in FORMAT.
                phase_set_to_write = (int64_t)phase.phase_set;
            } else {
                ++phase_declined;
            }
        }
    }

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name, out_variant);

    // PS last in FORMAT.
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

    // Merge near-identical called ALT alleles (vg call -L), turning 1/2 into 1/1. After
    // update_vcf_info, so the genotyper saw every candidate, and after flattening, so the surviving
    // allele's string, POS and REF are the same as without -L. The missing-allele fixup below runs
    // after it and is unaffected, since merging never empties alt. The GL layout depends on which
    // caller wrote the record, so it is passed in.
    const GLLayout gl_layout =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get())
            != nullptr
            ? GLLayout::Colexicographic
            : GLLayout::IMajor;
    const bool alleles_merged = merge_similar_alleles(graph, site_traversals, site_genotype,
                                                      sample_name, out_variant, gl_layout);
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
        // AD is Number=R, written by update_vcf_info for the alleles before this one was added, so it
        // gets an entry for the new allele. GL needs none: it is left out whenever the genotype has
        // a marker, as MISSING_ALLELE_MARKER is.
        auto ad_it = out_variant.samples[sample_name].find("AD");
        if (ad_it != out_variant.samples[sample_name].end()) {
            ad_it->second.push_back("0");
        }
    }

    // Tell descent, which runs next on this thread, that this snarl's traversals are available,
    // whether or not the record is added: a parent written as the reference has no line but still
    // has children to descend into.
    if (symbolic_manager != nullptr) {
        last_emit_valid = true;
    }

    // One record per difference block, where that changes the output. The site record above is
    // finished, so the blocks take every field they do not redefine from it. -1 means the site was
    // declined, and the site record below is written as it is.
    const int block_lines = emit_block_records(graph, snarl, called_traversals, genotype,
                                               ref_trav_idx, sample_name, out_variant,
                                               trav_to_allele, site_position_unflattened,
                                               gl_layout, genotype_snarls, alleles_merged);
    if (block_lines >= 0) {
        ++atomize_counters.split_sites;
        atomize_counters.split_lines += (size_t)block_lines;
        // There is no single line for this snarl, but the linkage model still needs to know
        // whether it has lines, as for a site record below: the mosaic accounts for every site
        // that does. The map is left empty, since each block numbers its own alleles.
        if (linkage_collector != nullptr) {
            linkage_collector->set_allele_map(record_key_of(snarl), vector<int>(),
                                              block_lines > 0);
        }
        return block_lines > 0;
    }

    // Whether this site wants a line. A pair of traversals differing from the reference only inside
    // child chains is written as allele 0 and leaves `alt` empty; such a site has no line, but its
    // children need it.
    const bool wants_line = genotype_snarls || !out_variant.alt.empty();
    bool added = false;
    if (wants_line) {
        added = add_variant(out_variant);
    } else if (include_nested) {
        // A site with nothing to report still knows where its children sit, so its reference
        // interval is kept for their RC, RS and RD, as Deconstructor::deconstruct_site does.
        suppressed_ref_info[omp_get_thread_num()][out_variant.id] =
            {out_variant.sequenceName, static_cast<size_t>(out_variant.position),
             out_variant.ref.length()};
    }
    // The linkage model gets the site whether or not it has a line. A parent written as the
    // reference still has two alleles, which differ only inside its children, and the children
    // need them to know which strand carries the chain. In VCF allele numbering such a parent is
    // 0/0; only in traversal space is it heterozygous.
    if (linkage_collector != nullptr) {
        // The site was recorded when it was genotyped. What remains is the traversal-to-VCF-allele
        // map, which depends on the alleles chosen just above, and whether a line was written.
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
    // True when the record has nothing to write, false only when add_variant refused a line; the
    // linkage pass depends on the difference.
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
        // oops, once of the above checks failed.  we tell caller we coudlnt find by hacking in a -1
        // coordinate.
        start_position = -1;
        end_position = -1;
    }

    if (backward) {
        return make_tuple(end_position, start_position, backward, out_end_step, out_start_step);
    } else {
        return make_tuple(start_position, end_position, backward, out_start_step, out_end_step);
    }
}

/// The 1-based position on the base path of the base `along_path` bases into `ref_path_name`, a
/// path that may name a subrange of its base path.
static int64_t base_path_position(const string& ref_path_name, int64_t along_path) {
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    const int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    return along_path + 1 + basepath_offset;
}

pair<string, int64_t> VCFOutputCaller::get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                                        int64_t ref_path_offset) const {
    const string basepath_name = Paths::strip_subrange(ref_path_name);
    const int64_t position = base_path_position(
        ref_path_name, get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_path_offset);
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

    // With min_allele_len 0 the decrement below would wrap max_flatten_len, a size_t, and the
    // backward pass would read past the end of an empty allele, such as a pure deletion's ALT.
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

string VCFOutputCaller::nesting_info_headers() {
    stringstream ss;
    ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree counting only ancestors whose record is on this record's own reference contig (0=top level for this contig)\">" << endl;
    ss << "##INFO=<ID=CH,Number=1,Type=Integer,Description=\"Nesting steps between VCF reference contigs: how many coordinate-system changes separate this record from a linear reference. Counted as the greater of the in-VCF ancestor hops and the record's own gref contig level, because counting only ancestors that happened to emit a record made a record on a gref fragment whose parent produced no line indistinguishable from one on the linear reference -- 29,843 of 41,669 off-reference records on a gref-covered chr20. So CH >= 1 no longer implies an in-VCF parent, and therefore no longer implies PS\">" << endl;
    ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    ss << "##INFO=<ID=RC,Number=1,Type=String,Description=\"CHROM of the topmost ancestor record in this VCF, or this record's own CHROM when it has none. On a gref fragment, where that own CHROM would be no use, the enclosing site is named even if it produced no record of its own; the tags are absent when there is no such site either\">" << endl;
    ss << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Start of the site named by RC: the POS of its record, or where the site begins when it produced none. A position on that contig, not a span of the snarl, so it can precede the sequence this record describes\">" << endl;
    ss << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"End of the site named by RC: RS plus the length of that site's REF allele\">" << endl;
    return ss.str();
}

string VCFOutputCaller::print_snarl(const HandleGraph* graph, const handle_t& snarl_start,
                                    const handle_t& snarl_end, bool in_brackets) const {
    return print_snarl(graph->get_id(snarl_start), graph->get_is_reverse(snarl_start),
                       graph->get_id(snarl_end), graph->get_is_reverse(snarl_end), in_brackets);
}
string VCFOutputCaller::print_snarl(const Snarl& snarl, bool in_brackets) const {
    return print_snarl(snarl.start().node_id(), snarl.start().backward(), snarl.end().node_id(),
                       snarl.end().backward(), in_brackets);
}
string VCFOutputCaller::print_flipped_snarl(const Snarl& snarl, bool in_brackets) const {
    return print_snarl(snarl.end().node_id(), !snarl.end().backward(), snarl.start().node_id(),
                       !snarl.start().backward(), in_brackets);
}
string VCFOutputCaller::print_snarl(nid_t start_node_id, bool start_backward, nid_t end_node_id,
                                    bool end_backward, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
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
    // Built in place rather than through a stringstream, and without a Snarl message for a
    // flipped or handle-given snarl: names are printed for every snarl of the graph when the VCF
    // is written, and for each ancestor of every record.
    string name;
    name.reserve(start_node.size() + end_node.size() + 4);
    if (in_brackets) {
        name += '(';
    }
    name += start_backward ? '<' : '>';
    name += start_node;
    name += end_backward ? '<' : '>';
    name += end_node;
    if (in_brackets) {
        name += ')';
    }
    return name;
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
            // path and handle oriented the same, we can just backtrack along the path to get
            // previous stuff
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

    // Merge the per-thread suppressed-site intervals collected during calling.  These are sites
    // that never reached the VCF, so pass 1 below cannot see them, but a record nested under one
    // has no other way to name a reference position.
    unordered_map<string, SuppressedRef> suppressed_ref;
    for (auto& buf : suppressed_ref_info) {
        for (auto& kv : buf) {
            suppressed_ref.emplace(kv.first, std::move(kv.second));
        }
        buf.clear();
        buf.rehash(0);
    }

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
    // Whether each interned contig is a synthetic gref fragment, by the same index.  Stored as a
    // bit per contig rather than looked up by name later, so the names are still stored once.
    // is_gref_name(), not is_gref_derived(): a gref copy of a real reference contig is a perfectly
    // good coordinate system -- it is what the whole VCF is deconstructed against -- and only the
    // "_<N>_alt" fragments are positions a reader cannot look up.
    vector<bool> chrom_is_gref_fragment;
    auto intern_chrom = [&](const string& chrom) -> uint32_t {
        auto result = chrom_index.emplace(chrom, (uint32_t)chrom_index.size());
        if (result.second) {
            chrom_is_gref_fragment.push_back(GrefCover::is_gref_name(chrom));
        }
        return result.first->second;
    };
    // One entry per snarl name.  A snarl ID is not unique -- a cyclic reference path that
    // traverses the same snarl twice emits two records with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa) -- but both occurrences are traversals of one
    // path, so they are on the same contig and it does not matter which one wins here.
    unordered_map<string, uint32_t> chrom_of_name;
    // What passes 1 and 2 read from each record: its site's name, CHROM, POS and REF length. The
    // records are decompressed once, in parallel, and the indexes are then filled in record
    // order, as reading the records one by one fills them.
    struct RecordFields {
        string name;
        string chrom;
        string pos;
        size_t ref_len = 0;
        bool top_level = false;
    };
    vector<vector<RecordFields>> record_fields(output_variants.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < output_variants.size(); ++b) {
        vector<RecordFields>& fields = record_fields[b];
        fields.reserve(output_variants[b].size());
        string output_variant_string;
        for (auto& output_variant_record : output_variants[b]) {
            output_variant_string.clear();
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            RecordFields f;
            f.name = block_site_name(toks[2]);
            f.ref_len = toks[3].length();
            f.chrom = std::move(toks[0]);
            f.pos = std::move(toks[1]);
            fields.push_back(std::move(f));
        }
    }
    for (const vector<RecordFields>& fields : record_fields) {
        for (const RecordFields& f : fields) {
            chrom_of_name.emplace(f.name, intern_chrom(f.chrom));
        }
    }

    // index the snarl tree by name
    //
    // Only the names of sites in the VCF are ever looked up, so only they are indexed. The
    // snarls are visited in the same order as for an index of every name, so a name that two
    // snarls print goes to the same one.
    unordered_map<string, const Snarl*> name_to_snarl;
    name_to_snarl.reserve(chrom_of_name.size());
    snarl_manager->for_each_snarl_preorder([&](const Snarl* snarl) {
            string snarl_name = print_snarl(*snarl);
            if (chrom_of_name.count(snarl_name) != 0) {
                name_to_snarl[std::move(snarl_name)] = snarl;
            }
            // also add a map from the flipped snarl (as call sometimes messes with orientation)
            string flipped_name = print_flipped_snarl(*snarl);
            if (chrom_of_name.count(flipped_name) != 0) {
                name_to_snarl[std::move(flipped_name)] = snarl;
            }
        });

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

    // Second pass through variants to extract ref info only for top-level snarls. Whether a
    // record's site is top level depends only on the indexes above, so the records are tested in
    // parallel; the ref info is then stored in record order.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < record_fields.size(); ++b) {
        for (RecordFields& f : record_fields[b]) {
            f.top_level = is_top_level(f.name);
        }
    }
    for (const vector<RecordFields>& fields : record_fields) {
        for (const RecordFields& f : fields) {
            if (f.top_level) {
                top_level_ref_info[f.name][make_pair(f.chrom, static_cast<size_t>(stoul(f.pos)))] =
                    f.ref_len;
            }
        }
    }
    vector<vector<RecordFields>>().swap(record_fields);

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
    // Returns: (contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
    //           suppressed_name)
    // ref_chrom_name is the topmost ancestor in the VCF that sits on a reference contig rather
    // than a gref one, and suppressed_name the topmost ancestor that was dropped for having no
    // variant.  Both feed the RC/RS/RD choice below; neither affects LV/CH/PS.
    function<tuple<size_t, size_t, string, string, string, string>(const string&, const string&)> get_nesting_tags =
        [&](const string& name, const string& my_chrom) {
        string parent_name;
        string ref_chrom_name;
        string suppressed_name;
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

            // Since it is possible that the snarl is actually flipped in the vcf, check for the
            // flipped version too
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
                // Not in the VCF.  If it was dropped for having no variant we still know where it
                // sits, and it may be the only ancestor that can give a reference position.
                auto sup_it = suppressed_ref.find(cur_name);
                if (sup_it == suppressed_ref.end()) {
                    sup_it = suppressed_ref.find(flipped_name);
                }
                if (sup_it != suppressed_ref.end()) {
                    suppressed_name = sup_it->first;
                }
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
            // ...and, separately, the topmost one actually on a reference contig.  An ancestor on
            // another gref contig can name a position, but not one a reader can look up in the
            // reference, so it is the weaker answer of the two.
            if (!chrom_is_gref_fragment[anc_chrom_id]) {
                ref_chrom_name = *hit;
            }
        }
        return make_tuple(contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                          suppressed_name);
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
            // Keyed by the site, so that a block record gets its site's tags.
            const string name = block_site_name(toks[2]);

            auto [contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                  suppressed_name] = get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig, so that a gRef fragment's
            // records start at level 0 on their own contig.
            //
            // CH counts the ancestors that have a record here, so it would be 0 for a record on a
            // gRef fragment whose enclosing site wrote no line. The contig's gRef level, which
            // equals the CH of every record on a fragment, is used as a floor. So CH >= 1 does not
            // imply a parent record in the VCF, or PS.
            size_t gref_level = 0;
            {
                auto it = gref_levels.find(toks[0]);
                if (it != gref_levels.end() && it->second > 0) {
                    gref_level = (size_t)it->second;
                }
            }
            string nesting_tags = ";LV=" + std::to_string(contig_level);
            nesting_tags += ";CH=" + std::to_string(max(contig_hops, gref_level));
            if (!parent_name.empty()) {
                // Not "if (lv != 0)": those were equivalent only while LV was the absolute
                // count.  A record can now legitimately be at LV=0 and still have a parent on
                // another contig, and it must keep PS -- vcfbub's rescue of the children of
                // popped bubbles is keyed on it.
                nesting_tags += ";PS=" + parent_name;
            }

            // Add RC, RS, RD tags: where to look this record up in the reference.
            //
            // Prefer, in order, the topmost ancestor in the VCF that is on a reference contig;
            // then the topmost ancestor in the VCF at all; then the topmost ancestor that was
            // dropped for having no variant.  The last is what rescues a gref fragment whose
            // parent snarl only the reference and its own gref copy span: the site is real and
            // has a reference interval, it just had nothing to report.
            //
            // If none of those exist there is no reference position to give, and the tags are
            // left off.  They used to fall back to this record's own contig and position, which
            // is not a reference coordinate at all -- on a gref contig it is a self-reference
            // that a reader cannot tell apart from the genuine case.
            const string* ref_source = nullptr;
            if (!ref_chrom_name.empty()) {
                ref_source = &ref_chrom_name;
            } else if (top_level_name != name) {
                ref_source = &top_level_name;
            }
            bool have_ref = true;
            RefInfo top_ref;
            if (ref_source == nullptr) {
                if (!GrefCover::is_gref_name(toks[0])) {
                    // Not on a gref fragment, so our own interval is already a position a reader
                    // can look up, and it is the narrower answer of the two.  Keep it rather than
                    // reach for an enclosing site that produced no record: doing that would
                    // repoint every such record on a reference contig at a site LV and CH say it
                    // has no ancestor in.  The records that need the reach are the fragments
                    // below, which have no usable coordinate of their own.
                    top_ref = {toks[0], static_cast<size_t>(stoul(toks[1])), toks[3].length()};
                } else {
                    auto sup_it = suppressed_name.empty() ? suppressed_ref.end()
                                                          : suppressed_ref.find(suppressed_name);
                    if (sup_it != suppressed_ref.end()) {
                        top_ref = {sup_it->second.chrom, sup_it->second.pos, sup_it->second.ref_len};
                    } else {
                        have_ref = false;
                    }
                }
            } else {
                const auto& candidates = top_level_ref_info.at(*ref_source);
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
            if (have_ref) {
                nesting_tags += ";RC=" + top_ref.chrom;
                nesting_tags += ";RS=" + std::to_string(top_ref.pos);
                nesting_tags += ";RD=" + std::to_string(top_ref.pos + top_ref.ref_len);
            }

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

/// Look a reference path up in one of the per-caller maps without inserting on a miss:
/// `operator[]` inserts, and these maps are read from worker threads.
static inline size_t ref_offset_of(const map<string, size_t>& offsets, const string& path) {
    auto it = offsets.find(path);
    return it != offsets.end() ? it->second : 0;
}
static inline int ref_ploidy_of(const map<string, int>& ploidies, const string& path) {
    auto it = ploidies.find(path);
    return it != ploidies.end() ? it->second : 0;
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
                                                            // note: because our traversal finder
                                                            // and support caller have different
                                                            // graphs, they can't share edge handles
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
        // `count` and `at`, not `operator[]`, which inserts on a missing key; these reads run on
        // worker threads.
        int ploidy = ploidy_at(path_name, get<0>(ref_interval),
                               ref_offset_of(ref_offsets, path_name),
                               ref_ploidy_of(ref_ploidies, path_name));
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

void FlowCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {
    GraphCaller::call_top_level_snarls(graph, recurse_type);
    if (show_progress) {
        report_descent_instrumentation();
    }
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

FlowCaller::TraversalNodeIndex FlowCaller::index_traversal_nodes(const SnarlTraversal& trav) {
    TraversalNodeIndex visits;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        visits[trav.visit(i).node_id()].push_back(i);
    }
    return visits;
}

int FlowCaller::crossings_of_child(const TraversalNodeIndex& visits, const Snarl& child) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    // Count crossings: an entry at one boundary followed by the other. Order matters: testing for
    // the two boundaries separately would count a traversal that touches both on unrelated
    // excursions, as find_child_traversal_set does. Only visits to the two boundary nodes can
    // change the state, so we walk their positions merged in ascending order. When both
    // boundaries are the same node, one list plays both roles.
    static const vector<int> none;
    auto s = visits.find(start);
    auto e = visits.find(end);
    const vector<int>& sp = (s == visits.end()) ? none : s->second;
    const vector<int>& ep = (e == visits.end() || end == start) ? none : e->second;
    int crossings = 0;
    nid_t open = 0;
    size_t i = 0, j = 0;
    while (i < sp.size() || j < ep.size()) {
        nid_t node;
        if (j >= ep.size() || (i < sp.size() && sp[i] <= ep[j])) {
            node = start;
            ++i;
        } else {
            node = end;
            ++j;
        }
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
        } else if (open != 0 && node == open) {
            ++crossings;
            open = 0;
        }
    }
    return crossings;
}


int FlowCaller::offset_of_child(const SnarlTraversal& trav, const Snarl& child) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    nid_t open = 0;
    int entry = -1;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        nid_t node = trav.visit(i).node_id();
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
            entry = i;
        } else if (open != 0 && node == open) {
            return entry;   // first complete crossing, entry side
        }
    }
    return -1;
}

int64_t FlowCaller::base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const {
    const int entry = offset_of_child(trav, child);
    if (entry < 0) {
        return -1;
    }
    int64_t bases = 0;
    for (int i = 0; i < entry && i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        bases += (int64_t)graph.get_length(graph.get_handle(trav.visit(i).node_id()));
    }
    return bases;
}

FlowCaller::ChildOffsets::ChildOffsets(const HandleGraph& graph, const SnarlTraversal& trav) {
    bases_before.assign((size_t)trav.visit_size() + 1, 0);
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        int64_t bases = 0;
        if (!visit.has_snarl()) {
            visits_of[visit.node_id()].push_back(i);
            bases = (int64_t)graph.get_length(graph.get_handle(visit.node_id()));
        }
        bases_before[i + 1] = bases_before[i] + bases;
    }
}

int64_t FlowCaller::ChildOffsets::base_offset(const Snarl& child) const {
    // `offset_of_child`'s rule: the entry is the first visit to either boundary node, and it counts
    // only if the other boundary node is visited after it.
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    auto start_visits = visits_of.find(start);
    auto end_visits = visits_of.find(end);
    const int none = numeric_limits<int>::max();
    const int first_start = start_visits == visits_of.end() ? none : start_visits->second.front();
    const int first_end = end_visits == visits_of.end() ? none : end_visits->second.front();
    const int entry = min(first_start, first_end);
    if (entry == none) {
        return -1;
    }
    const auto& closing = first_start <= first_end ? end_visits : start_visits;
    if (closing == visits_of.end()
        || std::upper_bound(closing->second.begin(), closing->second.end(), entry)
               == closing->second.end()) {
        return -1;
    }
    return bases_before[entry];
}

size_t FlowCaller::offset_along_genotype(
    const vector<SnarlTraversal>& travs, const vector<int>& genotype, const Snarl& child,
    unordered_map<const SnarlTraversal*, ChildOffsets>& offsets) const {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const SnarlTraversal* trav = &travs[allele];
        auto found = offsets.find(trav);
        if (found == offsets.end()) {
            found = offsets.emplace(trav, ChildOffsets(graph, *trav)).first;
        }
        const int64_t within = found->second.base_offset(child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

size_t FlowCaller::offset_along_genotype(const vector<SnarlTraversal>& travs,
                                         const vector<int>& genotype, const Snarl& child) const {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const int64_t within = base_offset_of_child(travs[allele], child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

uint64_t FlowCaller::child_crossing_mask(const vector<TraversalNodeIndex>& visits,
                                         const Snarl& child, bool* known) {
    if (known != nullptr) {
        *known = true;
    }
    // One bit per candidate traversal, not per VCF allele, since the linkage model chooses a site's
    // genotype as a pair of traversals.
    if (visits.size() > 64) {
        // Unknown rather than none: the mask cannot index this site's candidates.
        if (known != nullptr) {
            *known = false;
        }
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < visits.size(); ++i) {
        if (crossings_of_child(visits[i], child) > 0) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

int FlowCaller::child_ploidy(const vector<TraversalNodeIndex>& visits, const vector<int>& genotype,
                             const Snarl& child, int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)visits.size()) {
            continue;   // star or missing: that haplotype contributes no copy here
        }
        int crossings = crossings_of_child(visits[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;   // a cycle or tandem duplication; see the header comment
        }
        copies += crossings;
    }
    if (capped) {
        // Counted and reported once per run.
        ++descent_counters.child_multi_crossing;
    }
    return min(copies, cap);
}

void FlowCaller::set_stage_records(bool defer) {
    this->stage_records = defer;
    if (defer) {
        // Sized here rather than inside the parallel region that writes it.
        size_t threads = max((size_t)get_thread_count(), (size_t)omp_get_max_threads());
        // resize, not assign: PendingRecord holds a unique_ptr, so it cannot be copied.
        pending_records.clear();
        pending_records.resize(max(threads, (size_t)1));
        render_records.clear();
        render_records.resize(max(threads, (size_t)1));
    }
}

/// Total over the per-thread queues.
template <typename Queues>
static size_t total_queued(const Queues& queues) {
    size_t n = 0;
    for (const auto& queue : queues) {
        n += queue.size();
    }
    return n;
}

size_t FlowCaller::pending_record_count() const {
    return total_queued(pending_records);
}

size_t FlowCaller::render_record_count() const {
    return total_queued(render_records);
}

// The CallInfo is kept because update_vcf_info reads it when the record is rendered, to map the
// written alleles back to matrix columns, index GL and compute QUAL.
unique_ptr<FlowCaller::PendingRecord> FlowCaller::stage_render_record(
        const Snarl& snarl, const vector<int>& trav_genotype, int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info,
        const string& ref_path_name, int ref_offset, int ploidy) {
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
    rec->level = 0;
    rec->call_info = std::move(call_info);
    // `travs` is not moved here: descent runs after the emit and reads `travs` to find which
    // children the called alleles reach. The caller completes the record after descent.
    return rec;
}

/// The locus of a (path name, position) pair: the contig as the VCF names it, and the position,
/// held at 0 or above.
static FlowCaller::SiteLocus locus_of(pair<string, int64_t> pos_info) {
    const string locus = PathMetadata::parse_locus_name(pos_info.first);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        pos_info.first = locus;
    }
    return FlowCaller::SiteLocus{pos_info.first, (size_t)max((int64_t)0, pos_info.second)};
}

FlowCaller::SiteLocus FlowCaller::site_locus(const Snarl& snarl, const string& ref_path_name,
                                             int ref_offset) const {
    // The position before the record's alleles are trimmed, which can move POS.
    // `get_ref_position` names the base path, as in "CHM13#0#chr20".
    return locus_of(get_ref_position(graph, snarl, ref_path_name, ref_offset));
}

FlowCaller::SiteLocus FlowCaller::off_reference_site_locus(const string& ref_path_name,
                                                           int64_t stand_in_position) const {
    // Not `get_ref_position`: `get_ref_interval` asserts on a snarl the path does not pass
    // through.
    return locus_of(make_pair(ref_path_name, stand_in_position));
}

/// The quality inputs of the direct call `info`, which the linkage collector keeps for rewriting
/// the record if the model moves it. `caller` decides the GQ factor, since --no-share-quality and
/// --depth-quality are its settings.
static LinkageCollector::DirectQuality direct_quality_of(
    const SnarlCaller& caller, const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info) {
    const auto* rl_caller = dynamic_cast<const ReadLikelihoodSnarlCaller*>(&caller);
    return LinkageCollector::DirectQuality{
        .explained_share = info.explained_share,
        .gq_factor = rl_caller != nullptr ? rl_caller->gq_factor(info) : info.explained_share,
        .achievable_gap = info.achievable_gap,
    };
}

void FlowCaller::record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                            const vector<int>& trav_genotype,
                            const unique_ptr<SnarlCaller::CallInfo>& call_info, int ref_trav_idx,
                            const string& ref_path_name, int ref_offset,
                            bool no_reference, int64_t position_from_parent) {
    if (linkage_collector == nullptr) {
        return;
    }
    const auto* rl_info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
    if (rl_info == nullptr) {
        return;
    }
    // The same test the emitter uses, in traversal space: a genotype of one or two alleles, none of
    // them a missing or star marker. Haploid chains are included.
    const size_t site_ploidy = trav_genotype.size();
    if (site_ploidy != 1 && site_ploidy != 2) {
        return;
    }
    for (int allele : trav_genotype) {
        if (allele < 0) {
            return;
        }
    }
    const SiteLocus locus = no_reference
                                ? off_reference_site_locus(ref_path_name, position_from_parent)
                                : site_locus(snarl, ref_path_name, ref_offset);
    const int called_i = trav_genotype[0];
    const int called_j = site_ploidy > 1 ? trav_genotype[1] : called_i;
    // No allele map yet: the written alleles are chosen when the record is built, and
    // `set_allele_map` supplies the map then.
    static const vector<int> no_allele_map;
    linkage_collector->record(
        locus.contig, locus.position,
        rl_info->genotype_lls,
        panel_alleles(graph, travs),
        called_i, called_j, no_allele_map,
        record_key_of(snarl),
        direct_quality_of(snarl_caller, *rl_info), site_ploidy,
        (int64_t)snarl.start().node_id(), (int64_t)snarl.end().node_id(),
        // `nested` only when one copy of the chain is present, as for any other chain; a chain with
        // two copies joins its parent's diploid group.
        LinkageCollector::SiteContext{
            .nested = nested_context.one_copy,
            .parent_record_key = nested_context.parent_record_key,
            .parent_crossing = nested_context.parent_crossing,
            .level = current_level,
            .emitted = false,
            .unpositioned = no_reference,
            .chain_key = nested_context.chain_key,
            .freq_prior = site_freq_prior(travs, ref_trav_idx),
        });
}

double FlowCaller::site_freq_prior(const vector<SnarlTraversal>& travs, int ref_trav_idx) const {
    const LinkageModel::Params& params = linkage_collector->model_params();
    if (params.hp_prior <= 0.0) {
        return -1.0;
    }
    vector<string> alleles;
    alleles.reserve(travs.size());
    for (const SnarlTraversal& trav : travs) {
        alleles.push_back(trav_string(graph, trav));
    }
    const size_t ref = ref_trav_idx >= 0 ? (size_t)ref_trav_idx : (size_t)-1;
    return LinkageModel::run_length_site(alleles, params.hp_prior_run, ref) ? params.hp_prior : -1.0;
}

bool FlowCaller::snarl_is_leaf(const Snarl& snarl) const {
    // Through `manage`, not the address of this Snarl. `SnarlManager::record` casts a Snarl* to its
    // record, which is valid only for a Snarl the manager owns, and the Snarls here are copies.
    // `manage` throws for a snarl the manager does not own, as a nested chain reached by recursion
    // may be, so the call is guarded, and made only when --anchors-leaf-only needs the answer.
    try {
        const Snarl* managed = snarl_manager.manage(snarl);
        return managed != nullptr && snarl_manager.children_of(managed).empty();
    } catch (const std::runtime_error&) {
        // No answer, so treat it as a leaf rather than drop the site.
        return true;
    }
}

unordered_map<size_t, array<int, 3>> FlowCaller::chosen_snapshot() {
    // Each record's chosen pair and ploidy, keyed by record key.
    unordered_map<size_t, array<int, 3>> out;
    if (linkage_collector == nullptr) {
        return out;
    }
    for (const PendingRecord* recp : records_for_render()) {
        int a = -1, b = -1;
        size_t ploidy = 0;
        if (linkage_collector->chosen_traversals(recp->record_key, &a, &b, &ploidy)) {
            out[recp->record_key] = {a, b, (int)ploidy};
        }
    }
    return out;
}

size_t FlowCaller::snapshot_digest(const unordered_map<size_t, array<int, 3>>& snap) {
    // Independent of order, since the snapshot is a hash map: each record's contribution is
    // combined with a commutative mix.
    size_t acc = snap.size() * 1000003ULL;
    for (const auto& kv : snap) {
        size_t h = kv.first;
        h = h * 1000003ULL + (size_t)(kv.second[0] + 3);
        h = h * 1000003ULL + (size_t)(kv.second[1] + 3);
        h = h * 1000003ULL + (size_t)kv.second[2];
        acc ^= h + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
    }
    return acc;
}

size_t FlowCaller::chosen_changed(const unordered_map<size_t, array<int, 3>>& before) {
    const unordered_map<size_t, array<int, 3>> after = chosen_snapshot();
    size_t moved = 0;
    for (const auto& kv : after) {
        auto found = before.find(kv.first);
        if (found == before.end() || found->second != kv.second) {
            ++moved;
        }
    }
    // A record that had a chosen answer and now has none has changed too.
    for (const auto& kv : before) {
        if (after.count(kv.first) == 0) {
            ++moved;
        }
    }
    return moved;
}

const vector<int>& FlowCaller::cached_panel_alleles(PendingRecord& rec) {
    if (!rec.panel_cached) {
        rec.panel_cache = panel_alleles(graph, rec.travs);
        rec.panel_cached = true;
    }
    return rec.panel_cache;
}

vector<FlowCaller::PendingRecord*> FlowCaller::records_for_render(bool for_phasing) {
    vector<PendingRecord*> out;
    out.reserve(render_record_count() + deferred_pending.size());
    for (auto& queue : render_records) {
        for (PendingRecord& rec : queue) {
            out.push_back(&rec);
        }
    }
    for (PendingRecord& rec : deferred_pending) {
        // The same records the hand-off holds back, so that this matches what is rendered: a
        // dropped chain, which the parent's chosen genotype does not carry, and a
        // `reported_inline` one, which an enclosing block's ALT spells. Tested each time, since both
        // can change between linkage passes.
        if (rec.dropped || rec.reported_inline) {
            continue;
        }
        // A chain with no reference path has no REF or POS, so it is not rendered, but it is
        // genotyped, gets anchors, and has a meaningful strand, so phasing includes it.
        if (rec.no_reference && !for_phasing) {
            continue;
        }
        out.push_back(&rec);
    }
    return out;
}

size_t FlowCaller::cascade_nested_strands(vector<LinkageCollector::PhaseCall>& phased,
                                          const std::unordered_map<size_t, size_t>& phase_index,
                                          vector<NestedLink> links,
                                          const unordered_set<size_t>& flips) {
    // A nested site's `nested_strand` was set from its parent's chosen pair when the linkage pass
    // resolved its level, so swapping the parent leaves it naming the other strand. Sites are
    // visited top-down by level, so a parent is done before its children, and each inverts
    // its strand where the meaning of its parent's strand 0 changed:
    //   * under a diploid parent, strand 0 is the parent's first allele, so it changed if the
    //     parent was swapped;
    //   * under a haploid parent, strand 0 is the grandparent's, so it changed if the parent's own
    //     `nested_strand` inverted.
    std::stable_sort(links.begin(), links.end(), [](const NestedLink& a, const NestedLink& b) {
        return a.level < b.level;
    });
    std::unordered_map<size_t, bool> frame_flipped;
    frame_flipped.reserve(links.size() * 2);
    size_t moved = 0;
    for (const NestedLink& link : links) {
        const auto index = phase_index.find(link.key);
        if (index == phase_index.end()) {
            continue;
        }
        LinkageCollector::PhaseCall& pc = phased[index->second];
        bool parent_flipped = false;
        const auto at = frame_flipped.find(link.parent);
        if (at != frame_flipped.end()) {
            parent_flipped = at->second;
        }
        bool strand_moved = false;
        if (pc.nested_strand >= 0 && parent_flipped) {
            pc.nested_strand = pc.nested_strand == 0 ? 1 : 0;
            // The haplotype is held in the slot `nested_strand` names, and the other slot holds the
            // wildcard, which the mosaic reads as an empty strand.
            std::swap(pc.hap_first, pc.hap_second);
            strand_moved = true;
            ++moved;
        }
        frame_flipped[link.key] = pc.ploidy == 2 ? (flips.count(link.key) != 0) : strand_moved;
    }
    return moved;
}

void FlowCaller::apply_read_phasing() {
    if (!read_phasing || linkage_collector == nullptr || linkage_phased.empty()) {
        return;
    }
    // Reset, since re-genotyping calls this again on the new genotypes, and the report should
    // describe the phase the output carries.
    read_phasing_counters = ReadPhasingCounters();
    // Index the phasing by record key, the last one written winning, as in `build_render_phases`.
    std::unordered_map<size_t, size_t> phase_index;
    for (size_t i = 0; i < linkage_phased.size(); ++i) {
        phase_index[linkage_phased[i].record_key] = i;
    }

    // Kept in the member, since re-genotyping uses these sites.
    vector<PhaseSite>& sites = phase_sites;
    sites.clear();
    for (PendingRecord* recp : records_for_render(true)) {
        {
            PendingRecord& rec = *recp;
            const auto found = phase_index.find(rec.record_key);
            if (found == phase_index.end()) {
                continue;
            }
            const LinkageCollector::PhaseCall& pc = linkage_phased[found->second];
            if (pc.ploidy != 2 || pc.trav_first < 0 || pc.trav_second < 0
                || pc.trav_first == pc.trav_second) {
                // Homozygous, haploid, or unplaced: no two strands to order.
                continue;
            }
            const auto* info = dynamic_cast<
                const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
            if (info == nullptr) {
                continue;
            }
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            const size_t a0 = (size_t)pc.trav_first, a1 = (size_t)pc.trav_second;
            if (a0 >= pe->n_alleles || a1 >= pe->n_alleles) {
                // The PhaseCall names a traversal this site's matrix does not have. Skipped, since
                // reading another column would take the phase from the wrong allele.
                continue;
            }
            // Slot order is the PhaseCall's order, so slot 0 is strand 0, as for GT's first field and
            // the anchor file's slot column.
            PhaseSite site = reduce_to_pair(*pe, a0, a1);
            if (site.read_key.empty()) {
                continue;
            }
            site.record_key = rec.record_key;
            site.phase_set = phase_set_id(pc.contig, pc.phase_set);
            site.position = pc.position;
            sites.push_back(std::move(site));
        }
    }
    if (sites.empty()) {
        return;
    }

    phase_flips = read_phase_flips(sites, read_phasing_params, read_phasing_counters);
    const unordered_set<size_t>& flips = phase_flips;

    // Apply by swapping the chosen pair's order. The genotype is the same two traversals either
    // way, so no call changes, only which strand carries which allele. Nested sites are reordered
    // too. Under -A, block records spell the phase in their ALTs, so reordering a nested site can
    // change its GT's allele numbers.
    for (size_t key : flips) {
        const auto found = phase_index.find(key);
        if (found == phase_index.end()) {
            continue;
        }
        LinkageCollector::PhaseCall& pc = linkage_phased[found->second];
        std::swap(pc.trav_first, pc.trav_second);
        std::swap(pc.allele_first, pc.allele_second);
        std::swap(pc.hap_first, pc.hap_second);
    }

    // Carry the swaps down the nesting tree. Every recorded chain is linked, including one whose
    // line an enclosing block's ALT spells (`reported_inline`): it still has anchors, read from
    // its strand, and its children's strands depend on its own. A dropped chain is left out, since
    // the sample does not carry it or anything inside it. Read from the staged sites, since
    // between a linkage pass and the hand-off the nested records are in `deferred_pending`, not in
    // `render_records`.
    vector<NestedLink> links;
    auto link = [&](const PendingRecord& rec) {
        if (!rec.dropped && phase_index.count(rec.record_key) != 0) {
            links.push_back({rec.record_key, rec.parent_record_key, rec.level});
        }
    };
    for (const auto& queue : render_records) {
        for (const PendingRecord& rec : queue) {
            link(rec);
        }
    }
    for (const PendingRecord& rec : deferred_pending) {
        link(rec);
    }
    read_phasing_counters.strands_rederived +=
        cascade_nested_strands(linkage_phased, phase_index, std::move(links), flips);

    const ReadPhasingCounters& c = read_phasing_counters;
    cerr << "[vg call] read phasing: " << c.sites << " het sites, " << c.reliable
         << " reliable, " << c.chains << " blocks, " << c.breaks << " chain breaks ("
         << c.breaks_no_reads << " with no spanning read), " << c.hung
         << " sites hung off the chain (" << c.hung_no_reads << " with no read), " << c.flipped
         << " re-phased against the panel, " << c.strands_rederived
         << " nested strands carried with their parent" << endl;
    if (c.demoted_incoherent > 0) {
        cerr << "[vg call] read phasing: " << c.demoted_incoherent
             << " sites demoted for low phase coherence over " << c.coherence_rounds_run
             << " rounds"
             << (c.coherence_unconverged
                     ? " -- " + std::to_string(c.coherence_unconverged)
                           + " chains were STILL demoting at the round cap, so their chain is not"
                             " a coherent fixed point"
                     : " (every chain reached a coherent fixed point)")
             << endl;
    }
}

bool FlowCaller::apply_regenotyping() {
    if (!regenotype || linkage_collector == nullptr || linkage_phased.empty()) {
        return false;
    }
    // Reset the counters first, before `accumulate_lambda` fills the read counts, so that the report
    // describes this round. The calibration table and fitted temper are kept: they are set once, on
    // the first round.
    const double keep_temper = regenotype_counters.fitted_temper;
    const auto keep_abs = regenotype_counters.fit_abs_lambda;
    const auto keep_obs = regenotype_counters.fit_observed;
    const auto keep_pred = regenotype_counters.fit_predicted;
    const auto keep_n = regenotype_counters.fit_count;
    regenotype_counters = RegenotypeCounters();
    regenotype_counters.fitted_temper = keep_temper;
    regenotype_counters.fit_abs_lambda = keep_abs;
    regenotype_counters.fit_observed = keep_obs;
    regenotype_counters.fit_predicted = keep_pred;
    regenotype_counters.fit_count = keep_n;

    // Lambda over every site read phasing covered, in one pass, into a table keyed by read.
    LambdaTable lambda;
    accumulate_lambda(phase_sites, phase_flips, lambda, regenotype_counters);

    // Each site's phase set, the last PhaseCall written winning, as in `build_render_phases`. A
    // read's strand is usable only at sites of the phase set it was found in.
    unordered_map<size_t, size_t> site_phase_set;
    site_phase_set.reserve(linkage_phased.size() * 2);
    // And the allele the chain puts on strand 0 at each diploid site, against which the reads'
    // preferred order is reported.
    unordered_map<size_t, int> site_strand0;
    site_strand0.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        site_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
        site_strand0[pc.record_key] = pc.ploidy == 2 ? pc.trav_first : -1;
    }

    // Which strand of its parent each nested ploidy-1 chain sits on. `nested_strand` was set in the
    // linkage pass and corrected when its parent's pair was swapped, so it is in the same frame as
    // Lambda for reads of the chain's phase set: strand 0 of that phase set.
    unordered_map<size_t, int> haploid_strand;
    if (regenotype_params.haploid_include) {
        for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
            if (pc.ploidy == 1 && pc.nested_strand >= 0) {
                haploid_strand[pc.record_key] = pc.nested_strand == 0 ? 1 : -1;
            }
        }
    }

    double temper = regenotype_params.temper;
    double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
    if (temper < 0.0) {
        // The temper is fitted once, on the first round, and kept. It describes how reliable the
        // reads' summed strand log-odds are, not which genotypes are called, and fitting it again
        // each round would feed each round's result into the next fit.
        if (regenotype_counters.fitted_temper > 0.0) {
            temper = regenotype_counters.fitted_temper;
            ceiling = regenotype_counters.fitted_ceiling;
        } else {
            fit_calibration(phase_sites, phase_flips, lambda, regenotype_params, temper, ceiling,
                            regenotype_counters);
        }
    } else {
        regenotype_counters.fitted_temper = temper;
        regenotype_counters.fitted_ceiling = ceiling;
    }

    // Each site's own PhaseSite, so that its term can be subtracted from its reads' log-odds.
    unordered_map<size_t, const PhaseSite*> site_by_key;
    site_by_key.reserve(phase_sites.size() * 2);
    for (const PhaseSite& ps : phase_sites) {
        site_by_key[ps.record_key] = &ps;
    }

    ofstream ledger;
    const bool want_ledger = !regenotype_ledger.empty();
    if (want_ledger) {
        ledger.open(regenotype_ledger);
        if (!ledger) {
            cerr << "error [vg call]: cannot write --regeno-ledger " << regenotype_ledger << endl;
            exit(1);
        }
        ledger << "#regeno-ledger-version\t1" << endl;
        ledger << "#temper\t" << temper << endl;
        ledger << "#snarl\tcontig\tposition\tploidy\tcalled\tproposed\tdelta_ln\treads" << endl;
    }

    // Parallel over one flat list of records, strided across threads; `lambda`, `site_by_key` and
    // `phase_flips` are read only. Counters and ledger rows are kept per thread and merged
    // afterwards.
    const vector<PendingRecord*> all_records = records_for_render();
    const size_t n_queues = max<size_t>(1, render_records.size());
    // Only a read-likelihood caller makes the CallInfos corrected below.
    const auto* rl_caller = dynamic_cast<const ReadLikelihoodSnarlCaller*>(&snarl_caller);
    vector<RegenotypeCounters> thread_counters(n_queues);
    // Ledger rows are sorted before writing, since which thread handles a record depends on
    // scheduling.
    struct LedgerRow { string contig; size_t position; string snarl; string text; };
    vector<vector<LedgerRow>> thread_ledger(n_queues);
    vector<size_t> thread_moved(n_queues, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t qi = 0; qi < n_queues; ++qi) {
        RegenotypeCounters& counters = thread_counters[qi];
        unordered_map<uint64_t, double> own;
        size_t moved = 0;
        for (size_t ri = qi; ri < all_records.size(); ri += n_queues) {
            PendingRecord& rec = *all_records[ri];
            auto* info = dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (info == nullptr || rl_caller == nullptr) {
                continue;
            }
            // `converted` belongs to this iteration; nothing may point into it afterwards.
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            // This site's own term, or none: a homozygote has no PhaseSite and contributed nothing to
            // Lambda, so there is nothing to subtract, and it can be corrected into a heterozygote.
            own.clear();
            auto found_site = site_by_key.find(rec.record_key);
            if (found_site != site_by_key.end()) {
                site_own_log_odds(*found_site->second, phase_flips.count(rec.record_key) != 0,
                                  own);
            }

            // The likelihoods before correction, copied only when something reads them.
            map<vector<int>, double> before;
            if (want_ledger) {
                before = info->genotype_lls;
            }
            // At --regeno-passes 1 the correction is computed and reported, and nothing is kept.
            // `genotype_lls` is what GL and QUAL are written from, so correcting it in place would
            // change them while the genotypes stood still.
            const bool keep = regenotype_passes >= 2;
            map<vector<int>, double> scratch;
            if (keep) {
                // Correct the direct pass's likelihoods every round, not the previous round's: the first
                // round saves them, and later rounds restore them before correcting. GQ is
                // restored with them, so that a GQ recomputed in an earlier round does not outlive
                // the correction it was computed from.
                if (info->uncorrected_lls == nullptr) {
                    info->uncorrected_lls.reset(
                        new map<vector<int>, double>(info->genotype_lls));
                } else {
                    info->genotype_lls = *info->uncorrected_lls;
                    rl_caller->recompute_gq(*info);
                }
            } else {
                scratch = info->genotype_lls;
            }
            map<vector<int>, double>& target = keep ? info->genotype_lls : scratch;
            const auto ps = site_phase_set.find(rec.record_key);
            const size_t phase_set = ps != site_phase_set.end() ? ps->second : NO_PHASE_SET;
            const auto s0 = site_strand0.find(rec.record_key);
            const int strand0_allele = s0 != site_strand0.end() ? s0->second : -1;
            const auto hap = haploid_strand.find(rec.record_key);
            const bool site_moved =
                hap != haploid_strand.end()
                    ? haploid_inclusion_correction(*pe, lambda, phase_set, own, temper, ceiling,
                                                   hap->second, regenotype_params, target,
                                                   counters)
                    : phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                             ceiling, regenotype_params, target, counters);
            if (keep && site_moved) {
                // The correction changed the best genotype, so GQ is recomputed from the corrected
                // likelihoods, as the direct pass computes it. GQI and GQN are not: GQN's achievable gap
                // assumes the site's own mixture weights, not per-read ones.
                rl_caller->recompute_gq(*info);
            }
            // Both ploidies, as in the direct pass, since the linkage pass can move a chain from
            // ploidy 1 to 2. At ploidy 1 the correction is zero, but it is applied the same way.
            if (keep && info->alt_ploidy_info != nullptr) {
                auto& alt = *info->alt_ploidy_info;
                if (alt.uncorrected_lls == nullptr) {
                    alt.uncorrected_lls.reset(new map<vector<int>, double>(alt.genotype_lls));
                } else {
                    alt.genotype_lls = *alt.uncorrected_lls;
                    rl_caller->recompute_gq(alt);
                }
                RegenotypeCounters ignored;
                if (phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                           ceiling, regenotype_params, alt.genotype_lls,
                                           ignored)) {
                    rl_caller->recompute_gq(alt);
                }
            }
            if (!site_moved) {
                continue;
            }
            ++moved;
            if (want_ledger) {
                auto best_of = [](const map<vector<int>, double>& gl) {
                    const vector<int>* b = nullptr;
                    double v = -numeric_limits<double>::infinity();
                    for (const auto& kv : gl) {
                        if (kv.second > v) { v = kv.second; b = &kv.first; }
                    }
                    return std::make_pair(b, v);
                };
                auto spell = [](const vector<int>* g) {
                    string out;
                    if (g == nullptr) {
                        return string(".");
                    }
                    for (size_t i = 0; i < g->size(); ++i) {
                        out += (i ? "/" : "") + std::to_string((*g)[i]);
                    }
                    return out;
                };
                const auto a = best_of(before);
                const auto b = best_of(target);
                const string snarl_id = print_snarl(rec.snarl);
                std::ostringstream row;
                row << snarl_id << "\t" << rec.ref_path_name << "\t"
                    << rec.ref_offset << "\t" << rec.ploidy << "\t" << spell(a.first) << "\t"
                    << spell(b.first) << "\t" << (b.second - a.second) << "\t"
                    << pe->num_reads();
                thread_ledger[qi].push_back(
                    LedgerRow{rec.ref_path_name, (size_t)rec.ref_offset, snarl_id, row.str()});
            }
        }
        thread_moved[qi] = moved;
    }
    size_t moved = 0;
    vector<LedgerRow> rows;
    for (size_t qi = 0; qi < n_queues; ++qi) {
        merge_counters(thread_counters[qi], regenotype_counters);
        moved += thread_moved[qi];
        std::move(thread_ledger[qi].begin(), thread_ledger[qi].end(), std::back_inserter(rows));
    }
    if (ledger.is_open()) {
        // The snarl ID breaks ties between records at the same position.
        std::sort(rows.begin(), rows.end(), [](const LedgerRow& x, const LedgerRow& y) {
            if (x.contig != y.contig) return x.contig < y.contig;
            if (x.position != y.position) return x.position < y.position;
            return x.snarl < y.snarl;
        });
        for (const LedgerRow& row : rows) {
            ledger << row.text << endl;
        }
        ledger.close();
    }

    const RegenotypeCounters& c = regenotype_counters;
    cerr << "[vg call] re-genotyping: temper " << temper << ", ceiling " << ceiling << ", "
         << c.reads_with_lambda
         << " reads carry a strand log-odds (" << c.reads_singleton
         << " span one site, so are inert; " << c.reads_multi_phase_set << " span two blocks), "
         << c.sites_corrected << " of " << c.sites_considered << " sites corrected, "
         << c.sites_would_move << " would move (" << c.moved_hom_to_het << " hom->het, "
         << c.moved_het_to_hom << " het->hom, " << c.moved_het_to_het << " het->het), "
         << c.order_reversed << " where the reads prefer the other order" << endl;
    if (regenotype_params.haploid_include) {
        cerr << "[vg call] re-genotyping: " << c.haploid_sites
             << " nested haploid chains weighted by whether the reads belong to their strand, "
             << c.haploid_would_move << " would move" << endl;
    }
    if (show_progress && !c.fit_count.empty()) {
        // Only under --progress: the calibration table is a diagnostic.
        cerr << "[vg call] re-genotyping calibration, |Lambda| / observed / predicted / n:";
        for (size_t i = 0; i < c.fit_count.size(); ++i) {
            cerr << "  " << c.fit_abs_lambda[i] << " " << c.fit_observed[i] << " "
                 << c.fit_predicted[i] << " " << c.fit_count[i];
        }
        cerr << endl;
    }
    return moved > 0;
}

void FlowCaller::rerun_linkage_pass() {
    if (linkage_collector == nullptr) {
        return;
    }
    // Give the linkage model the corrected likelihoods, then run the linkage pass again in full, so that
    // every child is reassessed against its parent's new chosen pair, as on the first pass.
    const vector<PendingRecord*> records = records_for_render();
    // The loop below is serial, and most of its time would go to each record's first
    // `panel_alleles`, a GBWT lookup per allele. Each record's lookup is independent of the others',
    // so fill the caches in parallel first.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < records.size(); ++i) {
        cached_panel_alleles(*records[i]);
    }
    size_t rescored = 0, refused = 0;
    for (PendingRecord* recp : records) {
        PendingRecord& rec = *recp;
        const auto* info = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
            rec.call_info.get());
        if (info == nullptr || info->genotype_lls.empty()) {
            continue;
        }
        // The corrected best genotype, which becomes the entry's called pair.
        const vector<int>* best = nullptr;
        double best_ll = -numeric_limits<double>::infinity();
        for (const auto& kv : info->genotype_lls) {
            if (kv.second > best_ll) {
                best_ll = kv.second;
                best = &kv.first;
            }
        }
        if (best == nullptr || best->empty()) {
            continue;
        }
        const int called_i = (*best)[0];
        const int called_j = best->size() > 1 ? (*best)[1] : called_i;
        if (linkage_collector->rescore(rec.record_key, info->genotype_lls,
                                       cached_panel_alleles(rec), called_i, called_j)) {
            ++rescored;
        } else {
            // The key has no live entry, as for a chain this round has not reinstated, or the
            // corrected likelihoods cannot be compacted. A changed allele space is handled by
            // `rescore`.
            ++refused;
        }
    }
    cerr << "[vg call] re-genotyping: " << rescored << " sites re-scored into the layer, "
         << refused << " refused for want of a live entry or a compactable space" << endl;
    // The linkage pass again, in full.
    run_linkage_pass();
}

void FlowCaller::phase_and_regenotype() {
    // Round 1 ends with read phasing; its linkage pass has already run.
    apply_read_phasing();
    // Each later round re-genotypes from the phase: the current phase gives every read its strand
    // log-odds, the correction rescores every site from the direct pass's likelihoods, the linkage
    // pass chooses the genotypes from the result and reassesses every nested child, and read
    // phasing runs again on the new genotypes. Rounds stop when the correction moves no site's
    // direct call, or when the chosen genotypes stop changing, return to an earlier round's, or
    // reach --regeno-passes rounds. With --regeno-passes 1 the correction is only computed and
    // reported.
    if (regenotype && regenotype_passes >= 2) {
        // Every state the rounds have reached, so that a cycle is recognised. The rounds can cycle:
        // dropping and reinstating a subtree is a discrete change, and the phase is a chain whose
        // links move with the genotypes, so no single quantity must increase.
        vector<size_t> seen_states;
        for (size_t round = 2; round <= regenotype_passes; ++round) {
            const auto before = chosen_snapshot();
            if (round == 2) {
                // Round 1's genotypes.
                seen_states.push_back(snapshot_digest(before));
            }
            const bool calls_moved = apply_regenotyping();
            // Chosen even when the correction moved no direct call: it has already changed every
            // site's likelihoods in place, and GL is written from them, so the genotypes are
            // chosen from them too.
            rerun_linkage_pass();
            apply_read_phasing();
            const auto after = chosen_snapshot();
            const size_t moved = chosen_changed(before);
            cerr << "[vg call] re-genotyping round " << round << ": " << moved
                 << " chosen genotypes moved" << endl;
            if (!calls_moved) {
                cerr << "[vg call] re-genotyping: the correction moved no site's direct call;"
                     << " stopping after round " << round << endl;
                break;
            }
            if (moved == 0) {
                cerr << "[vg call] re-genotyping: converged after " << round << " rounds" << endl;
                break;
            }
            const size_t digest = snapshot_digest(after);
            for (size_t i = 0; i < seen_states.size(); ++i) {
                if (seen_states[i] == digest) {
                    cerr << "[vg call] re-genotyping: LIMIT CYCLE of period "
                         << (seen_states.size() - i) << ", entered at round " << (i + 1)
                         << ". The iteration does not converge and no round of a cycle is more"
                         << " the answer than another; stopping here and reporting it rather than"
                         << " presenting round " << round << " as a fixed point" << endl;
                    goto regeno_done;
                }
            }
            seen_states.push_back(digest);
            if (round == regenotype_passes) {
                if (regenotype_passes == 2) {
                    // The default number of passes.
                    cerr << "[vg call] re-genotyping: one correction round applied; the iteration"
                         << " was not run further (--regeno-passes)" << endl;
                } else {
                    cerr << "[vg call] re-genotyping: NOT CONVERGED and no repeated state seen --"
                         << " still moving " << moved << " genotypes at the round cap of "
                         << regenotype_passes << ". A cycle longer than the rounds run"
                         << " cannot be detected, so raise the cap before concluding there is"
                         << " none" << endl;
                }
            }
        }
    regeno_done:;
    } else {
        // One round: compute and report the correction, and keep nothing.
        apply_regenotyping();
    }
}

void FlowCaller::render_retained_records() {
    // Each read's strand log-odds, for the anchors collected during the render and the hand-off.
    // `phase_sites` and `phase_flips` are final here.
    build_render_lambda();
    // The phase, before any record is built, so that each record is phased as it is rendered. Also
    // before the hand-off, which collects anchors for the records that get no line
    // (`reported_inline` and `no_reference`) and reads `render_phases` to order them.
    build_render_phases();
    // Every linkage pass is done, so the records move to the render, once, which also keeps their
    // anchors from being collected twice.
    hand_off_deferred_records();
    if (render_records.empty()) {
        return;
    }
    // `nested_context` and `current_level` describe the snarl a direct pass thread is recording, and
    // only `record_site` and `call_snarl_internal` read them, neither of which the render calls. The
    // loop still clears them and restores them afterwards, so that it never runs under the context
    // the thread's last swept snarl left. The records are nested chains as well as top-level sites,
    // and each carries its own nesting in its `PendingRecord`.
    const size_t n_threads = render_records.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < n_threads; ++t) {
        NestedContext saved_ctx = nested_context;
        size_t saved_gen = current_level;
        nested_context = NestedContext();
        current_level = 0;
        for (PendingRecord& rec : render_records[t]) {
            // The chosen pair, not the direct pass's. The ALT list, whether a line is written at
            // all, QUAL, and the arity of AD, GL and GQI are all built from the genotype passed in,
            // so they agree with the call.
            vector<int> genotype = chosen_genotype_for(rec);
            // Before emit_variant, which passes the CallInfo on to update_vcf_info. The anchors are
            // collected in phase order, while `genotype` itself stays sorted, since emit_variant
            // builds the ALT list, AD, GL and QUAL from its order.
            collect_anchors_for_record(rec, genotype);
            emit_variant(graph, snarl_caller, rec.snarl, rec.travs, genotype, rec.ref_trav_idx,
                         rec.call_info, rec.ref_path_name, rec.ref_offset, genotype_snarls,
                         rec.ploidy);
        }
        nested_context = saved_ctx;
        current_level = saved_gen;
    }
    if (show_progress) {
        cerr << "[vg call] rendered " << render_record_count()
             << " retained records after the direct pass" << endl;
    }
}

void FlowCaller::run_linkage_pass() {
    if (!stage_records) {
        return;
    }
    // Descent already happened during the direct pass. What is left is to choose the chains' genotypes
    // in the order their ploidies depend on: a level's parents before its children. The reads are not
    // used.
    size_t levels = 0;
    if (linkage_collector != nullptr) {
        levels = linkage_collector->max_level();
    }
    // Merged from the per-thread queues on the first pass only. `deferred_pending` is a member,
    // since re-genotyping runs the linkage pass again.
    vector<PendingRecord>& pending = deferred_pending;
    pending.reserve(pending.size() + pending_record_count());
    for (auto& queue : pending_records) {
        std::move(queue.begin(), queue.end(), std::back_inserter(pending));
        queue.clear();
    }

    // On a later pass, everything a pass concludes is derived again from the chosen genotypes.
    // `dropped` is cleared, since a correction can move a parent onto an allele that crosses a
    // dropped chain; the level loop then drops again the chains still not crossed, and records
    // the others afresh.
    if (linkage_passes_run > 0) {
        for (PendingRecord& pr : pending) {
            pr.dropped = false;
        }
        // Appended to by every resolve, so cleared here.
        linkage_phased.clear();
        // Accumulated by every resolve too.
        linkage_changed = 0;
    }
    ++linkage_passes_run;

    // Parent record key -> indices of its pending children, so that dropping a chain can drop
    // everything under it. Built once, since `pending` does not grow during the linkage pass.
    unordered_map<size_t, vector<size_t>> children_of;
    children_of.reserve(pending.size() * 2);
    for (size_t i = 0; i < pending.size(); ++i) {
        children_of[pending[i].parent_record_key].push_back(i);
    }

    // Counters for the report.
    size_t revise_unrenderable = 0, pass_no_crossing = 0, pass_no_chosen = 0, pass_ploidy_unscored = 0;
    size_t pass_inline_rederived = 0;
    unordered_map<size_t, PendingRecord*> record_by_key;
    record_by_key.reserve((pending.size() + render_record_count()) * 2);
    // Each parent traversal's child offsets, for placing its chains. Keyed by address, which stays
    // valid for the pass as record_by_key's do, and the traversals do not change after the direct
    // pass.
    unordered_map<const SnarlTraversal*, ChildOffsets> child_offsets;
    for (PendingRecord& pr : pending) {
        record_by_key[pr.record_key] = &pr;
    }
    for (auto& queue : render_records) {
        for (PendingRecord& pr : queue) {
            record_by_key[pr.record_key] = &pr;
        }
    }
    // Drop a chain and its whole subtree: the chosen parent does not carry the chain, so the
    // sample has no copy of it or of anything inside it. Returns how many entries were retracted.
    // Iterative, over an explicit stack, since the depth depends on the data.
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

    size_t revised = 0, retracted = 0, gained = 0, crossing_unknown = 0, unspecifiable = 0;
    // One linkage pass per level, in order. `levels` is read again after each pass, since
    // a pass can add a chain at a deeper level, which must still be chosen.
    for (size_t gen = 0; gen <= levels; ++gen) {
        // The final pass has last=true and builds the phasing map and the mosaic from everything
        // accumulated. If it adds a deeper chain, the bound grows and a later pass rebuilds them.
        resolve_linkage_level(gen, gen == levels);

        // This level's parents are chosen, so each chain under one can be given the ploidy
        // its parent's chosen genotype implies before the chain's own level resolves. The
        // direct pass kept the answer at both ploidies, so this is a revision, not a new call.
        // Only the next level's parents are looked up, so only they are indexed; a key's last
        // PhaseCall wins, as it would in an index over every PhaseCall.
        unordered_set<size_t> next_parents;
        for (const PendingRecord& pr : pending) {
            if (pr.level == gen + 1) {
                next_parents.insert(pr.parent_record_key);
            }
        }
        unordered_map<size_t, const LinkageCollector::PhaseCall*> chosen;
        chosen.reserve(next_parents.size() * 2);
        for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
            if (next_parents.count(pc.record_key) != 0) {
                chosen[pc.record_key] = &pc;
            }
        }
        for (size_t i = 0; i < pending.size(); ++i) {
            PendingRecord& pr = pending[i];
            if (pr.level != gen + 1 || pr.dropped) {
                continue;
            }
            auto parent_record = record_by_key.find(pr.parent_record_key);
            if (linkage_collector != nullptr && parent_record != record_by_key.end()) {
                // Place the chain along the allele chosen for its parent, before this level's
                // linkage pass orders and spaces its sites by position. The parent's own offset
                // was placed in the previous level's iteration.
                const PendingRecord& par = *parent_record->second;
                const size_t offset =
                    par.chain_offset
                    + offset_along_genotype(par.travs, chosen_genotype_for(par), pr.snarl,
                                            child_offsets);
                if (offset != pr.chain_offset) {
                    if (pr.no_reference) {
                        pr.position_from_parent += (int64_t)offset - (int64_t)pr.chain_offset;
                        linkage_collector->set_position(
                            pr.record_key,
                            off_reference_site_locus(pr.ref_path_name, pr.position_from_parent)
                                .position);
                    }
                    pr.chain_offset = offset;
                }
            }
            if (!pr.crossing_known) {
                // The direct pass could not compute this chain's crossing mask, because its parent has
                // more candidate traversals than a 64-bit mask can hold. Left as it is and
                // counted, rather than read as "no allele crosses".
                ++crossing_unknown;
                continue;
            }
            if (pr.parent_crossing == 0) {
                // No candidate traversal of the parent crosses the chain, so no chosen genotype
                // can carry it: the sample has no copy, as at ploidy 0 below.
                ++pass_no_crossing;
                retracted += drop_subtree(i);
                continue;
            }
            // The chosen pair as traversals, which the crossing mask is indexed by, through
            // `LinkageCollector::relate_to_parent`, which `resolve_level` also uses to set
            // `nested_strand`.
            int chosen_first = -1, chosen_second = -1;
            bool have_pair = false;
            auto found = chosen.find(pr.parent_record_key);
            if (found != chosen.end()) {
                const LinkageCollector::PhaseCall& parent = *found->second;
                chosen_first = parent.trav_first;
                chosen_second = parent.ploidy == 2 ? parent.trav_second : -1;
                have_pair = true;
            } else if (linkage_collector != nullptr && parent_record != record_by_key.end()
                       && !parent_record->second->genotype.empty()) {
                // The linkage model gave the parent no PhaseCall, so the parent is rendered at its
                // own chosen genotype, which `chosen_genotype_for` reads.
                const vector<int> parent_genotype = chosen_genotype_for(*parent_record->second);
                chosen_first = parent_genotype[0];
                chosen_second = parent_genotype.size() > 1 ? parent_genotype[1] : -1;
                if (chosen_first < 0) {
                    std::swap(chosen_first, chosen_second);
                }
                have_pair = chosen_first >= 0;
            }
            if (!have_pair) {
                // Neither a PhaseCall nor a called allele of the parent can be read, and the chain
                // keeps the ploidy it was called at. Counted.
                ++pass_no_chosen;
                continue;
            }
            const LinkageCollector::Relation rel = LinkageCollector::relate_to_parent(
                pr.parent_crossing, chosen_first, chosen_second);
            int copies = (int)rel.copies;

            // How many copies of the chain the sample has, from the parent's chosen pair. Computed
            // again on every linkage pass, since the ploidy at descent came from the direct pass's
            // genotype.
            if (copies == 0) {
                // The sample has no copy of this chain, and everything inside it is missing too, so
                // the whole subtree is dropped, whether or not it had lines.
                retracted += drop_subtree(i);
                continue;
            }
            if (copies == pr.ploidy && linkage_collector != nullptr
                && linkage_collector->has_entry(pr.record_key)) {
                // The chain was called at the ploidy its parent's chosen genotype implies, so
                // nothing needs revising. `has_entry` matters: a chain that no called parent allele
                // reached in the direct pass is staged but not recorded, and falling through records it.
                continue;
            }

            // Leave alone a record that cannot be built: no traversals, or a genotype out of range.
            // emit_variant indexes `called_traversals[ref_trav_idx]` unchecked, so a record with a
            // reference path also needs a valid ref_trav_idx. A record with no reference path skips
            // that check, since it is never emitted.
            if (pr.travs.empty()
                || (!pr.no_reference
                    && (pr.ref_trav_idx < 0 || (size_t)pr.ref_trav_idx >= pr.travs.size()))) {
                ++revise_unrenderable;
                continue;
            }
            bool genotype_in_range = !pr.genotype.empty();
            for (int allele : pr.genotype) {
                if (allele >= 0 && (size_t)allele >= pr.travs.size()) {
                    genotype_in_range = false;
                }
            }
            if (!genotype_in_range) {
                ++revise_unrenderable;
                continue;
            }

            // Build the record at the ploidy the chosen parent implies, from the answers kept in the
            // direct pass; `alt_ploidy_info` holds the other ploidy's.
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
                // The two answers are exchanged, not one discarded, so that a chain can follow its
                // parent to either ploidy however often the linkage pass runs.
                use_genotype = rl->alt_ploidy_best;
                const vector<int> demoted_genotype = pr.genotype;
                unique_ptr<SnarlCaller::CallInfo> demoted = std::move(pr.call_info);
                // `rl` still points at it -- `demoted` owns what `pr.call_info` did.
                unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo> promoted(
                    rl->alt_ploidy_info.release());
                // The fields that do not depend on ploidy go with whichever answer is in front, since
                // the alternate does not copy them.
                promoted->anchor_evidence = std::move(rl->anchor_evidence);
                promoted->phase_evidence = std::move(rl->phase_evidence);
                // The replaced answer becomes the new alternate, with the genotype it was called at.
                promoted->alt_ploidy_best = demoted_genotype;
                promoted->alt_ploidy_info.reset(
                    static_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                        demoted.release()));
                use_info = std::move(promoted);
            } else {
                // No answer at that ploidy: the direct pass computed none, because the chain offers too
                // few traversals for a second genotype. The chain keeps a ploidy its chosen parent
                // contradicts, which is counted.
                ++pass_ploidy_unscored;
                continue;
            }
            // Whether this chain was already in the linkage model: whether it is being revised or
            // added, and whether there is an old entry to retract.
            const bool had_entry = linkage_collector != nullptr
                                   && linkage_collector->has_entry(pr.record_key);
            // Revise the staged site; the render builds its line once, at the end, from the chosen
            // genotype.
            pr.genotype = use_genotype;
            pr.ploidy = copies;
            if (use_info != nullptr) {
                pr.call_info = std::move(use_info);
            }
            const unique_ptr<SnarlCaller::CallInfo>& info = pr.call_info;

            const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* used =
                dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(info.get());
            if (used != nullptr) {
                // In traversal space, as `record_site` records it, so the linkage pass and the
                // direct pass describe a site the same way. No allele map yet, as in `record_site`.
                static const vector<int> no_allele_map;
                const vector<int>& trav_to_allele_vec = no_allele_map;
                const SiteLocus locus =
                    pr.no_reference
                        ? off_reference_site_locus(pr.ref_path_name, pr.position_from_parent)
                        : site_locus(pr.snarl, pr.ref_path_name, pr.ref_offset);
                int called_i = use_genotype.empty() ? -1 : use_genotype[0];
                int called_j = use_genotype.size() > 1 ? use_genotype[1] : called_i;
                const vector<int>& panel = cached_panel_alleles(pr);
                // Retract the old entry and record the site again, so that there is one way a site
                // enters the linkage model. Retract first: `live_index` returns the first live entry
                // for a key, so the new entry is the live one.
                if (had_entry) {
                    linkage_collector->retract(pr.record_key);
                }
                linkage_collector->record(
                    locus.contig, locus.position, used->genotype_lls, panel,
                    called_i, called_j, trav_to_allele_vec,
                    // The explained share the old entry carried, or 1.0 for a chain that had none.
                    // The quality inputs of the direct call recorded here, which the record's
                    // GQI and GL also come from.
                    pr.record_key, direct_quality_of(snarl_caller, *used),
                    (size_t)copies, pr.snarl.start().node_id(), pr.snarl.end().node_id(),
                    LinkageCollector::SiteContext{
                        .nested = copies == 1,
                        .parent_record_key = pr.parent_record_key,
                        .parent_crossing = pr.parent_crossing,
                        .level = pr.level,
                        .emitted = false,
                        .unpositioned = pr.no_reference,
                        .chain_key = pr.chain_key,
                        .freq_prior = site_freq_prior(pr.travs, pr.ref_trav_idx),
                    });
                if (!linkage_collector->has_entry(pr.record_key)) {
                    // `record` adds nothing for a site whose compact space it cannot describe: no called
                    // traversal, no likelihoods, or more than 127 alleles. The old entry stays
                    // retracted, and the record keeps its per-site call.
                    if (had_entry) {
                        ++unspecifiable;
                    }
                }
            }
            if (had_entry) {
                ++revised;
            } else {
                ++gained;
            }

            // This chain's chosen pair has changed, or the chain is new, so its children's crossing
            // masks are computed again, through `children_of`.
            const auto kids = children_of.find(pr.record_key);
            if (kids != children_of.end()) {
                // Once for this parent: see TraversalNodeIndex.
                vector<TraversalNodeIndex> pr_visits;
                pr_visits.reserve(pr.travs.size());
                for (const SnarlTraversal& t : pr.travs) {
                    pr_visits.push_back(index_traversal_nodes(t));
                }
                for (size_t ci : kids->second) {
                    PendingRecord& child = pending[ci];
                    bool known = true;
                    child.parent_crossing = child_crossing_mask(pr_visits, child.snarl, &known);
                    child.crossing_known = known;
                }
            }
        }
        if (linkage_collector != nullptr) {
            levels = max(levels, linkage_collector->max_level());
        }
    }

    // The exactly-once test, from the chosen genotypes the render builds each parent's blocks
    // from: `reported_inline` holds back a chain's line where an enclosing block's ALT spells it.
    // Without the linkage model every chosen genotype is the direct call the direct pass tested, so
    // there is nothing to redo.
    if (linkage_collector != nullptr && !children_of.empty()) {
        // Parents before their children, so that a chain inherits its parent's final flag.
        vector<pair<uint8_t, size_t>> parents;
        parents.reserve(children_of.size());
        for (const auto& kv : children_of) {
            auto parent = record_by_key.find(kv.first);
            if (parent != record_by_key.end()) {
                parents.emplace_back(parent->second->level, kv.first);
            }
        }
        sort(parents.begin(), parents.end());
        // Counted again from here, so that the report gives the chains held back now.
        atomize_counters.child_inlined = 0;
        for (const pair<uint8_t, size_t>& gk : parents) {
            const PendingRecord& parent = *record_by_key.at(gk.second);
            if (parent.dropped) {
                continue;   // its children were dropped with it
            }
            // The parts of the test that do not depend on the child, built once for this parent;
            // see VCFOutputCaller::ChainInlineContext.
            const ChainInlineContext ctx = build_chain_inline_context(
                parent.snarl, parent.travs, chosen_genotype_for(parent), parent.ref_trav_idx);
            for (size_t ci : children_of.at(gk.second)) {
                PendingRecord& child = pending[ci];
                if (child.dropped) {
                    continue;
                }
                const bool was = child.reported_inline;
                child.reported_inline =
                    parent.reported_inline || chain_reported_inline(ctx, child.snarl);
                if (was != child.reported_inline) {
                    ++pass_inline_rederived;
                }
            }
        }
    }
    if (show_progress) {
        // The bytes kept for the staged sites, counted by walking the objects. They are walked in
        // parallel; the totals are sums, so they do not depend on how the walk is split.
        size_t retained_bytes = 0, retained_visits = 0, retained_gls = 0;
        auto measure = [](const PendingRecord& rec, size_t& bytes, size_t& visits, size_t& gls) {
            bytes += sizeof(PendingRecord) + rec.ref_path_name.capacity()
                     + rec.genotype.capacity() * sizeof(int)
                     + rec.panel_cache.capacity() * sizeof(int);
            bytes += rec.travs.capacity() * sizeof(SnarlTraversal);
            for (const SnarlTraversal& t : rec.travs) {
                visits += (size_t)t.visit_size();
                bytes += (size_t)t.visit_size() * sizeof(Visit);
            }
            const auto* rl = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (rl != nullptr) {
                for (const auto& kv : rl->genotype_lls) {
                    ++gls;
                    bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                }
                if (rl->anchor_evidence != nullptr) {
                    bytes += rl->anchor_evidence->bytes();
                }
                if (rl->phase_evidence != nullptr) {
                    bytes += rl->phase_evidence->bytes();
                }
                // The parts re-genotyping adds.
                auto gl_bytes = [](const map<vector<int>, double>& gl) {
                    size_t n = 0;
                    for (const auto& kv : gl) {
                        n += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                    return n;
                };
                if (rl->uncorrected_lls != nullptr) {
                    bytes += gl_bytes(*rl->uncorrected_lls);
                }
                bytes += rl->scored_traversals.capacity() * sizeof(SnarlTraversal)
                         + rl->allele_support.capacity() * sizeof(double);
                if (rl->alt_ploidy_info != nullptr) {
                    // The alternate answer is kept too, with all its parts.
                    const auto& alt = *rl->alt_ploidy_info;
                    bytes += alt.scored_traversals.capacity() * sizeof(SnarlTraversal)
                             + alt.allele_support.capacity() * sizeof(double);
                    if (alt.uncorrected_lls != nullptr) {
                        bytes += gl_bytes(*alt.uncorrected_lls);
                    }
                }

                if (rl->alt_ploidy_info != nullptr) {
                    for (const auto& kv : rl->alt_ploidy_info->genotype_lls) {
                        ++gls;
                        bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                }
            }
        };
#pragma omp parallel reduction(+ : retained_bytes, retained_visits, retained_gls)
        {
            for (const auto& queue : render_records) {
#pragma omp for schedule(dynamic, 4096) nowait
                for (size_t r = 0; r < queue.size(); ++r) {
                    measure(queue[r], retained_bytes, retained_visits, retained_gls);
                }
            }
#pragma omp for schedule(dynamic, 4096) nowait
            for (size_t r = 0; r < pending.size(); ++r) {
                measure(pending[r], retained_bytes, retained_visits, retained_gls);
            }
        }
        // The read-phasing evidence. In the report below, the snarls the linkage pass will not revise
        // are the top-level ones and the children RecurseOnFail reaches without a ploidy override.
        for (const PhaseSite& ps : phase_sites) {
            retained_bytes += sizeof(PhaseSite) + ps.read_key.capacity() * sizeof(uint64_t)
                              + ps.q0.capacity() * sizeof(float) + ps.c.capacity() * sizeof(float);
        }
        retained_bytes += phase_flips.size() * (sizeof(size_t) + 16);
        cerr << "[vg call] retained for rendering: " << render_record_count()
             << " snarls the linkage pass will not revise, plus " << pending.size()
             << " nested chains; " << (retained_bytes / (1024.0 * 1024.0)) << " MB over "
             << retained_visits << " traversal visits and " << retained_gls
             << " genotype likelihoods" << endl;
        cerr << "[vg call] linkage pass exits: " << pass_no_crossing
             << " dropped because no parent candidate crosses them, " << pass_no_chosen
             << " whose parent's chosen pair could not be read, " << revise_unrenderable
             << " unrenderable so left unrevised, " << pass_ploidy_unscored
             << " stranded at a ploidy the direct pass never scored" << endl;
        if (pass_inline_rederived > 0) {
            cerr << "[vg call] linkage pass: " << pass_inline_rederived
                 << " children whose exactly-once suppression changed with their parent's"
                 << " chosen genotype" << endl;
        }
        cerr << "[vg call] single direct pass: " << pending.size() << " nested chains retained over "
             << (levels + 1) << " levels; " << revised << " revised, " << gained
             << " reachable only under the chosen parent, " << retracted << " retracted";
        if (crossing_unknown > 0) {
            cerr << ", " << crossing_unknown << " with a crossing mask the direct pass could not compute";
        }
        if (unspecifiable > 0) {
            cerr << ", " << unspecifiable << " dropped from the layer because the site's "
                 << "compact allele space could not be built";
        }
        cerr << endl;
    }

}

void FlowCaller::hand_off_deferred_records() {
    if (!stage_records) {
        return;
    }
    vector<PendingRecord>& pending = deferred_pending;
    // Hand every surviving chain to the render, so that nested and top-level records are written
    // in one place from their chosen genotypes. A dropped chain is not handed over, since the
    // sample has no copy of it. Spread over the queues, since the render is parallel over them.
    size_t next_queue = 0;
    size_t no_ref_unrendered = 0, inline_unrendered = 0;
    for (PendingRecord& pr : pending) {
        if (pr.dropped) {
            continue;
        }
        if (pr.reported_inline) {
            // An enclosing block's ALT already spells its variation, so it gets no line, but it still
            // gets anchors.
            collect_anchors_for_record(pr, chosen_genotype_for(pr));
            ++inline_unrendered;
            continue;
        }
        if (pr.no_reference) {
            // No reference path, so no REF or POS, and no line. Held back here, since the render calls
            // emit_variant for every record it is given. It still gets anchors, which are placed by
            // node ID.
            collect_anchors_for_record(pr, chosen_genotype_for(pr));
            ++no_ref_unrendered;
            continue;
        }
        if (render_records.empty()) {
            render_records.resize(max((size_t)1, (size_t)get_thread_count()));
        }
        render_records[next_queue % render_records.size()].push_back(std::move(pr));
        ++next_queue;
    }
    if (inline_unrendered > 0 && show_progress) {
        cerr << "[vg call] block emission: " << inline_unrendered
             << " chains genotyped, recorded and phased, but left unrendered because an enclosing"
             << " block's ALT already spells them out" << endl;
    }
    if (no_ref_unrendered > 0 && show_progress) {
        cerr << "[vg call] off-reference nested: " << no_ref_unrendered
             << " chains chosen and left unrendered, having no reference position to write" << endl;
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

    // Staged in the nested branch below and completed after descent, which reads `travs`, since
    // this record then takes ownership of them.
    unique_ptr<PendingRecord> pending_this;
    // The same, for a snarl the linkage pass will not revise.
    unique_ptr<PendingRecord> render_this;
    // Whether this call ran emit_variant, so that descent knows `last_emit_valid` describes this
    // snarl. A retained chain, or a GAF run, skips the emit.
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
    //
    // Only the flow traversal finder uses greedy_avg_flow, so the sum is computed only when there
    // is one.
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.deep_contents(&snarl, graph, false);
        if (snarl_contents.second.size() > max_snarl_edges) {
            // size cap needed as non-nested FlowCaller doesn't handle large snarls
            return false;
        }
        if (flow_trav_finder != nullptr) {
            size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
            size_t length = 0;
            for (auto i = snarl_contents.first.begin();
                 i != snarl_contents.first.end() && length < len_threshold; ++i) {
                length += graph.get_length(graph.get_handle(*i));
            }
            greedy_avg_flow = length > len_threshold;
        }
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
        // This test and the use_parent_interval test below must agree: otherwise get_ref_interval
        // would be called with the parent's reference path, which does not visit this snarl's
        // boundary nodes, and would assert.
        if ((parent_child_trav_sets == nullptr && !nested_context.no_reference)
            || parent_ref_path_name.empty()) {
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

    if (common_names.empty()) {
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
        // Left at -1 where no parent traversal set named a reference: travs[0] is the flow finder's
        // best-supported traversal, and making it REF would bias the genotyper's tie-break toward
        // it. The genotyper checks ref_trav_idx >= 0 before using it.
        if (ref_trav_idx < 0 && parent_child_trav_sets != nullptr) {
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

    // A ploidy from the parent overrides the contig's or the region BED's: it is the number of the
    // parent's called alleles that reach this child. The region's is still the number of the
    // sample's haplotypes here, which the depth term needs.
    const int region_ploidy = ploidy_at(ref_path_name, get<0>(ref_interval),
                                        ref_offset_of(ref_offsets, ref_path_name),
                                        ref_ploidy_of(ref_ploidies, ref_path_name));
    int ploidy = ploidy_override >= 0 ? ploidy_override : region_ploidy;

    // What both the parent-traversal-set branch and the top-level branch do with their genotype.
    // `trav_call_info` differs between them, so it is a parameter. `snarl` is captured by reference;
    // `flip_snarl` may already have rewritten it above.
    auto stage_or_emit = [&](unique_ptr<SnarlCaller::CallInfo>& trav_call_info) -> bool {
        bool added;
        if (!gaf_output) {
            // Staged, not emitted: `render_retained_records` writes it after the direct pass.
            // `added` stands in for emit_variant's return value, which here only gates recursion;
            // a staged site counts as added.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_trav_idx, ref_path_name,
                        ref_offset_of(ref_offsets, ref_path_name));
            render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                              ref_path_name, ref_offset_of(ref_offsets, ref_path_name), ploidy);
            added = render_this != nullptr;
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offset_of(ref_offsets, ref_path_name),
                                     genotype_snarls, ploidy);
            }
            emitted_this_call = true;
        } else {
            added = true;
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }
        return added;
    };


    if (traversals_only) {
        assert(gaf_output);
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offset_of(ref_offsets, ref_path_name));
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else if (parent_child_trav_sets != nullptr && !parent_child_trav_sets->empty()) {
        // Genotype using bounded search over traversal sets from parent
        // Each set contains traversals consistent with one parent allele
        ploidy = parent_child_trav_sets->size();

        // Track which set each traversal index belongs to (for phase consistency)
        // set_membership[i] = which parent allele set traversal i came from, or -1 if from finder
        vector<int> set_membership(travs.size(), -1);

        // Merge traversals from sets into travs, tracking membership

        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            const TraversalSet& tset = (*parent_child_trav_sets)[set_idx];

            if (tset.empty()) {
                // Empty set means parent allele doesn't traverse this child (star allele)
                continue;
            }

            // Add traversals from this set to travs (avoiding duplicates)
            // Keep track of indices for this set
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
            // Genotype at the ploidy that passes through the site, not at the parent's ploidy, since
            // a site only one strand reaches is not diploid. genotype() returns a sorted multiset, so
            // the alleles are then placed on the strands that pass through.
            int effective_ploidy = (int)traversing_sets.size();
            vector<int> called_alleles;
            ReadLikelihoodSnarlCaller::set_region_ploidy(region_ploidy);
            std::tie(called_alleles, trav_call_info) = snarl_caller.genotype(
                snarl, travs, ref_trav_idx, effective_ploidy, ref_path_name,
                make_pair(get<0>(ref_interval), get<1>(ref_interval)));
            ReadLikelihoodSnarlCaller::set_region_ploidy(0);

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
        } else {
            added = stage_or_emit(trav_call_info);
        }

        ret_val = trav_genotype.size() == ploidy && added;
    } else if (ploidy_override >= 0) {
        // A nested chain, reached by descent, at the ploidy its parent implied. Only a nested chain
        // can have its ploidy revised at the linkage pass, so only it needs the other ploidy's answer.
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(true);
        ReadLikelihoodSnarlCaller::set_region_ploidy(region_ploidy);
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(
            snarl, travs, ref_trav_idx, ploidy, ref_path_name,
            make_pair(get<0>(ref_interval), get<1>(ref_interval)));
        ReadLikelihoodSnarlCaller::set_region_ploidy(0);
        ReadLikelihoodSnarlCaller::set_want_alt_ploidy(false);

        const bool retain_only = nested_context.retain_only;
        // Whether this snarl's own boundaries are on no reference path, checked from the graph for
        // each snarl.
        const bool no_ref_position = use_parent_interval;

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
        bool added = true;
        if (no_ref_position) {
            // Genotyped and recorded, never written. Checked before retain_only, which does not
            // record. `added` is true, as for retain_only, since it gates descent into this chain's
            // children.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_trav_idx, ref_path_name,
                        ref_offset_of(ref_offsets, ref_path_name), /*no_reference*/ true,
                        // The parent's position, as `get_ref_position` gives it from the interval
                        // `use_parent_interval` set, plus the chain's offset along its parent, as
                        // `PendingRecord::position_from_parent` has it.
                        base_path_position(ref_path_name, get<0>(ref_interval)
                                                              + ref_offset_of(ref_offsets, ref_path_name))
                            + (int64_t)nested_context.parent_offset);
            ++descent_counters.no_ref_recorded;
            {
                int copies = 0;
                for (int a : trav_genotype) {
                    copies += (a >= 0);
                }
                descent_counters.no_ref_copies[copies < 3 ? copies : 2].fetch_add(1);
            }
            added = true;
        } else if (retain_only) {
            // No called parent allele reaches this chain, so nothing about it is written yet. It is
            // genotyped and kept, since the linkage model may still move the parent onto an allele
            // that reaches it.
            added = true;
        } else if (nested_context.reported_inline) {
            // An enclosing block's ALT already spells this chain, so it gets no line, but it is
            // genotyped and recorded, since its allele pair phases everything inside it. Checked
            // after retain_only, which does not record.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_trav_idx, ref_path_name,
                        ref_offset_of(ref_offsets, ref_path_name));
            added = true;
        } else if (!gaf_output) {
            // Recorded here rather than in emit_variant. A retained chain, on the path above, is
            // recorded only if the linkage pass later finds that the sample carries it.
            record_site(snarl, travs, trav_genotype, trav_call_info, ref_trav_idx, ref_path_name,
                        ref_offset_of(ref_offsets, ref_path_name));
            // Staged, not emitted, as at top level: the line is written after the linkage pass, from the
            // chosen genotype. `added` stands in for emit_variant's return value, which here only
            // gates recursion; a staged site counts as added.
            added = stage_records && !pending_records.empty();
            if (!added) {
                added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx,
                                     trav_call_info, ref_path_name, ref_offset_of(ref_offsets, ref_path_name),
                                     genotype_snarls, ploidy);
                emitted_this_call = true;
            }
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name,
                                                              ref_offset_of(ref_offsets, ref_path_name));
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx,
                             pos_info.first, pos_info.second, &support_finder);
        }

        // Stage the nested site without its traversals: descent below still reads `travs` to find
        // which children the called alleles reach, and they are moved in once descent is done.
        if (stage_records && !pending_records.empty()) {
            pending_this.reset(new PendingRecord());
            pending_this->snarl = snarl;
            pending_this->ref_path_name = ref_path_name;
            pending_this->ref_offset = ref_offset_of(ref_offsets, ref_path_name);
            pending_this->ref_trav_idx = ref_trav_idx;
            pending_this->genotype = trav_genotype;
            pending_this->ploidy = ploidy;
            pending_this->record_key = record_key_of(snarl);
            pending_this->parent_record_key = nested_context.parent_record_key;
            pending_this->parent_crossing = nested_context.parent_crossing;
            pending_this->chain_key = nested_context.chain_key;
            pending_this->no_reference = no_ref_position;
            pending_this->reported_inline = nested_context.reported_inline;
            pending_this->position_from_parent =
                no_ref_position
                    ? base_path_position(ref_path_name,
                                         get<0>(ref_interval) + ref_offset_of(ref_offsets, ref_path_name))
                          + (int64_t)nested_context.parent_offset
                    : 0;
            pending_this->chain_offset = nested_context.parent_offset;
            pending_this->crossing_known = nested_context.crossing_known;
            pending_this->level = (uint8_t)min(current_level, (size_t)255);
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
        added = stage_or_emit(trav_call_info);

        ret_val = trav_genotype.size() == ploidy && added;
    }

    // Nested calling: descend into each child the called alleles reach, at the ploidy they reach
    // it with.
    //
    // Descent does not depend on whether a line was written: a parent written as the reference
    // still has children to call. Children are genotyped independently, with no parent traversal
    // sets. Only a successful call descends, since a failed snarl has no genotype to take a child's
    // ploidy from. RecurseOnFail calls the children of a failed top-level snarl as top-level
    // snarls, but nothing does so for a failed nested snarl: its children are not called.
    if (ret_val && symbolic_manager != nullptr && !trav_genotype.empty() &&
        parent_child_trav_sets == nullptr) {
        const Snarl* managed_ptr = snarl_manager.into_which_snarl(snarl.start().node_id(),
                                                                  snarl.start().backward());
        if (managed_ptr != nullptr) {

            // The child-independent parts of the exactly-once test, built once for this snarl.
            const ChainInlineContext inline_ctx =
                build_chain_inline_context(snarl, travs, trav_genotype, ref_trav_idx);
            // Also once for this snarl: see TraversalNodeIndex.
            vector<TraversalNodeIndex> trav_visits;
            trav_visits.reserve(travs.size());
            for (const SnarlTraversal& t : travs) {
                trav_visits.push_back(index_traversal_nodes(t));
            }
            for (const Snarl* child : snarl_manager.children_of(managed_ptr)) {
                if (child == nullptr || snarl_manager.is_trivial(child, graph)) {
                    continue;
                }
                // A chain that no reference path passes through has no REF or POS for its records,
                // so it is skipped unless off-reference descent is on.
                bool child_off_reference = false;
                if (ref_trav_idx >= 0 && ref_trav_idx < (int)travs.size()) {
                    vector<int> ref_only(1, ref_trav_idx);
                    if (child_ploidy(trav_visits, ref_only, *child, 1) == 0) {
                        // With off-reference descent, such a chain is genotyped and recorded but has
                        // no line.
                        if (!off_reference_nesting) {
                            ++descent_counters.skipped_no_ref;
                            continue;
                        }
                        child_off_reference = true;
                        ++descent_counters.off_reference;
                    }
                }
                // Inherited: everything under a chain the reference does not cross is also off it.
                if (nested_context.no_reference) {
                    child_off_reference = true;
                }

                // The exactly-once test: under block emission, a chain that every called strand
                // crosses only inside a difference block is already spelled by that block's ALT. It
                // holds back the chain's line, not its descent, so the chain is still genotyped,
                // recorded and phased. Inherited by chains inside it. Does nothing when block
                // emission is off, or for a snarl whose projection has no symbols.
                bool child_reported_inline =
                    nested_context.reported_inline
                    || chain_reported_inline(inline_ctx, *child);

                int copies = child_ploidy(trav_visits, trav_genotype, *child, ploidy);
                bool retain_only = nested_context.retain_only;
                if (copies <= 0) {
                    // No called allele reaches it yet. Visited anyway, while this window's reads are
                    // in memory, since the linkage model may move the parent onto an allele that
                    // does reach it. Nothing about it is written unless the linkage pass says so.
                    ++descent_counters.skipped_no_copy;
                    if (!stage_records || linkage_collector == nullptr) {
                        // Without retention there is nothing to come back to. Without the linkage
                        // model nothing moves the parent after the direct pass, so the sample has no copy
                        // of this chain; the linkage pass, which has no chosen parent to read, would
                        // otherwise render it at the parent's ploidy.
                        continue;
                    }
                    retain_only = true;
                }

                // Saved and restored, since a child may descend further, and its own children must see
                // it as their parent.
                NestedContext saved = nested_context;
                nested_context.one_copy = (copies == 1);
                nested_context.parent_record_key = record_key_of(snarl);
                nested_context.retain_only = retain_only;
                nested_context.no_reference = child_off_reference;
                // Where this child starts along the first called allele that reaches it, added to
                // the offset of its parent. Only an off-reference chain uses it, but it is computed
                // for every chain, so that offsets add up down the tree.
                nested_context.parent_offset =
                    saved.parent_offset + offset_along_genotype(travs, trav_genotype, *child);
                nested_context.reported_inline = child_reported_inline;
                // The chain's identity, from its boundary nodes.
                {
                    const pair<nid_t, nid_t> cb = chain_bounds_of(child, snarl_manager);
                    nested_context.chain_key =
                        (size_t)((uint64_t)cb.first * 1000003ULL) ^ (size_t)(uint64_t)cb.second;
                }
                bool crossing_known = true;   // child_crossing_mask always sets it
                // The mask is over this snarl's own candidate traversals, which exist whether or not
                // a line was written.
                nested_context.parent_crossing =
                    child_crossing_mask(trav_visits, *child, &crossing_known);
                nested_context.crossing_known = crossing_known;
                size_t saved_level = current_level;
                current_level = saved_level + 1;
                ++g_descent_depth;
                if (g_descent_depth < 16) {
                    ++descent_counters.depth_hist[g_descent_depth];
                }
                // `copies` is zero only for a chain no called parent allele reaches, which is still
                // genotyped; it then takes the parent's ploidy, the most copies a child can have.
                // The other ploidy's answer is computed as well, so the linkage pass can change it
                // later.
                call_snarl_internal(*child, ref_path_name,
                                    make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                    nullptr, copies >= 1 ? copies : ploidy);
                --g_descent_depth;
                current_level = saved_level;
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

    // Descent above and the --top-down recursion, which builds each child's ChildTraversalSets
    // from `travs`, are done, so the staged site can take the traversals. At most one of these
    // is set.
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
                                   ref_offset_of(ref_offsets, ref_path_name),
                                   ref_ploidy_of(ref_ploidies, ref_path_name));
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
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        // and overlay
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
    // The interval the ploidy above was decided over. emit_snarl_recursive derives the ploidy
    // again from it, and indexes genotype_by_ploidy by that ploidy.
    record.ref_path_interval = make_pair((int64_t)gt_ref_interval.first,
                                         (int64_t)gt_ref_interval.second);
    record.ref_trav_idx = ref_trav_idx;

    // in the snarl graph, snarls a represented by a snarl end point and that's it.  here we fix up
    // the traversals to actually embed the snarls todo: should be able to avoid copy here!
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
            // indexed by it, so the record's own interval is used.
            ploidy = ploidy_at(record.ref_path_name, record.ref_path_interval.first,
                               ref_offset_of(ref_offsets, record.ref_path_name),
                               ref_ploidy_of(ref_ploidies, record.ref_path_name));
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
                                 ref_offset_of(ref_offsets, record.ref_path_name), genotype_snarls, ploidy, trav_to_flat_string);
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
                    // there's no call here. but we really want to emit something, so try picking
                    // the reference or first allele
                    if (record.ref_trav_idx >= 0) {
                        fallback_allele = record.ref_trav_idx;
                    } else if (!record.travs.empty()) {
                        fallback_allele = 0;
                    }
                }
                if (fallback_allele >= (int)record.travs.size()) {
                    flat_allele += "<...>";
                } else {
                    // todo: passing in a single ploidy simplisitic, would need to derive from the
                    // calls when reucrising in practice, the results will nearly the same but still
                    // needs fixing we try to get the allele from the genotype if possible, but
                    // fallback on the fallback_allele
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

