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
    cerr << "[vg call] descent skipped: " << g_descent_skipped_no_copy.load()
         << " children no called allele reaches, " << g_descent_skipped_no_ref.load()
         << " with no reference path through them" << endl;
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
    if (symbolic_manager != nullptr) {
        ss << "##INFO=<ID=NGT2,Number=1,Type=String,Description=\"For a nested site called at "
           << "ploidy 1, the genotype its own reads would give at ploidy 2, in this record's allele "
           << "numbering. Present only on nested haploid records. '.' means that genotype uses a "
           << "traversal this record does not carry, so reaching it would need a new ALT and not "
           << "just a different GT. Reported, not acted on: see FILTER=nested_diploid\">" << endl;
        // Nested descent chose each child's ploidy from its parent's pre-linkage genotype, and
        // linkage can then rewrite the parent. Where that changes which parent alleles cross the
        // child, the child record outlives the genotype that justified it. Flagged rather than
        // dropped: the read evidence at the child is real, and it is evidence against the parent's
        // new genotype as much as the parent is evidence against the child.
        ss << "##FILTER=<ID=nested_diploid,Description=\"Called at ploidy 1 because one parent "
           << "allele crossed this child chain, but after linkage both parent haplotypes cross it. "
           << "The locus is diploid here and this genotype names only one allele\">" << endl;
        ss << "##FILTER=<ID=nested_haploid,Description=\"Called at ploidy 2 because both parent "
           << "alleles crossed this child chain, but after linkage only one parent haplotype does. "
           << "The sample carries one copy of the chain and this genotype names two\">" << endl;
        ss << "##FILTER=<ID=nested_unreachable,Description=\"Called because a parent allele crossed "
           << "this child chain, but after linkage neither parent haplotype does. The sample carries "
           << "no copy of the chain under its own parent record, so this genotype has no haplotype "
           << "to sit on\">" << endl;
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
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
    if (gbwt != nullptr) {
        // One per thread, built here so the parallel region never allocates one.
        this->linkage_gbwt_cache.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            this->linkage_gbwt_cache.emplace_back(*gbwt);
        }
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

bool VCFOutputCaller::add_variant(vcflib::Variant& var) const {
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
    output_variants[omp_get_thread_num()].push_back(make_pair(make_pair(var.sequenceName, var.position), dest));
    return true;
}

void VCFOutputCaller::resolve_linkage() {
    if (!linkage_resolved) {
        resolve_linkage_generation(0, true);
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
    vector<LinkageCollector::NestedIncoherence> incoherent;
    // `linkage_phased` accumulates across generations rather than being replaced. The model needs
    // the earlier generations back: a nested site's strand is read off its parent's PhaseCall, and
    // a clamped site's phase is pinned to the pair already emitted for it.
    vector<LinkageCollector::Change> resolved =
        linkage_collector->resolve_generation(generation, last,
                                              emit_phasing ? &linkage_phased : nullptr,
                                              emit_phasing ? &incoherent : nullptr);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    linkage_seconds += seconds;
    linkage_changed += resolved.size();
    for (const LinkageCollector::Change& c : resolved) {
        linkage_changes[make_pair(c.contig, c.position)] = c;
    }
    for (const LinkageCollector::NestedIncoherence& ni : incoherent) {
        linkage_nested_filters[ni.record_key] = ni.kind;
    }
    if (!last) {
        // One line per intermediate generation, so a deferred-descent run shows its own shape:
        // how many sites each barrier settled and what it cost.
        cerr << "[vg call] linkage generation " << generation << ": "
             << linkage_collector->num_sites_at(generation) << " sites, "
             << resolved.size() << " genotypes changed, " << seconds << " s" << endl;
        return;
    }

    size_t unexplained = 0;
    size_t order_arbitrary = 0;
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        linkage_phasings[make_pair(pc.contig, pc.position)] = pc;
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
         << linkage_changed << " genotypes changed, " << linkage_seconds << " s" << endl;
    if (emit_phasing) {
        // The wildcard count is the honest caveat on a chromosome-length phase block: at
        // those sites the panel does not name a strand, so the phase either side of them
        // rests on the transition model alone.
        cerr << "[vg call] phasing: " << linkage_phased.size() << " sites phased, "
             << unexplained << " with a strand the panel does not explain" << endl;
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
        write_mosaic(linkage_phased);
    }
    if (!linkage_nested_filters.empty()) {
        size_t want_diploid = 0, want_haploid = 0, unreachable = 0;
        for (const auto& kv : linkage_nested_filters) {
            switch (kv.second) {
                case LinkageCollector::NestedIncoherence::WantsDiploid: ++want_diploid; break;
                case LinkageCollector::NestedIncoherence::WantsHaploid: ++want_haploid; break;
                default: ++unreachable; break;
            }
        }
        cerr << "[vg call] nested: " << linkage_nested_filters.size() << " records flagged, "
             << want_diploid << " nested_diploid, " << want_haploid << " nested_haploid, "
             << unreachable << " nested_unreachable" << endl;
    }
}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (include_nested) {
        update_nesting_info_tags(snarl_manager);
    }
    vector<pair<pair<string, size_t>, string>> all_variants;
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
    std::sort(all_variants.begin(), all_variants.end(), [](const pair<pair<string, size_t>, string>& v1,
                                                           const pair<pair<string, size_t>, string>& v2) {
            return v1.first.first < v2.first.first || (v1.first.first == v2.first.first && v1.first.second < v2.first.second);
        });
    // Phase two of the linkage pass. The records are already all here, compressed, with
    // (contig, position) uncompressed as the sort key -- so a change can be matched without ever
    // having kept the record itself, and only the records that actually move are re-parsed.
    resolve_linkage();

    for (const auto& v : all_variants) {
        string dest;
        int ret = zstdutil::DecompressString(v.second, dest);
        assert(ret == 0);
        if (!linkage_changes.empty()) {
            auto found = linkage_changes.find(v.first);
            if (found != linkage_changes.end()) {
                apply_linkage_change(dest, found->second);
            }
        }
        if (!linkage_phasings.empty()) {
            // After the change, never before: phasing has to describe the genotype that is
            // actually written out.
            auto found = linkage_phasings.find(v.first);
            if (found != linkage_phasings.end()) {
                apply_phasing(dest, found->second);
            }
        }
        if (!linkage_nested_filters.empty()) {
            // The record key is the hash of the ID column, which is how the linkage layer keyed it,
            // so the identity is recoverable from the line itself and nothing extra has to be
            // carried through the compressed buffer.
            size_t a = dest.find('\t');
            size_t b = a == string::npos ? string::npos : dest.find('\t', a + 1);
            size_t c = b == string::npos ? string::npos : dest.find('\t', b + 1);
            if (c != string::npos) {
                auto found = linkage_nested_filters.find(
                    std::hash<string>{}(dest.substr(b + 1, c - b - 1)));
                if (found != linkage_nested_filters.end()) {
                    apply_nested_filter(dest, found->second);
                }
            }
        }
        out_stream << dest << endl;
    }
}

void VCFOutputCaller::apply_linkage_change(string& line,
                                           const LinkageCollector::Change& change) const {
    // Edit the emitted line rather than rebuilding the record: everything else about it -- AD, DP,
    // GL, the traversals in AT -- is still the per-site truth and should not move.
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
        return;
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
        return;
    }
    size_t gt_field = keys.size();
    size_t gq_field = keys.size();
    size_t gqi_field = keys.size();
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "GT") {
            gt_field = i;
        } else if (keys[i] == "GQ") {
            gq_field = i;
        } else if (keys[i] == "GQI") {
            gqi_field = i;
        }
    }
    if (gt_field == keys.size()) {
        return;
    }

    // A haploid record's GT is a bare allele, not a pair, so its expected and replacement strings
    // have to be written that way. Building "i/j" regardless made the guard below reject every
    // haploid change: the linkage layer did the work, reported it in the progress line, and then
    // dropped all of it on the floor. chrY and non-pseudoautosomal chrX got no linkage correction
    // at all, silently -- and worse, the phasing and the mosaic are built from the *post*-linkage
    // genotypes, so the mosaic described genotypes the VCF did not contain.
    //
    // Ploidy is read off the record rather than off the change, because the record is what is
    // being patched. At ploidy 1 the collector sets both allele slots to the same value, so either
    // one is the allele.
    bool haploid_record = values[gt_field].find('/') == string::npos
                          && values[gt_field].find('|') == string::npos;

    // Guard. Two records can share a (contig, position) -- a nested site under -A, for one -- and
    // patching the wrong one would be silent. Only replace a genotype that is the one the per-site
    // model actually chose at the site this change came from.
    {
        string current = values[gt_field];
        std::replace(current.begin(), current.end(), '|', '/');
        if (haploid_record) {
            if (change.called_i != change.called_j
                || current != std::to_string(change.called_i)) {
                return;
            }
        } else {
            string expected = std::to_string(change.called_i) + "/"
                              + std::to_string(change.called_j);
            string expected_alt = std::to_string(change.called_j) + "/"
                                  + std::to_string(change.called_i);
            if (current != expected && current != expected_alt) {
                return;
            }
        }
    }

    values[gt_field] = haploid_record
                           ? std::to_string(change.allele_i)
                           : std::to_string(change.allele_i) + "/"
                                 + std::to_string(change.allele_j);
    if (gq_field != keys.size()) {
        // GQ becomes the phred-scaled complement of the posterior, which is what the quality means
        // once a posterior exists -- then discounted by the explained-read share, exactly as the
        // per-site GQ was. GQI is left alone: it is the per-site likelihood ratio, and keeping it
        // untouched is what makes the before/after pair comparable.
        //
        // Two corrections to the raw posterior, both measured over four datasets on AUC and on
        // false calls surviving at matched recall.
        //
        // The share discount applies here for the same reason it applies to the per-site GQ:
        // reads whose best allele lies outside the call enter every genotype's likelihood and
        // cancel, so neither the likelihood ratio nor the HMM posterior -- built from that same
        // emission -- can see them.
        //
        // The cap at GQI is the larger effect, and it is not a consistency nicety. The posterior
        // is computed under a strong prior (the panel frequency exponent defaults to 5), so
        // `1 - posterior` understates uncertainty -- and it does so exactly where the per-site
        // evidence was weakest, because that is where the linkage layer acts at all. Capping
        // re-anchors the reported quality to the read evidence: the linkage pass may lower
        // confidence and may not raise it above what the reads alone supported. Worth about
        // +0.003 of AUC and 1-2% fewer surviving false calls, against +0.0001 to +0.0009 for the
        // share discount alone.
        //
        // It also makes `GQ <= GQI` hold on every record. The share discount alone does not,
        // which was an early misreading: the posterior-based quality is not derived from GQI and
        // has no arithmetic relation to it, so scaling it by a factor in [0,1] cannot bound it.
        double q = change.posterior >= 1.0
                       ? 256.0 : -10.0 * log10(max(1.0 - change.posterior, 1e-26));
        q *= change.explained_share;
        if (gqi_field != keys.size()) {
            try {
                q = min(q, stod(values[gqi_field]));
            } catch (const std::exception&) {
                // GQI absent or unparsable: leave the discounted posterior uncapped rather
                // than dropping the change.
            }
        }
        values[gq_field] = std::to_string((int)min(256.0, max(0.0, round(q))));
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
    out << "#mosaic-version\t2\n";
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
        << "explain that strand there.\n";
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
    std::function<void(size_t, size_t, int, size_t)> emit_span =
        [&](size_t from, size_t to, int strand, size_t hap) {
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
                out << "*\t*";
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
        emit_span(hi, to, strand, hap);
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
                size_t prev = strand == 0 ? phasing[seg_start].hap_first
                                          : phasing[seg_start].hap_second;
                bool last = (t + 1 == j);
                bool changes = !last
                               && (strand == 0 ? phasing[t + 1].hap_first
                                               : phasing[t + 1].hap_second) != hap;
                (void)prev;
                if (last || changes) {
                    // One run of a single haplotype, but possibly several GBWT fragments of it.
                    // A segment must be walkable from one position, so the run is cut wherever the
                    // fragment under it changes: emit_span walks the sites, re-resolving whenever
                    // the position stops belonging to the same fragment, and emits one row per
                    // fragment. Without this a consumer following LF() from the segment's position
                    // would hit the endmarker partway and have no way to pick up the rest.
                    emit_span(seg_start, t, strand, hap);
                    seg_start = t + 1;
                }
            }
        }
        i = j;
    }
    cerr << "[vg call] mosaic: " << total_segments << " segments over " << phasing.size()
         << " sites, written to " << mosaic_path << endl;
}

void VCFOutputCaller::apply_nested_filter(string& line,
                                          LinkageCollector::NestedIncoherence::Kind kind) const {
    // Column 7 of the eight fixed fields. Edited in place for the same reason apply_linkage_change
    // is: re-parsing the record to a vcflib::Variant and reprinting it would reformat every float
    // in the line, which makes an unrelated diff on 90% of records.
    size_t start = 0;
    for (int field = 0; field < 6; ++field) {
        start = line.find('\t', start);
        if (start == string::npos) {
            return;
        }
        ++start;
    }
    size_t end = line.find('\t', start);
    if (end == string::npos) {
        return;
    }
    string tag;
    switch (kind) {
        case LinkageCollector::NestedIncoherence::WantsDiploid: tag = "nested_diploid"; break;
        case LinkageCollector::NestedIncoherence::WantsHaploid: tag = "nested_haploid"; break;
        default: tag = "nested_unreachable"; break;
    }
    // "." and "PASS" both mean "nothing to say", so the tag replaces them rather than joining them:
    // a FILTER of "PASS;nested_diploid" is a contradiction.
    string current = line.substr(start, end - start);
    if (current == "." || current == "PASS" || current.empty()) {
        line.replace(start, end - start, tag);
    } else if (current.find(tag) == string::npos) {
        line.replace(start, end - start, current + ";" + tag);
    }
}

void VCFOutputCaller::apply_phasing(string& line,
                                    const LinkageCollector::PhaseCall& phase) const {
    // Same line-editing approach as apply_linkage_change, and for the same reason: every other
    // field remains the per-site truth. Only GT's separator and order change, and PS is appended.
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
        return;
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
        return;
    }
    size_t gt_field = keys.size();
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "GT") {
            gt_field = i;
        }
    }
    if (gt_field == keys.size()) {
        return;
    }

    // The phased genotype must be a permutation of the unphased one. Two records can share a
    // position, and phasing the wrong one -- or phasing a genotype the record does not carry --
    // would silently emit a call that disagrees with its own likelihoods. Checked rather than
    // assumed, because the failure is invisible in the output.
    {
        string current = values[gt_field];
        std::replace(current.begin(), current.end(), '|', '/');
        size_t slash = current.find('/');
        if (slash == string::npos) {
            // Haploid: one allele, and it must be the one the path spells.
            if (phase.ploidy == 1 && current == std::to_string(phase.allele_first)) {
                slash = string::npos;   // fall through to the write below
            } else {
                return;
            }
        } else {
        string lhs = current.substr(0, slash);
        string rhs = current.substr(slash + 1);
        string want_a = std::to_string(phase.allele_first);
        string want_b = std::to_string(phase.allele_second);
        bool same = (lhs == want_a && rhs == want_b) || (lhs == want_b && rhs == want_a);
        if (!same) {
            return;
        }
        }
    }

    if (phase.ploidy == 1 && phase.nested_strand >= 0) {
        // A nested haploid site is one strand of a diploid locus, not a haploid locus: it is called
        // at ploidy 1 because the parent's *other* allele deletes the chain, so there is no sequence
        // on that strand to genotype. Written as a phased pair with the empty strand as ".", which
        // is the only place the VCF can carry which strand the allele is on -- a bare "a" names none,
        // and that is why the strand lived only in the mosaic and no phasing tool could read it.
        //
        // "." rather than "0": the other haplotype does not carry the reference sequence here, it
        // carries nothing. "*" would say the same thing and say it in ALT, but adding an allele
        // changes the arity of AD, GL and GQI on the record, and the strand is not known until
        // phasing, long after those were written.
        values[gt_field] = phase.nested_strand == 0
                           ? std::to_string(phase.allele_first) + "|."
                           : "." + ("|" + std::to_string(phase.allele_first));
    } else if (phase.ploidy == 1) {
        // A genuinely haploid locus -- a haploid contig, or a nested site with no strand to name.
        // One allele, no phase; only PS below is meaningful, and only as a block label. Writing
        // "a|a" would claim a homozygous diploid call.
        values[gt_field] = std::to_string(phase.allele_first);
    } else {
        values[gt_field] = std::to_string(phase.allele_first) + "|"
                           + std::to_string(phase.allele_second);
    }

    // PS identifies the phase block; phase is only comparable within one. Replaces any existing
    // value rather than appending a second copy, since write_variants can be called on a buffer
    // that has already been phased in a previous pass.
    size_t ps_field = keys.size();
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "PS") {
            ps_field = i;
        }
    }
    if (ps_field == keys.size()) {
        keys.push_back("PS");
        values.push_back(std::to_string(phase.phase_set));
    } else {
        values[ps_field] = std::to_string(phase.phase_set);
    }

    string fmt_keys, fmt_values;
    for (size_t i = 0; i < keys.size(); ++i) {
        if (i) {
            fmt_keys += ":";
            fmt_values += ":";
        }
        fmt_keys += keys[i];
        fmt_values += values[i];
    }
    fields[8] = fmt_keys;
    fields[9] = fmt_values;
    line.clear();
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) {
            line += "\t";
        }
        line += fields[i];
    }
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
                                            vcflib::Variant& out_variant) const {
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

    // GL is Number=G.  vg emits it i-major -- "for i; for j = i..n" -- which is not the VCF spec's
    // ordering for 3+ alleles, but the fold has to match what is actually written.  Take the max
    // over the old genotype classes mapping onto each new one: that is the max-marginal, i.e. the
    // merged allele scores as whichever of its members fit best.
    auto gl_it = sample.find("GL");
    if (gl_it != sample.end()) {
        size_t n_old = merge_to.size();
        auto gl_index = [](size_t i, size_t j, size_t n) { return i * n - (i * (i - 1)) / 2 + (j - i); };
        // The diploid layout is the only one that can occur: merging needs at least two distinct
        // called ALTs, site_genotype has one entry per ploidy, and PoissonSupportSnarlCaller --
        // whose update_vcf_info is the only writer of GL -- asserts in genotype that ploidy is
        // 1 or 2.  So n_old is always 3 here.
        assert(gl_it->second.size() == n_old * (n_old + 1) / 2);
        bool gl_usable = true;
        vector<double> folded((size_t)n_new * (n_new + 1) / 2,
                              -std::numeric_limits<double>::infinity());
        for (size_t i = 0; i < n_old && gl_usable; ++i) {
            for (size_t j = i; j < n_old && gl_usable; ++j) {
                double v = 0;
                gl_usable = parse_vcf_double(gl_it->second[gl_index(i, j, n_old)], v);
                if (!gl_usable) {
                    break;
                }
                size_t ni = new_index[i], nj = new_index[j];
                if (ni > nj) {
                    std::swap(ni, nj);
                }
                double& slot = folded[gl_index(ni, nj, (size_t)n_new)];
                slot = std::max(slot, v);
            }
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
            contigs.insert(output_variant_record.first.first);
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
thread_local VCFOutputCaller::EmittedAlleles VCFOutputCaller::last_emitted;

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
    last_emitted.recorded = false;

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

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name, out_variant);

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
    merge_similar_alleles(graph, site_traversals, site_genotype, sample_name, out_variant);
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
    }

    // What this nested haploid site would call at ploidy 2.
    //
    // Its ploidy came from the parent's genotype before linkage, and linkage can move the parent so
    // that both haplotypes cross the chain -- at which point the record names one allele where the
    // locus has two. Re-genotyping it needs the diploid likelihoods, so this reports what they say,
    // in the record's own allele numbering. "." means the diploid answer uses a traversal this
    // record does not carry, which is the case a genotype rewrite could not fix on its own.
    if (nested_context.active) {
        const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* rl_alt =
            dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
        if (rl_alt != nullptr && rl_alt->alt_ploidy_best.size() == 2) {
            string spelled;
            for (int t : rl_alt->alt_ploidy_best) {
                auto it = trav_to_allele.find(t);
                if (it == trav_to_allele.end()) {
                    spelled.clear();
                    break;
                }
                if (!spelled.empty()) {
                    spelled += "/";
                }
                spelled += std::to_string(it->second);
            }
            out_variant.info["NGT2"].push_back(spelled.empty() ? "." : spelled);
        }
    }

    if (genotype_snarls || !out_variant.alt.empty()) {
        bool added = add_variant(out_variant);
        if (added && linkage_collector != nullptr) {
            // Phase one of the linkage pass: keep the compact site, not the record. See
            // LinkageCollector for why the CallInfo cannot be retained instead.
            const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* rl_info =
                dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                    call_info.get());
            // Haploid chains are fed too. They used to be dropped here by a `size() == 2` guard,
            // which silently cost chrY and non-pseudoautosomal chrX both the linkage layer and any
            // mosaic -- about 5% of a genome, and the part where the mosaic is the whole answer.
            size_t site_ploidy = site_genotype.size();
            bool alleles_ok = (site_ploidy == 1 || site_ploidy == 2);
            for (size_t k = 0; k < site_genotype.size() && alleles_ok; ++k) {
                alleles_ok = site_genotype[k] >= 0;
            }
            if (rl_info != nullptr && alleles_ok) {
                size_t n_alleles = out_variant.alleles.size();
                size_t n_gt = site_ploidy == 1 ? n_alleles : n_alleles * (n_alleles + 1) / 2;
                vector<double> gls(n_gt, -numeric_limits<double>::infinity());
                for (const auto& entry : rl_info->genotype_lls) {
                    if (entry.first.size() != site_ploidy) {
                        continue;
                    }
                    if (site_ploidy == 1) {
                        // A haploid genotype is one allele, so the likelihood vector is indexed by
                        // allele rather than by the triangular pair index.
                        auto it = trav_to_allele.find(entry.first[0]);
                        if (it == trav_to_allele.end() || (size_t)it->second >= n_alleles) {
                            continue;
                        }
                        gls[(size_t)it->second] = entry.second;
                        continue;
                    }
                    // genotype_lls is keyed by traversal index. Using those as VCF allele indices
                    // fed the model a scrambled emission -- it changed 34% of confident genotypes
                    // against a 0.1% budget, which is what a wrong emission looks like.
                    auto a_it = trav_to_allele.find(entry.first[0]);
                    auto b_it = trav_to_allele.find(entry.first[1]);
                    if (a_it == trav_to_allele.end() || b_it == trav_to_allele.end()) {
                        continue;   // a traversal that never became a VCF allele
                    }
                    size_t i = (size_t)a_it->second, j = (size_t)b_it->second;
                    if (i >= n_alleles || j >= n_alleles) {
                        continue;
                    }
                    gls[LinkageModel::genotype_index(i, j)] = entry.second;
                }
                size_t gi = (size_t)site_genotype[0];
                size_t gj = site_ploidy == 2 ? (size_t)site_genotype[1] : gi;
                if (site_ploidy == 2 && gi > gj) {
                    std::swap(gi, gj);
                }
                linkage_collector->record(out_variant.sequenceName, out_variant.position,
                                          n_alleles, gls,
                                          // site_traversals, not called_traversals: the VCF
                                          // alleles come from the deduplicated list, so a
                                          // traversal index is not an allele index. Passing the
                                          // wrong one wrote past the end of the genotype vector
                                          // and corrupted the heap.
                                          panel_alleles(graph, site_traversals),
                                          gi, gj,
                                          // A key intrinsic to the site, so ordering cannot depend
                                          // on which thread got there first.
                                          std::hash<string>{}(out_variant.id),
                                          // So a rewritten GQ can carry the same discount the
                                          // per-site GQ did. See apply_linkage_change.
                                          rl_info->explained_share, site_ploidy,
                                          // Snarl boundaries, so the mosaic output can anchor on
                                          // node IDs rather than on reference positions.
                                          (int64_t)snarl.start().node_id(),
                                          (int64_t)snarl.end().node_id(),
                                          // Nested provenance, so the linkage pass can keep this
                                          // site out of the diploid chain runs, place it on the
                                          // parent's strand afterwards, and check its ploidy against
                                          // the genotype the parent finally ends up with.
                                          nested_context.active,
                                          nested_context.parent_record_key,
                                          nested_context.parent_slot,
                                          nested_context.parent_crossing,
                                          // Which resolve pass settles this site. Descent reads it
                                          // back off `last_emitted` to decide whether this snarl's
                                          // children can be visited now or must wait for linkage.
                                          current_generation);
                last_emitted.recorded = true;
            }
        }
        if (!added) {
            stringstream ss;
            ss << out_variant;
            cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":" << out_variant.position
                 << " with ID=" << out_variant.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
                 << VCFOutputCaller::max_vcf_line_length << endl;
        }
        return added;
    }
    return true;
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

uint64_t FlowCaller::child_crossing_mask(const vector<SnarlTraversal>& travs,
                                         const map<int, int>& trav_to_allele,
                                         const Snarl& child) {
    uint64_t mask = 0;
    for (const auto& kv : trav_to_allele) {
        if (kv.second < 0 || kv.second >= 64) {
            return 0;   // unknown, not none
        }
        if (kv.first < 0 || kv.first >= (int)travs.size()) {
            continue;
        }
        if (crossings_of_child(travs[kv.first], child) > 0) {
            mask |= (uint64_t)1 << kv.second;
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
    if (capped && show_progress) {
#pragma omp critical (cerr)
        cerr << "[vg call] --nested: an allele crosses child snarl " << start << "-" << end
             << " more than once; counted as one copy" << endl;
    }
    return min(copies, cap);
}

void FlowCaller::set_defer_nested_descent(bool defer) {
    this->defer_nested_descent = defer;
    if (defer) {
        // Sized once, here, rather than lazily inside the parallel region that writes it.
        size_t threads = max((size_t)get_thread_count(), (size_t)omp_get_max_threads());
        pending_descents.assign(max(threads, (size_t)1), {});
    }
}

size_t FlowCaller::pending_descent_count() const {
    size_t n = 0;
    for (const auto& queue : pending_descents) {
        n += queue.size();
    }
    return n;
}

size_t FlowCaller::descend_pending(size_t generation) {
    vector<PendingDescent> work;
    work.reserve(pending_descent_count());
    for (auto& queue : pending_descents) {
        std::move(queue.begin(), queue.end(), std::back_inserter(work));
        queue.clear();
    }
    if (work.empty()) {
        return 0;
    }

    // The parents' settled allele pairs, by record key. This is the phasing the VCF will carry, so
    // descent and the emitted parent record cannot disagree about which alleles the sample has.
    unordered_map<size_t, const LinkageCollector::PhaseCall*> settled;
    settled.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        settled[pc.record_key] = &pc;
    }

    // Node-ID order, so a read source that fetches by window still sees each window once in this
    // pass. Without it the pass is random access over the contig and every window is re-fetched.
    std::sort(work.begin(), work.end(), [](const PendingDescent& a, const PendingDescent& b) {
        return a.child->start().node_id() < b.child->start().node_id();
    });

    current_generation = generation;
    std::atomic<size_t> called(0);
    std::atomic<size_t> no_parent(0), unknown_mask(0), not_carried(0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)work.size(); ++i) {
        const PendingDescent& pd = work[i];
        auto found = settled.find(pd.parent_record_key);
        if (found == settled.end()) {
            ++no_parent;   // no phased parent record: nothing to take a ploidy from
            continue;
        }
        if (pd.crossing_mask == 0) {
            ++unknown_mask;   // descent could not express it, which is unknown rather than none
            continue;
        }
        const LinkageCollector::PhaseCall& parent = *found->second;
        bool first = parent.allele_first < 64
                     && ((pd.crossing_mask >> parent.allele_first) & 1);
        bool second = parent.ploidy == 2 && parent.allele_second < 64
                      && ((pd.crossing_mask >> parent.allele_second) & 1);
        int copies = min((int)first + (int)second, pd.parent_ploidy);
        if (copies <= 0) {
            ++not_carried;   // the settled genotype does not carry the chain, so no call belongs here
            continue;
        }
        nested_context.active = (copies == 1);
        nested_context.parent_record_key = pd.parent_record_key;
        // Which of the parent's strands carries the chain. On this side of the barrier it has to come
        // from the settled allele pair: the traversal-order slot the inline path records does not
        // exist here, the traversals having been dropped. The two conventions were measured against
        // the phased truth and are indistinguishable -- 1,655 switches against 1,661 on chr20 -- so
        // this is a forced choice that costs nothing.
        nested_context.parent_slot = first ? 0 : 1;
        nested_context.parent_crossing = pd.crossing_mask;
        call_snarl_internal(*pd.child, pd.ref_path_name,
                            make_pair((size_t)pd.ref_start, (size_t)pd.ref_end), nullptr, copies);
        nested_context = NestedContext();
        ++called;
    }
    if (show_progress) {
        cerr << "[vg call] deferred descent: generation " << generation << ", " << work.size()
             << " queued, " << called.load() << " called, " << not_carried.load()
             << " not carried by the settled parent genotype";
        if (no_parent.load() > 0 || unknown_mask.load() > 0) {
            cerr << ", " << no_parent.load() << " with no phased parent, " << unknown_mask.load()
                 << " with an unreadable mask";
        }
        cerr << endl;
    }
    return called.load();
}

void FlowCaller::run_deferred_descent() {
    if (!defer_nested_descent) {
        return;
    }
    size_t generation = 0;
    while (pending_descent_count() > 0) {
        resolve_linkage_generation(generation, false);
        ++generation;
        if (descend_pending(generation) == 0) {
            // Nothing became callable, so this generation has no sites and the queue those calls
            // would have filled is empty. The final pass below still owes the reporting.
            break;
        }
    }
    resolve_linkage_generation(generation, true);
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
            added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx, trav_call_info, ref_path_name,
                                 ref_offsets[ref_path_name], genotype_snarls, ploidy);
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
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
            added = emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx, trav_call_info, ref_path_name,
                                 ref_offsets[ref_path_name], genotype_snarls, ploidy);
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
            const EmittedAlleles parent_alleles = last_emitted;
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
            const bool defer = defer_nested_descent && emit_phasing && parent_alleles.recorded;
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
                if (defer) {
                    // Queued for every child the reference reaches, including the ones no called
                    // allele reaches. That is the point rather than an oversight: whether an allele
                    // crosses the chain is a question about the *settled* genotype, and asking it
                    // here rather than now is what recovers the 296 children per chromosome that
                    // linkage makes reachable and inline descent silently drops.
                    PendingDescent pd;
                    pd.child = child;
                    pd.ref_path_name = ref_path_name;
                    pd.ref_start = (int64_t)get<0>(ref_interval);
                    pd.ref_end = (int64_t)get<1>(ref_interval);
                    pd.parent_record_key = std::hash<string>{}(print_snarl(snarl, false));
                    pd.crossing_mask = parent_alleles.valid
                        ? child_crossing_mask(travs, parent_alleles.trav_to_allele, *child)
                        : 0;
                    pd.parent_ploidy = ploidy;
                    pending_descents[omp_get_thread_num()].push_back(std::move(pd));
                    continue;
                }
                int copies = child_ploidy(travs, trav_genotype, *child, ploidy);
                if (copies <= 0) {
                    // A child the *pre-linkage* genotype does not reach is dropped here and never
                    // looked at again, so nothing has counted how many of them the parent's final
                    // genotype does reach. Hand the crossing mask to the linkage pass, which is the
                    // only place that knows the final genotype, and let it answer.
                    ++g_descent_skipped_no_copy;
                    if (linkage_collector != nullptr && parent_alleles.valid) {
                        linkage_collector->record_skipped_child(
                            std::hash<string>{}(print_snarl(snarl, false)),
                            child_crossing_mask(travs, parent_alleles.trav_to_allele, *child));
                    }
                    continue;   // no called allele reaches it: a star allele, and nothing to call
                }
                // A haploid child hangs off exactly one parent slot, and which one decides the
                // strand its allele sits on once the parent is phased. Find it: the slot whose
                // called allele crosses the child.
                size_t slot = 0;
                if (copies == 1) {
                    for (size_t g = 0; g < trav_genotype.size(); ++g) {
                        vector<int> one(1, trav_genotype[g]);
                        if (child_ploidy(travs, one, *child, 1) == 1) {
                            slot = g;
                            break;
                        }
                    }
                }
                // Saved and restored rather than assigned: a child may descend further, and its own
                // children must see *it* as their parent, not this snarl.
                NestedContext saved = nested_context;
                nested_context.active = (copies == 1);
                nested_context.parent_record_key = std::hash<string>{}(print_snarl(snarl, false));
                nested_context.parent_slot = slot;
                // Which parent *alleles* cross this child, as opposed to which called traversal.
                // Linkage needs this to tell whether the ploidy the child is about to be genotyped
                // at survives the parent's final genotype: linkage can move the parent afterwards,
                // and a child called haploid because one allele crossed it is wrong if both do.
                // For every descended child, not only the haploid ones. A child called at ploidy 2
                // can be invalidated the other way -- linkage moves the parent to a genotype where
                // only one allele crosses, and the record then claims two haplotypes at a locus the
                // sample carries once. Nothing detected that while the mask was recorded only for
                // copies == 1, which is half the question the mask exists to answer.
                nested_context.parent_crossing =
                    parent_alleles.valid
                    ? child_crossing_mask(travs, parent_alleles.trav_to_allele, *child)
                    : 0;
                ++g_descent_depth;
                if (g_descent_depth < 16) {
                    ++g_descent_depth_hist[g_descent_depth];
                }
                call_snarl_internal(*child, ref_path_name,
                                    make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                    nullptr, copies);
                --g_descent_depth;
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

