#include "traversal_support.hpp"
#include "genotypekit.hpp"

//#define debug

namespace vg {

TraversalSupportFinder::TraversalSupportFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager) :
    graph(graph),
    snarl_manager(snarl_manager) {
}

TraversalSupportFinder::~TraversalSupportFinder() {
    
}

int64_t TraversalSupportFinder::get_edge_length(const edge_t& edge, const unordered_map<id_t, size_t>& ref_offsets) const {
    int len = -1;
    // use our reference traversal to try to come up with a deletion length for our edge
    // idea: if our edge corresponds to a huge deltion, it should be weighted accordingly
    auto s_it = ref_offsets.find(graph.get_id(edge.first));
    auto e_it = ref_offsets.find(graph.get_id(edge.second));
    if (s_it != ref_offsets.end() && e_it != ref_offsets.end()) {
        size_t start_offset = s_it->second;
        if (!graph.get_is_reverse(edge.first)) {
            start_offset += graph.get_length(edge.first);
        }
        size_t end_offset = e_it->second;
        if (graph.get_is_reverse(edge.second)) {
            end_offset += graph.get_length(edge.second);
        }
        if (start_offset > end_offset) {
            std::swap(start_offset, end_offset);
        }
        len = end_offset - start_offset;
    }
    return std::max(len, 1);
}

tuple<Support, Support, int> TraversalSupportFinder::get_child_support(const Snarl& snarl) const {
    // port over old functionality from support caller
    // todo: do we need to flag nodes as covered like it does?
    pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&snarl, graph, true);
    Support child_max_support;
    Support child_total_support;
    size_t child_size = 0;
    for (id_t node_id : contents.first) {
        Support child_support = get_avg_node_support(node_id);
        child_max_support = support_max(child_max_support, child_support);
        child_size += graph.get_length(graph.get_handle(node_id));
        child_total_support += child_support;
    }
    Support child_avg_support = child_total_support / child_size;
    // we always use child_max like the old support_caller.
    // this is the only way to get top-down recursion to work in many cases
    // todo: fix to use bottom up, get get support from actual traversals
    // every time!! 
    return std::tie(child_max_support, child_max_support, child_size);
}


Support TraversalSupportFinder::get_traversal_support(const SnarlTraversal& traversal) const {
    return get_traversal_set_support({traversal}, {}, {}, false, {}, false).at(0);
}

vector<Support> TraversalSupportFinder::get_traversal_genotype_support(const vector<SnarlTraversal>& traversals,
                                                                       const vector<int>& genotype,
                                                                       const set<int>& other_trav_subset,
                                                                       int ref_trav_idx) {
    set<int> tgt_trav_set(genotype.begin(), genotype.end());
    vector<int> tgt_travs(tgt_trav_set.begin(), tgt_trav_set.end());
    // get the support of just the alleles in the genotype, evenly splitting shared stuff
    vector<Support> allele_support = get_traversal_set_support(traversals, tgt_travs, tgt_trav_set, false, {}, true, ref_trav_idx);
    // get the support of everythin else, subtracting genotype supports, and splitting mutual supports
    vector<Support> other_support = get_traversal_set_support(traversals, tgt_travs, other_trav_subset, false, allele_support, true, ref_trav_idx);
    // combine the above two vectors
    for (int allele : tgt_travs) {
        other_support[allele] = allele_support[allele];
    }
    return other_support;
}

vector<Support> TraversalSupportFinder::get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                                  const vector<int>& shared_travs,
                                                                  const set<int>& tgt_travs,
                                                                  bool exclusive_only,
                                                                  const vector<Support>& exclusive_count_support,
                                                                  bool mutual_shared,
                                                                  int ref_trav_idx) const {

    // mutual_shared only makes sense when everything is shared
    assert(!mutual_shared || shared_travs.size() == traversals.size() || shared_travs.size() == tgt_travs.size());
    // exclusive_count_support corresponds to traversals
    assert(exclusive_count_support.empty() || exclusive_count_support.size() == traversals.size());
    
    // pass 1: how many times have we seen a node or edge
    unordered_map<id_t, double> node_counts;
    unordered_map<edge_t, double> edge_counts;
    map<Snarl, double> child_counts;

    for (auto trav_idx : shared_travs) {
        const SnarlTraversal& trav = traversals[trav_idx];
        for (int i = 0; i < trav.visit_size(); ++i) {
            const Visit& visit = trav.visit(i);
            double value = exclusive_count_support.empty() ? 1. : support_val(exclusive_count_support[trav_idx]);
            if (visit.node_id() != 0) {
                // Count the node once
                if (node_counts.count(visit.node_id())) {
                    node_counts[visit.node_id()] += value;
                } else {
                    node_counts[visit.node_id()] = value;
                }
            } else {
                // Count the child once
                if (child_counts.count(visit.snarl())) {
                    child_counts[visit.snarl()] += value;
                } else {
                    child_counts[visit.snarl()] = value;
                }
            }
            // note: there is no edge between adjacent snarls as they overlap
            // on their endpoints. 
            if (i > 0 && (trav.visit(i - 1).node_id() != 0 || trav.visit(i).node_id() != 0)) {
                edge_t edge = to_edge(graph, trav.visit(i - 1), visit);
                // Count the edge once
                if (edge_counts.count(edge)) {
                    edge_counts[edge] += value;
                } else {
                    edge_counts[edge] = value;
                }
            }
        }
    }

    // pass 1.5: get index for looking up deletion edge lengths (so far we aren't dependent
    // on having anything but a path handle graph, so we index on the fly)
    unordered_map<id_t, size_t> ref_offsets;
    if (ref_trav_idx >= 0) {
        ref_offsets = get_ref_offsets(traversals[ref_trav_idx]);
    }

    // pass 2: get the supports
    // we compute the various combinations of min/avg node/trav supports as we don't know which
    // we will need until all the sizes are known
    Support max_support;
    max_support.set_forward(numeric_limits<int>::max());
    vector<Support> min_supports_min(traversals.size(), max_support); // use min node support
    vector<Support> min_supports_avg(traversals.size(), max_support); // use avg node support
    vector<bool> has_support(traversals.size(), false);
    vector<Support> tot_supports_min(traversals.size()); // weighted by lengths, using min node support
    vector<Support> tot_supports_avg(traversals.size()); // weighted by lengths, using avg node support
    vector<int> tot_sizes(traversals.size(), 0); // to compute average from to_supports;
    vector<int> tot_sizes_all(traversals.size(), 0); // as above, but includes excluded lengths
    int max_trav_size = 0; // size of longest traversal

    bool count_end_nodes = false; // toggle to include snarl ends

    auto update_support = [&] (int trav_idx, const Support& min_support,
                               const Support& avg_support, int length, double share_count) {
        // keep track of overall size of longest traversal
        tot_sizes_all[trav_idx] += length;
        max_trav_size = std::max(tot_sizes_all[trav_idx], max_trav_size);

        // apply the scaling
        double denom_add = mutual_shared ? 0 : 1;
        double scale_factor = (exclusive_only && share_count > 0) ? 0. : 1. / (denom_add + share_count);
        
        // when looking at exclusive support, we don't normalize by skipped lengths
        if (scale_factor != 0 || !exclusive_only || !exclusive_count_support.empty()) {
            has_support[trav_idx] = true;
            Support scaled_support_min;
            Support scaled_support_avg;
            if (!exclusive_count_support.empty()) {
                scaled_support_min = min_support * scale_factor;
                scaled_support_avg = avg_support * scale_factor;
            } else {
                if (support_val(min_support) > share_count) {
                    scaled_support_min.set_forward(support_val(min_support) - share_count);
                }
                if (support_val(avg_support) > share_count) {
                    scaled_support_avg.set_forward(support_val(avg_support) - share_count);
                }                
            }
            tot_supports_min[trav_idx] += scaled_support_min;
            tot_supports_avg[trav_idx] += scaled_support_avg * length;
            tot_sizes[trav_idx] += length;
            min_supports_min[trav_idx] = support_min(min_supports_min[trav_idx], scaled_support_min);
            min_supports_avg[trav_idx] = support_min(min_supports_avg[trav_idx], scaled_support_avg * length);
        }
    };

    for (int trav_idx = 0; trav_idx < traversals.size(); ++trav_idx) {
        // target_set filter here
        if (!tgt_travs.empty() && !tgt_travs.count(trav_idx)) {
            continue;
        }
        const SnarlTraversal& trav = traversals[trav_idx];
        for (int visit_idx = 0; visit_idx < trav.visit_size(); ++visit_idx) {
            const Visit& visit = trav.visit(visit_idx);
            Support min_support;
            Support avg_support;
            int64_t length;
            int share_count = 0;

            if (visit.node_id() != 0) {
                // get the node support
                min_support = get_min_node_support(visit.node_id());
                avg_support = get_avg_node_support(visit.node_id());
                length = graph.get_length(graph.get_handle(visit.node_id()));
                if (node_counts.count(visit.node_id())) {
                    share_count = node_counts[visit.node_id()];
                } 
            } else {
                // get the child support
                tie(min_support, avg_support, length) = get_child_support(visit.snarl());
                if (child_counts.count(visit.snarl())) {
                    share_count = child_counts[visit.snarl()];
                } 
            }
            if (count_end_nodes || (visit_idx > 0 && visit_idx < trav.visit_size() - 1)) {
                update_support(trav_idx, min_support, avg_support, length, share_count);
            }
            share_count = 0;
            
            if (visit_idx > 0 && (trav.visit(visit_idx - 1).node_id() != 0 || trav.visit(visit_idx).node_id() != 0)) {
                // get the edge support
                edge_t edge = to_edge(graph, trav.visit(visit_idx - 1), visit);
                min_support = get_edge_support(edge);
                length = get_edge_length(edge, ref_offsets);
                if (edge_counts.count(edge)) {
                    share_count = edge_counts[edge];                    
                } 
                update_support(trav_idx, min_support, min_support, length, share_count);
            }
        }
    }

    // correct for case where no exclusive support found
    // or we're ignoring some traversals vg tgt_set interface
    for (int i = 0; i < min_supports_min.size(); ++i) {
        if (!has_support[i]) {
            min_supports_min[i] = Support();
            min_supports_avg[i] = Support();
        }
    }

    bool use_avg_trav_support = max_trav_size >= average_traversal_support_switch_threshold;
    bool use_avg_node_support = max_trav_size >= average_node_support_switch_threshold;

    if (use_avg_trav_support) {
        vector<Support>& tot_supports = use_avg_node_support ? tot_supports_avg : tot_supports_min;
        for (int i = 0; i < tot_supports.size(); ++i) {
            if (tot_sizes[i] > 0) {
                tot_supports[i] /= (double)tot_sizes[i];
            } else {
                tot_supports[i] = Support();
            }
        }
        return tot_supports;
    } else {
        return use_avg_node_support ? min_supports_avg : min_supports_min;
    }
}   

vector<int> TraversalSupportFinder::get_traversal_sizes(const vector<SnarlTraversal>& traversals) const {
    vector<int> sizes(traversals.size(), 0);
    for (int i = 0; i < traversals.size(); ++i) {
        for (int j = 0; j < traversals[i].visit_size(); ++j) {
            if (traversals[i].visit(j).node_id() != 0) {
                sizes[i] += graph.get_length(graph.get_handle(traversals[i].visit(j).node_id()));
            } else {
                // just summing up the snarl contents, which isn't a great heuristic but will
                // help in some cases
                pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(
                    snarl_manager.into_which_snarl(traversals[i].visit(j)), graph, true);
                for (id_t node_id : contents.first) {
                    sizes[i] += graph.get_length(graph.get_handle(node_id));
                }
            }
        }
    }
    return sizes;
    
}

size_t TraversalSupportFinder::get_average_traversal_support_switch_threshold() const {
    return average_traversal_support_switch_threshold;
}

unordered_map<id_t, size_t> TraversalSupportFinder::get_ref_offsets(const SnarlTraversal& ref_trav) const {
    unordered_map<id_t, size_t> ref_offsets;
    size_t offset = 0;
    for (int i = 0; i < ref_trav.visit_size(); ++i) {
        const Visit& visit = ref_trav.visit(i);
        if (visit.node_id() != 0) {
            if (visit.backward()) {
                offset += graph.get_length(graph.get_handle(visit.node_id()));
                ref_offsets[visit.node_id()] = offset;
            } else {
                ref_offsets[visit.node_id()] = offset;
                offset += graph.get_length(graph.get_handle(visit.node_id()));
            }
        }
    }
    return ref_offsets;
}

void TraversalSupportFinder::set_support_switch_threshold(size_t trav_thresh, size_t node_thresh) {
    average_traversal_support_switch_threshold = trav_thresh;
    average_node_support_switch_threshold = node_thresh;
}

PackedTraversalSupportFinder::PackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager) :
    TraversalSupportFinder(*dynamic_cast<const PathHandleGraph*>(packer.get_graph()), snarl_manager),
    packer(packer) {
}

PackedTraversalSupportFinder::~PackedTraversalSupportFinder() {
}

Support PackedTraversalSupportFinder::get_edge_support(const edge_t& edge) const {
    return get_edge_support(graph.get_id(edge.first), graph.get_is_reverse(edge.first),
                            graph.get_id(edge.second), graph.get_is_reverse(edge.second));
}

Support PackedTraversalSupportFinder::get_edge_support(id_t from, bool from_reverse,
                                                       id_t to, bool to_reverse) const {
    Edge proto_edge;
    proto_edge.set_from(from);
    proto_edge.set_from_start(from_reverse);
    proto_edge.set_to(to);
    proto_edge.set_to_end(to_reverse);
    Support support;
    support.set_forward(packer.edge_coverage(proto_edge));
    return support;
}

Support PackedTraversalSupportFinder::get_min_node_support(id_t node) const {
    Position pos;
    pos.set_node_id(node);
    size_t offset = packer.position_in_basis(pos);
    size_t coverage = packer.coverage_at_position(offset);
    size_t end_offset = offset + graph.get_length(graph.get_handle(node));
    for (int i = offset + 1; i < end_offset; ++i) {
        coverage = min(coverage, packer.coverage_at_position(i));
    }
    Support support;
    support.set_forward(coverage);
    return support;
}

Support PackedTraversalSupportFinder::get_avg_node_support(id_t node) const {
    Position pos;
    pos.set_node_id(node);
    size_t offset = packer.position_in_basis(pos);
    size_t coverage = 0;
    size_t length = graph.get_length(graph.get_handle(node));
    for (int i = 0; i < length; ++i) {
        coverage += packer.coverage_at_position(offset + i);
    }
    Support support;
    support.set_forward((double)coverage / (double)length);
    return support;
}


CachedPackedTraversalSupportFinder::CachedPackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager, size_t cache_size) :
    PackedTraversalSupportFinder(packer, snarl_manager) {
    size_t num_threads = get_thread_count();
    min_node_support_cache.resize(num_threads);
    avg_node_support_cache.resize(num_threads);
    edge_support_cache.resize(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        min_node_support_cache[i] = new LRUCache<nid_t, Support>(cache_size);
        avg_node_support_cache[i] = new LRUCache<nid_t, Support>(cache_size);
        edge_support_cache[i] = new LRUCache<edge_t, Support>(cache_size);
    }
}

CachedPackedTraversalSupportFinder::~CachedPackedTraversalSupportFinder() {
    for (size_t i = 0; i < min_node_support_cache.size(); ++i) {
        delete min_node_support_cache[i];
        delete avg_node_support_cache[i];
        delete edge_support_cache[i];
    }
}

Support CachedPackedTraversalSupportFinder::get_edge_support(id_t from, bool from_reverse,
                                                             id_t to, bool to_reverse) const {
    const HandleGraph* graph = packer.get_graph();
    edge_t edge = graph->edge_handle(graph->get_handle(from, from_reverse),
                                     graph->get_handle(to, to_reverse));
    
    auto& support_cache = *edge_support_cache[omp_get_thread_num()];
    pair<Support, bool> cached = support_cache.retrieve(edge);
    if (cached.second == true) {
        return cached.first;
    } else {
        Support support = PackedTraversalSupportFinder::get_edge_support(from, from_reverse, to, to_reverse);
        support_cache.put(edge, support);
        return support;
    }
}

Support CachedPackedTraversalSupportFinder::get_min_node_support(id_t node) const {
    auto& support_cache = *min_node_support_cache[omp_get_thread_num()];
    pair<Support, bool> cached = support_cache.retrieve(node);
    if (cached.second == true) {
        return cached.first;
    } else {
        Support support = PackedTraversalSupportFinder::get_min_node_support(node);
        support_cache.put(node, support);
        return support;
    }
}

Support CachedPackedTraversalSupportFinder::get_avg_node_support(id_t node) const {
    auto& support_cache = *avg_node_support_cache[omp_get_thread_num()];
    pair<Support, bool> cached = support_cache.retrieve(node);
    if (cached.second == true) {
        return cached.first;
    } else {
        Support support = PackedTraversalSupportFinder::get_avg_node_support(node);
        support_cache.put(node, support);
        return support;
    }
}


}
