#include "traversal_genotyper.hpp"
#include "genotypekit.hpp"

namespace vg {

SupportBasedTraversalGenotyper::SupportBasedTraversalGenotyper(const PathHandleGraph& graph,
                                                               bool use_avg_node_support,
                                                               bool use_avg_trav_support) :
    graph(graph), use_avg_node_support(use_avg_node_support), use_avg_trav_support(use_avg_trav_support) {
    if (use_avg_node_support) {
        get_node_support = [&] (id_t node) { return get_avg_node_support(node); };
    } else {
        get_node_support = [&] (id_t node) { return get_min_node_support(node); };
    }
}

SupportBasedTraversalGenotyper::~SupportBasedTraversalGenotyper() {
    
}


vector<int> SupportBasedTraversalGenotyper::genotype(const Snarl& snarl,
                                                     const vector<SnarlTraversal>& traversals,
                                                     int ref_trav_idx,
                                                     int ploidy) {
    return {0, 0};
}

Support SupportBasedTraversalGenotyper::get_traversal_support(const SnarlTraversal& traversal) {
    return get_traversal_set_support({traversal}, {0}, false).at(0);
}

vector<Support> SupportBasedTraversalGenotyper::get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                                          const vector<int>& trav_indxes,
                                                                          bool exclusive_only) {

    // pass 1: how many times have we seen a node or edge
    unordered_map<id_t, int> node_counts;
    unordered_map<edge_t, int> edge_counts;
    map<Snarl, int> child_counts;

    for (const SnarlTraversal& trav : traversals) {
        for (int i = 0; i < trav.visit_size(); ++i) {
            const Visit& visit = trav.visit(i);
            if (visit.node_id() != 0) {
                // Count the node once
                if (node_counts.count(visit.node_id())) {
                    node_counts[visit.node_id()] += 1;
                } else {
                    node_counts[visit.node_id()] = 1;
                }
            } else {
                // Count the child once
                if (child_counts.count(visit.snarl())) {
                    child_counts[visit.snarl()] += 1;
                } else {
                    child_counts[visit.snarl()] = 1;
                }
            }
            if (i > 0) {
                edge_t edge = to_edge(graph, trav.visit(i - 1), visit);
                // Count the edge once
                if (edge_counts.count(edge)) {
                    edge_counts[edge] += 1;
                } else {
                    edge_counts[edge] = 1;
                }
            }
        }
    }

    // pass 2: get the supports
    vector<Support> min_supports;
    vector<Support> tot_supports; // weighted by lengths
    vector<int> tot_sizes; // to compute average from to_supports;
    
    for (const SnarlTraversal& trav : traversals) {
        for (int i = 0; i < trav.visit_size(); ++i) {
            const Visit& visit = trav.visit(i);
            Support support;
            int64_t length;
            double scale_factor = 1.;

            if (visit.node_id() != 0) {
                // get the node support
                support = get_node_support(visit.node_id());
                length = graph.get_length(graph.get_handle(visit.node_id()));
                if (node_counts.count(visit.node_id())) {
                    scale_factor = 1. / (double)node_counts[visit.node_id()];
                }
            } else {
                // get the child support
                tie(support, length) = get_child_support(visit.snarl());
                if (child_counts.count(visit.snarl())) {
                    scale_factor = 1. / (double)child_counts[visit.snarl()];
                }
            }            
            if (i > 0) {
                // get the edge support
                edge_t edge = to_edge(graph, trav.visit(i - 1), visit);
                support = get_edge_support(edge);
                length = get_edge_length(edge);
                if (edge_counts.count(edge)) {
                    scale_factor = 1. / (double)edge_counts[edge];
                }
            }

            // apply the scaling
            if (exclusive_only && scale_factor < 1.) {
                scale_factor = 0.;
            }
            support = support * scale_factor;

            // update our support values for the traversal
            if (i == 0) {
                min_supports.push_back(support);
                tot_supports.push_back(support);
                tot_sizes.push_back(length);
            } else {
                min_supports.back() = support_min(min_supports.back(), support);
                tot_supports.back() += support;
                tot_sizes.back() += length;
            }
        }
    }

    if (use_avg_trav_support) {
        for (int i = 0; i < tot_supports.size(); ++i) {
            if (tot_sizes[i] > 0) {
                tot_supports[i] /= (double)tot_sizes[i];
            } else {
                tot_supports[i] = Support();
            }
        }
        return tot_supports;
    } else {
        return min_supports;
    }
}   

}
