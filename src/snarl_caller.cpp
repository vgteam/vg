#include "snarl_caller.hpp"
#include "genotypekit.hpp"

#define debug

namespace vg {

SnarlCaller::~SnarlCaller() {
}

SupportBasedSnarlCaller::SupportBasedSnarlCaller(const PathHandleGraph& graph,
                                                 bool use_avg_node_support,
                                                 bool use_avg_trav_support) :
    graph(graph), use_avg_node_support(use_avg_node_support), use_avg_trav_support(use_avg_trav_support) {
    if (use_avg_node_support) {
        get_node_support = [&] (id_t node) { return get_avg_node_support(node); };
    } else {
        get_node_support = [&] (id_t node) { return get_min_node_support(node); };
    }
}

SupportBasedSnarlCaller::~SupportBasedSnarlCaller() {
    
}


vector<int> SupportBasedSnarlCaller::genotype(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              int ref_trav_idx,
                                              int ploidy) {

    assert(ploidy == 2);

    // get the traversal sizes
    vector<int> traversal_sizes = get_traversal_sizes(traversals);

    // get the supports of each traversal independently
    vector<Support> supports = get_traversal_set_support(traversals, {}, false);
    int best_allele = get_best_support(supports, {});

    // we prune out traversals whose exclusive support (structure that is not shared with best traversal)
    // doesn't meet a certain cutoff
    vector<Support> secondary_exclusive_supports = get_traversal_set_support(traversals, {best_allele}, true);
    vector<int> skips = {best_allele};
    for (int i = 0; i < secondary_exclusive_supports.size(); ++i) {
        double bias = get_bias(traversal_sizes, i, best_allele, ref_trav_idx);
        if (support_val(secondary_exclusive_supports[i]) * bias <= support_val(supports[best_allele])) {
            skips.push_back(i);
        }
    }
    // get the supports of each traversal in light of best
    vector<Support> secondary_supports = get_traversal_set_support(traversals, {best_allele}, false);
    int second_best_allele = get_best_support(secondary_supports, {skips});

    // get the supports of each traversal in light of second best
    // for special case where we may call two alts, with each having less support than ref
    vector<Support> tertiary_supports;
    int third_best_allele = -1;
    if (second_best_allele != -1) {
        // prune out traversals whose exclusive support relative to second best doesn't pass cut
        vector<Support> tertiary_exclusive_supports = get_traversal_set_support(traversals, {second_best_allele}, true);
        skips = {best_allele, second_best_allele};
        for (int i = 0; i < tertiary_exclusive_supports.size(); ++i) {
            double bias = get_bias(traversal_sizes, i, second_best_allele, ref_trav_idx);
            if (support_val(tertiary_exclusive_supports[i]) * bias <= support_val(supports[second_best_allele])) {
                skips.push_back(i);
            }
        }
        tertiary_supports = get_traversal_set_support(traversals, {second_best_allele}, false);
        third_best_allele = get_best_support(tertiary_supports, skips);
    }


    // Decide if we're an indel by looking at the traversal sizes
    bool is_indel = traversal_sizes[best_allele] != traversal_sizes[ref_trav_idx] ||
        (second_best_allele != -1 && traversal_sizes[ref_trav_idx] != traversal_sizes[second_best_allele]);
    bool is_indel_ma_2 = (second_best_allele != -1 && traversal_sizes[ref_trav_idx] != traversal_sizes[second_best_allele]);
    bool is_indel_ma_3 = (third_best_allele != -1 && traversal_sizes[ref_trav_idx] != traversal_sizes[third_best_allele]);
    
    // Now make a genotype call at this site, up to the allowed copy number
    vector<int> genotype;
        
    // How much support do we have for the top two alleles?
    Support site_support = supports.at(best_allele);
    if(second_best_allele != -1) {
        site_support += supports.at(second_best_allele);
    }
    
    // Pull out the different supports. Some of them may be the same.
    Support best_support = supports.at(best_allele);
    Support second_best_support; // Defaults to 0
    if(second_best_allele != -1) {
        second_best_support = supports.at(second_best_allele);
    }
    Support third_best_support;
    if (third_best_allele != -1) {
        third_best_support = supports.at(third_best_allele);
    }
                    
    double bias_limit;
    if (best_allele == 0) {
        // Use ref bias limit
            
        // We decide closeness differently depending on whether best is ref
        // or not. In practice, we use this to slightly penalize homozygous
        // ref calls (by setting max_ref_het_bias higher than max_het_bias)
        // and rather make a less supported alt call instead.  This boost
        // max sensitivity, and because everything is homozygous ref by
        // default in VCF, any downstream filters will effectively reset
        // these calls back to homozygous ref. TODO: This shouldn't apply
        // when off the primary path!
        bias_limit = max_ref_het_bias;
    } else if (is_indel) {
        // This is an indel
        // Use indel bias limit
        bias_limit = max_indel_het_bias;
    } else {
        // Use normal het bias limit
        bias_limit = max_het_bias;
    }
        
#ifdef debug
    cerr << best_allele << ", " << best_support << " and "
         << second_best_allele << ", " << second_best_support << endl;
        
    if (support_val(second_best_support) > 0) {
        cerr << "Bias: (limit " << bias_limit  << "):"
             << support_val(best_support)/support_val(second_best_support) << endl;
    }
        
    cerr << bias_limit * support_val(second_best_support) << " vs "
         << support_val(best_support) << endl;
            
    cerr << total(second_best_support) << " vs " << min_total_support_for_call << endl;
#endif

    // Call 1/2 : REF-Alt1/Alt2 even if Alt2 has only third best support
    if (ploidy >= 2 &&
        best_allele == ref_trav_idx && 
        third_best_allele > 0 &&
        is_indel_ma_3 &&
        max_indel_ma_bias * support_val(third_best_support) >= support_val(best_support) &&
        total(second_best_support) > min_total_support_for_call &&
        total(third_best_support) > min_total_support_for_call) {
        // There's a second best allele and third best allele, and it's not too biased to call,
        // and both alleles exceed the minimum to call them present, and the
        // second-best and third-best alleles have enough support that it won't torpedo the
        // variant.
            
#ifdef debug
        cerr << "Call as second best/third best" << endl;
#endif
        // Say both are present
        genotype = {second_best_allele, third_best_allele};
    }
    else if (ploidy >= 2 &&
             second_best_allele != -1 &&
             bias_limit * support_val(second_best_support) >= support_val(best_support) &&
             total(best_support) > min_total_support_for_call &&
             total(second_best_support) > min_total_support_for_call) {
        // There's a second best allele, and it's not too biased to call,
        // and both alleles exceed the minimum to call them present, and the
        // second-best allele has enough support that it won't torpedo the
        // variant.
            
#ifdef debug
        cerr << "Call as best/second best" << endl;
#endif
            
        // Say both are present
        genotype = {best_allele, second_best_allele};
            
    } else if (ploidy >= 2 && total(best_support) > min_total_support_for_call) {
        // The second best allele isn't present or isn't good enough,
        // but the best allele has enough coverage that we can just call
        // two of it.
            
#ifdef debug
        cerr << "Call as best/best" << endl;
#endif
            
        // Say the best is present twice
        genotype = {best_allele, best_allele};

    } else if (ploidy >= 1 && total(best_support) > min_total_support_for_call) {
        // We're only supposed to have one copy, and the best allele is good enough to call
            
#ifdef debug
        cerr << "Call as best" << endl;
#endif
        genotype = {best_allele};
    } else {
        // Either coverage is too low, or we aren't allowed any copies.
        // We can't really call this as anything.
            
#ifdef debug
        cerr << "Do not call" << endl;
#endif
            
    }

    return genotype;
}

int64_t SupportBasedSnarlCaller::get_edge_length(const edge_t& edge) const {
    return 1;
}

pair<Support, int> SupportBasedSnarlCaller::get_child_support(const Snarl& snarl) const {
    return make_pair(Support(), 1);
}


Support SupportBasedSnarlCaller::get_traversal_support(const SnarlTraversal& traversal) const {
    return get_traversal_set_support({traversal}, {}, false).at(0);
}

vector<Support> SupportBasedSnarlCaller::get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                                   const vector<int>& shared_travs,
                                                                   bool exclusive_only) const {

    // pass 1: how many times have we seen a node or edge
    unordered_map<id_t, int> node_counts;
    unordered_map<edge_t, int> edge_counts;
    map<Snarl, int> child_counts;

    for (auto trav_idx : shared_travs) {
        const SnarlTraversal& trav = traversals[trav_idx];
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

    int trav_offset = 1; // change to 0 to consider snarl endpoints
    
    for (const SnarlTraversal& trav : traversals) {
        for (int i = trav_offset; i < trav.visit_size() - trav_offset; ++i) {
            const Visit& visit = trav.visit(i);
            Support support;
            int64_t length;
            int share_count = 0;

            if (visit.node_id() != 0) {
                // get the node support
                support = get_node_support(visit.node_id());
                length = graph.get_length(graph.get_handle(visit.node_id()));
                if (node_counts.count(visit.node_id())) {
                    share_count = node_counts[visit.node_id()];
                }
            } else {
                // get the child support
                tie(support, length) = get_child_support(visit.snarl());
                if (child_counts.count(visit.snarl())) {
                    share_count = child_counts[visit.snarl()];
                }
            }            
            if (i > 0) {
                // get the edge support
                edge_t edge = to_edge(graph, trav.visit(i - 1), visit);
                support = get_edge_support(edge);
                length = get_edge_length(edge);
                if (edge_counts.count(edge)) {
                    share_count = edge_counts[edge];
                }
            }

            // apply the scaling
            double scale_factor = (exclusive_only && share_count > 0) ? 0. : 1. / (1. + share_count);
            support = support * scale_factor;

            // update our support values for the traversal
            if (i == trav_offset) {
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

int SupportBasedSnarlCaller::get_best_support(const vector<Support>& supports, const vector<int>& skips) {
    int best_allele = -1;
    for(size_t i = 0; i < supports.size(); i++) {
        if(std::find(skips.begin(), skips.end(), i) == skips.end() && (
               best_allele == -1 || support_val(supports[best_allele]) <= support_val(supports[i]))) {
            best_allele = i;
        }
    }
    return best_allele;
}

vector<int> SupportBasedSnarlCaller::get_traversal_sizes(const vector<SnarlTraversal>& traversals) const {
    vector<int> sizes(traversals.size(), 0);
    for (int i = 0; i < traversals.size(); ++i) {
        for (int j = 0; j < traversals[i].visit_size(); ++j) {
            // Note: we are not counting nested snarls
            if (traversals[i].visit(j).node_id() != 0) {
                sizes[i] += graph.get_length(graph.get_handle(traversals[i].visit(j).node_id()));
            }
        }
    }
    return sizes;
    
}

double SupportBasedSnarlCaller::get_bias(const vector<int>& traversal_sizes, int trav_idx1, int trav_idx2, int ref_trav_idx) const {
    assert(trav_idx1 >= 0 && trav_idx2 >=0);

    return 2;
    
}


PackedSupportSnarlCaller::PackedSupportSnarlCaller(const Packer& packer, bool use_avg_node_support, bool use_avg_trav_support) :
    SupportBasedSnarlCaller(*packer.xgidx, use_avg_node_support, use_avg_trav_support),
    packer(packer) {
}

PackedSupportSnarlCaller::~PackedSupportSnarlCaller() {
}

Support PackedSupportSnarlCaller::get_edge_support(const edge_t& edge) const {
    Edge proto_edge;
    proto_edge.set_from(graph.get_id(edge.first));
    proto_edge.set_from_start(graph.get_is_reverse(edge.first));
    proto_edge.set_to(graph.get_id(edge.second));
    proto_edge.set_to_end(graph.get_is_reverse(edge.second));
    Support support;
    support.set_forward(packer.edge_coverage(proto_edge));
    return support;
}

Support PackedSupportSnarlCaller::get_min_node_support(id_t node) const {
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

Support PackedSupportSnarlCaller::get_avg_node_support(id_t node) const {
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


}
