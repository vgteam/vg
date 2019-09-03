#include "snarl_caller.hpp"
#include "genotypekit.hpp"

//#define debug

namespace vg {

SnarlCaller::~SnarlCaller() {
}

SupportBasedSnarlCaller::SupportBasedSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager) :
    graph(graph),
    snarl_manager(snarl_manager) {
}

SupportBasedSnarlCaller::~SupportBasedSnarlCaller() {
    
}


vector<int> SupportBasedSnarlCaller::genotype(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              int ref_trav_idx,
                                              int ploidy) {

#ifdef debug
    cerr << "Support calling site " << pb2json(snarl) << endl;
#endif

    // get the traversal sizes
    vector<int> traversal_sizes = get_traversal_sizes(traversals);

    // get the supports of each traversal independently
    vector<Support> supports = get_traversal_set_support(traversals, {}, false, ref_trav_idx);
    int best_allele = get_best_support(supports, {});

#ifdef debug
    for (int i = 0; i < traversals.size(); ++i) {
        cerr << "trav " << i << " size = " << traversal_sizes[i] << " support = " << support_val(supports[i]);
        if (i == ref_trav_idx) {
            cerr << " [Reference traversal]";
        }
        cerr << endl;
    }
#endif

    // we prune out traversals whose exclusive support (structure that is not shared with best traversal)
    // doesn't meet a certain cutoff
    vector<Support> secondary_exclusive_supports = get_traversal_set_support(traversals, {best_allele}, true, ref_trav_idx);    
    vector<int> skips = {best_allele};
    for (int i = 0; i < secondary_exclusive_supports.size(); ++i) {
        double bias = get_bias(traversal_sizes, i, best_allele, ref_trav_idx);
#ifdef debug
        cerr << "trav " << i << " exclusive support " << support_val(secondary_exclusive_supports[i])
             << " * bias " << bias << " vs " << support_val(supports[best_allele]) << endl;
#endif
        if (i != best_allele && support_val(secondary_exclusive_supports[i]) * bias <= support_val(supports[best_allele])) {
            skips.push_back(i);
        }
    }
    // get the supports of each traversal in light of best
    vector<Support> secondary_supports = get_traversal_set_support(traversals, {best_allele}, false, ref_trav_idx);
    int second_best_allele = get_best_support(secondary_supports, {skips});

    // get the supports of each traversal in light of second best
    // for special case where we may call two alts, with each having less support than ref
    vector<Support> tertiary_supports;
    int third_best_allele = -1;
    if (second_best_allele != -1) {
        // prune out traversals whose exclusive support relative to second best doesn't pass cut
        vector<Support> tertiary_exclusive_supports = get_traversal_set_support(traversals, {second_best_allele}, true, ref_trav_idx);
        skips.push_back(best_allele);
        skips.push_back(second_best_allele);
        for (int i = 0; i < tertiary_exclusive_supports.size(); ++i) {
            double bias = get_bias(traversal_sizes, i, second_best_allele, ref_trav_idx);
            if (support_val(tertiary_exclusive_supports[i]) * bias <= support_val(supports[second_best_allele])) {
                skips.push_back(i);
            }
        }
        tertiary_supports = get_traversal_set_support(traversals, {second_best_allele}, false, ref_trav_idx);
        third_best_allele = get_best_support(tertiary_supports, skips);
    }

    
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
                            
#ifdef debug
    cerr << "best allele=" << best_allele << ", best sup=" << best_support << " and "
         << "2nd_best_allele=" << second_best_allele << ", 2nd best sup=" << second_best_support << " and "
         << "3rd_best_allele=" << third_best_allele << ", 3rd best sup=" << third_best_support << endl;
        
    if (support_val(second_best_support) > 0) {
        cerr << "Bias: (limit " << get_bias(traversal_sizes, best_allele, second_best_allele, ref_trav_idx)  << "):"
             << support_val(best_support)/support_val(second_best_support) << endl;
    }
        
    cerr << get_bias(traversal_sizes, best_allele, second_best_allele, ref_trav_idx) * support_val(second_best_support) << " vs "
         << support_val(best_support) << endl;
            
    cerr << total(second_best_support) << " vs " << min_total_support_for_call << endl;
#endif

    // Single ploidy case when doing recursive genotyping.  Just return the best allele
    if (ploidy == 1) {
        return {best_allele};
    }
    // Call 1/2 : REF-Alt1/Alt2 even if Alt2 has only third best support
    else if (ploidy >= 2 &&
             third_best_allele > 0 &&
             best_allele == ref_trav_idx && 
             get_bias(traversal_sizes, second_best_allele, third_best_allele, ref_trav_idx) *
             support_val(third_best_support) >= support_val(best_support) &&
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
    // Call 1/2 : REF-Alt1/Alt2 even if Alt2 has only third best support (but ref is second best)
    else if (ploidy >= 2 &&
             third_best_allele > 0 &&
             second_best_allele == ref_trav_idx && 
             get_bias(traversal_sizes, best_allele, third_best_allele, ref_trav_idx) *
             support_val(third_best_support) >= support_val(second_best_support) &&
             total(best_support) > min_total_support_for_call &&
             total(third_best_support) > min_total_support_for_call) {
        // There's a second best allele and third best allele, and it's not too biased to call,
        // and both alleles exceed the minimum to call them present, and the
        // second-best and third-best alleles have enough support that it won't torpedo the
        // variant.
            
#ifdef debug
        cerr << "Call as second best/third best" << endl;
#endif
        // Say both are present
        genotype = {best_allele, third_best_allele};
    }    
    else if (ploidy >= 2 &&
             second_best_allele != -1 &&
             get_bias(traversal_sizes, best_allele, second_best_allele, ref_trav_idx) *
             support_val(second_best_support) >= support_val(best_support) &&
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

void SupportBasedSnarlCaller::update_vcf_info(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              const vector<int>& genotype,                                         
                                              const string& sample_name,
                                              vcflib::Variant& variant) {

    /// Compute the supports of the genotype
    /// (we're doing this a second time to fit the (unnecessarily modular?) interface where vcf
    ///  info is updated independently of calling)
    vector<SnarlTraversal> genotype_travs;
    map<int, int> allele_map; // map allele from "genotype" to position in "genotype_travs"/"allele_supports"
    int genotype_travs_ref_idx = -1;
    for (auto allele : genotype) {
        assert(allele < traversals.size());
        if (!allele_map.count(allele)) {
            genotype_travs.push_back(traversals[allele]);
            allele_map[allele] = genotype_travs.size() - 1;
        }
    }
    vector<int> shared_travs;
    if (genotype.size() >= 2 && genotype[0] != genotype[1]) {
        shared_travs.push_back(0);
    }
    int ref_trav_idx = allele_map.count(0) ? allele_map[0] : -1;
    vector<Support> allele_supports = get_traversal_set_support(genotype_travs, shared_travs, false, ref_trav_idx);

    // Set up the depth format field
    variant.format.push_back("DP");
    // And allelic depth
    variant.format.push_back("AD");
    // And the log likelihood from the assignment of reads among the
    // present alleles
    variant.format.push_back("XADL");
    // Also the alt allele depth
    variant.format.push_back("XAAD");

    // Compute the total support for all the alts that will be appearing
    Support total_support;
    // And total alt allele depth for the alt alleles
    Support alt_support;
    // Find the min total support of anything called
    double min_site_support = allele_supports.size() > 0 ? INFINITY : 0;

    if (!allele_supports.empty()) { //only add info if we made a call
        for (auto allele_pair : allele_map) {
            int genotype_allele = allele_pair.first;
            int allele_supports_idx = allele_pair.second;
            // For all the alleles we are using, look at the support.
            auto& support = allele_supports[allele_supports_idx];
                
            // Set up allele-specific stats for the allele
            variant.samples[sample_name]["AD"].push_back(std::to_string((int64_t)round(total(support))));
                
            // Sum up into total depth
            total_support += support;
                
            if (genotype_allele != 0) {
                // It's not the primary reference allele
                alt_support += support;
            }
            
            // Min all the total supports from the alleles called as present    
            min_site_support = min(min_site_support, total(support));
        }
    }
    
            
    // Find the binomial bias between the called alleles, if multiple were called.
    double ad_log_likelihood = INFINITY;
    if (allele_map.size() == 2) {
        int best_allele = genotype[0];
        int second_best_allele = genotype[1];
        // How many of the less common one do we have?
        size_t successes = round(total(allele_supports[allele_map[second_best_allele]]));
        // Out of how many chances
        size_t trials = successes + (size_t) round(total(allele_supports[allele_map[best_allele]]));
                
        assert(trials >= successes);
                
        // How weird is that?                
        ad_log_likelihood = binomial_cmf_ln(prob_to_logprob((real_t) 0.5), trials, successes);
                
        assert(!std::isnan(ad_log_likelihood));
                
        variant.samples[sample_name]["XADL"].push_back(std::to_string(ad_log_likelihood));
    } else {
        // No need to assign reads between two alleles
        variant.samples[sample_name]["XADL"].push_back(".");
    }

    // Set the variant's total depth            
    string depth_string = std::to_string((int64_t)round(total(total_support)));
    variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth
            
    // And for the sample
    variant.samples[sample_name]["DP"].push_back(depth_string);            
            
    // And its depth of non-0 alleles
    variant.samples[sample_name]["XAAD"].push_back(std::to_string((int64_t)round(total(alt_support))));

    // Set the total support of the min allele as the variant quality
    variant.quality = min_site_support;

    // Now do the filters
    variant.filter = "PASS";            
    if (min_site_support < min_mad_for_filter) {
        // Apply Min Allele Depth cutoff across all alleles (even ref)
        variant.filter = "lowad";
    } else if (min_ad_log_likelihood_for_filter != 0 &&
               ad_log_likelihood < min_ad_log_likelihood_for_filter) {
        // We have a het, but the assignment of reads between the two branches is just too weird
        variant.filter = "lowxadl";
    }
}

void SupportBasedSnarlCaller::update_vcf_header(string& header) const {
    header += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    header += "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
    header += "##FORMAT=<ID=XADL,Number=1,Type=Float,Description=\"Likelihood of allelic depths for called alleles\">\n";
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    header += "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">\n";
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    header += "##FILTER=<ID=lowad,Description=\"Variant does not meet minimum allele read support threshold of " +
        std::to_string(min_mad_for_filter) + "\">\n";
    header += "##FILTER=<ID=lowxadl,Description=\"Variant has AD log likelihood less than " +
        std::to_string(min_ad_log_likelihood_for_filter) + "\">\n";
}

int64_t SupportBasedSnarlCaller::get_edge_length(const edge_t& edge, const unordered_map<id_t, size_t>& ref_offsets) const {
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

tuple<Support, Support, int> SupportBasedSnarlCaller::get_child_support(const Snarl& snarl) const {
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
    return std::tie(child_max_support, child_avg_support, child_size);
}


Support SupportBasedSnarlCaller::get_traversal_support(const SnarlTraversal& traversal) const {
    return get_traversal_set_support({traversal}, {}, false).at(0);
}

vector<Support> SupportBasedSnarlCaller::get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                                   const vector<int>& shared_travs,
                                                                   bool exclusive_only,
                                                                   int ref_trav_idx) const {

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
            // note: there is no edge between adjacent snarls as they overlap
            // on their endpoints. 
            if (i > 0 && (trav.visit(i - 1).node_id() != 0 || trav.visit(i).node_id() != 0)) {
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
                               const Support& avg_support, int length, int share_count) {
        // keep track of overall size of longest traversal
        tot_sizes_all[trav_idx] += length;
        max_trav_size = std::max(tot_sizes_all[trav_idx], max_trav_size);
        
        // apply the scaling
        double scale_factor = (exclusive_only && share_count > 0) ? 0. : 1. / (1. + share_count);
        
        // when looking at exclusive support, we don't normalize by skipped lengths
        if (scale_factor != 0 || !exclusive_only) {
            has_support[trav_idx] = true;
            Support scaled_support_min = min_support * scale_factor * length;
            Support scaled_support_avg = avg_support * scale_factor * length;

            tot_supports_min[trav_idx] += scaled_support_min;
            tot_supports_avg[trav_idx] += scaled_support_avg;
            tot_sizes[trav_idx] += length;
            min_supports_min[trav_idx] = support_min(min_supports_min[trav_idx], scaled_support_min);
            min_supports_avg[trav_idx] = support_min(min_supports_avg[trav_idx], scaled_support_avg);
        }
    };

    for (int trav_idx = 0; trav_idx < traversals.size(); ++trav_idx) {
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

double SupportBasedSnarlCaller::get_bias(const vector<int>& traversal_sizes, int best_trav,
                                         int second_best_trav, int ref_trav_idx) const {
    bool is_indel = ((best_trav >= 0 && traversal_sizes[best_trav] != traversal_sizes[ref_trav_idx]) ||
                     (second_best_trav >=0 && traversal_sizes[second_best_trav] != traversal_sizes[ref_trav_idx]));

    double bias_limit = 1;

    if (best_trav >= 0 && second_best_trav >=0) {
        if (best_trav == ref_trav_idx) {
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
    }
    return bias_limit;
}

unordered_map<id_t, size_t> SupportBasedSnarlCaller::get_ref_offsets(const SnarlTraversal& ref_trav) const {
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

PackedSupportSnarlCaller::PackedSupportSnarlCaller(const Packer& packer, SnarlManager& snarl_manager) :
    SupportBasedSnarlCaller(*dynamic_cast<const PathHandleGraph*>(packer.graph), snarl_manager),
    packer(packer) {
}

PackedSupportSnarlCaller::~PackedSupportSnarlCaller() {
}

Support PackedSupportSnarlCaller::get_edge_support(const edge_t& edge) const {
    return get_edge_support(graph.get_id(edge.first), graph.get_is_reverse(edge.first),
                            graph.get_id(edge.second), graph.get_is_reverse(edge.second));
}

Support PackedSupportSnarlCaller::get_edge_support(id_t from, bool from_reverse,
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
