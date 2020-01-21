#include "snarl_caller.hpp"
#include "genotypekit.hpp"

//#define debug

namespace vg {

SnarlCaller::~SnarlCaller() {
}

function<bool(const SnarlTraversal&, int)> SnarlCaller::get_skip_allele_fn() const {
    // default implementation says don't skip anything
    return [](const SnarlTraversal&, int) { assert(false); return false; };
}

SupportBasedSnarlCaller::SupportBasedSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                                 TraversalSupportFinder& support_finder) :
    graph(graph),
    snarl_manager(snarl_manager),
    support_finder(support_finder) {

}

SupportBasedSnarlCaller::~SupportBasedSnarlCaller() {
    
}

void SupportBasedSnarlCaller::update_vcf_info(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              const vector<int>& genotype,
                                              const unique_ptr<CallInfo>& call_info,
                                              const string& sample_name,
                                              vcflib::Variant& variant) {
    
    
}

const TraversalSupportFinder& SupportBasedSnarlCaller::get_support_finder() const {
    return support_finder;
}

int SupportBasedSnarlCaller::get_min_total_support_for_call() const {
    return min_total_support_for_call;
}

void SupportBasedSnarlCaller::set_min_supports(double min_mad_for_call, double min_support_for_call, double min_site_support) {
    if (min_mad_for_call >= 0) {
        min_mad_for_filter = min_mad_for_call;
    }
    if (min_support_for_call >= 0) {
        min_total_support_for_call = min_support_for_call;
    }
    if (min_site_support >= 0) {
        min_site_depth = min_site_support;
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

function<bool(const SnarlTraversal&, int)> SupportBasedSnarlCaller::get_skip_allele_fn() const {
    // port over cutoff used in old support caller (there avg support used all the time, here
    // we use the same toggles as when genotyping)
    return [&](const SnarlTraversal& trav, int iteration) -> bool {
        return support_val(support_finder.get_traversal_support(trav)) < pow(2, iteration) * min_alt_path_support;
    };
}

RatioSupportSnarlCaller::RatioSupportSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                                 TraversalSupportFinder& support_finder) :
    SupportBasedSnarlCaller(graph, snarl_manager, support_finder)  {
}

RatioSupportSnarlCaller::~RatioSupportSnarlCaller() {
    
}

void RatioSupportSnarlCaller::set_het_bias(double het_bias, double ref_het_bias) {
    // want to move away from ugly hacks that treat the reference traversal differently,
    // so keep all these set the same
    if (het_bias >= 0) {
        max_het_bias = het_bias;
        max_ref_het_bias = het_bias;
        max_indel_het_bias = het_bias;
    }
    if (ref_het_bias >= 0) {
        max_ref_het_bias = ref_het_bias;
    }
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> RatioSupportSnarlCaller::genotype(const Snarl& snarl,
                                                                                       const vector<SnarlTraversal>& traversals,
                                                                                       int ref_trav_idx,
                                                                                       int ploidy,
                                                                                       const string& ref_path_name,
                                                                                       pair<size_t, size_t> ref_range) { 
    
#ifdef debug
    cerr << "Support calling site " << pb2json(snarl) << endl;
#endif

    // get the traversal sizes
    vector<int> traversal_sizes = support_finder.get_traversal_sizes(traversals);

    // get the supports of each traversal independently
    vector<Support> supports = support_finder.get_traversal_set_support(traversals, {}, {}, false, {}, false, ref_trav_idx);
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
    vector<Support> secondary_exclusive_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, {}, true, {}, false, ref_trav_idx);    
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
    vector<Support> secondary_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, {}, false, {}, false, ref_trav_idx);
    int second_best_allele = get_best_support(secondary_supports, {skips});

    // get the supports of each traversal in light of second best
    // for special case where we may call two alts, with each having less support than ref
    vector<Support> tertiary_supports;
    int third_best_allele = -1;
    if (second_best_allele != -1) {
        // prune out traversals whose exclusive support relative to second best doesn't pass cut
        vector<Support> tertiary_exclusive_supports = support_finder.get_traversal_set_support(traversals, {second_best_allele}, {}, true, {}, false, ref_trav_idx);
        skips.push_back(best_allele);
        skips.push_back(second_best_allele);
        for (int i = 0; i < tertiary_exclusive_supports.size(); ++i) {
            double bias = get_bias(traversal_sizes, i, second_best_allele, ref_trav_idx);
            if (support_val(tertiary_exclusive_supports[i]) * bias <= support_val(supports[second_best_allele])) {
                skips.push_back(i);
            }
        }
        tertiary_supports = support_finder.get_traversal_set_support(traversals, {second_best_allele}, {}, false, {}, false, ref_trav_idx);
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
        return make_pair(vector<int>(1, best_allele), unique_ptr<SnarlCaller::CallInfo>());
    }
    // Call 1/2 : REF-Alt1/Alt2 even if Alt2 has only third best support
    else if (ploidy >= 2 &&
             third_best_allele > 0 &&
             best_allele == ref_trav_idx &&
             max_ma_bias * 
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
             max_ma_bias * 
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

    // Todo: specify call_info to use new interface, then fix up update_vcf_info to read it,
    // and move common logic up to SupportBasedCaller if possible.
    return make_pair(genotype, unique_ptr<SnarlCaller::CallInfo>());
}

void RatioSupportSnarlCaller::update_vcf_info(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              const vector<int>& genotype,
                                              const unique_ptr<CallInfo>& call_info,
                                              const string& sample_name,
                                              vcflib::Variant& variant) {

    assert(traversals.size() == variant.alleles.size());

    set<int> called_allele_set(genotype.begin(), genotype.end());
    vector<int> shared_travs;
    if (called_allele_set.size() > 1) {
        shared_travs.push_back(genotype[0]);
    }
    // compute the support of our called alleles
    vector<Support> allele_supports = support_finder.get_traversal_genotype_support(traversals, genotype, {}, 0);
    
    // Compute the total support for all the alts that will be appearing
    Support total_support = std::accumulate(allele_supports.begin(), allele_supports.end(), Support());    
        
    // Set up the depth format field
    variant.format.push_back("DP");
    // And allelic depth
    variant.format.push_back("AD");
    // And the log likelihood from the assignment of reads among the
    // present alleles
    variant.format.push_back("XADL");
    // Also the alt allele depth
    variant.format.push_back("XAAD");
    
    // And total alt allele depth for the alt alleles
    Support alt_support;
    // Find the min total support of anything called
    double min_site_support = called_allele_set.size() > 0 ? INFINITY : 0;

    if (!allele_supports.empty()) { //only add info if we made a call
        for (int allele = 0; allele < traversals.size(); ++allele) {
            bool is_called = called_allele_set.count(allele);
            auto& support = allele_supports[allele];
            
            // Set up allele-specific stats for the allele
            variant.samples[sample_name]["AD"].push_back(std::to_string((int64_t)round(total(support))));
                
            if (allele != 0) {
                // It's not the primary reference allele
                alt_support += support;
            }
            
            // Min all the total supports from the alleles called as present
            if (is_called) {
                min_site_support = min(min_site_support, total(support));
            }
        }
    }
    
            
    // Find the binomial bias between the called alleles, if multiple were called.
    double ad_log_likelihood = INFINITY;
    if (called_allele_set.size() == 2) {
        int best_allele = genotype[0];
        int second_best_allele = genotype[1];
        // How many of the less common one do we have?
        size_t successes = round(total(allele_supports[second_best_allele]));
        // Out of how many chances
        size_t trials = successes + (size_t) round(total(allele_supports[best_allele]));
                
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

    // And store the minimum support just to be clear
    variant.format.push_back("MAD");
    variant.samples[sample_name]["MAD"].push_back(std::to_string((int)(min_site_support)));

    // Now do the filters
    variant.filter = "PASS";            
    if (min_site_support < min_mad_for_filter) {
        // Apply Min Allele Depth cutoff across all alleles (even ref)
        variant.filter = "lowad";
    } else if (min_ad_log_likelihood_for_filter != 0 &&
               ad_log_likelihood < min_ad_log_likelihood_for_filter) {
        // We have a het, but the assignment of reads between the two branches is just too weird
        variant.filter = "lowxadl";
    } else if ((int64_t)round(total(total_support)) < min_site_depth) {
        // we don't have enough support to want to make a call
        variant.filter = "lowdepth";
    }
}

void RatioSupportSnarlCaller::update_vcf_header(string& header) const {
    header += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    header += "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
    header += "##FORMAT=<ID=MAD,Number=1,Type=Integer,Description=\"Minimum site allele depth\">\n";
    header += "##FORMAT=<ID=XADL,Number=1,Type=Float,Description=\"Likelihood of allelic depths for called alleles\">\n";
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    header += "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">\n";
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    header += "##FILTER=<ID=lowad,Description=\"Variant does not meet minimum allele read support threshold of " +
        std::to_string(min_mad_for_filter) + "\">\n";
    header += "##FILTER=<ID=lowxadl,Description=\"Variant has AD log likelihood less than " +
        std::to_string(min_ad_log_likelihood_for_filter) + "\">\n";
    header += "##FILTER=<ID=lowdepth,Description=\"Variant has read depth less than " +
        std::to_string(min_site_depth) + "\">\n";    
}

double RatioSupportSnarlCaller::get_bias(const vector<int>& traversal_sizes, int best_trav,
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


PoissonSupportSnarlCaller::PoissonSupportSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                                     TraversalSupportFinder& support_finder,
                                                     const algorithms::BinnedDepthIndex& depth_index) :
    SupportBasedSnarlCaller(graph, snarl_manager, support_finder),
    depth_index(depth_index) {
    
}
    
PoissonSupportSnarlCaller::~PoissonSupportSnarlCaller() {
    
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> PoissonSupportSnarlCaller::genotype(const Snarl& snarl,
                                                                                         const vector<SnarlTraversal>& traversals,
                                                                                         int ref_trav_idx,
                                                                                         int ploidy,
                                                                                         const string& ref_path_name,
                                                                                         pair<size_t, size_t> ref_range) {
    
    
#ifdef debug
    cerr << "Poisson Support calling site " << pb2json(snarl)
         << " on path " << ref_path_name << ":" << ref_range.first << "-" << ref_range.second << endl;
#endif

    assert(ploidy == 2 || ploidy == 1);
    
    // get the traversal sizes
    vector<int> traversal_sizes = support_finder.get_traversal_sizes(traversals);

    // get the supports of each traversal independently
    vector<Support> supports = support_finder.get_traversal_set_support(traversals, {}, {}, false, {}, false, ref_trav_idx);

    // sort the traversals by support
    vector<int> ranked_traversals = rank_by_support(supports);
    size_t max_trav = std::min(top_k, (size_t)ranked_traversals.size());
    size_t max_sec_trav = std::min(top_m, (size_t)ranked_traversals.size());
    // take the top-m traversals in order to check against the top traversal
    set<int> top_traversals(ranked_traversals.begin(), ranked_traversals.begin() + max_sec_trav);

    // the candidate genotypes and their supports.  the numbers here are alleles as indexed in traversals[]
    set<vector<int>> candidates;
    // we always consider the reference allele

    // pre-filter out some alleles based on poor exclusive support
    set<int> skips;

    // consider each of the top 25 traversals as our top_traversal
    for (int i = 0; i < max_trav; ++i) {
        
        int best_allele = ranked_traversals[i];

        if (skips.count(best_allele)) {
            continue;
        }
        if (support_val(supports[best_allele]) < min_total_support_for_call && candidates.size() >= max_trav) {
            break;
        }

        if (ploidy == 1) {
            candidates.insert({best_allele});
        } else {
            assert(ploidy == 2);
        
            // we prune out traversals whose exclusive support (structure that is not shared with best traversal)
            // doesn't meet a certain cutoff
            vector<Support> secondary_exclusive_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, top_traversals, true, {}, false, ref_trav_idx);
            for (int j = 0; j < secondary_exclusive_supports.size(); ++j) {
                if (j != best_allele &&
                    support_val(secondary_exclusive_supports[j]) < min_total_support_for_call &&
                    support_val(secondary_exclusive_supports[j]) < support_val(supports[j])) {
                    skips.insert(j);
                }
            }

            // get the supports of each traversal in light of best
            vector<Support> secondary_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, top_traversals, false, {}, false, ref_trav_idx);
            vector<int> ranked_secondary_traversals = rank_by_support(secondary_supports);

            // add the homozygous genotype for our best allele
            candidates.insert({best_allele, best_allele});

            // now look at the top-k second-best traversals
            size_t sec_count = 0;
            for (int j = 0; j < ranked_secondary_traversals.size() && sec_count < top_k; ++j) {
                int second_best_allele = ranked_secondary_traversals[j];
                if (support_val(secondary_supports[second_best_allele]) < min_total_support_for_call && candidates.size() >= max_trav) {
                    break;
                }
                if (!skips.count(second_best_allele) && second_best_allele != best_allele) {
                    // canonical ordering for our set
                    candidates.insert({min(best_allele, second_best_allele), max(best_allele, second_best_allele)});
                    // also make sure we have our homozygous genotype for the second best allele
                    candidates.insert({second_best_allele, second_best_allele});
                    ++sec_count;
                }
            }
        }
    }

    // expected depth from our coverage
    auto depth_info = algorithms::get_depth_from_index(depth_index, ref_path_name, ref_range.first, ref_range.second);
    double exp_depth = depth_info.first;
    assert(!isnan(exp_depth));
    // variance/std-err can be nan when binsize < 2.  We just clamp it to 0
    double depth_err = depth_info.second ? !isnan(depth_info.second) : 0.;

    // genotype (log) likelihoods
    double best_genotype_likelihood = -numeric_limits<double>::max();
    double second_best_genotype_likelihood = -numeric_limits<double>::max();
    double total_likelihood = 0;
    vector<int> best_genotype;
    for (const auto& candidate : candidates) {
        double gl = genotype_likelihood(candidate, traversals, top_traversals, ref_trav_idx, exp_depth, depth_err);
        if (gl > best_genotype_likelihood) {
            second_best_genotype_likelihood = best_genotype_likelihood;
            best_genotype_likelihood = gl;
            best_genotype = candidate;
        } else if (gl > second_best_genotype_likelihood) {
            assert(gl <= best_genotype_likelihood);
            second_best_genotype_likelihood = gl;
        }
        total_likelihood = total_likelihood == 0 ? gl : add_log(total_likelihood, gl);
    }

    PoissonCallInfo* call_info = new PoissonCallInfo();

    call_info->posterior = 0;
    if (!candidates.empty()) {
        // compute the posterior from our likelihoods using a uniform prior
        call_info->posterior = best_genotype_likelihood - log(candidates.size()) - total_likelihood;  
    }
    
    // GQ computed as here https://gatk.broadinstitute.org/hc/en-us/articles/360035890451?id=11075
    // as difference between best and second best likelihoods
    call_info->gq = 0;
    if (!isnan(best_genotype_likelihood) && !isnan(second_best_genotype_likelihood)) {
        call_info->gq = logprob_to_phred(second_best_genotype_likelihood) - logprob_to_phred(best_genotype_likelihood);
    }

    call_info->expected_depth = exp_depth;

#ifdef debug
    cerr << " best genotype: "; for (auto a : best_genotype) {cerr << a <<",";} cerr << " gl=" << best_genotype_likelihood << endl;
#endif
    return make_pair(best_genotype, unique_ptr<SnarlCaller::CallInfo>(call_info));
}

double PoissonSupportSnarlCaller::genotype_likelihood(const vector<int>& genotype,
                                                      const vector<SnarlTraversal>& traversals,
                                                      const set<int>& trav_subset, 
                                                      int ref_trav_idx, double exp_depth, double depth_err) {
    
    assert(genotype.size() == 1 || genotype.size() == 2);

    // get the genotype support
    vector<Support> genotype_supports = support_finder.get_traversal_genotype_support(traversals, genotype, trav_subset, ref_trav_idx);

    // get the total support over the site
    Support total_site_support = std::accumulate(genotype_supports.begin(), genotype_supports.end(), Support());    

    // get the total support of traversals *not* in the genotype
    Support total_other_support;
    set<int> genotype_set(genotype.begin(), genotype.end());
    for (int i = 0; i < traversals.size(); ++i) {
        if (!genotype_set.count(i)) {
            total_other_support += genotype_supports[i];
        }
    }
 
    // split the homozygous support into two
    // from now on we'll treat it like two separate observations, each with half coverage
    vector<Support> fixed_genotype_supports = genotype_supports;
    if (std::equal(genotype.begin() + 1, genotype.end(), genotype.begin())) {
        for (int i = 0; i < genotype_supports.size(); ++i) {
            fixed_genotype_supports[i] = genotype_supports[i] / (double)genotype.size();
        }
    }
    
    // how many reads would we expect to not map to our genotype due to error
    // Note: The bin size is set quite a bit smaller than originally intended as it seems to
    // help nearly nevery benchmark.  But the small bin sizes means that depth_err, the
    // error from the binned coverage, is way too high and including it only causes trouble.
    // tldr: just use the baseline_mapping_error constant and forget about depth_err for now. 
    //double error_rate = std::min(0.05, depth_err + baseline_mapping_error);
    double error_rate = baseline_mapping_error;
    double other_poisson_lambda = error_rate * exp_depth; //support_val(total_site_support);

    // and our likelihood for the unmapped reads we see:
    double other_log_likelihood = poisson_prob_ln(std::round(support_val(total_other_support)), other_poisson_lambda);

    // how many reads do we expect for an allele?  we use the expected coverage and just
    // divide it out by the size of the genotype.  
    double allele_poisson_lambda = (exp_depth / (double)genotype.size()) * (1. - error_rate);

#ifdef debug
    cerr << "Computing prob of genotype: {";
    for (int i = 0; i < genotype.size(); ++i) {
        cerr << genotype[i] << ",";
    } 
    cerr << "}: tot_other_sup = " << total_other_support << " tot site sup = " << total_site_support 
         << " exp-depth = " << exp_depth << " depth-err = " << depth_err << " other-lambda = " << other_poisson_lambda
         << " allele-lambda " << allele_poisson_lambda << endl;
#endif
    
    // now we compute the likelihood of our genotype
    double alleles_log_likelihood = 0;
    for (int allele : genotype) {
        const Support& allele_support = fixed_genotype_supports[allele];
        double allele_ll = poisson_prob_ln(std::round(support_val(allele_support)), allele_poisson_lambda);
        alleles_log_likelihood += allele_ll;

#ifdef debug
        cerr << "  a[" << allele <<"]=" << " sup=" << genotype_supports[allele] << " fix-sup=" << allele_support
             << " prob " << allele_ll << endl;
#endif        
    }

#ifdef debug
    cerr  << " allele-log-prob " << alleles_log_likelihood << " other-log-prob " << other_log_likelihood
          << " total-prob " << (alleles_log_likelihood + other_log_likelihood) << endl;
#endif

    return alleles_log_likelihood + other_log_likelihood;
}

void PoissonSupportSnarlCaller::update_vcf_info(const Snarl& snarl,
                                                const vector<SnarlTraversal>& traversals,
                                                const vector<int>& genotype,
                                                const unique_ptr<CallInfo>& call_info,
                                                const string& sample_name,
                                                vcflib::Variant& variant) {

    assert(traversals.size() == variant.alleles.size());

    // get the traversal sizes
    vector<int> traversal_sizes = support_finder.get_traversal_sizes(traversals);

    // get the genotype support
    vector<Support> genotype_supports = support_finder.get_traversal_genotype_support(traversals, genotype, {}, 0);

    // get the genotype_

    // Get the depth of the site
    Support total_site_support = std::accumulate(genotype_supports.begin(), genotype_supports.end(), Support());    
    double total_site_depth = support_val(total_site_support);

    // Set the variant's total depth            
    string depth_string = std::to_string((int64_t)round(total_site_depth));
    variant.format.push_back("DP");
    variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth
            
    // And for the sample
    variant.samples[sample_name]["DP"].push_back(depth_string);            

    // get the allele depths
    variant.format.push_back("AD");

    set<int> genotype_set(genotype.begin(), genotype.end());
    double min_site_support = genotype.size() > 0 ? INFINITY : 0;

    // update the allele depths
    for (int i = 0; i < traversals.size(); ++i) {
        Support allele_support = genotype_supports[i];
        variant.samples[sample_name]["AD"].push_back(std::to_string((int64_t)round(support_val(allele_support))));
        if (genotype_set.count(i)) {
            // update the minimum support
            min_site_support = min(min_site_support, total(genotype_supports[i]));
        }
    }

    // get the genotype likelihoods
    vector<double> gen_likelihoods;
    double gen_likelihood;    
    variant.format.push_back("GL");

    // expected depth from our coverage.  we look at the reference-range from the snarl plus a bit of padding,
    // averaging over every depth bin this touches.  todo: adaptively compute nearby coverage without bins
    // (requires VCFGenotyper to be refactored to require a PathPositionIndex)
    size_t longest_traversal = *max_element(traversal_sizes.begin(), traversal_sizes.end());
    size_t padding = (depth_padding_factor * longest_traversal) / 2;
    pair<size_t, size_t> ref_range = make_pair(max((long)0, (long)(variant.position - padding)),
					       variant.position + variant.ref.length() + padding);
    auto depth_info = algorithms::get_depth_from_index(depth_index, variant.sequenceName, ref_range.first, ref_range.second);
    double exp_depth = depth_info.first;
    assert(!isnan(exp_depth));
    // variance/std-err can be nan when binsize < 2.  We just clamp it to 0
    double depth_err = depth_info.second ? !isnan(depth_info.second) : 0.;

    double total_likelihood = 0.;
    double ref_likelihood = 1.;
    double alt_likelihood = 0.;
    
    // assume ploidy 2
    for (int i = 0; i < traversals.size(); ++i) {
        for (int j = i; j < traversals.size(); ++j) {
            double gl = genotype_likelihood({i, j}, traversals, {}, 0, exp_depth, depth_err);
            gen_likelihoods.push_back(gl);
            if (vector<int>({i, j}) == genotype || vector<int>({j,i}) == genotype) {
                gen_likelihood = gl;
            }
            if (i == 0 && j == 0) {
                ref_likelihood = gl;
            } else {
                alt_likelihood = alt_likelihood == 0. ? gl : add_log(alt_likelihood, gl);
            }
            total_likelihood = total_likelihood == 0 ? gl : add_log(total_likelihood, gl);
            // convert from natural log to log10 by dividing by ln(10)
            variant.samples[sample_name]["GL"].push_back(std::to_string(gl / 2.30258));
        }
    }

    const SnarlCaller::CallInfo* s_call_info = call_info.get();
    const PoissonCallInfo* p_call_info = dynamic_cast<const PoissonCallInfo*>(call_info.get());
    variant.format.push_back("GQ");
    variant.samples[sample_name]["GQ"].push_back(std::to_string(min((int)256, max((int)0, (int)p_call_info->gq))));

    variant.format.push_back("GP");
    variant.samples[sample_name]["GP"].push_back(std::to_string(p_call_info->posterior));

    variant.format.push_back("XD");
    variant.samples[sample_name]["XD"].push_back(std::to_string(p_call_info->expected_depth));

    // The QUAL field is the probability that we have variation as a PHRED score (of wrongness)
    // We derive this from the posterior probability of the reference genotype.
    // But if it's a reference call, we take the total of all the alts
    variant.quality = 0;
    if (genotype.size() == 2) {
        // our flat prior and p[traversal coverage]
        double posterior = -log(gen_likelihoods.size()) - total_likelihood;
        if (genotype[0] != 0 || genotype[1] != 0) {
            posterior += ref_likelihood;
        } else {
            posterior += alt_likelihood;
        }
        variant.quality = logprob_to_phred(posterior);
    }

    // Minmum allele depth.  This historically has been our QUAL field for the sole reason
    // that it was better than anything else at making ROC curves
    variant.format.push_back("MAD");
    variant.samples[sample_name]["MAD"].push_back(std::to_string((int)(min_site_support)));

    // Now do the filters
    // todo: fix and share with other caller
    variant.filter = "PASS";            
    if (min_site_support < min_mad_for_filter) {
        // Apply Min Allele Depth cutoff across all alleles (even ref)
        variant.filter = "lowad";
    } else if ((int64_t)round(total_site_depth) < min_site_depth) {
        // we don't have enough support to want to make a call
        variant.filter = "lowdepth";
    }
}

void PoissonSupportSnarlCaller::update_vcf_header(string& header) const {
    header += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    header += "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
    header += "##FORMAT=<ID=MAD,Number=1,Type=Integer,Description=\"Minimum site allele depth\">\n";
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    header += "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n";
    header += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the Phred-scaled probability estimate of the called genotype\">\n";
    header += "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype Probability, the log-scaled posterior probability of the called genotype\">\n";
    header += "##FORMAT=<ID=XD,Number=1,Type=Float,Description=\"eXpected Depth, background coverage as used for the Poisson model\">\n";
    header += "##FILTER=<ID=lowad,Description=\"Variant does not meet minimum allele read support threshold of " +
        std::to_string(min_mad_for_filter) + "\">\n";
    header += "##FILTER=<ID=lowdepth,Description=\"Variant has read depth less than " +
        std::to_string(min_site_depth) + "\">\n";    
}

vector<int> PoissonSupportSnarlCaller::rank_by_support(const vector<Support>& supports) {
    vector<int> ranks(supports.size());
    for (int i = 0; i < supports.size(); ++i) {
        ranks[i] = i;
    }
    std::sort(ranks.begin(), ranks.end(), [&](int a, int b) {
            return support_val(supports[a]) > support_val(supports[b]);
        });
    return ranks;
}


}
