#include "snarl_caller.hpp"
#include "genotypekit.hpp"

//#define debug

namespace vg {

SnarlCaller::~SnarlCaller() {
}

function<bool(const SnarlTraversal&)> SnarlCaller::get_skip_allele_fn() const {
    // default implementation says don't skip anything
    return [](const SnarlTraversal&) { return false; };
}

SupportBasedSnarlCaller::SupportBasedSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                                 TraversalSupportFinder& support_finder) :
    graph(graph),
    snarl_manager(snarl_manager),
    support_finder(support_finder) {
}

SupportBasedSnarlCaller::~SupportBasedSnarlCaller() {
    
}

void SupportBasedSnarlCaller::set_het_bias(double het_bias, double ref_het_bias) {
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

vector<int> SupportBasedSnarlCaller::genotype(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              int ref_trav_idx,
                                              int ploidy) {

#ifdef debug
    cerr << "Support calling site " << pb2json(snarl) << endl;
#endif

    // get the traversal sizes
    vector<int> traversal_sizes = support_finder.get_traversal_sizes(traversals);

    // get the supports of each traversal independently
    vector<Support> supports = support_finder.get_traversal_set_support(traversals, {}, false, false, ref_trav_idx);
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
    vector<Support> secondary_exclusive_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, true, false, ref_trav_idx);    
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
    vector<Support> secondary_supports = support_finder.get_traversal_set_support(traversals, {best_allele}, false, false, ref_trav_idx);
    int second_best_allele = get_best_support(secondary_supports, {skips});

    // get the supports of each traversal in light of second best
    // for special case where we may call two alts, with each having less support than ref
    vector<Support> tertiary_supports;
    int third_best_allele = -1;
    if (second_best_allele != -1) {
        // prune out traversals whose exclusive support relative to second best doesn't pass cut
        vector<Support> tertiary_exclusive_supports = support_finder.get_traversal_set_support(traversals, {second_best_allele}, true, false, ref_trav_idx);
        skips.push_back(best_allele);
        skips.push_back(second_best_allele);
        for (int i = 0; i < tertiary_exclusive_supports.size(); ++i) {
            double bias = get_bias(traversal_sizes, i, second_best_allele, ref_trav_idx);
            if (support_val(tertiary_exclusive_supports[i]) * bias <= support_val(supports[second_best_allele])) {
                skips.push_back(i);
            }
        }
        tertiary_supports = support_finder.get_traversal_set_support(traversals, {second_best_allele}, false, false, ref_trav_idx);
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

    return genotype;
}

void SupportBasedSnarlCaller::update_vcf_info(const Snarl& snarl,
                                              const vector<SnarlTraversal>& traversals,
                                              const vector<int>& genotype,                                         
                                              const string& sample_name,
                                              vcflib::Variant& variant) {

    assert(traversals.size() == variant.alleles.size());

    set<int> called_allele_set(genotype.begin(), genotype.end());
    vector<int> shared_travs;
    if (called_allele_set.size() > 1) {
        shared_travs.push_back(genotype[0]);
    }
    // compute the support of our called alleles
    vector<Support> allele_supports = support_finder.get_traversal_set_support(traversals, shared_travs, false, false, 0);

    // get the support of our uncalled alleles, making sure to not include any called support
    // TODO: handle shared support within this set
    vector<Support> uncalled_supports = support_finder.get_traversal_set_support(traversals, genotype, false, true, 0);
        
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
    double min_site_support = called_allele_set.size() > 0 ? INFINITY : 0;

    if (!allele_supports.empty()) { //only add info if we made a call
        for (int allele = 0; allele < traversals.size(); ++allele) {
            bool is_called = called_allele_set.count(allele);
            auto& support = is_called ? allele_supports[allele] : uncalled_supports[allele];
            
            // Set up allele-specific stats for the allele
            variant.samples[sample_name]["AD"].push_back(std::to_string((int64_t)round(total(support))));
                
            // Sum up into total depth
            total_support += support;
                
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
    header += "##FILTER=<ID=lowdepth,Description=\"Variant has read depth less than " +
        std::to_string(min_site_depth) + "\">\n";    
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

function<bool(const SnarlTraversal&)> SupportBasedSnarlCaller::get_skip_allele_fn() const {
    // port over cutoff used in old support caller (there avg support used all the time, here
    // we use the same toggles as when genotyping)
    return [&](const SnarlTraversal& trav) -> bool {
        return support_val(support_finder.get_traversal_support(trav)) < min_alt_path_support;
    };
}


int SupportBasedSnarlCaller::get_min_total_support_for_call() const {
    return min_total_support_for_call;
}

const TraversalSupportFinder& SupportBasedSnarlCaller::get_support_finder() const {
    return support_finder;
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




}
