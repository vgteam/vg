#ifndef VG_SNARL_CALLER_HPP_INCLUDED
#define VG_SNARL_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"
#include "genotypekit.hpp"
#include "traversal_support.hpp"
#include "algorithms/coverage_depth.hpp"

namespace vg {

using namespace std;


/**
 * SnarlCaller: Given a list of traversals through a site, come up with a genotype
 * come up with a genotype
 */ 
class SnarlCaller {
public:
    virtual ~SnarlCaller();

    /// implementation-dependent metadata for calls that get paseed between genotype()
    /// and update_vcf_info().
    struct CallInfo {
        virtual ~CallInfo() = default;
    };

    /// Get the genotype of a site
    /// snarl : site
    /// traversals : all traversals to consider
    /// ref_trav_idx : index of reference path traversal in traversals (in case it needs special treatment)
    /// ref_path : the reference path associated with the snarl
    /// ref_range : the interval along the reference path (forward coordinates) spanned by snarl
    virtual pair<vector<int>, unique_ptr<CallInfo>> genotype(const Snarl& snarl,
                                                             const vector<SnarlTraversal>& traversals,
                                                             int ref_trav_idx,
                                                             int ploidy,
                                                             const string& ref_path_name,
                                                             pair<size_t, size_t> ref_range) = 0;
    
    /// Update INFO and FORMAT fields of the called variant
    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const unique_ptr<CallInfo>& call_info,
                                 const string& sample_name,
                                 vcflib::Variant& variant) = 0;

    /// Define any header fields needed by the above
    virtual void update_vcf_header(string& header) const = 0;

    /// Optional method used for pruning searches
    virtual function<bool(const SnarlTraversal&, int iteration)> get_skip_allele_fn() const;
};

/**
 * Interface for a caller that relies on a TraversalSupportFinder
 * and has a few very basic support-based cutoffs
 * Not every exciting but is currently required for the LegacySupportCaller
 * which needs this to interface with the RepresentativeTraversalFinder
 */ 
class SupportBasedSnarlCaller : public SnarlCaller {
public:
    SupportBasedSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                            TraversalSupportFinder& support_finder);

    virtual ~SupportBasedSnarlCaller();

    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const unique_ptr<CallInfo>& call_info,
                                 const string& sample_name,
                                 vcflib::Variant& variant);

    /// Set some of the parameters
    void set_min_supports(double min_mad_for_call, double min_support_for_call, double min_site_support);
    
    /// Get the traversal support finder
    TraversalSupportFinder& get_support_finder() const;

    /// Get the minimum total support for call
    virtual int get_min_total_support_for_call() const;

    /// Use min_alt_path_support threshold as cutoff
    virtual function<bool(const SnarlTraversal&, int iteration)> get_skip_allele_fn() const;

protected:

    /// Get the best support out of a list of supports, ignoring skips
    static int get_best_support(const vector<Support>& supports, const vector<int>& skips);

    /// Relic from old code
    static double support_val(const Support& support) { return total(support); };

    const PathHandleGraph& graph;

    SnarlManager& snarl_manager;    

    /// Get support from traversals
    TraversalSupportFinder& support_finder;
    
    /// What's the minimum integer number of reads that must support a call? We
    /// don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    int min_total_support_for_call = 1;
    /// what's the minimum ref or alt allele depth to give a PASS in the filter
    /// column? Also used as a min actual support for a second-best allele call
    size_t min_mad_for_filter = 1;
    /// what's the minimum total support (over all alleles) of the site to make
    /// a call
    size_t min_site_depth = 3;
    /// used only for pruning alleles in the VCFTraversalFinder:  minimum support
    /// of an allele's alt-path for it to be considered in the brute-force enumeration
    double min_alt_path_support = 0.5;
};


/**
 * Find the genotype of some traversals in a site using read support and
 * a bias ratio to tell heterozygous from homozygous
 */ 
class RatioSupportSnarlCaller : public SupportBasedSnarlCaller {
public:
    RatioSupportSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                            TraversalSupportFinder& support_finder);
    virtual ~RatioSupportSnarlCaller();

    /// Set some of the parameters
    void set_het_bias(double het_bias, double ref_het_bias = 0.);

    /// Get the genotype of a site
    virtual pair<vector<int>, unique_ptr<CallInfo>> genotype(const Snarl& snarl,
                                                             const vector<SnarlTraversal>& traversals,
                                                             int ref_trav_idx,
                                                             int ploidy,
                                                             const string& ref_path_name,
                                                             pair<size_t, size_t> ref_range);

    /// Update INFO and FORMAT fields of the called variant
    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const unique_ptr<CallInfo>& call_info,
                                 const string& sample_name,
                                 vcflib::Variant& variant);

    /// Define any header fields needed by the above
    virtual void update_vcf_header(string& header) const;

protected:

    /// Get the bias used to for comparing two traversals
    /// (It differrs heuristically depending whether they are alt/ref/het/hom/snp/indel
    ///  see tuning parameters below)
    double get_bias(const vector<int>& traversal_sizes, int best_trav,
                    int second_best_trav, int ref_trav_idx) const;

    /// get a map of the beginning of a node (in forward orientation) on a traversal
    /// used for up-weighting large deletion edges in complex snarls with average support
    unordered_map<id_t, size_t> get_ref_offsets(const SnarlTraversal& ref_trav) const;

    /// Tuning

    /// What fraction of the reads supporting an alt are we willing to discount?
    /// At 2, if twice the reads support one allele as the other, we'll call
    /// homozygous instead of heterozygous. At infinity, every call will be
    /// heterozygous if even one read supports each allele.
    double max_het_bias = 6;
    /// Like above, but applied to ref / alt ratio (instead of alt / ref)
    double max_ref_het_bias = 6;
    /// Like the max het bias, but applies to novel indels.
    double max_indel_het_bias = 6;
    /// Used for calling 1/2 calls.  If both alts (times this bias) are greater than
    /// the reference, the call is made.  set to 0 to deactivate.
    double max_ma_bias = 0;
    /// what's the min log likelihood for allele depth assignments to PASS?
    double min_ad_log_likelihood_for_filter = -9;    
};

/**
 * Find the genotype of some traversals in a site using read support 
 * and a Poisson model based on expected depth.  Inspired, in part,
 * by Paragraph, which uses a similar approach for genotyping break points
 * 
 **/ 
class PoissonSupportSnarlCaller : public SupportBasedSnarlCaller {
public:
    PoissonSupportSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                              TraversalSupportFinder& support_finder,
                              const algorithms::BinnedDepthIndex& depth_index,
                              bool use_mapq);
    virtual ~PoissonSupportSnarlCaller();

    struct PoissonCallInfo : public SnarlCaller::CallInfo {
        virtual ~PoissonCallInfo() = default;
        double gq;
        double posterior;
        double expected_depth;
        double depth_err;
        int max_trav_size;
    };

    /// Set some parameters
    void set_baseline_error(double small_variant_error, double large_variant_error);
    /// These are multipliers applied to the errors if the site has an insertion
    void set_insertion_bias(double insertion_threshold, double small_insertion_bias, double large_insertion_bias);    

    /// Get the genotype of a site
    virtual pair<vector<int>, unique_ptr<CallInfo>>  genotype(const Snarl& snarl,
                                                              const vector<SnarlTraversal>& traversals,
                                                              int ref_trav_idx,
                                                              int ploidy,
                                                              const string& ref_path_name,
                                                              pair<size_t, size_t> ref_range);
    
    /// Update INFO and FORMAT fields of the called variant
    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const unique_ptr<CallInfo>& call_info,
                                 const string& sample_name,
                                 vcflib::Variant& variant);

    /// Define any header fields needed by the above
    virtual void update_vcf_header(string& header) const;

protected:

    /// Compute likelihood of genotype as product of poisson probabilities
    /// P[allele1] * P[allle2] * P[uncalled alleles]
    /// Homozygous alleles are split into two, with half support each
    /// The (natural) logoarithm is returned
    /// If trav_subset is not empty, traversals outside that set (and genotype)
    /// will be ignored to save time
    double genotype_likelihood(const vector<int>& genotype,
                               const vector<SnarlTraversal>& traversals,
                               const set<int>& trav_subset,
                               const vector<int>& traversal_sizes,
                               const vector<double>& traversal_mapqs,
                               int ref_trav_idx, double exp_depth, double depth_err,
                               int max_trav_size, int ref_trav_size);

    /// Rank supports
    vector<int> rank_by_support(const vector<Support>& supports);

    /// Error rates are different for small and large variants, which depend
    /// more on base and mapping qualities respectively.  The switch threshold
    /// is in TraversalSupportFinder.  Error stats from the Packer object
    /// get added to these baselines when computing the scores. 
    
    /// Baseline error rate for larger variants
    double  baseline_error_large = 0.001;
    /// Baseline error rate for smaller variants
    double  baseline_error_small = 0.005;
    /// multiply error by this much in pressence of insertion
    double insertion_bias_large = 10.;
    double insertion_bias_small = 1.;
    /// a site is an insertion if one (supported)allele is this many times bigger than another
    /// unlike above, default comes from call_main.cpp (todo: straighten this out?)
    double insertion_threshold = 5.;

    /// Consider up to the top-k traversals (based on support) for genotyping
    size_t top_k = 20;
    /// Consider up to the tom-m secondary traversals (based on support) for each top traversal
    /// (so at most top_k * top_m considered)
    size_t top_m = 100;

    /// padding to apply wrt to longest traversal to snarl ranges when looking up binned depth
    double depth_padding_factor = 1.;
    
    /// Map path name to <mean, std_err> of depth coverage from the packer
    const algorithms::BinnedDepthIndex& depth_index;

    /// MAPQ information is available from the packer and we want to use it
    bool use_mapq;

};


// debug helpers
inline string to_string(const HandleGraph& graph, handle_t handle) {
    return std::to_string(graph.get_id(handle)) + ":" + std::to_string(graph.get_is_reverse(handle));
}
inline string to_string(const HandleGraph& graph, edge_t edge) {
    return to_string(graph, edge.first) + " -> " + to_string(graph, edge.second);
}

}
#endif
