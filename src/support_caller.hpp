#ifndef VG_SUPPORT_CALLER_HPP_INCLUDED
#define VG_SUPPORT_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "pileup.hpp"
#include "path_index.hpp"
#include "genotypekit.hpp"
#include "option.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

/**
 * SupportCaller: take an augmented graph from a Caller and produce actual calls in a
 * VCF.
 */
class SupportCaller : public Configurable {

public:

    /**
     * Set up to call with default parameters.
     */
    SupportCaller() = default;
    
    /**
     * We use this to represent a contig in the primary path, with its index and coverage info.
     */
    class PrimaryPath {
    public:
        /**
         * Index the given path in the given augmented graph, and compute all
         * the coverage bin information with the given bin size.
         */
        PrimaryPath(SupportAugmentedGraph& augmented, const string& ref_path_name, size_t ref_bin_size); 
    
        /**
         * Get the support at the bin appropriate for the given primary path
         * offset.
         */
        const Support& get_support_at(size_t primary_path_offset) const;
        
        /**
         * Get the index of the bin that the given path position falls in.
         */
        size_t get_bin_index(size_t primary_path_offset) const;
    
        /**
         * Get the bin with minimal coverage.
         */
        size_t get_min_bin() const;
    
        /**
         * Get the bin with maximal coverage.
         */
        size_t get_max_bin() const;
    
        /**
         * Get the support in the given bin.
         */
        const Support& get_bin(size_t bin) const;
        
        /**
         * Get the total number of bins that the path is divided into.
         */
        size_t get_total_bins() const;
        
        /**
         * Get the average support over the path.
         */
        Support get_average_support() const;
        
        /**
         * Get the average support over a collection of paths.
         */
        static Support get_average_support(const map<string, PrimaryPath>& paths);
        
        /**
         * Get the total support for the path.
         */
        Support get_total_support() const;
    
        /**
         * Get the PathIndex for this primary path.
         */
        PathIndex& get_index();
        
        /**
         * Get the PathIndex for this primary path.
         */
        const PathIndex& get_index() const;
        
        /**
         * Gets the path name we are representing.
         */
        const string& get_name() const;
        
    protected:
        /// How wide is each coverage bin along the path?
        size_t ref_bin_size;
        
        /// This holds the index for this path
        PathIndex index;
        
        /// This holds the name of the path
        string name;
        
        /// What's the expected in each bin along the path? Coverage gets split
        /// evenly over both strands.
        vector<Support> binned_support;
        
        /// Which bin has min support?
        size_t min_bin;
        /// Which bin has max support?
        size_t max_bin;
        
        /// What's the total Support over every bin?
        Support total_support;
    };

    // We'll fill this in with a PrimaryPath for every primary reference path
    // that is specified or detected.
    map<string, PrimaryPath> primary_paths;
    
    /**
     * Produce calls for the given annotated augmented graph. If a
     * pileup_filename is provided, the pileup is loaded again and used to add
     * comments describing variants
     */
    void call(SupportAugmentedGraph& augmented, SnarlManager& site_manager, string pileup_filename = "");

    /** 
     * Get the support and size for each traversal in a list. Discount support
     * of minus_traversal if it's specified.  Use average_support_switch_threshold and
     * use_average_support to decide whether to return min or avg supports.
     */
    tuple<vector<Support>, vector<size_t> > get_traversal_supports_and_sizes(
        SupportAugmentedGraph& augmented, SnarlManager& snarl_manager, const Snarl& site,
        const vector<SnarlTraversal>& traversals,
        const vector<const SnarlTraversal*>& minus_traversals = {});

    /**
     * Get the min support, total support, bp size (to divide total by for average
     * support), optionally special-casing the material used by another traversal. 
     * Material used by another traversal only makes half its coverage available to this traversal.
     * If ref_traversal specified (and not same as traversal), its contents will be forbidden from
     * boosting the average support to avoid the actual variation from getting
     * diluted by shared reference.
     */
    tuple<Support, Support, size_t> get_traversal_support(
        SupportAugmentedGraph& augmented, SnarlManager& snarl_manager, const Snarl* site,
        const SnarlTraversal& traversal, const vector<const SnarlTraversal*>& already_used = {},
        const SnarlTraversal* ref_traversal = nullptr);

    /** 
     * Get the edge supports of an inversion.  This is to be used for computing genotypes with average support,
     * as the normal supports will almost always return 0/1 due to the nodes having the same support.  If 
     * the alleles don't refer to a simple inversion, an empty vector is returned (and should be ignored)
     */
    vector<Support>  get_inversion_supports(
        SupportAugmentedGraph& augmented, SnarlManager& snarl_manager, const Snarl& site,
        const vector<SnarlTraversal>& traversals, const vector<size_t>& traversal_sizes,
        int best_allele, int second_best_allele);   

    /**
     * For the given snarl, find the reference traversal, the best traversal,
     * and the second-best traversal, recursively, if any exist. These
     * traversals will be fully filled in with nodes.
     *
     * Only snarls which are ultrabubbles can be called.
     *
     * Expects the given baseline support for a diploid call.
     *
     * Will not return more than 1 + copy_budget SnarlTraversals, and will
     * return less if some copies are called as having the same traversal.
     *
     * Does not deduplicate agains the ref traversal; it may be the same as the
     * best or second-best.
     *
     * Uses the given copy number allowance, and emits a Locus for this Snarl
     * and any child Snarls.
     *
     * If no path through the Snarl can be found, emits no Locus and returns no
     * SnarlTraversals.
     *
     */
    vector<SnarlTraversal> find_best_traversals(SupportAugmentedGraph& augmented,
        SnarlManager& snarl_manager, TraversalFinder* finder, const Snarl& site,
        const Support& baseline_support, size_t copy_budget,
        function<void(const Locus&, const Snarl*, const vcflib::Variant*)> emit_locus);

    /**
     * Instead of just emitting the locus, use the site info from the VCF
     * to output a locus for each variant in the VCF that touchces our site. 
     */
    void recall_locus(Locus& locus, const Snarl& site, vector<SnarlTraversal>& traversals,
                      vector<vector<int>>& trav_alleles, 
                      vector<vcflib::Variant*>& site_variants,
                      function<void(const Locus&, const Snarl*, const vcflib::Variant*)> emit_locus);

    /** This function emits the given variant on the given primary path, as
     * VCF. It needs to take the site as an argument because it may be
     * called for children of the site we're working on right now.
     */
    void emit_variant(map<string, string>& contig_names_by_path_name,
                      vcflib::VariantCallFile& vcf,
                      SupportAugmentedGraph& augmented,
                      Support& baseline_support,
                      Support& global_baseline_support, 
                      const Locus& locus, PrimaryPath& primary_path, const Snarl* site);

    /** Like emit_variant, but use the given vcf variant as a template and just 
     * compute the genotype and info */
    void emit_recall_variant(map<string, string>& contig_names_by_path_name,
                             vcflib::VariantCallFile& vcf,
                             SupportAugmentedGraph& augmented,
                             Support& baseline_support,
                             Support& global_baseline_support, 
                             const Locus& locus, PrimaryPath& primary_path, const Snarl* site,
                             const vcflib::Variant* recall_variant);

    /** add the info fields to a variant and actually emit it 
     * (used by both emit_variant and emit_recall_variant) */
    void add_variant_info_and_emit(vcflib::Variant& variant, SupportAugmentedGraph& augmented,
                                   const Locus& locus, const Genotype& genotype,
                                   int best_allele, int second_best_allele,
                                   const vector<int>& used_alleles,
                                   Support& baseline_support, Support& global_baseline_support);

    /**
     * Decide if the given SnarlTraversal is included in the original base graph
     * (true), or if it represents a novel variant (false).
     *
     * Looks at the nodes in the traversal, and sees if their calls are
     * CALL_REFERENCE or not.
     *
     * Handles single-edge traversals.
     *
     */
    bool is_reference(const SnarlTraversal& trav, AugmentedGraph& augmented);
    
    /**
     * Decide if the given Path is included in the original base graph (true) or
     * if it represents a novel variant (false).
     *
     * Looks at the nodes, and sees if their calls are CALL_REFERENCE or not.
     *
     * The path can't be empty; it has to be anchored to something (probably the
     * start and end of the snarl it came from).
     */
    bool is_reference(const Path& path, AugmentedGraph& augmented);
    
    /**
     * Find the primary path, if any, that the given site is threaded onto.
     *
     * TODO: can only work by brute-force search.
     */
    map<string, PrimaryPath>::iterator find_path(const Snarl& site);

    /**
     * Find the lenght of a deletion deletion edge along (the first found) primary path.  If no path found
     * or not a deletion edge, return 0.
     */
    size_t get_deletion_length(const NodeSide& end1, const NodeSide& end2, SupportAugmentedGraph& augmented);

    /** 
     * Get the amount of support.  Can use this function to toggle between unweighted (total from genotypekit)
     * and quality-weighted (support_quality below) in one place.
     */
    function<double(const Support&)> support_val;

    static double support_quality(const Support& support) {
        return support.quality();
    }
    
    // Option variables
    
    /// Should we output in VCF (true) or Protobuf Locus (false) format?
    Option<bool> convert_to_vcf{this, "no-vcf", "V", true,
        "output variants in binary Loci format instead of text VCF format"};
    /// How big should our output buffer be?
    size_t locus_buffer_size = 1000;
    
    /// What are the names of the reference paths, if any, in the graph?
    Option<vector<string>> ref_path_names{this, "ref", "r", {},
        "use the path with the given name as a reference path (can repeat)"};
    /// What name should we give each contig in the VCF file? Autodetected from
    /// path names if empty or too short.
    Option<vector<string>> contig_name_overrides{this, "contig", "c", {},
        "use the given name as the VCF name for the corresponding reference path (can repeat)"};
    /// What should the total sequence length reported in the VCF header be for
    /// each contig? Autodetected from path lengths if empty or too short.
    Option<vector<size_t>> length_overrides{this, "length", "l", {},
        "override total sequence length in VCF for the corresponding reference path (can repeat)"};
    /// What name should we use for the sample in the VCF file?
    Option<string> sample_name{this, "sample", "S", "SAMPLE",
        "name the sample in the VCF with the given name"};
    /// How far should we offset positions of variants?
    Option<int64_t> variant_offset{this, "offset", "o", 0,
        "offset variant positions by this amount in VCF"};
    /// How many nodes should we be willing to look at on our path back to the
    /// primary path? Keep in mind we need to look at all valid paths (and all
    /// combinations thereof) until we find a valid pair.
    Option<int64_t> max_search_depth{this, "max-search-depth", "D", 1000,
        "maximum depth for path search"};
    /// How many search states should we allow on the DFS stack when searching
    /// for traversals?
    Option<int64_t> max_search_width{this, "max-search-width", "wWmMsS", 1000,
        "maximum width for path search"};
    
    
    /// What fraction of average coverage should be the minimum to call a
    /// variant (or a single copy)? Default to 0 because vg call is still
    /// applying depth thresholding
    Option<double> min_fraction_for_call{this, "min-cov-frac", "F", 0,
        "min fraction of average coverage at which to call"};
    /// What fraction of the reads supporting an alt are we willing to discount?
    /// At 2, if twice the reads support one allele as the other, we'll call
    /// homozygous instead of heterozygous. At infinity, every call will be
    /// heterozygous if even one read supports each allele.
    Option<double> max_het_bias{this, "max-het-bias", "H", 10,
        "max imbalance factor to call heterozygous, alt major on SNPs"};
    /// Like above, but applied to ref / alt ratio (instead of alt / ref)
    Option<double> max_ref_het_bias{this, "max-ref-bias", "R", 4.5,
        "max imbalance factor to call heterozygous, ref major"};
    /// Like the max het bias, but applies to novel indels.
    Option<double> max_indel_het_bias{this, "max-indel-het-bias", "I", 3,
        "max imbalance factor to call heterozygous, alt major on indels"};
    /// Like the max het bias, but applies to multiallelic indels.
    Option<double> max_indel_ma_bias{this, "max-indel-ma-bias", "G", 6,
        "max imbalance factor between ref and alt2 to call 1/2 double alt on indels"};
    /// What's the minimum integer number of reads that must support a call? We
    /// don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    Option<size_t> min_total_support_for_call{this, "min-count", "n", 1, 
        "min total supporting read count to call a variant"};
    /// Bin size used for counting coverage along the reference path.  The
    /// bin coverage is used for computing the probability of an allele
    /// of a certain depth
    Option<size_t> ref_bin_size{this, "bin-size", "B", 250,
        "bin size used for counting coverage"};
    /// On some graphs, we can't get the coverage because it's split over
    /// parallel paths.  Allow overriding here
    Option<double> expected_coverage{this, "avg-coverage", "C", 0.0,
        "specify expected coverage (instead of computing on reference)"};
    /// Should we use average support instead of minimum support for our
    /// calculations?
    Option<bool> use_average_support{this, "use-avg-support", "u", false,
        "use average instead of minimum support"};
    /// Max traversal length threshold at which we switch from minimum support
    /// to average support (so we don't use average support on pairs of adjacent
    /// errors and miscall them, but we do use it on long runs of reference
    /// inside a deletion where the min support might not be representative.
    Option<size_t> average_support_switch_threshold{this, "use-avg-support-above", "uUaAtT", 100,
        "use average instead of minimum support for sites this long or longer"};
    
    /// What's the maximum number of bubble path combinations we can explore
    /// while finding one with maximum support?
    size_t max_bubble_paths = 100;
    /// what's the minimum ref or alt allele depth to give a PASS in the filter
    /// column? Also used as a min actual support for a second-best allele call
    Option<size_t> min_mad_for_filter{this, "min-mad", "E", 5,
        "min. ref/alt allele depth to PASS filter or be a second-best allele"};
    /// what's the maximum total depth to give a PASS in the filter column
    Option<size_t> max_dp_for_filter{this, "max-dp", "MmDdAaXxPp", 0,
        "max depth to PASS filter (0 for unlimited)"};
    /// what's the maximum total depth to give a PASS in the filter column, as a
    /// multiple of the global baseline coverage?
    Option<double> max_dp_multiple_for_filter{this, "max-dp-multiple", "MmDdAaXxPp", 0,
        "max portion of global expected depth to PASS filter (0 for unlimited)"};
    /// what's the maximum total depth to give a PASS in the filter column, as a
    /// multiple of the local baseline coverage?
    Option<double> max_local_dp_multiple_for_filter{this, "max-local-dp-multiple", "MmLlOoDdAaXxPp", 0,
        "max portion of local expected depth to PASS filter (0 for unlimited)"};
    /// what's the min log likelihood for allele depth assignments to PASS?
    Option<double> min_ad_log_likelihood_for_filter{this, "min-ad-log-likelihood", "MmAaDdLliI", -9.0,
        "min log likelihood for AD assignments to PASS filter (0 for unlimited)"};
        
    Option<bool> write_trivial_calls{this, "trival", "ivtTIRV", false,
        "write trivial vcf calls (ex 0/0 genotypes)"};
        
    /// Should we call on nodes/edges outside of snarls by coverage (true), or
    /// just assert that primary path things exist and off-path things don't
    /// (false)?
    Option<bool> call_other_by_coverage{this, "call-nodes-by-coverage", "cCoObB", false,
        "make calls on nodes/edges outside snarls by coverage"};

    /// Use total support count (true) instead of total support quality (false) when choosing
    /// top alleles and deciding gentypes based on the biases.  
    Option<bool> use_support_count{this, "use-support-count", "T", false,
        "use total support count instead of total support quality for selecting top alleles"};

    /// Path of supports file generated from the PileupAugmenter (via vg augment)
    Option<string> support_file_name{this, "support-file", "s", {},
            "path of file containing supports generated by vg augment -P -s"};

    /// Don't collapse the shared prefix and suffix of the different alleles
    /// in an output VCF line.  This is mostly for debugging.
    Option<bool> leave_shared_ends{this, "leave-shared-ends", "X", false,
            "don't collapse shared prefix and suffix of alleles in VCF output"};

    Option<int> max_inversion_size{this, "max-inv", "e", 1000,
            "maximum detectable inversion size in number of nodes"};

    Option<string> recall_vcf_filename{this, "recall-vcf", "f", "",
            "VCF to genotype against.  Must have been used to create input graph with vg construct -a"};

    Option<string> recall_ref_fasta_filename{this, "recall-fasta", "a", "",
            "Reference FASTA required for --recall-vcf in the presence of symbolic deletions or inversions"};

    Option<string> recall_ins_fasta_filename{this, "insertion-fasta", "Z", "",
            "Insertion FASTA required for --recall-vcf in the presence of symbolic insertions"};

    /// Path of pack file generated from vg pack
    Option<string> pack_file_name{this, "pack-file", "P", {}, 
            "path of pack file from vg pack"};

    Option<string> xg_file_name{this, "xg-file", "x", {},
            "path of xg file (required to read pack file with -P)"};

    /// structures to hold the recall vcf and fastas
    vcflib::VariantCallFile variant_file;
    unique_ptr<FastaReference> ref_fasta;
    unique_ptr<FastaReference> ins_fasta;
    /// minimum average base support on alt path for it to be considered
    double min_alt_path_support = 0.2;
    
    /// print warnings etc. to stderr
    bool verbose = false;

    /// inversion or deletion edges greater than this length with 0 support
    /// will clamp average support down to 0.  this is primarily to prevent
    /// FP inversions when using average support
    int max_unsupported_edge_size = 20;
    
};

}

#endif
