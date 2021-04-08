#ifndef VG_GENOTYPER_HPP_INCLUDED
#define VG_GENOTYPER_HPP_INCLUDED
// genotyper.hpp: defines the Genotyper, which is used to genotype from a graph
// and a collection of indexed alignments.

#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <vector>
#include <list>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "hash_map.hpp"
#include "types.hpp"
#include "genotypekit.hpp"
#include "statistics.hpp"

namespace vg {

using namespace std;

/**
 * Class to hold on to genotyping parameters and genotyping functions.
 */
class Genotyper {
public:

    // Represents an assertion that a read is consistent or not consistent with
    // an allele, including an affinity weight and a flag for strand
    struct Affinity {
        // Is the read consistent with this allele?
        bool consistent = false;
        // How consistent is the read with the allele, out of 1.
        double affinity = 0;
        // Is the read on the forward strand (false) or reverse strand (true)
        bool is_reverse = false;
        // What's the actual score (not necessarily normalized out of 1)?
        // We'll probably put per-base alignment score here
        double score = 0;    
        
        // What's the unnormalized log likelihood of the read given the allele?
        double likelihood_ln = 0;    
        
        // Have a default constructor
        Affinity() = default;
        
        // Have a useful constructor
        Affinity(double affinity, bool is_reverse) : consistent(affinity == 1), 
            affinity(affinity), is_reverse(is_reverse), score(affinity) {
            // Nothing to do
        }
        
    };


    // How many nodes max should we walk when checking if a path runs through a superbubble/snarl
    size_t max_path_search_steps = 100;
    
    // How long should we unfold graphs to?
    int unfold_max_length = 200;
    
    // How many steps of dagification should we do?
    int dagify_steps = 1;
    
    // What portion of the total affinity can be on one allele before we refuse
    // to call het and call homozygous instead.
    double max_het_bias = 4;
    
    // Should we use mapping quality when computing P(reads | genotype)?
    bool use_mapq = false;
    
    // Whould we do indel realignment, or should we use fast substring
    // affinities for everything?
    bool realign_indels = false;
    
    // If base qualities aren't available, what is the Phred-scale qualtiy of a
    // piece of sequence being correct?
    int default_sequence_quality = 15;
    
    // How many times must a path recur before we try aligning to it? Also, how
    // many times must a node in the graph be visited before we use it in indel
    // realignment for nearby indels? Note that the primary path counts as a
    // recurrence. TODO: novel inserts can't recur, and novel deletions can't be
    // filtered in this way.
    int min_recurrence = 2;
    
    // How much unique support must an alt have on each strand before we can
    // call it? Calls that fail this get dropped from the VCF.
    // TODO: use a VCF filter instead.
    int min_unique_per_strand = 2;
    
    // When we realign reads, what's the minimum per-base score for a read in
    // order to actually use it as supporting the thing we just aligned it to?
    double min_score_per_base = 0.90;
    
    // Now define the prior distribution on genotypes.
    
    // What should the prior probability of a snarl being diploid be?
    double diploid_prior_logprob = prob_to_logprob(0.999);
    // What should our prior on being heterozygous at a snarl be, given that it is diploid?
    double het_prior_logprob = prob_to_logprob(0.1);
    // What should the prior probability of a non-diploid snarl being haploid be?
    double haploid_prior_logprob = prob_to_logprob(0.5);
    // What should the prior probability of a non-diploid, non-haploid snarl
    // being completely deleted be?
    double deleted_prior_logprob = prob_to_logprob(0.1);
    // If not diploid, not haploid, and not deleted, the snarl is polyploid. We
    // model the number of additional copies over 2 with a geometric prior. What
    // should the parameter (success probability for stopping adding copies) be?
    double polyploid_prior_success_logprob = prob_to_logprob(0.5);

    // Provides a mechanism to translate back to the original graph
    Translator translator;
    
    
    // We need some data structures to track statistics about snarls and their
    // traversals. These data structures are only valid during the run method.
    // After that the snarl pointers they use may be left dangling.
    
    // This maps from length in the reference to number of snarls of that length.
    unordered_map<size_t, size_t> snarl_reference_length_histogram;
    
    // Which snarls have been traversed by at least one read?
    unordered_set<const Snarl*> snarl_traversals;
    
    // What snarls even exist?
    unordered_set<const Snarl*> all_snarls;
    
    // We need to have aligners in our genotyper, for realigning around indels.
    vector<Aligner> normal_aligners;
    vector<QualAdjAligner> quality_aligners;

    // Toggle traversal finder for testing
    enum TraversalAlg { Reads, Exhaustive, Representative, Adaptive };
    map<TraversalAlg, string> alg2name = {
        {Reads, "Reads"},
        {Exhaustive, "Exhaustive"},
        {Representative, "Representative"},
        {Adaptive, "Adaptive"}
    };
    TraversalAlg traversal_alg = TraversalAlg::Reads;

    // Show progress
    bool show_progress = false;

    /// Process and write output.
    /// Alignments must be embedded in the AugmentedGraph.
    void run(AugmentedGraph& graph,
             ostream& out,
             string ref_path_name,
             string contig_name = "",
             string sample_name = "",
             string augmented_file_name = "",
             bool output_vcf = false,
             bool output_json = false,
             int length_override = 0,
             int variant_offset = 0);
    
    /**
     * Given an Alignment and a Snarl, compute a phred score for the quality of
     * the alignment's bases within the snarl overall (not counting the start and
     * end nodes), which is supposed to be interpretable as the probability that
     * the call of the sequence is wrong (to the degree that it would no longer
     * support the alleles it appears to support).
     *
     * In practice we're just going to average the quality scores for all the
     * bases interior to the snarl (i.e. not counting the start and end nodes).
     *
     * If the alignment doesn't have base qualities, or no qualities are
     * available for bases internal to the snarl, returns a default value.
     */
    int alignment_qual_score(VG& graph, const Snarl* snarl, const Alignment& alignment);

    /**
     * Check if a mapping corresponds to the beginning or end of snarl by making
     * sure it crosses the given side in the expected direction. The handle
     * should be forward for the left side and reverse for the right side.
     */
    static bool mapping_enters_side(const Mapping& mapping, const handle_t& side, const HandleGraph* graph);
    static bool mapping_exits_side(const Mapping& mapping, const handle_t& side, const HandleGraph* graph);

    /**
     * Check if a snarl is small enough to be covered by reads (very conservative)
     */ 
    static bool is_snarl_smaller_than_reads(AugmentedGraph& augmented_graph,
                                            const Snarl* snarl,
                                            const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                            map<string, const Alignment*>& reads_by_name);
        
    /**
     * Get traversals of a snarl in one of several ways. 
     */
    vector<SnarlTraversal> get_snarl_traversals(AugmentedGraph& augmented_graph, SnarlManager& manager,
                                                map<string, const Alignment*>& reads_by_name,
                                                const Snarl* snarl,
                                                const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                                PathIndex* reference_index,
                                                TraversalAlg use_traversal_alg);
    
    /**
     * Get all the quality values in the alignment between the start and end
     * nodes of a snarl. Handles alignments that enter the snarl from the end, and
     * alignments that never make it through the snarl.
     *
     * If we run out of qualities, or qualities aren't present, returns no
     * qualities.
     *
     * If an alignment goes through the snarl multipe times, we get all the
     * qualities from when it is in the snarl.
     *
     * Does not return qualities on the start and end nodes. May return an empty
     * string.
     */
    string get_qualities_in_snarl(VG& graph, const Snarl* snarl, const Alignment& alignment);
    
    /**
     * Get the affinity of all the reads relevant to the superbubble to all the
     * paths through the superbubble.
     *
     * Affinity is a double out of 1.0. Higher is better.
     */ 
    map<const Alignment*, vector<Affinity>> get_affinities(AugmentedGraph& aug,
                                                           const map<string, const Alignment*>& reads_by_name,
                                                           const Snarl* snarl,
                                                           const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                                           const SnarlManager& manager,
                                                           const vector<SnarlTraversal>& superbubble_paths);
        
    /**
     * Get affinities as above but using only string comparison instead of
     * alignment. Affinities are 0 for mismatch and 1 for a perfect match.
     */
    map<const Alignment*, vector<Affinity>> get_affinities_fast(AugmentedGraph& aug,
                                                                const map<string, const Alignment*>& reads_by_name,
                                                                const Snarl* snarl,
                                                                const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                                                const SnarlManager& manager,
                                                                const vector<SnarlTraversal>& superbubble_paths,
                                                                bool allow_internal_alignments = false);
        
    /**
     * Compute annotated genotype from affinities and superbubble paths.
     * Needs access to the graph so it can chop up the alignments, which requires node sizes.
     */
    Locus genotype_snarl(VG& graph, const Snarl* snarl, const vector<SnarlTraversal>& superbubble_paths, const map<const Alignment*, vector<Affinity>>& affinities);
        
    /**
     * Compute the probability of the observed alignments given the genotype.
     *
     * Takes a genotype as a vector of allele numbers, and support data as a
     * collection of pairs of Alignments and vectors of bools marking whether
     * each alignment is consistent with each allele.
     *
     * Alignments should have had their quality values trimmed down to just the
     * part covering the superbubble.
     *
     * Returns a natural log likelihood.
     */
    double get_genotype_log_likelihood(VG& graph, const Snarl* snarl, const vector<int>& genotype, const vector<pair<const Alignment*, vector<Affinity>>>& alignment_consistency);
    
    /**
     * Compute the prior probability of the given genotype.
     *
     * Takes a genotype as a vector of allele numbers. It is not guaranteed that
     * allele 0 corresponds to any notion of primary reference-ness.
     *
     * Returns a natural log prior probability.
     *
     * TODO: add in strand bias
     */
    double get_genotype_log_prior(const vector<int>& genotype);
        
    /**
     * Make a VCFlib variant from a called Locus. Depends on an index of the
     * reference path we want to call against.
     *
     * Returns 0 or more variants we can articulate from the superbubble.
     * Sometimes if we can't make a variant for the superbubble against the
     * reference path, we'll emit 0 variants.
     */
    vector<vcflib::Variant> locus_to_variant(VG& graph, const Snarl* snarl,
                                             const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                                             const SnarlManager& manager,
                                             const PathIndex& index, vcflib::VariantCallFile& vcf,
                                             const Locus& locus,
                                             const string& sample_name = "SAMPLE");
    
    /**
     * Make a VCF header
     */
    void write_vcf_header(std::ostream& stream, const std::string& sample_name, const std::string& contig_name, size_t contig_size);
    
    /**
     * Start VCF output to a stream. Returns a VCFlib VariantCallFile that needs to be deleted.
     */
    vcflib::VariantCallFile* start_vcf(std::ostream& stream, const PathIndex& index, const string& sample_name, const string& contig_name, size_t contig_size);
    
    /**
     * Utility function for getting the reference bounds (start and past-end) of
     * a snarl with relation to a given reference index in the given graph.
     * Computes bounds of the variable region, not including the fixed start and
     * end node lengths. Also returns whether the reference path goes through
     * the snarl forwards (false) or backwards (true).
     */
    pair<pair<int64_t, int64_t>, bool> get_snarl_reference_bounds(const Snarl* snarl, const PathIndex& index,
        const HandleGraph* graph);
    
    /**
     * Tell the statistics tracking code that a snarl exists. We can do things
     * like count up the snarl length in the reference and so on. Called only
     * once per snarl, but may be called on multiple threads simultaneously.
     */
    void report_snarl(const Snarl* snarl, const SnarlManager& manager, const PathIndex* index, VG& graph,
                      PathIndex* reference_index);
    
    /**
     * Tell the statistics tracking code that a read traverses a snarl
     * completely. May be called multiple times for a given read and snarl, and
     * may be called in parallel.
     */
    void report_snarl_traversal(const Snarl* snarl, const SnarlManager& manager, VG& graph);

    /**
     * Print some information about affinities
     */
    void report_affinities(map<const Alignment*, vector<Genotyper::Affinity>>& affinities,
                           vector<SnarlTraversal>& paths, VG& graph);

    /**
     * Write the variant as output
     */
    void write_variant();
    
    /**
     * Print snarl statistics to the given stream.
     */
    void print_statistics(ostream& out);

    /**
     * Subset the graph
     */
    VG make_subset_graph(VG& graph, const string& ref_path_name,
                         map<string, const Alignment*>& reads_by_name); 
      

    /*
     *
     *
     *
     * Experimental alternative genotyper
     *
     *
     *
     *
     */
    
    // hash function for node traversals
    struct hash_node_traversal {
        size_t operator()(const NodeTraversal& node_traversal) const {
            return (size_t) 1099511628211ull * ((uintptr_t) node_traversal.node + node_traversal.backward * 16777619ull) + 14695981039346656037ull;
        }
    };
    // hash function for oriented edges (as pairs of node traversals)
    struct hash_oriented_edge {
        size_t operator()(const pair<const NodeTraversal, const NodeTraversal>& edge) const {
            hash_node_traversal hsh;
            return hsh(edge.first) ^ (hsh(edge.second) << 1);
        }
    };
    
    // hash function for ambiguous sets of alleles that alleles can match (implemented as sorted arrays of indices)
    struct hash_ambiguous_allele_set {
        size_t operator()(const vector<size_t>& ambiguous_set) const {
            size_t hash_value = 0;
            for (size_t i = 0; i < ambiguous_set.size(); i++) {
                hash_value += 1099511628211ull * ambiguous_set[i] + 14695981039346656037ull;
            }
            return hash_value;
        }
    };
    
    void edge_allele_labels(const VG& graph,
                            const Snarl* snarl,
                            const vector<list<NodeTraversal>>& superbubble_paths,
                            unordered_map<pair<NodeTraversal, NodeTraversal>,
                                          unordered_set<size_t>,
                                          hash_oriented_edge>* out_edge_allele_sets);
    void allele_ambiguity_log_probs(const VG& graph,
                                    const Snarl* snarl,
                                    const vector<list<NodeTraversal>>& superbubble_paths,
                                    const unordered_map<pair<NodeTraversal, NodeTraversal>,
                                                        unordered_set<size_t>,
                                                        hash_oriented_edge>& edge_allele_sets,
                                    vector<unordered_map<vector<size_t>,
                                                         double,
                                                         hash_ambiguous_allele_set>>* out_allele_ambiguity_probs);
};



}


#endif
