#ifndef VG_GENOTYPER_H
#define VG_GENOTYPER_H
// genotyper.hpp: defines the Genotyper, which is used to genotype from a graph
// and a collection of indexed alignments.

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <list>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"

namespace vg {

using namespace std;

/**
 * Holds indexes of the reference in a graph: position to node, node to position
 * and orientation, and the full reference string.
 */
struct ReferenceIndex {
    // Index from node ID to first position on the reference string and
    // orientation it occurs there.
    std::map<int64_t, std::pair<size_t, bool>> byId;
    
    // Index from start position on the reference to the oriented node that
    // begins there.  Some nodes may be backward (orientation true) at their
    // canonical reference positions. In this case, the last base of the node
    // occurs at the given position.
    std::map<size_t, vg::NodeTraversal> byStart;
    
    // The actual sequence of the reference.
    std::string sequence;
    
    // Make a ReferenceIndex from a path in a graph
    ReferenceIndex(VG& vg, string refPathName);
};

/**
 * Class to hold on to genotyping parameters and genotyping functions.
 */
class Genotyper {
public:

    // Represents a superbubble
    struct Site {
        NodeTraversal start;
        NodeTraversal end;
        set<id_t> contents;
    };

    // How many nodes max should we walk when checking if a path runs through a superbubble/site
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
    
    // If base qualities aren't available, what is the Phred-scale qualtiy of a
    // piece of sequence being correct?
    int default_sequence_quality = 15;
    
    // How mnay times must a path recur before we try aligning to it?
    // Note that the primary path counts as a recurrence.
    int min_recurrence = 2;
    
    // Remember that Phred scores can add to multiply just like normal log probs?
    // -10 log p1 + -10 log p2 = -10 (log p1 + log p2) = -10 log (p1 * p2)
    
    /**
     * Convert integer Phred score to probability.
     */
    double phred_to_prob(int phred);
     
    /**
     * Convert probability to integer Phred score.
     */
    int prob_to_phred(double prob);

    /**
     * Given an Alignment, compute a phred score for the quality of the
     * alignment's bases overall which is supposed to be interpretable as the
     * probability that the call of the sequence is wrong (to the degree that it
     * would no longer support the alleles it appears to support).
     *
     * In practice we're just going to average the quality scores for all the
     * bases.
     *
     * If the alignment doesn't have base qualities, returns 0.
     */
    int alignment_qual_score(const Alignment& alignment);

    /**
     * Unfold and dagify a graph, find the superbubbles, and then convert them
     * back to the space of the original graph.
     *
     * Returns a collection of Sites.
     */
    vector<Site> find_sites(VG& graph);

    /** 
     * Same as find_sites but use Cactus instead of Superbubbles.
     * This is more general and doesn't require DAGifcation etc., but we keep
     * both versions around for now for debugging and comparison
     */
    vector<Site> find_sites_with_cactus(VG& graph);
    
    /**
     * Given a path (which may run either direction through a site, or not touch
     * the ends at all), collect a list of NodeTraversals in order for the part
     * of the path that is inside the site, in the same orientation as the path.
     */
    list<NodeTraversal> get_traversal_of_site(VG& graph, const Site& site, const Path& path);
    
    /**
     * Make a list of NodeTraversals into the string they represent.
     */
    string traversals_to_string(const list<NodeTraversal>& path);
    
    /**
     * For the given site, emit all subpaths with unique sequences that run from
     * start to end, out of the paths in the graph.
     */
    vector<list<NodeTraversal>> get_paths_through_site(VG& graph, const Site& site);
    
    /**
     * Get all the quality values in the alignment between the start and end nodes
     * of a site. Handles alignments that enter the site from the end, and
     * alignments that never make it through the site.
     *
     * If we run out of qualities, or qualities aren't present, returns no
     * qualities.
     *
     * If an alignment goes through the site multipe times, we get all the qualities from when it is in the site.
     */
    string get_qualities_in_site(VG& graph, const Site& site, const Alignment& alignment);
    
    /**
     * Get the affinity of all the reads relevant to the superbubble to all the
     * paths through the superbubble.
     *
     * Affinity is a double out of 1.0. Higher is better.
     */ 
    map<Alignment*, vector<double>> get_affinities(VG& graph, const map<string, Alignment*>& reads_by_name,
        const Site& site,  const vector<list<NodeTraversal>>& superbubble_paths);
        
    /**
     * Get affinities as above but using only string comparison instead of
     * alignment. Affinities are 0 for mismatch and 1 for a perfect match.
     */
    map<Alignment*, vector<double>> get_affinities_fast(VG& graph, const map<string, Alignment*>& reads_by_name,
        const Site& site, const vector<list<NodeTraversal>>& superbubble_paths);
        
    /**
     * Compute annotated genotype from affinities and superbubble paths.
     * Needs access to the graph so it can chop up the alignments, which requires node sizes.
     */
    Genotype get_genotype(VG& graph, const Site& site, const vector<list<NodeTraversal>>& superbubble_paths, const map<Alignment*, vector<double>>& affinities);
        
    /**
     * Compute the probability of the observed alignments given the genotype.
     *
     * Takes a genotype as a vector of allele numbers, and support data as a
     * collection of pairs of Alignments and vectors of bools marking whether
     * each alignment is consistent with each allele.
     *
     * Alignments should have had their quality values trimmed down to just the
     * part covering the superbubble.
     */
    double get_genotype_probability(const vector<int>& genotype, const vector<pair<Alignment, vector<bool>>> alignment_consistency);
        
    /**
     * Make a VCFlib variant from a genotype. Depends on an index of the
     * reference path we want to call against.
     *
     * Returns 0 or more variants we can articulate from the superbubble.
     * Sometimes if we can't make a variant for the superbubble against the
     * reference path, we'll emit 0 variants.
     */
    vector<vcflib::Variant> genotype_to_variant(VG& graph, const ReferenceIndex& index, vcflib::VariantCallFile& vcf, const Genotype& genotype);
    
    /**
     * Make a VCF header
     */
    void write_vcf_header(std::ostream& stream, const std::string& sample_name, const std::string& contig_name, size_t contig_size);
    
    /**
     * Start VCF output to a stream. Returns a VCFlib VariantCallFile that needs to be deleted.
     */
    vcflib::VariantCallFile* start_vcf(std::ostream& stream, const ReferenceIndex& index, const string& ref_path_name);

};



}


#endif
