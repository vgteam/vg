#ifndef VG_VARIANT_ADDER_HPP_INCLUDED
#define VG_VARIANT_ADDER_HPP_INCLUDED

/** \file
 * variant_adder.hpp: defines a tool class used for adding variants from VCF
 * files into existing graphs.
 */

#include "vcf_buffer.hpp"
#include "path_index.hpp"
#include "handle.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include "progressive.hpp"
#include "name_mapper.hpp"
#include "graph_synchronizer.hpp"
#include "build_index.hpp"

namespace vg {

using namespace std;

/**
 * A tool class for adding variants to a VG graph. Integrated NameMapper
 * provides name translation for the VCF contigs.
 */
class VariantAdder : public NameMapper, public Progressive {

public:
    
    /**
     * Make a new VariantAdder to add variants to the given graph. Modifies the
     * graph in place.
     */
    VariantAdder(VG& graph);
    
    /**
     * Add in the variants from the given non-null VCF file. The file must be
     * freshly opened. The variants in the file must be sorted.
     *
     * May be called from multiple threads. Synchronizes internally on the
     * graph.
     */
    void add_variants(vcflib::VariantCallFile* vcf);
    
    /**
     * Align the given string to the given graph, between the given endpoints,
     * using the most appropriate alignment method, depending on the relative
     * sizes involved and whether a good alignment exists. max_span gives the
     * maximum length in the graph that we expect our string to possibly align
     * over (for cases of large deletions, where we might want to follow a long
     * path in the graph).
     *
     * The endpoints have to be heads/tails of the graph.
     *
     * Treats N/N substitutions as matches.
     *
     * TODO: now that we have a smart aligner that can synthesize deletions
     * without finding them with the banded global aligner, do we need max_span
     * anymore?
     *
     * Mostly exposed for testability.
     */
    Alignment smart_align(vg::VG& graph, pair<NodeSide, NodeSide> endpoints, const string& to_align, size_t max_span);
    
    /**
     * Turn any N/N substitutions in the given alignment against the given graph into matches.
     * Modifies the alignment in place.
     */
    void align_ns(vg::VG& graph, Alignment& aln);
    
    /// What's the maximum node size we should produce, and the size we should
    /// chop the input graph to? Since alt sequences are forced out to node
    /// boundaries, it makes sense for this to be small relative to
    /// whole_alignment_cutoff.
    size_t max_node_size = 32;
    
    /// How wide of a range in bases should we look for nearby variants in?
    size_t variant_range = 50;
    
    /// How much additional context should we try and add outside the radius of
    /// our group of variants we actually find?
    size_t flank_range = 100;
    
    /// Should we accept and ignore VCF contigs that we can't find in the graph?
    bool ignore_missing_contigs = false;
    
    /// What's the max radius on a variant we can have in order to use that
    /// variant as context for another main variant?
    size_t max_context_radius = 50;
    
    /// What's the cut-off for the graph's size or the alt's size in bp under
    /// which we can just use permissive banding and large band padding? If
    /// either is larger than this, we use the pinned-alignment-based do-each-
    /// end-and-splice mode.
    size_t whole_alignment_cutoff = 4096;
    
    /// When we're above that cutoff, what amount of band padding can we use
    /// looking for an existing version of our sequence?
    size_t large_alignment_band_padding = 30;
    
    /// When we're doing a restricted band padding alignment, how good does it
    /// have to be, as a fraction of the perfect match score for the whole
    /// context, in order to use it?
    double min_score_factor = 0.95;
    
    /// If the restricted band alignment doesn't find anything, we resort to
    /// pinned alignments from the ends and cutting and pasting together. How
    /// big should each pinned tail be?
    size_t pinned_tail_size = 200;
    
    /// We use this Aligner to hold the scoring parameters. It may be accessed
    /// by multiple threads at once.
    Aligner aligner;
    
    /// Sometimes, we have to make Mappers, for graphs too big to safely use our
    /// global banded aligner on. If we do that, what max edge crossing limit
    /// should we use for simplification?
    size_t edge_max = 0;
    /// What base kmer size should we use?
    size_t kmer_size = 16;
    /// What number of doubling steps should we use?
    size_t doubling_steps = 3;
    /// If nonzero, prune short subgraphs smaller than this before GCSA2-indexing
    size_t subgraph_prune = 0;
    
    // What's the cut-off above which we don't try thin banded aligners and do a
    // Mapper instead?
    size_t thin_alignment_cutoff = 10000;
    
    // What's the cutoff at which we decide to not even try the mapper-based
    // aligner (in query bases), and just resort to inserting duplicate
    // insertions?
    // TODO: This is disabled for now
    size_t mapper_alignment_cutoff = 0;
    
    /// We have code to skip large structural duplications, because aligners
    /// won't be able to distinguish the copies. TODO: we want to actually make
    /// them into cycles.
    bool skip_structural_duplications = true;
    
    /// Should we print out periodic updates about the variants we have
    /// processed?
    bool print_updates = false;
    
    const VG& get_graph() const;
    
protected:
    /// The graph we are modifying
    VG& graph;
    
    /// We keep a GraphSynchronizer so we can have multiple threads working on
    /// different parts of the same graph.
    GraphSynchronizer sync;
    
    /// We cache the set of valid path names, so we can detect/skip missing ones
    /// without locking the graph.
    set<string> path_names;
    
    /**
     * Get all the unique combinations of variant alts represented by actual
     * haplotypes. Arbitrarily phases unphased variants.
     *
     * Can (and should) take a WindowedVcfBuffer that owns the variants, and
     * from which cached pre-parsed genotypes can be extracted.
     *
     * Returns a set of vectors or one number per variant, giving the alt number
     * (starting with 0 for reference) that appears on the haplotype.
     *
     * TODO: ought to just take a collection of pre-barsed genotypes, but in an
     * efficient way (a vector of pointers to vectors of sample genotypes?)
     */
    set<vector<int>> get_unique_haplotypes(const vector<vcflib::Variant*>& variants, WindowedVcfBuffer* cache = nullptr) const;
    
    /**
     * Convert a haplotype on a list of variants into a string. The string will
     * run from the start of the first variant through the end of the last
     * variant.
     *
     * Can't be const because it relies on non-const operations on the
     * synchronizer.
     */
    string haplotype_to_string(const vector<int>& haplotype, const vector<vcflib::Variant*>& variants);
    
    /**
     * Get the radius of the variant around its center: the amount of sequence
     * that needs to be pulled out to make sure you have the ref and all the
     * alts, if they exist. This is just going to be twice the longest of the
     * ref and the alts.
     */
    static size_t get_radius(const vcflib::Variant& variant);
    
    /**
     * Get the center position of the given variant.
     */
    static size_t get_center(const vcflib::Variant& variant);
    
    /**
     * Get the center and radius around that center needed to extract everything
     * that might be involved in a group of variants.
     */
    static pair<size_t, size_t> get_center_and_radius(const vector<vcflib::Variant*>& variants);
    
    /**
     * Glom all the given variants into one vector, throwing out variants from
     * the before and after vectors that are too big to be in a context.
     */
    vector<vcflib::Variant*> filter_local_variants(const vector<vcflib::Variant*>& before,
        vcflib::Variant* variant, const vector<vcflib::Variant*>& after) const;
     
};

}

#endif

