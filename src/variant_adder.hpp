#ifndef VG_VARIANT_ADDER_HPP
#define VG_VARIANT_ADDER_HPP

/** \file
 * variant_adder.hpp: defines a tool class used for adding variants from VCF
 * files into existing graphs.
 */

#include "vcf_buffer.hpp"
#include "path_index.hpp"
#include "vg.hpp"
#include "progressive.hpp"
#include "name_mapper.hpp"
#include "graph_synchronizer.hpp"

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
    
    /// How wide of a range in bases should we look for nearby variants in?
    size_t variant_range = 50;
    
    /// How much additional context should we try and add outside the variants
    /// we grab? This sin't added to variant_range; it's counted from the outer
    /// ends of variants actually found.
    size_t flank_range = 20;
    
    /// Should we accept and ignore VCF contigs that we can't find in the graph?
    bool ignore_missing_contigs = false;
    
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
     */
    string haplotype_to_string(const vector<int>& haplotype, const vector<vcflib::Variant*>& variants);
 
};

}

#endif

