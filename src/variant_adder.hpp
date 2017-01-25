#ifndef VG_VARIANT_ADDER_HPP
#define VG_VARIANT_ADDER_HPP

/** \file
 * variant_adder.hpp: defines a tool class used for adding variants from VCF
 * files into existing graphs.
 */

#include "vcf_buffer.hpp"
#include "vg.hpp"

namespace vg {

using namespace std;

/**
 * A tool class for adding variants to a VG graph.
 */
class VariantAdder {

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
     * Each file of variants is added as a batch.
     */
    void add_variants(vcflib::VariantCallFile* vcf);
    
    /// How wide of a range in bases should we look for nearby variants in?
    size_t variant_range = 100;
    
protected:
    /// The graph we are modifying
    VG& graph;
    
    /**
     * Get all the unique combinations of variant alts represented by actual
     * haplotypes. Arbitrarily phases unphased variants.
     *
     * Returns a set of vectors or one number per variant, giving the alt number
     * (starting with 0 for reference) that appears on the haplotype.
     */
    set<vector<int>> get_unique_haplotypes(vector<vcflib::Variant*>& variants);

};

}

#endif

