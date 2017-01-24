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
    
protected:
    /// The graph we are modifying
    VG& graph;

};

}

#endif

