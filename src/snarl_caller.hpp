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
#include "traversal_finder.hpp"
#include "traversal_genotyper.hpp"

namespace vg {

using namespace std;

/**
 * SnarlCaller: Use the snarl decomposition to call snarls in a graph
 */
class SnarlCaller {
public:
    SnarlCaller(TraversalGenotyper& traversal_genotyper,
                SnarlManager& snarl_manager,
                ostream& out_stream);

    virtual ~SnarlCaller();

    /// Run call_snarl() on every top-level snarl in the manager.
    /// For any that return false, try the children, etc. (when recurse_on_fail true)
    /// Snarls are processed in parallel
    virtual void call_top_level_snarls(bool recurse_on_fail = true);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

protected:

    /// Write the vcf header
    virtual void header() = 0;
   
protected:

    /// Our snarls
    SnarlManager& snarl_manager;

    /// Where all output written
    ostream& out_stream;

    /// Our Genotyper
    TraversalGenotyper& traversal_genotyper;
};

/**
 * VCFGenotyper : Genotype variants in a given VCF file
 */
class VCFGenotyper : public SnarlCaller {
public:
    VCFGenotyper(TraversalGenotyper& traversal_genotyper,
                 SnarlManager& snarl_manager,
                 ostream& out_stream,
                 vcflib::VariantCallFile& variant_file);

    virtual ~VCFGenotyper();

    virtual bool call_snarl(const Snarl& snarl);

protected:
    
    virtual void header();

protected:

    /// input VCF to genotype, must have been loaded etc elsewhere
    vcflib::VariantCallFile& input_vcf;

    /// object that will take care of choosing our best traversals
    TraversalGenotyper& traversal_genotyper;

    /// output vcf
    vcflib::VariantCallFile output_vcf;

    /// traversal finder uses alt paths to map VCF alleles from input_vcf
    /// back to traversals in the snarl
    VCFTraversalFinder traversal_finder;
};


}

#endif
