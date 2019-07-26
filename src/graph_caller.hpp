#ifndef VG_GRAPH_CALLER_HPP_INCLUDED
#define VG_GRAPH_CALLER_HPP_INCLUDED

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
#include "snarl_caller.hpp"
#include "region.hpp"

namespace vg {

using namespace std;

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager,
                ostream& out_stream = cout);

    virtual ~GraphCaller();

    /// Run call_snarl() on every top-level snarl in the manager.
    /// For any that return false, try the children, etc. (when recurse_on_fail true)
    /// Snarls are processed in parallel
    virtual void call_top_level_snarls(bool recurse_on_fail = true);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs) const;
   
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// Where all output written
    ostream& out_stream;

};

/**
 * VCFGenotyper : Genotype variants in a given VCF file
 */
class VCFGenotyper : public GraphCaller {
public:
    VCFGenotyper(const PathHandleGraph& graph,
                 SnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 vcflib::VariantCallFile& variant_file,
                 const string& sample_name,
                 const vector<string>& ref_paths = {},
                 FastaReference* ref_fasta = nullptr,
                 FastaReference* ins_fasta = nullptr,
                 ostream& out_stream = cout);

    virtual ~VCFGenotyper();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs) const;

protected:


protected:

    /// the graph
    const PathHandleGraph& graph;

    /// input VCF to genotype, must have been loaded etc elsewhere
    vcflib::VariantCallFile& input_vcf;

    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;    

    /// traversal finder uses alt paths to map VCF alleles from input_vcf
    /// back to traversals in the snarl
    VCFTraversalFinder traversal_finder;

};


}

#endif
