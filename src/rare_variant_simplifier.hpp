#ifndef VG_RARE_VARIANT_SIMPLIFIER_HPP_INCLUDED
#define VG_RARE_VARIANT_SIMPLIFIER_HPP_INCLUDED


#include "progressive.hpp"
#include "vg.hpp"
#include "vcf_buffer.hpp"


/** \file 
 * Provides a class for simplifying graphs by removing rare variants.
 */
 
namespace vg {

using namespace std;

/**
 * A class that can be used to simplify a graph, by removing rare variants' alt
 * and ref paths and their exclusively-used nodes.
 */
class RareVariantSimplifier : public Progressive {

public:
    /// Make a simplifier that simplifies the given graph in place, using
    /// variants read using the given buffer.
    RareVariantSimplifier(MutablePathDeletableHandleGraph& graph, VcfBuffer& variant_source);
    
    /// Simplify the graph.
    void simplify();

    /// Keep variants at this total alt allele frequency or higher.
    double min_frequency_to_keep = 0;

    /// Keep variants with this total alt allele count or higher.
    /// AND'd with the frequency condition.
    size_t min_count_to_keep = 0;
     
protected:

    /// Holds a reference to the graph we're simplifying
    MutablePathDeletableHandleGraph& graph;

    /// Holds a reference to the variant buffer we are getting avriants from.
    VcfBuffer& variant_source;
};

}

#endif
