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
    RareVariantSimplifier(VG& graph, VcfBuffer& variant_source);
    
    /// Simplify the graph.
    void simplify();
     
protected:

    /// Holds a reference to the graph we're simplifying
    VG& graph;

    /// Holds a reference to the variant buffer we are getting avriants from.
    VcfBuffer& variant_source;
};

}

#endif
