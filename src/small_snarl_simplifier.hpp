#ifndef VG_SMALL_SNARL_SIMPLIFIER_HPP_INCLUDED
#define VG_SMALL_SNARL_SIMPLIFIER_HPP_INCLUDED


#include "progressive.hpp"
#include "vg.hpp"
#include <vg/vg.pb.h>
#include "traversal_finder.hpp"
#include "utility.hpp"
#include "path.hpp"
#include "feature_set.hpp"
#include "path_index.hpp"


/** \file 
 * Provides a class for simplifying graphs by popping bubbles.
 */
 
namespace vg {

using namespace std;

/**
 * A class that can be used to simplify a graph, by repeatedly popping leaf
 * bubbles under a certain size. Keeps graph paths and an optional set of BED-
 * like features up to date. TODO: doesn't handle path start and end positions
 * within nodes.
 */
class SmallSnarlSimplifier : public Progressive {

public:
    /// Make a simplifier that simplifies the given graph in place.
    SmallSnarlSimplifier(VG& graph);
    
    /// Simplify the graph by one step. Returns the number of nodes deleted and
    /// the number of edges deleted. Can be passed an iteration for its progress
    /// messages.
    pair<size_t, size_t> simplify_once(size_t iteration = 0);
    
    /// Simplify the graph until material stops being deleted or the maximum
    /// iteration count is reached.
    void simplify();
    
    /// What's the miniumum size of a bubble to keep, in involved bases?
    /// Everything smaller will get squished away.
    size_t min_size = 10;
    
    /// How many iterations of simplification should we allow in a simplify() call?
    size_t max_iterations = 10;
    
    /// Should we simplify bubbles where paths come in and leave through the
    /// enterance node (and delete those paths) (true)? Or should we leave those
    /// bubbles unsimplified?
    bool drop_hairpin_paths = false;
    
    /// If the user points this to a FeatureSet, that FeatureSet will get its
    /// features updated with changes to the graph as simplification proceeds.
    /// The user should load the features in and pull them out.
    /// TODO: Replace this with an on_path_edit event on this object that can be listened on.
    FeatureSet* features = nullptr;
    
protected:

    /// Holds a reference to the graph we're simplifying
    VG& graph;

    /// This keeps track of the sites to simplify
    SnarlManager site_manager;
    
    /// This is used to find traversals of those sites
    TrivialTraversalFinder traversal_finder;
    
    
};

}

#endif
