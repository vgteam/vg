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

namespace vg {

using namespace std;

/**
 * Class to hold on to genotyping parameters and genotyping functions.
 */
class Genotyper {
public:

    // How many nodes max should we walk when checking if a path runs through a superbubble/site
    size_t max_path_search_steps = 100;
    
    // We need a single aligner object to make aligning multiple reads efficient
    // (at least when using quality-score-aware alignment)
    QualAdjAligner aligner = QualAdjAligner();
    
    /**
     * For the superbubble/site between start and end in the given orientations,
     * emit all unique subpaths that run from start to end, out of the paths in
     * the graph.
     */
    vector<Path> get_paths_through_site(VG& graph, NodeTraversal start, NodeTraversal end);
    
    /**
     * Get the affinity of all the reads relevant to the superbubble to all the
     * paths through the superbubble. We need to know all the nodes involved in
     * the superbubble so that we can clip them and their edges out and replace
     * them with the paths in turn.
     */ 
    map<Alignment*, vector<int>> get_affinities(VG& graph, const map<string, Alignment*>& reads_by_name,
        vector<id_t>& superbubble_contents, vector<Path>& superbubble_paths);
        
    
    

};

}


#endif
