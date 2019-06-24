#ifndef VG_DECONSTRUCTOR_HPP_INCLUDED
#define VG_DECONSTRUCTOR_HPP_INCLUDED
#include <vector>
#include <string>
#include <sstream>
#include <ostream>
#include <sstream>
#include "genotypekit.hpp"
#include "path_index.hpp"
#include "Variant.h"
#include "path.hpp"
#include "vg.hpp"
#include "genotypekit.hpp"
#include "traversal_finder.hpp"
#include <vg/vg.pb.h>
#include "Fasta.h"

/** \file
* Deconstruct is getting rewritten.
* New functionality:
* -Detect superbubbles and bubbles
* -Fix command line interface.
* -harmonize on XG / raw graph (i.e. deprecate index)
* -Use unroll/DAGify if needed to avoid cycles

** Much of this is taken from Brankovic's
** "Linear-Time Superbubble Identification Algorithm for Genome Assembly"
*/
namespace vg{
using namespace std;

class Deconstructor{
public:

    Deconstructor();
    ~Deconstructor();

    // deconstruct the entire graph to cout
    void deconstruct(vector<string> refpaths, vg::VG* graph, SnarlManager* snarl_manager,
                     bool path_restricted_traversals,
                     const unordered_map<string, string>* path_to_sample = nullptr); 
    
private:

    // write a vcf record for the given site.  returns true if a record was written
    // (need to have a path going through the site)
    bool deconstruct_site(const Snarl* site);

    // convert traversals to strings.  returns mapping of traversal (offset in travs) to allele
    vector<int> get_alleles(vcflib::Variant& v, const vector<SnarlTraversal>& travs, int ref_path_idx, char prev_char);

    // write traversal path names as genotypes
    void get_genotypes(vcflib::Variant& v, const vector<string>& names, const vector<int>& trav_to_allele);

    // check to see if a snarl is too big to exhaustively traverse
    bool check_max_nodes(const Snarl* snarl);

    // get traversals from the exhaustive finder.  if they have nested visits, fill them in (exhaustively)
    // with node visits
    vector<SnarlTraversal> explicit_exhaustive_traversals(const Snarl* snarl);
    
    // output vcf object
    vcflib::VariantCallFile outvcf;
    
    // in memory path index for every specified reference path.
    map<string, PathIndex*> pindexes;

    // toggle between exhaustive and path restricted traversal finder
    bool path_restricted = false;

    // the graph
    VG* graph;

    // the snarl manager
    SnarlManager* snarl_manager;

    // the traversal finders. we always use a path traversal finder to get the reference path
    unique_ptr<PathRestrictedTraversalFinder> path_trav_finder;
    // we optionally use another (exhaustive for now) traversal finder if we don't want to rely on paths
    unique_ptr<TraversalFinder> trav_finder;

    // keep track of the non-ref paths as they will be our samples
    set<string> sample_names;

    // map the path name to the sample in the vcf
    const unordered_map<string, string>* path_to_sample;

    // upper limit of degree-2+ nodes for exhaustive traversal
    int max_nodes_for_exhaustive = 100;    
};

}
#endif
