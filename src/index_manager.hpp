#ifndef VG_INDEX_MANAGER_HPP_INCLUDED
#define VG_INDEX_MANAGER_HPP_INCLUDED

/**
 * \file index_manager.hpp: defines a system for managing index files using a basename and constructing missing ones
 */

#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <map>

#include <gbwt/gbwt.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/minimizer.h>

#include "handle.hpp"
#include "min_distance.hpp"

namespace vg {

using namespace std;

/**
 * Represents a set of indexes (including the actual graph) organized around a base name.
 * The base name is the name of the indexed FASTA file.
 */
class IndexManager {
public:
    /*
     * Make a new IndexManager with the given FASTA providing the basename, and
     * using the variants from the given VCF if indexes need to be constructed.
     */
    IndexManager(const string& fasta_filename, const string& vcf_filename = "");
    
    /**
     * Get the indexes that are used for mapping. If not available, they will be generated.
     */
    tuple<gbwtgraph::GBWTGraph*, gbwt::GBWT*, gbwtgraph::DefaultMinimizerIndex*, vg::MinimumDistanceIndex*> get_mapping_indexes();

protected:

    // Store the final mapping indexes
    unique_ptr<gbwtgraph::GBWTGraph> gbwtgraph;
    unique_ptr<gbwt::GBWT> gbwt;
    unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer;
    unique_ptr<vg::MinimumDistanceIndex> distance;

    // For some stages we need the full graph
    unique_ptr<handlegraph::MutablePathMutableHandleGraph> graph;
    
};

}

#endif
