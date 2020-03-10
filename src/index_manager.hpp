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
#include "snarls.hpp"
#include "progressive.hpp"

namespace vg {

using namespace std;

/**
 * Represents a set of indexes (including the actual graph) organized around a base name.
 * The base name is the name of the indexed FASTA file.
 *
 * We have an internal collection of shared_ptrs, one for each final or
 * intermediate index.
 *
 * We also have a bunch of internal ensure_whatever() methods that make sure
 * that the corresponding index is populated, either from disk if available, or
 * from the indexes it depends on, first calling the ensure method for each.
 */
class IndexManager : public Progressive {
public:
    /*
     * Make a new IndexManager with the given FASTA providing the basename, and
     * using the variants from the given VCF if indexes need to be constructed.
     *
     * The fasta must be .fa or .fa.gz
     */
    IndexManager(const string& fasta_filename = "", const string& vcf_filename = "");

    /// Set the FASTA filename (and thus the basename for looking for other indexes.
    void set_fasta_filename(const string& filename);

    /// Set the VCF filename
    void set_vcf_filename(const string& filename);


    // Functions to set a source filename override, and get an index, for each index type.

    /// Override the file to load the minimizer index from
    void set_minimizer_override(const string& filename);
    /// Get the minimizer index
    shared_ptr<gbwtgraph::DefaultMinimizerIndex> get_minimizer();

    /// Override the file to load the gbwtgraph index from
    void set_gbwtgraph_override(const string& filename);
    /// Get the gbwtgraph index
    shared_ptr<gbwtgraph::GBWTGraph> get_gbwtgraph();

    /// Override the file to load the gbwt index from
    void set_gbwt_override(const string& filename);
    /// Get the gbwt index
    shared_ptr<gbwt::GBWT> get_gbwt();

    /// Override the file to load the distance index from
    void set_distance_override(const string& filename);
    /// Get the gbwt index
    shared_ptr<vg::MinimumDistanceIndex> get_distance();

    /// Override the file to load the snarls from
    void set_snarls_override(const string& filename);
    /// Get the snarls
    shared_ptr<vg::SnarlManager> get_snarls();

    /// Override the file to load the graph from
    void set_graph_override(const string& filename);
    /// Get the graph
    shared_ptr<PathHandleGraph> get_graph();
    
    /**
     * Get the indexes that are used for mapping. If not available, they will be generated.
     */
    tuple<gbwtgraph::GBWTGraph*, gbwt::GBWT*, gbwtgraph::DefaultMinimizerIndex*, vg::MinimumDistanceIndex*> get_mapping_indexes();

    /// Minimizer kmer length to use when minimizer indexing
    size_t minimizer_k = 29;
    /// Minimizer window size to use when minimizer indexing
    size_t minimizer_w = 11;

protected:

    // Save the input FASTA filename
    string fasta_filename;
    // Save the input VCF name
    string vcf_filename;
    // And the basename of the FASTA, which is where we expect to find our indexes.
    string basename;

    // Override filenames for loading all the indexes from
    string minimizer_override;
    string gbwtgraph_override;
    string gbwt_override;
    string distance_override;
    string snarls_override;
    string graph_override;

    // Store the final mapping indexes
    shared_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer;
    // GBWTGraph needs to keep a reference to the used GBWT alive, so it is stored in a different type.
    pair<shared_ptr<gbwtgraph::GBWTGraph>, shared_ptr<gbwt::GBWT>> gbwtgraph;
    shared_ptr<gbwt::GBWT> gbwt;
    shared_ptr<vg::MinimumDistanceIndex> distance;
    
    // And then the intermediate types: snarls and base non-GBWT graph
    shared_ptr<SnarlManager> snarls;
    shared_ptr<PathHandleGraph> graph;
    
    // And the functions to fill them in if empty.
    
    /// Load the graph, or make it from the FASTA and VCF file names, and save it to disk.
    void ensure_graph();
    
    /// Load the snarls (including trivial snarls), or make them from the graph and save them to disk.
    void ensure_snarls();
    
    /// Load the distance index, or make it from the graph and the snarls and save it to disk.
    void ensure_distance();
    
    /// Load the GBWT, or make it from the graph and the VCF filename and save it to disk.
    void ensure_gbwt();
    
    /// Load the GBWTGraph, or make it from the GBWT and the base graph, and save it to disk.
    void ensure_gbwtgraph();
    
    /// Load the minimizer index, or make it from the GBWTGraph and the GBWT and save it to disk.
    void ensure_minimizer();

    /// We have a template to help us stamp out these ensure functions.
    /// We define it in the CPP since only we ever use it.
    template<typename IndexHolderType>
    void ensure(IndexHolderType& member, const string& filename_override, const string& extension,
        const function<void(istream&)>& load, const function<void(ostream&)>& make_and_save);


    
    /// Get the filename for the index file having the given extension.
    /// Extension should not include the dot.
    /// May or may not exist yet.
    string get_filename(const string& extension) const;

    
    
};

}

#endif
