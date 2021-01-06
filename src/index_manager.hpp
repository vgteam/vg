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
#include "haplotype_indexer.hpp"

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
 *
 * For each index type (foo), we need:
 *
 * 1. A protected member "foo" which is a shared_ptr to the index.
 * 2. A protected member "foo_override" which is a filename to load the index
 *    from, to override the basename-based logic.
 * 3. A protected method "ensure_foo()" which uses "ensure()" to load or create
 *    the index.
 * 4. A public method "set_foo_override()" to set the override filename.
 * 5. A public method "get_foo()" to fill it in (via "ensure_foo()") and return it.
 * 6. A public method "can_get_foo()" which returns true if it thinks get_foo()
 *    will be successful, and false (and prints a warning) if e.g. necessary
 *    files are missing.
 * 7. Optionally, an sub-struct in the config struct with parameters for the
 *    indexing algorithm.
 *
 * Note that the IndexManager may exit the process, instead of throwing an
 * exception, if something goes wrong during a "get_foo()" method, so you
 * should check "can_get_foo()" first.
 *
 * For each supported tool (bar), we need:
 * 
 * 1. A public method "get_all_for_bar()", to fill in and save all indexes the
 *    tool uses.
 * 2. A public methos "can_get_all_for_bar()", which returns true if it thinks
 *    get_all_for_bar() will be successful, and false (and prints a warning) if
 *    any required index has its can_get_foo() method fail.
 */
class IndexManager : public Progressive {
public:
    /**
     * Make a new IndexManager with the given FASTA providing the basename, and
     * using the variants from the given VCF if indexes need to be constructed.
     *
     * The fasta must be .fa or .fa.gz
     */
    IndexManager(const string& fasta_filename = "", const string& vcf_filename = "");
    
    /// For GBWTs, where should threads come from?
    // TODO: Not all of these are implemented here yet.
    enum thread_source_type { thread_source_default, thread_source_vcf, thread_source_paths, thread_source_gam, thread_source_gaf };
    
    /**
     * Configurations for each of the indexes to be generated.
     */
    struct {
        /// Configuration for the minimizer index
        struct {
            /// Minimizer kmer length to use when minimizer indexing
            size_t k = 29;
            /// Minimizer window size to use when minimizer indexing
            size_t w = 11;
            /// Syncmer smer length to use when minimizer indexing
            size_t s = 18;
        } minimizer;
        /// Configuration for GBWT (doubles as the haplotype indexer widget)
        struct GBWTConfig : public HaplotypeIndexer {
            /// Where should threads come from when generating the GBWT?
            thread_source_type thread_source = thread_source_default;
        } gbwt;
    } config;

    /// Set the FASTA filename (and thus the basename for looking for other indexes, if not already set).
    void set_fasta_filename(const string& filename);

    /// Set the VCF filename
    void set_vcf_filename(const string& filename);
    
    // Functions to get indexes for a particular tool
    
    /// Get all indexes needed for the Giraffe mapper (MinimizerMapper).
    /// Returns only the final indexes, not intermediates.
    tuple<
        shared_ptr<gbwtgraph::DefaultMinimizerIndex>,
        shared_ptr<gbwtgraph::GBWTGraph>,
        shared_ptr<gbwt::GBWT>,
        shared_ptr<vg::MinimumDistanceIndex>> get_all_for_giraffe();
    /// Returns true if the indexes needed for Giraffe are available or can be
    /// generated/loaded, and false otherwise.
    bool can_get_all_for_giraffe();


    // Functions to set a source filename override, and get an index, for each index type.

    /// Override the file to load the minimizer index from
    void set_minimizer_override(const string& filename);
    /// Get the minimizer index
    shared_ptr<gbwtgraph::DefaultMinimizerIndex> get_minimizer();
    /// Returns true if the minimizer index is available or can be generated/loaded, and false otherwise.
    bool can_get_minimizer() const;

    /// Override the file to load the gbwtgraph index from
    void set_gbwtgraph_override(const string& filename);
    /// Get the gbwtgraph index
    shared_ptr<gbwtgraph::GBWTGraph> get_gbwtgraph();
    /// Returns true if the GBWTGraph is available or can be generated/loaded, and false otherwise.
    bool can_get_gbwtgraph() const;

    /// Override the file to load the gbwt index from
    void set_gbwt_override(const string& filename);
    /// Get the gbwt index
    shared_ptr<gbwt::GBWT> get_gbwt();
    /// Returns true if the GBWT is available or can be generated/loaded, and false otherwise.
    bool can_get_gbwt() const;

    /// Override the file to load the distance index from
    void set_distance_override(const string& filename);
    /// Get the gbwt index
    shared_ptr<vg::MinimumDistanceIndex> get_distance();
    /// Returns true if the distance index is available or can be generated/loaded, and false otherwise.
    bool can_get_distance() const;

    /// Override the file to load the snarls from
    void set_snarls_override(const string& filename);
    /// Get the snarls
    shared_ptr<vg::SnarlManager> get_snarls();
    /// Returns true if the snarls are available or can be generated/loaded, and false otherwise.
    bool can_get_snarls() const;

    /// Override the file to load the graph from.
    /// Also sets index basename to be based on this graph file.
    void set_graph_override(const string& filename);
    /// Get the graph
    shared_ptr<PathHandleGraph> get_graph();
    /// Returns true if the graph is available or can be generated/loaded, and false otherwise.
    bool can_get_graph() const;
    
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
    /// Note that the ostream to make_and_save is only open if there is a
    /// basename and a file to write to.
    template<typename IndexHolderType>
    void ensure(IndexHolderType& member, const string& filename_override, const string& extension,
        const function<void(ifstream&)>& load, const function<void(ofstream&)>& make_and_save);
        
    
    /// We have a template for helping write the can_get functions. Defined in
    /// the CPP since only we use it.
    template<typename IndexHolderType>
    bool can_get(IndexHolderType& member, const string& filename_override, const string& extension,
        const function<bool(void)>& poll_dependencies) const;

    /// Get the filename for the index file having the given extension.
    /// Extension should not include the dot.
    /// May or may not exist yet.
    /// Returns "" if there is no basename.
    string get_filename(const string& extension) const;

    
    
};

}

#endif
