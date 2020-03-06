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
#include <gbwt/variants.h>
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


    /**
     * Generate a GBWT index.
     */
    static unique_ptr<gbwt::GBWT> make_gbwt(PathHandleGraph* xg_index, bool index_paths, const vector<string>& gam_file_names,  bool show_progress);

    /**
     * Parse a VCF file into the types needed for GBWT indexing.
     *
     * Takes a graph with alt paths embedded in it, and the corresponding VCF.
     * 
     * Path names in the graph are mapped to VCF contig names via path_to_vcf,
     * or used as-is if no entry there is found.
     *
     * If file_prefix is nonempty, a file for each contig is saved to
     * PREFIX_VCFCONTIG, and files for each batch of haplotypes are saved to
     * files named like PREFIX_VCFCONTIG_STARTSAMPLE_ENDSAMPLE. Otherwise, the
     * batch files are still saved, but to temporary files.
     *
     * Calls the callback with each contig's gbwt::VariantPaths, for each
     * gbwt::PhasingInformation batch of samples.
     */
    static void parse_vcf(const PathHandleGraph* graph, const string& vcf_filename, const map<string, string>& path_to_vcf, const string& batch_file_prefix, const function<void(const gbwt::VariantPaths&, const gbwt::PhasingInformation&)>& handle_contig_haplotype_batch);

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
