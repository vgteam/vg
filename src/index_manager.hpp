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
#include "progressive.hpp"

namespace vg {

using namespace std;

/**
 * Allows indexing haplotypes, either to pre-parsed haplotype files or to a GBWT.
 */
class HaplotypeIndexer : public Progressive {
public:

    /// Print a warning if variants in the VCF can't be found in the graph
    bool warn_on_missing_variants = true;
    /// Track the number of variants in the phasing VCF that aren't found in the graph
    /// TODO: Make atomic?
    size_t found_missing_variants = 0; 
    /// Only report up to this many of them
    size_t max_missing_variant_warnings = 10; 

    /// Path names in the graph are mapped to VCF contig names via path_to_vcf,
    /// or used as-is if no entry there is found.
    map<string, string> path_to_vcf;
    
    /// Use graph path names instead of VCF path names when composing variant
    /// alt paths.
    bool rename_variants = true;
    
    /// If batch_file_prefix is nonempty, a file for each contig is saved to
    /// PREFIX_VCFCONTIG, and files for each batch of haplotypes are saved to
    /// files named like PREFIX_VCFCONTIG_STARTSAMPLE_ENDSAMPLE. Otherwise, the
    /// batch files are still saved, but to temporary files.
    string batch_file_prefix = "";
    
    /// If set to true, store paths from the graph alognside haplotype threads
    /// from the VCF, if any.
    bool index_paths = false;

    /// Phase homozygous unphased variants
    bool phase_homozygous = true;
    
    /// Arbitrarily phase all unphased variants
    bool force_phasing = false;
    
    /// Join together overlapping haplotypes
    bool discard_overlaps = false;

    /// Number of samples to process together in a haplotype batch.
    size_t samples_in_batch = 200;
    
    /// Range of VCF samples to process (first to past-last).
    pair<size_t, size_t> sample_range = pair<size_t, size_t>(0, numeric_limits<size_t>::max());
    
    /// Region restrictions for contigs, in VCF name space, as 0-based
    /// exclusive-end ranges.
    map<string, pair<size_t, size_t>> regions;

    /**
     * Parse a VCF file into the types needed for GBWT indexing.
     *
     * Takes a graph, a map of alt paths by name that will be extracted from
     * the graph if not populated, a vector of contigs in the graph to process,
     * in order, and the corresponding VCF file, already open. Sample parsing
     * on the VCF file should be turned off.
     *
     * Uses the given vector of sample names, which must be pre-populated with
     * any other samples (such as "ref") already in a GBWTBuilder that you are
     * using with the results of this function. Sample names from the VCF will
     * be added.
     *
     * Calls the callback with the contig number, each contig's
     * gbwt::VariantPaths, for each gbwt::PhasingInformation batch of samples.
     * The gbwt::PhasingInformation is not const because the GBWT library needs
     * to modify it in order to generate haplotypes from it efficiently.
     *
     * Doesn't create threads for embedded graph paths itself.
     *
     * Returns the number of haplotypes created (2 per sample) This number will
     * need to be adjusted if any samples' haplotypes are filtered out later.
     * This function ignores any sample filters and processes the entire VCF.
     */
    size_t parse_vcf(const PathHandleGraph* graph, map<string, Path>& alt_paths, const vector<path_handle_t>& contigs, vcflib::VariantCallFile& variant_file, std::vector<std::string>& sample_names, const function<void(size_t, const gbwt::VariantPaths&, gbwt::PhasingInformation&)>& handle_contig_haplotype_batch);
    
    /**
     * Generate a GBWT index.
     */
    static unique_ptr<gbwt::GBWT> make_gbwt(const PathHandleGraph* graph, bool index_paths, const string& vcf_filename, const vector<string>& gam_filenames);
};

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
