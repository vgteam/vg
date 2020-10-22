#ifndef VG_HAPLOTYPE_INDEXER_HPP_INCLUDED
#define VG_HAPLOTYPE_INDEXER_HPP_INCLUDED

/**
 * \file haplotype_indexer.hpp: defines how to construct GBWT indexes and VCF parses.
 */

#include <string>
#include <utility>
#include <map>
#include <memory>
#include <unordered_set>
#include <vector>

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/variants.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/minimizer.h>


#include "handle.hpp"
#include "progressive.hpp"

namespace vg {

using namespace std;

/**
 * Allows indexing haplotypes, either to pre-parsed haplotype files or to a GBWT.
 */
class HaplotypeIndexer : public Progressive {
public:

    /// Treat the embedded paths in the graph as samples. By default,
    /// the paths are interpreted as contigs.
    bool paths_as_samples = false;

    /// Print a warning if variants in the VCF can't be found in the graph
    bool warn_on_missing_variants = true;

    /// Only report up to this many of them
    size_t max_missing_variant_warnings = 10; 

    /// Path names in the graph are mapped to VCF contig names via path_to_vcf,
    /// or used as-is if no entry there is found.
    std::map<std::string, std::string> path_to_vcf;
    
    /// Use graph path names instead of VCF path names when composing variant
    /// alt paths.
    bool rename_variants = true;
    
    /// If batch_file_prefix is nonempty, a file for each contig is saved to
    /// PREFIX_VCFCONTIG, and files for each batch of haplotypes are saved to
    /// files named like PREFIX_VCFCONTIG_STARTSAMPLE_ENDSAMPLE. Otherwise, the
    /// batch files are still saved, but to temporary files.
    std::string batch_file_prefix = "";

    /// Phase homozygous unphased variants
    bool phase_homozygous = true;
    
    /// Arbitrarily phase all unphased variants
    bool force_phasing = false;
    
    /// Join together overlapping haplotypes
    bool discard_overlaps = false;

    /// Number of samples to process together in a haplotype batch.
    size_t samples_in_batch = 200;
    
    /// Size of the GBWT buffer in millions of nodes
    size_t gbwt_buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION;
    
    /// Interval at which to sample for GBWT locate
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
    
    /// Range of VCF samples to process (first to past-last).
    std::pair<size_t, size_t> sample_range = std::pair<size_t, size_t>(0, std::numeric_limits<size_t>::max());
    
    /// Region restrictions for contigs, in VCF name space, as 0-based
    /// exclusive-end ranges.
    std::map<std::string, std::pair<size_t, size_t>> regions;
    
    /// Excluded VCF sample names, for which threads will not be generated.
    /// Ignored during VCF parsing.
    std::unordered_set<std::string> excluded_samples;
    
    /// Perform initialization of backing libraries
    HaplotypeIndexer();

    /**
     * Parse a VCF file into the types needed for GBWT indexing.
     *
     * Takes a graph, a vector of contigs in the graph to process, in order,
     * and the corresponding VCF file, already open. Sample parsing on the VCF
     * file should be turned off.
     *
     * Inserts the sample names from the VCF into sample_names.
     *
     * Calls the callback serially with the contig number, each contig's
     * gbwt::VariantPaths, for each gbwt::PhasingInformation batch of samples.
     * The gbwt::PhasingInformation is not const because the GBWT library needs
     * to modify it in order to generate haplotypes from it efficiently.
     *
     * If batch_file_prefix is set on the object, also dumps VCF parse
     * information.
     *
     * If needed, this function can delete the graph to save memory.
     *
     * Returns the number of haplotypes created (2 per sample) This number will
     * need to be adjusted if any samples' haplotypes are filtered out later.
     * This function ignores any sample filters and processes the entire
     * sample range.
     */
    size_t parse_vcf(PathHandleGraph* graph, const std::vector<path_handle_t>& contigs,
        vcflib::VariantCallFile& variant_file, std::vector<std::string>& sample_names,
        const function<void(size_t, const gbwt::VariantPaths&, gbwt::PhasingInformation&)>& handle_contig_haplotype_batch,
        bool delete_graph) const;
    
    /**
     * Build a GBWT from the haplotypes in the given VCF file.
     *
     * Respects excluded_samples and does not produce threads for them.
     *
     * If needed, this function can delete the graph to save memory.
     *
     * TODO: We copy the file name, as vcflib requires a non-const name.
     */
    std::unique_ptr<gbwt::DynamicGBWT> build_gbwt(PathHandleGraph* graph, std::string vcf_filename,
        bool delete_graph) const;

    /**
     * Build a GBWT from the haplotypes in the given VCF file, but only
     * for the specified paths.
     *
     * Respects excluded_samples and does not produce threads for them.
     *
     * If needed, this function can delete the graph to save memory.
     *
     * TODO: We copy the file name, as vcflib requires a non-const name.
     */
    std::unique_ptr<gbwt::DynamicGBWT> build_gbwt(PathHandleGraph* graph, std::string vcf_filename,
        const std::vector<path_handle_t>& path_handles, bool delete_graph) const;

    /**
     * Build a GBWT from the embedded non-alt paths in the graph. Use
     * paths_as_samples to choose whether we treat the paths as contigs or
     * samples.
     */
    std::unique_ptr<gbwt::DynamicGBWT> build_gbwt(const PathHandleGraph* graph) const;

    /**
     * Build a GBWT from the alignments. Each distinct alignment name becomes
     * a sample in the GBWT metadata. If there are multiple alignments with
     * the same name, the corresponding GBWT path names will have the same
     * sample identifier but different values in the count field.
     *
     * aln_format can be "GAM" or "GAF"
     */
    std::unique_ptr<gbwt::DynamicGBWT> build_gbwt(const PathHandleGraph* graph,
        const std::vector<std::string>& aln_filenames, const std::string& aln_format) const;
};

}

#endif
