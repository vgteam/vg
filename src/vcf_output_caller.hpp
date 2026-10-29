#ifndef VG_VCF_OUTPUT_CALLER_HPP_INCLUDED
#define VG_VCF_OUTPUT_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "../artifacts/vcf_record.hpp"
#include "../artifacts/site_tree.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/// Return the stretch of the reference that the REF of `record` covers, 1-based and end-inclusive,
/// as VCF POS is.
Region record_ref_region(const vcflib::Variant& record);

/**
 * Writer of a sorted VCF, which variant callers inherit from or hold.
 *
 * Callers add records from any thread, and the writer buffers them compressed. `write_variants`
 * adds the nesting INFO tags, sorts the records and writes them. A caller builds the record of
 * each site with a `SiteEmitter`, and adds its header lines with `add_format_header` and
 * `add_info_header`.
 */
class VCFOutputCaller {
public:
    /**
     * The sort key of a buffered VCF record.
     *
     * Several records can share a contig and POS, such as a nested site and its parent. The sort
     * is not stable, so the output order is the same on every run only if the records at one POS
     * differ in ID or block number.
     */
    struct BufferedRecordKey {
        string contig;
        size_t position = 0;
        /// The ID column.
        string id;
        /// The number of the difference block of its snarl that the record holds, or 0 for a
        /// record of a whole snarl. Two blocks of one snarl share an ID and can share a POS.
        size_t block = 0;

        /// Return true if this key sorts before `other`, comparing contig, then POS, then ID, then
        /// block number.
        bool operator<(const BufferedRecordKey& other) const;
    };

    explicit VCFOutputCaller(const string& sample_name);

    virtual ~VCFOutputCaller();

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const;

    /**
     * Add a record to the output.
     *
     * A hidden record is not written. The writer keeps only its ID and the stretch of the
     * reference its REF covers, so that a record nested inside it can still be placed on the
     * reference when the nesting tags are computed.
     *
     * @param block The number of the difference block of its site that the record holds, or 0
     * for a record of a whole site. Two blocks of one site share an ID and can share a POS.
     * @return false, without adding the record, if its line is longer than
     * `max_vcf_line_length`.
     */
    bool add_variant(vcflib::Variant& var, size_t block = 0, bool hidden = false) const;

    /**
     * Write the buffered records, sorted, taking the nesting INFO tags from `snarl_manager`.
     *
     * This empties the buffer, so it can be called once. `snarl_manager` may be null if nested
     * output is off.
     */
    void write_variants(ostream& out_stream, const SnarlManager* snarl_manager = nullptr);

    /**
     * Write the buffered records, sorted, taking the nesting INFO tags from the tree of sites.
     *
     * This empties the buffer, so it can be called once. `sites` may be null if nested output is
     * off.
     */
    void write_variants(ostream& out_stream, const SiteTree* sites);

    /// Run vcffixup from vcflib
    void vcf_fixup(vcflib::Variant& var) const;

    /// Add a translation map
    void set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation);

    /// Assume writing nested snarls is enabled
    void set_nested(bool nested);

    /// Set the gRef nesting level of each gRef contig, by contig name.
    ///
    /// A record on one of these contigs gets at least its level as INFO/CH.
    void set_gref_levels(map<string, int> levels);

    /// Enable post-genotyping merging of near-identical called ALT alleles, so that a 1/2 call of
    /// two effectively-identical alleles collapses to 1/1 with a single ALT.
    ///
    /// Two alleles merge when their similarity, a modified length-weighted Jaccard, is at least
    /// `threshold`. A threshold of 1.0 turns merging off.
    void set_allele_merge(double threshold, int64_t min_len);

    /// The set of reference contigs that actually have a record.  Reads the sort keys of the
    /// output buffer, so it costs nothing (no decompression) and does not need the snarl tree.
    /// Only meaningful once calling is finished and before write_variants() drains the buffer.
    unordered_set<string> get_output_contigs() const;

    /// Remove ##contig lines whose ID is not in keep, leaving every other line alone.
    /// A reference contig that produced no record is not worth declaring: with a gref cover
    /// most contigs are fragments, and on a human chromosome a third of them carry nothing.
    string prune_header_contigs(const string& header, const unordered_set<string>& keep) const;

    /// Add a FORMAT line, without its newline, to the header that `vcf_header` writes, after the
    /// nesting INFO lines.
    void add_format_header(string line);

    /// Add an INFO line, without its newline, to the header that `vcf_header` writes, after the
    /// AT line.
    void add_info_header(string line);

    /**
     * Add the record for a called site to the output buffer.
     *
     * The record has the called traversals as its alleles, each spelled as the sequence of its
     * traversal, and the snarl caller's INFO and FORMAT fields.
     *
     * @return false only if add_variant refused a line the site wanted.
     */
    bool emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int64_t ref_offset, bool genotype_snarls, int ploidy);

    /**
     * Return the VCF header of a caller that genotypes with `snarl_caller`.
     *
     * The header has the base lines, GT, the lines of `snarl_caller`, FILTER, SAMPLE and the
     * column line. We also open the output VCF with it.
     */
    string snarl_caller_vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides,
                                   const SnarlCaller& snarl_caller) const;

    /// Return the VCF header of a caller whose genotyper adds its own lines to the header with
    /// `add_genotyper_lines`.
    string genotyper_vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                const vector<size_t>& contig_length_overrides,
                                const function<void(string& header)>& add_genotyper_lines) const;

    /// Get the node translation given to `set_translation`, or null.
    const NodeTranslation* translation() const;

    /// Get the name of the sample the records are for.
    const string& sample_name() const;

    /// A record longer than this is refused.
    static const int64_t max_vcf_line_length = 2000000000;

protected:

    /// The FORMAT lines a caller added to the header.
    vector<string> format_header_lines;
    /// The INFO lines a caller added to the header.
    vector<string> info_header_lines;

    // The nesting INFO headers (LV/CH/PS/RC/RS/RD), for both vg call and vg deconstruct.
    //
    // One definition on purpose.  These used to be written out verbatim in two places, and
    // drifted: 54bfd0f2d corrected the CH description in graph_caller.cpp while deconstructor.cpp
    // -- the copy deconstruct actually emits -- kept the text that commit's own message called
    // false, so the released VCFs carried the wrong one.
    static string nesting_info_headers();

    /// Add the nesting INFO tags (LV, CH, PS, RC, RS, RD) to every buffered record, from the tree
    /// of sites.
    void update_nesting_info_tags(const SiteTree& sites);
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample;

    /// The buffered records, one buffer per thread, each a sort key and the zstd-compressed line.
    /// A vcflib::Variant is too big to keep in memory for every record.
    mutable vector<vector<pair<BufferedRecordKey, string>>> output_variants;


    // optional node translation to apply to snarl names in variant IDs
    const NodeTranslation* node_translation;

    // need to write LV/PS info tags
    bool include_nested;
    /// The gRef nesting level of each gRef contig, by contig name. Empty without a gRef cover.
    map<string, int> gref_levels;

    // post-genotyping ALT merging (vg call -L / --cluster-min-len).  Deliberately NOT named
    // cluster_threshold / cluster_min_allele_len: Deconstructor derives from this class and already
    // declares both for its own pre-allele-string clustering, and -Wshadow is silent when a derived
    // member shadows a base one.
    double allele_merge_threshold = 1.0;
    int64_t allele_merge_min_len = 0;


private:
    /// The reference interval of each hidden record, by its ID, one table per thread like
    /// `output_variants`.
    ///
    /// A hidden record is one of a site that had no variant to report. A record nested under such
    /// a site may have no other reference coordinate to point its RC, RS and RD at. This is
    /// common in gref graphs, where the parent of an island of non-reference sequence is often a
    /// large snarl that only the reference and its own gref copy span.
    mutable vector<unordered_map<string, Region>> suppressed_ref_info;
};

}

#endif
