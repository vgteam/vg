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
#include "ploidy_regions.hpp"
#include "vcf_genotype_likelihoods.hpp"
#include "vcf_record.hpp"
#include "site_tree.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/**
 * Writer of a sorted VCF, which variant callers inherit from or hold.
 *
 * Callers add records from any thread, and the writer buffers them compressed. `write_variants`
 * adds the nesting INFO tags, sorts the records and writes them. A caller that holds a writer,
 * such as MultiPassCaller, passes its own callbacks to `emit_variant` and adds its header lines
 * with `set_header_steps`.
 */
class VCFOutputCaller {
public:
    /**
     * The sort key of a buffered VCF record: contig, POS, ID, then block number.
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
    };

    /// Order buffered record keys by contig, POS, ID, then block number.
    static bool buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b);

    explicit VCFOutputCaller(const string& sample_name);

    virtual ~VCFOutputCaller();

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const;

    /**
     * Add a record to the output buffer, as difference block `block` of its site.
     *
     * @return false, without adding the record, if its line is longer than
     * `max_vcf_line_length`.
     */
    bool add_variant(vcflib::Variant& var, size_t block = 0) const;

    /// Set the ploidy of the regions of the reference that `regions` covers.
    void set_ploidy_regions(PloidyRegions regions) { ploidy_regions = std::move(regions); }

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
    /// two effectively-identical alleles collapses to 1/1 with a single ALT.  Uses the same
    /// similarity metric and the same core-length gate as "vg deconstruct -L/--cluster-min-len" (a
    /// length-weighted Jaccard, except that a pure deletion is scored against the site -- see
    /// weighted_traversal_similarity).  The gate is applied to the alleles each tool emits, and
    /// those sets differ, so the two can disagree at a given site:
    /// similarity is >= threshold to merge, and min_len > 0 restricts merging to sites whose
    /// core length reaches min_len bp (see allele_core_length).
    /// A threshold of 1.0 (the default) disables merging entirely.
    void set_allele_merge(double threshold, int64_t min_len);

    /// The set of reference contigs that actually have a record.  Reads the sort keys of the
    /// output buffer, so it costs nothing (no decompression) and does not need the snarl tree.
    /// Only meaningful once calling is finished and before write_variants() drains the buffer.
    unordered_set<string> get_output_contigs() const;

    /// Remove ##contig lines whose ID is not in keep, leaving every other line alone.
    /// A reference contig that produced no record is not worth declaring: with a gref cover
    /// most contigs are fragments, and on a human chromosome a third of them carry nothing.
    string prune_header_contigs(const string& header, const unordered_set<string>& keep) const;

    /// Callbacks that a caller adds to the writing of a site's record. Each is empty unless the
    /// caller needs it.
    struct SiteRecordSteps {
        /// The ID of the site's records, or null to print the ID from the snarl's bounds.
        const string* id = nullptr;
        /// Return whether called traversal `trav` is written as the reference allele, because it
        /// takes the reference traversal's route through the site.
        function<bool(const Snarl& site, const vector<SnarlTraversal>& travs, int trav,
                      int ref_trav_idx)> same_as_reference;
        /// Count the site, before its record is built.
        function<void(const PathPositionHandleGraph& graph, const Snarl& site,
                      const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                      int ref_trav_idx)> count_site;
        /// Phase the record's genotype, as `SiteHooks::phase` does.
        function<int64_t(const Snarl& site, const vector<int>& site_genotype,
                         const map<int, int>& trav_to_allele, string& gt)> phase;
        /// Return the order in which the snarl caller wrote the GL of a call.
        function<GLLayout(const SnarlCaller::CallInfo* call_info)> gl_layout;
        /// Write the site as several records. Return the number of lines written, or -1 to have
        /// the site record written instead.
        function<int(const PathPositionHandleGraph& graph, const Snarl& site,
                     const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                     int ref_trav_idx, const SiteRecord& record, GLLayout gl_layout,
                     bool genotype_snarls)> write_blocks;
        /// Change the finished site record just before it is added to the output.
        function<void(vcflib::Variant& record)> finish_record;
        /// Add the caller's INFO and FORMAT fields, as `SiteHooks::fill_info` does, in place of
        /// the snarl caller's `update_vcf_info`.
        function<void(const vector<int>& site_trav, const vector<int>& site_genotype,
                      vcflib::Variant& variant)> fill_info;
        /// Take note of a site whose records are done. `trav_to_allele` gives the VCF allele of
        /// each traversal in the genotype, out of `traversal_count` traversals. `has_line` says
        /// whether the site got a line.
        function<void(const Snarl& site, const map<int, int>& trav_to_allele,
                      size_t traversal_count, bool has_line)> site_filed;
    };

    /// Callbacks that add a caller's lines to the VCF header. Each is empty unless the caller
    /// needs it.
    struct HeaderSteps {
        /// Return FORMAT lines, which go after the nesting INFO lines.
        function<string()> format_header;
        /// Return INFO lines, which go after the AT line.
        function<string()> info_header;
    };

    /// Add the lines from `steps` to the header that `vcf_header` writes.
    void set_header_steps(HeaderSteps steps) { header_steps = std::move(steps); }

    /**
     * Add the record for a called site to the output buffer.
     *
     * We build the record with build_site_record, from the called traversals and the snarl
     * caller's INFO and FORMAT fields, and run the callbacks in `steps` along the way.
     *
     * @param trav_to_string Spells an allele. When null, an allele is spelled as the sequence of
     * its traversal.
     * @return false only if add_variant refused a line the site wanted.
     */
    bool emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int64_t ref_offset, bool genotype_snarls, int ploidy,
                      const SiteRecordSteps& steps,
                      function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string = nullptr);

    /// Add the record for a called site to the output buffer, with no added callbacks.
    bool emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int64_t ref_offset, bool genotype_snarls, int ploidy,
                      const function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)>& trav_to_string = nullptr);

    /**
     * Return the VCF header of a caller that genotypes with `snarl_caller`.
     *
     * The header has the base lines, GT, the lines of `snarl_caller`, FILTER, SAMPLE and the
     * column line. We also open the output VCF with it.
     */
    string snarl_caller_vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides,
                                   const SnarlCaller& snarl_caller) const;

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// The same as above, but print the snarl as if its orientation has been flipped
    string print_flipped_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// Print a snarl from its two boundary visits, as the other print_snarl overloads do.
    string print_snarl(nid_t start_id, bool start_backward, nid_t end_id, bool end_backward,
                       bool in_brackets) const;

    /// convert a traversal into an allele string
    string trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const;

    /// Get the node translation given to `set_translation`, or null.
    const unordered_map<nid_t, pair<string, size_t>>* get_translation() const { return translation; }

protected:

    /// add a traversal to the VCF info field in the format of a GFA W-line or GAF path
    void add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele,
                                 const Traversal& trav, bool reversed, bool one_based) const;
    /// legacy version of above
    void add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav, bool reversed, bool one_based) const;
    
    
    /// Return the core length of a variant, as `vg::allele_core_length` defines it.
    static int64_t allele_core_length(const vector<string>& alleles) {
        return vg::allele_core_length(alleles);
    }

    /// The callbacks that add a caller's lines to the header.
    HeaderSteps header_steps;

    /// Return the options build_site_record takes from this caller.
    RecordOptions record_options() const;

    /// get the interval of a snarl from our reference path using the PathPositionHandleGraph interface
    /// the bool is true if the snarl's backward on the path
    /// first returned value -1 if no traversal found 
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(const PathPositionHandleGraph& graph, const Snarl& snarl,
                                                                                 const string& ref_path_name) const;

    /// used for making gaf traversal names
    pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                           int64_t ref_path_offset) const;

    /// clean up the alleles to not share common prefixes / suffixes
    /// if len_override given, just do that many bases without thinking
    void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const;

    // The nesting INFO headers (LV/CH/PS/RC/RS/RD), for both vg call and vg deconstruct.
    //
    // One definition on purpose.  These used to be written out verbatim in two places, and
    // drifted: 54bfd0f2d corrected the CH description in graph_caller.cpp while deconstructor.cpp
    // -- the copy deconstruct actually emits -- kept the text that commit's own message called
    // false, so the released VCFs carried the wrong one.
    static string nesting_info_headers();

    /// do the opposite of print_snarl
    /// So a string that looks like AACT(>12<17)TTT would invoke the callback three times with
    /// ("AACT", Snarl), ("", Snarl(12,-17)), ("TTT", Snarl(12,-17))
    /// The parameters are to be treated as unions:  A sequence fragment if non-empty, otherwise a snarl
    void scan_snarl(const string& allele_string, const function<void(const string&, Snarl&)>& callback) const;

    /// Add the nesting INFO tags (LV, CH, PS, RC, RS, RD) to every buffered record, from the tree
    /// of sites.
    void update_nesting_info_tags(const SiteTree& sites);
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;

    /// The buffered records, one buffer per thread, each a sort key and the zstd-compressed line.
    /// A vcflib::Variant is too big to keep in memory for every record.
    mutable vector<vector<pair<BufferedRecordKey, string>>> output_variants;

    /// Reference interval of a site that was visited but not emitted, because every traversal
    /// through it was the reference (or absent) and so it had no variant to report.  Such a site
    /// is invisible to the RC/RS/RD walk, which only sees sites that reached the VCF, and a record
    /// nested under one would otherwise have no reference coordinate to point at.  Common in gref
    /// graphs, where the parent of an island of non-reference sequence is often a large snarl that
    /// only the reference and its own gref copy span.
    ///
    /// Keyed by snarl name as print_snarl() spells it, which is how record IDs and chrom_of_name
    /// are keyed too.  One buffer per thread, like output_variants, merged in
    /// update_nesting_info_tags().
    struct SuppressedRef {
        string chrom;
        size_t pos;
        size_t ref_len;
    };
    mutable vector<unordered_map<string, SuppressedRef>> suppressed_ref_info;

    /// print up to this many uncalled alleles when doing ref-genotpes in -a mode
    size_t max_uncalled_alleles = 5;

    /// The ploidy of regions of the reference, which overrides the ploidy of their paths. Empty
    /// unless the run gave a BED file of them.
    PloidyRegions ploidy_regions;

    // optional node translation to apply to snarl names in variant IDs
    const unordered_map<nid_t, pair<string, size_t>>* translation;

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

    // prevent giant variants
    static const int64_t max_vcf_line_length = 2000000000;
};

}

#endif
