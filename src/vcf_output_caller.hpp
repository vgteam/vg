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
#include "../artifacts/ploidy_regions.hpp"
#include "../artifacts/vcf_record.hpp"
#include "../artifacts/site_tree.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/**
 * Writer of a sorted VCF, which variant callers inherit from or hold.
 *
 * Callers add records from any thread, and the writer buffers them compressed. `write_variants`
 * adds the nesting INFO tags, sorts the records and writes them. A caller that holds a writer,
 * such as MultiPassCaller, writes each site with its own subclass of `SiteEmitter` and adds its
 * header lines with `add_format_header` and `add_info_header`.
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
     * Add a record to the output buffer, as difference block `block` of its site.
     *
     * @return false, without adding the record, if its line is longer than
     * `max_vcf_line_length`.
     */
    bool add_variant(vcflib::Variant& var, size_t block = 0) const;

    /// Set the ploidy of the regions of the reference that `regions` covers.
    void set_ploidy_regions(PloidyRegions regions);

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

    /**
     * A strategy for emitting the records of individual sites. Subclass it to customize how a
     * site is emitted.
     */
    class SiteEmitter {
    public:
        /// Make an emitter that adds records to `output`, which must outlive it.
        ///
        /// @param genotype_snarls Write the record even for a reference call, with alleles that
        /// were not called added, and trim the alleles only as far as the site's boundary nodes.
        SiteEmitter(VCFOutputCaller& output, bool genotype_snarls);
        virtual ~SiteEmitter() = default;

        /**
         * Build the record of `site` and add it to the output.
         *
         * A record with no ALT gets no line unless `genotype_snarls` is set. For a site with
         * no line, keep its reference interval for its children's nesting tags.
         *
         * @return false only if add_variant refused a line the site wanted.
         */
        bool emit(const PathPositionHandleGraph& graph, const CalledSite& site);

    protected:
        // The alleles, which every subclass gives.

        /// Return the sequence written for traversal `trav` where it fills entry `genotype_index`
        /// of the genotype.
        virtual string spell(int trav, int genotype_index) = 0;

        /// Return the visits of traversal `trav`.
        virtual vector<NodeVisit> visits(int trav) = 0;

        // Decisions about the record.

        /// Return whether traversal `trav`, which is not the reference traversal, is written as
        /// the reference allele, because it takes the reference allele's route through the site.
        /// By default no other traversal is.
        virtual bool same_as_reference(int trav, int ref_trav_idx);

        /// Phase the genotype of the record, and return the phase set to write as PS, or -1 for
        /// none.
        ///
        /// `site_genotype` is in VCF allele numbers, and `trav_to_allele` gives the VCF allele of
        /// each traversal in the genotype. It may replace the GT string `gt`. By default the
        /// genotype stays unphased, with no PS.
        virtual int64_t phase(const vector<int>& site_genotype, const map<int, int>& trav_to_allele,
                              string& gt);

        /// Return the order in which `fill_info` writes GL. By default it is i-major.
        virtual GLLayout gl_layout();

        /// Write the site as several records, and return the number of lines written, or -1 to
        /// have the site record written instead. By default the site record is written.
        virtual int write_blocks(const PathPositionHandleGraph& graph, const vector<int>& genotype,
                                 int ref_trav_idx, const SiteRecord& record, GLLayout gl_layout);

        // Hooks, which do nothing by default.

        /// Count the site, before its record is built.
        virtual void count_site(const vector<int>& genotype, int ref_trav_idx);

        /// Add the caller's INFO and FORMAT fields to `variant`.
        ///
        /// `site_trav` gives the traversal of each VCF allele, or -1 for "*". The record has its GT
        /// and not yet its PS.
        virtual void fill_info(const vector<int>& site_trav, const vector<int>& site_genotype,
                               vcflib::Variant& variant);

        /// Change the finished site record just before it is added to the output.
        virtual void finish_record(vcflib::Variant& record);

        /// Take note of a site whose records are done. `trav_to_allele` gives the VCF allele of
        /// each called allele in the genotype, out of `traversal_count` alleles. `has_line` says
        /// whether the site got a line. A site written as blocks gives an empty map.
        virtual void site_filed(const map<int, int>& trav_to_allele, size_t traversal_count,
                                bool has_line);

        /// The output the records go to.
        VCFOutputCaller& output;
        /// Whether to write the record even for a reference call, with alleles that were not
        /// called added, and to trim the alleles only as far as the site's boundary nodes.
        const bool genotype_snarls;

    private:
        /// Build the record for `site`, with GL in order `layout`. Traversals with the same
        /// sequence become one allele. Where the genotype has an allele marker, the record has
        /// "*" or ".", and a genotype with "." and no ALT gets "*" as its ALT, so that the record
        /// is valid VCF.
        SiteRecord build_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                                GLLayout layout);

        /// Number the reference allele and the alleles of the called genotype of `site`, whose
        /// REF is already in `record`. Traversals that spell the same sequence share one allele,
        /// and the reference allele is number 0.
        ///
        /// Fill the empty `allele_to_gt` with the VCF allele number of each allele sequence, and
        /// the empty `site_trav` with the traversal each VCF allele is spelled from, or -1 for
        /// "*". Fill in `record.genotype` and `record.trav_to_allele` to match.
        void number_called_alleles(const CalledSite& site, SiteRecord& record,
                                   map<string, int>& allele_to_gt, vector<int>& site_trav);

        /// Add up to the output's `max_uncalled_alleles` alleles of `site` that were not called,
        /// in order of their sequences, to `allele_to_gt` and `site_trav`.
        void add_uncalled_alleles(const CalledSite& site, map<string, int>& allele_to_gt,
                                  vector<int>& site_trav);

        /// Add a built site record to the output, or let `write_blocks` write it as blocks, and
        /// tell `site_filed` how the site was filed.
        ///
        /// @return false only if add_variant refused a line the site wanted.
        bool file_record(const PathPositionHandleGraph& graph, SiteRecord& record,
                         const CalledSite& site, GLLayout gl_layout);
    };

    /// Add a FORMAT line, without its newline, to the header that `vcf_header` writes, after the
    /// nesting INFO lines.
    void add_format_header(string line);

    /// Add an INFO line, without its newline, to the header that `vcf_header` writes, after the
    /// AT line.
    void add_info_header(string line);

    /**
     * Add the record for a called site to the output buffer.
     *
     * We build the record with a SiteEmitter, from the called traversals and the snarl caller's
     * INFO and FORMAT fields.
     *
     * @param trav_to_string Spells an allele. When null, an allele is spelled as the sequence of
     * its traversal.
     * @return false only if add_variant refused a line the site wanted.
     */
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

    /// Return the VCF header of a caller whose genotyper adds its lines to the header with
    /// `add_genotyper_lines`, as `snarl_caller_vcf_header` does with a snarl caller's.
    string genotyper_vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                const vector<size_t>& contig_length_overrides,
                                const function<void(string& header)>& add_genotyper_lines) const;

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// The same as above, but print the snarl as if its orientation has been flipped
    string print_flipped_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// Print a snarl from its two boundary visits, as >3435<12222, or as (>3435<12222) with
    /// `in_brackets`.
    string print_snarl(nid_t start_id, bool start_backward, nid_t end_id, bool end_backward,
                       bool in_brackets = false) const;

    /// convert a traversal into an allele string
    string trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const;

    /// Get the node translation given to `set_translation`, or null.
    const unordered_map<nid_t, pair<string, size_t>>* get_translation() const;

protected:

    /// add a traversal to the VCF info field in the format of a GFA W-line or GAF path
    void add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele,
                                 const Traversal& trav, bool reversed) const;
    
    
    /// The writer of a site called by a SnarlCaller, which `emit_variant` uses.
    class SnarlTraversalEmitter;

    /// The FORMAT lines a caller added to the header.
    vector<string> format_header_lines;
    /// The INFO lines a caller added to the header.
    vector<string> info_header_lines;

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

    /// The ploidy of regions of the reference, which overrides the ploidy of their paths.
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
