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

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/// Special marker value for star alleles in genotype vectors.
/// A star allele (*) represents a haplotype that spans a nested site in the
/// parent but doesn't have a defined traversal at the child level.
constexpr int STAR_ALLELE_MARKER = -2;

/// Special marker value for missing alleles in genotype vectors.
/// Used when a parent allele doesn't traverse a child snarl and star_allele
/// mode is disabled. Outputs as '.' in VCF to maintain consistent ploidy.
constexpr int MISSING_ALLELE_MARKER = -1;

/**
 * Helper class that vcf writers can inherit from to for some common code to output sorted VCF
 */
class VCFOutputCaller {
public:
    VCFOutputCaller(const string& sample_name);

    virtual ~VCFOutputCaller();

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const;

    /// Add a variant to our buffer
    /// Returns false if the variant line length exceeds VCFOutputCaller::max_vcf_line_length
    bool add_variant(vcflib::Variant& var) const;

    /// Sort then write variants in the buffer
    /// snarl_manager needed if include_nested is true
    void write_variants(ostream& out_stream, const SnarlManager* snarl_manager = nullptr);

    /// Run vcffixup from vcflib
    void vcf_fixup(vcflib::Variant& var) const;

    /// Add a translation map
    void set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation);

    /// Assume writing nested snarls is enabled
    void set_nested(bool nested);

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

protected:

    /// add a traversal to the VCF info field in the format of a GFA W-line or GAF path
    void add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele,
                                 const Traversal& trav, bool reversed, bool one_based) const;
    /// legacy version of above
    void add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav, bool reversed, bool one_based) const;
    
    
    /// convert a traversal into an allele string
    string trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const;

    /// Convert a SnarlTraversal to the handle vector the clustering code works on.  Returns false
    /// (leaving out_trav unspecified) if the traversal cannot be represented: fewer than two visits
    /// (the "*" placeholder pushed for a star allele), or a visit carrying a child Snarl rather
    /// than a node, which NestedFlowCaller produces via SnarlGraph::embed_snarl.  (LegacyCaller
    /// expands its children into node visits in top_down_genotype, so it never reaches here.)
    static bool snarl_traversal_to_handles(const HandleGraph& graph, const SnarlTraversal& trav,
                                           Traversal& out_trav);

    /// The CORE LENGTH of a variant: the length of the longest allele after stripping the prefix
    /// and the suffix that every non-"*" allele shares.  This is the single definition of "how big
    /// is this variant" behind --cluster-min-len in BOTH vg call and vg deconstruct.  It is
    /// invariant to how much shared flanking context a caller keeps in its allele strings, which is
    /// the point: vg call flattens down to an anchor base while vg deconstruct emits the whole
    /// snarl interior, so a raw string length answers differently for the same variant.
    /// Consequences, all intended:
    ///   - the anchor base flatten_common_allele_ends must leave on every indel is a shared prefix,
    ///     so it is stripped: a 49bp indel measures 49, not 50.
    ///   - REF participates, so a pure deletion measures the deleted length.  A maximum over ALTs
    ///     alone measures 1 for a deletion of any size.
    ///   - "*" is a marker, not sequence, so it is excluded from both the affixes and the maximum.
    ///     That also neutralizes flatten_common_allele_ends being a no-op whenever a "*" is
    ///     present -- without -a because min_allele_len becomes 1 and max_flatten_len decrements to
    ///     0, and with -a because "*" matches no base at the first offset compared.  Either way the
    ///     un-flattened boundary sequence is common to every real allele, so it is stripped here.
    /// Note this measures the SPAN of the variant, not the size of any one event inside it: a
    /// haplotype differing from the reference at two bases 59bp apart has a core length of 60.
    static int64_t allele_core_length(const vector<string>& alleles);

    /// Merge near-identical called ALT alleles in an already-populated variant.  Must run AFTER
    /// SnarlCaller::update_vcf_info and after flatten_common_allele_ends, so that the genotyper and
    /// the allele-flattening both see the full pre-merge allele set: merging earlier drops the
    /// absorbed allele's reads from AD/DP and from the Poisson caller's total_other_support term.
    /// Rewrites the allele-indexed fields (alleles/alt, AT, AD, GL, GT, MAD) and records what was
    /// merged in the MAT info field.  Returns true if anything merged.
    bool merge_similar_alleles(const PathPositionHandleGraph& graph,
                               const vector<SnarlTraversal>& site_traversals,
                               vector<int>& site_genotype,
                               const string& sample_name,
                               vcflib::Variant& out_variant) const;

    /// print a vcf variant
    /// return value is taken from add_variant (see above)
    bool emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                      function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string = nullptr);

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

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    // The nesting INFO headers (LV/CH/PS/RC/RS/RD), for both vg call and vg deconstruct.
    //
    // One definition on purpose.  These used to be written out verbatim in two places, and
    // drifted: 54bfd0f2d corrected the CH description in graph_caller.cpp while deconstructor.cpp
    // -- the copy deconstruct actually emits -- kept the text that commit's own message called
    // false, so the released VCFs carried the wrong one.
    static string nesting_info_headers();

    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// The same as above, but print the snarl as if its orientation has been flipped
    string print_flipped_snarl(const Snarl& snarl, bool in_brackets = false) const;

    /// do the opposite of above
    /// So a string that looks like AACT(>12<17)TTT would invoke the callback three times with
    /// ("AACT", Snarl), ("", Snarl(12,-17)), ("TTT", Snarl(12,-17))
    /// The parameters are to be treated as unions:  A sequence fragment if non-empty, otherwise a snarl
    void scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const;

    // update the PS and LV tags in the output buffer (called in write_variants if include_nested is true)
    void update_nesting_info_tags(const SnarlManager* snarl_manager);
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;

    /// output buffers (1/thread) (for sorting)
    /// variants stored as strings (and position key pairs) because vcflib::Variant in-memory struct so huge
    mutable vector<vector<pair<pair<string, size_t>, string>>> output_variants;

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

    // optional node translation to apply to snarl names in variant IDs
    const unordered_map<nid_t, pair<string, size_t>>* translation;

    // need to write LV/PS info tags
    bool include_nested;

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
