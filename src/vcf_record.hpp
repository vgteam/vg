#ifndef VG_VCF_RECORD_HPP_INCLUDED
#define VG_VCF_RECORD_HPP_INCLUDED

/** \file
 * Building VCF records for sites in a graph.
 *
 * A site is a region of the graph between two boundary nodes, and these functions know it only by
 * those two nodes, given as handles oriented into the site. An allele is one way through the site,
 * called a traversal. The functions read a site's traversals through callbacks, so they serve any
 * caller, whatever structure it keeps its sites and traversals in.
 */

#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Variant.h"
#include "handle.hpp"
#include "traversal_finder.hpp"
#include "vcf_genotype_likelihoods.hpp"

namespace vg {

using namespace std;

/// Special marker value for star alleles in genotype vectors.
/// A star allele (*) represents a haplotype that spans a nested site in the
/// parent but doesn't have a defined traversal at the child level.
constexpr int STAR_ALLELE_MARKER = -2;

/// Special marker value for missing alleles in genotype vectors.
/// Used when a parent allele doesn't traverse a child snarl and star_allele
/// mode is disabled. Outputs as '.' in VCF to maintain consistent ploidy.
constexpr int MISSING_ALLELE_MARKER = -1;

/// For each node ID of a graph, the name of the node it was made from, and the node's offset in
/// it. Records then name nodes by these names instead of by ID.
using NodeTranslation = unordered_map<nid_t, pair<string, size_t>>;

/// A visit to a node, as the node's ID and whether the visit reads the node backward. A visit to a
/// child site rather than to a node has node ID 0.
using NodeVisit = pair<nid_t, bool>;

/// The 1-based position on the base path of the base `along_path` bases into `ref_path_name`, a
/// path that may name a subrange of its base path.
int64_t base_path_position(const string& ref_path_name, int64_t along_path);

/// The ID of a site's record: each boundary node as ">" or "<" (forward or backward) and its ID, or
/// its name under `translation` if that is not null. With `in_brackets`, the ID is put in
/// parentheses, as when it stands for a child site inside an allele.
string site_name(nid_t start_id, bool start_backward, nid_t end_id, bool end_backward,
                 const NodeTranslation* translation, bool in_brackets = false);

/// The same, for the site with bounds `start` and `end`.
string site_name(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                 const NodeTranslation* translation, bool in_brackets = false);

/// The key that identifies a site's record from one pass of a caller to the next: a hash of the
/// record's ID.
size_t record_key_of(const string& site_id);

/// The same, for the site with bounds `start` and `end`.
size_t record_key_of(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                     const NodeTranslation* translation);

/// Where the site with bounds `start` and `end` lies on reference path `ref_path_name`: the
/// positions of its two boundary nodes on the path, in path order, whether the site runs backward
/// along the path, and the steps at those two positions. The positions are -1 if the path does not
/// pass through the site from one boundary node to the other.
tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(
    const PathPositionHandleGraph& graph, const handle_t& start, const handle_t& end,
    const string& ref_path_name);

/// The base path of `ref_path_name`, and the 1-based position on it where the site with bounds
/// `start` and `end` begins, `ref_path_offset` bases further on.
pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const handle_t& start,
                                       const handle_t& end, const string& ref_path_name,
                                       int64_t ref_path_offset);

/// The visits of a traversal, with node ID 0 for a visit to a child snarl.
vector<NodeVisit> visits_of(const SnarlTraversal& trav);

/// Write the nodes an allele visits into INFO/AT for allele number `allele`, as a walk in the form
/// of a GFA W-line or a GAF path, read in reverse if `reversed` is set. Nodes are named by
/// `translation` if it is not null, and consecutive visits to pieces of one translated node are
/// written once. An allele that visits nothing is written as ".".
void add_allele_path_to_info(vcflib::Variant& v, int allele, const vector<NodeVisit>& visits,
                             bool reversed, const NodeTranslation* translation);

/// The visits of an allele as handles, for clustering. False, leaving `walk` unspecified, if the
/// allele cannot be given as handles: when it has fewer than two visits, as the placeholder for a
/// "*" allele has, or a visit to a child site, which has no single handle.
bool visits_to_walk(const HandleGraph& graph, const vector<NodeVisit>& visits, Traversal& walk);

/// The core length of a variant: the length of the longest allele after stripping the prefix and
/// the suffix that every non-"*" allele shares. vg call and vg deconstruct both use it to decide
/// whether a variant is long enough to have its alleles merged. It does not depend on how much
/// shared flanking sequence a caller keeps in its alleles, which differs between the two. It follows
/// that:
///   - the anchor base that flatten_common_allele_ends leaves on every indel is a shared prefix,
///     so it is stripped: a 49bp indel measures 49, not 50.
///   - REF takes part, so a pure deletion measures the deleted length. A maximum over the ALTs
///     alone would measure 1 for a deletion of any size.
///   - "*" is a marker, not sequence, so it is left out of both the affixes and the maximum.
///     flatten_common_allele_ends trims nothing from a record with a "*", and the boundary
///     sequence it leaves there is common to every real allele, so it is stripped here.
/// It measures the span of the variant, not the size of any one event inside it: a haplotype
/// differing from the reference at two bases 59bp apart has a core length of 60.
int64_t allele_core_length(const vector<string>& alleles);

/// Merge near-identical called ALT alleles in an already populated variant. Must run after the
/// caller's INFO and FORMAT fields are written and after flatten_common_allele_ends, so that both
/// see every allele: merging earlier would drop the absorbed allele's reads from AD and DP.
/// Rewrites the allele-indexed fields (alleles/alt, AT, AD, GL, GT, MAD) and records the merge in
/// INFO/MAT. Returns true if anything merged.
///
/// `allele_walk` gives allele `a` of the `allele_count` alleles as handles, or returns false for an
/// allele that is to be left alone. Alleles are merged when their walks are at least `threshold`
/// similar, and only in a variant whose core length is at least `min_len`. `gl_layout` is the order
/// in which the caller wrote GL, which cannot be recovered from the record.
bool merge_similar_alleles(const PathPositionHandleGraph& graph, size_t allele_count,
                           const function<bool(size_t allele, Traversal& walk)>& allele_walk,
                           vector<int>& site_genotype, const string& sample_name,
                           vcflib::Variant& out_variant, GLLayout gl_layout, double threshold,
                           int64_t min_len);

/// Trim the sequence every allele shares from the start of the alleles, or from the end if
/// `backward` is set, leaving at least one base in each. If `len_override` is not 0, trim that many
/// bases instead, without comparing them.
void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override);

/// What stays the same for every record a caller writes.
struct RecordOptions {
    /// The sample the FORMAT fields are written for.
    string sample_name;
    /// Names for node IDs, or null to write the IDs.
    const NodeTranslation* translation;
    /// The most alleles that were not called that are added to a reference call written anyway.
    size_t max_uncalled_alleles;
    /// Similar alleles are merged when their walks are at least this similar; 1 merges none.
    double allele_merge_threshold;
    /// Alleles are merged only in a record whose core length is at least this.
    int64_t allele_merge_min_len;
};

/// The site a record is built for, and its call.
struct SiteToWrite {
    /// The site's bounds: its two boundary nodes, oriented into the site.
    handle_t start;
    handle_t end;
    /// The reference path the record is placed on, which may name a subrange of its base path.
    const string& ref_path_name;
    /// Added to positions on the reference path.
    int ref_offset;
    /// The called genotype, as numbers of the site's traversals or as the allele markers. Empty
    /// when the site has no call.
    const vector<int>& genotype;
    /// The traversal that is the reference allele.
    int ref_trav_idx;
    /// How many traversals the site has.
    size_t traversal_count;
    /// How many "." the genotype of a site with no call has.
    int ploidy;
    /// Write the record even for a reference call, with alleles that were not called added, and
    /// trim the alleles only as far as the site's boundary nodes.
    bool genotype_snarls;
};

/// A site's traversals, numbered from 0, which build_site_record reads through these functions.
struct SiteAlleles {
    /// The sequence written for traversal `trav` where it fills entry `genotype_index` of the
    /// genotype. A caller that writes a child site inside an allele spells it from the child's own
    /// call, which can differ between the entries of the genotype.
    function<string(int trav, int genotype_index)> spell;
    /// The visits of traversal `trav`.
    function<vector<NodeVisit>(int trav)> visits;
};

/// Steps a caller adds to building a record. Any of the functions may be left empty.
struct SiteHooks {
    /// Whether traversal `trav` takes the reference traversal's route through the site, differing
    /// only inside child sites, so that it is written as the reference allele. When empty, only the
    /// reference traversal is.
    function<bool(int trav)> same_as_reference;
    /// Phases the genotype. It is given the genotype in VCF allele numbers and the VCF allele of
    /// each traversal in the genotype. It may replace the GT string, and returns the phase set to
    /// write as PS, or -1 for none.
    function<int64_t(const vector<int>& site_genotype, const map<int, int>& trav_to_allele,
                     string& gt)> phase;
    /// Adds the caller's INFO and FORMAT fields, after GT and before PS. `site_trav` gives the
    /// traversal of each VCF allele, or -1 for "*".
    function<void(const vector<int>& site_trav, const vector<int>& site_genotype,
                  vcflib::Variant& variant)> fill_info;
    /// The order in which `fill_info` writes GL.
    GLLayout gl_layout;
};

/// A record from build_site_record, with what its caller needs to file it.
struct SiteRecord {
    vcflib::Variant variant;
    /// The genotype in VCF allele numbers, with MISSING_ALLELE_MARKER for ".".
    vector<int> genotype;
    /// The VCF allele of each traversal in the genotype.
    map<int, int> trav_to_allele;
    /// The record's position before its alleles were trimmed: that of the first base of the
    /// reference allele.
    int64_t unflattened_position;
    /// Whether similar alleles were merged.
    bool alleles_merged;
};

/// A nested ploidy-1 genotype: one allele on a named strand, with "." on the other, since the other
/// strand carries nothing here, its parent allele having deleted the chain. Shared by a site record
/// and its block records.
string nested_strand_genotype(int allele, int strand);

/// Build the record for a site. Traversals with the same sequence become one allele. Where the
/// genotype has an allele marker, the record has "*" or ".", and a genotype with "." and no ALT
/// gets "*" as its ALT, so that the record is valid VCF.
SiteRecord build_site_record(const PathPositionHandleGraph& graph, const SiteToWrite& site,
                             const SiteAlleles& alleles, const SiteHooks& hooks,
                             const RecordOptions& options);

}

#endif
