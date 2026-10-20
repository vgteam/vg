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

/// The genotype entry for a star allele, written as "*".
///
/// It stands for a strand that visits the parent site but has no traversal of this child site.
constexpr int STAR_ALLELE_MARKER = -2;

/// The genotype entry for a missing allele, written as ".".
///
/// It stands for a strand whose parent allele does not visit this child site, when star alleles
/// are off. Writing "." keeps the genotype at its ploidy.
constexpr int MISSING_ALLELE_MARKER = -1;

/// For each node ID of a graph, the name of the node it was made from, and the node's offset in
/// it. Records then name nodes by these names instead of by ID.
using NodeTranslation = unordered_map<nid_t, pair<string, size_t>>;

/// A visit to a node, as the node's ID and whether the visit reads the node backward. A visit to a
/// child site rather than to a node has node ID 0.
using NodeVisit = pair<nid_t, bool>;

/// Return the 1-based position on the base path of the base `along_path` bases into
/// `ref_path_name`, a path that may name a subrange of its base path.
int64_t base_path_position(const string& ref_path_name, int64_t along_path);

/// Return the ID of the record of a site.
///
/// The ID gives each boundary node as ">" or "<", for forward or reverse, and its ID, or its name
/// under `translation` if that is not null. With `in_brackets`, the ID is put in parentheses, as
/// when it stands for a child site inside an allele.
string site_name(nid_t start_id, bool start_backward, nid_t end_id, bool end_backward,
                 const NodeTranslation* translation, bool in_brackets = false);

/// Return the ID of the record of the site with bounds `start` and `end`.
string site_name(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                 const NodeTranslation* translation, bool in_brackets = false);

/// Return the key that identifies the record of a site from one pass of a caller to the next.
///
/// The key is a hash of the ID of the record.
size_t record_key_of(const string& site_id);

/// Return the key of the record of the site with bounds `start` and `end`.
size_t record_key_of(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                     const NodeTranslation* translation);

/// Find where the site with bounds `start` and `end` lies on reference path `ref_path_name`.
///
/// Return the positions of its two boundary nodes on the path, in path order, whether the site
/// runs backward along the path, and the steps at those two positions. The positions are -1 if
/// the path does not go through the site from one boundary node to the other.
tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(
    const PathPositionHandleGraph& graph, const handle_t& start, const handle_t& end,
    const string& ref_path_name);

/// Return the base path of `ref_path_name`, and the 1-based position on it where the site with
/// bounds `start` and `end` begins, moved `ref_path_offset` bases further on.
pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const handle_t& start,
                                       const handle_t& end, const string& ref_path_name,
                                       int64_t ref_path_offset);

/// Return the visits of a traversal, with node ID 0 for a visit to a child snarl.
vector<NodeVisit> visits_of(const SnarlTraversal& trav);

/// Return the visits of a walk.
vector<NodeVisit> visits_of(const HandleGraph& graph, const Traversal& walk);

/// Write the nodes an allele visits into INFO/AT for allele number `allele`, as a walk in the form
/// of a GFA W-line or a GAF path, read in reverse if `reversed` is set. Nodes are named by
/// `translation` if it is not null, and consecutive visits to pieces of one translated node are
/// written once. An allele that visits nothing is written as ".".
void add_allele_path_to_info(vcflib::Variant& v, int allele, const vector<NodeVisit>& visits,
                             bool reversed, const NodeTranslation* translation);

/// Put the visits of an allele into `walk` as handles, for clustering.
///
/// Return false, leaving `walk` unspecified, if the allele has fewer than two visits, as the
/// placeholder for a "*" allele has, or has a visit to a child site, which has no single handle.
bool visits_to_walk(const HandleGraph& graph, const vector<NodeVisit>& visits, Traversal& walk);

/**
 * Return the core length of a variant.
 *
 * That is the length of the longest allele after stripping the prefix and the suffix that every
 * allele other than "*" shares. vg call and vg deconstruct both use it to decide whether a variant
 * is long enough to have its alleles merged. It does not depend on how much shared flanking
 * sequence a caller keeps in its alleles, which differs between the two.
 *
 * The anchor base that flatten_common_allele_ends leaves on an indel is shared, so a 49bp indel
 * measures 49. REF counts, so a pure deletion measures the deleted length. "*" is not sequence,
 * so it counts toward neither the affixes nor the maximum.
 *
 * It measures the span of the variant: a haplotype that differs from the reference at two bases
 * 59bp apart has a core length of 60.
 */
int64_t allele_core_length(const vector<string>& alleles);

/**
 * Merge near-identical called ALT alleles in a finished variant. Return true if any merged.
 *
 * Run it after the caller's INFO and FORMAT fields are written, and after
 * flatten_common_allele_ends. Merging earlier would drop the reads of the absorbed allele from AD
 * and DP. It rewrites the fields with an entry per allele (the alleles, AT, AD, GL, GT and MAD)
 * and records the merge in INFO/MAT.
 *
 * `allele_walk` gives allele `a` of the `allele_count` alleles as handles, or returns false for an
 * allele that is to be left alone. Alleles are merged when their walks are at least `threshold`
 * similar, and only in a variant whose core length is at least `min_len`. `gl_layout` is the order
 * in which the caller wrote GL, which cannot be recovered from the record.
 */
bool merge_similar_alleles(const PathPositionHandleGraph& graph, size_t allele_count,
                           const function<bool(size_t allele, Traversal& walk)>& allele_walk,
                           vector<int>& site_genotype, const string& sample_name,
                           vcflib::Variant& out_variant, GLLayout gl_layout, double threshold,
                           int64_t min_len);

/// Trim the sequence every allele shares from the start of the alleles, or from the end if
/// `backward` is set, leaving at least one base in each. If `len_override` is not 0, trim at most
/// that many bases.
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
    int64_t ref_offset;
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
    /// The site's ID, written as the record's. When null, the ID is printed from the bounds
    /// (`site_name`).
    const string* id = nullptr;
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

/// Return the ID of the site a record belongs to.
///
/// That is the ID of a block record without the "_<index>" that tells the block records of a site
/// apart, or the ID of any other record unchanged.
string block_site_name(const string& id);

/// Return the phased GT of a nested ploidy-1 site: `allele` on strand `strand`, and "." on the other.
///
/// The other strand carries nothing here, since its parent allele deletes the chain.
string nested_strand_genotype(int allele, int strand);

/// Build the record for a site. Traversals with the same sequence become one allele. Where the
/// genotype has an allele marker, the record has "*" or ".", and a genotype with "." and no ALT
/// gets "*" as its ALT, so that the record is valid VCF.
SiteRecord build_site_record(const PathPositionHandleGraph& graph, const SiteToWrite& site,
                             const SiteAlleles& alleles, const SiteHooks& hooks,
                             const RecordOptions& options);

}

#endif
