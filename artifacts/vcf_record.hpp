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
#include "../src/handle.hpp"
#include "../src/traversal_finder.hpp"
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
 * `allele_visits` gives the visits of each allele, empty for "*". An allele whose visits make
 * no walk (see visits_to_walk) is left alone. Alleles are merged when their walks are at least `threshold`
 * similar, and only in a variant whose core length is at least `min_len`. `gl_layout` is the order
 * in which the caller wrote GL, which cannot be recovered from the record.
 */
bool merge_similar_alleles(const PathPositionHandleGraph& graph,
                           const vector<vector<NodeVisit>>& allele_visits,
                           vector<int>& site_genotype, const string& sample_name,
                           vcflib::Variant& out_variant, GLLayout gl_layout, double threshold,
                           int64_t min_len);

/// Trim the sequence every allele shares from the start of the alleles, or from the end if
/// `backward` is set, leaving at least one base in each. If `len_override` is not 0, trim at most
/// that many bases.
void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override);

/// A called site, as `SiteEmitter` writes it.
struct CalledSite;

/// A record built for a site, with what its caller needs to file it.
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

/// Return the GT text of an unphased genotype in VCF allele numbers, with "." for
/// MISSING_ALLELE_MARKER.
string unphased_gt_text(const vector<int>& site_genotype);

/// Return the GT text of a genotype with no call: `ploidy` "." separated by "/".
string no_call_gt_text(int ploidy);

/// Set the contig, position and ID of a record for `site`.
void place_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                  const NodeTranslation* translation, SiteRecord& record);

/// Add "*" as the ALT of a record whose genotype has "." and no ALT, so that it is valid VCF.
///
/// This happens where a strand of the parent does not visit a nested site.
void add_star_for_missing(vcflib::Variant& out_variant, const vector<int>& site_genotype,
                          const string& sample_name);

/// Return the ID of the site a record belongs to.
///
/// That is the ID of a block record without the "_<index>" that tells the block records of a site
/// apart, or the ID of any other record unchanged.
string block_site_name(const string& id);

/// Return the phased GT of a nested ploidy-1 site: `allele` on strand `strand`, and "." on the other.
///
/// The other strand carries nothing here, since its parent allele deletes the chain.
string nested_strand_genotype(int allele, int strand);

}

#endif
