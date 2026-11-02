#ifndef VG_SITE_EMITTER_HPP_INCLUDED
#define VG_SITE_EMITTER_HPP_INCLUDED

/** \file
 * Building the VCF record of a called site and adding it to a sorted VCF.
 */

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Variant.h"
#include "handle.hpp"
#include "sorted_vcf_writer.hpp"
#include "../artifacts/vcf_genotype_likelihoods.hpp"
#include "../artifacts/vcf_record.hpp"

namespace vg {

using namespace std;

/// Represents a called site for translation into VCF.
struct CalledSite {
    /// The site's boundary nodes, oriented into the site.
    handle_t start;
    handle_t end;
    /// The reference path the record is placed on, which may name a subrange of its base path.
    const string& ref_path_name;
    /// The offset added to positions on the reference path.
    int64_t ref_offset;
    /// The called genotype, as numbers of the site's traversals or as the allele markers. Empty
    /// when the site has no call.
    const vector<int>& genotype;
    /// The traversal that is the reference allele.
    int ref_trav_idx;
    /// How many traversals the site has. The allele source can spell any of them.
    size_t traversal_count;
    /// The ploidy at the site.
    int ploidy;
    /// The record's ID, or null to name the record automatically.
    const string* id = nullptr;
};

/// Spell traversal `allele` of a site: its sequence, and its visits. The visits are empty for a
/// traversal with no walk through the graph.
using allele_source_t = function<pair<string, vector<NodeVisit>>(size_t allele)>;

/**
 * Write the genotyper's fields into a record, given which traversal each VCF allele is spelled
 * from (-1 for "*") and the genotype in VCF allele numbers. It is expected to fill:
 *
 * - AD, the read support of each allele.
 * - GL, the genotype likelihoods, in the emitter's `gl_layout` order.
 * - DP and GQ.
 * - Any other INFO or FORMAT field the genotyper reports, and FILTER.
 */
using field_writer_t = function<void(const vector<int>& allele_traversals,
                                     const vector<int>& genotype, vcflib::Variant& record)>;

/**
 * An emitter that accepts CalledSites and produces VCF records into a SortedVCFWriter.
 *
 * It can be subclassed to customize its VCF record generation behavior.
 */
class SiteEmitter {
public:
    /// Make an emitter that writes into `output`, which must outlive it, and declare the AT
    /// field in the header of `output`.
    explicit SiteEmitter(SortedVCFWriter& output);
    virtual ~SiteEmitter() = default;

    /// Write a record even for a reference call, with up to `max_uncalled_alleles` alleles that
    /// were not called added, and trim the alleles only as far as the site's boundary nodes.
    bool genotype_snarls = false;
    /// Merge called ALT alleles whose similarity, a modified length-weighted Jaccard, is at least
    /// this. 1.0 merges none.
    double allele_merge_threshold = 1.0;
    /// Merge alleles only in a record whose core length is at least this.
    int64_t allele_merge_min_len = 0;
    /// The order that will be used in GL fields.
    GLLayout gl_layout = GLLayout::IMajor;

    /// How many alleles that were not called a reference call written under `genotype_snarls`
    /// gets at most.
    static constexpr size_t max_uncalled_alleles = 5;

    /**
     * Write the record of `site`.
     *
     * Traversals that spell the same sequence become one allele. Where the genotype has an
     * allele marker, the record has "*" or ".", and a genotype with "." and no ALT gets "*" as
     * its ALT, so that the record is valid VCF. A record with no ALT is hidden unless
     * `genotype_snarls` is set.
     *
     * @param alleles Spells the site's traversals.
     * @param write_fields Writes the genotyper's fields once the alleles are numbered.
     * @return false only if the writer refused a line the site wanted.
     */
    bool emit(const PathPositionHandleGraph& graph, const CalledSite& site,
              const allele_source_t& alleles, const field_writer_t& write_fields);

protected:
    /// Return whether traversal `trav`, which is not the reference traversal, should be treated
    /// as the same as the reference allele. By default no other traversal should be.
    virtual bool same_as_reference(const CalledSite& site, int trav);

    /// Decide how the genotype should be phased. Return the phase set it should be in, or -1
    /// for none, and replace the GT string `gt` with the phased one if it should be phased.
    /// `genotype` is in VCF allele numbers, and `trav_to_allele` gives the VCF allele of each
    /// traversal in it. By default the genotype should stay unphased, in no phase set.
    virtual int64_t phase(const CalledSite& site, const vector<int>& genotype,
                          const map<int, int>& trav_to_allele, string& gt);

    /// Decide whether the site should be split into several records. If it should, write them
    /// and return the number of lines written. If it should not, write nothing and return -1,
    /// and the site gets one record. The default implementation never splits a site.
    virtual int write_blocks(const PathPositionHandleGraph& graph, const CalledSite& site,
                             const SiteRecord& record);

    /// Change the finished record just before it is written. By default it is left alone.
    virtual void finish_record(const CalledSite& site, vcflib::Variant& record);

    /**
     * Handle a site whose records have all been sent to the writer. By default do nothing.
     *
     * @param trav_to_allele If the site got one record, the VCF allele of each called traversal
     * in it. If `write_blocks` split the site, empty, since each record numbers its own alleles.
     * @param has_line If the site got one record, whether that record was written rather than
     * hidden or refused. If the site was split, whether any of its records was written.
     */
    virtual void on_site_done(const CalledSite& site, const map<int, int>& trav_to_allele,
                              bool has_line);

    /// The output the records go to.
    SortedVCFWriter& output;

private:
    /// The spelling of each traversal that `emit` has asked the allele source for, by traversal.
    using Spellings = map<size_t, pair<string, vector<NodeVisit>>>;

    /// Return the spelling of traversal `trav`, asking `alleles` for it only the first time.
    static const pair<string, vector<NodeVisit>>& spelling(const allele_source_t& alleles,
                                                           Spellings& spelled, size_t trav);

    /**
     * Number the reference allele and the alleles of the called genotype of `site`, whose REF is
     * already in `record`. Traversals that spell the same sequence share one allele, and the
     * reference allele is number 0.
     *
     * Fill the empty `allele_to_gt` with the VCF allele number of each allele sequence, and the
     * empty `site_trav` with the traversal each VCF allele is spelled from, or -1 for "*". Fill in
     * `record.genotype` and `record.trav_to_allele` to match.
     */
    void number_called_alleles(const CalledSite& site, const allele_source_t& alleles,
                               Spellings& spelled, SiteRecord& record,
                               map<string, int>& allele_to_gt, vector<int>& site_trav);

    /// Add up to `max_uncalled_alleles` alleles of `site` that were not called, in order of their
    /// sequences, to `allele_to_gt` and `site_trav`.
    static void add_uncalled_alleles(const CalledSite& site, const allele_source_t& alleles,
                                     Spellings& spelled, map<string, int>& allele_to_gt,
                                     vector<int>& site_trav);

    /// Build the record of `site`, up to the point where it could be split into blocks.
    SiteRecord build_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                            const allele_source_t& alleles, const field_writer_t& write_fields);

    /// Send a built record to the writer, unless `write_blocks` splits the site, and report the
    /// site to `on_site_done`.
    ///
    /// @return false only if the writer refused a line the site wanted.
    bool file_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                     SiteRecord& record);
};

}

#endif
