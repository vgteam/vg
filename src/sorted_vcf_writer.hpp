#ifndef VG_SORTED_VCF_WRITER_HPP_INCLUDED
#define VG_SORTED_VCF_WRITER_HPP_INCLUDED

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "Variant.h"
#include "handle.hpp"
#include "region.hpp"
#include "../artifacts/vcf_record.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/// Return the stretch of the reference that the REF of `record` covers, 1-based and end-inclusive,
/// as VCF POS is.
Region record_ref_region(const vcflib::Variant& record);

/**
 * A sorted VCF, built up from records and written out at the end.
 *
 * Records can be added from any thread. This class is responsible for buffering them until the
 * end, giving them the nesting INFO tags, sorting them, and writing them.
 */
class SortedVCFWriter {
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

    explicit SortedVCFWriter(const string& sample_name);

    /**
     * Add a record to the output. The record is attached to the writer's VCF file, which gives it
     * the writer's sample columns, and is then kept as its text.
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
    bool add_variant(vcflib::Variant& var, size_t block = 0, bool hidden = false);

    /// A contig of the reference, for a ##contig line of the header.
    struct ReferenceContig {
        /// The name that the CHROM column of the records gives the contig.
        string name;
        size_t length = 0;
    };

    /**
     * Set the contigs that the header declares, in order.
     *
     * @param only_with_records Declare only the contigs that have records.
     */
    void set_contigs(vector<ReferenceContig> declared, bool only_with_records = false);

    /// Add one ##INFO line, without its newline, to the header. The added INFO lines follow the
    /// nesting lines, in the order they were added.
    void add_info_line(string line);
    /// Add one ##FORMAT line, without its newline, to the header. The added FORMAT lines follow
    /// the GT line, in the order they were added.
    void add_format_line(string line);
    /// Add one ##FILTER line, without its newline, to the header. The added FILTER lines follow
    /// the PASS line, in the order they were added.
    void add_filter_line(string line);

    /**
     * Write the header, then the records, sorted, and empty the buffer. Call it once, after
     * every record and header line has been added.
     *
     * With nesting on, each record first gets the nesting INFO tags, which say where its site
     * lies in the snarl tree `sites` and where to find it on the reference. A record is matched
     * to its site by its ID, so `record_id_of_site` must give each snarl's ID as its records
     * carry it. When nesting is off, `sites` may be null.
     *
     * @param graph The graph of `sites`, which the default record IDs read node IDs from. It
     * may be null when `record_id_of_site` is given or nesting is off.
     * @param record_id_of_site Returns the record ID of a snarl of `sites`, read the way the
     * net handle traverses it. When empty, a snarl is named from its boundary nodes.
     */
    void write_variants(ostream& out, const SnarlDecomposition* sites = nullptr,
                        const HandleGraph* graph = nullptr,
                        const function<string(const net_handle_t&)>& record_id_of_site = {});

    /// Write the nesting INFO tags (LV, CH, PS, RC, RS, RD) to every record.
    void set_nested(bool nested);
    /// Set the gRef nesting level of each gRef contig, by contig name.
    ///
    /// A record on one of these contigs gets at least its level as INFO/CH.
    void set_gref_levels(map<string, int> levels);

    /// Name nodes in record IDs and AT by `translation` instead of by their IDs.
    void set_translation(const NodeTranslation* translation);
    const NodeTranslation* translation() const;

    const string& sample_name() const;

    /// Write the sample columns of the records for `names`, in order, instead of for the one
    /// sample the writer was made for.
    void set_sample_names(vector<string> names);

    /// Run vcflib's vcffixup on a record, filling AC, AN and the like.
    void vcf_fixup(vcflib::Variant& var) const;

    /// The nesting INFO header lines, for both vg call and vg deconstruct.
    static string nesting_info_headers();

    /// A record longer than this is refused.
    static const int64_t max_vcf_line_length = 2000000000;

private:
    /**
     * Return the text of the whole header, each line ending in a newline. If `keep_contigs`
     * is not null, only the contigs in it get ##contig lines.
     *
     * The header starts with the fileformat line and the ##contig lines. The INFO, FORMAT and
     * FILTER lines follow, grouped by type. Each group starts with the lines that the writer
     * always writes, and then has the added lines. The ##SAMPLE lines and the column line come
     * last.
     */
    string vcf_header(const unordered_set<string>* keep_contigs) const;

    /// Return the reference contigs that have any records in the buffer.
    unordered_set<string> get_output_contigs() const;

    /// Add the nesting INFO tags to every buffered record, from the snarl tree `sites`, whose
    /// sites have the record IDs `record_id_of_site` gives.
    void update_nesting_info_tags(const SnarlDecomposition& sites,
                                  const function<string(const net_handle_t&)>& record_id_of_site);

    /// The name of the sample the writer was made for.
    string sample;

    /**
     * The VCF file that the records belong to. Its sample names are the samples whose columns
     * the records have, in order.
     *
     * We take a deliberate shortcut with this file. vcflib expects a file to be opened with its
     * header before any of its records exist. Here the records are printed to text as they are
     * added, and the header is known only at the end. So while records are added, the file holds
     * only the sample names. A record attached to the file copies them, and printing the record
     * reads only those copies. The shortcut is safe because nothing reads the INFO or FORMAT
     * types of the file through a record before the header is in place: a record is kept only as
     * text once it is added. The file is opened with the real header just before anything is
     * written.
     */
    vcflib::VariantCallFile output_vcf;

    /// The contigs for the ##contig lines, in order.
    vector<ReferenceContig> contigs;
    /// Whether the header declares only the contigs that have records.
    bool only_contigs_with_records = false;

    /// The buffered records, one buffer per thread, each a sort key and the zstd-compressed line.
    /// A vcflib::Variant is too big to keep in memory for every record.
    vector<vector<pair<BufferedRecordKey, string>>> output_variants;

    /// The reference interval of each hidden record, by its ID, one table per thread like
    /// `output_variants`.
    ///
    /// A hidden record is one of a site that had no variant to report. A record nested under such
    /// a site may have no other reference coordinate to point its RC, RS and RD at. This is
    /// common in gref graphs, where the parent of an island of non-reference sequence is often a
    /// large snarl that only the reference and its own gref copy span.
    vector<unordered_map<string, Region>> suppressed_ref_info;

    /// The ##INFO lines added to the header.
    vector<string> info_lines;
    /// The ##FORMAT lines added to the header.
    vector<string> format_lines;
    /// The ##FILTER lines added to the header.
    vector<string> filter_lines;

    /// The node translation that names nodes, or null to name them by ID.
    const NodeTranslation* node_translation = nullptr;

    /// Whether to write the nesting INFO tags.
    bool include_nested = false;
    /// The gRef nesting level of each gRef contig, by contig name. Empty without a gRef cover.
    map<string, int> gref_levels;
};

}

#endif
