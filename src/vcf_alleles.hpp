#ifndef VG_VCF_ALLELES_HPP_INCLUDED
#define VG_VCF_ALLELES_HPP_INCLUDED

/**
 * \file vcf_alleles.hpp
 *
 * The ALT alleles of a VCF of graph variants, as walks in the graph.
 *
 * Tools that describe variants in a graph record the walk of each allele. vg
 * deconstruct writes the walk of every allele of a record in INFO/AT.
 * PanGenie's panel preparation splits sites into biallelic variants named
 * CHROM-POS-TYPE-WALK-LENGTH, with the ALT walk in the name, and PanGenie keeps
 * that name as the ID of each call.
 *
 * The first and last nodes of each walk lie on the reference, on either side
 * of the variant. For vg deconstruct they are the boundary nodes of the
 * variant's snarl. The allele is the walk without them. It spells the ALT
 * sequence without its padding base, and with any bases that a tool trimmed
 * from both REF and ALT.
 */

#include "walk_interval.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <array>
#include <functional>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace vg {

/// Return true if htslib reads the file as VCF or BCF.
bool is_variant_file(const std::string& path);

/**
 * A cursor over the lines of a VCF file, plain or compressed.
 *
 * On a record line, it splits off the eight fixed columns, from CHROM to INFO.
 * It leaves the sample columns as they are, so a VCF with many samples reads
 * quickly and can be copied exactly.
 */
class VcfCursor {
public:
    /// The fixed columns of a record line, in order.
    enum Column { CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FIXED_COLUMNS };

    /**
     * Open the VCF and decompress it with the given number of threads.
     *
     * @throws std::runtime_error if the file cannot be opened or is BCF.
     */
    VcfCursor(const std::string& path, int threads);
    ~VcfCursor();
    VcfCursor(const VcfCursor&) = delete;
    VcfCursor& operator=(const VcfCursor&) = delete;

    /**
     * Move to the next line, skipping empty ones.
     *
     * @returns false at the end of the file.
     * @throws std::runtime_error on a read error, or on a record line with
     * fewer than eight columns.
     */
    bool next();

    /// Return true if the line is a header line, which starts with #.
    bool at_header() const;
    /// Return the line, without its line break.
    std::string_view line() const;
    /// Return a fixed column of a record line.
    std::string_view column(Column column) const;
    /// Return the rest of a record line after its INFO column, from the tab
    /// before FORMAT on. It is empty in a VCF without samples.
    std::string_view after_info() const;
    /// Return the ALT alleles of a record line. ALT "." has none.
    std::vector<std::string_view> alt_alleles() const;
    /// Return the value of an INFO key on a record line, or nothing if the
    /// record lacks the key. A flag has an empty value.
    std::optional<std::string_view> info_value(std::string_view key) const;
    /// Return the GT of each sample on a record line, in sample order. Return
    /// none if the record has no samples, or if FORMAT does not start with GT.
    std::vector<std::string_view> genotypes() const;
    /// Return the place of the record, as CHROM:POS.
    std::string locus() const;

private:
    htsFile* file = nullptr;
    /// The line, read by htslib.
    kstring_t buffer = KS_INITIALIZE;
    /// The [start, end) offsets of each fixed column in a record line.
    std::array<std::pair<size_t, size_t>, FIXED_COLUMNS> columns;
};

/// An ALT allele of a VCF record, as a walk in the graph.
struct VcfAllele {
    /// The ID of the record, followed by :N for the Nth ALT allele of a record
    /// with several.
    std::string name;
    /// The index of the allele's record among the records of the VCF.
    size_t record = 0;
    /// The walk of the allele, from the node before the variant to the node
    /// after it. Empty if the VCF does not give it.
    std::vector<handlegraph::handle_t> walk;
    /// The bases of the walk without its first and last nodes, if any nodes
    /// are left. A deletion's walk has none.
    std::optional<WalkInterval> interval;
    /// True if the sequence of the interval holds the ALT column, apart from a
    /// padding base. A walk from another graph usually does not.
    bool holds_alt = true;
};

/// The place of a VCF record on its reference.
struct VcfSite {
    std::string chrom;
    /// The 1-based position of the first base of REF.
    uint64_t pos = 0;
    uint64_t ref_length = 0;
    /// The index of the record's first ALT allele among the ALT alleles of the
    /// VCF. The record's other ALT alleles follow it.
    size_t first_allele = 0;
    size_t alt_count = 0;
};

/// The records of a VCF of graph variants, with their ALT alleles as walks.
struct VcfCalls {
    /// The sample names, in column order.
    std::vector<std::string> samples;
    /// The place of each record, in file order.
    std::vector<VcfSite> sites;
    /// The ALT alleles of all the records, in file order.
    std::vector<VcfAllele> alleles;
};

/**
 * Read the records of a VCF of graph variants, and the walks of their ALT
 * alleles.
 *
 * Each allele's walk comes from INFO/AT. A record without AT and with one ALT
 * allele takes the walk from its ID if that is a PanGenie variant ID.
 *
 * @throws std::runtime_error on a malformed VCF, on a walk that is not in the
 * graph, or on a walk with fewer than two nodes.
 */
VcfCalls read_vcf_calls(const std::string& path, const gbwtgraph::GBWTGraph& graph, int threads);

/// An INFO key with one value for each ALT allele of a VCF (Number=A).
struct VcfAlleleInfo {
    std::string key;
    /// The ##INFO line that declares the key.
    std::string header;
    /// The value for each ALT allele of the VCF, in file order. A "." is a
    /// missing value.
    std::vector<std::string> values;
};

/**
 * Copy a VCF to out with INFO keys added.
 *
 * The keys replace any values and declarations that the VCF has for them. A
 * record gets no value for a key if all its alleles have "." there.
 *
 * @param header Header lines to add, each unless the VCF already declares its
 * INFO key, FORMAT key or contig.
 * @param listener Called with the index of each record and a cursor on it,
 * after the record is written.
 * @throws std::runtime_error if the VCF does not have one ALT allele for each
 * value.
 */
void write_vcf_with_allele_info(const std::string& path, std::ostream& out, const std::vector<VcfAlleleInfo>& info,
                                const std::vector<std::string>& header,
                                const std::function<void(size_t, const VcfCursor&)>& listener, int threads);

}

#endif
