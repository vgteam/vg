#ifndef VG_NESTED_VARIANTS_HPP_INCLUDED
#define VG_NESTED_VARIANTS_HPP_INCLUDED

/**
 * \file nested_variants.hpp
 *
 * Variants that called alleles carry along gref alt contigs, broken out as
 * records on the contigs.
 *
 * A gref alt contig is one haplotype's copy of sequence that the reference
 * lacks, such as an insertion. Another haplotype's ALT allele can run along
 * the contig and differ from it, such as an insertion with SNPs in it.
 * Comparing the walk of the allele with the contig breaks the SNPs out as
 * biallelic variants on the contig, as PanGenie's panel preparation breaks
 * variants out of walks by comparing them with the reference. The calls of
 * the alleles that read the contig give the variants their genotypes.
 */

#include "gref.hpp"
#include "vcf_alleles.hpp"

#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace vg {

/**
 * A biallelic variant on a gref alt contig.
 *
 * Its alleles are in the contig's orientation. As in VCF, an empty allele is
 * padded with the contig base before it, or with the base after it at the
 * start of the contig. Bases that REF and ALT share at either end are trimmed,
 * down to one, as PanGenie's panel preparation trims them.
 */
struct NestedVariant {
    /// The gref path of the contig.
    uint32_t contig = 0;
    /// The 1-based position of the first base of REF on the contig.
    uint64_t pos = 0;
    std::string ref;
    std::string alt;
    /// The walks of REF and ALT, from the node before the variant to the node
    /// after it.
    std::vector<handlegraph::handle_t> ref_walk;
    std::vector<handlegraph::handle_t> alt_walk;

    /// Return the 0-based contig coordinate of the first base of REF.
    uint64_t ref_start() const;
    /// Return the 0-based contig coordinate after the last base of REF.
    uint64_t ref_end() const;
    /// Return SNV, INS, DEL or COMPLEX, as PanGenie classifies variants.
    std::string type() const;
    /// Return the ID that PanGenie gives a variant, CONTIG-POS-TYPE-WALK-LENGTH,
    /// with the ALT walk.
    std::string id(const std::string& contig_name, const gbwtgraph::GBWTGraph& graph) const;
};

/// A stretch of a walk that reads a gref alt contig in order, in one orientation.
struct ContigPass {
    uint32_t contig = 0;
    /// The contig bases [start, end) that the stretch reads. It reads through
    /// to an end of the contig if it has a variant there.
    uint64_t start = 0;
    uint64_t end = 0;
    /// The steps of the walk at the first and last nodes that the stretch
    /// visits on the contig.
    size_t first_step = 0;
    size_t last_step = 0;
    /// The variants of the stretch relative to the contig, in contig order.
    std::vector<NestedVariant> variants;
};

/**
 * A decomposer that breaks walks into passes along gref alt contigs, and finds
 * the variants of each pass relative to its contig.
 *
 * Between two consecutive visits of a pass to its contig, the nodes of the
 * walk replace those of the contig, as a variant. A pass that leaves the
 * contig before its end has a variant there if a later node of the walk
 * follows the contig's last node in the graph: the walk's nodes up to that
 * node replace the contig's tail. Otherwise the walk really leaves the contig,
 * and the pass ends at its last visit. The start of the contig works the same
 * way, with a node that precedes the contig's first node.
 *
 * The reference lies outside every contig, so a walk that visits it has left
 * them all. Its next visit to a contig starts a new pass, and a variant at the
 * start or end of a contig cannot run across the reference. Only the node that
 * closes such a variant can be a reference node, next to where the contig
 * branches off. A walk on the reference that shares a node with a contig, such
 * as a common SNP allele, therefore reads only that node.
 */
class ContigDecomposer {
public:
    /// Make a decomposer for the gref paths of the index that are alt contigs.
    ContigDecomposer(const gbwtgraph::GBWTGraph& graph, const GrefIndex& gref);

    /// Return true if the gref path is an alt contig, not a copy of a
    /// reference path.
    bool is_contig(uint32_t path) const;
    /// Return the passes of the walk along alt contigs, in the order of their
    /// first visits.
    std::vector<ContigPass> decompose(const std::vector<handlegraph::handle_t>& walk) const;

private:
    /// A step of a walk onto a node of an alt contig.
    struct Visit {
        size_t step;
        uint32_t contig;
        /// True if the walk reads the contig in reverse.
        bool is_reverse;
        /// The contig bases [start, end) on the node.
        uint64_t start;
        uint64_t end;
    };

    /// A walk read in one orientation, with its steps numbered that way.
    class OrientedWalk;

    /// Return the pass made of the visits, which read one contig in order.
    ContigPass make_pass(const std::vector<handlegraph::handle_t>& walk, std::vector<Visit> visits) const;
    /**
     * Add to the pass the variant in which the walk's nodes strictly between
     * steps left and right replace the contig's bases [ref_start, ref_end).
     *
     * The walk must read the contig forward. A walk that spells the same
     * bases on other nodes has no variant there.
     */
    void add_variant(ContigPass& pass, const OrientedWalk& walk, uint64_t ref_start, uint64_t ref_end, size_t left,
                     size_t right) const;

    const gbwtgraph::GBWTGraph& graph;
    /// Finds and spells the contigs, and finds where walks go past their ends.
    const GrefIndex& gref;
};

/**
 * A catalog of the variants that the ALT alleles of a VCF carry along gref alt
 * contigs.
 *
 * It lists each variant once, on its contig, however many alleles carry it.
 * The alleles that read a contig are its readers, and the catalog keeps the
 * passes of each reader along the contig.
 */
class NestedVariantCatalog {
public:
    /// A pass of a reader along a contig.
    struct Reading {
        /// The contig bases [start, end) that the pass reads.
        uint64_t start;
        uint64_t end;
        /// The variants that the pass has, as indices into variants().
        std::vector<size_t> variants;
    };

    /// An allele that reads a contig, with its passes along it.
    struct Reader {
        /// The index of the allele in the calls.
        size_t allele;
        std::vector<Reading> passes;
    };

    /// A gref alt contig that alleles read.
    struct Contig {
        uint32_t path;
        /// The name of the contig in VCF, as vg deconstruct names it.
        std::string name;
        uint64_t length;
        /// The variants on the contig, as indices into variants(), in VCF
        /// order.
        std::vector<size_t> variants;
        std::vector<Reader> readers;
        /// One more than the number of contigs that enclose this one in the
        /// walks that read it, as in vg deconstruct's CH tag.
        size_t level;
        /// The CHROM, first POS and last REF end of the calls that read the
        /// contig, as in vg deconstruct's RC, RS and RD tags.
        std::string site_chrom;
        uint64_t site_start;
        uint64_t site_end;
        /// The index of the last record with an allele that reads the contig.
        size_t last_record;
    };

    /// Decompose the walks of the ALT alleles of the calls.
    NestedVariantCatalog(const VcfCalls& calls, const ContigDecomposer& decomposer, const GrefIndex& gref,
                         const gbwtgraph::GBWTGraph& graph);

    const std::vector<NestedVariant>& variants() const { return variant_list; }
    /// Return the ID of a variant, as PanGenie makes IDs.
    const std::string& id(size_t variant) const { return ids[variant]; }
    /// Return the contigs, in the order of their last records.
    const std::vector<Contig>& contigs() const { return contig_list; }
    /// Return the contigs that an allele reads, as pairs of an index into
    /// contigs() and an index into its readers.
    const std::vector<std::pair<size_t, size_t>>& readings(size_t allele) const { return allele_readings[allele]; }

private:
    /// The index of each variant in variants(), by its contig, position and alleles.
    using VariantKeys = std::map<std::tuple<uint32_t, uint64_t, std::string, std::string>, size_t>;

    /// Return the index of the allele among the readers of the contig, adding
    /// it as a reader if it is not one yet.
    size_t reader_of(size_t allele, size_t contig);
    /// Return the index of the variant in variants(), adding it to the catalog
    /// and to the variants of the contig if it is new.
    size_t list_variant(const NestedVariant& variant, size_t contig, VariantKeys& keys,
                        const gbwtgraph::GBWTGraph& graph);
    /// Put the contigs in the order of their last records.
    void sort_contigs();

    std::vector<NestedVariant> variant_list;
    std::vector<std::string> ids;
    std::vector<Contig> contig_list;
    /// For each allele of the calls, the contigs that it reads.
    std::vector<std::vector<std::pair<size_t, size_t>>> allele_readings;
};

/**
 * A genotyper that writes the variants of a catalog as VCF records, with
 * genotypes from the calls of the alleles that read their contigs.
 *
 * Each copy of a called ALT allele is a walk that the sample carries. A walk
 * that has a variant gives the sample a 1 there, and a walk that reads the
 * variant's part of the contig without it gives a 0, even if it has another
 * variant there. A sample with no such walk gets "." for the haplotype, as in
 * vg deconstruct's records on alt contigs, so allele frequencies are among the
 * carriers of the contig. Phased calls give phased genotypes. A sample with
 * more walks there than haplotypes gets a missing genotype.
 *
 * A genotype is only as good as the calls it comes from, so its GQ is the
 * lowest GQ of the calls whose walks give it, if they all have one. A record's
 * PS names the top-level snarl of its contig, the site that the contig is an
 * allele of.
 *
 * Records go out a contig at a time, as soon as no later call reads the
 * contig, so the genotyper holds the walks of a few contigs at most.
 */
class NestedGenotyper {
public:
    /**
     * Make a genotyper for the calls that writes to out.
     *
     * @param extra_info Extra INFO for the record of each variant, starting
     * with ";" if not empty.
     * @param top_level_snarls The top-level snarl of each contig, by its name
     * in the VCF, for INFO/PS. Records on contigs that it leaves out get no PS.
     */
    NestedGenotyper(const NestedVariantCatalog& catalog, const VcfCalls& calls, std::vector<std::string> extra_info,
                    const std::unordered_map<std::string, std::string>& top_level_snarls,
                    const gbwtgraph::GBWTGraph& graph, std::ostream& out);

    /// Return the header lines that declare the INFO keys, the FORMAT key and
    /// the contigs of the records. Contigs without variants are left out.
    std::vector<std::string> header_lines() const;
    /// Take the genotypes of a call, then write the records of the contigs
    /// that no later call reads.
    void add_call(size_t record, const VcfCursor& cursor);
    /// Return the number of genotypes set missing because the sample carries
    /// more walks there than it has haplotypes.
    size_t conflicts() const { return conflict_count; }
    /// Return the number of GT alleles that named an ALT allele the record does
    /// not have, which are read as REF.
    size_t alleles_past_end() const { return past_end_count; }

private:
    /// A walk that a sample carries along a contig.
    struct Carried {
        /// The index of its allele among the readers of the contig.
        size_t reader;
        /// The haplotype that carries it, or -1 if the call is not phased.
        int haplotype;
        /// The GQ of the call, or -1 if it has none.
        int quality;
    };

    /// A sample's genotype at a variant.
    class Genotype;
    /// The sample columns of a record, with the counts for its INFO.
    class SampleColumns;

    /**
     * Add the walks of the ALT alleles that a sample's genotype calls to the
     * contigs that they read.
     *
     * @param quality The GQ of the call, or -1 if it has none.
     * @returns The number of haplotypes in the genotype.
     */
    size_t carry_genotype(const VcfSite& site, size_t sample, std::string_view genotype, int quality);
    /**
     * Add the walk of an ALT allele that a sample carries to each contig that
     * the allele reads.
     *
     * @param haplotype The haplotype that carries the walk, or -1 if the call
     * is not phased.
     * @param quality The GQ of the call, or -1 if it has none.
     */
    void carry(size_t sample, size_t allele, int haplotype, int quality);
    /// Write the records of the variants on a contig, and forget its walks.
    void write_contig(size_t contig);
    /// Write the record of a variant on the contig, with the sample columns.
    void write_record(const NestedVariantCatalog::Contig& contig, size_t variant, const SampleColumns& columns);

    const NestedVariantCatalog& catalog;
    const VcfCalls& calls;
    /// Extra INFO for the record of each variant.
    std::vector<std::string> extra_info;
    const std::unordered_map<std::string, std::string>& top_level_snarls;
    const gbwtgraph::GBWTGraph& graph;
    std::ostream& out;
    /// For each contig, the walks that each sample carries along it. Empty
    /// until a call reads the contig.
    std::vector<std::vector<std::vector<Carried>>> carried;
    /// The number of haplotypes of each sample, from its genotypes.
    std::vector<size_t> ploidy;
    /// The index of the next contig to write.
    size_t next_contig = 0;
    size_t conflict_count = 0;
    size_t past_end_count = 0;
};

}

#endif
