/**
 * \file nested_variants.cpp
 * Variants on gref alt contigs, see nested_variants.hpp.
 */

#include "nested_variants.hpp"
#include "sorted_vcf_writer.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>
#include <tuple>
#include <unordered_map>

#include <omp.h>

namespace vg {

using handlegraph::handle_t;
using handlegraph::PathMetadata;
using gbwtgraph::GBWTGraph;

//------------------------------------------------------------------------------

uint64_t NestedVariant::ref_start() const {
    return pos - 1;
}

uint64_t NestedVariant::ref_end() const {
    return pos - 1 + ref.size();
}

std::string NestedVariant::type() const {
    if (ref.size() == 1) {
        return alt.size() == 1 ? "SNV" : "INS";
    }
    return alt.size() == 1 ? "DEL" : "COMPLEX";
}

std::string NestedVariant::id(const std::string& contig_name, const GBWTGraph& graph) const {
    // A complex variant is as long as its longer allele. Another is as long as
    // what it adds or removes, and a SNV is 1 long.
    size_t length = std::max<size_t>(ref.size() > alt.size() ? ref.size() - alt.size() : alt.size() - ref.size(), 1);
    if (ref.size() > 1 && alt.size() > 1) {
        length = std::max(ref.size(), alt.size());
    }
    return contig_name + "-" + std::to_string(pos) + "-" + type() + "-" + walk_text(alt_walk, graph) + "-" +
           std::to_string(length);
}

//------------------------------------------------------------------------------

ContigDecomposer::ContigDecomposer(const GBWTGraph& graph, const GrefIndex& gref) : graph(graph), gref(gref) {
}

bool ContigDecomposer::is_contig(uint32_t path) const {
    return path < gref.path_count() && gref.is_alt_contig(path);
}

std::vector<ContigPass> ContigDecomposer::decompose(const std::vector<handle_t>& walk) const {
    // Group the walk's visits to each contig into passes. A visit continues
    // the contig's latest pass if it reads the contig the same way, further on.
    std::vector<std::vector<Visit>> groups;
    std::unordered_map<uint32_t, size_t> latest;
    for (size_t step = 0; step < walk.size(); step++) {
        const GrefIndex::Position& place = gref.find(graph.get_id(walk[step]));
        if (place.path == GrefIndex::NO_PATH) {
            continue;
        }
        if (!is_contig(place.path)) {
            // The walk is back on the reference, so it has left every contig.
            latest.clear();
            continue;
        }
        Visit visit { step, place.path, graph.get_is_reverse(walk[step]) != place.is_reverse, place.offset,
                      place.offset + graph.get_length(walk[step]) };
        if (auto found = latest.find(visit.contig); found != latest.end()) {
            const Visit& last = groups[found->second].back();
            if (last.is_reverse == visit.is_reverse &&
                (visit.is_reverse ? visit.end <= last.start : visit.start >= last.end)) {
                groups[found->second].push_back(visit);
                continue;
            }
        }
        latest[visit.contig] = groups.size();
        groups.push_back({ visit });
    }

    std::vector<ContigPass> passes;
    passes.reserve(groups.size());
    for (std::vector<Visit>& visits : groups) {
        passes.push_back(make_pass(walk, std::move(visits)));
    }
    return passes;
}

/**
 * A walk read in one orientation, with its steps numbered that way.
 *
 * Read in reverse, the walk's last step comes first, with its handle flipped.
 */
class ContigDecomposer::OrientedWalk {
public:
    OrientedWalk(const GBWTGraph& graph, const std::vector<handle_t>& walk, bool is_reverse) :
        graph(graph), walk(walk), is_reverse(is_reverse) {
    }

    /// Return the handle that this reading visits at step i.
    handle_t at(size_t i) const {
        return is_reverse ? graph.flip(walk[to_stored(i)]) : walk[i];
    }
    /// Return the number in the stored walk of step i of this reading.
    size_t to_stored(size_t i) const {
        return is_reverse ? walk.size() - 1 - i : i;
    }
    /// Return the number in this reading of step i of the stored walk.
    size_t from_stored(size_t i) const {
        return is_reverse ? walk.size() - 1 - i : i;
    }

private:
    const GBWTGraph& graph;
    const std::vector<handle_t>& walk;
    bool is_reverse;
};

namespace {

/// Trim the bases that REF and ALT share off the end, then off the start,
/// keeping one base in each, and move POS past the bases trimmed off the start.
void trim_shared_bases(NestedVariant& variant) {
    std::string& ref = variant.ref;
    std::string& alt = variant.alt;
    while (ref.size() >= 2 && alt.size() >= 2 && ref.back() == alt.back()) {
        ref.pop_back();
        alt.pop_back();
    }
    size_t shared = 0;
    while (ref.size() - shared >= 2 && alt.size() - shared >= 2 && ref[shared] == alt[shared]) {
        shared++;
    }
    ref.erase(0, shared);
    alt.erase(0, shared);
    variant.pos += shared;
}

}

ContigPass ContigDecomposer::make_pass(const std::vector<handle_t>& walk, std::vector<Visit> visits) const {
    ContigPass pass;
    pass.contig = visits.front().contig;
    pass.first_step = visits.front().step;
    pass.last_step = visits.back().step;
    uint64_t length = gref.path_length(pass.contig);

    // Read the walk in the contig's orientation, with steps numbered that way.
    bool is_reverse = visits.front().is_reverse;
    OrientedWalk oriented(graph, walk, is_reverse);
    if (is_reverse) {
        std::reverse(visits.begin(), visits.end());
        for (Visit& visit : visits) {
            visit.step = oriented.from_stored(visit.step);
        }
    }
    // Ask the index where the walk comes past the start of the contig, and
    // where it rejoins the graph past the end. The index searches the walk as
    // stored, so a pass that reads the contig in reverse asks about the other
    // end, and gets back a step to renumber the contig's way.
    auto step_before_start = [&](size_t i) {
        std::optional<size_t> found =
            is_reverse ? gref.step_past_contig_end(walk, oriented.to_stored(i), pass.contig, true)
                       : gref.step_before_contig_start(walk, i, pass.contig, false);
        return found ? std::optional<size_t>(oriented.from_stored(*found)) : std::nullopt;
    };
    auto step_past_end = [&](size_t i) {
        std::optional<size_t> found =
            is_reverse ? gref.step_before_contig_start(walk, oriented.to_stored(i), pass.contig, true)
                       : gref.step_past_contig_end(walk, i, pass.contig, false);
        return found ? std::optional<size_t>(oriented.from_stored(*found)) : std::nullopt;
    };

    // A pass that joins the contig after its start replaces the contig's head
    // if the walk comes past the start.
    const Visit& head = visits.front();
    pass.start = head.start;
    std::optional<size_t> before;
    if (head.start > 0) {
        before = step_before_start(head.step);
    }
    if (before.has_value()) {
        add_variant(pass, oriented, 0, head.start, *before, head.step);
        pass.start = 0;
    }
    // Between two visits, the walk's nodes replace the contig's, unless the
    // two nodes follow each other on both.
    for (size_t k = 0; k + 1 < visits.size(); k++) {
        const Visit& left = visits[k];
        const Visit& right = visits[k + 1];
        if (left.end < right.start || right.step > left.step + 1) {
            add_variant(pass, oriented, left.end, right.start, left.step, right.step);
        }
    }
    // A pass that leaves the contig before its end replaces the contig's tail
    // if the walk rejoins the graph past the end.
    const Visit& tail = visits.back();
    pass.end = tail.end;
    std::optional<size_t> past;
    if (tail.end < length) {
        past = step_past_end(tail.step);
    }
    if (past.has_value()) {
        add_variant(pass, oriented, tail.end, length, tail.step, *past);
        pass.end = length;
    }
    return pass;
}

void ContigDecomposer::add_variant(ContigPass& pass, const OrientedWalk& walk, uint64_t ref_start, uint64_t ref_end,
                                   size_t left, size_t right) const {
    NestedVariant variant;
    variant.contig = pass.contig;
    std::string ref_bases;
    std::string alt_bases;
    variant.ref_walk.push_back(walk.at(left));
    if (ref_start < ref_end) {
        WalkInterval replaced = gref.interval(pass.contig, ref_start, ref_end);
        ref_bases = replaced.sequence(graph);
        variant.ref_walk.insert(variant.ref_walk.end(), replaced.walk.begin(), replaced.walk.end());
    }
    variant.ref_walk.push_back(walk.at(right));
    variant.alt_walk.push_back(walk.at(left));
    for (size_t i = left + 1; i < right; i++) {
        variant.alt_walk.push_back(walk.at(i));
        alt_bases += graph.get_sequence(walk.at(i));
    }
    variant.alt_walk.push_back(walk.at(right));
    if (ref_bases == alt_bases) {
        // The walk spells the contig on other nodes. It has no variant
        // here, but it still reads this part of the contig.
        return;
    }

    // Pad with the base before, or after at the start of the contig.
    if (ref_start > 0) {
        char before = gref.base(pass.contig, ref_start - 1);
        variant.ref = before + ref_bases;
        variant.alt = before + alt_bases;
        variant.pos = ref_start;
    } else {
        char after = gref.base(pass.contig, ref_end);
        variant.ref = ref_bases + after;
        variant.alt = alt_bases + after;
        variant.pos = 1;
    }
    trim_shared_bases(variant);
    pass.variants.push_back(std::move(variant));
}

//------------------------------------------------------------------------------

namespace {

/// Return the level of the contig that a pass of a walk reads. It is one more
/// than the number of the walk's passes along other contigs that enclose the pass.
size_t nesting_level(const std::vector<ContigPass>& passes, const ContigPass& pass) {
    auto enclosing = std::count_if(passes.begin(), passes.end(), [&](const ContigPass& other) {
        return other.contig != pass.contig && other.first_step < pass.first_step && pass.last_step < other.last_step;
    });
    return 1 + static_cast<size_t>(enclosing);
}

}

NestedVariantCatalog::NestedVariantCatalog(const VcfCalls& calls, const ContigDecomposer& decomposer,
                                           const GrefIndex& gref, const GBWTGraph& graph) :
    allele_readings(calls.alleles.size()) {
    // Decompose the walks in parallel, then list each variant once, by its
    // contig, position and alleles.
    std::vector<std::vector<ContigPass>> passes(calls.alleles.size());
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t allele = 0; allele < calls.alleles.size(); allele++) {
        if (!calls.alleles[allele].walk.empty()) {
            passes[allele] = decomposer.decompose(calls.alleles[allele].walk);
        }
    }

    std::unordered_map<uint32_t, size_t> contig_of_path;
    VariantKeys variant_keys;
    for (size_t allele = 0; allele < passes.size(); allele++) {
        size_t record = calls.alleles[allele].record;
        const VcfSite& site = calls.sites[record];
        for (const ContigPass& pass : passes[allele]) {
            size_t level = nesting_level(passes[allele], pass);
            auto [found, is_new] = contig_of_path.try_emplace(pass.contig, contig_list.size());
            if (is_new) {
                Contig& contig = contig_list.emplace_back();
                contig.path = pass.contig;
                contig.name = PathMetadata::parse_locus_name(gref.path_name(pass.contig));
                contig.length = gref.path_length(pass.contig);
                contig.level = level;
                contig.site_chrom = site.chrom;
                contig.site_start = site.pos;
                contig.site_end = site.pos + site.ref_length;
                contig.last_record = record;
            }
            size_t contig_index = found->second;
            Contig& contig = contig_list[contig_index];
            contig.level = std::min(contig.level, level);
            contig.site_start = std::min(contig.site_start, site.pos);
            contig.site_end = std::max(contig.site_end, site.pos + site.ref_length);
            contig.last_record = std::max(contig.last_record, record);

            // A walk can pass along a contig more than once, as one reader.
            Reading& reading = contig.readers[reader_of(allele, contig_index)].passes.emplace_back();
            reading.start = pass.start;
            reading.end = pass.end;
            for (const NestedVariant& variant : pass.variants) {
                reading.variants.push_back(list_variant(variant, contig_index, variant_keys, graph));
            }
        }
    }

    for (Contig& contig : contig_list) {
        std::sort(contig.variants.begin(), contig.variants.end(), [&](size_t a, size_t b) {
            const NestedVariant& x = variant_list[a];
            const NestedVariant& y = variant_list[b];
            return std::tie(x.pos, x.ref, x.alt) < std::tie(y.pos, y.ref, y.alt);
        });
    }
    sort_contigs();
}

size_t NestedVariantCatalog::reader_of(size_t allele, size_t contig) {
    std::vector<std::pair<size_t, size_t>>& readings = allele_readings[allele];
    if (auto same = std::find_if(readings.begin(), readings.end(), [&](const std::pair<size_t, size_t>& reading) {
            return reading.first == contig;
        });
        same != readings.end()) {
        return same->second;
    }
    std::vector<Reader>& readers = contig_list[contig].readers;
    readings.emplace_back(contig, readers.size());
    readers.push_back({ allele, {} });
    return readings.back().second;
}

size_t NestedVariantCatalog::list_variant(const NestedVariant& variant, size_t contig, VariantKeys& keys,
                                          const GBWTGraph& graph) {
    auto [entry, is_new] = keys.try_emplace(std::make_tuple(variant.contig, variant.pos, variant.ref, variant.alt),
                                            variant_list.size());
    if (is_new) {
        variant_list.push_back(variant);
        ids.push_back(variant.id(contig_list[contig].name, graph));
        contig_list[contig].variants.push_back(entry->second);
    }
    return entry->second;
}

void NestedVariantCatalog::sort_contigs() {
    std::vector<size_t> order(contig_list.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return contig_list[a].last_record < contig_list[b].last_record;
    });
    std::vector<size_t> rank(order.size());
    std::vector<Contig> sorted;
    sorted.reserve(order.size());
    for (size_t i = 0; i < order.size(); i++) {
        rank[order[i]] = i;
        sorted.push_back(std::move(contig_list[order[i]]));
    }
    contig_list = std::move(sorted);
    for (std::vector<std::pair<size_t, size_t>>& readings : allele_readings) {
        for (auto& [contig, reader] : readings) {
            contig = rank[contig];
        }
    }
}

//------------------------------------------------------------------------------

NestedGenotyper::NestedGenotyper(const NestedVariantCatalog& catalog, const VcfCalls& calls,
                                 std::vector<std::string> extra_info,
                                 const std::unordered_map<std::string, std::string>& top_level_snarls,
                                 const GBWTGraph& graph, std::ostream& out) :
    catalog(catalog), calls(calls), extra_info(std::move(extra_info)), top_level_snarls(top_level_snarls), graph(graph),
    out(out),
    carried(catalog.contigs().size()), ploidy(calls.samples.size(), 2) {
}

std::vector<std::string> NestedGenotyper::header_lines() const {
    std::vector<std::string> lines {
        "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">",
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
        "##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in "
        "the graph with conflicting alleles\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality: phred-scaled probability that the "
        "genotype is wrong. On a record vg lift breaks out of called alleles, the lowest GQ of the calls that give "
        "the genotype.\">"
    };
    std::istringstream nesting(SortedVCFWriter::nesting_info_headers());
    for (std::string line; std::getline(nesting, line);) {
        lines.push_back(line);
    }
    for (const NestedVariantCatalog::Contig& contig : catalog.contigs()) {
        if (!contig.variants.empty()) {
            lines.push_back("##contig=<ID=" + contig.name + ",length=" + std::to_string(contig.length) + ">");
        }
    }
    return lines;
}

namespace {

/// Return the number in a GQ field, or -1 if it holds none.
int parse_quality(std::string_view text) {
    int parsed = 0;
    auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
    return result.ec == std::errc() && result.ptr == text.data() + text.size() ? parsed : -1;
}

/// Return the number of the ALT allele that an allele of a GT field calls, or
/// 0 if it calls REF or is missing. A number above the record's ALT count is
/// returned as it stands, for the caller to count and read as REF.
size_t alt_number(std::string_view allele) {
    size_t number = 0;
    auto parsed = std::from_chars(allele.data(), allele.data() + allele.size(), number);
    return parsed.ec == std::errc() ? number : 0;
}

/**
 * A table of the allele that each reader of a contig has at each variant on
 * the contig.
 *
 * A reader has 1, for ALT, at a variant that one of its passes has. It has 0,
 * for REF, at a variant whose part of the contig a pass reads without the
 * variant. It has -1, for no allele, at the rest.
 */
class ReaderAlleleTable {
public:
    ReaderAlleleTable(const NestedVariantCatalog::Contig& contig, const std::vector<NestedVariant>& variants) {
        // The index of each variant on the contig, by its index in the catalog.
        std::unordered_map<size_t, size_t> local;
        for (size_t i = 0; i < contig.variants.size(); i++) {
            local[contig.variants[i]] = i;
        }
        for (const NestedVariantCatalog::Reader& reader : contig.readers) {
            rows.push_back(reader_alleles(reader, contig, variants, local));
        }
    }

    /// Return the allele that the reader has at the contig's ith variant.
    int8_t allele(size_t reader, size_t i) const {
        return rows[reader][i];
    }

private:
    /// Return the alleles of the reader, in the order of the contig's variants.
    static std::vector<int8_t> reader_alleles(const NestedVariantCatalog::Reader& reader,
                                              const NestedVariantCatalog::Contig& contig,
                                              const std::vector<NestedVariant>& variants,
                                              const std::unordered_map<size_t, size_t>& local) {
        std::vector<int8_t> row(contig.variants.size(), -1);
        for (const NestedVariantCatalog::Reading& reading : reader.passes) {
            for (size_t i = 0; i < contig.variants.size(); i++) {
                const NestedVariant& variant = variants[contig.variants[i]];
                if (reading.start <= variant.ref_start() && variant.ref_end() <= reading.end) {
                    row[i] = std::max<int8_t>(row[i], 0);
                }
            }
            for (size_t variant : reading.variants) {
                row[local.at(variant)] = 1;
            }
        }
        return row;
    }

    /// The alleles of each reader, in the order of the contig's variants.
    std::vector<std::vector<int8_t>> rows;
};

}

/**
 * A sample's genotype at a variant, from the walks that it carries along the
 * variant's contig.
 *
 * Each walk that has an allele at the variant fills a haplotype with it. The
 * walks of a phased sample fill their own haplotypes, and those of an
 * unphased one fill the haplotypes in order of allele. More walks than
 * haplotypes, or two on one, conflict, and leave every haplotype missing. The
 * genotype's GQ is the lowest GQ of the calls of its walks.
 */
class NestedGenotyper::Genotype {
public:
    /// Make a sample's genotype at the contig's ith variant. The sample carries
    /// the walks along the contig, and has ploidy haplotypes.
    Genotype(const std::vector<Carried>& walks, const ReaderAlleleTable& table, size_t i, size_t ploidy) :
        alleles(std::max<size_t>(ploidy, 1), '.') {
        // The allele of each walk that has one here, with its haplotype.
        std::vector<std::pair<int, int>> values;
        bool phased = true;
        for (const Carried& walk : walks) {
            if (int8_t allele = table.allele(walk.reader, i); allele >= 0) {
                values.emplace_back(walk.haplotype, allele);
                phased = phased && walk.haplotype >= 0;
                lowest_quality = std::min(lowest_quality, walk.quality);
            }
        }
        has_walks = !values.empty();
        conflict = values.size() > alleles.size();
        if (!conflict && phased) {
            place_phased(values);
        } else if (!conflict) {
            place_unphased(values);
        }
        if (conflict) {
            std::fill(alleles.begin(), alleles.end(), '.');
        }
        separator = phased && has_walks && !conflict ? '|' : '/';
    }

    /// Return true if the sample's walks conflict.
    bool conflicts() const {
        return conflict;
    }
    /// Return the number of haplotypes with an allele.
    size_t called() const {
        return alleles.size() - std::count(alleles.begin(), alleles.end(), '.');
    }
    /// Return the number of haplotypes with the ALT allele.
    size_t alt_called() const {
        return std::count(alleles.begin(), alleles.end(), '1');
    }
    /// Return the GQ, or -1 if it is unknown.
    int quality() const {
        return !has_walks || conflict || lowest_quality < 0 ? -1 : lowest_quality;
    }
    /// Append the genotype as GT, such as 0|1.
    void append_gt(std::string& text) const {
        for (size_t h = 0; h < alleles.size(); h++) {
            if (h > 0) {
                text += separator;
            }
            text += alleles[h];
        }
    }

private:
    /// Fill the haplotype of each walk with its allele. Conflict if a walk's
    /// haplotype is past the last one, or already filled.
    void place_phased(const std::vector<std::pair<int, int>>& values) {
        for (const auto& [haplotype, allele] : values) {
            if (size_t(haplotype) >= alleles.size() || alleles[haplotype] != '.') {
                conflict = true;
                return;
            }
            alleles[haplotype] = char('0' + allele);
        }
    }
    /// Fill the first haplotypes with the walks' alleles, REF first.
    void place_unphased(std::vector<std::pair<int, int>>& values) {
        std::sort(values.begin(), values.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.second < b.second;
        });
        for (size_t h = 0; h < values.size(); h++) {
            alleles[h] = char('0' + values[h].second);
        }
    }

    /// The allele of each haplotype, as '0', '1', or '.' if missing.
    std::vector<char> alleles;
    /// The separator of the alleles in GT, which is '|' if phased and '/' if not.
    char separator = '/';
    /// True if some walk has an allele here.
    bool has_walks = false;
    bool conflict = false;
    /// The lowest GQ of the calls of the walks, which is -1 for a call without one.
    int lowest_quality = std::numeric_limits<int>::max();
};

/**
 * The sample columns of a record, with the counts for its INFO.
 *
 * It takes the genotypes of the samples in column order. It keeps the columns
 * with GT alone, and with GT and GQ. A record whose calls have no GQ gets the
 * first.
 */
class NestedGenotyper::SampleColumns {
public:
    /// Add the genotype of the next sample.
    void add(const Genotype& genotype, std::string_view sample) {
        size_t sample_start = genotypes.size();
        genotypes += '\t';
        genotype.append_gt(genotypes);
        allele_count += genotype.called();
        alt_count += genotype.alt_called();
        if (genotype.called() > 0) {
            with_data++;
        }
        if (genotype.conflicts()) {
            conflicted.append(conflicted.empty() ? "" : ",").append(sample);
        }
        genotypes_with_quality.append(genotypes, sample_start).append(":");
        if (genotype.quality() >= 0) {
            genotypes_with_quality += std::to_string(genotype.quality());
            has_quality = true;
        } else {
            genotypes_with_quality += '.';
        }
    }

    /// Write INFO/AC, AN, AF and NS, each after a ";". Leave out AF if no
    /// haplotype has an allele.
    void write_counts(std::ostream& out) const {
        out << ";AC=" << alt_count << ";AN=" << allele_count;
        if (allele_count > 0) {
            out << ";AF=" << static_cast<double>(alt_count) / static_cast<double>(allele_count);
        }
        out << ";NS=" << with_data;
    }
    /// Return the names of the samples whose walks conflict, separated by commas.
    const std::string& conflicted_samples() const {
        return conflicted;
    }
    /// Write FORMAT and the sample columns, each after a tab.
    void write(std::ostream& out) const {
        out << (has_quality ? "\tGT:GQ" : "\tGT") << (has_quality ? genotypes_with_quality : genotypes);
    }

private:
    /// The sample columns with GT alone, and with GT and GQ.
    std::string genotypes;
    std::string genotypes_with_quality;
    /// True if some genotype has a GQ.
    bool has_quality = false;
    /// The number of haplotypes with an allele, and with ALT.
    size_t allele_count = 0;
    size_t alt_count = 0;
    /// The number of samples with an allele.
    size_t with_data = 0;
    std::string conflicted;
};

void NestedGenotyper::add_call(size_t record, const VcfCursor& cursor) {
    const VcfSite& site = calls.sites[record];
    bool reads_contig = false;
    for (size_t i = 0; i < site.alt_count; i++) {
        reads_contig |= !catalog.readings(site.first_allele + i).empty();
    }
    std::vector<std::string_view> genotypes;
    std::vector<std::string_view> qualities;
    if (reads_contig) {
        genotypes = cursor.genotypes();
        qualities = cursor.sample_values("GQ");
    }
    for (size_t sample = 0; sample < genotypes.size() && sample < calls.samples.size(); sample++) {
        int quality = sample < qualities.size() ? parse_quality(qualities[sample]) : -1;
        ploidy[sample] = carry_genotype(site, sample, genotypes[sample], quality);
    }

    const std::vector<NestedVariantCatalog::Contig>& contigs = catalog.contigs();
    while (next_contig < contigs.size() && contigs[next_contig].last_record <= record) {
        write_contig(next_contig++);
    }
}

size_t NestedGenotyper::carry_genotype(const VcfSite& site, size_t sample, std::string_view genotype, int quality) {
    bool phased = genotype.find('/') == std::string_view::npos;
    size_t haplotype = 0;
    size_t start = 0;
    while (start <= genotype.size()) {
        size_t end = std::min(genotype.find_first_of("/|", start), genotype.size());
        // A REF or missing allele adds no walk.
        size_t number = alt_number(genotype.substr(start, end - start));
        if (number > site.alt_count) {
            // The GT names an allele the record does not have, so there is no
            // walk to take. Reading it as REF is a guess, so say how often.
            past_end_count++;
            number = 0;
        }
        if (number > 0) {
            carry(sample, site.first_allele + number - 1, phased ? int(haplotype) : -1, quality);
        }
        start = end + 1;
        haplotype++;
    }
    return haplotype;
}

void NestedGenotyper::carry(size_t sample, size_t allele, int haplotype, int quality) {
    for (const auto& [contig, reader] : catalog.readings(allele)) {
        if (carried[contig].empty()) {
            carried[contig].resize(calls.samples.size());
        }
        carried[contig][sample].push_back({ reader, haplotype, quality });
    }
}

void NestedGenotyper::write_contig(size_t contig_index) {
    const NestedVariantCatalog::Contig& contig = catalog.contigs()[contig_index];
    ReaderAlleleTable table(contig, catalog.variants());
    // The walks that each sample carries along the contig. The list is empty
    // if no call carried a reader of the contig.
    const std::vector<std::vector<Carried>>& walks = carried[contig_index];
    const std::vector<Carried> no_walks;
    for (size_t i = 0; i < contig.variants.size(); i++) {
        SampleColumns columns;
        for (size_t sample = 0; sample < calls.samples.size(); sample++) {
            Genotype genotype(walks.empty() ? no_walks : walks[sample], table, i, ploidy[sample]);
            if (genotype.conflicts()) {
                conflict_count++;
            }
            columns.add(genotype, calls.samples[sample]);
        }
        write_record(contig, contig.variants[i], columns);
    }
    carried[contig_index] = {};
}

void NestedGenotyper::write_record(const NestedVariantCatalog::Contig& contig, size_t variant_index,
                                   const SampleColumns& columns) {
    const NestedVariant& variant = catalog.variants()[variant_index];
    bool has_samples = !calls.samples.empty();
    out << contig.name << '\t' << variant.pos << '\t' << catalog.id(variant_index) << '\t' << variant.ref << '\t'
        << variant.alt << "\t.\t.\tAT=" << walk_text(variant.ref_walk, graph) << ','
        << walk_text(variant.alt_walk, graph);
    if (has_samples) {
        columns.write_counts(out);
    }
    out << ";LV=0;CH=" << contig.level;
    if (auto snarl = top_level_snarls.find(contig.name); snarl != top_level_snarls.end()) {
        out << ";PS=" << snarl->second;
    }
    out << ";RC=" << contig.site_chrom << ";RS=" << contig.site_start << ";RD=" << contig.site_end;
    if (!columns.conflicted_samples().empty()) {
        out << ";CONFLICT=" << columns.conflicted_samples();
    }
    out << extra_info[variant_index];
    if (has_samples) {
        columns.write(out);
    }
    out << '\n';
}

}
