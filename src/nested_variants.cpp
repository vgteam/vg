/**
 * \file nested_variants.cpp
 * Variants on gref alt contigs, see nested_variants.hpp.
 */

#include "nested_variants.hpp"
#include "vcf_output_caller.hpp"

#include <algorithm>
#include <charconv>
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

ContigDecomposer::ContigDecomposer(const GBWTGraph& graph, const GrefIndex& gref) :
    graph(graph), gref(gref), contig_ends(gref.path_count()) {
    // Alt contigs have fragment names, which end in _<N>_alt.
    for (uint32_t path = 0; path < gref.path_count(); path++) {
        uint64_t length = gref.path_length(path);
        if (GrefCover::is_gref_name(gref.path_name(path)) && length > 0) {
            contig_ends[path] = std::make_pair(gref.interval(path, 0, 1).walk.front(),
                                               gref.interval(path, length - 1, length).walk.front());
        }
    }
}

bool ContigDecomposer::is_contig(uint32_t path) const {
    return path < contig_ends.size() && contig_ends[path].has_value();
}

std::vector<ContigPass> ContigDecomposer::decompose(const std::vector<handle_t>& walk) const {
    // Group the walk's visits to each contig into passes. A visit continues
    // the contig's latest pass if it reads the contig the same way, further on.
    std::vector<std::vector<Visit>> groups;
    std::unordered_map<uint32_t, size_t> latest;
    for (size_t step = 0; step < walk.size(); step++) {
        const GrefIndex::Position& place = gref.find(graph.get_id(walk[step]));
        if (place.path == GrefIndex::NO_PATH || !is_contig(place.path)) {
            continue;
        }
        Visit visit { step, place.path, graph.get_is_reverse(walk[step]) != place.is_reverse, place.offset,
                      place.offset + graph.get_length(walk[step]) };
        auto found = latest.find(visit.contig);
        if (found != latest.end()) {
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

ContigPass ContigDecomposer::make_pass(const std::vector<handle_t>& walk, std::vector<Visit> visits) const {
    ContigPass pass;
    pass.contig = visits.front().contig;
    pass.first_step = visits.front().step;
    pass.last_step = visits.back().step;
    const auto& [first_node, last_node] = *contig_ends[pass.contig];
    uint64_t length = gref.path_length(pass.contig);

    // Read the walk in the contig's orientation, with steps numbered that way.
    bool is_reverse = visits.front().is_reverse;
    size_t n = walk.size();
    auto step = [&](size_t i) {
        return is_reverse ? graph.flip(walk[n - 1 - i]) : walk[i];
    };
    if (is_reverse) {
        std::reverse(visits.begin(), visits.end());
        for (Visit& visit : visits) {
            visit.step = n - 1 - visit.step;
        }
    }
    auto on_contig = [&](size_t i) {
        return gref.find(graph.get_id(step(i))).path == pass.contig;
    };

    // The walk's nodes strictly between two steps replace the contig's bases
    // [ref_start, ref_end).
    auto add_variant = [&](uint64_t ref_start, uint64_t ref_end, size_t left, size_t right) {
        NestedVariant variant;
        variant.contig = pass.contig;
        std::string ref_bases;
        std::string alt_bases;
        variant.ref_walk.push_back(step(left));
        if (ref_start < ref_end) {
            WalkInterval replaced = gref.interval(pass.contig, ref_start, ref_end);
            ref_bases = replaced.sequence(graph);
            variant.ref_walk.insert(variant.ref_walk.end(), replaced.walk.begin(), replaced.walk.end());
        }
        variant.ref_walk.push_back(step(right));
        variant.alt_walk.push_back(step(left));
        for (size_t i = left + 1; i < right; i++) {
            variant.alt_walk.push_back(step(i));
            alt_bases += graph.get_sequence(step(i));
        }
        variant.alt_walk.push_back(step(right));
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
        // Trim shared bases off the end, then off the start, keeping one in each.
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
        pass.variants.push_back(std::move(variant));
    };

    // A pass that joins the contig after its start replaces the contig's head
    // if an earlier node precedes the contig's first node.
    const Visit& head = visits.front();
    pass.start = head.start;
    if (head.start > 0) {
        for (size_t i = head.step; i-- > 0 && !on_contig(i);) {
            if (graph.has_edge(step(i), first_node)) {
                add_variant(0, head.start, i, head.step);
                pass.start = 0;
                break;
            }
        }
    }
    // Between two visits, the walk's nodes replace the contig's, unless the
    // two nodes follow each other on both.
    for (size_t k = 0; k + 1 < visits.size(); k++) {
        const Visit& left = visits[k];
        const Visit& right = visits[k + 1];
        if (left.end < right.start || right.step > left.step + 1) {
            add_variant(left.end, right.start, left.step, right.step);
        }
    }
    // A pass that leaves the contig before its end replaces the contig's tail
    // if a later node follows the contig's last node.
    const Visit& tail = visits.back();
    pass.end = tail.end;
    if (tail.end < length) {
        for (size_t i = tail.step + 1; i < n && !on_contig(i); i++) {
            if (graph.has_edge(last_node, step(i))) {
                add_variant(tail.end, length, tail.step, i);
                pass.end = length;
                break;
            }
        }
    }
    return pass;
}

//------------------------------------------------------------------------------

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
    std::map<std::tuple<uint32_t, uint64_t, std::string, std::string>, size_t> variant_of_key;
    for (size_t allele = 0; allele < passes.size(); allele++) {
        const VcfSite& site = calls.sites[calls.alleles[allele].record];
        for (const ContigPass& pass : passes[allele]) {
            // Each pass of the walk along another contig, before and after
            // this one, puts the contig a level deeper.
            size_t level = 1;
            for (const ContigPass& other : passes[allele]) {
                level += (other.contig != pass.contig && other.first_step < pass.first_step &&
                          pass.last_step < other.last_step);
            }

            auto [found, is_new] = contig_of_path.emplace(pass.contig, contig_list.size());
            if (is_new) {
                Contig& contig = contig_list.emplace_back();
                contig.path = pass.contig;
                contig.name = PathMetadata::parse_locus_name(gref.path_name(pass.contig));
                contig.length = gref.path_length(pass.contig);
                contig.level = level;
                contig.site_chrom = site.chrom;
                contig.site_start = site.pos;
                contig.site_end = site.pos + site.ref_length;
                contig.last_record = calls.alleles[allele].record;
            }
            size_t contig_index = found->second;
            Contig& contig = contig_list[contig_index];
            contig.level = std::min(contig.level, level);
            contig.site_start = std::min(contig.site_start, site.pos);
            contig.site_end = std::max(contig.site_end, site.pos + site.ref_length);
            contig.last_record = std::max(contig.last_record, calls.alleles[allele].record);

            // A walk can pass along a contig more than once, as one reader.
            std::vector<std::pair<size_t, size_t>>& readings = allele_readings[allele];
            auto same = std::find_if(readings.begin(), readings.end(), [&](const std::pair<size_t, size_t>& r) {
                return r.first == contig_index;
            });
            size_t reader = same == readings.end() ? contig.readers.size() : same->second;
            if (same == readings.end()) {
                readings.emplace_back(contig_index, reader);
                contig.readers.push_back({ allele, {} });
            }
            Reading& reading = contig.readers[reader].passes.emplace_back();
            reading.start = pass.start;
            reading.end = pass.end;
            for (const NestedVariant& variant : pass.variants) {
                auto key = std::make_tuple(variant.contig, variant.pos, variant.ref, variant.alt);
                auto [entry, is_new_variant] = variant_of_key.emplace(key, variant_list.size());
                if (is_new_variant) {
                    variant_list.push_back(variant);
                    ids.push_back(variant.id(contig.name, graph));
                    contig.variants.push_back(entry->second);
                }
                reading.variants.push_back(entry->second);
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

    // Put the contigs in the order in which their last records come.
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
        for (std::pair<size_t, size_t>& reading : readings) {
            reading.first = rank[reading.first];
        }
    }
}

//------------------------------------------------------------------------------

NestedGenotyper::NestedGenotyper(const NestedVariantCatalog& catalog, const VcfCalls& calls,
                                 std::vector<std::string> info, const GBWTGraph& graph, std::ostream& out) :
    catalog(catalog), calls(calls), info(std::move(info)), graph(graph), out(out),
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
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    };
    std::istringstream nesting(VCFOutputCaller::nesting_info_headers());
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

void NestedGenotyper::add_call(size_t record, const VcfCursor& vcf) {
    const VcfSite& site = calls.sites[record];
    bool reads_contig = false;
    for (size_t i = 0; i < site.alt_count; i++) {
        reads_contig |= !catalog.readings(site.first_allele + i).empty();
    }
    std::vector<std::string_view> genotypes;
    if (reads_contig) {
        genotypes = vcf.genotypes();
    }
    for (size_t sample = 0; sample < genotypes.size() && sample < calls.samples.size(); sample++) {
        // Each called ALT allele adds its walk to the contigs that it reads.
        std::string_view genotype = genotypes[sample];
        bool phased = genotype.find('/') == std::string_view::npos;
        size_t haplotype = 0;
        for (size_t start = 0; start <= genotype.size(); haplotype++) {
            size_t end = std::min(genotype.find_first_of("/|", start), genotype.size());
            std::string_view allele = genotype.substr(start, end - start);
            start = end + 1;
            size_t number = 0;
            auto parsed = std::from_chars(allele.data(), allele.data() + allele.size(), number);
            if (parsed.ec != std::errc() || number == 0 || number > site.alt_count) {
                // A REF or missing allele adds no walk.
                continue;
            }
            for (const auto& [contig, reader] : catalog.readings(site.first_allele + number - 1)) {
                if (carried[contig].empty()) {
                    carried[contig].resize(calls.samples.size());
                }
                carried[contig][sample].push_back({ reader, phased ? int(haplotype) : -1 });
            }
        }
        ploidy[sample] = haplotype;
    }

    const std::vector<NestedVariantCatalog::Contig>& contigs = catalog.contigs();
    while (next_contig < contigs.size() && contigs[next_contig].last_record <= record) {
        write_contig(next_contig++);
    }
}

void NestedGenotyper::write_contig(size_t contig_index) {
    const NestedVariantCatalog::Contig& contig = catalog.contigs()[contig_index];
    const std::vector<NestedVariant>& variants = catalog.variants();
    size_t variant_count = contig.variants.size();

    // What each reader says about each variant: 1 if it has the variant, 0 if
    // it reads the variant's part of the contig without it, -1 otherwise.
    std::unordered_map<size_t, size_t> local;
    for (size_t i = 0; i < variant_count; i++) {
        local[contig.variants[i]] = i;
    }
    std::vector<int8_t> says(contig.readers.size() * variant_count, -1);
    for (size_t reader = 0; reader < contig.readers.size(); reader++) {
        int8_t* row = says.data() + reader * variant_count;
        for (const NestedVariantCatalog::Reading& reading : contig.readers[reader].passes) {
            for (size_t i = 0; i < variant_count; i++) {
                const NestedVariant& variant = variants[contig.variants[i]];
                if (reading.start <= variant.ref_start() && variant.ref_end() <= reading.end) {
                    row[i] = std::max<int8_t>(row[i], 0);
                }
            }
            for (size_t variant : reading.variants) {
                row[local.at(variant)] = 1;
            }
        }
    }

    const std::vector<std::vector<Carried>>& walks = carried[contig_index];
    size_t sample_count = calls.samples.size();
    std::vector<std::pair<int, int>> values;
    std::string genotypes;
    for (size_t i = 0; i < variant_count; i++) {
        size_t variant_index = contig.variants[i];
        const NestedVariant& variant = variants[variant_index];
        size_t allele_count = 0;
        size_t alt_count = 0;
        size_t with_data = 0;
        std::string conflicted;
        genotypes.clear();
        for (size_t sample = 0; sample < sample_count; sample++) {
            // The values of the walks that the sample carries here, by haplotype.
            values.clear();
            bool phased = true;
            if (!walks.empty()) {
                for (const Carried& walk : walks[sample]) {
                    int8_t value = says[walk.reader * variant_count + i];
                    if (value >= 0) {
                        values.emplace_back(walk.haplotype, value);
                        phased &= walk.haplotype >= 0;
                    }
                }
            }
            // The walks of a phased sample fill their own haplotypes, and those
            // of an unphased one fill the haplotypes in order of value. More
            // walks than haplotypes, or two on one, conflict.
            size_t haplotypes = std::max<size_t>(ploidy[sample], 1);
            std::vector<char> slots(haplotypes, '.');
            bool conflict = values.size() > haplotypes;
            if (phased) {
                for (const auto& [haplotype, value] : values) {
                    conflict |= (size_t(haplotype) >= haplotypes || slots[haplotype] != '.');
                    if (!conflict) {
                        slots[haplotype] = char('0' + value);
                    }
                }
            } else if (!conflict) {
                std::sort(values.begin(), values.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                    return a.second < b.second;
                });
                for (size_t h = 0; h < values.size(); h++) {
                    slots[h] = char('0' + values[h].second);
                }
            }
            if (conflict) {
                std::fill(slots.begin(), slots.end(), '.');
                conflict_count++;
                conflicted.append(conflicted.empty() ? "" : ",").append(calls.samples[sample]);
            }
            genotypes += '\t';
            for (size_t h = 0; h < haplotypes; h++) {
                if (h > 0) {
                    genotypes += (phased && !values.empty() && !conflict) ? '|' : '/';
                }
                genotypes += slots[h];
                allele_count += slots[h] != '.';
                alt_count += slots[h] == '1';
            }
            with_data += std::any_of(slots.begin(), slots.end(), [](char c) { return c != '.'; });
        }

        out << contig.name << '\t' << variant.pos << '\t' << catalog.id(variant_index) << '\t' << variant.ref << '\t'
            << variant.alt << "\t.\t.\tAT=" << walk_text(variant.ref_walk, graph) << ','
            << walk_text(variant.alt_walk, graph);
        if (sample_count > 0) {
            out << ";AC=" << alt_count << ";AN=" << allele_count;
            if (allele_count > 0) {
                out << ";AF=" << double(alt_count) / allele_count;
            }
            out << ";NS=" << with_data;
        }
        out << ";LV=0;CH=" << contig.level << ";RC=" << contig.site_chrom << ";RS=" << contig.site_start
            << ";RD=" << contig.site_end;
        if (!conflicted.empty()) {
            out << ";CONFLICT=" << conflicted;
        }
        out << info[variant_index];
        if (sample_count > 0) {
            out << "\tGT" << genotypes;
        }
        out << '\n';
    }
    carried[contig_index] = {};
}

}
