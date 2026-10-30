/**
 * \file vcf_alleles.cpp
 * VCF alleles as graph walks, see vcf_alleles.hpp.
 */

#include "vcf_alleles.hpp"

#include <algorithm>
#include <cctype>
#include <ostream>
#include <stdexcept>
#include <unordered_set>

namespace vg {

using handlegraph::handle_t;
using gbwtgraph::GBWTGraph;

namespace {

/// Split the text at each separator.
std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    size_t start = 0;
    for (size_t end = text.find(separator); end != std::string_view::npos; end = text.find(separator, start)) {
        parts.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    parts.push_back(text.substr(start));
    return parts;
}

/**
 * Return the ALT walk in a PanGenie variant ID, or nothing if the ID is not one.
 *
 * The ID is CHROM-POS-TYPE-WALK-LENGTH, where TYPE is SNV, INS, DEL or COMPLEX.
 */
std::optional<std::string_view> pangenie_walk(std::string_view id) {
    // The walk and the length hold no dashes, so we find the fields from the
    // end. The chromosome name may hold dashes.
    size_t length_dash = id.rfind('-');
    if (length_dash == std::string_view::npos || length_dash == 0) {
        return std::nullopt;
    }
    size_t walk_dash = id.rfind('-', length_dash - 1);
    if (walk_dash == std::string_view::npos || walk_dash == 0) {
        return std::nullopt;
    }
    size_t type_dash = id.rfind('-', walk_dash - 1);
    if (type_dash == std::string_view::npos) {
        return std::nullopt;
    }
    std::string_view type = id.substr(type_dash + 1, walk_dash - type_dash - 1);
    std::string_view walk = id.substr(walk_dash + 1, length_dash - walk_dash - 1);
    std::string_view length = id.substr(length_dash + 1);
    bool is_type = (type == "SNV" || type == "INS" || type == "DEL" || type == "COMPLEX");
    bool is_walk = (!walk.empty() && (walk[0] == '>' || walk[0] == '<') &&
                    walk.find_first_not_of("<>0123456789") == std::string_view::npos);
    bool is_length = (!length.empty() && length.find_first_not_of("0123456789") == std::string_view::npos);
    if (is_type && is_walk && is_length) {
        return walk;
    }
    return std::nullopt;
}

/// Return true if the sequence contains the bases, ignoring case.
bool contains_bases(std::string_view sequence, std::string_view bases) {
    auto same = [](char x, char y) {
        return std::toupper(static_cast<unsigned char>(x)) == std::toupper(static_cast<unsigned char>(y));
    };
    return bases.empty() || std::search(sequence.begin(), sequence.end(), bases.begin(), bases.end(), same) != sequence.end();
}

/**
 * Return true if the sequence of an allele's walk holds its ALT column.
 *
 * VCF puts a padding base before an indel, or after one at the start of a
 * contig, and the walk lacks it. Tools such as PanGenie's panel preparation
 * also trim the bases that ALT shares with REF from the columns, but keep the
 * whole nodes in the walk.
 */
bool holds_allele(std::string_view sequence, std::string_view alt) {
    if (contains_bases(sequence, alt)) {
        return true;
    }
    return !alt.empty() &&
           (contains_bases(sequence, alt.substr(1)) || contains_bases(sequence, alt.substr(0, alt.size() - 1)));
}

}

//------------------------------------------------------------------------------

bool is_variant_file(const std::string& path) {
    htsFile* file = hts_open(path.c_str(), "r");
    if (file == nullptr) {
        return false;
    }
    bool result = (hts_get_format(file)->category == variant_data);
    hts_close(file);
    return result;
}

VcfCursor::VcfCursor(const std::string& path, int threads) {
    file = hts_open(path.c_str(), "r");
    if (file == nullptr) {
        throw std::runtime_error("could not open " + path);
    }
    const htsFormat* format = hts_get_format(file);
    if (format->format == bcf) {
        hts_close(file);
        throw std::runtime_error(path + " is BCF, which vg reads only as VCF");
    }
    if (format->compression == bgzf && threads > 1) {
        hts_set_threads(file, threads);
    }
}

VcfCursor::~VcfCursor() {
    ks_free(&buffer);
    hts_close(file);
}

bool VcfCursor::next() {
    do {
        int status = hts_getline(file, '\n', &buffer);
        if (status == -1) {
            return false;
        }
        if (status < -1) {
            throw std::runtime_error("could not read the VCF");
        }
        if (buffer.l > 0 && buffer.s[buffer.l - 1] == '\r') {
            buffer.s[--buffer.l] = '\0';
        }
    } while (buffer.l == 0);
    if (at_header()) {
        return true;
    }
    std::string_view text = line();
    size_t start = 0;
    for (size_t i = 0; i < FIXED_COLUMNS; i++) {
        size_t end = text.find('\t', start);
        if (end == std::string_view::npos) {
            if (i + 1 < FIXED_COLUMNS) {
                throw std::runtime_error("VCF record with fewer than 8 columns: " + std::string(text.substr(0, 100)));
            }
            end = text.size();
        }
        columns[i] = { start, end };
        start = end + 1;
    }
    return true;
}

bool VcfCursor::at_header() const {
    return buffer.l > 0 && buffer.s[0] == '#';
}

std::string_view VcfCursor::line() const {
    return std::string_view(buffer.s, buffer.l);
}

std::string_view VcfCursor::column(Column column) const {
    return line().substr(columns[column].first, columns[column].second - columns[column].first);
}

std::string_view VcfCursor::after_info() const {
    return line().substr(columns[INFO].second);
}

std::vector<std::string_view> VcfCursor::alt_alleles() const {
    std::string_view alt = column(ALT);
    if (alt == ".") {
        return {};
    }
    return split(alt, ',');
}

std::optional<std::string_view> VcfCursor::info_value(std::string_view key) const {
    std::string_view info = column(INFO);
    if (info == ".") {
        return std::nullopt;
    }
    for (std::string_view entry : split(info, ';')) {
        size_t equals = entry.find('=');
        if (entry.substr(0, equals) == key) {
            return equals == std::string_view::npos ? std::string_view() : entry.substr(equals + 1);
        }
    }
    return std::nullopt;
}

std::vector<std::string_view> VcfCursor::genotypes() const {
    std::string_view rest = after_info();
    // The FORMAT column, and then the sample columns, each after a tab.
    std::vector<std::string_view> columns = split(rest.substr(std::min<size_t>(rest.size(), 1)), '\t');
    if (rest.empty() || columns.front().substr(0, columns.front().find(':')) != "GT") {
        return {};
    }
    std::vector<std::string_view> result;
    result.reserve(columns.size() - 1);
    for (size_t i = 1; i < columns.size(); i++) {
        result.push_back(columns[i].substr(0, columns[i].find(':')));
    }
    return result;
}

std::string VcfCursor::locus() const {
    return std::string(column(CHROM)) + ":" + std::string(column(POS));
}

//------------------------------------------------------------------------------

VcfCalls read_vcf_calls(const std::string& path, const GBWTGraph& graph, int threads) {
    VcfCalls calls;
    std::vector<VcfAllele>& alleles = calls.alleles;
    VcfCursor vcf(path, threads);
    while (vcf.next()) {
        if (vcf.at_header()) {
            // The column names line ends with the sample names, after FORMAT.
            std::vector<std::string_view> names = split(vcf.line(), '\t');
            if (vcf.line().substr(0, 6) == "#CHROM" && names.size() > 9) {
                calls.samples.assign(names.begin() + 9, names.end());
            }
            continue;
        }
        std::vector<std::string_view> alts = vcf.alt_alleles();
        auto fail = [&](const std::string& problem) {
            throw std::runtime_error("VCF record at " + vcf.locus() + ": " + problem);
        };
        VcfSite& site = calls.sites.emplace_back();
        site.chrom = vcf.column(VcfCursor::CHROM);
        try {
            site.pos = std::stoull(std::string(vcf.column(VcfCursor::POS)));
        } catch (const std::exception&) {
            fail("bad POS");
        }
        site.ref_length = vcf.column(VcfCursor::REF).size();
        site.first_allele = alleles.size();
        site.alt_count = alts.size();

        // Find the walk of each ALT allele. INFO/AT starts with the REF walk.
        std::vector<std::string_view> walks;
        if (std::optional<std::string_view> at = vcf.info_value("AT")) {
            walks = split(*at, ',');
            if (walks.size() != alts.size() + 1) {
                fail("INFO/AT has " + std::to_string(walks.size()) + " walks for " + std::to_string(alts.size() + 1) +
                     " alleles");
            }
            walks.erase(walks.begin());
        } else if (alts.size() == 1) {
            if (std::optional<std::string_view> walk = pangenie_walk(vcf.column(VcfCursor::ID))) {
                walks.push_back(*walk);
            }
        }

        for (size_t i = 0; i < alts.size(); i++) {
            VcfAllele& allele = alleles.emplace_back();
            allele.name = vcf.column(VcfCursor::ID);
            if (alts.size() > 1) {
                allele.name += ":" + std::to_string(i + 1);
            }
            allele.record = calls.sites.size() - 1;
            if (i >= walks.size() || walks[i] == ".") {
                continue;
            }
            try {
                allele.walk = parse_walk(walks[i], graph);
            } catch (const std::runtime_error& e) {
                fail(e.what());
            }
            if (allele.walk.size() < 2) {
                fail("allele walk " + std::string(walks[i]) + " has fewer than two nodes");
            }
            // The allele is the walk without its first and last nodes, which
            // lie on the reference on either side of the variant.
            if (allele.walk.size() > 2) {
                WalkInterval interval;
                interval.walk.assign(allele.walk.begin() + 1, allele.walk.end() - 1);
                interval.length = interval.walk_length(graph);
                allele.holds_alt = holds_allele(interval.sequence(graph), alts[i]);
                allele.interval = std::move(interval);
            }
        }
    }
    return calls;
}

void write_vcf_with_allele_info(const std::string& path, std::ostream& out, const std::vector<VcfAlleleInfo>& info,
                                const std::vector<std::string>& header,
                                const std::function<void(size_t, const VcfCursor&)>& listener, int threads) {
    size_t allele_count = info.empty() ? 0 : info.front().values.size();
    auto is_added = [&](std::string_view key) {
        return std::any_of(info.begin(), info.end(), [&](const VcfAlleleInfo& added) { return added.key == key; });
    };
    // The INFO keys, FORMAT keys and contigs that the VCF declares, as the
    // start of their declarations, up to the ID.
    std::unordered_set<std::string> declared;
    auto declares = [](std::string_view line) {
        size_t id = line.find("=<ID=");
        return line.substr(0, 2) == "##" && id != std::string_view::npos ? line.substr(0, line.find_first_of(",>", id))
                                                                        : std::string_view();
    };

    VcfCursor vcf(path, threads);
    // Index of the record's first ALT allele among all the ALT alleles of the VCF.
    size_t first = 0;
    size_t record = 0;
    std::string fields;
    while (vcf.next()) {
        std::string_view line = vcf.line();
        if (vcf.at_header()) {
            // Declare the added keys just before the column names, in place of
            // any old declarations, and then the other added lines that the
            // VCF does not have yet.
            std::string_view declaration = declares(line);
            if (declaration.substr(0, 11) == "##INFO=<ID=" && is_added(declaration.substr(11))) {
                continue;
            }
            declared.emplace(declaration);
            if (line.substr(0, 6) == "#CHROM") {
                for (const std::string& added : header) {
                    if (!declared.count(std::string(declares(added)))) {
                        out << added << '\n';
                    }
                }
                for (const VcfAlleleInfo& added : info) {
                    out << added.header << '\n';
                }
            }
            out << line << '\n';
            continue;
        }

        size_t count = vcf.alt_alleles().size();
        if (first + count > allele_count) {
            throw std::runtime_error("VCF record at " + vcf.locus() + " has more ALT alleles than when it was read");
        }
        // Keep the old INFO fields other than the added keys, and add those.
        fields.clear();
        std::string_view old_info = vcf.column(VcfCursor::INFO);
        if (old_info != ".") {
            for (std::string_view entry : split(old_info, ';')) {
                if (!is_added(entry.substr(0, entry.find('=')))) {
                    fields.append(fields.empty() ? "" : ";").append(entry);
                }
            }
        }
        for (const VcfAlleleInfo& added : info) {
            auto begin = added.values.begin() + first;
            auto end = begin + count;
            if (std::all_of(begin, end, [](const std::string& value) { return value == "."; })) {
                continue;
            }
            fields.append(fields.empty() ? "" : ";").append(added.key).append("=");
            for (auto it = begin; it != end; ++it) {
                fields.append(it == begin ? "" : ",").append(*it);
            }
        }
        if (fields.empty()) {
            fields = ".";
        }
        size_t info_start = vcf.column(VcfCursor::INFO).data() - line.data();
        out << line.substr(0, info_start) << fields << vcf.after_info() << '\n';
        first += count;
        if (listener) {
            listener(record, vcf);
        }
        record++;
    }
    if (first != allele_count) {
        throw std::runtime_error("the VCF has fewer ALT alleles than when it was read");
    }
}

}
