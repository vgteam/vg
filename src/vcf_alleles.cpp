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

#include <sys/stat.h>

namespace vg {

using handlegraph::handle_t;
using gbwtgraph::GBWTGraph;

namespace {

/// Split the text at each separator.
std::vector<std::string_view> split(std::string_view text, char separator) {
    std::vector<std::string_view> parts;
    size_t start = 0;
    size_t end = text.find(separator);
    while (end != std::string_view::npos) {
        parts.push_back(text.substr(start, end - start));
        start = end + 1;
        end = text.find(separator, start);
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
    return is_type && is_walk && is_length ? std::optional<std::string_view>(walk) : std::nullopt;
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
bool holds_allele(std::string_view sequence, std::string_view ref, std::string_view alt) {
    if (contains_bases(sequence, alt)) {
        return true;
    }
    if (alt.size() <= 1) {
        // A one-base ALT against a REF of the same length is a substitution, which
        // the search above has tested against the walk. Against a longer REF it is
        // taken as untestable rather than tested: it is usually VCF's padding base,
        // left in the column by trimming that took the allele's own sequence out,
        // and where the ALT is a placeholder such as * there is no sequence to test
        // at all.
        return ref.size() > alt.size();
    }
    // One end of the ALT is the padding base, which lies outside the walk.
    return contains_bases(sequence, alt.substr(1)) ||
           contains_bases(sequence, alt.substr(0, alt.size() - 1));
}

/// Return the value of the given field, counting from 0, in a sample column
/// whose fields are separated by colons, or nothing if it has fewer fields.
std::string_view sample_field(std::string_view column, size_t field) {
    std::string_view value = column;
    for (size_t skipped = 0; skipped < field; skipped++) {
        size_t colon = value.find(':');
        if (colon == std::string_view::npos) {
            return {};
        }
        value = value.substr(colon + 1);
    }
    return value.substr(0, value.find(':'));
}

}

//------------------------------------------------------------------------------

bool is_variant_file(const std::string& path) {
    // Probe a regular file only. Opening a FIFO or a process substitution reads
    // its first block into htslib's format-detection buffer, and closing the
    // file throws that block away, so whoever opens the path next starts in the
    // middle. A VCF has to be a regular file here in any case, since the lift
    // reads it twice.
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) {
        return false;
    }
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
    std::vector<std::string_view> format_and_samples = split(rest.substr(std::min<size_t>(rest.size(), 1)), '\t');
    if (std::string_view format = format_and_samples.front(); rest.empty() || format.substr(0, format.find(':')) != "GT") {
        return {};
    }
    std::vector<std::string_view> result;
    result.reserve(format_and_samples.size() - 1);
    for (size_t i = 1; i < format_and_samples.size(); i++) {
        result.push_back(format_and_samples[i].substr(0, format_and_samples[i].find(':')));
    }
    return result;
}

std::vector<std::string_view> VcfCursor::sample_values(std::string_view key) const {
    std::string_view rest = after_info();
    if (rest.empty()) {
        return {};
    }
    std::vector<std::string_view> format_and_samples = split(rest.substr(1), '\t');
    std::vector<std::string_view> keys = split(format_and_samples.front(), ':');
    size_t field = std::find(keys.begin(), keys.end(), key) - keys.begin();
    if (field == keys.size()) {
        return {};
    }
    std::vector<std::string_view> result;
    result.reserve(format_and_samples.size() - 1);
    for (size_t i = 1; i < format_and_samples.size(); i++) {
        result.push_back(sample_field(format_and_samples[i], field));
    }
    return result;
}

std::string VcfCursor::locus() const {
    return std::string(column(CHROM)) + ":" + std::string(column(POS));
}

//------------------------------------------------------------------------------

namespace {

/// Throw an error about the record at the cursor.
[[noreturn]] void fail_record(const VcfCursor& cursor, const std::string& problem) {
    throw std::runtime_error("VCF record at " + cursor.locus() + ": " + problem);
}

/**
 * Return the walk of each ALT allele of the record at the cursor, as text.
 *
 * The walks come from INFO/AT, which starts with the REF walk, or else from
 * the ID of a record with one ALT allele, if that is a PanGenie variant ID. A
 * record with neither has none.
 *
 * @throws std::runtime_error if INFO/AT does not have a walk for each allele.
 */
std::vector<std::string_view> alt_walks(const VcfCursor& cursor, size_t alt_count) {
    std::vector<std::string_view> walks;
    if (std::optional<std::string_view> at = cursor.info_value("AT")) {
        walks = split(*at, ',');
        if (walks.size() != alt_count + 1) {
            fail_record(cursor, "INFO/AT has " + std::to_string(walks.size()) + " walks for " +
                                std::to_string(alt_count + 1) + " alleles");
        }
        walks.erase(walks.begin());
    } else if (alt_count == 1) {
        if (std::optional<std::string_view> walk = pangenie_walk(cursor.column(VcfCursor::ID))) {
            walks.push_back(*walk);
        }
    }
    return walks;
}

/**
 * Give an ALT allele the walk in the text, and the interval of its bases if
 * the walk has nodes between its first and last.
 *
 * @param alt The ALT column of the allele.
 * @throws std::runtime_error if the walk is not in the graph, or has fewer
 * than two nodes.
 */
void add_walk(VcfAllele& allele, std::string_view walk, std::string_view alt, const VcfCursor& cursor,
              const GBWTGraph& graph) {
    try {
        allele.walk = parse_walk(walk, graph);
    } catch (const std::runtime_error& e) {
        fail_record(cursor, e.what());
    }
    if (allele.walk.size() < 2) {
        fail_record(cursor, "allele walk " + std::string(walk) + " has fewer than two nodes");
    }
    // The allele is the walk without its first and last nodes, which lie on
    // the reference on either side of the variant.
    if (allele.walk.size() > 2) {
        WalkInterval interval;
        interval.walk.assign(allele.walk.begin() + 1, allele.walk.end() - 1);
        interval.length = interval.walk_length(graph);
        allele.holds_alt = holds_allele(interval.sequence(graph), cursor.column(VcfCursor::REF), alt);
        allele.interval = std::move(interval);
    }
}

}

VcfCalls read_vcf_calls(const std::string& path, const GBWTGraph& graph, int threads) {
    VcfCalls calls;
    std::vector<VcfAllele>& alleles = calls.alleles;
    VcfCursor cursor(path, threads);
    while (cursor.next()) {
        if (cursor.at_header()) {
            // The column names line ends with the sample names, after FORMAT.
            if (std::vector<std::string_view> names = split(cursor.line(), '\t');
                cursor.line().starts_with("#CHROM") && names.size() > 9) {
                calls.samples.assign(names.begin() + 9, names.end());
            }
            continue;
        }
        std::vector<std::string_view> alts = cursor.alt_alleles();
        VcfSite& site = calls.sites.emplace_back();
        site.chrom = cursor.column(VcfCursor::CHROM);
        try {
            site.pos = std::stoull(std::string(cursor.column(VcfCursor::POS)));
        } catch (const std::exception&) {
            fail_record(cursor, "bad POS");
        }
        site.ref_length = cursor.column(VcfCursor::REF).size();
        site.first_allele = alleles.size();
        site.alt_count = alts.size();

        std::vector<std::string_view> walks = alt_walks(cursor, alts.size());
        for (size_t i = 0; i < alts.size(); i++) {
            VcfAllele& allele = alleles.emplace_back();
            allele.name = cursor.column(VcfCursor::ID);
            if (alts.size() > 1) {
                allele.name += ":" + std::to_string(i + 1);
            }
            allele.record = calls.sites.size() - 1;
            if (i < walks.size() && walks[i] != ".") {
                add_walk(allele, walks[i], alts[i], cursor, graph);
            }
        }
    }
    return calls;
}

namespace {

/// Return the start of the declaration on a header line, up to the end of its
/// ID, or nothing if the line declares nothing with an ID.
std::string_view declaration_of(std::string_view line) {
    size_t id = line.find("=<ID=");
    return line.starts_with("##") && id != std::string_view::npos ? line.substr(0, line.find_first_of(",>", id))
                                                                   : std::string_view();
}

/**
 * An editor that adds INFO keys with a value for each ALT allele to the lines
 * of a VCF.
 *
 * In the header, it declares the added keys just before the column names, in
 * place of any old declarations, along with the added header lines that the
 * VCF lacks. In each record, it replaces any old values of the keys with the
 * values of the record's ALT alleles.
 */
class InfoEditor {
public:
    /// Make an editor that adds the INFO keys and the header lines.
    InfoEditor(const std::vector<VcfAlleleInfo>& info, const std::vector<std::string>& header) :
        info(info), header(header) {
    }

    /// Write a header line. Write the added declarations before the column
    /// names line.
    void write_header_line(std::string_view line, std::ostream& out) {
        std::string_view declaration = declaration_of(line);
        if (declaration.starts_with("##INFO=<ID=") && is_added(declaration.substr(11))) {
            return;
        }
        declared.emplace(declaration);
        if (line.starts_with("#CHROM")) {
            for (const std::string& added : header) {
                if (!declared.contains(std::string(declaration_of(added)))) {
                    out << added << '\n';
                }
            }
            for (const VcfAlleleInfo& added : info) {
                out << added.header << '\n';
            }
        }
        out << line << '\n';
    }

    /**
     * Set fields to the INFO column of a record with the added keys.
     *
     * It keeps the old INFO fields other than the added keys. The record's
     * ALT alleles are alleles first to first + count - 1 of the VCF.
     */
    void edit_info(std::string_view old_info, size_t first, size_t count, std::string& fields) const {
        fields.clear();
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
            if (std::all_of(begin, end, [](std::string_view value) { return value == "."; })) {
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
    }

private:
    /// Return true if the editor adds the INFO key.
    bool is_added(std::string_view key) const {
        return std::any_of(info.begin(), info.end(), [&](const VcfAlleleInfo& added) { return added.key == key; });
    }

    const std::vector<VcfAlleleInfo>& info;
    const std::vector<std::string>& header;
    /// The INFO keys, FORMAT keys and contigs that the VCF declares, as the
    /// start of their declarations, up to the ID.
    std::unordered_set<std::string> declared;
};

}

void write_vcf_with_allele_info(const std::string& path, std::ostream& out, const std::vector<VcfAlleleInfo>& info,
                                const std::vector<std::string>& header,
                                const std::function<void(size_t, const VcfCursor&)>& listener, int threads) {
    size_t allele_count = info.empty() ? 0 : info.front().values.size();
    InfoEditor editor(info, header);
    VcfCursor cursor(path, threads);
    // Index of the record's first ALT allele among all the ALT alleles of the VCF.
    size_t first = 0;
    size_t record = 0;
    std::string fields;
    while (cursor.next()) {
        std::string_view line = cursor.line();
        if (cursor.at_header()) {
            editor.write_header_line(line, out);
            continue;
        }
        size_t count = cursor.alt_alleles().size();
        if (first + count > allele_count) {
            throw std::runtime_error("VCF record at " + cursor.locus() + " has more ALT alleles than when it was read");
        }
        editor.edit_info(cursor.column(VcfCursor::INFO), first, count, fields);
        size_t info_start = cursor.column(VcfCursor::INFO).data() - line.data();
        out << line.substr(0, info_start) << fields << cursor.after_info() << '\n';
        first += count;
        if (listener) {
            listener(record, cursor);
        }
        record++;
    }
    if (first != allele_count) {
        throw std::runtime_error("the VCF has fewer ALT alleles than when it was read");
    }
}

}
