#pragma once

/**
 * \file gafkluge.hpp
 *
 * Barebones header-only standalone GAF parser reads GAF lines into and out of structs
 * Named in honour of gfakluge
 */

#include <string>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <map>
#include <functional>

namespace gafkluge {


/**
 * We allow pretty much any field to be set as "*" in the GAF, which gets mapped to a -1
 * in the numeric fields (as there are no valid negative values)
 */
static const std::string missing_string = "*";
static const int64_t missing_int = -1;
inline std::string int_to_string(int64_t i) {
    return i == missing_int ? missing_string : std::to_string(i);
}
inline int64_t string_to_int(const std::string& s) {
    return s == missing_string ? missing_int : std::stol(s);
}
inline bool is_missing(const std::string& s) {
    return s == missing_string;
}
inline bool is_missing(int64_t i) {
    return i == missing_int;
}

/**
 * One step of a GAF path
 */
struct GafStep {

    std::string name;          // Either a path name or segment/node name (see above)
    bool is_reverse;           // In reverse orientation ('<' in GAF)
    bool is_stable;            // True if it's a stable path name (as opposed to segment/node name)
    bool is_interval;          // True if it's an interval of a stable path (false if it's the whole path)
    int64_t start;             // 0-based start (inclusive). only defined if is_stable and is_interval are true
    int64_t end;               // 0-based end (inclusive). only defined if is_stable and is_interval are true
};

/**
 * One line of GAF as described here: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
 */
struct GafRecord {

    std::string query_name;      // Query sequence name
    int64_t query_length;        // Query sequence length
    int64_t query_start;         // 0-based, closed
    int64_t query_end;           // 0-based, open
    int64_t path_length;
    int64_t path_start;          // Start position on the path (0-based)
    int64_t path_end;            // End position on the path (0-based)
    int64_t matches;             // Number of residue matches
    int64_t block_length;        // Alignment block length
    int32_t mapq;                // Mapping quality (0-255; 255 for missing)
    char strand;                 // strand relative to the path + or -
    std::vector<GafStep> path;   // the path

    // Map a tag name to its type and value
    // ex: "de:f:0.2183" in the GAF would appear as opt_fields["de"] = ("f", "0.2183")
    std::map<std::string, std::pair<std::string, std::string>>  opt_fields;

    // Init everything to missing
    GafRecord() : query_length(missing_int), query_start(missing_int), query_end(missing_int),
                  strand(missing_string[0]), path_length(missing_int), path_start(missing_int),
                  path_end(missing_int), matches(missing_int), block_length(missing_int), mapq(255) {}
};

/**
 * Parse a single GAF record
 */ 
inline void parse_gaf_record(const std::string& gaf_line, GafRecord& gaf_record) {

    std::istringstream in(gaf_line);
    std::string buffer;
    
    int col = 1;
    
    std::function<void(void)> scan_column = [&]() {
        getline(in, buffer, '\t');
        if (!in || buffer.empty()) {
            throw std::runtime_error("Error parsing GAF column " + std::to_string(col));
        }
        ++col;
    };
    
    scan_column();
    gaf_record.query_name = std::move(buffer);

    scan_column();
    gaf_record.query_length = string_to_int(buffer);

    scan_column();
    gaf_record.query_start = string_to_int(buffer);

    scan_column();
    gaf_record.query_end = string_to_int(buffer);

    scan_column();
    if (buffer == "-" || buffer == missing_string || buffer == "+") {
        gaf_record.strand = buffer[0];
    } else {
        throw std::runtime_error("Error parsing GAF strand: " + buffer);
    }

    gaf_record.path.clear();
    scan_column();
    if (buffer[0] == '<' || buffer[0] == '>') {
        // our path is a list of oriented segments or intervales
        size_t pos = 0;
        size_t next;
        do {
            GafStep step;
            pos = buffer.find_first_of("><", pos);
            next = buffer.find_first_of("><", pos + 1);
            std::string step_token = buffer.substr(pos, next);
            size_t colon = step_token.find_first_of(':');
            step.is_reverse = step_token[0] == '<';
            if (colon == std::string::npos) {
                // no colon, we interpret the step as a segID
                step.name = step_token.substr(1);
                step.is_stable = false;
                step.is_interval = false;
            } else {
                // colon, we interpret the step as a stable path interval
                step.name = step_token.substr(1, colon);
                step.is_stable = true;
                step.is_interval = true;
                size_t dash = step_token.find_first_of('-', colon);
                if (dash == std::string::npos) {
                    throw std::runtime_error("Error parsing GAF range of " + step_token);
                }
                step.start = std::stol(step_token.substr(colon + 1, dash));
                step.end = std::stol(step_token.substr(dash + 1));
            }
            gaf_record.path.push_back(step);
            pos = next;
        } while (next != std::string::npos);
    } else if (buffer != "*") {
        // our path is a stable path name
        gaf_record.path.resize(1);
        gaf_record.path[0].name = buffer;
        gaf_record.path[0].is_reverse = false;
        gaf_record.path[0].is_stable = true;
        gaf_record.path[0].is_interval = false;
    }

    scan_column();
    gaf_record.path_length = string_to_int(buffer);
    
    scan_column();
    gaf_record.path_start = string_to_int(buffer);
    
    scan_column();
    gaf_record.path_end = string_to_int(buffer);

    scan_column();
    gaf_record.matches = string_to_int(buffer);
    
    scan_column();
    gaf_record.block_length = string_to_int(buffer);

    scan_column();
    if (buffer == missing_string) {
        gaf_record.mapq = missing_int;
    } else {
        gaf_record.mapq = std::stoi(buffer);
        if (gaf_record.mapq >= 255) {
            gaf_record.mapq = missing_int;
        }
    }

    gaf_record.opt_fields.clear();
    do {
        getline(in, buffer, '\t');
        if (in && !buffer.empty()) {
            size_t col1 = buffer.find_first_of(':');
            size_t col2 = buffer.find_first_of(':', col1 + 1);
            if (buffer.length() < 5 || col1 == std::string::npos || col2 == std::string::npos) {
                throw std::runtime_error("Unable to parse optional tag " + buffer);
            }
            std::string tag = buffer.substr(0, col1);
            std::string type = buffer.substr(col1 + 1, col2 - col1 - 1);
            std::string val = buffer.substr(col2 + 1);
            if (gaf_record.opt_fields.count(tag)) {
                throw std::runtime_error("Duplicate optional field found: " + tag);
            }
            gaf_record.opt_fields[tag] = make_pair(type, val);
        }
    } while (bool(in));

}

/*
 * Visit each CS cigar record as a string.  CS cigars are described here: 
 * https://github.com/lh3/minimap2#the-cs-optional-tag
 */
inline void for_each_cs(const GafRecord& gaf_record, std::function<void(const std::string&)> fn) {
    if (gaf_record.opt_fields.count("cs")) {
        const std::string& cs_cigar = gaf_record.opt_fields.find("cs")->second.second;
        size_t next;
        for (size_t co = 0; co != std::string::npos; co = next) {
            next = cs_cigar.find_first_of(":*-+", co + 1);
            fn(cs_cigar.substr(co, next == std::string::npos ? std::string::npos : next - co));
        }
    }
}
    
/*
 * Write a GAF Step to a stream
 */
inline std::ostream& operator<<(std::ostream& os, const gafkluge::GafStep& gaf_step) {
    if (!gaf_step.is_stable || gaf_step.is_interval) {
        os << (gaf_step.is_reverse ? "<" : ">");
    }
    os << gaf_step.name;
    if (gaf_step.is_interval) {
        os << gaf_step.start << "-" << gaf_step.end;
    }
    return os;
}

/**
 * Write a GAF record to a stream
 */
inline std::ostream& operator<<(std::ostream& os, const gafkluge::GafRecord& gaf_record) {

    os << (gaf_record.query_name.empty() ? gafkluge::missing_string : gaf_record.query_name)  << "\t"
       << gafkluge::int_to_string(gaf_record.query_length) << "\t"
       << gafkluge::int_to_string(gaf_record.query_start) << "\t"
       << gafkluge::int_to_string(gaf_record.query_end) << "\t"
       << gaf_record.strand << "\t";

    if (gaf_record.path.empty()) {
        os << gafkluge::missing_string << "\t"
           << gafkluge::missing_string << "\t"
           << gafkluge::missing_string << "\t"
           << gafkluge::missing_string << "\t"
           << gafkluge::missing_string << "\t"
           << gafkluge::missing_string << "\t";
    } else {
        for (const gafkluge::GafStep& step : gaf_record.path) {
            os << step;
        }
        
        os << "\t"
           << gafkluge::int_to_string(gaf_record.path_length) << "\t"
           << gafkluge::int_to_string(gaf_record.path_start) << "\t"
           << gafkluge::int_to_string(gaf_record.path_end) << "\t"
           << gafkluge::int_to_string(gaf_record.matches) << "\t"
           << gafkluge::int_to_string(gaf_record.block_length) << "\t";
    }
    
    os << (gaf_record.mapq == gafkluge::missing_int ? 255 : gaf_record.mapq);

    for (const auto& it : gaf_record.opt_fields) {
        os << "\t" << it.first << ":" << it.second.first << ":" << it.second.second;
    }

    return os;
}

} // namesapce gafkluge

