#include "gfaz_to_handle.hpp"

#include "../path.hpp"

#include <GFAz/core/codec/serialization.hpp>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <tuple>

namespace vg {
namespace algorithms {
using gfaz::GFAZ_MAGIC;
using gfaz::NodeId;

bool GFAzParser::has_magic(const char* buffer, size_t length) {
    if (length < sizeof(uint32_t)) {
        return false;
    }
    uint32_t magic = 0;
    memcpy(&magic, buffer, sizeof(uint32_t));
    return magic == GFAZ_MAGIC;
}

bool GFAzParser::looks_like_gfaz(const string& filename) {
    if (filename.empty() || filename == "-") {
        return false;
    }
    ifstream in(filename, ios::binary);
    if (!in) {
        return false;
    }
    char buffer[sizeof(uint32_t)];
    in.read(buffer, sizeof(buffer));
    return has_magic(buffer, static_cast<size_t>(in.gcount()));
}

void GFAzParser::prepare_decode() {
    pending_rgfa_visits.clear();
    rgfa_path_ranks.clear();

    auto& map_info = id_map();
    map_info.numeric_mode = true;
    map_info.direct_numeric_lookup = true;
    map_info.max_id = 0;
}

void GFAzParser::finish_decode() {
    sort(pending_rgfa_visits.begin(), pending_rgfa_visits.end(),
         [](const PendingRGFAVisit& a, const PendingRGFAVisit& b) {
             if (a.path_name != b.path_name) {
                 return a.path_name < b.path_name;
             }
             if (a.offset != b.offset) {
                 return a.offset < b.offset;
             }
             return a.id < b.id;
         });

    for (const auto& visit : pending_rgfa_visits) {
        handle_duplicate_paths('S', [&]() {
            for (auto& listener : rgfa_listeners) {
                listener(visit.id, visit.offset, visit.length,
                         visit.path_name, visit.path_rank);
            }
        });
    }
}

void GFAzParser::on_header(const string& header_line) {
    auto parsed = GFAParser::parse_h(header_line);
    for (auto& listener : header_listeners) {
        listener(get<0>(parsed));
    }
}

void GFAzParser::on_segment(uint32_t id, const string& sequence,
                            const gfaz::GfaTagList& tags) {
    auto& map_info = id_map();
    map_info.max_id = max(map_info.max_id, static_cast<nid_t>(id));

    for (auto& listener : node_listeners) {
        listener(static_cast<nid_t>(id), as_chars(sequence), tags);
    }

    if (max_rgfa_rank < 0 || tags.size() < 3) {
        return;
    }

    string path_name;
    int64_t offset = 0;
    int64_t path_rank = -1;
    if (!decode_rgfa_tags(tags, &path_name, &offset, &path_rank) ||
        path_rank > max_rgfa_rank) {
        return;
    }
    if (path_rank < 0) {
        throw GFAFormatError("rGFA path " + path_name +
                             " has negative rank " + to_string(path_rank));
    }

    auto found = rgfa_path_ranks.find(path_name);
    if (found == rgfa_path_ranks.end()) {
        rgfa_path_ranks.emplace(path_name, path_rank);
    } else if (found->second != path_rank) {
        throw GFAFormatError("rGFA path " + path_name +
                             " has conflicting ranks " +
                             to_string(path_rank) + " and " +
                             to_string(found->second));
    }

    pending_rgfa_visits.push_back({static_cast<nid_t>(id), offset,
                                   sequence.size(), path_name, path_rank});
}

void GFAzParser::on_link(uint32_t from_id, bool from_is_reverse,
                         uint32_t to_id, bool to_is_reverse,
                         const string& overlap,
                         const gfaz::GfaTagList& tags) {
    for (auto& listener : edge_listeners) {
        listener(static_cast<nid_t>(from_id), from_is_reverse,
                 static_cast<nid_t>(to_id), to_is_reverse,
                 as_chars(overlap), tags);
    }
}

GFAParser::visit_source_t
GFAzParser::make_visit_source(const vector<NodeId>& visits) {
    return [visits](const visit_iteratee_t& visit_step) {
        if (visits.empty()) {
            visit_step(-1, 0, false);
            return;
        }
        for (size_t i = 0; i < visits.size(); ++i) {
            const NodeId visit = visits[i];
            const bool is_reverse = visit < 0;
            const nid_t node_id = is_reverse
                                      ? -static_cast<nid_t>(visit)
                                      : static_cast<nid_t>(visit);
            if (!visit_step(i, node_id, is_reverse)) {
                return;
            }
        }
    };
}

void GFAzParser::on_path(const string& name, const vector<NodeId>& visits,
                         const string&, const gfaz::GfaTagList& tags) {
    auto visit_source = make_visit_source(visits);
    handle_duplicate_paths('P', [&]() {
        for (auto& listener : path_listeners) {
            listener(name, visit_source, tags);
        }
    });
}

void GFAzParser::on_walk(const string& sample_name, uint32_t haplotype,
                         const string& contig_name, int64_t sequence_start,
                         int64_t sequence_end, const vector<NodeId>& visits) {
    auto visit_source = make_visit_source(visits);
    static const tag_list_t no_tags;
    subrange_t subrange = PathMetadata::NO_SUBRANGE;
    if (sequence_start >= 0) {
        subrange.first = sequence_start;
        if (sequence_end >= 0) {
            subrange.second = sequence_end;
        }
    }

    handle_duplicate_paths('W', [&]() {
        for (auto& listener : walk_listeners) {
            listener(sample_name, haplotype, contig_name, subrange,
                     visit_source, no_tags);
        }
    });
}

void GFAzParser::handle_duplicate_paths(
    char line_type, const function<void(void)>& callback) {
    try {
        callback();
    } catch (GFADuplicatePathError& e) {
        if (stop_on_duplicate_paths) {
            throw;
        }
#pragma omp critical (cerr)
        cerr << "warning:[GFAzParser] Skipping GFA " << line_type
             << " line: " << e.what() << endl;
    }
}

void GFAzParser::parse(istream& in) {
    prepare_decode();
    gfaz::decode_gfa_records(in, *this);
    finish_decode();
}

void GFAzParser::parse(const string& filename) {
    if (filename == "-") {
        parse(cin);
        return;
    }
    prepare_decode();
    gfaz::decode_gfa_records(filename, *this);
    finish_decode();
}

// Machinery for helping call the parser

/// Helper function to make a parser of the right type from a filename.
static unique_ptr<GFAParser> make_gfa_family_parser_for_file(const string& filename) {
    if (filename != "-" && GFAzParser::looks_like_gfaz(filename)) {
        return make_unique<GFAzParser>();
    }
    return make_unique<GFATextParser>();
}

void load_gfa_or_gfaz_to_handle_graph(const string& filename,
                                      MutableHandleGraph* graph,
                                      GFAIDMapInfo* translation) {
    auto parser = make_gfa_family_parser_for_file(filename);
    attach_parser(*parser, translation);
    attach_parser(*parser, graph);
    parser->parse(filename);
}

void load_gfa_or_gfaz_to_handle_graph(const string& filename,
                                      MutableHandleGraph* graph,
                                      const string& translation_filename) {
    GFAIDMapInfo id_map_info;
    load_gfa_or_gfaz_to_handle_graph(filename, graph, &id_map_info);
    id_map_info.write_gfa_translation(translation_filename);
}

void load_gfa_or_gfaz_to_path_handle_graph(
    const string& filename, MutablePathMutableHandleGraph* graph,
    GFAIDMapInfo* translation, int64_t max_rgfa_rank,
    unordered_set<PathSense>* ignore_sense) {
    auto parser = make_gfa_family_parser_for_file(filename);
    attach_parser(*parser, translation);
    attach_parser(*parser, graph, max_rgfa_rank, ignore_sense);
    parser->parse(filename);
}

void load_gfa_or_gfaz_to_path_handle_graph(
    const string& filename, MutablePathMutableHandleGraph* graph,
    int64_t max_rgfa_rank, const string& translation_filename,
    unordered_set<PathSense>* ignore_sense) {
    GFAIDMapInfo id_map_info;
    load_gfa_or_gfaz_to_path_handle_graph(filename, graph, &id_map_info,
                                          max_rgfa_rank, ignore_sense);
    id_map_info.write_gfa_translation(translation_filename);
}

} // namespace algorithms
} // namespace vg
