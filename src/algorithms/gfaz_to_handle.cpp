#include "gfa_to_handle.hpp"

#include "../path.hpp"

#include <GFAz/codec/codec.hpp>
#include <GFAz/codec/serialization.hpp>
#include <GFAz/io/gfa_write_utils.hpp>

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <tuple>

namespace vg {
namespace algorithms {
using gfaz::gfa_write_utils::SequenceOffsets;
using gfaz::gfa_write_utils::build_offsets;
using gfaz::gfa_write_utils::decode_rules;
using gfaz::gfa_write_utils::decompress_optional_column;
using gfaz::gfa_write_utils::decompress_string_column;
using gfaz::CompressedData;
using gfaz::LinkData;
using gfaz::NodeId;
using gfaz::OptionalFieldColumn;
using gfaz::GFAZ_MAGIC;
using gfaz::deserialize_compressed_data;
namespace Codec = gfaz::Codec;

struct StreamingGFAZPaths {
  string header_line;
  size_t node_count = 0;
  vector<size_t> node_lengths;
  vector<OptionalFieldColumn> segment_optional_fields;

  LinkData links;

  vector<int32_t> rules_first;
  vector<int32_t> rules_second;
  vector<int32_t> paths_flat;
  vector<int32_t> walks_flat;
  SequenceOffsets path_offsets;
  SequenceOffsets walk_offsets;
  SequenceOffsets original_path_offsets;
  SequenceOffsets original_walk_offsets;

  vector<string> path_names;
  vector<string> walk_sample_ids;
  vector<uint32_t> walk_hap_indices;
  vector<string> walk_seq_ids;
  vector<int64_t> walk_seq_starts;
  vector<int64_t> walk_seq_ends;

  uint32_t min_rule_id = 0;
};

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

static CompressedData load_gfaz_compressed(const string &filename) {
  if (filename == "-") {
    throw invalid_argument("GFAZ input from stdin (-) is not supported");
  }
  return deserialize_compressed_data(filename);
}

static void expand_rule(uint32_t rule_id, bool reverse,
                        const vector<int32_t> &first,
                        const vector<int32_t> &second, uint32_t min_id,
                        uint32_t max_id, vector<NodeId> &out) {
  const uint32_t idx = rule_id - min_id;
  const int32_t a = first[idx];
  const int32_t b = second[idx];

  if (!reverse) {
    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_id && abs_a < max_id) {
      expand_rule(abs_a, a < 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(a);
    }

    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_id && abs_b < max_id) {
      expand_rule(abs_b, b < 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(b);
    }
  } else {
    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_id && abs_b < max_id) {
      expand_rule(abs_b, b >= 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(-b);
    }

    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_id && abs_a < max_id) {
      expand_rule(abs_a, a >= 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(-a);
    }
  }
}

static vector<NodeId> decode_sequence_at_index(
    const vector<int32_t> &flat, const SequenceOffsets &compressed_offsets,
    const SequenceOffsets &original_offsets, size_t index,
    const vector<int32_t> &rules_first, const vector<int32_t> &rules_second,
    uint32_t min_rule_id, int delta_round) {
  if (index + 1 >= compressed_offsets.size()) {
    throw out_of_range("GFAZ traversal index out of range");
  }

  const size_t start = compressed_offsets[index];
  const size_t end = compressed_offsets[index + 1];
  if (end > flat.size()) {
    throw runtime_error("GFAZ traversal block is truncated");
  }

  const uint32_t max_rule_id =
      min_rule_id + static_cast<uint32_t>(rules_first.size());
  const size_t original_length =
      (index + 1 < original_offsets.size())
          ? (original_offsets[index + 1] - original_offsets[index])
          : (end - start);

  vector<NodeId> decoded;
  decoded.reserve(original_length);
  for (size_t pos = start; pos < end; ++pos) {
    const NodeId node = flat[pos];
    const uint32_t abs_id = static_cast<uint32_t>(std::abs(node));
    if (abs_id >= min_rule_id && abs_id < max_rule_id) {
      expand_rule(abs_id, node < 0, rules_first, rules_second, min_rule_id,
                  max_rule_id, decoded);
    } else {
      decoded.push_back(node);
    }
  }

  vector<vector<NodeId>> seqs(1);
  seqs[0] = std::move(decoded);
  for (int round = 0; round < delta_round; ++round) {
    Codec::inverse_delta_transform(seqs);
  }
  return std::move(seqs[0]);
}

static void dispatch_rgfa_visits(
    size_t node_count, const vector<OptionalFieldColumn>& segment_optional_fields,
    const vector<size_t>& node_lengths, int64_t max_rgfa_rank,
    const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback);


static GFAParser::chars_t as_chars(const string& value) {
  return make_pair(value.begin(), value.end());
}

static GFAParser::tag_list_t collect_optional_tags_for_row(
    const vector<OptionalFieldColumn>& columns,
    vector<size_t>& string_offsets,
    vector<size_t>& byte_offsets,
    size_t row_index) {
  GFAParser::tag_list_t tags;
  tags.reserve(columns.size());
  for (size_t i = 0; i < columns.size(); ++i) {
    const auto& col = columns[i];
    switch (col.type) {
      case 'A':
        if (row_index < col.char_values.size()) {
          tags.push_back(col.tag + ":A:" + string(1, col.char_values[row_index]));
        }
        break;
      case 'i':
        if (row_index < col.int_values.size()) {
          tags.push_back(col.tag + ":i:" + to_string(col.int_values[row_index]));
        }
        break;
      case 'f':
        if (row_index < col.float_values.size()) {
          tags.push_back(col.tag + ":f:" + to_string(col.float_values[row_index]));
        }
        break;
      case 'Z':
      case 'J':
      case 'H':
        if (row_index < col.string_lengths.size()) {
          uint32_t len = col.string_lengths[row_index];
          string value = col.concatenated_strings.substr(string_offsets[i], len);
          string_offsets[i] += len;
          tags.push_back(col.tag + ":" + string(1, col.type) + ":" + value);
        }
        break;
      case 'B':
        if (row_index < col.b_lengths.size() && row_index < col.b_subtypes.size()) {
          uint32_t len = col.b_lengths[row_index];
          size_t offset = byte_offsets[i];
          byte_offsets[i] += len;
          string value = col.tag + ":B:" + string(1, col.b_subtypes[row_index]);
          for (uint32_t j = 0; j < len && offset + j < col.b_concat_bytes.size(); ++j) {
            value.push_back(',');
            value += to_string(col.b_concat_bytes[offset + j]);
          }
          tags.push_back(std::move(value));
        }
        break;
      default:
        break;
    }
  }
  return tags;
}

static GFAParser::visit_iteratee_t make_gfaz_visit_iteratee(const vector<NodeId>& visits) {
  return [visits](const GFAParser::visit_step_t& visit_step) {
    if (visits.empty()) {
      visit_step(-1, 0, std::string_view(), false);
      return;
    }
    // GFAz encodes node IDs numerically, so we hand them to the listener
    // directly and skip the listener's name-based lookup.
    for (size_t i = 0; i < visits.size(); ++i) {
      NodeId visit = visits[i];
      bool is_reverse = visit < 0;
      nid_t node_id = is_reverse ? -static_cast<nid_t>(visit) : static_cast<nid_t>(visit);
      if (!visit_step(i, node_id, std::string_view(), is_reverse)) {
        return;
      }
    }
  };
}

template<class Callback>
static void handle_duplicate_paths(GFAParser& parser, char line_type, Callback&& callback) {
  try {
    callback();
  } catch (GFADuplicatePathError& e) {
    if (parser.stop_on_duplicate_paths) {
      throw;
    }
    #pragma omp critical (cerr)
    cerr << "warning:[GFAzParser] Skipping GFA " << line_type << " line: "
         << e.what() << endl;
  }
}

static void dispatch_rgfa_visits(size_t node_count,
                                 const vector<OptionalFieldColumn>& segment_optional_fields,
                                 const vector<size_t>& node_lengths,
                                 int64_t max_rgfa_rank,
                                 const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback) {
  if (max_rgfa_rank < 0) {
    return;
  }

  const OptionalFieldColumn* sn_col = nullptr;
  const OptionalFieldColumn* so_col = nullptr;
  const OptionalFieldColumn* sr_col = nullptr;
  for (const auto& col : segment_optional_fields) {
    if (col.tag == "SN" && col.type == 'Z') {
      sn_col = &col;
    } else if (col.tag == "SO" && col.type == 'i') {
      so_col = &col;
    } else if (col.tag == "SR" && col.type == 'i') {
      sr_col = &col;
    }
  }
  if (!sn_col || !so_col || !sr_col) {
    return;
  }

  vector<size_t> sn_offsets(sn_col->string_lengths.size() + 1, 0);
  for (size_t i = 0; i < sn_col->string_lengths.size(); ++i) {
    sn_offsets[i + 1] = sn_offsets[i] + sn_col->string_lengths[i];
  }

  struct RGFAVisit {
    int64_t offset;
    nid_t node_id;
    size_t length;
  };
  struct RGFAPath {
    int64_t rank = 0;
    bool rank_set = false;
    vector<RGFAVisit> visits;
  };
  unordered_map<string, RGFAPath> by_path;
  const int64_t missing_i64 = numeric_limits<int64_t>::min();

  for (nid_t node_id = 1; static_cast<size_t>(node_id) <= node_count; ++node_id) {
    size_t idx = node_id - 1;
    if (idx >= sn_col->string_lengths.size() ||
        idx >= so_col->int_values.size() ||
        idx >= sr_col->int_values.size()) {
      continue;
    }
    uint32_t sn_len = sn_col->string_lengths[idx];
    int64_t so = so_col->int_values[idx];
    int64_t sr = sr_col->int_values[idx];
    if (sn_len == 0 || so == missing_i64 || sr == missing_i64 || sr > max_rgfa_rank) {
      continue;
    }

    string path_name = sn_col->concatenated_strings.substr(sn_offsets[idx], sn_len);
    auto& path_info = by_path[path_name];
    if (path_info.rank_set && path_info.rank != sr) {
      throw GFAFormatError("rGFA path " + path_name +
                           " has conflicting ranks " + to_string(sr) + " and " +
                           to_string(path_info.rank));
    }
    path_info.rank = sr;
    path_info.rank_set = true;
    size_t length = (static_cast<size_t>(node_id) < node_lengths.size()) ? node_lengths[node_id] : 0;
    path_info.visits.push_back({so, node_id, length});
  }

  for (auto& kv : by_path) {
    auto& path_name = kv.first;
    auto& path_info = kv.second;
    sort(path_info.visits.begin(), path_info.visits.end(), [](const RGFAVisit& a, const RGFAVisit& b) {
      return a.offset < b.offset;
    });
    for (const auto& visit : path_info.visits) {
      callback(visit.node_id, visit.offset, visit.length, path_name, path_info.rank);
    }
  }
}

void GFAzParser::parse(istream& in) {
  (void)in;
  throw invalid_argument("GFAZ input from streams is not supported directly");
}

void GFAzParser::parse(const string& filename) {
  CompressedData compressed = load_gfaz_compressed(filename);
  StreamingGFAZPaths gfaz_paths;

  {
    string segment_sequences =
        Codec::zstd_decompress_string(compressed.segment_sequences_zstd);
    vector<uint32_t> segment_lengths =
        Codec::zstd_decompress_uint32_vector(compressed.segment_seq_lengths_zstd);
    gfaz_paths.header_line = compressed.header_line;
    gfaz_paths.node_count = segment_lengths.size();
    gfaz_paths.node_lengths.resize(gfaz_paths.node_count + 1, 0);

    auto& map_info = this->id_map();
    map_info.numeric_mode = true;
    map_info.direct_numeric_lookup = true;
    map_info.max_id = gfaz_paths.node_count;
    map_info.name_to_id->clear();
    map_info.id_to_name.reset();

    gfaz_paths.segment_optional_fields.reserve(compressed.segment_optional_fields_zstd.size());
    for (const auto& column : compressed.segment_optional_fields_zstd) {
      gfaz_paths.segment_optional_fields.push_back(decompress_optional_column(column));
    }

    if (!gfaz_paths.header_line.empty() && gfaz_paths.header_line[0] == 'H') {
      auto parsed = GFAParser::parse_h(gfaz_paths.header_line);
      for (auto& listener : this->header_listeners) {
        listener(get<0>(parsed));
      }
    }

    vector<size_t> optional_string_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    vector<size_t> optional_byte_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    size_t segment_seq_offset = 0;
    for (size_t i = 0; i < segment_lengths.size(); ++i) {
      nid_t id = static_cast<nid_t>(i + 1);
      uint32_t len = segment_lengths[i];
      if (segment_seq_offset + len > segment_sequences.size()) {
        throw runtime_error("GFAZ segment sequence column is truncated");
      }
      string sequence = segment_sequences.substr(segment_seq_offset, len);
      auto tags = collect_optional_tags_for_row(gfaz_paths.segment_optional_fields,
                                                optional_string_offsets,
                                                optional_byte_offsets,
                                                i);
      for (auto& listener : this->node_listeners) {
        listener(id, as_chars(sequence), tags);
      }
      gfaz_paths.node_lengths[id] = len;
      segment_seq_offset += len;
    }
  }

  gfaz_paths.links.from_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_from_ids_zstd, compressed.num_links);
  gfaz_paths.links.to_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_to_ids_zstd, compressed.num_links);
  gfaz_paths.links.from_orients = Codec::decompress_orientations(
      compressed.link_from_orients_zstd, compressed.num_links);
  gfaz_paths.links.to_orients = Codec::decompress_orientations(
      compressed.link_to_orients_zstd, compressed.num_links);
  gfaz_paths.links.overlap_nums =
      Codec::zstd_decompress_uint32_vector(compressed.link_overlap_nums_zstd);
  gfaz_paths.links.overlap_ops =
      Codec::zstd_decompress_char_vector(compressed.link_overlap_ops_zstd);

  for (size_t i = 0; i < gfaz_paths.links.from_ids.size(); ++i) {
    string overlap;
    if (i < gfaz_paths.links.overlap_ops.size() && gfaz_paths.links.overlap_ops[i] != '\0') {
      overlap = to_string(i < gfaz_paths.links.overlap_nums.size() ? gfaz_paths.links.overlap_nums[i] : 0) +
                gfaz_paths.links.overlap_ops[i];
    }
    for (auto& listener : this->edge_listeners) {
      listener(gfaz_paths.links.from_ids[i],
               i < gfaz_paths.links.from_orients.size() ? gfaz_paths.links.from_orients[i] == '-' : false,
               gfaz_paths.links.to_ids[i],
               i < gfaz_paths.links.to_orients.size() ? gfaz_paths.links.to_orients[i] == '-' : false,
               as_chars(overlap),
               GFAParser::tag_list_t());
    }
  }

  gfaz_paths.path_names =
      decompress_string_column(compressed.names_zstd, compressed.name_lengths_zstd);
  gfaz_paths.walk_sample_ids = decompress_string_column(
      compressed.walk_sample_ids_zstd, compressed.walk_sample_id_lengths_zstd);
  gfaz_paths.walk_hap_indices =
      Codec::zstd_decompress_uint32_vector(compressed.walk_hap_indices_zstd);
  gfaz_paths.walk_seq_ids = decompress_string_column(
      compressed.walk_seq_ids_zstd, compressed.walk_seq_id_lengths_zstd);
  gfaz_paths.walk_seq_starts = Codec::decompress_varint_int64(
      compressed.walk_seq_starts_zstd, compressed.walk_lengths.size());
  gfaz_paths.walk_seq_ends = Codec::decompress_varint_int64(
      compressed.walk_seq_ends_zstd, compressed.walk_lengths.size());

  auto decoded_rules = decode_rules(compressed);
  gfaz_paths.rules_first = std::move(decoded_rules.first);
  gfaz_paths.rules_second = std::move(decoded_rules.second);
  gfaz_paths.min_rule_id = compressed.min_rule_id();

  if (!compressed.paths_zstd.payload.empty()) {
    gfaz_paths.paths_flat = Codec::zstd_decompress_int32_vector(compressed.paths_zstd);
  }
  if (!compressed.walks_zstd.payload.empty()) {
    gfaz_paths.walks_flat = Codec::zstd_decompress_int32_vector(compressed.walks_zstd);
  }
  gfaz_paths.path_offsets = build_offsets(compressed.sequence_lengths);
  gfaz_paths.walk_offsets = build_offsets(compressed.walk_lengths);
  gfaz_paths.original_path_offsets = build_offsets(compressed.original_path_lengths);
  gfaz_paths.original_walk_offsets = build_offsets(compressed.original_walk_lengths);

  size_t encoded_path_count =
      gfaz_paths.path_offsets.empty() ? 0 : gfaz_paths.path_offsets.size() - 1;
  size_t path_count = min(gfaz_paths.path_names.size(), encoded_path_count);
  for (size_t i = 0; i < path_count; ++i) {
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.paths_flat, gfaz_paths.path_offsets,
        gfaz_paths.original_path_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, compressed.delta_round);
    auto visit_iteratee = make_gfaz_visit_iteratee(visits);
    handle_duplicate_paths(*this, 'P', [&]() {
      for (auto& listener : this->path_listeners) {
        listener(gfaz_paths.path_names[i], visit_iteratee, GFAParser::tag_list_t());
      }
    });
  }

  size_t walk_count = gfaz_paths.walk_offsets.empty() ? 0 : gfaz_paths.walk_offsets.size() - 1;
  for (size_t i = 0; i < walk_count; ++i) {
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.walks_flat, gfaz_paths.walk_offsets,
        gfaz_paths.original_walk_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, compressed.delta_round);
    auto visit_iteratee = make_gfaz_visit_iteratee(visits);
    string sample_name = i < gfaz_paths.walk_sample_ids.size() ? gfaz_paths.walk_sample_ids[i] : "*";
    int64_t haplotype = i < gfaz_paths.walk_hap_indices.size() ? gfaz_paths.walk_hap_indices[i] : 0;
    string contig_name = i < gfaz_paths.walk_seq_ids.size() ? gfaz_paths.walk_seq_ids[i] : PathMetadata::NO_LOCUS_NAME;
    int64_t start = i < gfaz_paths.walk_seq_starts.size() ? gfaz_paths.walk_seq_starts[i] : PathMetadata::NO_END_POSITION;
    int64_t end = i < gfaz_paths.walk_seq_ends.size() ? gfaz_paths.walk_seq_ends[i] : PathMetadata::NO_END_POSITION;
    subrange_t subrange = PathMetadata::NO_SUBRANGE;
    if (start != PathMetadata::NO_END_POSITION) {
      subrange.first = start;
      if (end != PathMetadata::NO_END_POSITION) {
        subrange.second = end;
      }
    }
    handle_duplicate_paths(*this, 'W', [&]() {
      for (auto& listener : this->walk_listeners) {
        listener(sample_name, haplotype, contig_name, subrange, visit_iteratee, GFAParser::tag_list_t());
      }
    });
  }

  dispatch_rgfa_visits(gfaz_paths.node_count, gfaz_paths.segment_optional_fields,
                       gfaz_paths.node_lengths, this->max_rgfa_rank,
                       [&](nid_t id, int64_t offset, size_t length, const string& path_name, int64_t path_rank) {
                         handle_duplicate_paths(*this, 'S', [&]() {
                           for (auto& listener : this->rgfa_listeners) {
                             listener(id, offset, length, path_name, path_rank);
                           }
                         });
                       });
}

} // namespace algorithms
} // namespace vg
