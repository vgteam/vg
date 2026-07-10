#include "gfaz_to_handle.hpp"

#include "../path.hpp"

#include <GFAz/core/codec/codec.hpp>
#include <GFAz/core/codec/serialization.hpp>
#include <GFAz/compress/io/gfa_write_utils.hpp>

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

// TODO: The GFAz library doesn't really implement a way to access a GFAz file
// without getting all involved in the format. So we implement our own GFAz
// interpretation logic which we keep penned in in dedicated classes.

/**
 * Class which manages decoding GFAz path information.
 *
 * Contains information needed to generate the rGFA visits (which includes node length information and tags), and information used to store and decode the path compression rules.
 *
 * Also used to hold that information during the main decode process.
 * TODO: Split that responsibility out?
 */
struct GFAzPathDecoder {
  /**
   * Holds the length of each graph node, by 1-based ID.
   * TODO: This wastes space at index 0.
   */
  vector<size_t> node_lengths;
  /**
   * Get the number of nodes in the graph.
   */
  inline size_t node_count() const {
    return node_lengths.size() - 1;
  }
  
  vector<OptionalFieldColumn> segment_optional_fields;

  std::pair<vector<int32_t>, vector<int32_t>> rules;
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

  /**
   * Expand a rule for encoded visited node IDs.
   */
  void expand_rule(uint32_t rule_id, bool reverse,
                   uint32_t max_id, vector<NodeId> &out) const;

  /**
   * Decode a sequence of visited node IDs.
   */
  vector<NodeId> decode_sequence_at_index(size_t index, int delta_round) const;
  
  /**
   * Call the callback with rGFA path visits, in order along each path, for
   * rGFA paths a the given rank and below.
   */
  void for_each_rgfa_visit(
    int64_t max_rgfa_rank,
    const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback
  ) const;
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

void GFAzPathDecoder::expand_rule(uint32_t rule_id, bool reverse,
                                     uint32_t max_rule_id, vector<NodeId> &out) const {
  const uint32_t idx = rule_id - min_rule_id;
  const int32_t a = rules.first[idx];
  const int32_t b = rules.second[idx];

  if (!reverse) {
    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_rule_id && abs_a < max_rule_id) {
      expand_rule(abs_a, a < 0, max_rule_id, out);
    } else {
      out.push_back(a);
    }

    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_rule_id && abs_b < max_rule_id) {
      expand_rule(abs_b, b < 0, max_rule_id, out);
    } else {
      out.push_back(b);
    }
  } else {
    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_rule_id && abs_b < max_rule_id) {
      expand_rule(abs_b, b >= 0, max_rule_id, out);
    } else {
      out.push_back(-b);
    }

    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_rule_id && abs_a < max_rule_id) {
      expand_rule(abs_a, a >= 0, max_rule_id, out);
    } else {
      out.push_back(-a);
    }
  }
}

vector<NodeId> GFAzPathDecoder::decode_sequence_at_index(
    size_t index,
    int delta_round
) const {
  if (index + 1 >= path_offsets.size()) {
    throw out_of_range("GFAZ traversal index out of range");
  }

  const size_t start = path_offsets[index];
  const size_t end = path_offsets[index + 1];
  if (end > paths_flat.size()) {
    throw runtime_error("GFAZ traversal block is truncated");
  }

  const uint32_t max_rule_id =
      min_rule_id + static_cast<uint32_t>(rules.first.size());
  const size_t original_length =
      (index + 1 < original_path_offsets.size())
          ? (original_path_offsets[index + 1] - original_path_offsets[index])
          : (end - start);

  vector<NodeId> decoded;
  decoded.reserve(original_length);
  for (size_t pos = start; pos < end; ++pos) {
    const NodeId node = paths_flat[pos];
    const uint32_t abs_id = static_cast<uint32_t>(std::abs(node));
    if (abs_id >= min_rule_id && abs_id < max_rule_id) {
      expand_rule(abs_id, node < 0, max_rule_id, decoded);
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

/**
 * Decode the optional fields AKA tags for a sequence line/row.
 *
 * string_offsets and byte_offsets are decoder state that need to have been
 * updated by decoding all prior rows. They should have one field per entry in
 * columns and should start at 0.
 *
 * Must be called in order; does not provide random access.
 *
 * TODO: Make this a class or a function with a callback instead of something
 * that needs to be called exactly right with the right scratch state.
 *
 * This is a static helper because if we put it on GFAzParser we'd need to
 * include more GFAz headers.
 */
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

/// Translate a vector of GFAz node IDs into a function that produces vg-style
/// rank, id, orientation visit calls. This is a static helper because if we
/// put it on GFAzParser we'd need to include more GFAz headers.
static GFAParser::visit_source_t make_gfaz_visit_source(const vector<NodeId>& visits) {
  return [visits](const GFAParser::visit_iteratee_t& visit_step) {
    if (visits.empty()) {
      visit_step(-1, 0, false);
      return;
    }
    for (size_t i = 0; i < visits.size(); ++i) {
      NodeId visit = visits[i];
      bool is_reverse = visit < 0;
      nid_t node_id = is_reverse ? -static_cast<nid_t>(visit) : static_cast<nid_t>(visit);
      if (!visit_step(i, node_id, is_reverse)) {
        return;
      }
    }
  };
}

void GFAzParser::handle_duplicate_paths(char line_type, const std::function<void(void)>& callback) {
  try {
    callback();
  } catch (GFADuplicatePathError& e) {
    if (stop_on_duplicate_paths) {
      throw;
    }
    #pragma omp critical (cerr)
    cerr << "warning:[GFAzParser] Skipping GFA " << line_type << " line: "
         << e.what() << endl;
  }
}

void GFAzPathDecoder::for_each_rgfa_visit(
  int64_t max_rgfa_rank,
  const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback
) const {
  if (max_rgfa_rank < 0) {
    // We don't need to find any rGFA paths actually
    return;
  }

  // Find the tag columns with the rGFA data
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
    // We don't have all the columns we need to have any rGFA paths
    return;
  }

  // Partition the path name string into the individual path names per tagged node
  vector<size_t> sn_offsets(sn_col->string_lengths.size() + 1, 0);
  for (size_t i = 0; i < sn_col->string_lengths.size(); ++i) {
    sn_offsets[i + 1] = sn_offsets[i] + sn_col->string_lengths[i];
  }

  // We need to produce the path visits in sorted order, so we need to store
  // them.

  /// This represents a visit along an rGFA path, so they can be sorted by
  /// offset.
  struct RGFAVisit {
    int64_t offset;
    nid_t node_id;
    size_t length;
  };
  /// This represents an rGFA path.
  ///
  struct RGFAPath {
    /// Rank of the rGFA path, or -1 for unset.
    /// 
    int64_t rank = -1;
    vector<RGFAVisit> visits;
  };
  // This holds all the rGFA paths by name.
  unordered_map<string, RGFAPath> by_path;

  for (nid_t node_id = 1; static_cast<size_t>(node_id) <= node_count(); ++node_id) {
    // Go through all the nodes and make rGFA visits for them.
    // Determine the GFAz node index from the ID.
    // TODO: is there a better way to write this loop?
    size_t idx = node_id - 1;
    
    if (idx >= sn_col->string_lengths.size() ||
        idx >= so_col->int_values.size() ||
        idx >= sr_col->int_values.size()) {
      // This node is past the range of nodes that can have all the values we
      // need.
      continue;
    }
    // Check the node's tag values.
    uint32_t sn_len = sn_col->string_lengths[idx];
    int64_t so = so_col->int_values[idx];
    int64_t sr = sr_col->int_values[idx];
    if (
      sn_len == 0 ||
      so == numeric_limits<int64_t>::min() ||
      sr == numeric_limits<int64_t>::min() ||
      sr > max_rgfa_rank
    ) {
      // This node doesn't have all the tags we need.
      // TODO: If it has some but not all, should we error?
      continue;
    }
    
    // Unpack the actual rGFA path name
    string path_name = sn_col->concatenated_strings.substr(sn_offsets[idx], sn_len);
    
    if (sr < 0) {
      // Real ranks can't be negative
      throw GFAFormatError("rGFA path " + path_name +
                           " has negative rank " + to_string(sr));
    }

    auto& path_info = by_path[path_name];
    if (path_info.rank != -1 && path_info.rank != sr) {
      throw GFAFormatError("rGFA path " + path_name +
                           " has conflicting ranks " + to_string(sr) + " and " +
                           to_string(path_info.rank));
    }
    path_info.rank = sr;
    // Get the node's length.
    // We know we loaded all the node lengths, so if it's not there, something is very wrong.
    size_t length = node_lengths.at(node_id);
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
  if (filename == "-") {
    // TODO: Handle parsing from standard input when the library can.
    throw invalid_argument("GFAZ input from stdin (-) is not supported");
  }

  // Load the compressed data.
  CompressedData compressed = deserialize_compressed_data(filename);

  if (!compressed.header_line.empty() && compressed.header_line[0] == 'H') {
    // We got a header, so parse it.
    auto parsed = GFAParser::parse_h(compressed.header_line);
    for (auto& listener : this->header_listeners) {
      // And show the listeners.
      listener(get<0>(parsed));
    }
  }
  // TODO: What if there's no header? Shouldn't that be an error?
 
  // We will populate this with information needed to decode paths.
  // Some of that information we also use more generally.
  GFAzPathDecoder gfaz_paths;

  {
    // Do the nodes

    // Temporarily unpack the sequence data
    string segment_sequences =
        Codec::zstd_decompress_string(compressed.segment_sequences_zstd);
    vector<uint32_t> segment_lengths =
        Codec::zstd_decompress_uint32_vector(compressed.segment_seq_lengths_zstd);
    
    // Reserve space for all the nodes' lengths, plus an empty slot at ID 0.
    gfaz_paths.node_lengths.resize(segment_lengths.size() + 1, 0);
  
    // Configure the ID map to pass through IDs.
    auto& map_info = this->id_map();
    map_info.numeric_mode = true;
    map_info.direct_numeric_lookup = true;
    map_info.max_id = gfaz_paths.node_count();
    // We know the ID map is empty when parsing starts so we don't have to clear it.

    // Unpack all the optional fields (node tags)
    gfaz_paths.segment_optional_fields.reserve(compressed.segment_optional_fields_zstd.size());
    for (const auto& column : compressed.segment_optional_fields_zstd) {
      gfaz_paths.segment_optional_fields.push_back(decompress_optional_column(column));
    }
  
    // Set up decoder scratch for the tags
    vector<size_t> optional_string_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    vector<size_t> optional_byte_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    // And for the segment sequences
    size_t segment_seq_offset = 0;
    for (size_t i = 0; i < segment_lengths.size(); ++i) {
      // Get the ID of this segment
      nid_t id = static_cast<nid_t>(i + 1);
      
      // Load its sequence
      uint32_t len = segment_lengths[i];
      if (segment_seq_offset + len > segment_sequences.size()) {
        throw runtime_error("GFAZ segment sequence column is truncated");
      }
      string sequence = segment_sequences.substr(segment_seq_offset, len);

      // Generate a string representation for its tags
      auto tags = collect_optional_tags_for_row(gfaz_paths.segment_optional_fields,
                                                optional_string_offsets,
                                                optional_byte_offsets,
                                                i);

      for (auto& listener : this->node_listeners) {
        // Show it to the listeners
        listener(id, as_chars(sequence), tags);
      }

      // Remember its length because we need node lengths later when doing rGFA visits.
      // TODO: Why not just generarte and hold onto rGFA visits now???
      gfaz_paths.node_lengths[id] = len;

      // The next node's sequence is after this one's.
      segment_seq_offset += len;
    }
  }
  

  {
    // Do the edges

    LinkData links;

    // Temporarily unpack the link information
    links.from_ids = Codec::decompress_delta_varint_uint32(
        compressed.link_from_ids_zstd, compressed.num_links);
    links.to_ids = Codec::decompress_delta_varint_uint32(
        compressed.link_to_ids_zstd, compressed.num_links);
    links.from_orients = Codec::decompress_orientations(
        compressed.link_from_orients_zstd, compressed.num_links);
    links.to_orients = Codec::decompress_orientations(
        compressed.link_to_orients_zstd, compressed.num_links);
    links.overlap_nums =
        Codec::zstd_decompress_uint32_vector(compressed.link_overlap_nums_zstd);
    links.overlap_ops =
        Codec::zstd_decompress_char_vector(compressed.link_overlap_ops_zstd);

    for (size_t i = 0; i < links.from_ids.size(); ++i) {
      string overlap;
      if (i < links.overlap_ops.size() && links.overlap_ops[i] != '\0') {
        // Reconstruct the overlap string
        overlap = to_string(i < links.overlap_nums.size() ? links.overlap_nums[i] : 0) +
                  links.overlap_ops[i];
      }
      for (auto& listener : this->edge_listeners) {
        // Show the link to the listener
        listener(links.from_ids[i],
                 i < links.from_orients.size() ? links.from_orients[i] == '-' : false,
                 links.to_ids[i],
                 i < links.to_orients.size() ? links.to_orients[i] == '-' : false,
                 as_chars(overlap),
                 GFAParser::tag_list_t());
      }
    }
  }
  
  // Unpack the rest of the path information
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
  // Including the rule information
  gfaz_paths.rules = decode_rules(compressed);
  gfaz_paths.min_rule_id = compressed.min_rule_id();
  // And the optional flat information
  if (!compressed.paths_zstd.payload.empty()) {
    gfaz_paths.paths_flat = Codec::zstd_decompress_int32_vector(compressed.paths_zstd);
  }
  if (!compressed.walks_zstd.payload.empty()) {
    gfaz_paths.walks_flat = Codec::zstd_decompress_int32_vector(compressed.walks_zstd);
  }

  // Build the offset information
  gfaz_paths.path_offsets = build_offsets(compressed.sequence_lengths);
  gfaz_paths.walk_offsets = build_offsets(compressed.walk_lengths);
  gfaz_paths.original_path_offsets = build_offsets(compressed.original_path_lengths);
  gfaz_paths.original_walk_offsets = build_offsets(compressed.original_walk_lengths);

  // Figure out how many paths there are supposed to be.
  size_t encoded_path_count =
      gfaz_paths.path_offsets.empty() ? 0 : gfaz_paths.path_offsets.size() - 1;
  // TODO: Shouldn't it be an error if we have more or fewer paths than there are supposed to be?
  size_t path_count = min(gfaz_paths.path_names.size(), encoded_path_count);
  for (size_t i = 0; i < path_count; ++i) {
    // Decode each path
    vector<NodeId> visits = gfaz_paths.decode_sequence_at_index(i, compressed.delta_round);
    // Make a visit source to translate it into vg terms
    auto visit_source = make_gfaz_visit_source(visits);
    handle_duplicate_paths('P', [&]() {
      for (auto& listener : this->path_listeners) {
        // And show it to all the listeners
        listener(gfaz_paths.path_names[i], visit_source, GFAParser::tag_list_t());
      }
    });
  }
  
  // Figure out how many walks there are
  size_t walk_count = gfaz_paths.walk_offsets.empty() ? 0 : gfaz_paths.walk_offsets.size() - 1;
  // TODO: We handle having fewer paths than there should be, but what about walks?
  for (size_t i = 0; i < walk_count; ++i) {
    // Decode each walk
    vector<NodeId> visits = gfaz_paths.decode_sequence_at_index(i, compressed.delta_round);
    // Make a visit source to translate it into vg terms
    auto visit_source = make_gfaz_visit_source(visits);

    // Decode the metadata
    string sample_name = i < gfaz_paths.walk_sample_ids.size() ? gfaz_paths.walk_sample_ids[i] : "*";
    int64_t haplotype = i < gfaz_paths.walk_hap_indices.size() ? gfaz_paths.walk_hap_indices[i] : 0;
    string contig_name = i < gfaz_paths.walk_seq_ids.size() ? gfaz_paths.walk_seq_ids[i] : PathMetadata::NO_LOCUS_NAME;
    subrange_t subrange = PathMetadata::NO_SUBRANGE;
    if (i < gfaz_paths.walk_seq_starts.size()) {
      subrange.first = gfaz_paths.walk_seq_starts[i];
      if (i < gfaz_paths.walk_seq_ends.size()) {
        subrange.second = gfaz_paths.walk_seq_ends[i];
      }
    }

    handle_duplicate_paths('W', [&]() {
      for (auto& listener : this->walk_listeners) {
        listener(sample_name, haplotype, contig_name, subrange, visit_source, GFAParser::tag_list_t());
      }
    });
  }

  // Now do the rGFA paths.
  gfaz_paths.for_each_rgfa_visit(
    this->max_rgfa_rank,
    [&](nid_t id, int64_t offset, size_t length, const string& path_name, int64_t path_rank) {
      // When we get an rGFA path visit
      handle_duplicate_paths('S', [&]() {
        for (auto& listener : this->rgfa_listeners) {
          // Tell all the listeners about it
          listener(id, offset, length, path_name, path_rank);
        }
      });
    }
  );
}


/// Helper function to make a parser of the right type from a filename.
static unique_ptr<GFAParser> make_gfa_family_parser_for_file(const string& filename) {
    if (filename != "-" && GFAzParser::looks_like_gfaz(filename)) {
        return make_unique<GFAzParser>();
    }
    return make_unique<GFATextParser>();
}

void load_gfa_or_gfaz_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                                      GFAIDMapInfo* translation) {
    auto parser = make_gfa_family_parser_for_file(filename);
    attach_parser(*parser, translation);
    attach_parser(*parser, graph);
    parser->parse(filename);
}

void load_gfa_or_gfaz_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                                      const string& translation_filename) {

    GFAIDMapInfo id_map_info;
    load_gfa_or_gfaz_to_handle_graph(filename, graph, &id_map_info);
    id_map_info.write_gfa_translation(translation_filename);
}

void load_gfa_or_gfaz_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                                           GFAIDMapInfo* translation, int64_t max_rgfa_rank,
                                           unordered_set<PathSense>* ignore_sense) {
    auto parser = make_gfa_family_parser_for_file(filename);
    attach_parser(*parser, translation);
    attach_parser(*parser, graph, max_rgfa_rank, ignore_sense);
    parser->parse(filename);
}

void load_gfa_or_gfaz_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                                           int64_t max_rgfa_rank, const string& translation_filename,
                                           unordered_set<PathSense>* ignore_sense) {

    GFAIDMapInfo id_map_info;
    load_gfa_or_gfaz_to_path_handle_graph(filename, graph, &id_map_info, max_rgfa_rank, ignore_sense);
    id_map_info.write_gfa_translation(translation_filename);
}

} // namespace algorithms
} // namespace vg
