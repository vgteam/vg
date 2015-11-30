#ifndef INDEX_H
#define INDEX_H

#include <iostream>
#include <exception>
#include <sstream>
#include <climits>

#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/options.h"
#include "rocksdb/write_batch.h"
#include "rocksdb/memtablerep.h"
#include "rocksdb/statistics.h"
#include "rocksdb/cache.h"
#include "rocksdb/slice_transform.h"
#include "rocksdb/table.h"
#include "rocksdb/filter_policy.h"

#include "json2pb.h"
#include "vg.hpp"
#include "hash_map.hpp"

namespace vg {

#ifdef __APPLE__
#include <machine/endian.h>
#include <libkern/OSByteOrder.h>

#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)

#endif

/*

  Cache our variant graph in a database (rocksdb-backed) which enables us to quickly:
  1) obtain specific nodes and edges from a large graph
  2) search nodes and edges by kmers that they contain or overlap them
  3) index the kmers of the graph
  4) store paths and determine the relative locations of nodes and edges in them

  Each of these functions uses a different subset of the namespace. Our key format is:

  +=\x00 is our 'start' separator
  -=\xff is our 'end' separator --- this makes it easy to do range queries

  ids are stored as raw int64_t

  bools are stored a '0' or '1', not as sizeof(bool) bytes, since sizeof(bool) can vary.

  Note that all the graph keys have a node ID and then a "type" character.

  Also note that "pos" in path-related keys is the base-pair coordinate along the path, not the rank of the node.

  Note that we store the edge data for self loops twice.

  // key                                // value
  --------------------------------------------------------------
  +m+metadata_key                       value // various information about the table
  +g+node_id+n                          node [vg::Node]
  +g+node_id+s+other_id+backward        edge [vg::Edge] if node_id <= other_id, else null. edge is on start
  +g+node_id+e+other_id+backward        edge [vg::Edge] if node_id <= other_id, else null. edge is on end
  +g+node_id+p+path_id+pos+backward     mapping [vg::Mapping]
  +k+kmer+node_id                       position of kmer in node [int32_t]
  +p+path_id+pos+backward+node_id       mapping [vg::Mapping]
  +s+node_id+offset                     mapping [vg::Mapping] // mapping-only "side" against one node
  +a+node_id+offset                     alignment [vg::Alignment]

 */

class Index {

public:

    Index(void);
    Index(string& name);
    ~Index(void);

    rocksdb::Options GetOptions(void);
    void open(const std::string& dir, bool read_only = false);
    void open_read_only(string& dir);
    void open_for_write(string& dir);
    void open_for_bulk_load(string& dir);

    void reset_options(void);
    void flush(void);
    void compact(void);
    void close(void);

    string name;

    char start_sep;
    char end_sep;
    int threads;

    rocksdb::DB* db;
    bool is_open;
    bool use_snappy;
    rocksdb::Options db_options;
    rocksdb::WriteOptions write_options;
    rocksdb::ColumnFamilyOptions column_family_options;
    bool bulk_load;
    bool mem_env;
    size_t block_cache_size;

    void load_graph(VG& graph);
    void dump(std::ostream& out);
    void for_all(std::function<void(string&, string&)> lambda);
    void for_range(string& key_start, string& key_end,
                   std::function<void(string&, string&)> lambda);

    void put_node(const Node* node);
    void put_edge(const Edge* edge);
    void batch_node(const Node* node, rocksdb::WriteBatch& batch);
    void batch_edge(const Edge* edge, rocksdb::WriteBatch& batch);
    // Put a kmer that starts at the given index in the given node in the index.
    // The index only stores the kmers that are on the forward strand at their
    // start positions. The aligner is responsible for searching both strands of
    // any query string.
    void put_kmer(const string& kmer,
                  const int64_t id,
                  const int32_t pos);
    void batch_kmer(const string& kmer,
                    const int64_t id,
                    const int32_t pos,
                    rocksdb::WriteBatch& batch);
    void put_metadata(const string& tag, const string& data);
    void put_node_path(int64_t node_id, int64_t path_id, int64_t path_pos, bool backward, const Mapping& mapping);
    void put_path_position(int64_t path_id, int64_t path_pos, bool backward, int64_t node_id, const Mapping& mapping);
    void put_mapping(const Mapping& mapping);
    void put_alignment(const Alignment& alignment);

    rocksdb::Status get_node(int64_t id, Node& node);
    // Takes the nodes and orientations and gets the Edge object with any associated edge data.
    rocksdb::Status get_edge(int64_t from, bool from_start, int64_t to, bool to_end, Edge& edge);
    rocksdb::Status get_metadata(const string& key, string& data);
    // Gets information about the first time the given node appears in the given
    // path, and returns the number of times it appears.
    int get_node_path(int64_t node_id, int64_t path_id, int64_t& path_pos, bool& backward, Mapping& mapping);
    void get_mappings(int64_t node_id, vector<Mapping>& mappings);
    void get_alignments(int64_t node_id, vector<Alignment>& alignments);

    // obtain the key corresponding to each entity
    const string key_for_node(int64_t id);
    const string key_for_edge_on_start(int64_t node_id, int64_t other, bool backward);
    const string key_for_edge_on_end(int64_t node_id, int64_t other, bool backward);
    const string key_prefix_for_edges_on_node_start(int64_t node);
    const string key_prefix_for_edges_on_node_end(int64_t node);
    const string key_for_kmer(const string& kmer, int64_t id);
    const string key_prefix_for_kmer(const string& kmer);
    const string key_for_metadata(const string& tag);
    const string key_for_path_position(int64_t path_id, int64_t path_pos, bool backward, int64_t node_id);
    const string key_for_node_path_position(int64_t node_id, int64_t path_id, int64_t path_pos, bool backward);
    const string key_prefix_for_node_path(int64_t node_id, int64_t path_id);
    const string key_for_mapping_prefix(int64_t node_id);
    const string key_for_mapping(const Mapping& mapping);
    const string key_for_alignment_prefix(int64_t node_id);
    const string key_for_alignment(const Alignment& alignment);

    // deserialize a key/value pair
    void parse_node(const string& key, const string& value, int64_t& id, Node& node);
    // Parse an edge from any of the three kinds of edge keys. For the key types
    // that don't actually store the Edge object, this really constructs a new
    // Edge which won't have the data payload and which might have from and to
    // swapped, but which is equivalent to the actual edge. Populates id1 and id2 with the from and to nodes, and Edge with the actual edge. Populates type with 's' for on-start keys, 'e' for on-end keys, or 'n' for "normal" two-ID edge keys.
    void parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge);
    // We have an overload that doesn't actually fill in an Edge and just looks at the key.
    void parse_edge(const string& key, char& type, int64_t& node_id, int64_t& other_id, bool& backward);
    void parse_kmer(const string& key, const string& value, string& kmer, int64_t& id, int32_t& pos);
    void parse_node_path(const string& key, const string& value,
                         int64_t& node_id, int64_t& path_id, int64_t& path_pos, bool& backward, Mapping& mapping);
    void parse_path_position(const string& key, const string& value,
                             int64_t& path_id, int64_t& path_pos, bool& backward, int64_t& node_id, Mapping& mapping);
    void parse_mapping(const string& key, const string& value, int64_t& node_id, string& hash, Mapping& mapping);
    void parse_alignment(const string& key, const string& value, int64_t& node_id, string& hash, Alignment& alignment);

    // for dumping graph state/ inspection
    string entry_to_string(const string& key, const string& value);
    string graph_entry_to_string(const string& key, const string& value);
    string kmer_entry_to_string(const string& key, const string& value);
    string position_entry_to_string(const string& key, const string& value);
    string metadata_entry_to_string(const string& key, const string& value);
    string node_path_to_string(const string& key, const string& value);
    string path_position_to_string(const string& key, const string& value);
    string mapping_entry_to_string(const string& key, const string& value);
    string alignment_entry_to_string(const string& key, const string& value);

    // accessors, traversal, context
    void get_context(int64_t id, VG& graph);
    // Augment the given graph with the nodes referenced by orphan edges, and
    // all the edges of those nodes, repeatedly for the given number of steps.
    void expand_context(VG& graph, int steps);
    // Add all the elements in the given range to the given graph, if they aren't in it already.
    void get_range(int64_t from_id, int64_t to_id, VG& graph);
    void for_graph_range(int64_t from_id, int64_t to_id, function<void(string&, string&)> lambda);
    void get_connected_nodes(VG& graph);
    // Get the edges on the end of the given node
    void get_edges_on_end(int64_t node, vector<Edge>& edges);
    // Get the edges on the start of the given node
    void get_edges_on_start(int64_t node, vector<Edge>& edges);
    // Get the IDs and orientations of the nodes to the right of the given oriented node
    void get_nodes_next(int64_t node, bool backward, vector<pair<int64_t, bool>>& destinations);
    // Get the IDs and orientations of the nodes to the left of the given oriented node
    void get_nodes_prev(int64_t node, bool backward, vector<pair<int64_t, bool>>& destinations);
    // Get the specified region in bases (start inclusive, end exclusive) along
    // the named path, and poipulate the given graph with it. Also gets dangling
    // edges not on the path.
    void get_path(VG& graph, const string& name, int64_t start, int64_t end);
    // TODO: unimplemented. Supposed to get the position of a node in a path, or
    // relative to a path (using a BFS to find the nearest previous path node)
    // if not actually in the path.
    void node_path_position(int64_t id, string& path_name, int64_t& position, bool& backward, int64_t& offset);

    // Given a node ID and orientation, and the ID of a path, fill in path_pos
    // with the position along the path of the nearest node left of the node
    // specified that is on that path. Fill in relative_orientation with the
    // orientation of the node on the path relative to the specified orientation
    // of the starting node. Returns the path taken by the breadth-first search,
    // and the ID and orientation of the node on the target path that was
    // reached.
    pair<list<pair<int64_t, bool>>, pair<int64_t, bool>>
    get_nearest_node_prev_path_member(int64_t node_id, bool backward, int64_t path_id,
                                      int64_t& path_pos, bool& relative_orientation,
                                      int max_steps = 4);
    // Given a node ID and orientation, and the ID of a path, fill in path_pos
    // with the position along the path of the nearest node *right* of the node
    // specified that is on that path. Fill in relative_orientation with the
    // orientation of the node on the path relative to the specified orientation
    // of the starting node. Returns the path taken by the breadth-first search,
    // and the ID and orientation of the node on the target path that was
    // reached.
    pair<list<pair<int64_t, bool>>, pair<int64_t, bool>>
    get_nearest_node_next_path_member(int64_t node_id, bool backward, int64_t path_id,
                                      int64_t& path_pos, bool& relative_orientation,
                                      int max_steps = 4);
    // Get the relative position, in both directions, of the given orientation of the given node along the given path.
    bool get_node_path_relative_position(int64_t node_id, bool backward, int64_t path_id,
                                         list<pair<int64_t, bool>>& path_prev, int64_t& prev_pos, bool& prev_orientation,
                                         list<pair<int64_t, bool>>& path_next, int64_t& next_pos, bool& next_orientation);
    // Get a Mapping for this node relative to the given path. The mapping will point to the given node in the given orientation.
    Mapping path_relative_mapping(int64_t node_id, bool backward, int64_t path_id,
                                  list<pair<int64_t, bool>>& path_prev, int64_t& prev_pos, bool& prev_orientation,
                                  list<pair<int64_t, bool>>& path_next, int64_t& next_pos, bool& next_orientation);

    bool surject_alignment(const Alignment& source,
                           set<string>& path_names,
                           Alignment& surjection,
                           string& path_name,
                           int64_t& path_pos,
                           int window = 5);
    // Populates layout with path start and end nodes (and orientations),
    // indexed by path names, and lengths with path lengths indexed by path
    // names.
    void path_layout(map<string, pair<pair<int64_t, bool>, pair<int64_t, bool>> >& layout,
                     map<string, int64_t>& lengths);
    pair<int64_t, bool> path_first_node(int64_t path_id);
    pair<int64_t, bool> path_last_node(int64_t path_id, int64_t& path_length);

    // kmers
    void get_kmer_subgraph(const string& kmer, VG& graph);
    // This is in bytes, and is often 0 for things that occur only once.
    uint64_t approx_size_of_kmer_matches(const string& kmer);
    // This is in bytes, and is often 0 for things that occur only once.
    void approx_sizes_of_kmer_matches(const vector<string>& kmers, vector<uint64_t>& sizes);
    // Run the given function on all the keys and values in the database describing instances of the given kmer.
    void for_kmer_range(const string& kmer, function<void(string&, string&)> lambda);
    // In the given map by node ID, fill in the vector with the offsets in that node at which the given kmer starts.
    void get_kmer_positions(const string& kmer, map<int64_t, vector<int32_t> >& positions);
    // In the given map by kmer, fill in the vector with the node IDs and offsets at which the given kmer starts.
    void get_kmer_positions(const string& kmer, map<string, vector<pair<int64_t, int32_t> > >& positions);
    void prune_kmers(int max_kb_on_disk);

    void remember_kmer_size(int size);
    set<int> stored_kmer_sizes(void);
    void store_batch(map<string, string>& items);
    //void store_positions(VG& graph, std::map<long, Node*>& node_path, std::map<long, Edge*>& edge_path);

    // once we have indexed the kmers, we can get the nodes and edges matching
    void kmer_matches(std::string& kmer, std::set<int64_t>& node_ids, std::set<int64_t>& edge_ids);

    // find lowest key with given kmer string (empty string returned if not found)
    // does not check reverse_complement
    string first_kmer_key(const string& kmer);
    // compare kmers with other index: count the number of unique kmers (taking into account strand)
    // in this index that are found in other, and the number not found.  return <#found, #not found> pair
    pair<int64_t, int64_t> compare_kmers(Index& other);

    // paths
    int64_t get_max_path_id(void);
    void put_max_path_id(int64_t id);
    int64_t new_path_id(const string& name);
    string path_name_prefix(const string& name);
    string path_id_prefix(int64_t id);
    void put_path_id_to_name(int64_t id, const string& name);
    void put_path_name_to_id(int64_t id, const string& name);
    string get_path_name(int64_t id);
    int64_t get_path_id(const string& name);
    void load_paths(VG& graph);
    void store_paths(VG& graph); // of graph
    void store_path(VG& graph, Path& path); // path of graph
    map<string, int64_t> paths_by_id(void);

    // alignments and mappings
    void for_each_mapping(function<void(const Mapping&)> lambda);
    void for_each_alignment(function<void(const Alignment&)> lambda);

    // what table is the key in
    char graph_key_type(const string& key);

};

class indexOpenException: public exception
{
    virtual const char* what() const throw()
    {
        return "unable to open variant graph index";
    }
};

class keyNotFoundException: public exception
{
    virtual const char* what() const throw()
    {
        return "unable to find key in index";
    }
};

}

#endif
