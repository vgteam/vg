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

#include "pb2json.h"
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

  // key                     // value
  --------------------------------------------------------------
  +m+metadata_key            value // various information about the table
  +g+node_id                 node [vg::Node]
  +g+from_id+f+to_id         edge [vg::Edge]
  +g+to_id+t+from_id         null // already stored under from_id+to_id, but this provides reverse index
  +g+node_id+p+path_id+pos   mapping [vg::Mapping]
  +k+kmer+node_id            position of kmer in node [int32_t]
  +p+path_id+pos+node_id     mapping [vg::Mapping]

 */

class Index {

public:

    Index(void);
    Index(string& name);
    ~Index(void);

    rocksdb::Options GetOptions(void);
    void open(const std::string& dir, bool read_only);
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
    void put_kmer(const string& kmer,
                  const int64_t id,
                  const int32_t pos);
    void batch_kmer(const string& kmer,
                    const int64_t id,
                    const int32_t pos,
                    rocksdb::WriteBatch& batch);
    void put_metadata(const string& tag, const string& data);
    void put_node_path(int64_t node_id, int64_t path_id, int64_t path_pos, const Mapping& mapping);
    void put_path_position(int64_t path_id, int64_t path_pos, int64_t node_id, const Mapping& mapping);

    rocksdb::Status get_node(int64_t id, Node& node);
    rocksdb::Status get_edge(int64_t from, int64_t to, Edge& edge);
    rocksdb::Status get_metadata(const string& key, string& data);
    int get_node_path(int64_t node_id, int64_t path_id, int64_t& path_pos, Mapping& mapping);

    // obtain the key corresponding to each entity
    const string key_for_node(int64_t id);
    const string key_for_edge_from_to(int64_t from, int64_t to);
    const string key_for_edge_to_from(int64_t to, int64_t from);
    const string key_prefix_for_edges_from_node(int64_t from);
    const string key_prefix_for_edges_to_node(int64_t to);
    const string key_for_kmer(const string& kmer, int64_t id);
    const string key_prefix_for_kmer(const string& kmer);
    const string key_for_metadata(const string& tag);
    const string key_for_path_position(int64_t path_id, int64_t path_pos, int64_t node_id);
    const string key_for_node_path_position(int64_t node_id, int64_t path_id, int64_t path_pos);
    const string key_prefix_for_node_path(int64_t node_id, int64_t path_id);

    // deserialize a key/value pair
    void parse_node(const string& key, const string& value, int64_t& id, Node& node);
    void parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge);
    void parse_kmer(const string& key, const string& value, string& kmer, int64_t& id, int32_t& pos);
    void parse_node_path(const string& key, const string& value,
                         int64_t& node_id, int64_t& path_id, int64_t& path_pos, Mapping& mapping);
    void parse_path_position(const string& key, const string& value,
                             int64_t& path_id, int64_t& path_pos, int64_t& node_id, Mapping& mapping);

    // for dumping graph state/ inspection
    string entry_to_string(const string& key, const string& value);
    string graph_entry_to_string(const string& key, const string& value);
    string kmer_entry_to_string(const string& key, const string& value);
    string position_entry_to_string(const string& key, const string& value);
    string metadata_entry_to_string(const string& key, const string& value);
    string node_path_to_string(const string& key, const string& value);
    string path_position_to_string(const string& key, const string& value);

    // accessors, traversal, context
    void get_context(int64_t id, VG& graph);
    void expand_context(VG& graph, int steps);
    void get_range(int64_t from_id, int64_t to_id, VG& graph);
    void for_graph_range(int64_t from_id, int64_t to_id, function<void(string&, string&)> lambda);
    void get_connected_nodes(VG& graph);
    void get_edges_of(int64_t id, vector<Edge>& edges);
    void get_edges_from(int64_t from, vector<Edge>& edges);
    void get_edges_to(int64_t to, vector<Edge>& edges);
    void get_path(VG& graph, const string& name, int64_t start, int64_t end);
    void node_path_position(int64_t id, string& path_name, int64_t& position, int64_t& offset);
    pair<list<int64_t>, int64_t> get_nearest_node_prev_path_member(int64_t node_id, int64_t path_id, int64_t& path_pos, int max_steps = 4);
    pair<list<int64_t>, int64_t> get_nearest_node_next_path_member(int64_t node_id, int64_t path_id, int64_t& path_pos, int max_steps = 4);
    bool get_node_path_relative_position(int64_t node_id, int64_t path_id,
                                         list<int64_t>& path_prev, int64_t& prev_pos,
                                         list<int64_t>& path_next, int64_t& next_pos);
    Mapping path_relative_mapping(int64_t node_id, int64_t path_id,
                                  list<int64_t>& path_prev, int64_t& prev_pos,
                                  list<int64_t>& path_next, int64_t& next_pos);
    bool project_path(const Path& source, string path_name, Alignment& projection, int window = 5);
    map<string, pair<int64_t, int64_t> > path_layout(void);
    int64_t path_first_node(int64_t path_id);
    int64_t path_last_node(int64_t path_id);

    // kmers
    void get_kmer_subgraph(const string& kmer, VG& graph);
    uint64_t approx_size_of_kmer_matches(const string& kmer);
    void approx_sizes_of_kmer_matches(const vector<string>& kmers, vector<uint64_t>& sizes);
    void for_kmer_range(const string& kmer, function<void(string&, string&)> lambda);
    void get_kmer_positions(const string& kmer, map<int64_t, vector<int32_t> >& positions);
    void get_kmer_positions(const string& kmer, map<string, vector<pair<int64_t, int32_t> > >& positions);
    void prune_kmers(int max_kb_on_disk);

    void remember_kmer_size(int size);
    set<int> stored_kmer_sizes(void);
    void store_batch(map<string, string>& items);
    //void store_positions(VG& graph, std::map<long, Node*>& node_path, std::map<long, Edge*>& edge_path);

    // once we have indexed the kmers, we can get the nodes and edges matching
    void kmer_matches(std::string& kmer, std::set<int64_t>& node_ids, std::set<int64_t>& edge_ids);

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

    // what table is the key in
    char graph_key_type(string& key);

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
