#ifndef INDEX_H
#define INDEX_H

#include <iostream>
#include <exception>
#include <sstream>
#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/options.h"
#include "rocksdb/write_batch.h"
#include "rocksdb/memtablerep.h"
#include "rocksdb/statistics.h"
#include "pb2json.h"
#include "vg.hpp"
#include "hash_map.hpp"

namespace vg {

/*

  Cache our variant graph in a database (rocksdb-backed) which enables us to quickly:
  1) obtain specific nodes and edges from a large graph
  2) search nodes and edges by kmers that they contain or overlap them
  3) use a positional index to quickly build a small portion of the overall

  Each of these functions uses a different subset of the namespace. Our key format is:

  +=\x00 is our 'start' separator
  -=\xff is our 'end' separator --- this makes it easy to do range queries
  ids are stored as raw int64_t

  +m+metadata_key       value // various information about the table
  +g+node_id            node [vg::Node]
  +g+from_id+f+to_id    edge [vg::Edge]
  +g+to_id+t+from_id    null // already stored under from_id+to_id, but this provides reverse index
  +k+kmer+id            position of kmer in node
  +p+position           position overlaps [protobuf]

 */

class Index {

public:

    Index(void);
    Index(string& name);
    ~Index(void);

    void prepare_for_bulk_load(void);
    void open(void);
    void open_read_only(void);
    void open(string& dir);
    void open_read_only(string& dir);
    void reset_options(void);
    void flush(void);
    void compact(void);
    void close(void);

    string name;

    char start_sep;
    char end_sep;

    rocksdb::DB* db;
    rocksdb::Options options;
    rocksdb::WriteOptions write_options;

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

    rocksdb::Status get_node(int64_t id, Node& node);
    rocksdb::Status get_edge(int64_t from, int64_t to, Edge& edge);
    
    const string key_for_node(int64_t id);
    const string key_for_edge_from_to(int64_t from, int64_t to);
    const string key_for_edge_to_from(int64_t to, int64_t from);
    const string key_prefix_for_edges_from_node(int64_t from);
    const string key_prefix_for_edges_to_node(int64_t to);
    const string key_for_kmer(const string& kmer, int64_t id);
    const string key_prefix_for_kmer(const string& kmer);
    const string key_for_metadata(const string& tag);

    void parse_node(const string& key, const string& value, int64_t& id, Node& node);
    void parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge);
    void parse_kmer(const string& key, const string& value, string& kmer, int64_t& id, int32_t& pos);

    void get_context(int64_t id, VG& graph);
    void expand_context(VG& graph, int steps);
    void get_connected_nodes(VG& graph);
    void get_edges_of(int64_t id, vector<Edge>& edges);
    void get_edges_from(int64_t from, vector<Edge>& edges);
    void get_edges_to(int64_t to, vector<Edge>& edges);
    void get_kmer_subgraph(const string& kmer, VG& graph);
    uint64_t approx_size_of_kmer_matches(const string& kmer);
    void for_kmer_range(const string& kmer, function<void(string&, string&)> lambda);
    void get_kmer_positions(const string& kmer, map<int64_t, vector<int32_t> >& positions);
    void get_kmer_positions(const string& kmer, map<int64_t, vector<pair<string, int32_t> > >& positions);

    // for dumping graph state/ inspection
    string entry_to_string(const string& key, const string& value);
    string graph_entry_to_string(const string& key, const string& value);
    string kmer_entry_to_string(const string& key, const string& value);
    string position_entry_to_string(const string& key, const string& value);
    string metadata_entry_to_string(const string& key, const string& value);

    void remember_kmer_size(int size);
    set<int> stored_kmer_sizes(void);
    void store_batch(map<string, string>& items);
    //void store_positions(VG& graph, std::map<long, Node*>& node_path, std::map<long, Edge*>& edge_path);

    // once we have indexed the kmers, we can get the nodes and edges matching
    void kmer_matches(std::string& kmer, std::set<int64_t>& node_ids, std::set<int64_t>& edge_ids);

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
