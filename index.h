#ifndef INDEX_H
#define INDEX_H

#include <iostream>
#include <exception>
#include <sstream>
#include "leveldb/db.h"
#include "pb2json.h"
#include "vg.h"
#include "portable_endian.h"

namespace vg {

/*

  Cache our variant graph in a database (leveldb-backed) which enables us to quickly:
  1) obtain specific nodes and edges from a large graph
  2) search nodes and edges by kmers that they contain or overlap them
  3) use a positional index to quickly build a small portion of the overall

  Each of these functions uses a different subset of the namespace. Our key format is:

  +=\xff is our default separtor
  (-=\x00 also has use in some cases?)
  ids are stored as raw int64_t

  +m+metadata_key       value // various information about the table
  +g+node_id            node [protobuf]
  +g+from_id+f+to_id    edge [protobuf]
  +g+to_id+t+from_id    null // already stored under from_id+to_id, but this provides reverse index
  +k+kmer               kmer hits [protobuf]
  +p+position           position overlaps [protobuf]

 */

class Index {

public:

    Index(string& name);
    ~Index(void);

    char start_sep;
    char end_sep;

    leveldb::DB* db;
    leveldb::Options options;

    void load_graph(VariantGraph& graph);
    void dump(std::ostream& out);

    void put_node(const Node& node);
    void put_edge(const Edge& edge);

    leveldb::Status get_node(int64_t id, Node& node);
    leveldb::Status get_edge(int64_t from, int64_t to, Edge& edge);
    
    const string key_for_node(int64_t id);
    const string key_for_edge_from_to(int64_t from, int64_t to);
    const string key_for_edge_to_from(int64_t to, int64_t from);

    void parse_node(const string& key, const string& value, int64_t& id, Node& node);
    void parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge);

    void get_context(int64_t id, VariantGraph& graph);
    void get_edges_of(int64_t id, vector<Edge>& edges);
    void get_edges_from(int64_t from, vector<Edge>& edges);
    void get_edges_to(int64_t to, vector<Edge>& edges);

    // for dumping graph state/ inspection
    string entry_to_string(const string& key, const string& value);
    string graph_entry_to_string(const string& key, const string& value);
    string kmer_entry_to_string(const string& key, const string& value);
    string position_entry_to_string(const string& key, const string& value);
    string metadata_entry_to_string(const string& key, const string& value);

    void index_kmers(VariantGraph& graph, int kmer_size = 15);
    void index_positions(VariantGraph& graph, std::map<long, Node*>& node_path, std::map<long, Edge*>& edge_path);

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
