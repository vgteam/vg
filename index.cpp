#include "index.h"

namespace vg {

using namespace std;

Index::Index(string& name) {
    options.create_if_missing = true;
    options.error_if_exists = true;
    leveldb::Status status = leveldb::DB::Open(options, name, &db);
    if (!status.ok()) {
        throw indexOpenException();
    }
}

Index::~Index(void) {
    delete db;
}

const string Index::key_for_node(int64_t id) {
    string key = '\xff' + "n" + '\xff';
    key.resize(key.size() + sizeof(int64_t));
    memcpy(&id, key.c_str()+key.size(), sizeof(int64_t));
    return key;
}

const string Index::key_for_edge(int64_t from, int64_t to) {
    string key = '\xff' + "e" + '\xff';
    key.resize(key.size() + 2*sizeof(int64_t) + 1);
    memcpy(&from, key.c_str()+key.size(), sizeof(int64_t));
    key.at(key.size() + sizeof(int64_t)) = '\xff';
    memcpy(&to, key.c_str()+key.size() + sizeof(int64_t) + 1, sizeof(int64_t));
    return key;
}

void Index::put_node(const Node& node) {
    string data;
    node.SerializeToString(&data);
    string key = key_for_node(node.id());
    db->Put(leveldb::WriteOptions(), key, data);
}

void Index::put_edge(const Edge& edge) {
    string data;
    edge.SerializeToString(&data);
    string key = key_for_edge(edge.from(), edge.to());
    db->Put(leveldb::WriteOptions(), key, data);
}

void Index::load_graph(VariantGraph& graph) {
    Graph& g = graph.graph;
    for (int i = 0; i < g.nodes_size(); ++i) {
        put_node(g.nodes(i));
    }
    for (int i = 0; i < g.edges_size(); ++i) {
        put_edge(g.edges(i));
    }
}

void index_kmers(VariantGraph& graph, int kmer_size = 15) {
// 
}

void index_positions(VariantGraph& graph, map<long, Node*>& node_path, map<long, Edge*>& edge_path) {

}

}
