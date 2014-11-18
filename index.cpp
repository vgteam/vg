#include "index.h"

namespace vg {

using namespace std;

Index::Index(string& name) {
    options.create_if_missing = true;
    //options.error_if_exists = true;
    leveldb::Status status = leveldb::DB::Open(options, name, &db);
    if (!status.ok()) {
        throw indexOpenException();
    }
}

Index::~Index(void) {
    delete db;
}

const string Index::key_for_node(int64_t id) {
    string key;
    key.resize(2*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    memcpy((void*)(k + sizeof(char)*2), &id, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_from_to(int64_t from, int64_t to) {
    string key;
    key.resize(5*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    memcpy((void*)(k + sizeof(char)*2), &from, sizeof(int64_t));
    k[2 + sizeof(int64_t)] = '\xff';
    k[2 + sizeof(int64_t) + 1] = 'f';
    k[2 + sizeof(int64_t) + 2] = '\xff';
    memcpy((void*)(k + sizeof(char)*2 + sizeof(int64_t) + 3*sizeof(char)), &to, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_to_from(int64_t to, int64_t from) {
    string key;
    key.resize(5*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    memcpy((void*)(k + sizeof(char)*2), &from, sizeof(int64_t));
    k[2 + sizeof(int64_t)] = '\xff';
    k[2 + sizeof(int64_t) + 1] = 't';
    k[2 + sizeof(int64_t) + 2] = '\xff';
    memcpy((void*)(k + sizeof(char)*2 + sizeof(int64_t) + 3*sizeof(char)), &to, sizeof(int64_t));
    return key;
}

void Index::dump(ostream& out) {
    leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        
        out << it->key().ToString() << ": "  << it->value().ToString() << endl;
    }
    assert(it->status().ok());  // Check for any errors found during the scan
    delete it;
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
    db->Put(leveldb::WriteOptions(), key_for_edge_from_to(edge.from(), edge.to()), data);
    // only store in from_to key
    string null_data;
    db->Put(leveldb::WriteOptions(), key_for_edge_to_from(edge.to(), edge.from()), null_data);
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

leveldb::Status Index::get_node(int64_t id, Node& node) {
    string value;
    leveldb::Status s = db->Get(leveldb::ReadOptions(), key_for_node(id), &value);
    if (s.ok()) {
        node.ParseFromString(value);
    }
    return s;
}

leveldb::Status Index::get_edge(int64_t from, int64_t to, Edge& edge) {
    string value;
    leveldb::Status s = db->Get(leveldb::ReadOptions(), key_for_edge_from_to(from, to), &value);
    if (s.ok()) {
        edge.ParseFromString(value);
    }
    return s;
}

void index_kmers(VariantGraph& graph, int kmer_size = 15) {
// 
}

void index_positions(VariantGraph& graph, map<long, Node*>& node_path, map<long, Edge*>& edge_path) {

}

}
