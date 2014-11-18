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
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    k[2] = '\xff';
    memcpy((void*)(k + sizeof(char)*3), &id, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_from_to(int64_t from, int64_t to) {
    string key;
    key.resize(5*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    k[2] = '\xff';
    memcpy((void*)(k + sizeof(char)*3), &from, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = '\xff';
    k[3 + sizeof(int64_t) + 1] = 'f';
    k[3 + sizeof(int64_t) + 2] = '\xff';
    memcpy((void*)(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char)), &to, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_to_from(int64_t to, int64_t from) {
    string key;
    key.resize(5*sizeof(char) + 3*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = '\xff';
    k[1] = 'g'; // graph elements
    k[2] = '\xff';
    memcpy((void*)(k + sizeof(char)*3), &from, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = '\xff';
    k[3 + sizeof(int64_t) + 1] = 't';
    k[3 + sizeof(int64_t) + 2] = '\xff';
    memcpy((void*)(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char)), &to, sizeof(int64_t));
    return key;
}

string Index::entry_to_string(const string& key, const string& value) {
    char type = key[1];
    switch (type) {
    case 'g':
        return graph_entry_to_string(key, value);
        break;
    case 'k':
        return kmer_entry_to_string(key, value);
        break;
    case 'p':
        return position_entry_to_string(key, value);
        break;
    case 'm':
        return metadata_entry_to_string(key, value);
        break;
    default:
        break;
    }
}

string Index::graph_entry_to_string(const string& key, const string& value) {
    // do we have a node or edge?
    stringstream s;
    const char* k = key.c_str();
    if (key.size() == (3*sizeof(char) + sizeof(int64_t))) {
        // it's a node
        int64_t id;
        memcpy((void*) &id, (k + 3*sizeof(char)), sizeof(int64_t));
        Node node;
        node.ParseFromString(value);
        char *json2 = pb2json(node);
        s << "{\"key\":\"+g+" << id << "\", \"value\":"<<json2 << "}";
        free(json2);
    } else {
        int64_t id1, id2;
        memcpy(&id1, (k + 3*sizeof(char)), sizeof(int64_t));
        memcpy(&id2, (k + 3*sizeof(char)+sizeof(int64_t)+3*sizeof(char)), sizeof(int64_t));
        char type = k[3*sizeof(char)+sizeof(int64_t)+1*sizeof(char)];
        if (type == 'f') { // from
            Edge edge;
            edge.ParseFromString(value);
            char *json2 = pb2json(edge);
            s << "{\"key\":\"+g+" << id1 << "+f+" << id2 << "\", \"value\":"<<json2 << "}";
            free(json2);
        } else {
            s << "{\"key\":\"+g+" << id1 << "+t+" << id2 << "\", \"value\":null}";
        }
    }
    return s.str();
}

string Index::kmer_entry_to_string(const string& key, const string& value) {
}

string Index::position_entry_to_string(const string& key, const string& value) {
}

string Index::metadata_entry_to_string(const string& key, const string& value) {
}


void Index::dump(ostream& out) {
    leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        out << entry_to_string(it->key().ToString(), it->value().ToString()) << endl;
        //out << it->key().ToString() << ": "  << it->value().ToString() << endl;
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
