#include "index.h"

namespace vg {

using namespace std;

Index::Index(string& name) {
    start_sep = '\x00';
    end_sep = '\xff';
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

// todo: replace with union / struct
const string Index::key_for_node(int64_t id) {
    string key;
    id = htobe64(id);
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy((void*)(k + sizeof(char)*3), &id, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_from_to(int64_t from, int64_t to) {
    // reverse endianness for sorting
    to = htobe64(to);
    from = htobe64(from);
    string key;
    key.resize(6*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy((void*)(k + sizeof(char)*3), &from, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[3 + sizeof(int64_t) + 1] = 'f';
    k[3 + sizeof(int64_t) + 2] = start_sep;
    memcpy((void*)(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char)), &to, sizeof(int64_t));
    return key;
}

const string Index::key_for_edge_to_from(int64_t to, int64_t from) {
    // reverse endianness for sorting
    to = htobe64(to);
    from = htobe64(from);
    string key;
    key.resize(5*sizeof(char) + 3*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy((void*)(k + sizeof(char)*3), &to, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[3 + sizeof(int64_t) + 1] = 't';
    k[3 + sizeof(int64_t) + 2] = start_sep;
    memcpy((void*)(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char)), &from, sizeof(int64_t));
    return key;
}

const string Index::key_for_kmer(const string& kmer) {
    string key;
    key.resize(3*sizeof(char) + kmer.size());
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'k'; // kmers
    k[2] = start_sep;
    memcpy((char*)(k + sizeof(char)*3), (char*)kmer.c_str(), kmer.size());
    return key;
}

char Index::graph_key_type(string& key) {
    if (key.size() == (3*sizeof(char) + sizeof(int64_t))) return 'n';
    return key.c_str()[4*sizeof(char) + sizeof(int64_t)];
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

void Index::parse_node(const string& key, const string& value, int64_t& id, Node& node) {
    const char* k = key.c_str();
    memcpy((void*) &id, (k + 3*sizeof(char)), sizeof(int64_t));
    id = be64toh(id);
    node.ParseFromString(value);
}

void Index::parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge) {
    const char* k = key.c_str();
    memcpy(&id1, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&id2, (k + 3*sizeof(char)+sizeof(int64_t)+3*sizeof(char)), sizeof(int64_t));
    id1 = be64toh(id1);
    id2 = be64toh(id2);
    type = k[3*sizeof(char)+sizeof(int64_t)+1*sizeof(char)];
    if (type == 'f') {
        edge.ParseFromString(value);
    }
}

string Index::graph_entry_to_string(const string& key, const string& value) {
    // do we have a node or edge?
    stringstream s;
    switch (graph_key_type((string&)key)) {
    case 'n': {
        // it's a node
        int64_t id;
        Node node;
        parse_node(key, value, id, node);
        char *json = pb2json(node);
        s << "{\"key\":\"+g+" << id << "+n\", \"value\":"<<json << "}";
        free(json);
    } break;
    case 'f': {
        Edge edge;
        int64_t id1, id2;
        char type;
        parse_edge(key, value, type, id1, id2, edge);
        char *json = pb2json(edge);
        s << "{\"key\":\"+g+" << id1 << "+f+" << id2 << "\", \"value\":"<<json << "}";
        free(json);
    } break;
    case 't': {
        Edge edge;
        int64_t id1, id2;
        char type;
        parse_edge(key, value, type, id1, id2, edge);
        get_edge(id2, id1, edge);
        char *json = pb2json(edge);
        s << "{\"key\":\"+g+" << id1 << "+t+" << id2 << "\", \"value\":"<<json << "}";
        free(json);
    } break;
    }
    return s.str();
}

void Index::parse_kmer(const string& key, const string& value, string& kmer, Matches& matches) {
    const char* k = key.c_str();
    kmer = string(k+3*sizeof(char));
    matches.ParseFromString(value);
}

string Index::kmer_entry_to_string(const string& key, const string& value) {
    stringstream s;
    Matches matches;
    string kmer;
    parse_kmer(key, value, kmer, matches);
    char *json = pb2json(matches);
    s << "{\"key\":\"+k+" << kmer << "\", \"value\":"<<json << "}";
    free(json);
    return s.str();
}

string Index::position_entry_to_string(const string& key, const string& value) {
}

string Index::metadata_entry_to_string(const string& key, const string& value) {
}


void Index::dump(ostream& out) {
    leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        out << entry_to_string(it->key().ToString(), it->value().ToString()) << endl;
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
    for (int i = 0; i < g.node_size(); ++i) {
        put_node(g.node(i));
    }
    for (int i = 0; i < g.edge_size(); ++i) {
        put_edge(g.edge(i));
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

void Index::expand_context(VariantGraph& graph, int steps = 1) {
    Graph& g = graph.graph; // ugly
    for (int step = 0; step < steps; ++step) {
        vector<int64_t> ids;
        for (int i = 0; i < g.node_size(); ++i) {
            Node* node = g.mutable_node(i);
            ids.push_back(node->id());
        }
        for (vector<int64_t>::iterator id = ids.begin(); id != ids.end(); ++id) {
            get_context(*id, graph);
        }
    }
}

void Index::get_context(int64_t id, VariantGraph& graph) {
    leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
    string key_start = key_for_node(id).substr(0,3+sizeof(int64_t));
    leveldb::Slice start = leveldb::Slice(key_start);
    string key_end = key_start+end_sep;
    leveldb::Slice end = leveldb::Slice(key_end);
    // TODO ... uhhh return a valid graph
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string s = it->key().ToString();
        char keyt = graph_key_type(s);
        switch (keyt) {
        case 'n': {
            Node node;
            node.ParseFromString(it->value().ToString());
            graph.add_node(node);
        } break;
        case 'f': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            graph.add_edge(edge);
            // get to node
            Node node;
            get_node(id2, node);
            graph.add_node(node);
        } break;
        case 't': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            get_edge(id2, id1, edge);
            graph.add_edge(edge);
            // get from node
            Node node;
            get_node(id2, node);
            graph.add_node(node);
        } break;
        default:
            cerr << "vg::Index unrecognized key type " << keyt << endl;
            exit(1);
            break;
        }
    }
    delete it;
}

void Index::get_kmer_subgraph(const string& kmer, VariantGraph& graph) {
    string value;
    leveldb::Status s = db->Get(leveldb::ReadOptions(), key_for_kmer(kmer), &value);
    //if (!s.ok()) cerr << "read of kmer " << kmer << " is not OK" << endl;
    Matches matches;
    matches.ParseFromString(value);
    for (int i = 0; i < matches.match_size(); ++i) {
        Match* match = matches.mutable_match(i);
        get_context(match->node_id(), graph);
    }
}

void Index::get_edges_from(int64_t from, vector<Edge>& edges) {
}

void Index::get_edges_to(int64_t to, vector<Edge>& edges) {
}

void Index::put_kmer(const string& kmer, const Matches& matches) {
    string data;
    matches.SerializeToString(&data);
    string key = key_for_kmer(kmer);
    leveldb::Status s = db->Put(leveldb::WriteOptions(), key, data);
    if (!s.ok()) cerr << "put failed" << endl;
}

void Index::batch_kmer(const string& kmer, const Matches& matches, leveldb::WriteBatch& batch) {
    string data;
    matches.SerializeToString(&data);
    string key = key_for_kmer(kmer);
    batch.Put(key, data);
}

void Index::populate_matches(Matches& matches, map<Node*, int>& kmer_node_pos) {
    for (map<Node*, int>::iterator m = kmer_node_pos.begin(); m != kmer_node_pos.end(); ++m) {
        Node* n = m->first;
        int pos = m->second;
        Match* match = matches.add_match();
        match->set_node_id(n->id());
        match->set_position(pos);
    }
}

void Index::store_kmers(map<string, map<Node*, int> >& kmer_map) {
    leveldb::WriteBatch batch;
    for (map<string, map<Node*, int> >::iterator k = kmer_map.begin(); k != kmer_map.end(); ++k) {
        const string& kmer = k->first;
        map<Node*, int>& kmer_node_pos = k->second;
        Matches matches;
        populate_matches(matches, kmer_node_pos);
        batch_kmer(kmer, matches, batch);
    }
    leveldb::Status s = db->Write(leveldb::WriteOptions(), &batch);
    if (!s.ok()) cerr << "an error occurred while inserting kmers" << endl;
}

void index_positions(VariantGraph& graph, map<long, Node*>& node_path, map<long, Edge*>& edge_path) {

}

}
