#include "index.hpp"

namespace vg {

using namespace std;

Index::Index(string& dir) : name(dir) {
    start_sep = '\x00';
    end_sep = '\xff';
    options.create_if_missing = true;
    //options.env->SetBackgroundThreads(omp_get_num_procs());
    //options.compression = rocksdb::kBZip2Compression;
    options.compression = rocksdb::kZlibCompression;
    options.compaction_style = rocksdb::kCompactionStyleLevel;
    options.IncreaseParallelism(omp_get_num_procs());
    options.write_buffer_size = 1024*1024*16; // 16mb
    //options.error_if_exists = true;
    rocksdb::Status status = rocksdb::DB::Open(options, name, &db);
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
    memcpy(k + sizeof(char)*3, &id, sizeof(int64_t));
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
    memcpy(k + sizeof(char)*3, &from, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[3 + sizeof(int64_t) + 1] = 'f';
    k[3 + sizeof(int64_t) + 2] = start_sep;
    memcpy(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char), &to, sizeof(int64_t));
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
    memcpy(k + sizeof(char)*3, &to, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[3 + sizeof(int64_t) + 1] = 't';
    k[3 + sizeof(int64_t) + 2] = start_sep;
    memcpy(k + sizeof(char)*3 + sizeof(int64_t) + 3*sizeof(char), &from, sizeof(int64_t));
    return key;
}

const string Index::key_for_kmer(const string& kmer, int64_t id) {
    id = htobe64(id);
    string key;
    key.resize(4*sizeof(char) + kmer.size() + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'k'; // kmers
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, kmer.c_str(), kmer.size());
    k[sizeof(char)*3 + kmer.size()] = start_sep;
    memcpy(k + sizeof(char)*4 + kmer.size(), &id, sizeof(int64_t));
    return key;
}

const string Index::key_prefix_for_kmer(const string& kmer) {
    string key;
    key.resize(3*sizeof(char) + kmer.size());
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'k'; // kmers
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, kmer.c_str(), kmer.size());
    return key;
}

const string Index::key_for_metadata(const string& tag) {
    string key;
    key.resize(3*sizeof(char) + tag.size());
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'm'; // metadata
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, tag.c_str(), tag.size());
    return key;
}

const string Index::key_prefix_for_edges_from_node(int64_t from) {
    string key = key_for_node(from);
    key.resize(key.size() + 2);
    key[key.size() - 2] = start_sep;
    key[key.size() - 1] = 'f';
    return key;
}

const string Index::key_prefix_for_edges_to_node(int64_t to) {
    string key = key_for_node(to);
    key.resize(key.size() + 2);
    key[key.size() - 2] = start_sep;
    key[key.size() - 1] = 't';
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
    memcpy(&id, (k + 3*sizeof(char)), sizeof(int64_t));
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

void Index::parse_kmer(const string& key, const string& value, string& kmer, int64_t& id, int32_t& pos) {
    const char* k = key.c_str();
    kmer = string(k+3*sizeof(char));
    memcpy(&id, k+4*sizeof(char)+kmer.size(), sizeof(int64_t));
    id = be64toh(id);
    memcpy(&pos, (char*)value.c_str(), sizeof(int32_t));
}

string Index::kmer_entry_to_string(const string& key, const string& value) {
    stringstream s;
    int64_t id;
    int32_t pos;
    string kmer;
    parse_kmer(key, value, kmer, id, pos);
    s << "{\"key\":\"+k+" << kmer << "+" << id << "\", \"value\":"<< pos << "}";
    return s.str();
}

string Index::position_entry_to_string(const string& key, const string& value) {
}

string Index::metadata_entry_to_string(const string& key, const string& value) {
    stringstream s;
    s << "{\"key\":\"" << "+" << key[1] << "+" << key.substr(2) << "\", \"value\":\""<< value << "\"}";
    return s.str();
}


void Index::dump(ostream& out) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        out << entry_to_string(it->key().ToString(), it->value().ToString()) << endl;
    }
    assert(it->status().ok());  // Check for any errors found during the scan
    delete it;
}

void Index::put_node(const Node* node) {
    string data;
    node->SerializeToString(&data);
    string key = key_for_node(node->id());
    db->Put(rocksdb::WriteOptions(), key, data);
}

void Index::put_edge(const Edge* edge) {
    string data;
    edge->SerializeToString(&data);
    db->Put(rocksdb::WriteOptions(), key_for_edge_from_to(edge->from(), edge->to()), data);
    // only store in from_to key
    string null_data;
    db->Put(rocksdb::WriteOptions(), key_for_edge_to_from(edge->to(), edge->from()), null_data);
}

void Index::put_metadata(const string& tag, const string& data) {
    string key = key_for_metadata(tag);
    db->Put(rocksdb::WriteOptions(), key, data);
}

void Index::load_graph(VG& graph) {
    graph.for_each_node_parallel([this](Node* n) { put_node(n); });
    graph.for_each_edge_parallel([this](Edge* e) { put_edge(e); });
}

rocksdb::Status Index::get_node(int64_t id, Node& node) {
    string value;
    rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key_for_node(id), &value);
    if (s.ok()) {
        node.ParseFromString(value);
    }
    return s;
}

rocksdb::Status Index::get_edge(int64_t from, int64_t to, Edge& edge) {
    string value;
    rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key_for_edge_from_to(from, to), &value);
    if (s.ok()) {
        edge.ParseFromString(value);
    }
    return s;
}

void Index::expand_context(VG& graph, int steps = 1) {
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

void Index::get_context(int64_t id, VG& graph) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    string key_start = key_for_node(id).substr(0,3+sizeof(int64_t));
    rocksdb::Slice start = rocksdb::Slice(key_start);
    string key_end = key_start+end_sep;
    rocksdb::Slice end = rocksdb::Slice(key_end);
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

void Index::get_kmer_subgraph(const string& kmer, VG& graph) {
    // get the nodes in the kmer subgraph
    auto add_node_matching_kmer = [&graph, this](string& key, string& value) {
        int64_t id;
        string kmer;
        int32_t pos;
        parse_kmer(key, value, kmer, id, pos);
        Node node;
        get_node(id, node);
        graph.add_node(node);
    };
    string start = key_prefix_for_kmer(kmer);
    string end = start + end_sep;
    start = start + start_sep;
    // apply to the range matching the kmer in the db
    for_range(start, end, add_node_matching_kmer);

    auto add_edges_from_index = [&graph, this](Node* n) {
        vector<Edge> edges;
        get_edges_from(n->id(), edges);
        get_edges_to(n->id(), edges);
        for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
            if (graph.has_node(e->to()) && graph.has_node(e->from())) {
                graph.add_edge(*e);
            }
        }
    };
    // add the edges between the matching nodes from the index
    graph.for_each_node(add_edges_from_index);
}

void Index::get_edges_from(int64_t from, vector<Edge>& edges) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    string key_start = key_prefix_for_edges_from_node(from);
    rocksdb::Slice start = rocksdb::Slice(key_start);
    string key_end = key_start+end_sep;
    rocksdb::Slice end = rocksdb::Slice(key_end);
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string s = it->key().ToString();
        char keyt = graph_key_type(s);
        switch (keyt) {
        case 'f': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            edges.push_back(edge);
        } break;
        default:
            // there should only be edges from here
            cerr << keyt << endl;
            assert(false);
            break;
        }
    }
}

void Index::get_edges_to(int64_t to, vector<Edge>& edges) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    string key_start = key_prefix_for_edges_to_node(to);
    rocksdb::Slice start = rocksdb::Slice(key_start);
    string key_end = key_start+end_sep;
    rocksdb::Slice end = rocksdb::Slice(key_end);
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string s = it->key().ToString();
        char keyt = graph_key_type(s);
        switch (keyt) {
        case 't': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            get_edge(id2, id1, edge);
            edges.push_back(edge);
        } break;
        default:
            // there should only be edges from here
            cerr << keyt << endl;
            assert(false);
            break;
        }
    }
}

void Index::put_kmer(const string& kmer,
                     const int64_t id,
                     const int32_t pos) {
    string key = key_for_kmer(kmer, id);
    string data(sizeof(int32_t), '\0');
    memcpy((char*)data.c_str(), &pos, sizeof(int32_t));
    rocksdb::Status s = db->Put(rocksdb::WriteOptions(), key, data);
    if (!s.ok()) { cerr << "put of " << kmer << " " << id << "@" << pos << " failed" << endl; exit(1); }
}

void Index::batch_kmer(const string& kmer,
                       const int64_t id,
                       const int32_t pos,
                       rocksdb::WriteBatch& batch) {
    string key = key_for_kmer(kmer, id);
    string data(sizeof(int32_t), '\0');
    memcpy((char*) data.c_str(), &pos, sizeof(int32_t));
    batch.Put(key, data);
}

void Index::populate_matches(Matches& matches, hash_map<Node*, int>& kmer_node_pos) {
    for (hash_map<Node*, int>::iterator m = kmer_node_pos.begin(); m != kmer_node_pos.end(); ++m) {
        Node* n = m->first;
        int pos = m->second;
        Match* match = matches.add_match();
        match->set_node_id(n->id());
        match->set_position(pos);
    }
}

void Index::store_kmers(string_hash_map<string, hash_map<Node*, int> >& kmer_map) {
    rocksdb::WriteBatch batch;
    for (string_hash_map<string, hash_map<Node*, int> >::iterator k = kmer_map.begin();
         k != kmer_map.end(); ++k) {
        const string& kmer = k->first;
        hash_map<Node*, int>& kmer_node_pos = k->second;
        for (auto kv : kmer_node_pos) {
            batch_kmer(kmer, kv.first->id(), kv.second, batch);
        }
    }
    rocksdb::Status s = db->Write(rocksdb::WriteOptions(), &batch);
    if (!s.ok()) cerr << "an error occurred while inserting kmers" << endl;
}

void Index::for_range(string& key_start, string& key_end,
                      std::function<void(string&, string&)> lambda) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    rocksdb::Slice start = rocksdb::Slice(key_start);
    rocksdb::Slice end = rocksdb::Slice(key_end);
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string key = it->key().ToString();
        string value = it->value().ToString();
        lambda(key, value);
    }
}

void Index::remember_kmer_size(int size) {
    stringstream s;
    s << "k=" << size;
    put_metadata(s.str(), "");
}

set<int> Index::stored_kmer_sizes(void) {
    set<int> sizes;
    auto lambda = [&sizes](string& key, string& value) {
        sizes.insert(atoi(key.substr(5).c_str()));
    };
    string start = key_for_metadata("k=");
    string end = start + end_sep;
    start = start + start_sep;
    for_range(start, end, lambda);
    return sizes;
}


void index_positions(VG& graph, map<long, Node*>& node_path, map<long, Edge*>& edge_path) {

}

}
