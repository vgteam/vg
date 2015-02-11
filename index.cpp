#include "index.hpp"

namespace vg {

using namespace std;

Index::Index(void) {

    start_sep = '\x00';
    end_sep = '\xff';
    write_options = rocksdb::WriteOptions();

    threads = 1;
#pragma omp parallel
    {
#pragma omp master
        threads = omp_get_num_threads();
    }

}

// from https://github.com/facebook/rocksdb/blob/master/utilities/spatialdb/spatial_db.cc#L660-L716

rocksdb::Options Index::GetOptions(void) {

    rocksdb::Options options;

    options.create_if_missing = true;
    options.max_open_files = 100000;
    options.compression = rocksdb::kSnappyCompression;
    options.compaction_style = rocksdb::kCompactionStyleLevel;
    options.IncreaseParallelism(threads);
    options.max_background_compactions = threads;
    options.max_background_flushes = threads;

    options.num_levels = 2;
    options.target_file_size_base = (long) 64 * 1024 * 1024 * 1024; // 64G

    if (bulk_load) {
        options.PrepareForBulkLoad();
        options.write_buffer_size = 1024 * 1024 * 256;
        //options.target_file_size_base = 1024 * 1024 * 512;
        //options.target_file_size_base = (long) 64 * 1024 * 1024 * 1024; // 64G
        options.IncreaseParallelism(threads);
        options.max_background_compactions = threads;
        options.max_background_flushes = threads;
        options.max_write_buffer_number = threads;
        options.compaction_style = rocksdb::kCompactionStyleNone;
        options.memtable_factory.reset(new rocksdb::VectorRepFactory(1000));
    }

    options.compression_per_level.resize(
        options.num_levels);
    for (int i = 0; i < options.num_levels; ++i) {
        if (i == 0) {
            options.compression_per_level[i] = rocksdb::kSnappyCompression;
        } else {
            options.compression_per_level[i] = rocksdb::kZlibCompression;
        }
    }

    return options;
}

void Index::open(const std::string& dir, bool read_only) {

    name = dir;
    db_options = GetOptions();

    rocksdb::Status s;
    if (read_only) {
        s = rocksdb::DB::OpenForReadOnly(db_options, name, &db);
    } else {
        s = rocksdb::DB::Open(db_options, name, &db);
    }
    if (!s.ok()) {
        throw indexOpenException();
    }

}

void Index::open_read_only(string& dir) {
    bulk_load = false;
    open(dir, true);
}

void Index::open_for_write(string& dir) {
    bulk_load = false;
    open(dir, false);
}

void Index::open_for_bulk_load(string& dir) {
    bulk_load = true;
    open(dir, false);
}

Index::~Index(void) {
    close();
}

void Index::close(void) {
    flush();
    delete db;
}

void Index::flush(void) {
    db->Flush(rocksdb::FlushOptions());
}

void Index::compact(void) {
    db->CompactRange(NULL, NULL);
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
        //edge.ParseFromString(value);
        edge.set_from(id1);
        edge.set_to(id2);
    } else {
        // XXX big hack
        // we don't need to properly parse the edge until we stash data there
        // if we do stash data there, it's probably best to avoid the second lookup and duplicate the edge
        // or normalize things and store the edges in their own subset of the index
        edge.set_from(id2);
        edge.set_to(id1);
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
        //get_edge(id2, id1, edge);
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
    db->Put(write_options, key, data);
}

void Index::batch_node(const Node* node, rocksdb::WriteBatch& batch) {
    string data;
    node->SerializeToString(&data);
    string key = key_for_node(node->id());
    batch.Put(key, data);
}

void Index::put_edge(const Edge* edge) {
    string data;
    edge->SerializeToString(&data);
    db->Put(write_options, key_for_edge_from_to(edge->from(), edge->to()), data);
    // only store in from_to key
    string null_data;
    db->Put(write_options, key_for_edge_to_from(edge->to(), edge->from()), null_data);
}

void Index::batch_edge(const Edge* edge, rocksdb::WriteBatch& batch) {
    string data;
    edge->SerializeToString(&data);
    db->Put(write_options, key_for_edge_from_to(edge->from(), edge->to()), data);
    // only store in from_to key
    string null_data;
    batch.Put(key_for_edge_to_from(edge->to(), edge->from()), null_data);
}

void Index::put_metadata(const string& tag, const string& data) {
    string key = key_for_metadata(tag);
    db->Put(write_options, key, data);
}

void Index::load_graph(VG& graph) {
    // a bit of a hack--- the logging only works with for_each_*parallel
    // also the high parallelism may be causing issues
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    //omp_set_num_threads(1);
    graph.create_progress("indexing nodes of " + graph.name, graph.graph.node_size());
    vector<rocksdb::WriteBatch> batch;
    batch.resize(thread_count);
    graph.for_each_node_parallel([this, &batch](Node* n) { batch_node(n, batch.at(omp_get_thread_num())); });
    graph.destroy_progress();
    graph.create_progress("indexing edges of " + graph.name, graph.graph.edge_size());
    graph.for_each_edge_parallel([this, &batch](Edge* e) { batch_edge(e, batch.at(omp_get_thread_num())); });
    for (auto& b : batch) {
        rocksdb::Status s = db->Write(rocksdb::WriteOptions(), &b);
    }
    graph.destroy_progress();
    omp_set_num_threads(thread_count);
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
    for (int step = 0; step < steps; ++step) {
        set<int64_t> ids;
        graph.for_each_edge([this, &graph, &ids](Edge* edge) {
                if (!graph.has_node(edge->from())) {
                    ids.insert(edge->from());
                }
                if (!graph.has_node(edge->to())) {
                    ids.insert(edge->to());
                }
            });
        for (auto id : ids) {
            get_context(id, graph);
        }
    }
}

void Index::get_connected_nodes(VG& graph) {
    graph.for_each_edge([this, &graph](Edge* edge) {
            if (!graph.has_node(edge->from())) {
                Node node;
                get_node(edge->from(), node);
                graph.add_node(node);
            }
            if (!graph.has_node(edge->to())) {
                Node node;
                get_node(edge->to(), node);
                graph.add_node(node);
            }
        });
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
        } break;
        case 't': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            // avoid a second lookup
            // probably we should index these twice and pay the penalty on *write* rather than read
            //get_edge(id2, id1, edge);
            graph.add_edge(edge);

        } break;
        default:
            cerr << "vg::Index unrecognized key type " << keyt << endl;
            exit(1);
            break;
        }
    }
    delete it;
}

void Index::get_range(int64_t from_id, int64_t to_id, VG& graph) {
    auto handle_entry = [this, &graph](string& key, string& value) {
        char keyt = graph_key_type(key);
        switch (keyt) {
        case 'n': {
            Node node;
            node.ParseFromString(value);
            graph.add_node(node);
        } break;
        case 'f': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(key, value, type, id1, id2, edge);
            graph.add_edge(edge);
        } break;
        case 't': {
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(key, value, type, id1, id2, edge);
            graph.add_edge(edge);

        } break;
        default:
            cerr << "vg::Index unrecognized key type " << keyt << endl;
            exit(1);
            break;
        }
    };
    for_graph_range(from_id, to_id, handle_entry);
}

void Index::get_kmer_subgraph(const string& kmer, VG& graph) {
    // get the nodes in the kmer subgraph
    for_kmer_range(kmer, [&graph, this](string& key, string& value) {
            int64_t id;
            string kmer;
            int32_t pos;
            parse_kmer(key, value, kmer, id, pos);
            get_context(id, graph);
        });
}

void Index::get_kmer_positions(const string& kmer, map<int64_t, vector<int32_t> >& positions) {
    for_kmer_range(kmer, [&positions, this](string& key, string& value) {
            int64_t id;
            string kmer;
            int32_t pos;
            parse_kmer(key, value, kmer, id, pos);
            positions[id].push_back(pos);
        });
}

void Index::get_kmer_positions(const string& kmer, map<string, vector<pair<int64_t, int32_t> > >& positions) {
    for_kmer_range(kmer, [&positions, this](string& key, string& value) {
            int64_t id;
            string kmer;
            int32_t pos;
            parse_kmer(key, value, kmer, id, pos);
            positions[kmer].push_back(make_pair(id, pos));
        });
}

void Index::for_kmer_range(const string& kmer, function<void(string&, string&)> lambda) {
    string start = key_prefix_for_kmer(kmer);
    string end = start + end_sep;
    start = start + start_sep;
    // apply to the range matching the kmer in the db
    for_range(start, end, lambda);
}

void Index::for_graph_range(int64_t from_id, int64_t to_id, function<void(string&, string&)> lambda) {
    string start = key_for_node(from_id);
    string end = key_for_node(to_id+1);
    // apply to the range matching the kmer in the db
    for_range(start, end, lambda);
}

uint64_t Index::approx_size_of_kmer_matches(const string& kmer) {
    uint64_t size;
    string start = key_prefix_for_kmer(kmer);
    string end = start + end_sep;
    rocksdb::Range range = rocksdb::Range(start, end);
    db->GetApproximateSizes(&range, 1, &size);
    return size;
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
    rocksdb::Status s = db->Put(write_options, key, data);
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

void Index::store_batch(map<string, string>& items) {
    rocksdb::WriteBatch batch;
    for (auto& i : items) {
        const string& k = i.first;
        const string& v = i.second;
        batch.Put(k, v);
    }
    rocksdb::Status s = db->Write(write_options, &batch);
    if (!s.ok()) cerr << "an error occurred while inserting items" << endl;
}

void Index::for_all(std::function<void(string&, string&)> lambda) {
    string start(1, start_sep);
    string end(1, end_sep);
    for_range(start, end, lambda);
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

// todo, get range estimated size

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
