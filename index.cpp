#include "index.hpp"

namespace vg {

using namespace std;

Index::Index(void) {

    start_sep = '\x00';
    end_sep = '\xff';
    write_options = rocksdb::WriteOptions();
    mem_env = false;
    block_cache_size = 1024 * 1024 * 1024; // 1GB

    threads = 1;
#pragma omp parallel
    {
#pragma omp master
        threads = omp_get_num_threads();
    }

}

rocksdb::Options Index::GetOptions(void) {

    rocksdb::Options options;

    if (mem_env) {
        options.env = rocksdb::NewMemEnv(options.env);
    }

    options.create_if_missing = true;
    options.max_open_files = 100000;
    options.compression = rocksdb::kSnappyCompression;
    options.compaction_style = rocksdb::kCompactionStyleLevel;
    // we are unlikely to reach either of these limits
    options.IncreaseParallelism(threads);
    options.max_background_flushes = threads;
    options.max_background_compactions = threads;

    options.num_levels = 2;
    options.target_file_size_base = (long) 1024 * 1024 * 512; // 512MB
    options.write_buffer_size = 1024 * 1024 * 256;

    // doesn't work this way
    rocksdb::BlockBasedTableOptions topt;
    topt.filter_policy.reset(rocksdb::NewBloomFilterPolicy(16, true));
    //topt.block_cache = rocksdb::NewLRUCache(block_cache_size, 32);
    //topt.no_block_cache = true;
    options.table_factory.reset(NewBlockBasedTableFactory(topt));

    if (bulk_load) {
        options.PrepareForBulkLoad();
        options.max_write_buffer_number = threads;
        options.max_background_flushes = threads;
        options.max_background_compactions = threads;
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
        //s = rocksdb::DB::Open(db_options, name, &db);
        s = rocksdb::DB::OpenForReadOnly(db_options, name, &db);
    } else {
        s = rocksdb::DB::Open(db_options, name, &db);
    }
    if (!s.ok()) {
        throw indexOpenException();
    }
    is_open = true;

}

void Index::open_read_only(string& dir) {
    bulk_load = false;
    //mem_env = true;
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
    if (is_open) {
        close();
    }
}

void Index::close(void) {
    flush();
    delete db;
    is_open = false;
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
    key.resize(5*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &id, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[3 + sizeof(int64_t) + 1] = 'n';
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

const string Index::key_for_node_path_position(int64_t node_id, int64_t path_id, int64_t path_pos) {
    node_id = htobe64(node_id);
    path_id = htobe64(path_id);
    path_pos = htobe64(path_pos);
    string key;
    key.resize(7*sizeof(char) + 3*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[4 + sizeof(int64_t)] = 'p';
    k[5 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*6 + sizeof(int64_t), &path_id, sizeof(int64_t));
    k[6 + 2*sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*7 + 2*sizeof(int64_t), &path_pos, sizeof(int64_t));
    return key;
}

const string Index::key_prefix_for_node_path(int64_t node_id, int64_t path_id) {
    node_id = htobe64(node_id);
    path_id = htobe64(path_id);
    string key;
    key.resize(6*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[4 + sizeof(int64_t)] = 'p';
    k[5 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*6 + sizeof(int64_t), &path_id, sizeof(int64_t));
    return key;
}

const string Index::key_for_path_position(int64_t path_id, int64_t path_pos, int64_t node_id) {
    node_id = htobe64(node_id);
    path_id = htobe64(path_id);
    path_pos = htobe64(path_pos);
    string key;
    key.resize(5*sizeof(char) + 3*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'p'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &path_id, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*4 + sizeof(int64_t), &path_pos, sizeof(int64_t));
    k[4 + 2*sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*5 + 2*sizeof(int64_t), &node_id, sizeof(int64_t));
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
    string key = key_for_edge_from_to(from, 0);
    return key.substr(0, key.size()-sizeof(int64_t));
}

const string Index::key_prefix_for_edges_to_node(int64_t to) {
    string key = key_for_edge_to_from(to, 0);
    return key.substr(0, key.size()-sizeof(int64_t));
}

char Index::graph_key_type(string& key) {
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
        return path_position_to_string(key, value);
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
    case 'p': {
        s << node_path_to_string(key, value);
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

void Index::parse_node_path(const string& key, const string& value,
                            int64_t& node_id, int64_t& path_id, int64_t& path_pos, Mapping& mapping) {
    const char* k = key.c_str();
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&path_id, (k + 6*sizeof(char)+sizeof(int64_t)), sizeof(int64_t));
    memcpy(&path_pos, (k + 7*sizeof(char)+2*sizeof(int64_t)), sizeof(int64_t));
    node_id = be64toh(node_id);
    path_id = be64toh(path_id);
    path_pos = be64toh(path_pos);
    mapping.ParseFromString(value);
}

void Index::parse_path_position(const string& key, const string& value,
                                int64_t& path_id, int64_t& path_pos, int64_t& node_id, Mapping& mapping) {
    const char* k = key.c_str();
    memcpy(&path_id, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&path_pos, (k + 4*sizeof(char)+sizeof(int64_t)), sizeof(int64_t));
    memcpy(&node_id, (k + 5*sizeof(char)+2*sizeof(int64_t)), sizeof(int64_t));
    node_id = be64toh(node_id);
    path_id = be64toh(path_id);
    path_pos = be64toh(path_pos);
    mapping.ParseFromString(value);
}

string Index::node_path_to_string(const string& key, const string& value) {
    Mapping mapping;
    int64_t node_id, path_id, path_pos;
    parse_node_path(key, value, node_id, path_id, path_pos, mapping);
    stringstream s;
    char *json = pb2json(mapping);
    s << "{\"key\":\"+g+" << node_id << "+p+" << path_id << "+" << path_pos << "\", \"value\":"<<json << "}";
    free(json);
    return s.str();
}

string Index::path_position_to_string(const string& key, const string& value) {
    Mapping mapping;
    int64_t node_id, path_id, path_pos;
    parse_path_position(key, value, path_id, path_pos, node_id, mapping);
    stringstream s;
    char *json = pb2json(mapping);
    s << "{\"key\":\"+p+" << path_id << "+" << path_pos << "+" << node_id << "\", \"value\":"<<json << "}";
    free(json);
    return s.str();
}

string Index::metadata_entry_to_string(const string& key, const string& value) {
    stringstream s;
    string prefix = key.substr(3);
    string val = value;
    if (prefix == "max_path_id"
        || prefix.substr(0,9) == "path_name") {
        stringstream v;
        int64_t id;
        memcpy(&id, (char*)value.c_str(), sizeof(int64_t));
        v << id;
        val = v.str();
    } else if (prefix.substr(0,7) == "path_id") {
        stringstream v;
        int64_t id;
        memcpy(&id, ((char*)prefix.c_str())+7, sizeof(int64_t));
        v << id;
        prefix = prefix.substr(0,7) + "+" + v.str();
    }
    s << "{\"key\":\"" << "+" << key[1] << "+" << prefix << "\", \"value\":\""<< val << "\"}";
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

void Index::put_node_path(int64_t node_id, int64_t path_id, int64_t path_pos, const Mapping& mapping) {
    string data;
    mapping.SerializeToString(&data);
    db->Put(write_options, key_for_node_path_position(node_id, path_id, path_pos), data);
}

void Index::put_path_position(int64_t path_id, int64_t path_pos, int64_t node_id, const Mapping& mapping) {
    string data;
    mapping.SerializeToString(&data);
    db->Put(write_options, key_for_path_position(path_id, path_pos, node_id), data);
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
    omp_set_num_threads(1);
    graph.create_progress("indexing nodes of " + graph.name, graph.graph.node_size());
    rocksdb::WriteBatch batch;
    graph.for_each_node_parallel([this, &batch](Node* n) { batch_node(n, batch); });
    graph.destroy_progress();
    graph.create_progress("indexing edges of " + graph.name, graph.graph.edge_size());
    graph.for_each_edge_parallel([this, &batch](Edge* e) { batch_edge(e, batch); });
    rocksdb::Status s = db->Write(rocksdb::WriteOptions(), &batch);
    omp_set_num_threads(thread_count);
}

void Index::load_paths(VG& graph) {
    graph.destroy_progress();
    graph.create_progress("indexing paths of " + graph.name, graph.paths._paths->size());
    store_paths(graph);
    graph.destroy_progress();
}

int64_t Index::get_max_path_id(void) {
    string data;
    int64_t id;
    rocksdb::Status s = get_metadata("max_path_id", data);
    if (!s.ok()) {
        id = 0;
        put_max_path_id(id);
    } else {
        memcpy(&id, data.c_str(), sizeof(int64_t));
    }
    return id;
}

void Index::put_max_path_id(int64_t id) {
    string data;
    data.resize(sizeof(int64_t));
    memcpy((char*)data.c_str(), &id, sizeof(int64_t));
    put_metadata("max_path_id", data);
}

int64_t Index::new_path_id(const string& path_name) {
    int64_t max_id = get_max_path_id();
    int64_t new_id = max_id + 1;
    put_max_path_id(new_id);
    put_path_id_to_name(new_id, path_name);
    put_path_name_to_id(new_id, path_name);
    return new_id;
}

string Index::path_name_prefix(const string& name) {
    return "path_name" + start_sep + name;
}

string Index::path_id_prefix(int64_t id) {
    string prefix = "path_id" + start_sep;
    size_t prefix_size = prefix.size();
    prefix.resize(prefix.size() + sizeof(int64_t));
    memcpy((char*)prefix.c_str() + prefix_size, &id, sizeof(int64_t));
    return prefix;
}

void Index::put_path_id_to_name(int64_t id, const string& name) {
    put_metadata(path_id_prefix(id), name);
}

void Index::put_path_name_to_id(int64_t id, const string& name) {
    string data;
    data.resize(sizeof(int64_t));
    memcpy((char*)data.c_str(), &id, sizeof(int64_t));
    put_metadata(path_name_prefix(name), data);
}

string Index::get_path_name(int64_t id) {
    string data;
    get_metadata(path_id_prefix(id), data);
    return data;
}

int64_t Index::get_path_id(const string& name) {
    string data;
    get_metadata(path_name_prefix(name), data);
    int64_t id = 0;
    memcpy(&id, (char*)data.c_str(), sizeof(int64_t));
    return id;
}

void Index::store_paths(VG& graph) {
    for (auto& path : *graph.paths._paths) {
        store_path(graph, path);
    }
}

void Index::store_path(VG& graph, Path& path) {
    // get a new path id
    // if there is no name, cry
    if (!path.has_name()) {
        cerr << "[vg::Index] error, path has no name" << endl;
        exit(1);
    }
    // check if the path name/id mapping already exists
    int64_t path_id;
    path_id = get_path_id(path.name());
    // if it doesn't, create it
    if (!path_id) {
        path_id = new_path_id(path.name());
    }
    // keep track of position
    int64_t path_pos = 0;
    // for each node in the path
    for (int64_t i = 0; i < path.mapping_size(); ++i) {

        const Mapping& mapping = path.mapping(i);
        // put an entry in the path table
        put_path_position(path_id, path_pos, mapping.node_id(), mapping);
        // put an entry in the graph table
        put_node_path(mapping.node_id(), path_id, path_pos, mapping);

        // get the node, to find the size of this step
        Node node;
        get_node(mapping.node_id(), node);
        // TODO use the cigar... if there is one
        path_pos += node.sequence().size();

        graph.update_progress(graph.progress_count+1);
    }
}

rocksdb::Status Index::get_metadata(const string& key, string& data) {
    rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key_for_metadata(key), &data);
    return s;
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

int Index::get_node_path(int64_t node_id, int64_t path_id, int64_t& path_pos, Mapping& mapping) {
    string value;
    string key = key_prefix_for_node_path(node_id, path_id);
    string start = key + start_sep;
    string end = key + end_sep;
    // NB: uses the first position in the range
    // apply to the range matching the kmer in the db
    int count = 0;
    for_range(start, end, [this, &count, &node_id, &path_id, &path_pos, &mapping](string& key, string& value) {
            if (count == 0) {
                parse_node_path(key, value,
                                node_id, path_id,
                                path_pos, mapping);
            }
            ++count;
        });
    return count;
}

pair<list<int64_t>, int64_t> Index::get_nearest_node_prev_path_member(int64_t node_id, int64_t path_id, int64_t& path_pos, int max_steps) {
    list<int64_t> nullpath;
    list<int64_t> bpath;
    map<list<int64_t>, Node> nq;
    { // handle this node
        bpath.push_front(node_id);
        Node& node = nq[bpath];
        get_node(node_id, node);
        Mapping mapping;
        if (get_node_path(node_id, path_id, path_pos, mapping) > 0) {
            return make_pair(bpath, node_id);
        }
        Node n = node;
        nq.clear();
        nq[nullpath] = n;
    }
    // BFS back
    int steps_back = 0;
    while (steps_back++ < max_steps) {
        map<list<int64_t>, Node> cq;
        for (auto& n : nq) {
            Node& node = n.second;
            const list<int64_t>& path = n.first;
            vector<Edge> e_to;
            get_edges_to(node.id(), e_to);
            for (auto& edge : e_to) {
                int64_t id = edge.from();
                list<int64_t> npath = path;
                npath.push_front(id);
                Node& node = cq[npath];
                get_node(id, node);
                Mapping mapping;
                if (get_node_path(id, path_id, path_pos, mapping) > 0) {
                    path_pos += node.sequence().size();
                    return make_pair(npath, id);
                }
            }
        }
        nq = cq;
    }
    return make_pair(nullpath, 0);
}

pair<list<int64_t>, int64_t> Index::get_nearest_node_next_path_member(int64_t node_id, int64_t path_id, int64_t& path_pos, int max_steps) {
    list<int64_t> nullpath;
    list<int64_t> bpath;
    map<list<int64_t>, Node> nq;
    { // handle this node
        bpath.push_back(node_id);
        Node& node = nq[bpath];
        get_node(node_id, node);
        Mapping mapping;
        if (get_node_path(node_id, path_id, path_pos, mapping) > 0) {
            return make_pair(bpath, node_id);
        }
        Node n = node;
        nq.clear();
        nq[nullpath] = n;
    }
    // BFS back
    int steps_back = 0;
    while (steps_back++ < max_steps) {
        map<list<int64_t>, Node> cq;
        for (auto& n : nq) {
            Node& node = n.second;
            const list<int64_t>& path = n.first;
            vector<Edge> e_from;
            get_edges_from(node.id(), e_from);
            for (auto& edge : e_from) {
                int64_t id = edge.to();
                list<int64_t> npath = path;
                npath.push_back(id);
                Node& node = cq[npath];
                get_node(id, node);
                Mapping mapping;
                if (get_node_path(id, path_id, path_pos, mapping) > 0) {
                    return make_pair(npath, id);
                }
            }
        }
        nq = cq;
    }
    return make_pair(nullpath, 0);
}

bool Index::get_node_path_relative_position(int64_t node_id, int64_t path_id,
                                            list<int64_t>& path_prev, int64_t& prev_pos,
                                            list<int64_t>& path_next, int64_t& next_pos) {
    // scan the range before the node
    // start with our node, and walk back BFS until we find a node with a path
    // are any parents part of the path?
    list<int64_t> nullpath;

    auto null_pair = make_pair(nullpath, (int64_t)0);
    auto to_path_prev = get_nearest_node_prev_path_member(node_id, path_id, prev_pos);
    if (to_path_prev == null_pair) {
        cerr << "no to path" << endl;
        return false;
    } else {
        path_prev = to_path_prev.first;
    }

    auto to_path_next = get_nearest_node_next_path_member(node_id, path_id, next_pos);
    if (to_path_next == null_pair) {
        cerr << "no from path" << endl;
        return false;
    } else {
        path_next = to_path_next.first;
    }

    return true;
}

Mapping Index::path_relative_mapping(int64_t node_id, int64_t path_id,
                                     list<int64_t>& path_prev, int64_t& prev_pos,
                                     list<int64_t>& path_next, int64_t& next_pos) {
    Mapping mapping;
    mapping.set_node_id(node_id);
    if (get_node_path_relative_position(node_id, path_id,
                                        path_prev, prev_pos, path_next, next_pos)) {
        Edit* edit = mapping.add_edit();
        Node node; get_node(node_id, node);
        bool in_path = path_prev.back() == node_id && path_next.front() == node_id;
        int32_t to_length = node.sequence().size();
        int32_t from_length = in_path ? to_length : next_pos - prev_pos;
        if (from_length == to_length) {
            edit->set_from_length(from_length);
        } else {
            edit->set_from_length(from_length);
            edit->set_to_length(to_length);
        }
    }
    return mapping;
}

// transform the path into a path relative to another path (defined by path_name)
// source -> surjection (in path_name coordinate space)
// the product is equivalent to a pairwise alignment between this path and the other

// new approach
// get path sequence
// get graph component overlapping path
// removing elements which aren't in the path of interest
// realign to this graph
// cross fingers

            // what we need:
            /*
                         const string& refseq,
                         const int32_t refpos,
                         const string& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         const int32_t tlen);
            */
bool Index::surject_alignment(const Alignment& source,
                              set<string>& path_names,
                              Alignment& surjection,
                              string& path_name,
                              int64_t& path_pos,
                              int window) {
    VG graph;
    // get start and end nodes in path
    // get range between +/- window
    if (!source.has_path()) {
        return false;
    }
    int64_t from_id = source.path().mapping(0).node_id() - window;
    int64_t to_id = source.path().mapping(source.path().mapping_size()-1).node_id() + window;
    get_range(max((int64_t)0, from_id), to_id, graph);
    graph.remove_orphan_edges();
    //string source_seq = graph.path_sequence(source);
    // which path(s) did we keep?
    set<string> kept_paths;
    graph.keep_paths(path_names, kept_paths);
    //graph.serialize_to_file("surjection.vg"); // debugging
    //surjection.set_sequence(source.sequence());
    surjection = source;
    surjection.clear_path();
    graph.align(surjection);
    if (surjection.path().mapping_size() > 0 && kept_paths.size() == 1) {
        // determine the paths of the node we mapped into
        //  ... get the id of the first node, get the pahs of it
        assert(kept_paths.size() == 1);
        path_name = *kept_paths.begin();

        int64_t path_id = get_path_id(path_name);
        int64_t hit_id = surjection.path().mapping(0).node_id();
        int64_t pos = surjection.path().position();
        // we pick up positional information using the index
        int64_t prev_pos=0, next_pos=0;
        list<int64_t> path_prev, path_next;
        get_node_path_relative_position(hit_id, path_id, path_prev, prev_pos, path_next, next_pos);
        // as a surjection this node must be in the path
        // the nearest place we map should be the head of this node
        assert(path_prev.size()==1 && hit_id == path_prev.front());
        // we then add the position of this node against the path to our offset in the node
        // to get the final chrom+position for the surjection
        path_pos = prev_pos + pos;
        // we need the cigar, but this comes from a function on the alignment itself
        return true;
    } else {
        return false;
    }
}

map<string, int64_t> Index::paths_by_id(void) {
    map<string, int64_t> byid;
    string start = key_for_metadata(path_id_prefix(0));
    start = start.substr(0, start.size()-sizeof(int64_t));
    string end = start + end_sep;
    for_range(start, end, [this, &byid](string& key, string& value) {
            int64_t& id = byid[value];
            memcpy(&id, (void*)(key.c_str() + 10*sizeof(char)), sizeof(int64_t));
        });
    return byid;
}

int64_t Index::path_first_node(int64_t path_id) {
    string k = key_for_path_position(path_id, 0, 0);
    k = k.substr(0, 4 + sizeof(int64_t));
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    rocksdb::Slice start = rocksdb::Slice(k);
    rocksdb::Slice end = rocksdb::Slice(k+end_sep);
    int64_t node_id = 0;
    it->Seek(start);
    if (it->Valid()) {
        string key = it->key().ToString();
        string value = it->value().ToString();
        int64_t path_id2, path_pos; Mapping mapping;
        parse_path_position(key, value, path_id2, path_pos, node_id, mapping);
    }
    delete it;
    return node_id;
}

int64_t Index::path_last_node(int64_t path_id, int64_t& path_length) {
    // we aim to seek to the first item in the next path, then step back
    string key_start = key_for_path_position(path_id, 0, 0);
    string key_end = key_for_path_position(path_id+1, 0, 0);
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    //rocksdb::Slice start = rocksdb::Slice(key_start);
    rocksdb::Slice end = rocksdb::Slice(key_end);
    int64_t node_id = 0;
    it->Seek(end);
    it->Prev();
    // horrible hack
    if (!it->Valid()) it->SeekToLast(); // XXXX
    if (it->Valid()) {
        string key = it->key().ToString();
        string value = it->value().ToString();
        int64_t path_id2, path_pos; Mapping mapping;
        parse_path_position(key, value, path_id2, path_pos, node_id, mapping);
        Node node; get_node(node_id, node);
        path_length = path_pos + node.sequence().size();
    }
    delete it;
    return node_id;
}

void Index::path_layout(map<string, pair<int64_t, int64_t> >& layout,
                        map<string, int64_t>& lengths) {
    map<string, int64_t> pbyid = paths_by_id();
    // for each path
    for (auto& p : pbyid) {
        // find the start and end nodes
        int64_t path_length;
        layout[p.first] = make_pair(path_first_node(p.second),
                                    path_last_node(p.second, path_length));
        lengths[p.first] = path_length;
    }
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
        case 'p': {
            int64_t node_id, path_id, path_pos;
            Mapping mapping;
            parse_node_path(it->key().ToString(), it->value().ToString(),
                            node_id, path_id, path_pos, mapping);
            graph.paths.append_mapping(get_path_name(path_id), mapping);
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
        case 'p': {
            int64_t node_id, path_id, path_pos;
            Mapping mapping;
            parse_node_path(key, value,
                            node_id, path_id, path_pos, mapping);
            graph.paths.append_mapping(get_path_name(path_id), mapping);
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

void Index::approx_sizes_of_kmer_matches(const vector<string>& kmers, vector<uint64_t>& sizes) {
    sizes.resize(kmers.size());
    vector<rocksdb::Range> ranges;
    for (auto& kmer : kmers) {
        string start = key_prefix_for_kmer(kmer);
        string end = start + end_sep;
        ranges.push_back(rocksdb::Range(start, end));
    }
    db->GetApproximateSizes(&ranges[0], kmers.size(), &sizes[0]);
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

void Index::get_path(VG& graph, const string& name, int64_t start, int64_t end) {
    // picks up the specified range in the given path
    if (start < 0 && end < 0) {
        start = 0; end = LONG_MAX;
    }
    int64_t path_id = get_path_id(name);
    string key_start = key_for_path_position(path_id, start, 0);
    string key_end = key_for_path_position(path_id, end, 0);
    for_range(key_start, key_end, [this, &graph](string& key, string& data) {
            Mapping mapping;
            int64_t path_id, path_pos, node_id;
            parse_path_position(key, data,
                                path_id, path_pos,
                                node_id, mapping);
            get_context(node_id, graph);
        });
    // scan the path record in the db to find included nodes
    // get these and drop them into the graph
}

void node_path_position(int64_t id, string& path_name, int64_t& position, int64_t& offset) {
    // if we are in the path, trivial
    // if not, run a BFS back to the nearest node in the path
    
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
    delete it;
}

// todo, get range estimated size

void Index::prune_kmers(int max_kb_on_disk) {
    string start = key_prefix_for_kmer("");
    string end = start + end_sep;
    for_range(start, end, [this, max_kb_on_disk](string& key, string& value) {
            string kmer;
            int64_t id;
            int32_t pos;
            parse_kmer(key, value, kmer, id, pos);
            if (approx_size_of_kmer_matches(kmer) > max_kb_on_disk) {
                //cerr << "pruning kmer " << kmer << endl;
                db->Delete(write_options, key);
            }
        });
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
