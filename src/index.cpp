#include "index.hpp"

namespace vg {

using namespace std;

Index::Index(void) {

    start_sep = '\x00';
    end_sep = '\xff';
    write_options = rocksdb::WriteOptions();
    mem_env = false;
    use_snappy = false;
    // We haven't opened the index yet. We don't get false by default on all platforms.
    is_open = false;
    db = nullptr;
    //block_cache_size = 1024 * 1024 * 10; // 10MB
    rng.seed(time(NULL));

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
    options.max_open_files = -1;
    options.compression = rocksdb::kSnappyCompression;
    options.compaction_style = rocksdb::kCompactionStyleLevel;
    // we are unlikely to reach either of these limits
    options.IncreaseParallelism(threads);
    options.max_background_flushes = threads;
    options.max_background_compactions = threads;

    options.num_levels = 2;
    options.target_file_size_base = (long) 1024 * 1024 * 512; // ~512MB (bigger in practice)
    options.write_buffer_size = 1024 * 1024 * 256; // ~256MB

    // doesn't work this way
    rocksdb::BlockBasedTableOptions topt;
    topt.filter_policy.reset(rocksdb::NewBloomFilterPolicy(10, true));
    topt.block_cache = rocksdb::NewLRUCache(512 * 1024 * 1024, 7);
    topt.no_block_cache = true;
    options.table_factory.reset(NewBlockBasedTableFactory(topt));
    options.table_cache_numshardbits = 7;
    options.allow_mmap_reads = true;
    options.allow_mmap_writes = false;

    if (bulk_load) {
        options.PrepareForBulkLoad();
        options.max_write_buffer_number = threads;
        options.max_background_flushes = threads;
        options.max_background_compactions = threads;
        options.compaction_style = rocksdb::kCompactionStyleNone;
        options.memtable_factory.reset(new rocksdb::VectorRepFactory(1000));
    }

    options.compression_per_level.resize(options.num_levels);
    for (int i = 0; i < options.num_levels; ++i) {
        if (i == 0 || use_snappy == true) {
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
        throw indexOpenException("can't open " + dir);
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

const string Index::key_for_edge_on_start(int64_t node, int64_t other, bool backward) {
    // reverse endianness for sorting
    node = htobe64(node);
    other = htobe64(other);
    string key;
    key.resize(8*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[4 + sizeof(int64_t)] = 's'; // edge on start
    k[5 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*6 + sizeof(int64_t), &other, sizeof(int64_t));
    k[6 + 2*sizeof(int64_t)] = start_sep;
    k[7 + 2*sizeof(int64_t)] = backward ? '1' : '0';
    return key;
}

const string Index::key_for_edge_on_end(int64_t node, int64_t other, bool backward) {
    // reverse endianness for sorting
    node = htobe64(node);
    other = htobe64(other);
    string key;
    key.resize(8*sizeof(char) + 2*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'g'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    k[4 + sizeof(int64_t)] = 'e'; // edge on end
    k[5 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*6 + sizeof(int64_t), &other, sizeof(int64_t));
    k[6 + 2*sizeof(int64_t)] = start_sep;
    k[7 + 2*sizeof(int64_t)] = backward ? '1' : '0';
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

const string Index::key_for_node_path_position(int64_t node_id, int64_t path_id, int64_t path_pos, bool backward) {
    node_id = htobe64(node_id);
    path_id = htobe64(path_id);
    path_pos = htobe64(path_pos);
    string key;
    key.resize(9*sizeof(char) + 3*sizeof(int64_t));
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
    k[7 + 3*sizeof(int64_t)] = start_sep;
    k[8 + 3*sizeof(int64_t)] = backward ? '1' : '0';
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

const string Index::key_for_path_position(int64_t path_id, int64_t path_pos, bool backward, int64_t node_id) {
    node_id = htobe64(node_id);
    path_id = htobe64(path_id);
    path_pos = htobe64(path_pos);
    string key;
    key.resize(7*sizeof(char) + 3*sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'p'; // graph elements
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &path_id, sizeof(int64_t));
    k[3 + sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*4 + sizeof(int64_t), &path_pos, sizeof(int64_t));
    k[4 + 2*sizeof(int64_t)] = start_sep;
    k[5 + 2*sizeof(int64_t)] = backward ? '1' : '0';
    k[6 + 2*sizeof(int64_t)] = start_sep;
    memcpy(k + sizeof(char)*7 + 2*sizeof(int64_t), &node_id, sizeof(int64_t));
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

const string Index::key_for_mapping_prefix(int64_t node_id) {
    node_id = htobe64(node_id);
    string key;
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 's'; // mappings (~sides)
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    return key;
}

const string Index::key_for_mapping(const Mapping& mapping) {
    const string prefix = key_for_mapping_prefix(mapping.position().node_id());
    // use first 8 chars of sha1sum of object; space is 16^8 = 4294967296
    uniform_int_distribution<int> dist(0, 1e8);
    stringstream t;
    t << dist(rng);
    string data = t.str();
    //mapping.SerializeToString(&data);
    const string hash = sha1head(data, 8);
    string key = prefix;
    key.resize(prefix.size() + sizeof(char) + hash.size());
    char* k = (char*) key.c_str();
    k[prefix.size()] = start_sep;
    memcpy(k + prefix.size() + sizeof(char), hash.c_str(), hash.size());
    return key;
}

const string Index::key_for_alignment_prefix(int64_t node_id) {
    node_id = htobe64(node_id);
    string key;
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'a'; // alignments
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    return key;
}

const string Index::key_for_alignment(const Alignment& alignment) {
    const string prefix = key_for_alignment_prefix(alignment.path().mapping(0).position().node_id());
    // use first 8 chars of sha1sum of object; space is 16^8 = 4294967296
    // maybe this shouldn't be a hash, but a random nonce
    uniform_int_distribution<int> dist(0, 1e8);
    stringstream t;
    t << dist(rng);
    string data = t.str();
    //alignment.SerializeToString(&data);
    const string hash = sha1head(data, 8);
    string key = prefix;
    key.resize(prefix.size() + sizeof(char) + hash.size());
    char* k = (char*) key.c_str();
    k[prefix.size()] = start_sep;
    memcpy(k + prefix.size() + sizeof(char), hash.c_str(), hash.size());
    return key;
}

const string Index::key_for_base(int64_t aln_id) {
    string key;
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'b'; // base-alignments
    k[2] = start_sep;
    aln_id = htobe64(aln_id);
    memcpy(k + sizeof(char)*3, &aln_id, sizeof(int64_t));
    return key;
}

const string Index::key_prefix_for_traversal(int64_t node_id) {
    string key;
    key.resize(3*sizeof(char) + sizeof(int64_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 't'; // traversals
    k[2] = start_sep;
    node_id = htobe64(node_id);
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    return key;
}

const string Index::key_for_traversal(int64_t aln_id, const Mapping& mapping) {
    int64_t node_id = mapping.position().node_id();
    node_id = htobe64(node_id);
    string key;
    key.resize(3*sizeof(char) + 2*sizeof(int64_t) + sizeof(int16_t));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 't'; // traversals
    k[2] = start_sep;
    memcpy(k + sizeof(char)*3, &node_id, sizeof(int64_t));
    aln_id = htobe64(aln_id);
    memcpy(k + sizeof(char)*3+sizeof(int64_t), &aln_id, sizeof(int64_t));
    int16_t rank = mapping.rank() * (mapping.position().is_reverse() ? -1 : 1);
    memcpy(k + sizeof(char)*3+sizeof(int64_t)*2, &rank, sizeof(int16_t));
    return key;
}

const string Index::key_prefix_for_edges_on_node_start(int64_t node) {
    string key = key_for_edge_on_start(node, 0, false);
    return key.substr(0, key.size()-sizeof(int64_t)-2*sizeof(char));
}

const string Index::key_prefix_for_edges_on_node_end(int64_t node) {
    string key = key_for_edge_on_end(node, 0, false);
    return key.substr(0, key.size()-sizeof(int64_t)-2*sizeof(char));
}

char Index::graph_key_type(const string& key) {
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
    case 's':
        return mapping_entry_to_string(key, value);
        break;
    case 'a':
        return alignment_entry_to_string(key, value);
        break;
    case 'b':
        return base_entry_to_string(key, value);
        break;
    case 't':
        return traversal_entry_to_string(key, value);
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

void Index::parse_edge(const string& key, char& type, int64_t& node_id, int64_t& other_id, bool& backward) {
    // Parse the edge just out of the key
    const char* k = key.c_str();

    // Work out what type the key is ('s' or 'e' depending on if it's on the first node's start or end).
    type = graph_key_type(key);

    // Get the node IDs involved.
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&other_id, (k + 6*sizeof(char)) + sizeof(int64_t), sizeof(int64_t));
    node_id = be64toh(node_id);
    other_id = be64toh(other_id);

    // Is the relative orientation forward ('0') or backward ('1')?
    char backward_char;
    memcpy(&backward_char, (k + 7*sizeof(char)) + 2*sizeof(int64_t), sizeof(char));
    backward = backward_char == '1';
}

void Index::parse_edge(const string& key, const string& value, char& type, int64_t& id1, int64_t& id2, Edge& edge) {
    // We can take either of the two edge keys:
    // +g+node_id+s+other_id+backward
    // +g+node_id+e+other_id+backward


    if(value.size() > 0) {
        // We can just deserialize the edge.
        edge.ParseFromString(value);

        // But we still need to fill in our output parameters
        type = graph_key_type(key);
        id1 = edge.from();
        id2 = edge.to();

    } else {
        // We have to synthesize an edge.

        // Get what we can from the key. Arbitrarily say this node is the from.
        bool backward;
        parse_edge(key, type, id1, id2, backward);

        // Work out if the edge should be from the start
        bool from_start = type == 's';
        // And if it should be to the end. We attach to the end of the other
        // node when we attached to the start of this node and we want to be
        // forward, or when we attached to the end of this node and we want to
        // be backward. That works out to: XOR(on start, should be backward).
        bool to_end = from_start != backward;

        if(from_start && to_end) {
            // If we got that it should be both, we can replace it with the
            // normal end to start edge going the other way.
            swap(id1, id2);
            from_start = to_end = false;
        }

        // Build the edge
        edge.set_from(id1);
        edge.set_to(id2);
        edge.set_from_start(from_start);
        edge.set_to_end(to_end);

        // TODO: get the edge data somehow in these cases instead of making up edges.
    }
}

string Index::graph_entry_to_string(const string& key, const string& value) {
    // do we have a node or edge?
    stringstream s;
    switch (graph_key_type(key)) {
    case 'n': {
        // it's a node
        int64_t id;
        Node node;
        parse_node(key, value, id, node);
        s << "{\"key\":\"+g+" << id << "+n\", \"value\":"<< pb2json(node) << "}";
    } break;
    case 's': {
        Edge edge;
        int64_t id1, id2;
        char type;
        bool backward;
        if(value.size() > 0) {
            edge.ParseFromString(value);
        }
        parse_edge(key, type, id1, id2, backward);
        s << "{\"key\":\"+g+" << id1 << "+s+" << id2 << "+" << (backward ? '1' : '0')
          << "\", \"value\":"<< (value.size() > 0 ? pb2json(edge) : "") << "}";
    } break;
    case 'e': {
        Edge edge;
        int64_t id1, id2;
        char type;
        bool backward;
        if(value.size() > 0) {
            edge.ParseFromString(value);
        }
        parse_edge(key, type, id1, id2, backward);
        s << "{\"key\":\"+g+" << id1 << "+e+" << id2 << "+" << (backward ? '1' : '0')
          << "\", \"value\":"<< (value.size() > 0 ? pb2json(edge) : "") << "}";
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
                            int64_t& node_id, int64_t& path_id, int64_t& path_pos, bool& backward, Mapping& mapping) {
    const char* k = key.c_str();
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&path_id, (k + 6*sizeof(char)+sizeof(int64_t)), sizeof(int64_t));
    memcpy(&path_pos, (k + 7*sizeof(char)+2*sizeof(int64_t)), sizeof(int64_t));
    backward = (k[8 + 3*sizeof(int64_t)] == '1');
    node_id = be64toh(node_id);
    path_id = be64toh(path_id);
    path_pos = be64toh(path_pos);
    mapping.ParseFromString(value);
}

void Index::parse_path_position(const string& key, const string& value,
                                int64_t& path_id, int64_t& path_pos, bool& backward, int64_t& node_id, Mapping& mapping) {
    const char* k = key.c_str();
    memcpy(&path_id, (k + 3*sizeof(char)), sizeof(int64_t));
    memcpy(&path_pos, (k + 4*sizeof(char)+sizeof(int64_t)), sizeof(int64_t));
    backward = (k[5 + 2*sizeof(int64_t)] == '1');
    memcpy(&node_id, (k + 7*sizeof(char)+2*sizeof(int64_t)), sizeof(int64_t));
    node_id = be64toh(node_id);
    path_id = be64toh(path_id);
    path_pos = be64toh(path_pos);
    mapping.ParseFromString(value);
}

void Index::parse_mapping(const string& key, const string& value, int64_t& node_id, string& hash, Mapping& mapping) {
    const char* k = key.c_str();
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    hash.resize(8);
    memcpy((char*)hash.c_str(), (k + 4*sizeof(char) + sizeof(int64_t)), 8*sizeof(char));
    node_id = be64toh(node_id);
    mapping.ParseFromString(value);
}

void Index::parse_alignment(const string& key, const string& value, int64_t& node_id, string& hash, Alignment& alignment) {
    const char* k = key.c_str();
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    hash.resize(8);
    memcpy((char*)hash.c_str(), (k + 4*sizeof(char) + sizeof(int64_t)), 8*sizeof(char));
    node_id = be64toh(node_id);
    alignment.ParseFromString(value);
}

void Index::parse_base(const string& key, const string& value, int64_t& aln_id, Alignment& alignment) {
    const char* k = key.c_str();
    memcpy(&aln_id, (k + 3*sizeof(char)), sizeof(int64_t));
    aln_id = be64toh(aln_id);
    alignment.ParseFromString(value);
}

void Index::parse_traversal(const string& key, const string& value, int64_t& node_id, int16_t& rank, bool& backward, int64_t& aln_id) {
    const char* k = key.c_str();
    memcpy(&node_id, (k + 3*sizeof(char)), sizeof(int64_t));
    node_id = be64toh(node_id);
    memcpy(&aln_id, (k + 3*sizeof(char)+sizeof(int64_t)), sizeof(int64_t));
    aln_id = be64toh(aln_id);
    memcpy(&rank, (k + 3*sizeof(char) + 2*sizeof(int64_t)), sizeof(int16_t));
    if (rank < 0) { backward = true; } else { backward = false; }
    rank = abs(rank);
}

string Index::node_path_to_string(const string& key, const string& value) {
    Mapping mapping;
    int64_t node_id, path_id, path_pos;
    bool backward;
    parse_node_path(key, value, node_id, path_id, path_pos, backward, mapping);
    stringstream s;
    s << "{\"key\":\"+g+" << node_id << "+p+" << path_id << "+" << path_pos << "+" << (backward ? '1' : '0')
      << "\", \"value\":"<< pb2json(mapping) << "}";
    return s.str();
}

string Index::path_position_to_string(const string& key, const string& value) {
    Mapping mapping;
    int64_t node_id, path_id, path_pos;
    bool backward;
    parse_path_position(key, value, path_id, path_pos, backward, node_id, mapping);
    stringstream s;
    s << "{\"key\":\"+p+" << path_id << "+" << path_pos << "+" << (backward ? '1' : '0') << "+" << node_id
      << "\", \"value\":"<< pb2json(mapping) << "}";
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

string Index::mapping_entry_to_string(const string& key, const string& value) {
    Mapping mapping;
    int64_t node_id;
    string hash;
    parse_mapping(key, value, node_id, hash, mapping);
    stringstream s;
    s << "{\"key\":\"+s+" << node_id << "+" << hash << "\", \"value\":"<< pb2json(mapping) << "}";
    return s.str();
}

string Index::alignment_entry_to_string(const string& key, const string& value) {
    Alignment alignment;
    int64_t node_id;
    string hash;
    parse_alignment(key, value, node_id, hash, alignment);
    stringstream s;
    s << "{\"key\":\"+a+" << node_id << "+" << hash << "\", \"value\":"<< pb2json(alignment) << "}";
    return s.str();
}

string Index::base_entry_to_string(const string& key, const string& value) {
    Alignment alignment;
    int64_t aln_id;
    parse_base(key, value, aln_id, alignment);
    stringstream s;
    s << "{\"key\":\"+b+" << aln_id << "\", \"value\":"<< pb2json(alignment) << "}";
    return s.str();
}

string Index::traversal_entry_to_string(const string& key, const string& value) {
    int64_t node_id;
    bool backward;
    int64_t aln_id;
    int16_t rank;
    parse_traversal(key, value, node_id, rank, backward, aln_id);
    stringstream s;
    s << "{\"key\":\"+t+" << node_id << (backward?"r":"f") << rank << "+" << aln_id << "\", \"value\":\"\"}";
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
    // At least one edge key will hold the serialized edge data
    string data;
    edge->SerializeToString(&data);

    // One will probably hold an empty string, unless this is a self loop somehow.
    string null_data;

    // only store serialized edge in the key linking the edge to the smaller
    // node. If the two node IDs are equal, store in both keys (which might just actually be one key).
    string& from_data = (edge->from() <= edge->to()) ? data : null_data;
    string& to_data = (edge->to() <= edge->from()) ? data : null_data;

    // Is the edge reversing relative node orientation?
    bool backward = (edge->from_start() != edge->to_end());

    if(edge->from_start()) {
        // On the from node, we're on the start
        db->Put(write_options, key_for_edge_on_start(edge->from(), edge->to(), backward), from_data);
    } else {
        // On the from node, we're on the end
        db->Put(write_options, key_for_edge_on_end(edge->from(), edge->to(), backward), from_data);
    }

    if(edge->to_end()) {
        // On the to node, we're on the end
        db->Put(write_options, key_for_edge_on_end(edge->to(), edge->from(), backward), to_data);
    } else {
        // On the to node, we're on the start
        db->Put(write_options, key_for_edge_on_start(edge->to(), edge->from(), backward), to_data);
    }
}

void Index::batch_edge(const Edge* edge, rocksdb::WriteBatch& batch) {
    // At least one edge key will hold the serialized edge data
    string data;
    edge->SerializeToString(&data);

    // One will probably hold an empty string, unless this is a self loop somehow.
    string null_data;

    // only store serialized edge in the key linking the edge to the smaller
    // node. If the two node IDs are equal, store in both keys (which might just actually be one key).
    string& from_data = (edge->from() <= edge->to()) ? data : null_data;
    string& to_data = (edge->to() <= edge->from()) ? data : null_data;

    // Is the edge reversing relative node orientation?
    bool backward = (edge->from_start() != edge->to_end());

    if(edge->from_start()) {
        // On the from node, we're on the start
        batch.Put(key_for_edge_on_start(edge->from(), edge->to(), backward), from_data);
    } else {
        // On the from node, we're on the end
        batch.Put(key_for_edge_on_end(edge->from(), edge->to(), backward), from_data);
    }

    if(edge->to_end()) {
        // On the to node, we're on the end
        batch.Put(key_for_edge_on_end(edge->to(), edge->from(), backward), to_data);
    } else {
        // On the to node, we're on the start
        batch.Put(key_for_edge_on_start(edge->to(), edge->from(), backward), to_data);
    }
}

void Index::put_metadata(const string& tag, const string& data) {
    string key = key_for_metadata(tag);
    db->Put(write_options, key, data);
}

void Index::put_node_path(int64_t node_id, int64_t path_id, int64_t path_pos, bool backward, const Mapping& mapping) {
    string data;
    mapping.SerializeToString(&data);
    db->Put(write_options, key_for_node_path_position(node_id, path_id, path_pos, backward), data);
}

void Index::put_path_position(int64_t path_id, int64_t path_pos, bool backward, int64_t node_id, const Mapping& mapping) {
    string data;
    mapping.SerializeToString(&data);
    db->Put(write_options, key_for_path_position(path_id, path_pos, backward, node_id), data);
}

void Index::put_mapping(const Mapping& mapping) {
    string data;
    mapping.SerializeToString(&data);
    db->Put(write_options, key_for_mapping(mapping), data);
}

void Index::put_alignment(const Alignment& alignment) {
    string data;
    alignment.SerializeToString(&data);
    db->Put(write_options, key_for_alignment(alignment), data);
}

void Index::put_base(int64_t aln_id, const Alignment& alignment) {
    string data;
    alignment.SerializeToString(&data);
    db->Put(write_options, key_for_base(aln_id), data);
}

void Index::put_traversal(int64_t aln_id, const Mapping& mapping) {
    string data; // empty data
    db->Put(write_options, key_for_traversal(aln_id, mapping), data);
}

void Index::cross_alignment(int64_t aln_id, const Alignment& alignment) {
    put_base(aln_id, alignment);
    if (alignment.has_path()) {
        auto& path = alignment.path();
        for (int i = 0; i < path.mapping_size(); ++i) {
            put_traversal(aln_id, path.mapping(i));
        }
    }
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
    graph.preload_progress("indexing nodes of " + graph.name);
    rocksdb::WriteBatch batch;
    graph.for_each_node_parallel([this, &batch](Node* n) { batch_node(n, batch); });
    graph.preload_progress("indexing edges of " + graph.name);
    graph.for_each_edge_parallel([this, &batch](Edge* e) { batch_edge(e, batch); });
    rocksdb::Status s = db->Write(write_options, &batch);
    omp_set_num_threads(thread_count);
}

void Index::load_paths(VG& graph) {
    graph.create_progress("indexing paths of " + graph.name, graph.paths._paths.size());
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
    // TODO: reraise errors other than NotFound...
    if (get_metadata(path_id_prefix(id), data).ok()) {
        return data;
    }
    return string();
}

int64_t Index::get_path_id(const string& name) {
    string data;
    int64_t id = 0;
    // TODO: reraise errors other than NotFound...
    if (get_metadata(path_name_prefix(name), data).ok()) {
        memcpy(&id, (char*)data.c_str(), sizeof(int64_t));
    }
    return id;
}

void Index::store_paths(VG& graph) {
    function<void(const Path&)> lambda = [this, &graph](const Path& path) {
        store_path(graph, path);
    };
    graph.paths.for_each(lambda);
}

void Index::store_path(VG& graph, const Path& path) {
    // get a new path id
    // if there is no name, cry
    if (path.name().empty()) {
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
        put_path_position(path_id, path_pos, mapping.position().is_reverse(), mapping.position().node_id(), mapping);
        // put an entry in the graph table
        put_node_path(mapping.position().node_id(), path_id, path_pos, mapping.position().is_reverse(), mapping);

        // get the node, to find the size of this step
        Node node;
        get_node(mapping.position().node_id(), node);
        // TODO use the cigar... if there is one
        path_pos += node.sequence().size();

        graph.increment_progress();
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

rocksdb::Status Index::get_edge(int64_t from, bool from_start, int64_t to, bool to_end, Edge& edge) {
    // Are we looking for a reversing edge?
    bool backward = from_start != to_end;

    // What key do we need to look up to get the edge data?
    string key;

    // TODO: restructure keys so we don't need to do so much figuring to work out what to look up.
    if(from < to) {
        // We will find the edge data on the record for its attachment to the from node.
        if(from_start) {
            key = key_for_edge_on_start(from, to, backward);
        } else {
            key = key_for_edge_on_end(from, to, backward);
        }
    } else {
        // We will find the edge data on the record for its attachment to the to node.
        if(to_end) {
            key = key_for_edge_on_end(to, from, backward);
        } else {
            key = key_for_edge_on_start(to, from, backward);
        }
    }

    string value;
    rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key, &value);
    if (s.ok()) {
        edge.ParseFromString(value);
    }
    return s;
}

void Index::get_mappings(int64_t node_id, vector<Mapping>& mappings) {
    string start = key_for_mapping_prefix(node_id);
    string end = start + end_sep;
    for_range(start, end, [this, &mappings](string& key, string& value) {
            mappings.emplace_back();
            Mapping& mapping = mappings.back();
            mapping.ParseFromString(value);
        });
}

void Index::get_alignments(int64_t node_id, vector<Alignment>& alignments) {
    string start = key_for_alignment_prefix(node_id);
    string end = start + end_sep;
    for_range(start, end, [this, &alignments](string& key, string& value) {
            alignments.emplace_back();
            Alignment& alignment = alignments.back();
            alignment.ParseFromString(value);
        });
}

void Index::get_alignments(int64_t id1, int64_t id2, vector<Alignment>& alignments) {
    string start = key_for_alignment_prefix(id1);
    string end = key_for_alignment_prefix(id2) + end_sep;
    for_range(start, end, [this, &alignments](string& key, string& value) {
            alignments.emplace_back();
            Alignment& alignment = alignments.back();
            alignment.ParseFromString(value);
        });
}

void Index::for_alignment_in_range(int64_t id1, int64_t id2, std::function<void(const Alignment&)> lambda) {
    string start = key_for_alignment_prefix(id1);
    string end = key_for_alignment_prefix(id2) + end_sep;
    for_range(start, end, [this, &lambda](string& key, string& value) {
            Alignment alignment;
            alignment.ParseFromString(value);
            lambda(alignment);
        });
}

void Index::for_alignment_to_nodes(const vector<int64_t>& ids, std::function<void(const Alignment&)> lambda) {
    set<int64_t> aln_ids;
    for (auto id : ids) {
        string start = key_prefix_for_traversal(id);
        string end = start + end_sep;
        for_range(start, end, [this, &lambda, &aln_ids](string& key, string& value) {
                // parse the alignment id out
                int64_t node_id;
                int16_t rank;
                bool backward;
                int64_t aln_id;
                parse_traversal(key, value, node_id, rank, backward, aln_id);
                aln_ids.insert(aln_id);
            });
    }
    for_base_alignments(aln_ids, lambda);
}

void Index::for_base_alignments(const set<int64_t>& aln_ids, std::function<void(const Alignment&)> lambda) {
    for (auto id : aln_ids) {
        string start = key_for_base(id);
        string end = start + end_sep;
        for_range(start, end, [this, &lambda](string& key, string& value) {
                Alignment alignment;
                int64_t aln_id;
                parse_base(key, value, aln_id, alignment);
                lambda(alignment);
            });
    }
}

int Index::get_node_path(int64_t node_id, int64_t path_id, int64_t& path_pos, bool& backward, Mapping& mapping) {
    string value;
    string key = key_prefix_for_node_path(node_id, path_id);
    string start = key + start_sep;
    string end = key + end_sep;
    // NB: uses the first position in the range
    // apply to the range matching the kmer in the db
    int count = 0;
    for_range(start, end, [this, &count, &node_id, &path_id, &path_pos, &backward, &mapping](string& key, string& value) {
            if (count == 0) {
                parse_node_path(key, value,
                                node_id, path_id,
                                path_pos, backward, mapping);
            }
            ++count;
        });
    return count;
}

pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> Index::get_nearest_node_prev_path_member(
    int64_t node_id, bool backward, int64_t path_id, int64_t& path_pos, bool& relative_orientation, int max_steps) {

    list<pair<int64_t, bool>> nullpath;
    list<pair<int64_t, bool>> bpath;

    // Keeps a list of oriented nodes we can reach, by path taken to reach them
    map<list<pair<int64_t, bool>>, pair<Node, bool>> nq;

    { // handle this node
        // Put this node on the path
        bpath.push_front(make_pair(node_id, backward));
        // Load the node
        Node& node = nq[bpath].first;
        get_node(node_id, node);
        nq[bpath].second = backward;


        Mapping mapping;
        bool backward_on_path;
        if (get_node_path(node_id, path_id, path_pos, backward_on_path, mapping) > 0) {
            // This node is on the target path.

            // Report if we were looking at it backward relative to the path
            relative_orientation = (backward != backward_on_path);

            // Return a search path of just this node, and its ID. We inclue it
            // in the search path (due to being the end of the search) even
            // though it wouldn't normally be included (due to being the start
            // of the search).
            return make_pair(bpath, bpath.front());
        }

        // Otherwise, say we're at this node after taking the empty path.
        Node n = node;
        nq.clear();
        nq[nullpath] = make_pair(n, backward);
    }

    // BFS back
    int steps_back = 0;
    while (steps_back++ < max_steps) {
        // We're going to extend all the paths and populate this.
        map<list<pair<int64_t, bool>>, pair<Node, bool>> cq;
        for (auto& n : nq) {
            // Unpack the entry
            Node& node = n.second.first;
            bool orientation = n.second.second;
            const list<pair<int64_t, bool>>& path = n.first;

            // Look off the left side of this oriented node, and get the oriented nodes you find there.
            vector<pair<int64_t, bool>> destinations;
            get_nodes_prev(node.id(), orientation, destinations);

            for(auto& destination : destinations) {
                int64_t id = destination.first;

                // Extend the path on the left with this destination
                list<pair<int64_t, bool>> npath = path;
                npath.push_front(destination);

                // Fill in the Node object and orientation you can reach via this path
                Node& node = cq[npath].first;
                get_node(id, node);
                cq[npath].second = destination.second;

                Mapping mapping;
                bool backward_on_path;
                if (get_node_path(id, path_id, path_pos, backward_on_path, mapping) > 0) {
                    // This node we just reached is on the path

                    // Report if we were looking at it backward relative to the path
                    relative_orientation = (destination.second != backward_on_path);

                    if(!relative_orientation) {
                        // The right side of this oriented node, which we reached, comes later in the path.
                        path_pos += node.sequence().size();
                    }
                    return make_pair(npath, npath.front());
                }
            }
        }
        // Advance to the next search stage.
        nq = cq;
    }

    // If we get here, we failed to find a path
    relative_orientation = false;
    return make_pair(nullpath, make_pair(0, false));
}

pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> Index::get_nearest_node_next_path_member(
    int64_t node_id, bool backward, int64_t path_id, int64_t& path_pos, bool& relative_orientation, int max_steps) {

    list<pair<int64_t, bool>> nullpath;
    list<pair<int64_t, bool>> bpath;

    // Keeps a list of oriented nodes we can reach, by path taken to reach them
    map<list<pair<int64_t, bool>>, pair<Node, bool>> nq;

    { // handle this node
        // Put this node on the path
        bpath.push_back(make_pair(node_id, backward));
        // Load the node
        Node& node = nq[bpath].first;
        get_node(node_id, node);
        nq[bpath].second = backward;


        Mapping mapping;
        bool backward_on_path;
        if (get_node_path(node_id, path_id, path_pos, backward_on_path, mapping) > 0) {
            // This node is on the target path.

            // Report if we were looking at it backward relative to the path
            relative_orientation = (backward != backward_on_path);

            // Return a search path of just this node, and its ID. We inclue it
            // in the search path (due to being the end of the search) even
            // though it wouldn't normally be included (due to being the start
            // of the search).
            return make_pair(bpath, bpath.back());
        }

        // Otherwise, say we're at this node after taking the empty path.
        Node n = node;
        nq.clear();
        nq[nullpath] = make_pair(n, backward);
    }

    // BFS forward
    int steps_forward = 0;
    while (steps_forward++ < max_steps) {
        // We're going to extend all the paths and populate this.
        map<list<pair<int64_t, bool>>, pair<Node, bool>> cq;
        for (auto& n : nq) {
            // Unpack the entry
            Node& node = n.second.first;
            bool orientation = n.second.second;
            const list<pair<int64_t, bool>>& path = n.first;

            // Look off the right side of this oriented node, and get the oriented nodes you find there.
            vector<pair<int64_t, bool>> destinations;
            get_nodes_next(node.id(), orientation, destinations);

            for(auto& destination : destinations) {
                int64_t id = destination.first;

                // Extend the path on the left with this destination
                list<pair<int64_t, bool>> npath = path;
                npath.push_back(destination);

                // Fill in the Node object and orientation you can reach via this path
                Node& node = cq[npath].first;
                get_node(id, node);
                cq[npath].second = destination.second;

                Mapping mapping;
                bool backward_on_path;
                if (get_node_path(id, path_id, path_pos, backward_on_path, mapping) > 0) {
                    // This node we just reached is on the path

                    // Report if we were looking at it backward relative to the path
                    relative_orientation = (destination.second != backward_on_path);

                    if(relative_orientation) {
                        // The *left* side of this oriented node, which we reached, comes later in the path.
                        path_pos += node.sequence().size();
                    }
                    return make_pair(npath, npath.back());
                }
            }
        }
        // Advance to the next search stage.
        nq = cq;
    }

    // If we get here, we failed to find a path
    relative_orientation = false;
    return make_pair(nullpath, make_pair(0, false));
}

bool Index::get_node_path_relative_position(int64_t node_id, bool backward, int64_t path_id,
                                            list<pair<int64_t, bool>>& path_prev, int64_t& prev_pos, bool& prev_orientation,
                                            list<pair<int64_t, bool>>& path_next, int64_t& next_pos, bool& next_orientation) {
    // scan the range before the node
    // start with our node, and walk back BFS until we find a node with a path
    // are any parents part of the path?

    list<pair<int64_t, bool>> nullpath;
    auto null_pair = make_pair(nullpath, make_pair((int64_t)0, false));

    auto to_path_prev = get_nearest_node_prev_path_member(node_id, backward, path_id, prev_pos, prev_orientation);
    if (to_path_prev == null_pair) {
        cerr << "no to path" << endl;
        return false;
    } else {
        path_prev = to_path_prev.first;
    }

    auto to_path_next = get_nearest_node_next_path_member(node_id, backward, path_id, next_pos, next_orientation);
    if (to_path_next == null_pair) {
        cerr << "no from path" << endl;
        return false;
    } else {
        path_next = to_path_next.first;
    }

    if(next_orientation != prev_orientation) {
        // TODO: this will only happen if cycles are possible, but it's not clear how to handle it.
        cerr << "meets path in different orientations from different ends" << endl;
        return false;
    }

    return true;
}

Mapping Index::path_relative_mapping(int64_t node_id, bool backward, int64_t path_id,
                                     list<pair<int64_t, bool>>& path_prev, int64_t& prev_pos, bool& prev_orientation,
                                     list<pair<int64_t, bool>>& path_next, int64_t& next_pos, bool& next_orientation) {
    Mapping mapping;
    // TODO: shouldn't this point to the node(s?) we're changing, not the one we changed to?
    // TODO: yes it should, but it does not...
    mapping.mutable_position()->set_node_id(node_id);
    mapping.mutable_position()->set_is_reverse(backward);
    // what about offset?
    // TODO I assume this condition works for now, but I do so with a whole salt lick of salt.
    if (get_node_path_relative_position(node_id, backward, path_id,
                                        path_prev, prev_pos, prev_orientation, path_next, next_pos, next_orientation)) {
        // We found a way to the path.

        Edit* edit = mapping.add_edit();
        Node node; get_node(node_id, node);
        // See if we're actually just on that path.
        bool in_path = path_prev.back().first == node_id && path_next.front().first == node_id;
        int32_t to_length = node.sequence().size();
        // The length we replace is either our length if we're on the path, or
        // the distance between where we meet the path on our right and our left
        // otherwise. We need to account for being backwards relative to path coordinates though.
        int32_t from_length = in_path ? to_length : max(next_pos, prev_pos) - min(next_pos, prev_pos);
        //Case to_len == from_len: we're working with a SNP or an exact match. Kinda the base case.
        if (from_length == to_length) {
            edit->set_from_length(from_length);
            edit->set_to_length(to_length);
            // TODO set sequence
            // TODO path_prev is a lost of <node_id, bool> pairs
            // so accessing the front of it is like accessing the nearest node on the ref path
            int64_t p_id = (path_prev.front()).first;
            Node p; get_node(p_id, p);
            Node n; get_node((path_next.front()).first, n);
            string seq = "";
            // Now that we have the alternate node and its neighbors on the path, get
            // nodes in the path that are across from the alt node
            // (i.e. between the next and prev nodes on path).
            vector<pair<int64_t, bool>> level_nodes_prev;
            get_nodes_next(p.id(), (path_prev.front()).second, level_nodes_prev);
            Mapping m;
            int i;
            for (i = 0; i < level_nodes_prev.size(); i++){
                pair<int64_t, bool> n_id_and_backward = level_nodes_prev[i];
                //TODO I'm suspicious the position argument here is incorrect. We should be using the position
                //of the level node.
                if (get_node_path(n_id_and_backward.first, path_id, prev_pos, n_id_and_backward.second, m) <= 0){
                    Node n_id_node; get_node(n_id_and_backward.first, n_id_node);
                    seq = n_id_node.sequence();
                }
            }
            edit->set_sequence(seq);
        } else {
            edit->set_from_length(from_length);
            edit->set_to_length(to_length);
            //TODO set sequence
            edit->set_sequence(node.sequence());
        }
    }
    return mapping;
}

// transform the path into a path relative to another path (defined by path_name)
// source -> surjection (in path_name coordinate space)
// the product is equivalent to a pairwise alignment between this path and the other

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

pair<int64_t, bool> Index::path_first_node(int64_t path_id) {
    string k = key_for_path_position(path_id, 0, false, 0);
    k = k.substr(0, 4 + sizeof(int64_t));
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    rocksdb::Slice start = rocksdb::Slice(k);
    rocksdb::Slice end = rocksdb::Slice(k+end_sep);
    int64_t node_id = 0;
    bool backward;
    it->Seek(start);
    if (it->Valid()) {
        string key = it->key().ToString();
        string value = it->value().ToString();
        int64_t path_id2, path_pos; Mapping mapping;
        parse_path_position(key, value, path_id2, path_pos, backward, node_id, mapping);
    }
    delete it;
    return make_pair(node_id, backward);
}

pair<int64_t, bool> Index::path_last_node(int64_t path_id, int64_t& path_length) {
    // we aim to seek to the first item in the next path, then step back
    string key_start = key_for_path_position(path_id, 0, false, 0);
    string key_end = key_for_path_position(path_id+1, 0, false, 0);
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    //rocksdb::Slice start = rocksdb::Slice(key_start);
    rocksdb::Slice end = rocksdb::Slice(key_end);
    int64_t node_id = 0;
    bool backward;
    it->Seek(end);
    if (it->Valid()) {
        it->Prev();
    }
    else {
        it->SeekToLast();
    }
    if (it->Valid()) {
        string key = it->key().ToString();
        string value = it->value().ToString();
        int64_t path_id2, path_pos; Mapping mapping;
        parse_path_position(key, value, path_id2, path_pos, backward, node_id, mapping);
        Node node; get_node(node_id, node);
        path_length = path_pos + node.sequence().size();
    }
    delete it;
    return make_pair(node_id, backward);
}

void Index::path_layout(map<string, pair<pair<int64_t, bool>, pair<int64_t, bool>> >& layout,
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

void Index::for_each_alignment(function<void(const Alignment&)> lambda) {
    string key;
    key.resize(2*sizeof(char));
    char* k = (char*) key.c_str();
    k[0] = start_sep;
    k[1] = 'a'; // alignments
    string start = key;
    string end = start + end_sep;
    for_range(start, end, [this, &lambda](string& key, string& value) {
            Alignment alignment;
            alignment.ParseFromString(value);
            lambda(alignment);
        });
}

void Index::for_each_mapping(function<void(const Mapping&)> lambda) {
    string start = start_sep + "s";
    string end = start + end_sep;
    for_range(start, end, [this, &lambda](string& key, string& value) {
            Mapping mapping;
            mapping.ParseFromString(value);
            lambda(mapping);
        });
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
        // TODO: optimize this to only look at newly added edges on subsequent steps.
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
            // Key describes the node
            Node node;
            node.ParseFromString(it->value().ToString());
            graph.add_node(node);
        } break;
        case 's': {
            // Key describes an edge on the start of the node
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(it->key().ToString(), it->value().ToString(), type, id1, id2, edge);
            graph.add_edge(edge);
        } break;
        case 'e': {
            // Key describes an edge on the end of the node
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
            // Key describes a path membership
            int64_t node_id, path_id, path_pos;
            Mapping mapping;
            bool backward;
            parse_node_path(it->key().ToString(), it->value().ToString(),
                            node_id, path_id, path_pos, backward, mapping);
            // We don't need to pass backward here since it's included in the Mapping object.
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
            // Key describes a node
            Node node;
            node.ParseFromString(value);
            graph.add_node(node);
        } break;
        case 's': {
            // Key describes an edge on the start of a node
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(key, value, type, id1, id2, edge);
            graph.add_edge(edge);
        } break;
        case 'e': {
            // Key describes an edge on the end of a node
            Edge edge;
            int64_t id1, id2;
            char type;
            parse_edge(key, value, type, id1, id2, edge);
            graph.add_edge(edge);
        } break;
        case 'p': {
            // Key describes a path membership
            int64_t node_id, path_id, path_pos;
            Mapping mapping;
            bool backward;
            parse_node_path(key, value,
                            node_id, path_id, path_pos, backward, mapping);
            // We don't need to pass backward here since it's included in the Mapping object.
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
    // We can't rely on edge keys coming after their node keys, so we need to
    // trim off the trailing "+n" from the first key, so we get all the edges.
    string start = key_for_node(from_id).substr(0,3+sizeof(int64_t));
    // Similarly, we need to stop before the edges attached to this other node,
    // if there are any.
    string end = key_for_node(to_id+1).substr(0,3+sizeof(int64_t));
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

void Index::get_edges_on_start(int64_t node_id, vector<Edge>& edges) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    string key_start = key_prefix_for_edges_on_node_start(node_id);
    rocksdb::Slice start = rocksdb::Slice(key_start);
    string key_end = key_start+end_sep;
    rocksdb::Slice end = rocksdb::Slice(key_end);
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string s = it->key().ToString();
        char keyt = graph_key_type(s);
        switch (keyt) {
        case 's': {
            // Parse out the edge
            Edge edge;
            int64_t id1, id2;
            char type;
            bool backward;
            parse_edge(it->key().ToString(), type, id1, id2, backward);

            // TODO: If we can know we don't really need the edge metadata, we could stop here and save a lookup.

            // What's the other node involved in this edge?
            assert(node_id == id1);
            int64_t other_id = id2;

            if(other_id < node_id) {
                // The edge metadata wasn't stored here. We need to look it up.

                // What's the key for the other end of this edge? If this edge
                // is reversing, then it's on the other node's start, too.
                // Otherwise it's on the other node's end.
                string other_key = backward ?
                    key_for_edge_on_start(other_id, node_id, backward) :
                    key_for_edge_on_end(other_id, node_id, backward);

                // Load up that key
                string value;
                rocksdb::Status status = db->Get(rocksdb::ReadOptions(), other_key, &value);
                if (status.ok()) {
                    edge.ParseFromString(value);
                } else {
                    cerr << entry_to_string(s, "") << " looking for " << entry_to_string(other_key, "") << endl;
                    throw std::runtime_error("Could not find other end of edge on start");
                }
            }
            edges.push_back(edge);
        } break;
        default:
            // there should only be edges on the start
            cerr << keyt << endl;
            assert(false);
            break;
        }
    }
}

void Index::get_edges_on_end(int64_t node_id, vector<Edge>& edges) {
    rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    string key_start = key_prefix_for_edges_on_node_end(node_id);
    rocksdb::Slice start = rocksdb::Slice(key_start);
    string key_end = key_start+end_sep;
    rocksdb::Slice end = rocksdb::Slice(key_end);
    for (it->Seek(start);
         it->Valid() && it->key().ToString() < key_end;
         it->Next()) {
        string s = it->key().ToString();
        char keyt = graph_key_type(s);
        switch (keyt) {
        case 'e': {
            Edge edge;
            int64_t id1, id2;
            char type;
            bool backward;
            parse_edge(it->key().ToString(), type, id1, id2, backward);
            // TODO: If we can know we don't really need the edge metadata, we could stop here and save a lookup.

            // What's the other node involved in this edge?
            assert(node_id == id1);
            int64_t other_id = id2;

            if(other_id < node_id) {
                // The edge metadata wasn't stored here. We need to look it up.

                // What's the key for the other end of this edge? If this edge
                // is reversing, then it's on the other node's end, too.
                // Otherwise it's on the other node's start.
                string other_key = backward ?
                    key_for_edge_on_end(other_id, node_id, backward) :
                    key_for_edge_on_start(other_id, node_id, backward);

                // Load up that key
                string value;
                rocksdb::Status status = db->Get(rocksdb::ReadOptions(), other_key, &value);
                if (status.ok()) {
                    edge.ParseFromString(value);
                } else {
                    cerr << entry_to_string(s, "") << " looking for " << entry_to_string(other_key, "") << endl;
                    throw std::runtime_error("Could not find other end of edge on end");
                }
            } else {
                // We have the whole Edge right here
                edge.ParseFromString(it->value().ToString());
            }
            edges.push_back(edge);
        } break;
        default:
            // there should only be edges on the end
            cerr << keyt << endl;
            assert(false);
            break;
        }
    }
}

void Index::get_nodes_next(int64_t node, bool backward, vector<pair<int64_t, bool>>& destinations) {

    // Get all the edges off the appropriate side of the node.
    vector<Edge> edges_to_follow;
    if(backward) {
        // "next" = right = start
        get_edges_on_start(node, edges_to_follow);
    } else {
        // "next" = right = end
        get_edges_on_end(node, edges_to_follow);
    }

    for(Edge& e : edges_to_follow) {
        // Get the other node involved in the edge
        int64_t other_node = (e.to() == node ? e.from() : e.to());

        // Work out if this is a reversing edge
        bool reversing_edge = e.from_start() != e.to_end();

        // Put in the other node ID and the relative orientation, which is our
        // orientation, only reversed if we crossed a reversing edge.
        destinations.emplace_back(other_node, backward != reversing_edge);
    }
}

void Index::get_nodes_prev(int64_t node, bool backward, vector<pair<int64_t, bool>>& destinations) {
    // TODO: combine with get_nodes_next, since they're basically the same code.

    // Get all the edges off the appropriate side of the node.
    vector<Edge> edges_to_follow;
    if(backward) {
        // "prev" = left = end
        get_edges_on_end(node, edges_to_follow);
    } else {
        // "prev" = left = start
        get_edges_on_start(node, edges_to_follow);
    }

    for(Edge& e : edges_to_follow) {
        // Get the other node involved in the edge
        int64_t other_node = (e.to() == node ? e.from() : e.to());

        // Work out if this is a reversing edge
        bool reversing_edge = e.from_start() != e.to_end();

        // Put in the other node ID and the relative orientation, which is our
        // orientation, only reversed if we crossed a reversing edge.
        destinations.emplace_back(other_node, backward != reversing_edge);
    }

}

void Index::get_path(VG& graph, const string& name, int64_t start, int64_t end) {
    // picks up the specified range in the given path
    if (start < 0 && end < 0) {
        start = 0; end = LONG_MAX;
    }
    int64_t path_id = get_path_id(name);
    string key_start = key_for_path_position(path_id, start, false, 0);
    // This is deliberately before any key we would get for the actual end, because the end is exclusive.
    string key_end = key_for_path_position(path_id, end, false, 0);

    for_range(key_start, key_end, [this, &graph](string& key, string& data) {
            Mapping mapping;
            int64_t path_id, path_pos, node_id;
            bool backward;
            parse_path_position(key, data,
                                path_id, path_pos, backward,
                                node_id, mapping);
            get_context(node_id, graph);
        });
    // scan the path record in the db to find included nodes
    // get these and drop them into the graph
}


void node_path_position(int64_t id, string& path_name, int64_t& position, bool backward, int64_t& offset) {
    // if we are in the path, trivial
    // if not, run a BFS back to the nearest node in the path
    // if (get_node_path_relative_position()){
    //   Node n; // find the path start.
    //   // iterate over nodes in path until you arrive at this one,
    //   // and sum all sequence lengths from said traversal.
    // }
    // else{
    //   // Not on the path, so get previous member.
    //   get_nearest_node_prev_path_member();
    //   // iterate over nodes in path until you arrive at this one,
    //   // and sum all sequence lengths from said traversal.
    //   // add this node's sequence as well IE add offset
    // }

    throw runtime_error("node_path_position not yet implemented");

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
    // TODO: support orientation here.
}

string Index::first_kmer_key(const string& kmer) {
    string found_key;
    function<void(string&, string&)> lambda = [&found_key](string& key, string& value) {
        if (found_key.empty()) {
            found_key = key;
        }
    };
    string first_key = key_for_kmer(kmer, 0);
    string last_key = key_for_kmer(kmer, numeric_limits<int64_t>::max());
    for_range(first_key, last_key, lambda);
    return found_key;
}

pair<int64_t, int64_t> Index::compare_kmers(Index& other) {
    int64_t outFound = 0;
    int64_t outNotFound = 0;
    string prev_kmer;

    function<void(string&, string&)> lambda = [&](string& key, string& value) {
        if (key[1] == 'k') {
            int64_t id;
            int32_t pos;
            string kmer;
            parse_kmer(key, value, kmer, id, pos);

            // only visit first kmer when multiple occurances with dif. ids in a row
            if (kmer != prev_kmer) {

                string remk = reverse_complement(kmer);
                string remk_key = first_kmer_key(remk);

                // only visit canonical strand (ie lexicographic less than reverse comp)
                if (remk_key.empty() || key < remk_key) {

                    // put together a key range that will find all matches to kmer (i think)
                    string first_key = other.key_for_kmer(kmer, 0);
                    string last_key = other.key_for_kmer(kmer, numeric_limits<int64_t>::max());

                    // search other index
                    bool found = false;
                    function<void(string&, string&)> lambda1 = [&found](string& key, string& value) {
                        found = true;
                    };
                    other.for_range(first_key, last_key, lambda1);

                    // wasn't found in other, try reverse complement
                    if (!found) {
                        first_key = other.key_for_kmer(remk, 0);
                        last_key = other.key_for_kmer(remk,  numeric_limits<int64_t>::max());
                        other.for_range(first_key, last_key, lambda1);
                    }

                    // update stats
                    if (found) {
                        ++outFound;
                    } else {
                        ++outNotFound;
                    }
                }
            }
            swap(kmer, prev_kmer);
        }
    };

    // skip things that aren't kmers
    string first_search_key(1, start_sep);
    first_search_key += 'k';
    string last_search_key(1, start_sep);
    last_search_key += 'l';

    for_range(first_search_key, last_search_key, lambda);

    return pair<int64_t, int64_t>(outFound, outNotFound);
}

}
