#include "xg.hpp"
#include "stream.hpp"

#include <bitset>
#include <arpa/inet.h>

//#define VERBOSE_DEBUG
//#define debug_algorithms

namespace xg {

id_t side_id(const side_t& side) {
    return abs(side);
}

bool side_is_end(const side_t& side) {
    return side < 0;
}

side_t make_side(id_t id, bool is_end) {
    return !is_end ? id : -1 * id;
}

id_t trav_id(const trav_t& trav) {
    return abs(trav.first);
}

bool trav_is_rev(const trav_t& trav) {
    return trav.first < 0;
}

int32_t trav_rank(const trav_t& trav) {
    return trav.second;
}

trav_t make_trav(id_t id, bool is_rev, int32_t rank) {
    return make_pair(!is_rev ? id : -1 * id, rank);
}

int dna3bit(char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    default:
        return 4;
    }
}

char revdna3bit(int i) {
    switch (i) {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    default:
        return 'N';
    }
}

const XG::destination_t XG::BS_SEPARATOR = 1;
const XG::destination_t XG::BS_NULL = 0;

XG::XG(istream& in)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    load(in);
}

XG::XG(Graph& graph)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    from_graph(graph);
}

XG::XG(function<void(function<void(Graph&)>)> get_chunks)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    from_callback(get_chunks);
}

XG::~XG(void) {
    // Clean up any created XGPaths
    while (!paths.empty()) {
        delete paths.back();
        paths.pop_back();
    }
}

void XG::load(istream& in) {

    if (!in.good()) {
        throw XGFormatError("Index file does not exist or index stream cannot be read");
    }

    // Version 0 is the last XG format without an explicit version specifier.
    // If we find a version specifier we will up this.
    uint32_t file_version = 0;

    // We need to look for the magic value
    char buffer;
    in.get(buffer);
    if (buffer == 'X') {
        in.get(buffer);
        if (buffer == 'G') {
            // We found the magic value!
            
            // Don't put it back, but the next 4 bytes are a version number.
            in.read((char*) &file_version, sizeof(file_version));
            // Make sure to convert from network to host byte order
            file_version = ntohl(file_version);
            
        } else {
            // Put back both characters
            in.unget();
            in.unget();
        }        
    } else {
        // Put back the one character
        in.unget();
    }
    
    if (file_version > MAX_INPUT_VERSION) {
        // This XG file is from the distant future.
        throw XGFormatError("XG index file version " + to_string(file_version) +
            " is too new to interpret (max version = " + to_string(MAX_INPUT_VERSION) + ")");
    }
    
    try {
        
        ////////////////////////////////////////////////////////////////////////
        // DO NOT CHANGE THIS CODE without creating a new XG version:
        // 1. Increment OUTPUT_VERSION to a new integer.
        // 2. Change the serialization code.
        // 3. Add a new case here with new deserialization code.
        // 4. Up MAX_INPUT_VERSION to allow your new version to be read.
        ////////////////////////////////////////////////////////////////////////
        switch (file_version) {
        
        case 5:
            {
                sdsl::read_member(seq_length, in);
                sdsl::read_member(node_count, in);
                sdsl::read_member(edge_count, in);
                sdsl::read_member(path_count, in);
                size_t entity_count = node_count + edge_count;
                //cerr << sequence_length << ", " << node_count << ", " << edge_count << endl;
                sdsl::read_member(min_id, in);
                sdsl::read_member(max_id, in);

                i_iv.load(in);
                r_iv.load(in);

                g_iv.load(in);
                g_bv.load(in);
                g_bv_rank.load(in, &g_bv);
                g_bv_select.load(in, &g_bv);

                s_iv.load(in);
                s_bv.load(in);
                s_bv_rank.load(in, &s_bv);
                s_bv_select.load(in, &s_bv);

                tn_csa.load(in);
                tn_cbv.load(in);
                tn_cbv_rank.load(in, &tn_cbv);
                tn_cbv_select.load(in, &tn_cbv);
                tin_civ.load(in);
                tio_civ.load(in);
                side_thread_wt.load(in);
                
                pn_iv.load(in);
                pn_csa.load(in);
                pn_bv.load(in);
                pn_bv_rank.load(in, &pn_bv);
                pn_bv_select.load(in, &pn_bv);
                pi_iv.load(in);
                sdsl::read_member(path_count, in);
                for (size_t i = 0; i < path_count; ++i) {
                    auto path = new XGPath;
                    path->load(in);
                    paths.push_back(path);
                }
                np_iv.load(in);
                np_bv.load(in);
                np_bv_rank.load(in, &np_bv);
                np_bv_select.load(in, &np_bv);
                
                h_civ.load(in);
                ts_civ.load(in);

                // Load all the B_s arrays for sides.
                // Baking required before serialization.
                deserialize(bs_single_array, in);
            }
            break;
        default:
            throw XGFormatError("Unimplemented XG format version: " + to_string(file_version));
        }
    } catch (const XGFormatError& e) {
        // Pass XGFormatErrors through
        throw e;
    } catch (const bad_alloc& e) {
        // We get std::bad_alloc generally if we try to read arbitrary data as an xg index.
        throw XGFormatError("XG input data not in XG version " + to_string(file_version) + " format (" + e.what() + ")");
    } catch (const exception& e) {
        // Other things will get re-thrown with a hint.
        cerr << "error [xg]: Unexpected error parsing XG data. Is it in version " << file_version << " XG format?" << endl;
        throw e;
    }

}

void XGPath::load(istream& in) {
    nodes.load(in);
    nodes_rank.load(in, &nodes);
    nodes_select.load(in, &nodes);
    ids.load(in);
    directions.load(in);
    ranks.load(in);
    positions.load(in);
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
}

size_t XGPath::serialize(std::ostream& out,
                         sdsl::structure_tree_node* v,
                         std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += nodes.serialize(out, child, "path_membership_" + name);
    written += nodes_rank.serialize(out, child, "path_membership_rank_" + name);
    written += nodes_select.serialize(out, child, "path_membership_select_" + name);
    written += ids.serialize(out, child, "path_node_ids_" + name);
    written += directions.serialize(out, child, "path_node_directions_" + name);
    written += ranks.serialize(out, child, "path_mapping_ranks_" + name);
    written += positions.serialize(out, child, "path_node_offsets_" + name);
    written += offsets.serialize(out, child, "path_node_starts_" + name);
    written += offsets_rank.serialize(out, child, "path_node_starts_rank_" + name);
    written += offsets_select.serialize(out, child, "path_node_starts_select_" + name);
    
    sdsl::structure_tree::add_size(child, written);
    
    return written;
}

XGPath::XGPath(const string& path_name,
               const vector<trav_t>& path,
               size_t node_count,
               XG& graph,
               size_t* unique_member_count_out) {

    // nodes in the path
    bit_vector nodes_bv;
    util::assign(nodes_bv, bit_vector(node_count));
    // node ids, the literal path
    int_vector<> ids_iv;
    util::assign(ids_iv, int_vector<>(path.size()));
    // directions of traversal (typically forward, but we allow backwards
    bit_vector directions_bv;
    util::assign(directions_bv, bit_vector(path.size()));
    // node positions in path
    util::assign(positions, int_vector<>(path.size()));
    // mapping ranks in path
    util::assign(ranks, int_vector<>(path.size()));

    size_t path_off = 0;
    size_t members_off = 0;
    size_t positions_off = 0;
    size_t path_length = 0;

    // determine total length
    for (size_t i = 0; i < path.size(); ++i) {
        auto node_id = trav_id(path[i]);
        path_length += graph.node_length(node_id);
        ids_iv[i] = node_id;
        // we will explode if the node isn't in the graph
    }

    // make the bitvector for path offsets
    util::assign(offsets, bit_vector(path_length));
    set<int64_t> uniq_nodes;
    //cerr << "path " << path_name << " has " << path.size() << endl;
    for (size_t i = 0; i < path.size(); ++i) {
        //cerr << i << endl;
        auto& trav = path[i];
        auto node_id = trav_id(trav);
        bool is_reverse = trav_is_rev(trav);
        //cerr << node_id << endl;
        // record node
        nodes_bv[graph.id_to_rank(node_id)-1] = 1;
        // record direction of passage through node
        directions_bv[i] = is_reverse;
        // and the external rank of the mapping
        ranks[i] = trav_rank(trav);
        // we've seen another entity
        uniq_nodes.insert(node_id);
        // and record node offset in path
        positions[positions_off++] = path_off;
        // record position of node
        offsets[path_off] = 1;
        // and update the offset counter
        path_off += graph.node_length(node_id);
    }
    //cerr << uniq_nodes.size() << " vs " << path.size() << endl;
    if(unique_member_count_out) {
        // set member count as the unique entities that are in the path
        // We don't need it but our caller might
        *unique_member_count_out = uniq_nodes.size();
    }
    // compress path membership vectors
    util::assign(nodes, nodes_bv);
    util::assign(nodes_rank, rrr_vector<>::rank_1_type(&nodes));
    util::assign(nodes_select, rrr_vector<>::select_1_type(&nodes));
    // and traversal information
    util::assign(directions, sd_vector<>(directions_bv));
    // handle entity lookup structure (wavelet tree)
    util::bit_compress(ids_iv);
    construct_im(ids, ids_iv);
    // bit compress the positional offset info
    util::bit_compress(positions);
    // bit compress mapping ranks
    util::bit_compress(ranks);

    // and set up rank/select dictionary on them
    util::assign(offsets_rank, rank_support_v<1>(&offsets));
    util::assign(offsets_select, bit_vector::select_1_type(&offsets));
}

Mapping XGPath::mapping(size_t offset) const {
    // TODO actually store the "real" mapping
    Mapping m;
    // store the starting position and series of edits
    m.mutable_position()->set_node_id(ids[offset]);
    m.mutable_position()->set_is_reverse(directions[offset]);
    m.set_rank(ranks[offset]);
    return m;
}

size_t XG::serialize(ostream& out, sdsl::structure_tree_node* s, std::string name) {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;
    
    // Do the magic number
    out << "XG";
    written += 2;
    
    // And the version number
    int32_t version_buffer = htonl(OUTPUT_VERSION);
    out.write((char*) &version_buffer, sizeof(version_buffer));
    written += sizeof(version_buffer) / sizeof(char);
    
    ////////////////////////////////////////////////////////////////////////
    // DO NOT CHANGE THIS CODE without creating a new XG version:
    // 1. Increment OUTPUT_VERSION to a new integer.
    // 2. Add your new serialization code.
    // 3. Add a new case for your new version to XG::load()
    // 4. Up MAX_INPUT_VERSION to allow your new version to be read.
    ////////////////////////////////////////////////////////////////////////

    written += sdsl::write_member(s_iv.size(), out, child, "sequence_length");
    written += sdsl::write_member(i_iv.size(), out, child, "node_count");
    written += sdsl::write_member(edge_count, out, child, "edge_count");
    written += sdsl::write_member(path_count, out, child, "path_count");
    written += sdsl::write_member(min_id, out, child, "min_id");
    written += sdsl::write_member(max_id, out, child, "max_id");

    written += i_iv.serialize(out, child, "id_rank_vector");
    written += r_iv.serialize(out, child, "rank_id_vector");

    written += g_iv.serialize(out, child, "graph_vector");
    written += g_bv.serialize(out, child, "graph_bit_vector");
    written += g_bv_rank.serialize(out, child, "graph_bit_vector_rank");
    written += g_bv_select.serialize(out, child, "graph_bit_vector_select");
    
    written += s_iv.serialize(out, child, "seq_vector");
    written += s_bv.serialize(out, child, "seq_node_starts");
    written += s_bv_rank.serialize(out, child, "seq_node_starts_rank");
    written += s_bv_select.serialize(out, child, "seq_node_starts_select");

    // save the thread name index
    written += tn_csa.serialize(out, child, "thread_name_csa");
    written += tn_cbv.serialize(out, child, "thread_name_cbv");
    written += tn_cbv_rank.serialize(out, child, "thread_name_cbv_rank");
    written += tn_cbv_select.serialize(out, child, "thread_name_cbv_select");
    written += tin_civ.serialize(out, child, "thread_start_node_civ");
    written += tio_civ.serialize(out, child, "thread_start_offset_civ");
    written += side_thread_wt.serialize(out, child, "side_thread_wt");

    // Treat the paths as their own node
    size_t paths_written = 0;
    auto paths_child = sdsl::structure_tree::add_child(child, "paths", sdsl::util::class_name(*this));

    paths_written += pn_iv.serialize(out, paths_child, "path_names");
    paths_written += pn_csa.serialize(out, paths_child, "path_names_csa");
    paths_written += pn_bv.serialize(out, paths_child, "path_names_starts");
    paths_written += pn_bv_rank.serialize(out, paths_child, "path_names_starts_rank");
    paths_written += pn_bv_select.serialize(out, paths_child, "path_names_starts_select");
    paths_written += pi_iv.serialize(out, paths_child, "path_ids");
    paths_written += sdsl::write_member(paths.size(), out, paths_child, "path_count");    
    for (size_t i = 0; i < paths.size(); i++) {
        XGPath* path = paths[i];
        paths_written += path->serialize(out, paths_child, "path:" + path_name(i + 1));
    }
    
    paths_written += np_iv.serialize(out, paths_child, "node_path_mapping");
    paths_written += np_bv.serialize(out, paths_child, "node_path_mapping_starts");
    paths_written += np_bv_rank.serialize(out, paths_child, "node_path_mapping_starts_rank");
    paths_written += np_bv_select.serialize(out, paths_child, "node_path_mapping_starts_select");
    
    sdsl::structure_tree::add_size(paths_child, paths_written);
    written += paths_written;

    // Treat the threads as their own node.
    // This will mess up any sort of average size stats, but it will also be useful.
    auto threads_child = sdsl::structure_tree::add_child(child, "threads", sdsl::util::class_name(*this));
    
    size_t threads_written = 0;
    threads_written += h_civ.serialize(out, threads_child, "thread_usage_count");
    threads_written += ts_civ.serialize(out, threads_child, "thread_start_count");
    // Stick all the B_s arrays in together. Must be baked.
    threads_written += xg::serialize(bs_single_array, out, threads_child, "bs_single_array");
    
    sdsl::structure_tree::add_size(threads_child, threads_written);
    written += threads_written;

    sdsl::structure_tree::add_size(child, written);
    return written;
    
}

void XG::from_stream(istream& in, bool validate_graph, bool print_graph,
    bool store_threads, bool is_sorted_dag) {

    from_callback([&](function<void(Graph&)> handle_chunk) {
        // TODO: should I be bandying about function references instead of
        // function objects here?
        stream::for_each(in, handle_chunk);
    }, validate_graph, print_graph, store_threads, is_sorted_dag);
}

void XG::from_graph(Graph& graph, bool validate_graph, bool print_graph,
    bool store_threads, bool is_sorted_dag) {

    from_callback([&](function<void(Graph&)> handle_chunk) {
        // There's only one chunk in this case.
        handle_chunk(graph);
    }, validate_graph, print_graph, store_threads, is_sorted_dag);

}

void XG::from_callback(function<void(function<void(Graph&)>)> get_chunks, 
    bool validate_graph, bool print_graph, bool store_threads, bool is_sorted_dag) {

    // temporaries for construction
    vector<pair<id_t, string> > node_label;
    // need to store node sides
    unordered_map<side_t, vector<side_t> > from_to;
    unordered_map<side_t, vector<side_t> > to_from;
    map<string, vector<trav_t> > path_nodes;

    // This takes in graph chunks and adds them into our temporary storage.
    function<void(Graph&)> lambda = [this,
                                     &node_label,
                                     &from_to,
                                     &to_from,
                                     &path_nodes](Graph& graph) {

        for (int64_t i = 0; i < graph.node_size(); ++i) {
            const Node& n = graph.node(i);
            node_label.push_back(make_pair(n.id(), n.sequence()));
        }
        for (int64_t i = 0; i < graph.edge_size(); ++i) {
            // Canonicalize every edge, so only canonical edges are in the index.
            Edge e = canonicalize(graph.edge(i));
            bool new_from = from_to.find(make_side(e.from(), e.from_start())) == from_to.end();
            bool new_edge = true;
            if (!new_from) {
                side_t to_add = make_side(e.to(), e.to_end());
                for (auto& side : from_to[make_side(e.from(), e.from_start())]) {
                    if (side == to_add) { new_edge = false; break; }
                }
            }
            if (new_edge) {
                ++edge_count;
                from_to[make_side(e.from(), e.from_start())].push_back(make_side(e.to(), e.to_end()));
                to_from[make_side(e.to(), e.to_end())].push_back(make_side(e.from(), e.from_start()));
            }
        }

        // Print out all the paths in the graph we are loading
        for (int64_t i = 0; i < graph.path_size(); ++i) {
            const Path& p = graph.path(i);
            const string& name = p.name();
#ifdef VERBOSE_DEBUG
            cerr << "Path " << name << ": ";
#endif
            vector<trav_t>& path = path_nodes[name];
            for (int64_t j = 0; j < p.mapping_size(); ++j) {
                const Mapping& m = p.mapping(j);
                path.push_back(make_trav(m.position().node_id(), m.position().is_reverse(), m.rank()));
#ifdef VERBOSE_DEBUG
                cerr << m.position().node_id() * 2 + m.position().is_reverse() << "; ";
#endif
            }
#ifdef VERBOSE_DEBUG
            cerr << endl;
#endif

        }
    };

    // Get all the chunks via the callback, and have them called back to us.
    // The other end handles figuring out how much to loop.
    get_chunks(lambda);

    // sort the node labels and remove any duplicates
    std::sort(node_label.begin(), node_label.end());
    node_label.erase(std::unique(node_label.begin(), node_label.end()), node_label.end());
    for (auto& p : node_label) {
        ++node_count;
        seq_length += p.second.size();
    }

    path_count = path_nodes.size();
    
    // sort the paths using mapping rank
    // and remove duplicates
    for (auto& p : path_nodes) {
        vector<trav_t>& path = path_nodes[p.first];
        if (!std::is_sorted(path.begin(), path.end(),
                            [](const trav_t& m1, const trav_t& m2) { return trav_rank(m1) < trav_rank(m2); })) {
            cerr << "[xg] warning: path " << p.first << " is not in sorted order by rank" << endl;
            std::sort(path.begin(), path.end(),
                      [](const trav_t& m1, const trav_t& m2) { return trav_rank(m1) < trav_rank(m2); });
        }
        auto last_unique = std::unique(path.begin(), path.end(),
                                       [](const trav_t& m1, const trav_t& m2) {
                                           return trav_rank(m1) == trav_rank(m2);
                                       });
        
        if (last_unique != path.end()) {
            cerr << "[xg] error: path " << p.first << " contains duplicate node ranks" << endl;
            exit(1);
        }
    }

    build(node_label, from_to, to_from, path_nodes, validate_graph, print_graph,
        store_threads, is_sorted_dag);
    
}

void XG::build(vector<pair<id_t, string> >& node_label,
               unordered_map<side_t, vector<side_t> >& from_to,
               unordered_map<side_t, vector<side_t> >& to_from,
               map<string, vector<trav_t> >& path_nodes,
               bool validate_graph,
               bool print_graph,
               bool store_threads,
               bool is_sorted_dag) {

    size_t entity_count = node_count + edge_count;

#ifdef VERBOSE_DEBUG
    cerr << "graph has " << seq_length << "bp in sequence, "
         << node_count << " nodes, "
         << edge_count << " edges, and "
         << path_count << " paths "
         << "for a total of " << entity_count << " entities" << endl;
#endif

    // for mapping of ids to ranks using a vector rather than wavelet tree
    min_id = node_label.begin()->first;
    max_id = node_label.rbegin()->first;
    
    // set up our compressed representation
    util::assign(s_iv, int_vector<>(seq_length, 0, 3));
    util::assign(s_bv, bit_vector(seq_length));
    util::assign(i_iv, int_vector<>(node_count));
    util::assign(r_iv, int_vector<>(max_id-min_id+1)); // note possibly discontiguous
    
    // for each node in the sequence
    // concatenate the labels into the s_iv
#ifdef VERBOSE_DEBUG
    cerr << "storing node labels" << endl;
#endif
    size_t i = 0; // insertion point
    size_t r = 1;
    for (auto& p : node_label) {
        int64_t id = p.first;
        const string& l = p.second;
        s_bv[i] = 1; // record node start
        i_iv[r-1] = id;
        // store ids to rank mapping
        r_iv[id-min_id] = r;
        ++r;
        for (auto c : l) {
            s_iv[i++] = dna3bit(c); // store sequence
        }
    }
    // keep only if we need to validate the graph
    if (!validate_graph) node_label.clear();

    // we have to process all the nodes before we do the edges
    // because we need to ensure full coverage of node space

    util::bit_compress(i_iv);
    util::bit_compress(r_iv);

    /*
#ifdef VERBOSE_DEBUG    
    cerr << "storing forward edges and adjacency table" << endl;
#endif
    size_t f_itr = 0;
    size_t j_itr = 0; // edge adjacency pointer
    for (size_t k = 0; k < node_count; ++k) {
        int64_t f_id = i_iv[k];
        size_t f_rank = k+1;
        f_iv[f_itr] = f_rank;
        f_bv[f_itr] = 1;
        ++f_itr;
        for (auto end : { false, true }) {
            if (from_to.find(make_side(f_id, end)) != from_to.end()) {
                auto t_side_itr = from_to.find(make_side(f_id, end));
                if (t_side_itr != from_to.end()) {
                    for (auto& t_side : t_side_itr->second) {
                        size_t t_rank = id_to_rank(side_id(t_side));
                        // store link
                        f_iv[f_itr] = t_rank;
                        f_bv[f_itr] = 0;
                        // store side for start of edge
                        f_from_start_bv[f_itr] = end;
                        f_to_end_bv[f_itr] = side_is_end(t_side);
                        ++f_itr;
                    }
                }
            }
        }
    }

    // compress the forward direction side information
    util::assign(f_from_start_cbv, sd_vector<>(f_from_start_bv));
    util::assign(f_to_end_cbv, sd_vector<>(f_to_end_bv));
    
    //assert(e_iv.size() == edge_count*3);
#ifdef VERBOSE_DEBUG
    cerr << "storing reverse edges" << endl;
#endif

    size_t t_itr = 0;
    for (size_t k = 0; k < node_count; ++k) {
        //cerr << k << endl;
        int64_t t_id = i_iv[k];
        size_t t_rank = k+1;
        t_iv[t_itr] = t_rank;
        t_bv[t_itr] = 1;
        ++t_itr;
        for (auto end : { false, true }) {
            if (to_from.find(make_side(t_id, end)) != to_from.end()) {
                auto f_side_itr = to_from.find(make_side(t_id, end));
                if (f_side_itr != to_from.end()) {
                    for (auto& f_side : f_side_itr->second) {
                        size_t f_rank = id_to_rank(side_id(f_side));
                        // store link
                        t_iv[t_itr] = f_rank;
                        t_bv[t_itr] = 0;
                        // store side for end of edge
                        t_to_end_bv[t_itr] = end;
                        t_from_start_bv[t_itr] = side_is_end(f_side);
                        ++t_itr;
                    }
                }
            }
        }
    }

    // compress the reverse direction side information
    util::assign(t_to_end_cbv, sd_vector<>(t_to_end_bv));
    util::assign(t_from_start_cbv, sd_vector<>(t_from_start_bv));
    */

    // to label the paths we'll need to compress and index our vectors
    util::bit_compress(s_iv);
    util::assign(s_bv_rank, rank_support_v<1>(&s_bv));
    util::assign(s_bv_select, bit_vector::select_1_type(&s_bv));
    
    // now that we've set up our sequence indexes, we can build the locally traversable graph storage
    // calculate g_iv size
    size_t g_iv_size =
        node_count * 5 // record headers
        + edge_count * 4; // edges (stored twice) plus edge types
    util::assign(g_iv, int_vector<>(g_iv_size));
    util::assign(g_bv, bit_vector(g_iv_size));
    int64_t g = 0; // pointer into g_iv and g_bv
    // for each node
    for (int64_t i = 0; i < i_iv.size(); ++i) {
        Node n = node(i_iv[i]);
        // now build up the record
        g_bv[g] = 1; // mark record start for later query
        g_iv[g++] = n.id(); // save id
        g_iv[g++] = node_start(n.id());
        g_iv[g++] = n.sequence().size(); // sequence length
        size_t to_edge_count = 0;
        size_t from_edge_count = 0;
        size_t to_edge_count_idx = g++;
        size_t from_edge_count_idx = g++;
        // write the edges in id-based format
        // we will next convert these into relative format
        for (auto end : { false, true }) {
            auto& to_sides = to_from[make_side(n.id(), end)];
            for (auto& e : to_sides) {
                g_iv[g++] = side_id(e);
                g_iv[g++] = edge_type(side_is_end(e), end);
                ++to_edge_count;
            }
        }
        g_iv[to_edge_count_idx] = to_edge_count;
        for (auto end : { false, true }) {
            auto& from_sides = from_to[make_side(n.id(), end)];
            for (auto& e : from_sides) {
                g_iv[g++] = side_id(e);
                g_iv[g++] = edge_type(end, side_is_end(e));
                ++from_edge_count;
            }
        }
        g_iv[from_edge_count_idx] = from_edge_count;
    }

    // set up rank and select supports on g_bv so we can locate nodes in g_iv
    util::assign(g_bv_rank, rank_support_v<1>(&g_bv));
    util::assign(g_bv_select, bit_vector::select_1_type(&g_bv));

    // convert the edges in g_iv to relativistic form
    for (int64_t i = 0; i < i_iv.size(); ++i) {
        int64_t id = i_iv[i];
        // find the start of the node's record in g_iv
        int64_t g = g_bv_select(id_to_rank(id));
        // get to the edges to
        int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
        int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
        int64_t t = g + G_NODE_HEADER_LENGTH;
        int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
        for (int64_t j = t; j < f; ) {
            g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
            j += 2;
        }
        for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
            g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
            j += 2;
        }
    }

    util::bit_compress(g_iv);

#if GPBWT_MODE == MODE_SDSL
    // We have one B_s array for every side, but the first 2 numbers for sides
    // are unused. But max node rank is inclusive, so it evens out...
    bs_arrays.resize(max_node_rank() * 2);
#endif

#ifdef VERBOSE_DEBUG
    cerr << "storing paths" << endl;
#endif
    // paths
    string path_names;
    size_t path_node_count = 0; // count of node path memberships
    for (auto& pathpair : path_nodes) {
        // add path name
        const string& path_name = pathpair.first;
        //cerr << path_name << endl;
        path_names += start_marker + path_name + end_marker;
        // The path constructor helpfully counts unique path members for us
        size_t unique_member_count;
        XGPath* path = new XGPath(path_name, pathpair.second, node_count, *this, &unique_member_count);
        paths.push_back(path);
        path_node_count += unique_member_count;
    }

    // handle path names
    util::assign(pn_iv, int_vector<>(path_names.size()));
    util::assign(pn_bv, bit_vector(path_names.size()));
    // now record path name starts
    for (size_t i = 0; i < path_names.size(); ++i) {
        pn_iv[i] = path_names[i];
        if (path_names[i] == start_marker) {
            pn_bv[i] = 1; // register name start
        }
    }
    util::assign(pn_bv_rank, rank_support_v<1>(&pn_bv));
    util::assign(pn_bv_select, bit_vector::select_1_type(&pn_bv));
    
    //util::bit_compress(pn_iv);
    string path_name_file = "@pathnames.iv";
    store_to_file((const char*)path_names.c_str(), path_name_file);
    construct(pn_csa, path_name_file, 1);

    // node -> paths
    util::assign(np_iv, int_vector<>(path_node_count+node_count));
    util::assign(np_bv, bit_vector(path_node_count+node_count));
    size_t np_off = 0;
    for (size_t i = 0; i < node_count; ++i) {
        np_bv[np_off] = 1;
        np_iv[np_off] = 0; // null so we can detect entities with no path membership
        ++np_off;
        for (size_t j = 0; j < paths.size(); ++j) {
            if (paths[j]->nodes[i] == 1) {
                np_iv[np_off++] = j+1;
            }
        }
    }

    util::bit_compress(np_iv);
    //cerr << ep_off << " " << path_entities << " " << entity_count << endl;
    assert(np_off <= path_node_count+node_count);
    util::assign(np_bv_rank, rank_support_v<1>(&np_bv));
    util::assign(np_bv_select, bit_vector::select_1_type(&np_bv));

    if(store_threads) {

// Prepare empty vectors for path indexing
#ifdef VERBOSE_DEBUG
        cerr << "creating empty succinct thread store" << endl;
#endif

        util::assign(h_iv, int_vector<>(g_iv.size() * 2, 0));
        util::assign(ts_iv, int_vector<>((node_count + 1) * 2, 0));
        
#ifdef VERBOSE_DEBUG
        cerr << "storing threads" << endl;
#endif
    
        // If we're a sorted DAG we'll batch up the paths and use a batch
        // insert.
        vector<thread_t> batch;
        vector<string> batch_names;
    
        // Just store all the paths that are all perfect mappings as threads.
        // We end up converting *back* into thread_t objects.
        for (auto& pathpair : path_nodes) {
            thread_t reconstructed;
            
            // Grab the trav_ts, which are now sorted by rank
            for (auto& m : pathpair.second) {
                // Convert the mapping to a ThreadMapping
                // trav_ts are already rank sorted and deduplicated.
                ThreadMapping mapping = {trav_id(m), trav_is_rev(m)};
                reconstructed.push_back(mapping);
            }
            
#if GPBWT_MODE == MODE_SDSL
            if(is_sorted_dag) {
                // Save for a batch insert
                batch.push_back(reconstructed);
                batch_names.push_back(pathpair.first);
            }
            // TODO: else case!
#elif GPBWT_MODE == MODE_DYNAMIC
            // Insert the thread right now
            insert_thread(reconstructed, pathpair.first);
#endif
            
        }
        
#if GPBWT_MODE == MODE_SDSL
        if(is_sorted_dag) {
            // Do the batch insert
            insert_threads_into_dag(batch, batch_names);
        }
        // TODO: else case!
#endif
        util::assign(h_civ, h_iv);
        util::assign(ts_civ, ts_iv);
    }
    

#ifdef DEBUG_CONSTRUCTION
    cerr << "|g_iv| = " << size_in_mega_bytes(g_iv) << endl;
    cerr << "|g_bv| = " << size_in_mega_bytes(g_bv) << endl;
    cerr << "|s_iv| = " << size_in_mega_bytes(s_iv) << endl;
    cerr << "|f_iv| = " << size_in_mega_bytes(f_iv) << endl;
    cerr << "|t_iv| = " << size_in_mega_bytes(t_iv) << endl;

    cerr << "|f_from_start_cbv| = " << size_in_mega_bytes(f_from_start_cbv) << endl;
    cerr << "|t_to_end_cbv| = " << size_in_mega_bytes(t_to_end_cbv) << endl;

    cerr << "|f_bv| = " << size_in_mega_bytes(f_bv) << endl;
    cerr << "|t_bv| = " << size_in_mega_bytes(t_bv) << endl;

    cerr << "|i_iv| = " << size_in_mega_bytes(i_iv) << endl;
    //cerr << "|i_wt| = " << size_in_mega_bytes(i_wt) << endl;

    cerr << "|s_bv| = " << size_in_mega_bytes(s_bv) << endl;
    
    cerr << "|h_civ| = " << size_in_mega_bytes(h_civ) << endl;
    cerr << "|ts_civ| = " << size_in_mega_bytes(ts_civ) << endl;

    long double paths_mb_size = 0;
    cerr << "|pn_iv| = " << size_in_mega_bytes(pn_iv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_iv);
    cerr << "|pn_csa| = " << size_in_mega_bytes(pn_csa) << endl;
    paths_mb_size += size_in_mega_bytes(pn_csa);
    cerr << "|pn_bv| = " << size_in_mega_bytes(pn_bv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_bv);
    paths_mb_size += size_in_mega_bytes(pn_bv_rank);
    paths_mb_size += size_in_mega_bytes(pn_bv_select);
    paths_mb_size += size_in_mega_bytes(pi_iv);
    cerr << "|np_iv| = " << size_in_mega_bytes(np_iv) << endl;
    paths_mb_size += size_in_mega_bytes(np_iv);
    cerr << "|np_bv| = " << size_in_mega_bytes(np_bv) << endl;
    paths_mb_size += size_in_mega_bytes(np_bv);
    paths_mb_size += size_in_mega_bytes(np_bv_rank);
    paths_mb_size += size_in_mega_bytes(np_bv_select);
    cerr << "total paths size " << paths_mb_size << endl;
    // TODO you are missing the rest of the paths size in xg::paths
    // but this fragment should be factored into a function anyway
    
    cerr << "total size [MB] = " << (
        size_in_mega_bytes(s_iv)
        + size_in_mega_bytes(f_iv)
        + size_in_mega_bytes(t_iv)
        //+ size_in_mega_bytes(s_bv)
        + size_in_mega_bytes(f_bv)
        + size_in_mega_bytes(t_bv)
        + size_in_mega_bytes(i_iv)
        //+ size_in_mega_bytes(i_wt)
        + size_in_mega_bytes(s_bv)
        + size_in_mega_bytes(h_civ)
        + size_in_mega_bytes(ts_civ)
        // TODO: add in size of the bs_arrays in a loop
        + paths_mb_size
        ) << endl;

#endif

    if (print_graph) {
        cerr << "printing graph" << endl;
        // we have to print the relativistic graph manually because the default sdsl printer assumes unsigned integers are stored in it
        for (size_t i = 0; i < g_iv.size(); ++i) {
            cerr << (int64_t)g_iv[i] << " ";
        } cerr << endl;
        for (int64_t i = 0; i < i_iv.size(); ++i) {
            int64_t id = i_iv[i];
            // find the start of the node's record in g_iv
            int64_t g = g_bv_select(id_to_rank(id));
            // get to the edges to
            int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
            int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
            int sequence_size = g_iv[g+G_NODE_LENGTH_OFFSET];
            size_t seq_start = g_iv[g+G_NODE_SEQ_START_OFFSET];
            cerr << id << " ";
            for (int64_t j = seq_start; j < seq_start+sequence_size; ++j) {
                cerr << revdna3bit(s_iv[j]);
            } cerr << " : ";
            int64_t t = g + G_NODE_HEADER_LENGTH;
            int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
            cerr << " from ";
            for (int64_t j = t; j < f; ) {
                cerr << i_iv[g_bv_rank(g+g_iv[j])] << " ";
                j += 2;
            }
            cerr << "to ";
            for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
                cerr << i_iv[g_bv_rank(g+g_iv[j])] << " ";
                j += 2;
            }
            cerr << endl;
        }
        cerr << s_iv << endl;
        for (size_t i = 0; i < s_iv.size(); ++i) {
            cerr << revdna3bit(s_iv[i]);
        } cerr << endl;
        cerr << s_bv << endl;
        cerr << i_iv << endl;
        cerr << "paths (" << paths.size() << ")" << endl;
        for (size_t i = 0; i < paths.size(); i++) {
            // Go through paths by number, so we can determine rank
            XGPath* path = paths[i];
            
            cerr << path_name(i + 1) << endl;
            cerr << path->nodes << endl;
            // manually print IDs because simplified wavelet tree doesn't support ostream for some reason
            for (size_t j = 0; j + 1 < path->ids.size(); j++) {
                cerr << path->ids[j] << " ";
            }
            if (path->ids.size() > 0) {
                cerr << path->ids[path->ids.size() - 1];
            }
            cerr << endl;
            cerr << path->ranks << endl;
            cerr << path->directions << endl;
            cerr << path->positions << endl;
            cerr << path->offsets << endl;
        }
        cerr << np_bv << endl;
        cerr << np_iv << endl;
    }

    if (validate_graph) {
        cerr << "validating graph sequence" << endl;
        int max_id = s_bv_rank(s_bv.size());
        for (auto& p : node_label) {
            int64_t id = p.first;
            const string& l = p.second;
            //size_t rank = node_rank[id];
            size_t rank = id_to_rank(id);
            //cerr << rank << endl;
            // find the node in the array
            //cerr << "id = " << id << " rank = " << s_bv_select(rank) << endl;
            // this should be true given how we constructed things
            if (rank != s_bv_rank(s_bv_select(rank)+1)) {
                cerr << rank << " != " << s_bv_rank(s_bv_select(rank)+1) << " for node " << id << endl;
                assert(false);
            }
            // get the sequence from the s_iv
            string s = node_sequence(id);

            string ltmp, stmp;
            if (l.size() != s.size()) {
                cerr << l << " != " << endl << s << endl << " for node " << id << endl;
                assert(false);
            } else {
                int j = 0;
                for (auto c : l) {
                    if (dna3bit(c) != dna3bit(s[j++])) {
                        cerr << l << " != " << endl << s << endl << " for node " << id << endl;
                        assert(false);
                    }
                }
            }
        }
        node_label.clear();

        // -1 here seems weird
        // what?
        /*
        cerr << "validating forward edge table" << endl;
        for (size_t j = 0; j < f_iv.size()-1; ++j) {
            if (f_bv[j] == 1) continue;
            // from id == rank
            size_t fid = i_iv[f_bv_rank(j)-1];
            // to id == f_cbv[j]
            size_t tid = i_iv[f_iv[j]-1];
            bool from_start = f_from_start_bv[j];
            // get the to_end
            bool to_end = false;
            for (auto& side : from_to[make_side(fid, from_start)]) {
                if (side_id(side) == tid) {
                    to_end = side_is_end(side);
                }
            }
            bool has_edge = false;
            for (auto& side : from_to[make_side(fid, from_start)]) {
                if (side == make_side(tid, to_end)) { has_edge = true; break; }
            }
            if (!has_edge) {
                cerr << "could not find edge (f) "
                     << fid << (from_start ? "+" : "-")
                     << " -> "
                     << tid << (to_end ? "+" : "-")
                     << endl;
                assert(false);
            }
        }

        cerr << "validating reverse edge table" << endl;
        for (size_t j = 0; j < t_iv.size()-1; ++j) {
            //cerr << j << endl;
            if (t_bv[j] == 1) continue;
            // from id == rank
            size_t tid = i_iv[t_bv_rank(j)-1];
            // to id == f_cbv[j]
            size_t fid = i_iv[t_iv[j]-1];
            //cerr << tid << " " << fid << endl;

            bool to_end = t_to_end_bv[j];
            // get the to_end
            bool from_start = false;
            for (auto& side : to_from[make_side(tid, to_end)]) {
                if (side_id(side) == fid) {
                    from_start = side_is_end(side);
                }
            }
            bool has_edge = false;
            for (auto& side : to_from[make_side(tid, to_end)]) {
                if (side == make_side(fid, from_start)) { has_edge = true; break; }
            }
            if (!has_edge) {
                cerr << "could not find edge (t) "
                     << fid << (from_start ? "+" : "-")
                     << " -> "
                     << tid << (to_end ? "+" : "-")
                     << endl;
                assert(false);
            }
        }
        */

        /*
        cerr << "validating paths" << endl;
        for (auto& pathpair : path_nodes) {
            const string& name = pathpair.first;
            auto& path = pathpair.second;
            size_t prank = path_rank(name);
            //cerr << path_name(prank) << endl;
            assert(path_name(prank) == name);
            rrr_vector<>& pe_bv = paths[prank-1]->nodes;
            int_vector<>& pp_iv = paths[prank-1]->positions;
            sd_vector<>& dir_bv = paths[prank-1]->directions;
            // check each entity in the nodes is present
            // and check node reported at the positions in it
            size_t pos = 0;
            size_t in_path = 0;
            for (auto& m : path) {
                int64_t id = trav_id(m);
                bool rev = trav_is_rev(m);
                // todo rank
                assert(pe_bv[node_rank_as_entity(id)-1]);
                assert(dir_bv[in_path] == rev);
                Node n = node(id);
                //cerr << id << " in " << name << endl;
                auto p = position_in_path(id, name);
                assert(std::find(p.begin(), p.end(), pos) != p.end());
                for (size_t k = 0; k < n.sequence().size(); ++k) {
                    //cerr << "id " << id << " ==? " << node_at_path_position(name, pos+k) << endl;
                    assert(id == node_at_path_position(name, pos+k));
                    assert(id == mapping_at_path_position(name, pos+k).position().node_id());
                }
                pos += n.sequence().size();
                ++in_path;
            }
            //cerr << path_name << " rank = " << prank << endl;
            // check membership now for each entity in the path
        }
        */
        
#if GPBWT_MODE == MODE_SDSL
        if(store_threads && is_sorted_dag) {
#elif GPBWT_MODE == MODE_DYNAMIC
        if(store_threads) {
#endif
        
            cerr << "validating threads" << endl;
            
            // How many thread orientations are in the index?
            size_t threads_found = 0;
            // And how many shoukd we have inserted?
            size_t threads_expected = 0;
            list<thread_t> threads;
            for (auto& t : extract_threads(false)) {
                for (auto& k : t.second) threads.push_back(k);
            }
            for (auto& t : extract_threads(true)) {
                for (auto& k : t.second) threads.push_back(k);
            }
            for(auto thread : threads) {
#ifdef VERBOSE_DEBUG
                cerr << "Thread: ";
                for(size_t i = 0; i < thread.size(); i++) {
                    ThreadMapping mapping = thread[i];
                    cerr << mapping.node_id * 2 + mapping.is_reverse << "; ";
                }
                cerr << endl;
#endif
                // Make sure we can search all the threads we find present in the index
                assert(count_matches(thread) > 0);
                
                // Flip the thread around
                reverse(thread.begin(), thread.end());
                for(auto& mapping : thread) {
                    mapping.is_reverse = !mapping.is_reverse;
                }
                
                // We need to be able to find it backwards as well
                assert(count_matches(thread) > 0);
                
                threads_found++;
            }
            
            for (auto& pathpair : path_nodes) {
                Path reconstructed;
                
                // Grab the name
                reconstructed.set_name(pathpair.first);
                
                // This path should have been inserted. Look for it.
                assert(count_matches(reconstructed) > 0);
                
                threads_expected += 2;
                
            }
            
            // Make sure we have the right number of threads.
            assert(threads_found == threads_expected);
        }

        cerr << "graph ok" << endl;
    }
}

const uint64_t* XG::sequence_data(void) const {
    return s_iv.data();
}

const size_t XG::sequence_bit_size(void) const {
    return s_iv.bit_size();
}

bool XG::has_node(int64_t id) const {
    return id_to_rank(id) != 0;
}

Node XG::node(int64_t id) const {
    Node n;
    n.set_id(id);
    //cerr << omp_get_thread_num() << " looks for " << id << endl;
    n.set_sequence(node_sequence(id));
    return n;
}

string XG::node_sequence(int64_t id) const {
    size_t rank = id_to_rank(id);
    assert(rank != 0); // We can crash if we try to look up rank 0.
    size_t start = s_bv_select(rank);
    size_t end = rank == node_count ? s_bv.size() : s_bv_select(rank+1);
    string s; s.resize(end-start);
    for (size_t i = start; i < s_bv.size() && i < end; ++i) {
        s[i-start] = revdna3bit(s_iv[i]);
    }
    return s;
}

size_t XG::node_length(int64_t id) const {
    size_t rank = id_to_rank(id);
    size_t start = s_bv_select(rank);
    size_t end = rank == node_count ? s_bv.size() : s_bv_select(rank+1);
    return end-start;
}

char XG::pos_char(int64_t id, bool is_rev, size_t off) const {
    assert(off < node_length(id));
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t pos = s_bv_select(rank) + off;
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return c;
    } else {
        size_t rank = id_to_rank(id);
        size_t pos = s_bv_select(rank+1) - (off+1);
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return reverse_complement(c);
    }
}

string XG::pos_substr(int64_t id, bool is_rev, size_t off, size_t len) const {
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t start = s_bv_select(rank) + off;
        assert(start < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t end;
        if (!len) {
            end = s_bv_select(rank+1);
        } else {
            end = min(start + len, (size_t)s_bv_select(rank+1));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_bv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return s;
    } else {
        size_t rank = id_to_rank(id);
        size_t end = s_bv_select(rank+1) - off;
        assert(end < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t start;
        if (len > end || !len) {
            start = s_bv_select(rank);
        } else {
            start = max(end - len, (size_t)s_bv_select(rank));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_bv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return reverse_complement(s);
    }
}

size_t XG::id_to_rank(int64_t id) const {
    size_t x = id-min_id;
    if (x < 0 || x >= r_iv.size()) return 0;
    return r_iv[id-min_id];
}

int64_t XG::rank_to_id(size_t rank) const {
    if(rank == 0) {
        cerr << "[xg] error: Request for id of rank 0" << endl;
        assert(false);
    }
    if(rank > i_iv.size()) {
        cerr << "[xg] error: Request for id of rank " << rank << "/" << i_iv.size() << endl;
        assert(false);
    }
    return i_iv[rank-1];
}

int XG::edge_type(bool from_start, bool to_end) const {
    if (from_start && to_end) {
        return 4;
    } else if (from_start) {
        return 3;
    } else if (to_end) {
        return 2;
    } else {
        return 1;
    }
}

int XG::edge_type(const Edge& edge) const {
    if (edge.from_start() && edge.to_end()) {
        return 4;
    } else if (edge.from_start()) {
        return 3;
    } else if (edge.to_end()) {
        return 2;
    } else {
        return 1;
    }
}

Edge XG::edge_from_encoding(int64_t from, int64_t to, int type) const {
    Edge edge;
    edge.set_from(from);
    edge.set_to(to);
    switch (type) {
    case 1:
        break;
    case 2:
        edge.set_to_end(true);
        break;
    case 3:
        edge.set_from_start(true);
        break;
    case 4:
        edge.set_from_start(true);
        edge.set_to_end(true);
        break;
    default:
        assert(false);
        break;
    }
    return edge;
}

void XG::idify_graph(Graph& graph) const {
    // map into the id space; offsets in gciv contain the actual ids
    for (auto& node : *graph.mutable_node()) {
        node.set_id(g_iv[node.id()+G_NODE_ID_OFFSET]);
    }
    for (auto& edge : *graph.mutable_edge()) {
        edge.set_from(g_iv[edge.from()+G_NODE_ID_OFFSET]);
        edge.set_to(g_iv[edge.to()+G_NODE_ID_OFFSET]);
    }
}

Graph XG::node_subgraph_id(int64_t id) const {
    Graph graph = node_subgraph_g(g_bv_select(id_to_rank(id)));
    idify_graph(graph);
    return graph;
}

/// returns the graph with graph offsets rather than ids on edges and nodes
Graph XG::node_subgraph_g(int64_t g) const {
    Graph graph;
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    int sequence_size = g_iv[g+G_NODE_LENGTH_OFFSET];
    size_t seq_start = g_iv[g+G_NODE_SEQ_START_OFFSET];
    string sequence; sequence.resize(sequence_size);
    int i = 0;
    for (int64_t j = seq_start; j < seq_start+sequence_size; ++j, ++i) {
        sequence[i] = revdna3bit(s_iv[j]);
    }
    Node* node = graph.add_node();
    node->set_sequence(sequence);
    node->set_id(g);
    int64_t t = g + G_NODE_HEADER_LENGTH;
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    for (int64_t j = t; j < f; ) {
        int64_t from = g+g_iv[j++];
        int type = g_iv[j++];
        *graph.add_edge() = edge_from_encoding(from, g, type);
    }
    for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
        int64_t to = g+g_iv[j++];
        int type = g_iv[j++];
        *graph.add_edge() = edge_from_encoding(g, to, type);
    }
    return graph;
}

Graph XG::graph_context_id(const pos_t& pos, int64_t length) const {
    pos_t g = pos;
    vg::get_id(g) = g_bv_select(id_to_rank(id(pos)));
    Graph graph = graph_context_g(g, length);
    idify_graph(graph);
    return graph;
}

Graph XG::graph_context_g(const pos_t& g_pos, int64_t length) const {
    // helper
    auto is_inverting = [](const Edge& e) {
            return !(e.from_start() == e.to_end())
            && (e.from_start() || e.to_end());
    };
    // walk the graph from this position forward
    // adding the nodes we run into to the graph
    Graph graph;
    set<pos_t> seen;
    set<pos_t> nexts;
    nexts.insert(g_pos);
    int64_t distance = -offset(g_pos); // don't count what we won't traverse
    while (!nexts.empty()) {
        set<pos_t> todo;
        int nextd = 0;
        for (auto& next : nexts) {
            if (!seen.count(next)) {
                seen.insert(next);
                // add the node and its edges to the graph
                Graph node_graph = node_subgraph_g(id(next));
                graph.MergeFrom(node_graph);
                const Node& node = node_graph.node(0);
                nextd = nextd == 0 ? node.sequence().size() : min(nextd, (int)node.sequence().size());
                // where to next
                // look at the next positions we could reach
                pos_t pos = next;
                if (!is_rev(pos)) {
                    // we are on the forward strand, the next things from this node come off the end
                    for (auto& edge : node_graph.edge()) {
                        if((edge.to() == id(pos) && edge.to_end()) || (edge.from() == id(pos) && !edge.from_start())) {
                            id_t nid = (edge.from() == id(pos) ?
                                        edge.to()
                                        : edge.from());
                            todo.insert(make_pos_t(nid, is_inverting(edge), 0));
                        }
                    }
                } else {
                    // we are on the reverse strand, the next things from this node come off the start
                    for (auto& edge : node_graph.edge()) {
                        if((edge.to() == id(pos) && !edge.to_end()) || (edge.from() == id(pos) && edge.from_start())) {
                            id_t nid = (edge.to() == id(pos) ?
                                        edge.from()
                                        : edge.to());
                            todo.insert(make_pos_t(nid, !is_inverting(edge), 0));
                        }
                    }
                }
            }
        }
        distance += nextd;
        if (distance > length) {
            break;
        }
        nexts = todo;
    }
    return graph;
}

handle_t XG::get_handle(const id_t& node_id, bool is_reverse) const {
    // Handles will be g vector index with is_reverse in the high bit
    
    // Where in the g vector do we need to be
    size_t handle = g_bv_select(id_to_rank(node_id));
    // And set the high bit if it's reverse
    if (is_reverse) handle |= HIGH_BIT;
    return as_handle(handle);
}

id_t XG::get_id(const handle_t& handle) const {
    // Go get the g offset and then look up the noder ID
    return g_iv[(as_integer(handle) & LOW_BITS) + G_NODE_ID_OFFSET];
}

bool XG::get_is_reverse(const handle_t& handle) const {
    return as_integer(handle) & HIGH_BIT;
}

handle_t XG::flip(const handle_t& handle) const {
    return as_handle(as_integer(handle) ^ HIGH_BIT);
}

size_t XG::get_length(const handle_t& handle) const {
    return g_iv[(as_integer(handle) & LOW_BITS) + G_NODE_LENGTH_OFFSET];
}

string XG::get_sequence(const handle_t& handle) const {
    
    // Figure out how big it should be
    size_t sequence_size = get_length(handle);
    // Allocate the sequence string
    string sequence(sequence_size, '\0');
    // Extract the node record start
    size_t g = as_integer(handle) & LOW_BITS;
    // Figure out where the sequence starts
    size_t sequence_start = g_iv[g + G_NODE_SEQ_START_OFFSET];
    for (int64_t i = 0; i < sequence_size; i++) {
        // Blit the sequence out
        sequence[i] = revdna3bit(s_iv[sequence_start + i]);
    }
    
    if (as_integer(handle) & HIGH_BIT) {
        return reverse_complement(sequence);
    } else {
        return sequence;
    }
}

bool XG::edge_filter(int type, bool is_to, bool want_left, bool is_reverse) const {
    // Return true if we want an edge of the given type, where we are the from
    // or to node (according to is_to), when we are looking off the right or
    // left side of the node (according to want_left), and when the node is
    // forward or reverse (accoridng to is_reverse).
    
    // Edge type encoding:
    // 1: end to start
    // 2: end to end
    // 3: start to start
    // 4: start to end
    
    // First compute what we want looking off the right of a node in the forward direction.
    bool wanted = !is_to && (type == 1 || type == 2) || is_to && (type == 2 || type == 4);
    
    // We computed whether we wanted it assuming we were looking off the right. The complement is what we want looking off the left.
    wanted = wanted != want_left;
    
    // We computed whether we wanted ot assuming we were in the forward orientation. The complement is what we want in the reverse orientation.
    wanted = wanted != is_reverse;
    
    return wanted;
}

bool XG::do_edges(const size_t& g, const size_t& start, const size_t& count, bool is_to,
    bool want_left, bool is_reverse, const function<bool(const handle_t&)>& iteratee) const {
    
    // OK go over all those edges
    for (size_t i = 0; i < count; i++) {
        // What edge type is the edge?
        int type = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_TYPE_OFFSET];
        
        // Make sure we got a valid edge type and we haven't wandered off into non-edge data.
        assert(type >= 0);
        assert(type <= 3);
        
        if (edge_filter(type, is_to, want_left, is_reverse)) {
            
            // What's the offset to the other node?
            int64_t offset = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_OFFSET_OFFSET];
            
            // Make sure we haven't gone off the rails into non-edge data.
            assert((int64_t) g + offset >= 0);
            assert(g + offset < g_iv.size());
            
            // Should we invert?
            // We only invert if we cross an end to end edge. Or a start to start edge
            bool new_reverse = is_reverse != (type == 2 || type == 3);
            
            // Compose the handle for where we are going
            handle_t next_handle = as_handle((g + offset) | (new_reverse ? HIGH_BIT : 0));
            
            // We want this edge
            
            if (!iteratee(next_handle)) {
                // Stop iterating
                return false;
            }
        }
        else {
            // TODO: delete this after using it to debug
            int64_t offset = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_OFFSET_OFFSET];
            bool new_reverse = is_reverse != (type == 2 || type == 3);
            handle_t next_handle = as_handle((g + offset) | (new_reverse ? HIGH_BIT : 0));
        }
    }
    // Iteratee didn't stop us.
    return true;
}

void XG::follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {

    // Unpack the handle
    size_t g = as_integer(handle) & LOW_BITS;
    bool is_reverse = get_is_reverse(handle);

    // How many edges are there of each type?
    size_t edges_to_count = g_iv[g + G_NODE_TO_COUNT_OFFSET];
    size_t edges_from_count = g_iv[g + G_NODE_FROM_COUNT_OFFSET];
    
    // Where does each edge run start?
    size_t to_start = g + G_NODE_HEADER_LENGTH;
    size_t from_start = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    
    // We will look for all the edges on the appropriate side, which means we have to check the from and to edges
    if (do_edges(g, to_start, edges_to_count, true, go_left, is_reverse, iteratee)) {
        // All the edges where we're to were accepted, so do the edges where we're from
        do_edges(g, from_start, edges_from_count, false, go_left, is_reverse, iteratee);
    }
}

void XG::for_each_handle(const function<bool(const handle_t&)>& iteratee) const {
    // How big is the g vector entry size we are on?
    size_t entry_size = 0;
    for (size_t g = 0; g < g_iv.size(); g += entry_size) {
        // Make the handle
        // We need to make sure our index won't set the orientation bit.
        assert((g & (~LOW_BITS)) == 0);
        
        // Just make it into a handle; we're always forward.
        handle_t handle = as_handle(g);
        
        // Run the iteratee
        if (!iteratee(handle)) {
            // The iteratee is bored and wants to stop.
            return;
        }
    
        // How many edges are there of each type on this record?
        size_t edges_to_count = g_iv[g + G_NODE_TO_COUNT_OFFSET];
        size_t edges_from_count = g_iv[g + G_NODE_FROM_COUNT_OFFSET];
        
        // This record is the header plus all the edge records it contains
        entry_size = G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * (edges_to_count + edges_from_count);
    }
}

size_t XG::node_size() const {
    return this->node_count;
}

vector<Edge> XG::edges_of(int64_t id) const {
    size_t g = g_bv_select(id_to_rank(id));
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    int64_t t = g + G_NODE_HEADER_LENGTH;
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    vector<Edge> edges;
    for (int64_t j = t; j < f; ) {
        int64_t from = g+g_iv[j++];
        int type = g_iv[j++];
        edges.push_back(edge_from_encoding(from, g, type));
    }
    for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
        int64_t to = g+g_iv[j++];
        int type = g_iv[j++];
        edges.push_back(edge_from_encoding(g, to, type));
    }
    for (auto& edge : edges) { 
        edge.set_from(g_iv[edge.from()+G_NODE_ID_OFFSET]);
        edge.set_to(g_iv[edge.to()+G_NODE_ID_OFFSET]);
    }
    return edges;
}

vector<Edge> XG::edges_to(int64_t id) const {
    size_t g = g_bv_select(id_to_rank(id));
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    int64_t t = g + G_NODE_HEADER_LENGTH;
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    vector<Edge> edges;
    for (int64_t j = t; j < f; ) {
        int64_t from = g+g_iv[j++];
        int type = g_iv[j++];
        edges.push_back(edge_from_encoding(from, g, type));
    }
    for (auto& edge : edges) { 
        edge.set_from(g_iv[edge.from()+G_NODE_ID_OFFSET]);
        edge.set_to(g_iv[edge.to()+G_NODE_ID_OFFSET]);
    }
    return edges;
}

vector<Edge> XG::edges_from(int64_t id) const {
    size_t g = g_bv_select(id_to_rank(id));
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    //int64_t t = g + G_NODE_HEADER_LENGTH;
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    int64_t e = f + G_EDGE_LENGTH * edges_from_count;
    vector<Edge> edges;
    for (int64_t j = f; j < e; ) {
        int64_t to = g+g_iv[j++];
        int type = g_iv[j++];
        edges.push_back(edge_from_encoding(g, to, type));
    }
    for (auto& edge : edges) { 
        edge.set_from(g_iv[edge.from()+G_NODE_ID_OFFSET]);
        edge.set_to(g_iv[edge.to()+G_NODE_ID_OFFSET]);
    }
    return edges;
}

vector<Edge> XG::edges_on_start(int64_t id) const {
    vector<Edge> edges;
    for (auto& edge : edges_of(id)) {
        if((edge.to() == id && !edge.to_end()) || (edge.from() == id && edge.from_start())) {
            edges.push_back(edge);
        }
    }
    return edges;
}

vector<Edge> XG::edges_on_end(int64_t id) const {
    vector<Edge> edges;
    for (auto& edge : edges_of(id)) {
        if((edge.to() == id && edge.to_end()) || (edge.from() == id && !edge.from_start())) {
            edges.push_back(edge);
        }
    }
    return edges;
}

size_t XG::max_node_rank(void) const {
    return s_bv_rank(s_bv.size());
}

int64_t XG::node_at_seq_pos(size_t pos) const {
    return rank_to_id(s_bv_rank(pos));
}

size_t XG::node_start(int64_t id) const {
    return s_bv_select(id_to_rank(id));
}

size_t XG::max_path_rank(void) const {
    //cerr << pn_bv << endl;
    //cerr << "..." << pn_bv_rank(pn_bv.size()) << endl;
    return pn_bv_rank(pn_bv.size());
}

// snoop through the forward table to check if the edge exists
bool XG::has_edge(int64_t id1, bool from_start, int64_t id2, bool to_end) const {
    Edge query;
    query.set_from(id1);
    query.set_from_start(from_start);
    query.set_to(id2);
    query.set_to_end(to_end);
    query = canonicalize(query);
    Graph node_graph = node_subgraph_id(id1);
    for (auto& edge : node_graph.edge()) {
        if (edge.from() == query.from()
            && edge.from_start() == query.from_start()
            && edge.to() == query.to()
            && edge.to_end() == query.to_end()) {
            return true;
        }
    }
    return false;
}

bool XG::has_edge(const Edge& edge) const {
    auto fixed = canonicalize(edge);
    return has_edge(fixed.from(), fixed.from_start(), fixed.to(), fixed.to_end());
}

size_t XG::node_graph_idx(int64_t id) const {
    return g_bv_select(id_to_rank(id));
}

size_t XG::edge_graph_idx(const Edge& edge_in) const {
    auto edge = canonicalize(edge_in);
    int64_t id = edge.from();
    size_t g = g_bv_select(id_to_rank(id));
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    //int64_t t = g + G_NODE_HEADER_LENGTH;
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    int64_t e = f + G_EDGE_LENGTH * edges_from_count;
    vector<Edge> edges;
    int i = 1;
    for (int64_t j = f; j < e; ++i) {
        int64_t to = g+g_iv[j++];
        int type = g_iv[j++];
        Edge curr = edge_from_encoding(g, to, type);
        curr.set_from(g_iv[curr.from()+G_NODE_ID_OFFSET]);
        curr.set_to(g_iv[curr.to()+G_NODE_ID_OFFSET]);
        if (curr.to() == edge.to()
            && curr.from_start() == edge.from_start()
            && curr.to_end() == edge.to_end()) {
            return g + i;
        }
    }
    assert(false);
    return 0;
}

Edge XG::canonicalize(const Edge& edge) const {
    // An edge is canonical if it is not doubly reversing and, if it is singly
    // reversing, the lower side comes first.
    // XG can't actually handle doubly reversing edges, I think.

    if ((edge.from_start() && edge.to_end()) || ((edge.from_start() || edge.to_end()) && edge.from() > edge.to())) {
        // Doubly reversing or (singly reversing and in the wrong direction)
        return make_edge(edge.to(), !edge.to_end(), edge.from(), !edge.from_start());
    } else {
        return edge;
    }
}

Path XG::path(const string& name) const {
    // Extract a whole path by name
    
    // First find the XGPath we're using to store it.
    const XGPath& xgpath = *(paths[path_rank(name)-1]);
    
    // Make a new path to fill in
    Path to_return;
    // Fill in the name
    to_return.set_name(name);
    
    // There's one ID entry per node visit    
    size_t total_nodes = xgpath.ids.size();
    
    for(size_t i = 0; i < total_nodes; i++) {
        // For everything on the XGPath, put a Mapping on the real path.
        Mapping* m = to_return.add_mapping();
        *m = xgpath.mapping(i);
        // Add one full length match edit, because the XGPath doesn't know how
        // to make it.
        Edit* e = m->add_edit();
        e->set_from_length(node_length(m->position().node_id()));
        e->set_to_length(e->from_length());
    }
    
    return to_return;
    
}

size_t XG::path_rank(const string& name) const {
    // find the name in the csa
    string query = start_marker + name + end_marker;
    auto occs = locate(pn_csa, query);
    if (occs.size() > 1) {
        cerr << "multiple hits for " << query << endl;
        assert(false);
    }
    if(occs.size() == 0) {
        // This path does not exist. Give back 0, which can never be a real path
        // rank.
        return 0;
    }
    //cerr << "path named " << name << " is at " << occs[0] << endl;
    return pn_bv_rank(occs[0])+1; // step past '#'
}

string XG::path_name(size_t rank) const {
    size_t start = pn_bv_select(rank)+1; // step past '#'
    size_t end = rank == path_count ? pn_iv.size() : pn_bv_select(rank+1);
    end -= 1;  // step before '$'
    string name; name.resize(end-start);
    for (size_t i = start; i < end; ++i) {
        name[i-start] = pn_iv[i];
    }
    return name;
}

bool XG::path_contains_node(const string& name, int64_t id) const {
    return 1 == paths[path_rank(name)-1]->nodes[id_to_rank(id)-1];
}

vector<size_t> XG::paths_of_node(int64_t id) const {
    size_t off = np_bv_select(id_to_rank(id));
    assert(np_bv[off++]);
    vector<size_t> path_ranks;
    while (off < np_bv.size() ? np_bv[off] == 0 : false) {
        path_ranks.push_back(np_iv[off++]);
    }
    return path_ranks;
}

vector<pair<size_t, vector<pair<size_t, bool>>>> XG::oriented_paths_of_node(int64_t id) const {
    vector<pair<size_t, vector<pair<size_t, bool>>>> path_occurrences;
    for (size_t path_rank : paths_of_node(id)) {
        path_occurrences.emplace_back(path_rank,  vector<pair<size_t, bool>>());
        auto& occurrences = path_occurrences.back().second;
        for (size_t i : node_ranks_in_path(id, path_rank)) {
            occurrences.emplace_back(i, paths[path_rank-1]->directions[i]);
        }
    }
    return path_occurrences;
}
    
map<string, vector<Mapping>> XG::node_mappings(int64_t id) const {
    map<string, vector<Mapping>> mappings;
    // for each time the node crosses the path
    for (auto i : paths_of_node(id)) {
        // get the path name
        string name = path_name(i);
        // get reference to the offset of the mapping in the path
        // to get the direction and (stored) rank
        for (auto j : node_ranks_in_path(id, name)) {
            // nb: path rank is 1-based, path index is 0-based
            mappings[name].push_back(paths[i-1]->mapping(j));
        }
    }
    return mappings;
}

void XG::neighborhood(int64_t id, size_t dist, Graph& g, bool use_steps) const {
    *g.add_node() = node(id);
    expand_context(g, dist, true, use_steps);
}

void XG::expand_context(Graph& g, size_t dist, bool add_paths, bool use_steps,
                        bool expand_forward, bool expand_backward,
                        int64_t until_node) const {
    if (use_steps) {
        expand_context_by_steps(g, dist, add_paths, expand_forward, expand_backward, until_node);
    } else {
        expand_context_by_length(g, dist, add_paths, expand_forward, expand_backward, until_node);
    }
}

void XG::expand_context_by_steps(Graph& g, size_t steps, bool add_paths,
                                 bool expand_forward, bool expand_backward,
                                 int64_t until_node) const {
    map<int64_t, Node*> nodes;
    map<pair<side_t, side_t>, Edge*> edges;
    set<int64_t> to_visit;
    // start with the nodes in the graph
    for (size_t i = 0; i < g.node_size(); ++i) {
        to_visit.insert(g.node(i).id());
        // handles the single-node case: we should still get the paths
        Node* np = g.mutable_node(i);
        nodes[np->id()] = np;
    }
    for (size_t i = 0; i < g.edge_size(); ++i) {
        auto& edge = g.edge(i);
        to_visit.insert(edge.from());
        to_visit.insert(edge.to());
        edges[make_pair(make_side(edge.from(), edge.from_start()),
                        make_side(edge.to(), edge.to_end()))] = g.mutable_edge(i);
    }
    // and expand
    for (size_t i = 0; i < steps; ++i) {
        set<int64_t> to_visit_next;
        for (auto id : to_visit) {
            // build out the graph
            // if we have nodes we haven't seeen
            if (nodes.find(id) == nodes.end()) {
                Node* np = g.add_node();
                nodes[id] = np;
                *np = node(id);
            }
            vector<Edge> edges_todo;
            if (expand_forward && expand_backward) {
                edges_todo = edges_of(id);
            } else if (expand_forward) {
                edges_todo = edges_from(id);
            } else if (expand_backward) {
                edges_todo = edges_to(id);
            } else {
                cerr << "[xg] error: Requested neither forward no backward context expansion" << endl;
                exit(1);
            }
            for (auto& edge : edges_todo) {
                auto sides = make_pair(make_side(edge.from(), edge.from_start()),
                                       make_side(edge.to(), edge.to_end()));
                if (edges.find(sides) == edges.end()) {
                    Edge* ep = g.add_edge(); *ep = edge;
                    edges[sides] = ep;
                }
                if (edge.from() == id) {
                    to_visit_next.insert(edge.to());
                } else {
                    to_visit_next.insert(edge.from());
                }
            }
            if (until_node != 0 && nodes.find(until_node) != nodes.end()) {
                break;
            }
        }
        to_visit = to_visit_next;
    }
    // then add connected nodes that we have edges to but didn't pull in yet.
    // These are the nodes reached on the last step; we won't follow their edges
    // to new noded.
    set<int64_t> last_step_nodes;
    for (auto& e : edges) {
        auto& edge = e.second;
        // get missing nodes
        int64_t f = edge->from();
        if (nodes.find(f) == nodes.end()) {
            Node* np = g.add_node();
            nodes[f] = np;
            *np = node(f);
            last_step_nodes.insert(f);
        }
        int64_t t = edge->to();
        if (nodes.find(t) == nodes.end()) {
            Node* np = g.add_node();
            nodes[t] = np;
            *np = node(t);
            last_step_nodes.insert(t);
        }
    }
    // We do need to find edges that connect the nodes we just grabbed on the
    // last step. Otherwise we'll produce something that isn't really a useful
    // subgraph, because there might be edges connecting the nodes you have that
    // you don't see.
    for(auto& n : last_step_nodes) {
        for (auto& edge : edges_from(n)) {
            if(last_step_nodes.count(edge.to())) {
                // This edge connects two nodes that were added on the last
                // step, and so wouldn't have been found by the main loop.
                
                // Determine if it's been seen already (somehow).
                // TODO: it probably shouldn't have been, unless it's a self loop or something.
                auto sides = make_pair(make_side(edge.from(), edge.from_start()),
                                       make_side(edge.to(), edge.to_end()));
                if (edges.find(sides) == edges.end()) {
                    // If it isn't there already, add it to the graph
                    Edge* ep = g.add_edge(); *ep = edge;
                    edges[sides] = ep;
                }
            }
        }       
    }
    // Edges between the last step nodes and other nodes will have already been
    // pulled in, on the step when those other nodes were processed by the main
    // loop.
    if (add_paths) {
        add_paths_to_graph(nodes, g);
    }
}

void XG::expand_context_by_length(Graph& g, size_t length, bool add_paths,
                                  bool expand_forward, bool expand_backward,
                                  int64_t until_node) const {

    // map node_id --> min-distance-to-left-side, min-distance-to-right-side
    // these distances include the length of the node in the table. 
    map<int64_t, pair<int64_t, int64_t> > node_table;
    // nodes and edges in graph, so we don't duplicate when we add to protobuf
    map<int64_t, Node*> nodes;
    map<pair<side_t, side_t>, Edge*> edges;
    // bfs queue (id, enter-on-left-size)
    queue<int64_t> to_visit;

    // add starting graph with distance 0
    for (size_t i = 0; i < g.node_size(); ++i) {
        Node* np = g.mutable_node(i);
        node_table[np->id()] = pair<int64_t, int64_t>(0, 0);
        nodes[np->id()] = np;
        to_visit.push(np->id());
    }

    // add starting edges
    for (size_t i = 0; i < g.edge_size(); ++i) {
        auto& edge = g.edge(i);
        edges[make_pair(make_side(edge.from(), edge.from_start()),
                        make_side(edge.to(), edge.to_end()))] = g.mutable_edge(i);
    }

    // expand outward breadth-first
    while (!to_visit.empty() && (until_node == 0 || nodes.find(until_node) == nodes.end())) {
        int64_t id = to_visit.front();
        to_visit.pop();
        pair<int64_t, int64_t> dists = node_table[id];
        if (dists.first < length || dists.second < length) {
            vector<Edge> edges_todo;
            if (expand_forward && expand_backward) {
                edges_todo = edges_of(id);
            } else if (expand_forward) {
                edges_todo = edges_from(id);
            } else if (expand_backward) {
                edges_todo = edges_to(id);
            } else {
                cerr << "[xg] error: Requested neither forward no backward context expansion" << endl;
                exit(1);
            }
            for (auto& edge : edges_todo) {
                // update distance table with other end of edge
                function<void(int64_t, bool, bool)> lambda = [&](
                    int64_t other, bool from_start, bool to_end) {

                    int64_t dist = !from_start ? dists.first : dists.second;
                    Node other_node = node(other);
                    int64_t other_dist = dist + other_node.sequence().size();
                    if (dist < length) {
                        auto it = node_table.find(other);
                        bool updated = false;
                        if (it == node_table.end()) {
                            auto entry = make_pair(numeric_limits<int64_t>::max(),
                                                   numeric_limits<int64_t>::max());
                            it = node_table.insert(make_pair(other, entry)).first;
                            updated = true;
                        }
                        if (!to_end && other_dist < it->second.first) {
                            updated = true;
                            node_table[other].first = other_dist;
                        } else if (to_end && other_dist < it->second.second) {
                            updated = true;
                            node_table[other].second = other_dist;
                        }
                        // create the other node
                        if (nodes.find(other) == nodes.end()) {
                            Node* np = g.add_node();
                            nodes[other] = np;
                            *np = other_node;
                        }
                        // create all links back to graph, so as not to break paths
                        for (auto& other_edge : edges_of(other)) {
                            auto sides = make_pair(make_side(other_edge.from(),
                                                             other_edge.from_start()),
                                                   make_side(other_edge.to(),
                                                             other_edge.to_end()));
                            int64_t other_from = other_edge.from() == other ? other_edge.to() :
                                other_edge.from();
                            if (nodes.find(other_from) != nodes.end() &&
                                edges.find(sides) == edges.end()) {
                                Edge* ep = g.add_edge(); *ep = other_edge;
                                edges[sides] = ep;
                            }
                        }
                        // revisit the other node
                        if (updated) {
                            // this may be overly conservative (bumping any updated node)
                            to_visit.push(other);
                        }
                    }
                };
                // we can actually do two updates if we have a self loop, hence no else below
                if (edge.from() == id) {
                    lambda(edge.to(), edge.from_start(), edge.to_end());
                }
                if (edge.to() == id) {
                    lambda(edge.from(), !edge.to_end(), !edge.from_start());
                }
            }
        }
    }

    if (add_paths) {
        add_paths_to_graph(nodes, g);
    }
}
    
// if the graph ids partially ordered, this works no prob
// otherwise... owch
// the paths become disordered due to traversal of the node ids in order
void XG::add_paths_to_graph(map<int64_t, Node*>& nodes, Graph& g) const {
    // map from path name to (map from mapping rank to mapping)
    map<string, map<size_t, Mapping>> paths;
    // mappings without 
    map<string, vector<Mapping>> unplaced;
    // use:
    //size_t node_position_in_path(int64_t id, const string& name) const;

    // pick up current paths in the graph
    for (size_t i = 0; i < g.path_size(); ++i) {
        auto& path = g.path(i);
        for (size_t j = 0; j < path.mapping_size(); ++j) {
            auto& m = path.mapping(j);
            if (m.rank()) {
                paths[path.name()][m.rank()] = m;
            } else {
                unplaced[path.name()].push_back(m);
            }                
        }
    }
    // do the same for the mappings in the list of nodes
    for (auto& n : nodes) {
        auto& id = n.first;
        auto& node = n.second;
        for (auto& n : node_mappings(id)) {
            auto& name = n.first;
            for (auto& m : n.second) {
                if (m.rank()) {
                    paths[name][m.rank()] = m;
                } else {
                    unplaced[name].push_back(m);
                }
            }
        }
    }
    // rebuild graph's paths
    // NB: mapping ranks allow us to remove this bit
    // only adding what we haven't seen before
    g.clear_path();
    for (auto& p : paths) {
        auto& name = p.first;
        auto& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        for (auto& n : mappings) {
            *path->add_mapping() = n.second;
        }
        if (unplaced.find(name) != unplaced.end()) {
            auto& unp = unplaced[name];
            for (auto& m : unp) {
                *path->add_mapping() = m;
            }
        }
    }
}

void XG::get_id_range(int64_t id1, int64_t id2, Graph& g) const {
    id1 = max(min_id, id1);
    id2 = min(max_id, id2);
    for (auto i = id1; i <= id2; ++i) {
        if(id_to_rank(i) != 0) { 
            // We actually have a node with that ID.
            *g.add_node() = node(i);
        }
    }
}

// walk forward in id space, collecting nodes, until at least length bases covered
// (or end of graph reached).  if forward is false, do go backward
void XG::get_id_range_by_length(int64_t id, int64_t length, Graph& g, bool forward) const {
    // find out first base of node's position in the sequence vector
    size_t rank = id_to_rank(id);
    size_t start = s_bv_select(rank);
    size_t end;
    // jump by length, checking to make sure we stay in bounds
    if (forward) {
        end = s_bv_rank(min(s_bv.size() - 1, start + node_length(id) + length));
    } else {
        end = s_bv_rank(1 + max((int64_t)0, (int64_t)(start  - length)));
    }
    // convert back to id
    int64_t id2 = rank_to_id(end);

    // get the id range
    if (!forward) {
        swap(id, id2);
    }
    get_id_range(id, id2, g);
}


/*
void XG::get_connected_nodes(Graph& g) {
}
*/

size_t XG::path_length(const string& name) const {
    return paths[path_rank(name)-1]->offsets.size();
}

size_t XG::path_length(size_t rank) const {
    return paths[rank-1]->offsets.size();
}

pair<int64_t, vector<size_t> > XG::nearest_path_node(int64_t id, int max_steps) const {
    set<int64_t> todo;
    set<int64_t> seen;
    todo.insert(id);
    int i = 0;
    while (!todo.empty() && i++ < max_steps) {
        set<int64_t> next;
        for (auto& id : todo) {
            if (seen.count(id)) continue;
            seen.insert(id);
            vector<size_t> path_ids = paths_of_node(id);
            if (!path_ids.empty()) {
                // if we found a node on a path, return
                return make_pair(id, path_ids);
            } else {
                for (auto& edge : edges_of(id)) {
                    next.insert(edge.from());
                    next.insert(edge.to());
                }
            }
        }
        todo = next;
    }
    return make_pair(id, vector<size_t>());
}

// like above, but find minumum over list of paths.  if names is empty, do all paths
// don't actually take strict minumum over all paths.  rather, prefer paths that
// contain the nodes when possible. 
int64_t XG::min_approx_path_distance(int64_t id1, int64_t id2) const {

    int64_t min_distance = numeric_limits<int64_t>::max();
    pair<int64_t, vector<size_t> > near1 = nearest_path_node(id1);
    pair<int64_t, vector<size_t> > near2 = nearest_path_node(id2);
    if (near1.second.size() || near2.second.size()) {
        map<int, size_t> paths;
        for (auto& i : near1.second) paths[i]++;
        for (auto& i : near2.second) paths[i]++;
        for (auto& i : paths) {
            if (i.second < 2) continue;
            auto name = path_name(i.first);
            for (auto& p1 : position_in_path(near1.first, name)) {
                for (auto& p2 : position_in_path(near2.first, name)) {
                    int64_t distance = abs((int64_t)p1 - (int64_t)p2);
                    min_distance = min(distance, min_distance);
                }
            }
        }
    }
    return min_distance;
}
    
void XG::memoized_oriented_paths_of_node(int64_t id,
                                         vector<pair<size_t, vector<pair<size_t, bool>>>>& local_paths_var,
                                         vector<pair<size_t, vector<pair<size_t, bool>>>>*& paths_of_node_ptr_out,
                                         unordered_map<int64_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo) const {
    if (paths_of_node_memo) {
        auto iter = paths_of_node_memo->find(id);
        if (iter != paths_of_node_memo->end()) {
            paths_of_node_ptr_out = &(iter->second);
        }
        else {
            (*paths_of_node_memo)[id] = oriented_paths_of_node(id);
            paths_of_node_ptr_out = &((*paths_of_node_memo)[id]);
        }
    }
    else {
        local_paths_var = oriented_paths_of_node(id);
        paths_of_node_ptr_out = &local_paths_var;
    }
}
    
handle_t XG::memoized_get_handle(int64_t id, bool rev, unordered_map<pair<int64_t, bool>, handle_t>* handle_memo) const {
    if (handle_memo) {
        auto iter = handle_memo->find(make_pair(id, rev));
        if (iter != handle_memo->end()) {
            return iter->second;
        }
        else {
            handle_t handle = get_handle(id, rev);
            (*handle_memo)[make_pair(id, rev)] = handle;
            return handle;
        }
    }
    else {
        return get_handle(id, rev);
    }
}
    
int64_t XG::closest_shared_path_oriented_distance(int64_t id1, size_t offset1, bool rev1,
                                                  int64_t id2, size_t offset2, bool rev2,
                                                  size_t max_search_dist,
                                                  unordered_map<id_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo,
                                                  unordered_map<pair<int64_t, bool>, handle_t>* handle_memo) const {
    
#ifdef debug_algorithms
    cerr << "[XG] estimating oriented distance between " << id1 << "[" << offset1 << "]" << (rev1 ? "-" : "+") << " and " << id2 << "[" << offset2 << "]" << (rev2 ? "-" : "+") << " with max search distance of " << max_search_dist << endl;
#endif
    // maps of oriented paths to (node id, strand, oriented distance) tuples
    unordered_map<pair<size_t, bool>, tuple<int64_t, bool, int64_t>> path_dists_1;
    unordered_map<pair<size_t, bool>, tuple<int64_t, bool, int64_t>> path_dists_2;
    
    unordered_set<pair<size_t, bool>> shared_paths;
    
    // a local variable to store local results if no memo is provided
    vector<pair<size_t, vector<pair<size_t, bool>>>> local_paths_var1, local_paths_var2;
    // a pointer that will be set to point to valid results of oriented_paths_of_node
    vector<pair<size_t, vector<pair<size_t, bool>>>>* curr_node_path_occurrences;
    
    // ensure that the paths of the start nodes are added, even if their ends are too far away
    // from the positions for the search to explore
    // TODO: this leaves the ambiguity that a node might occur multiple times on the same path, in which case
    // the tie for closest traversal to the path is broken arbitrarily
    memoized_oriented_paths_of_node(id1, local_paths_var1, curr_node_path_occurrences, paths_of_node_memo);
    for (pair<size_t, vector<pair<size_t, bool>>>& oriented_occurrences : *curr_node_path_occurrences) {
        for (const pair<size_t, bool>& occurrence : oriented_occurrences.second) {
            pair<size_t, bool> path_occurrence(oriented_occurrences.first, occurrence.second != rev1);
            path_dists_1[path_occurrence] = make_tuple(id1, rev1, -((int64_t) offset1));
            
#ifdef debug_algorithms
            cerr << "[XG] first position " << id1 << "[" << offset1 << "]" << (rev1 ? "-" : "+") << " has an initial path occurrence on " << path_occurrence.first << (path_occurrence.second ? "-" : "+") << endl;
#endif
        }
    }
    memoized_oriented_paths_of_node(id2, local_paths_var1, curr_node_path_occurrences, paths_of_node_memo);
    for (pair<size_t, vector<pair<size_t, bool>>>& oriented_occurrences : *curr_node_path_occurrences) {
        for (const pair<size_t, bool>& occurrence : oriented_occurrences.second) {
            pair<size_t, bool> path_occurrence(oriented_occurrences.first, occurrence.second != rev2);
            path_dists_2[path_occurrence] = make_tuple(id2, rev2, -((int64_t) offset2));
            
#ifdef debug_algorithms
            cerr << "[XG] second position " << id2 << "[" << offset2 << "]" << (rev2 ? "-" : "+") << " has an initial path occurrence on " << path_occurrence.first << (path_occurrence.second ? "-" : "+") << endl;
#endif
            
            if (path_dists_1.count(path_occurrence)) {
#ifdef debug_algorithms
                cerr << "[XG] this occurrence is on a shared path" << endl;
#endif
                shared_paths.insert(path_occurrence);
            }
        }
    }
    
    // a local struct for traversals that are ordered by search distance and keep
    // track of whether we're searching to the left or right
    struct Traversal {
        Traversal() {}
        Traversal(int64_t dist, handle_t handle, bool search_left) :
                 dist(dist), handle(handle), search_left(search_left) {}
        handle_t handle;
        int64_t dist;
        bool search_left;
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    // if we already found shared paths on the start nodes, don't search anymore
    if (shared_paths.empty() && max_search_dist > 0) {
#ifdef debug_algorithms
        cerr << "[XG] no shared paths detected, beginning traversals" << endl;
#endif
        // priority queues over traversals
        // distance is measure at the end of the node, so it's actually the distance
        // to the next nodes we will traverse to
        // there is a separate queue for each of the positions
        priority_queue<Traversal> queue_1, queue_2;
        unordered_set<handle_t> queued_1, queued_2;
        
#ifdef debug_algorithms
        cerr << "[XG] getting handles for search starts" << endl;
#endif
        // get handles to the starting positions
        handle_t handle1 = memoized_get_handle(id1, rev1, handle_memo);
        handle_t handle2 = memoized_get_handle(id2, rev2, handle_memo);
        
#ifdef debug_algorithms
        cerr << "[XG] initializing queues" << endl;
#endif

        // add in the initial traversals in both directions from both start nodes
        queue_1.emplace((int64_t) offset1 - (int64_t) get_length(handle1), handle1, true);
        queue_1.emplace(-((int64_t) offset1), handle1, false);
        queue_2.emplace((int64_t) offset2 - (int64_t) get_length(handle2), handle2, true);
        queue_2.emplace(-((int64_t) offset2), handle2, false);
        
        // mark the initial traversals as queued
        queued_1.insert(handle1);
        queued_2.insert(handle2);
        
        while (!(queue_1.empty() && queue_2.empty()) && shared_paths.empty()) {
            // get the queue that has the next shortest path
#ifdef debug_algorithms
            cerr << "[XG] choosing queue for next traversal" << endl;
#endif
            auto curr_queue = &queue_1;
            auto curr_queued = &queued_1;
            auto curr_path_dists = &path_dists_1;
            auto other_path_dists = &path_dists_2;
            if (queue_1.empty() ? true : (queue_2.empty() ? false : queue_1.top().dist > queue_2.top().dist)) {
                curr_queue = &queue_2;
                curr_queued = &queued_2;
                std::swap(curr_path_dists, other_path_dists);
            }
            
            Traversal trav = curr_queue->top();
            curr_queue->pop();
            
#ifdef debug_algorithms
            cerr << "[XG] traversing " << get_id(trav.handle) << (get_is_reverse(trav.handle) ? "-" : "+") << " in " << (trav.search_left ? "leftward" : "rightward") << " direction at distance " << trav.dist << endl;
#endif
            
            // don't look any further if the next closest traversal is beyond the maximum distance
            if ( trav.dist > (int64_t) max_search_dist) {
                break;
            }
            
            int64_t dist = trav.dist + get_length(trav.handle);
            
            int64_t trav_id = get_id(trav.handle);
            bool trav_is_rev = get_is_reverse(trav.handle);
            
            if ((trav_id != id1 || trav_is_rev != rev1) && (trav_id != id2 || trav_is_rev != rev2)) {
                // this is not one of the start positions, so it might have new paths on it
                memoized_oriented_paths_of_node(trav_id, local_paths_var1, curr_node_path_occurrences, paths_of_node_memo);
                for (pair<size_t, vector<pair<size_t, bool>>>& oriented_occurrences : *curr_node_path_occurrences) {
                    for (const pair<size_t, bool>& occurrence : oriented_occurrences.second) {
#ifdef debug_algorithms
                        cerr << "\tnode is on path " << oriented_occurrences.first << " in " << (occurrence.second ? "reverse" : "forward") << " orientation" << endl;
#endif
                        pair<size_t, bool> path_orientation(oriented_occurrences.first, occurrence.second != trav_is_rev);
                        if (!curr_path_dists->count(path_orientation)) {
                            // record the oriented distance to the forward beginning of the node, relative to the start traversal
                            (*curr_path_dists)[path_orientation] = make_tuple(trav_id, trav_is_rev,
                                                                              trav.search_left ? -dist : trav.dist);
                            // have we found nodes that share a path yet?
                            if (other_path_dists->count(path_orientation)) {
                                shared_paths.insert(path_orientation);
                            }
                        }
                    }
                }
            }
            
            function<bool(const handle_t& next)> enqueue_next = [&](const handle_t& next) {
                if (!curr_queued->count(next)) {
#ifdef debug_algorithms
                    cerr << "\tfollowing edge to " << get_id(next) << (get_is_reverse(trav.handle) ? "-" : "+") << " at dist " << dist << endl;
#endif
                    curr_queued->emplace(next);
                    curr_queue->emplace(dist, next, trav.search_left);
                }
                return true;
            };
            
            follow_edges(trav.handle, trav.search_left, enqueue_next);
        }
    }
    
#ifdef debug_algorithms
    cerr << "[XG] found a shared path or exhausted search distance" << endl;
#endif
    
    // we will look for minimum absolute distance, so set it to the max to begin
    int64_t approx_dist = std::numeric_limits<int64_t>::max();
    for (const pair<size_t, bool>& oriented_path : shared_paths) {
#ifdef debug_algorithms
        cerr << "[XG] estimating distance with shared path " << oriented_path.first << (oriented_path.second ? "-" : "+") << endl;
#endif
        XGPath& path = *paths[oriented_path.first - 1];
        auto& node_trav_1 = path_dists_1[oriented_path];
        auto& node_trav_2 = path_dists_2[oriented_path];
        
        // the net distance searched between the two points to get to nodes on the path
        int64_t relative_offset = get<2>(node_trav_1) - get<2>(node_trav_2);
        
#ifdef debug_algorithms
        cerr << "[XG] search offset adds up to " << relative_offset << endl;
#endif
        
        memoized_oriented_paths_of_node(get<0>(node_trav_1), local_paths_var1, curr_node_path_occurrences, paths_of_node_memo);
        vector<pair<size_t, vector<pair<size_t, bool>>>>& path_occurrences_1 = *curr_node_path_occurrences;
        memoized_oriented_paths_of_node(get<0>(node_trav_2), local_paths_var2, curr_node_path_occurrences, paths_of_node_memo);
        vector<pair<size_t, vector<pair<size_t, bool>>>>& path_occurrences_2 = *curr_node_path_occurrences;
        
        // get the records corresponding to this shared path
        size_t k = 0;
        while (path_occurrences_1[k].first != oriented_path.first) {
            k++;
        }
        auto& oriented_occurrences_1 = path_occurrences_1[k].second;
        k = 0;
        while (path_occurrences_2[k].first != oriented_path.first) {
            k++;
        }
        auto& oriented_occurrences_2 = path_occurrences_2[k].second;

        // check each pair of node traversals for orientation consistency
        for (size_t i = 0; i < oriented_occurrences_1.size(); i++) {
            // is the first traversal oriented correctly on the path?
            if ((get<1>(node_trav_1) != oriented_occurrences_1[i].second) != oriented_path.second) {
                continue;
            }
            for (size_t j = 0; j < oriented_occurrences_2.size(); j++) {
                // is the second traversal oriented correctly on the path?
                if ((get<1>(node_trav_2) != oriented_occurrences_2[j].second) != oriented_path.second) {
                    continue;
                }
                
                // add up the distance between these nodes in the path plus the distance traversed to them
                // n.b. if the node traversals occur on the reverse strand of the path, the traversal
                // distances represent searches to the opposite sides of the node, so the interval
                // is offset by one node
#ifdef debug_algorithms
                cerr << "[XG] found consistently ordered traversals at ranks " << oriented_occurrences_1[i].first << " and " << oriented_occurrences_2[j].first << endl;
#endif
                int64_t interval_dist = relative_offset;
                if (oriented_path.second) {
                    interval_dist += (path.positions[oriented_occurrences_1[i].first] + node_length(get<0>(node_trav_1))) - (path.positions[oriented_occurrences_2[j].first] + node_length(get<0>(node_trav_2)));
                }
                else {
                    interval_dist += path.positions[oriented_occurrences_2[j].first] - path.positions[oriented_occurrences_1[i].first];
                }
                
#ifdef debug_algorithms
                cerr << "[XG] estimating distance on path " << oriented_path.first << (oriented_path.second ? "-" : "+") << " at " << interval_dist << endl;
#endif
                // find the minimum absolute distance, but retain signing
                if (abs(interval_dist) < abs(approx_dist)) {
                    approx_dist = interval_dist;
                }
            }
        }
    }
    
    return approx_dist;
}
    
vector<tuple<int64_t, bool, size_t>> XG::jump_along_closest_path(int64_t id, bool is_rev, size_t offset, int64_t jump_dist, size_t max_search_dist,
                                                                 unordered_map<id_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo,
                                                                 unordered_map<pair<int64_t, bool>, handle_t>* handle_memo) const {
    
#ifdef debug_algorithms
    cerr << "[XG] jumping " << jump_dist << " from position " << id << (is_rev ? " rev:" : " fwd:") << offset << " with a max search dist of " << max_search_dist << endl;
#endif
    
    vector<tuple<int64_t, bool, size_t>> to_return;
    
    // a local variable to store local results if no memo is provided
    vector<pair<size_t, vector<pair<size_t, bool>>>> local_paths_var;
    // a pointer that will be set to point to valid results of oriented_paths_of_node
    vector<pair<size_t, vector<pair<size_t, bool>>>>* curr_node_path_occurrences;
    
    
    // a local struct for traversals that are ordered by search distance and keep
    // track of whether we're searching to the left or right
    struct Traversal {
        Traversal() {}
        Traversal(int64_t dist, handle_t handle, bool search_left) :
        dist(dist), handle(handle), search_left(search_left) {}
        handle_t handle;
        int64_t dist;
        bool search_left;
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    function<bool(handle_t,int64_t,bool)> look_for_jumpable_paths = [&](handle_t handle, int64_t search_dist, bool search_left) {
        
        bool found_jumpable_path = false;
        
        int64_t trav_id = get_id(handle);
        bool trav_is_rev = get_is_reverse(handle);
        
#ifdef debug_algorithms
        cerr << "[XG] checking for jumpable paths for " << trav_id << (trav_is_rev ? "-" : "+") << " at search dist " << search_dist << " from searching " << (search_left ? "leftwards" : "rightwards") << endl;
#endif
        
        memoized_oriented_paths_of_node(trav_id, local_paths_var, curr_node_path_occurrences, paths_of_node_memo);
        
        for (pair<size_t, vector<pair<size_t, bool>>>& oriented_occurrences : *curr_node_path_occurrences) {
            
            const XGPath& path = *paths[oriented_occurrences.first - 1];
            
            for (const pair<size_t, bool>& occurrence : oriented_occurrences.second) {
#ifdef debug_algorithms
                cerr << "\tnode is on path " << oriented_occurrences.first << " in " << (occurrence.second ? "reverse" : "forward") << " orientation at rank " << occurrence.first << " in path of size " << path.positions.size() << endl;
#endif
                bool path_rev_strand = occurrence.second != trav_is_rev;
                int64_t dist_to_jump;
                if (path_rev_strand && search_left) {
                    dist_to_jump = -jump_dist - search_dist;
                }
                else if (path_rev_strand) {
                    dist_to_jump = -jump_dist + search_dist + get_length(handle);
                }
                else if (search_left) {
                    dist_to_jump = jump_dist + search_dist + get_length(handle);
                }
                else {
                    dist_to_jump = jump_dist - search_dist;
                }
                
                int64_t target_path_pos = path.positions[occurrence.first] + dist_to_jump;
                
#ifdef debug_algorithms
                cerr << "\tjump position should be at offset " << target_path_pos << endl;
#endif
                
                if (target_path_pos >= 0 && target_path_pos < path.offsets.size()) {
                    size_t jump_rank = path.offsets_rank(target_path_pos + 1) - 1;
#ifdef debug_algorithms
                    cerr << "\tthis position is found at path index " << jump_rank << endl;
#endif
                    
                    size_t node_offset = path.positions[jump_rank];
                    int64_t node_id = path.ids[jump_rank];
                    bool node_is_rev = path.directions[jump_rank];
                    
                    size_t offset = path_rev_strand ? node_offset + node_length(node_id) - target_path_pos : target_path_pos - node_offset;
#ifdef debug_algorithms
                    cerr << "\tthe node at this index is " << node_id << (node_is_rev ? "-" : "+") << " at path position " << node_offset << " with additional offset " << offset << endl;
#endif
                    to_return.emplace_back(node_id, node_is_rev != path_rev_strand, offset);
                    
                    found_jumpable_path = true;
                }
            }
        }
        
        return found_jumpable_path;
    };
    
    priority_queue<Traversal> queue;
    unordered_set<handle_t> traversed;
    
    handle_t handle = memoized_get_handle(id, is_rev, handle_memo);
    
    // add in the initial traversals in both directions from the start position
    queue.emplace(offset, handle, true);
    queue.emplace(get_length(handle) - offset, handle, false);
    
    // mark the initial node as traversed
    traversed.insert(handle);
    
    bool found_jumpable_path = look_for_jumpable_paths(handle, -offset, false);
    while (!queue.empty() && !found_jumpable_path) {
        // get the queue that has the next shortest path
        
        Traversal trav = queue.top();
        queue.pop();
        
#ifdef debug_algorithms
        cerr << "[XG] traversing " << get_id(trav.handle) << (get_is_reverse(trav.handle) ? "-" : "+") << " in " << (trav.search_left ? "leftward" : "rightward") << " direction at distance " << trav.dist << endl;
#endif
        
        function<bool(const handle_t& next)> check_next = [&](const handle_t& next) {
#ifdef debug_algorithms
            cerr << "\tfollowing edge to " << get_id(next) << (get_is_reverse(next) ? "-" : "+") << " at dist " << trav.dist << endl;
#endif
            
            if (!traversed.count(next)) {
                found_jumpable_path = look_for_jumpable_paths(next, trav.dist, trav.search_left);
                
                int64_t dist_thru = trav.dist + get_length(next);
                
                if (dist_thru <= (int64_t) max_search_dist) {
                    queue.emplace(dist_thru, next, trav.search_left);
                }
                traversed.emplace(next);
            }
            return !found_jumpable_path;
        };
        
        follow_edges(trav.handle, trav.search_left, check_next);
    }
    
    return std::move(to_return);
}

void XG::for_path_range(const string& name, int64_t start, int64_t stop,
                        function<void(int64_t)> lambda, bool is_rev) const {

    // what is the node at the start, and at the end
    auto& path = *paths[path_rank(name)-1];
    size_t plen = path.offsets.size();
    if (start > plen) return; // no overlap with path
    // careful not to exceed the path length
    if (stop >= plen) stop = plen-1;
    if (is_rev) {
        start = plen - start;
        stop = plen - stop;
    }
    size_t pr1 = path.offsets_rank(start+1)-1;
    size_t pr2 = path.offsets_rank(stop+1)-1;

    // Grab the IDs visited in order along the path
    auto& pi_wt = path.ids;
    for (size_t i = pr1; i <= pr2; ++i) {
        // For all the visits along this section of path, grab the node being visited and all its edges.
        int64_t id = pi_wt[i];
        lambda(id);
    }
}

void XG::get_path_range(const string& name, int64_t start, int64_t stop, Graph& g, bool is_rev) const {

    set<int64_t> nodes;
    set<pair<side_t, side_t> > edges;

    for_path_range(name, start, stop, [&](int64_t id) {
            nodes.insert(id);
            for (auto& e : edges_from(id)) {
                edges.insert(make_pair(make_side(e.from(), e.from_start()), make_side(e.to(), e.to_end())));
            }
            for (auto& e : edges_to(id)) {
                edges.insert(make_pair(make_side(e.from(), e.from_start()), make_side(e.to(), e.to_end())));
            }
        }, is_rev);

    for (auto& n : nodes) {
        *g.add_node() = node(n);
    }
    
    map<string, Path*> local_paths;
    for (auto& n : nodes) {
        for (auto& m : node_mappings(n)) {
            if (local_paths.find(m.first) == local_paths.end()) {
                Path* p = g.add_path();
                local_paths[m.first] = p;
                p->set_name(m.first);
            }
            Path* new_path = local_paths[m.first];
            // TODO output mapping direction
            //if () { m.second.is_reverse(true); }
            for (auto& n : m.second) {
                *new_path->add_mapping() = n;
            }
        }
    }
    for (auto& e : edges) {
        Edge edge;
        edge.set_from(side_id(e.first));
        edge.set_from_start(!side_is_end(e.first));
        edge.set_to(side_id(e.second));
        edge.set_to_end(side_is_end(e.second));
        *g.add_edge() = edge;
    }
}

size_t XG::node_occs_in_path(int64_t id, const string& name) const {
    return node_occs_in_path(id, path_rank(name));
}

size_t XG::node_occs_in_path(int64_t id, size_t rank) const {
    size_t p = rank-1;
    auto& pi_wt = paths[p]->ids;
    return pi_wt.rank(pi_wt.size(), id);
}

vector<size_t> XG::node_ranks_in_path(int64_t id, const string& name) const {
    return node_ranks_in_path(id, path_rank(name));
}

vector<size_t> XG::node_ranks_in_path(int64_t id, size_t rank) const {
    vector<size_t> ranks;
    size_t p = rank-1;
    size_t occs = node_occs_in_path(id, rank);
    for (size_t i = 1; i <= occs; ++i) {
        ranks.push_back(paths[p]->ids.select(i, id));
    }
    return ranks;
}

vector<size_t> XG::position_in_path(int64_t id, const string& name) const {
    return position_in_path(id, path_rank(name));
}

vector<size_t> XG::position_in_path(int64_t id, size_t rank) const {
    auto& path = *paths[rank-1];
    vector<size_t> pos_in_path;
    for (auto i : node_ranks_in_path(id, rank)) {
        pos_in_path.push_back(path.positions[i]);
    }
    return pos_in_path;
}

map<string, vector<size_t> > XG::position_in_paths(int64_t id, bool is_rev, size_t offset) const {
    map<string, vector<size_t> > positions;
    for (auto& prank : paths_of_node(id)) {
        auto& path = *paths[prank-1];
        auto& pos_in_path = positions[path_name(prank)];
        for (auto i : node_ranks_in_path(id, prank)) {
            size_t pos = offset + (is_rev ?
                                   path_length(prank) - path.positions[i] - node_length(id)
                                   : path.positions[i]);
            pos_in_path.push_back(pos);
        }
    }
    return positions;
}

map<string, vector<size_t> > XG::distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                                   int64_t id2, bool is_rev2, size_t offset2) const {
    auto pos1 = position_in_paths(id1, is_rev1, offset1);
    auto pos2 = position_in_paths(id2, is_rev2, offset2);
    map<string, vector<size_t> > dist;
    // position in a path is undefined in inversion
    if (is_rev1 != is_rev2) {
        return dist;
    }
    for (auto& c1 : pos1) {
        auto c2 = pos2.find(c1.first);
        if (c2 != pos2.end()) {
            auto& d = dist[c1.first];
            // distances are a cross of the points
            for (auto o1 : c1.second) {
                for (auto o2 : c2->second) {
                    d.push_back(o1-o2);
                }
            }
        }
    }
    return dist;
}

int64_t XG::min_distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                  int64_t id2, bool is_rev2, size_t offset2) const {
    auto dist = distance_in_paths(id1, is_rev1, offset1,
                                  id2, is_rev2, offset2);
    size_t min_dist = std::numeric_limits<size_t>::max();
    for (auto& c : dist) {
        for (auto& o : c.second) {
            if (o <  min_dist) {
                min_dist = o;
            }
        }
    }
    return min_dist;
}

int64_t XG::node_at_path_position(const string& name, size_t pos) const {
    size_t p = path_rank(name)-1;
    return paths[p]->ids[paths[p]->offsets_rank(pos+1)-1];
}

Mapping XG::mapping_at_path_position(const string& name, size_t pos) const {
    size_t p = path_rank(name)-1;
    return paths[p]->mapping(paths[p]->offsets_rank(pos+1)-1);
}

size_t XG::node_start_at_path_position(const string& name, size_t pos) const {
    size_t p = path_rank(name)-1;
    size_t position_rank = paths[p]->offsets_rank(pos+1);
    return paths[p]->offsets_select(position_rank);
}

Alignment XG::target_alignment(const string& name, size_t pos1, size_t pos2, const string& feature) const {
    Alignment aln;
    const XGPath& path = *paths[path_rank(name)-1];
    size_t first_node_start = path.offsets_select(path.offsets_rank(pos1+1));
    int64_t trim_start = pos1 - first_node_start;
    {
        Mapping* first_mapping = aln.mutable_path()->add_mapping();
        *first_mapping = mapping_at_path_position(name, pos1);
        Edit* e = first_mapping->add_edit();
        e->set_to_length(node_length(first_mapping->position().node_id())-trim_start);
        e->set_from_length(e->to_length());
    }
    // get p to point to the next step (or past it, if we're a feature on a single node)
    int64_t p = (int64_t)pos1 + node_length(aln.path().mapping(0).position().node_id()) - trim_start;
    while (p < pos2) {
        Mapping m = mapping_at_path_position(name, p);
        Edit* e = m.add_edit();
        e->set_to_length(node_length(m.position().node_id()));
        e->set_from_length(e->to_length());
        *aln.mutable_path()->add_mapping() = m;
        p += mapping_from_length(aln.path().mapping(aln.path().mapping_size()-1));
    }
    // trim to the target
    int64_t trim_end = p - pos2;
    if (trim_start) {
        *aln.mutable_path() = cut_path(aln.path(), trim_start).second;
    }
    if (trim_end) {
        *aln.mutable_path() = cut_path(aln.path(), path_from_length(aln.path()) - trim_end).first;
    }
    aln.set_name(feature);
    return aln;
}

Mapping new_mapping(const string& name, int64_t id, size_t rank, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.mutable_position()->set_is_reverse(is_reverse);
    m.set_rank(rank);
    return m;
}

void parse_region(const string& target, string& name, int64_t& start, int64_t& end) {
    start = -1;
    end = -1;
    size_t foundFirstColon = target.find(":");
    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        name = target;
    } else {
        name = target.substr(0, foundFirstColon);
	    size_t foundRangeDash = target.find("-", foundFirstColon);
        if (foundRangeDash == string::npos) {
            start = atoi(target.substr(foundFirstColon + 1).c_str());
            end = start;
        } else {
            start = atoi(target.substr(foundFirstColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
            end = atoi(target.substr(foundRangeDash + 1).c_str());
        }
    }
}

void to_text(ostream& out, Graph& graph) {
    out << "H" << "\t" << "HVN:Z:1.0" << endl;
    for (size_t i = 0; i < graph.node_size(); ++i) {
        auto& node = graph.node(i);
        out << "S" << "\t" << node.id() << "\t" << node.sequence() << endl;
    }
    for (size_t i = 0; i < graph.path_size(); ++i) {
        auto& path = graph.path(i);
        for (size_t j = 0; j < path.mapping_size(); ++j) {
            auto& mapping = path.mapping(i);
            string orientation = mapping.position().is_reverse() ? "-" : "+";
            out << "P" << "\t" << mapping.position().node_id() << "\t" << path.name() << "\t"
                << mapping.rank() << "\t" << orientation << "\n";
        }
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        auto& edge = graph.edge(i);
        out << "L" << "\t" << edge.from() << "\t"
            << (edge.from_start() ? "-" : "+") << "\t"
            << edge.to() << "\t"
            << (edge.to_end() ? "-" : "+") << endl;
    }
}

int64_t XG::where_to(int64_t current_side, int64_t visit_offset, int64_t new_side) const {
    // Given that we were at visit_offset on the current side, where will we be
    // on the new side? 
    
    // Work out where we're going as a node and orientation
    int64_t new_node_id = rank_to_id(new_side / 2);
    bool new_node_is_reverse = new_side % 2;
    
    // Work out what node and orientation we came from
    int64_t old_node_id = rank_to_id(current_side / 2);
    bool old_node_is_reverse = current_side % 2;
    
    // Work out what edges are going into the place we're going into.
    vector<Edge> edges = new_node_is_reverse ? edges_on_end(new_node_id) : edges_on_start(new_node_id);
    
    // Look at the edges we could have taken next
    vector<Edge> edges_out = old_node_is_reverse ? edges_on_start(old_node_id) : edges_on_end(old_node_id);
    
    return where_to(current_side, visit_offset, new_side, edges, edges_out);
}

int64_t XG::node_height(XG::ThreadMapping node) const {
  return h_civ[node_graph_idx(node.node_id) * 2 + node.is_reverse];
}

int64_t XG::where_to(int64_t current_side, int64_t visit_offset, int64_t new_side, vector<Edge>& edges_into_new, vector<Edge>& edges_out_of_old) const {
    // Given that we were at visit_offset on the current side, where will we be
    // on the new side?

    // What will the new visit offset be?
    int64_t new_visit_offset = 0;

    // Work out where we're going as a node and orientation
    int64_t new_node_id = rank_to_id(new_side / 2);
    bool new_node_is_reverse = new_side % 2;

    // Work out what edges are going into the place we're going into.

    // Work out what node and orientation we came from
    int64_t old_node_id = rank_to_id(current_side / 2);
    bool old_node_is_reverse = current_side % 2;

    // What edge are we following
    Edge edge_taken = make_edge(old_node_id, old_node_is_reverse, new_node_id, new_node_is_reverse);

    // Make sure we find it
    bool edge_found = false;

    for(auto& edge : edges_into_new) {
        // Look at every edge in order.

        if(edges_equivalent(edge, edge_taken)) {
            // If we found the edge we're taking, break.
            edge_found = true;
            break;
        }

        // Otherwise add in the threads on this edge to the offset

        // Get the orientation number for this edge (which is *not* the edge we
        // are following and so doesn't involve old_node_anything). We only
        // treat it as reverse if the reverse direction is the only direction we
        // can take to get here.
        int64_t edge_orientation_number = (edge_graph_idx(edge) * 2) + arrive_by_reverse(edge, new_node_id, new_node_is_reverse);

        int64_t contribution = h_civ[edge_orientation_number];
#ifdef VERBOSE_DEBUG
        cerr << contribution << " (from prev edge " << edge.from() << (edge.from_start() ? "L" : "R") << "-" << edge.to() << (edge.to_end() ? "R" : "L") << " at " << edge_orientation_number << ") + ";
#endif
        new_visit_offset += contribution;
    }

    assert(edge_found);

    // What edge out of all the edges we can take are we taking?
    int64_t edge_taken_index = -1;

    // Look at the edges we could have taken next

    for(int64_t i = 0; i < edges_out_of_old.size(); i++) {
        if(edges_equivalent(edges_out_of_old[i], edge_taken)) {
            // i is the index of the edge we took, of the edges available to us.
            edge_taken_index = i;
            break;
        }
    }

    assert(edge_taken_index != -1);

    // Get the rank in B_s[] for our current side of our visit offset among
    // B_s[] entries pointing to the new node and add that in. Make sure to +2
    // to account for the nulls and separators.
    int64_t contribution = bs_rank(current_side, visit_offset, edge_taken_index + 2);
#ifdef VERBOSE_DEBUG
    cerr << contribution << " (via this edge) + ";
#endif
    new_visit_offset += contribution;

    // Get the number of threads starting at the new side and add that in.
    contribution = ts_civ[new_side];
#ifdef VERBOSE_DEBUG
    cerr << contribution << " (starting here) = ";
#endif
    new_visit_offset += contribution;
#ifdef VERBOSE_DEBUG
    cerr << new_visit_offset << endl;
#endif

    // Now we know where it actually ends up: after all the threads that start,
    // all the threads that come in via earlier edges, and all the previous
    // threads going there that come via this edge.
    return new_visit_offset;
}

XG::thread_t XG::extract_thread(xg::XG::ThreadMapping node, int64_t offset = 0, int64_t max_length = 0) {
  thread_t path;
  int64_t side = (node.node_id)*2 + node.is_reverse;
  bool continue_search = true;
  while(continue_search) {
    xg::XG::ThreadMapping m = {rank_to_id(side / 2), (bool) (side % 2)};
    path.push_back(m);
    // Work out where we go:
    // What edge of the available edges do we take?
    int64_t edge_index = bs_get(side, offset);
    assert(edge_index != BS_SEPARATOR);
    if(edge_index == BS_NULL) {
        // Path ends here.
        break;
    } else {
        // If we've got an edge, convert to an actual edge index
        edge_index -= 2;
    }
    // We also should not have negative edges.
    assert(edge_index >= 0);
    // Look at the edges we could have taken next
    vector<Edge> edges_out = side % 2 ? edges_on_start(rank_to_id(side / 2)) :
          edges_on_end(rank_to_id(side / 2));
    assert(edge_index < edges_out.size());
    Edge& taken = edges_out[edge_index];
    // Follow the edge
    int64_t other_node = taken.from() == rank_to_id(side / 2) ? taken.to() :
          taken.from();
    bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
    // Get the side
    int64_t other_side = id_to_rank(other_node) * 2 + other_orientation;
    // Go there with where_to
    offset = where_to(side, offset, other_side);
    side = other_side;
    // Continue the process from this new side
    if(max_length != 0) {
      if(path.size() >= max_length) {
        continue_search = false;
      }
    }
  }
  return path;
}

void XG::insert_threads_into_dag(const vector<thread_t>& t, const vector<string>& names) {

    util::assign(h_iv, int_vector<>(g_iv.size() * 2, 0));
    util::assign(ts_iv, int_vector<>((node_count + 1) * 2, 0));
    // Store the names
    for (auto& name : names) {
        names_str.append("$" + name);
    }
#ifdef VERBOSE_DEBUG
    cerr << "Compressing thread names..." << endl;
#endif
    tn_bake();

    // store the sides in order of their addition to the threads
    int_vector<> sides_ordered_by_thread_id(t.size()*2); // fwd and reverse

    // store the start+offset for each thread and its reverse complement
    int_vector<> tin_iv(t.size()*2+2);
    int_vector<> tio_iv(t.size()*2+2);
    int thread_count = 0;
    
    auto emit_destinations = [&](int64_t node_id, bool is_reverse, vector<size_t> destinations) {
        // We have to take this destination vector and store it in whatever B_s
        // storage we are using.
        
        int64_t node_side = id_to_rank(node_id) * 2 + is_reverse;
        
        // Copy all the destinations into the succinct B_s storage
        bs_set(node_side, destinations);
        
        // Set the number of total visits to this side.
        h_iv[node_graph_idx(node_id) * 2 + is_reverse] = destinations.size();
    
#ifdef VERBOSE_DEBUG
        cerr << "Found " << destinations.size() << " visits total to node " << node_id << (is_reverse ? "-" : "+") << endl;
#endif
    };
    
    auto emit_edge_traversal = [&](int64_t node_id, bool from_start, int64_t next_node_id, bool to_end) {
        // We have to store the fact that we traversed this edge in the specified direction in our succinct storage. 
        
        // Find the edge as it actually appears in the graph.
        Edge canonical = canonicalize(make_edge(node_id, from_start, next_node_id, to_end));
        
        // We're departing along this edge, so our orientation cares about
        // whether we have to take the edge forward or backward when departing.
        auto edge_idx = edge_graph_idx(canonical);
        assert(edge_idx != numeric_limits<size_t>::max()); // We must actually have the edge
        int64_t edge_orientation_number = edge_idx * 2 +
            depart_by_reverse(canonical, node_id, from_start);

#ifdef VERBOSE_DEBUG
        cerr << "We need to add 1 to the traversals of oriented edge " << edge_orientation_number << endl;
#endif
                
        // Increment the count for the edge
        h_iv[edge_orientation_number]++; 
    };
    
    auto emit_thread_start = [&](int64_t node_id, bool is_reverse) {
        // Record that an (orientation of) a thread starts at this node in this
        // orientation. We have to update our thread start succinct data
        // structure.
        
        int64_t node_side = id_rev_to_side(node_id, is_reverse);
        
        // Say we start one more thread on this side.
        ts_iv[node_side]++;
        // Record the side
        sides_ordered_by_thread_id[thread_count] = node_side;
        // Sum the number of threads
        thread_count++;

#ifdef VERBOSE_DEBUG
        cerr << "A thread starts at " << node_side << " for " <<  node_id << (is_reverse ? "-" : "+") << endl;
#endif
        
    };

    // We want to go through and insert running forward through the DAG, and
    // then again backward through the DAG.
    auto insert_in_direction = [&](bool insert_reverse) {

        // First sort out the thread numbers by the node they start at.
        // We know all the threads go the same direction through each node.
        map<int64_t, list<size_t>> thread_numbers_by_start_node;
        
        for(size_t i = 0; i < t.size(); i++) {
            if(t[i].size() > 0) {
                // Do we start with the first or last mapping in the thread?
                size_t thread_start = insert_reverse ? t[i].size() - 1 : 0;
                auto& mapping = t[i][thread_start];
                thread_numbers_by_start_node[mapping.node_id].push_back(i);
                // we know the mapping node id and rank
                // so we can construct the start position for this entity
                int k = 2*(i+1) + insert_reverse;
                tin_iv[k] = mapping.node_id;
                // nb: the rank of this thread among those starting at this side
                tio_iv[k] = thread_numbers_by_start_node[mapping.node_id].size()-1;
                // Say a thread starts here, going in the orientation determined
                // by how the node is visited and how we're traversing the path.
                emit_thread_start(mapping.node_id, mapping.is_reverse != insert_reverse);
            }
        }
        
        // We have this message-passing architecture, where we send groups of
        // threads along edges to destination nodes. This records, by edge rank of
        // the traversed edge (with 0 meaning starting there), the group of threads
        // coming in along that edge, and the offset in each thread that the visit
        // to the node is at. These are messages passed along the edge from the
        // earlier node to the later node (since we know threads follow a DAG).
        map<size_t, list<pair<size_t, size_t>>> edge_to_ordered_threads;
        
        for(size_t node_rank = (insert_reverse ? max_node_rank() : 1);
            node_rank != (insert_reverse ? 0 : max_node_rank() + 1);
            node_rank += (insert_reverse ? -1 : 1)) {
            // Then we start at the first node in the DAG
            
            int64_t node_id = rank_to_id(node_rank);
            
#ifdef VERBOSE_DEBUG
            if(node_id % 10000 == 1) {
                cerr << "Processing node " << node_id << endl;
            }
#endif
            
            // We order the thread visits starting there, and then all the threads
            // coming in from other places, ordered by edge traversed.
            // Stores a pair of thread number and mapping index in the thread.
            list<pair<size_t, size_t>> threads_visiting;
            
            if(thread_numbers_by_start_node.count(node_id)) {
                // Grab the threads starting here
                for(size_t thread_number : thread_numbers_by_start_node.at(node_id)) {
                    // For every thread that starts here, say it visits here
                    // with its first mapping (0 for forward inserts, last one
                    // for reverse inserts).
                    threads_visiting.emplace_back(thread_number, insert_reverse ? t[thread_number].size() - 1 : 0);
                }
                thread_numbers_by_start_node.erase(node_id);
            }
            
            
            for(Edge& in_edge : edges_of(node_id)) {
                // Look at all the edges on the node. Messages will only exist on
                // the incoming ones.
                auto edge_idx = edge_graph_idx(in_edge);
                if(edge_to_ordered_threads.count(edge_idx)) {
                    // We have messages coming along this edge on our start
                    
                    // These threads come in next. Splice them in because they
                    // already have the right mapping indices.
                    threads_visiting.splice(threads_visiting.end(), edge_to_ordered_threads[edge_idx]);
                    edge_to_ordered_threads.erase(edge_idx);
                }
            }
            
            if(threads_visiting.empty()) {
                // Nothing visits here, so there's no cool succinct data structures
                // to generate.
                continue;
            }
            
            
            // Some threads visit here! Determine our orientation from the first
            // and assume it applies to all threads visiting us.
            auto& first_visit = threads_visiting.front();
            bool node_is_reverse = t[first_visit.first][first_visit.second].is_reverse;
            // When we're inserting threads backwards, we need to treat forward
            // nodes as reverse and visa versa, to leave the correct side.
            node_is_reverse = node_is_reverse != insert_reverse;
            
            // Now we have all the threads coming through this node, and we know
            // which way they are going.

            // Make a map from outgoing edge rank to the B_s array number (0 for
            // stop here, 1 reserved as a separator, and 2 through n corresponding
            // to outgoing edges in order) for this node's outgoing side.
            map<size_t, size_t> edge_idx_to_local_edge_number;
            auto outgoing_edges = node_is_reverse ? edges_on_start(node_id) : edges_on_end(node_id);
            for(size_t i = 0; i < outgoing_edges.size(); i++) {
                size_t edge_idx = edge_graph_idx(outgoing_edges[i]);
                edge_idx_to_local_edge_number[edge_idx] = i + 2;
            }
            
            // Make a vector we'll fill in with all the B array values (0 for stop,
            // 2 + edge number for outgoing edge)
            vector<size_t> destinations;
            
            for(auto& visit : threads_visiting) {
                // Now go through all the path visits, fill in the edge numbers (or 0
                // for stop) they go to, and stick visits to the next mappings in the
                // correct message lists.
                if(insert_reverse ? (visit.second != 0) : (visit.second + 1 < t[visit.first].size())) {
                    // This visit continues on from here
                    
                    // Make a visit to the next mapping on the path
                    auto next_visit = visit;
                    next_visit.second += insert_reverse ? -1 : 1;
                    
                    // Work out what node that is, and what orientation
                    auto& next_mapping = t[next_visit.first][next_visit.second];
                    int64_t next_node_id = next_mapping.node_id;
                    bool next_is_reverse = next_mapping.is_reverse != insert_reverse;
                    
                    // Figure out the rank of the edge we need to take to get
                    // there. Note that we need to make sure we can handle going
                    // forward and backward over edges.
                    Edge canonical = canonicalize(make_edge(node_id, node_is_reverse, next_node_id, next_is_reverse));
                    size_t next_edge_idx = edge_graph_idx(canonical);
                    
                    // Look up what local edge number that edge gets and say we follow it.
                    destinations.push_back(edge_idx_to_local_edge_number.at(next_edge_idx));
                    
                    // Send the new mapping along the edge after all the other ones
                    // we've sent along the edge
                    edge_to_ordered_threads[next_edge_idx].push_back(next_visit);
                    
                    // Say we traverse an edge going from this node in this
                    // orientation to that node in that orientation.
                    emit_edge_traversal(node_id, node_is_reverse, next_node_id, next_is_reverse);
                    
                } else {
                    // This visit ends here
#ifdef VERBOSE_DEBUG
                    if(insert_reverse) {
                        cerr << "A thread ends here because " << visit.second << " is 0 " << endl;
                    } else {
                        cerr << "A thread ends here because " << visit.second + 1 << " is >= " << t[visit.first].size() << endl;
                    }
#endif
                    destinations.push_back(BS_NULL);
                }
            }
            
            // Emit the destinations array for the node. Store it in whatever
            // sort of succinct storage we are using...
            // We need to send along the side (false for left, true for right)
            emit_destinations(node_id, node_is_reverse, destinations);
            
            // We repeat through all nodes until done.
        }
    
        // OK now we have gone through the whole set of everything and inserted
        // in this direction.
    };
    
    // Actually call the inserts
#ifdef VERBOSE_DEBUG
    cerr << "Inserting threads forwards..." << endl;
#endif
    insert_in_direction(false);
#ifdef VERBOSE_DEBUG
    cerr << "Inserting threads backwards..." << endl;
#endif
    insert_in_direction(true);
    
    // Actually build the B_s arrays for rank and select.
#ifdef VERBOSE_DEBUG
    cerr << "Creating final compressed array..." << endl;
#endif
    bs_bake();
    // compress the visit counts
    util::assign(h_civ, h_iv);
    util::assign(ts_civ, ts_iv);
    // compress the starts for the threads
    util::assign(tin_civ, int_vector<>(tin_iv));
    util::assign(tio_civ, int_vector<>(tio_iv));
    util::bit_compress(sides_ordered_by_thread_id);
    // and build up the side wt
    construct_im(side_thread_wt, sides_ordered_by_thread_id);
}

void XG::insert_thread(const thread_t& t, const string& name) {
    // We're going to insert this thread
    
    auto insert_thread_forward = [&](const thread_t& thread) {
    
#ifdef VERBOSE_DEBUG
        cerr << "Inserting thread with " << thread.size() << " mappings" << endl;
#endif
        // Where does the current visit fall on its node? On the first node we
        // arbitrarily decide to be first of all the threads starting there.
        // TODO: Make sure that we actually end up ordering starts based on path
        // name or something later.
        int64_t visit_offset = 0;
        for(size_t i = 0; i < thread.size(); i++) {
            // For each visit to a node...
        
            // What side are we visiting?
            int64_t node_id = thread[i].node_id;
            bool node_is_reverse = thread[i].is_reverse;
            int64_t node_side = id_to_rank(node_id) * 2 + node_is_reverse;
            
#ifdef VERBOSE_DEBUG
            cerr << "Visit side " << node_side << endl;
#endif

            // Where are we going next?
            
            if(i == thread.size() - 1) {
                // This is the last visit. Send us off to null
                
#ifdef VERBOSE_DEBUG
                cerr << "End the thread." << endl;
#endif
                // Stick a new entry in the B array at the place where it belongs.
                bs_insert(node_side, visit_offset, BS_NULL);
            } else {
                // This is not the last visit. Send us off to the next place, and update the count on the edge.
                
                // Work out where we're actually going next
                int64_t next_id = thread[i + 1].node_id;
                bool next_is_reverse = thread[i + 1].is_reverse;
                int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;

                // What edge do we take to get there? We're going to search for
                // the actual edge articulated in the official orientation.
                Edge edge_wanted = make_edge(node_id, node_is_reverse, next_id, next_is_reverse);
                
                // And what is its index in the edges we can take?
                int64_t edge_taken_index = -1;
                
                // Look at the edges we could have taken
                vector<Edge> edges_out = node_is_reverse ? edges_on_start(node_id) : edges_on_end(node_id);
                for(int64_t j = 0; j < edges_out.size(); j++) {
                    if(edge_taken_index == -1 && edges_equivalent(edges_out[j], edge_wanted)) {
                        // i is the index of the edge we took, of the edges available to us.
                        edge_taken_index = j;
                    }
#ifdef VERBOSE_DEBUG
                    cerr << "Edge #" << j << ": "  << edges_out[j].from() << (edges_out[j].from_start() ? "L" : "R") << "-" << edges_out[j].to() << (edges_out[j].to_end() ? "R" : "L") << endl;
                    
                    // Work out what its number is for the orientation it goes
                    // out in. We know we have this edge, so we just see if we
                    // have to depart in reveres and if so take the edge in
                    // reverse.
                    int64_t edge_orientation_number = edge_graph_idx(edges_out[j]) * 2 + depart_by_reverse(edges_out[j], node_id, node_is_reverse);
                    
                    cerr << "\tOrientation rank: " << edge_orientation_number << endl;
#endif
                }
                
                if(edge_taken_index == -1) {
                    cerr << "[xg] error: step " << i << " of thread: "
                        << edge_wanted.from() << (edge_wanted.from_start() ? "L" : "R") << "-" 
                        << edge_wanted.to() << (edge_wanted.to_end() ? "R" : "L") << " does not exist" << endl;
                    cerr << "[xg] error: Possibilities: ";
                    for(Edge& edge_out : edges_out) {
                        cerr << edge_out.from() << (edge_out.from_start() ? "L" : "R") << "-" 
                            << edge_out.to() << (edge_out.to_end() ? "R" : "L") << " ";
                    }
                    cerr << endl;
                    cerr << "[xg] error: On start: ";
                    for(Edge& edge_out : edges_on_start(node_id)) {
                        cerr << edge_out.from() << (edge_out.from_start() ? "L" : "R") << "-" 
                            << edge_out.to() << (edge_out.to_end() ? "R" : "L") << " ";
                    }
                    cerr << endl;
                    cerr << "[xg] error: On end: ";
                    for(Edge& edge_out : edges_on_end(node_id)) {
                        cerr << edge_out.from() << (edge_out.from_start() ? "L" : "R") << "-" 
                            << edge_out.to() << (edge_out.to_end() ? "R" : "L") << " ";
                    }
                    cerr << endl;
                    cerr << "[xg] error: Both: ";
                    for(Edge& edge_out : edges_of(node_id)) {
                        cerr << edge_out.from() << (edge_out.from_start() ? "L" : "R") << "-" 
                            << edge_out.to() << (edge_out.to_end() ? "R" : "L") << " ";
                    }
                    cerr << endl;
                    cerr << "[xg] error: has_edge: " << has_edge(edge_wanted.from(), edge_wanted.from_start(), edge_wanted.to(), edge_wanted.to_end()) << endl;
                    assert(false);
                }
                
                assert(edge_taken_index != -1);

#ifdef VERBOSE_DEBUG
                cerr << "Proceed to " << next_side << " via edge #" << edge_taken_index << "/" << edges_out.size() << endl;
#endif
            
                // Make a nice reference to the edge we're taking in its real orientation.
                auto& edge_taken = edges_out[edge_taken_index];
            
                // Stick a new entry in the B array at the place where it belongs.
                // Make sure to +2 to leave room in the number space for the separators and null destinations.
                bs_insert(node_side, visit_offset, edge_taken_index + 2);
                
                // Update the usage count for the edge going form here to the next node
                // Make sure that edge storage direction is correct.
                // Which orientation of an edge are we crossing?
                int64_t edge_orientation_number = edge_graph_idx(edge_taken) * 2 + depart_by_reverse(edge_taken, node_id, node_is_reverse);
                
#ifdef VERBOSE_DEBUG
                cerr << "We need to add 1 to the traversals of oriented edge " << edge_orientation_number << endl;
#endif
                
                // Increment the count for the edge
                h_iv[edge_orientation_number]++;
                
#ifdef VERBOSE_DEBUG
                cerr << "h = [";
                for(int64_t k = 0; k < h_iv.size(); k++) {
                    cerr << h_iv[k] << ", ";
                }
                cerr << "]" << endl;
#endif
                
                // Now where do we go to on the next visit?
                visit_offset = where_to(node_side, visit_offset, next_side);
                
#ifdef VERBOSE_DEBUG

                cerr << "Offset " << visit_offset << " at side " << next_side << endl;
#endif
                
            }
            
#ifdef VERBOSE_DEBUG
            
            cerr << "Node " << node_id << " orientation " << node_is_reverse <<
                " has rank " <<
                (node_graph_idx(node_id) * 2 + node_is_reverse) <<
                " of " << h_civ.size() << endl;
#endif
            
            // Increment the usage count of the node in this orientation
            h_iv[node_graph_idx(node_id) * 2 + node_is_reverse]++;
            
            if(i == 0) {
                // The thread starts here
                ts_iv[node_side]++;
            }
        }
        
    };
    
    // We need a simple reverse that works only for perfect match paths
    auto simple_reverse = [&](const thread_t& thread) {
        // Make a reversed version
        thread_t reversed;
        
        // TODO: give it a reversed name or something
        
        for(size_t i = thread.size(); i != (size_t) -1; i--) { 
            // Copy the mappings from back to front, flipping the is_reverse on their positions.
            ThreadMapping reversing = thread[thread.size() - 1 - i];
            reversing.is_reverse = !reversing.is_reverse;
            reversed.push_back(reversing);
        }
        
        return reversed;        
    };
    
    // Insert forward
    insert_thread_forward(t);
    
    // Insert reverse
    insert_thread_forward(simple_reverse(t));
    
    // Save the name
    names_str.append("$" + name);
}

auto XG::extract_threads_matching(const string& pattern, bool reverse) const -> map<string, list<thread_t>> {

    map<string, list<thread_t> > found;

    // get the set of threads that match the given pattern
    // for each thread, get its start position
    vector<int64_t> threads = threads_named_starting(pattern);
    for (auto& id : threads) {
        // For every thread starting there
        // make a new thread
        thread_t path;

        // Get the name of this thread
        string path_name = thread_name(id);

        // Start the side at i and the offset at j
        auto p = thread_start(id, reverse);
        int64_t side = id_rev_to_side(p.first, reverse);
        int64_t offset = p.second;

        while(true) {
                
            // Unpack the side into a node traversal
            ThreadMapping m = {rank_to_id(side/2), (bool) (side % 2)};

            // Add the mapping to the thread
            path.push_back(m);
                
            // Work out where we go
                
            // What edge of the available edges do we take?
            int64_t edge_index = bs_get(side, offset);
                
            // If we find a separator, we're very broken.
            assert(edge_index != BS_SEPARATOR);
                
            if(edge_index == BS_NULL) {
                // Path ends here.
                break;
            } else {
                // Convert to an actual edge index
                edge_index -= 2;
            }
                
            // We also should not have negative edges.
            assert(edge_index >= 0);
                
            // Look at the edges we could have taken next
            vector<Edge> edges_out = side % 2 ? edges_on_start(rank_to_id(side / 2)) : edges_on_end(rank_to_id(side / 2));
                
            assert(edge_index < edges_out.size());
                
            Edge& taken = edges_out[edge_index];
                
            // Follow the edge
            int64_t other_node = taken.from() == rank_to_id(side / 2) ? taken.to() : taken.from();
            bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
                
            // Get the side 
            int64_t other_side = id_to_rank(other_node) * 2 + other_orientation;
                
            // Go there with where_to
            offset = where_to(side, offset, other_side);
            side = other_side;
    
        }
            
        found[path_name].push_back(path);
            
    }
    
    return found;
    
}

auto XG::extract_threads(bool extract_reverse) const -> map<string, list<thread_t>> {

    // Fill in a map of lists of paths found by name
    map<string, list<thread_t> > found;

#ifdef VERBOSE_DEBUG
    cerr << "Extracting threads" << endl;
#endif
    // We consider "reverse" threads to be those that start on the higher-
    // numbered sides of their nodes.
    // We know sides 0 and 1 are unused, so the smallest side is 2.
    int64_t begin = !extract_reverse ? 2 : 3;
    int64_t end = !extract_reverse ? ts_civ.size()-1 : ts_civ.size();
    for(int64_t i = begin; i < end; i+=2) {
        // For every other real side
    
#ifdef VERBOSE_DEBUG
        cerr << ts_civ[i] << " threads start at side " << i << endl;
#endif
        // For each side
        if(ts_civ[i] == 0) {
            // Skip it if no threads start at it
            continue;
        }

        for(int64_t j = 0; j < ts_civ[i]; j++) {
            // For every thread starting there
      
#ifdef debug    
            cerr << "Extracting thread " << j << endl;
#endif
            
            // make a new thread
            thread_t path;
            
            // Start the side at i and the offset at j
            int64_t side = i;
            int64_t offset = j;

            // Get the name of this thread
            string path_name = thread_name(thread_starting_at(side, offset));

            while(true) {
                
                // Unpack the side into a node traversal
                ThreadMapping m = {rank_to_id(side / 2), (bool) (side % 2)};
                
                // Add the mapping to the thread
                path.push_back(m);
                
#ifdef VERBOSE_DEBUG
                cerr << "At side " << side << endl;
                
#endif
                // Work out where we go
                
                // What edge of the available edges do we take?
                int64_t edge_index = bs_get(side, offset);
                
                // If we find a separator, we're very broken.
                assert(edge_index != BS_SEPARATOR);
                
                if(edge_index == BS_NULL) {
                    // Path ends here.
                    break;
                } else {
                    // Convert to an actual edge index
                    edge_index -= 2;
                }
                
#ifdef VERBOSE_DEBUG
                cerr << "Taking edge #" << edge_index << " from " << side << endl;
#endif

                // We also should not have negative edges.
                assert(edge_index >= 0);
                
                // Look at the edges we could have taken next
                vector<Edge> edges_out = side % 2 ? edges_on_start(rank_to_id(side / 2)) : edges_on_end(rank_to_id(side / 2));
                
                assert(edge_index < edges_out.size());
                
                Edge& taken = edges_out[edge_index];
                
#ifdef VERBOSE_DEBUG
                cerr << edges_out.size() << " edges possible." << endl;
#endif
                // Follow the edge
                int64_t other_node = taken.from() == rank_to_id(side / 2) ? taken.to() : taken.from();
                bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
                
                // Get the side 
                int64_t other_side = id_to_rank(other_node) * 2 + other_orientation;
                
#ifdef VERBOSE_DEBUG
                cerr << "Go to side " << other_side << endl;
#endif
                
                // Go there with where_to
                offset = where_to(side, offset, other_side);
                side = other_side;
    
            }
            
            found[path_name].push_back(path);
            
        }
    }
    
    return found;
}

XG::destination_t XG::bs_get(int64_t side, int64_t offset) const {
#if GPBWT_MODE == MODE_SDSL
    if(!bs_arrays.empty()) {
        // We still have per-side arrays
        return bs_arrays.at(side - 2)[offset];
    } else {
        // We have a single big array
#ifdef VERBOSE_DEBUG
        cerr << "Range " << side << " has its separator at " << bs_single_array.select(side, BS_SEPARATOR) << endl;
        cerr << "Offset " << offset << " puts us at " << bs_single_array.select(side, BS_SEPARATOR) + 1 + offset << endl;
#endif
        return bs_single_array[bs_single_array.select(side, BS_SEPARATOR) + 1 + offset];
    }
#elif GPBWT_MODE == MODE_DYNAMIC
    // Start after the separator for the side and go offset from there.
    
    auto& bs_single_array = const_cast<rank_select_int_vector&>(this->bs_single_array);
    
    return bs_single_array.at(bs_single_array.select(side - 2, BS_SEPARATOR) + 1 + offset);
#endif
    
}

size_t XG::bs_rank(int64_t side, int64_t offset, destination_t value) const {
#if GPBWT_MODE == MODE_SDSL
    if(!bs_arrays.empty()) {
        throw runtime_error("No rank support until bs_bake() is called!");
    } else {
        size_t range_start = bs_single_array.select(side, BS_SEPARATOR) + 1;
        return bs_single_array.rank(range_start + offset, value) - bs_single_array.rank(range_start, value);
    }
#elif GPBWT_MODE == MODE_DYNAMIC

    auto& bs_single_array = const_cast<rank_select_int_vector&>(this->bs_single_array);

    // Where does the B_s[] range for the side we're interested in start?
    int64_t bs_start = bs_single_array.select(side - 2, BS_SEPARATOR) + 1;
    
    // Get the rank difference between the start and the start plus the offset.
    return bs_single_array.rank(bs_start + offset, value) - bs_single_array.rank(bs_start, value);
#endif
}

void XG::bs_set(int64_t side, vector<destination_t> new_array) {
#if GPBWT_MODE == MODE_SDSL
    // We always know bs_arrays will be big enough.
    
    // Turn the new array into a string of bytes.
    // TODO: none of the destinations can be 255 or greater!
    bs_arrays.at(side - 2) = string(new_array.size(), 0);
    copy(new_array.begin(), new_array.end(), bs_arrays.at(side - 2).begin());
    
#ifdef VERBOSE_DEBUG
    cerr << "B_s for " << side << ": ";
    for(auto entry : bs_arrays.at(side - 2)) { 
        cerr << to_string(entry);
    }
    cerr << endl;
#endif
#elif GPBWT_MODE == MODE_DYNAMIC
    auto current_separators = bs_single_array.rank(bs_single_array.size(), BS_SEPARATOR);
    while(current_separators <= side - 2) {
        // Tack on separators until we have enough
        bs_single_array.insert(bs_single_array.size(), BS_SEPARATOR);
        current_separators++;
    }

    // Where does the block we want start?
    size_t this_range_start = bs_single_array.select(side - 2, BS_SEPARATOR) + 1;
    
    // Where is the first spot not in the range for this side?
    int64_t this_range_past_end = (side - 2 == bs_single_array.rank(bs_single_array.size(), BS_SEPARATOR) - 1 ?
        bs_single_array.size() : bs_single_array.select(side - 2 + 1, BS_SEPARATOR));
    
    if(this_range_start != this_range_past_end) {
        // We can't overwrite! Just explode.
        throw runtime_error("B_s overwrite not supported");
    }
    
    size_t bs_insert_index = this_range_start;
    
    for(auto destination : new_array) {
        // Blit everything into the B_s array
        bs_single_array.insert(bs_insert_index, destination);
        bs_insert_index++;
    }
#endif
}

void XG::bs_insert(int64_t side, int64_t offset, destination_t value) {
#if GPBWT_MODE == MODE_SDSL
    // This is a pretty slow insert. Use set instead.

    auto& array_to_expand = bs_arrays.at(side - 2);
    
    // Stick one copy of the new entry in at the right position.
    array_to_expand.insert(offset, 1, value);
#elif GPBWT_MODE == MODE_DYNAMIC
     // Find the place to put it in the correct side's B_s and insert
     bs_single_array.insert(bs_single_array.select(side - 2, BS_SEPARATOR) + 1 + offset, value);
#endif
}

void XG::bs_bake() {
#if GPBWT_MODE == MODE_SDSL
    // First pass: determine required size
    size_t total_visits = 1;
    for(auto& bs_array : bs_arrays) {
        total_visits += 1; // For the separator
        total_visits += bs_array.size();
    }

#ifdef VERBOSE_DEBUG
    cerr << "Allocating giant B_s array of " << total_visits << " bytes..." << endl;
#endif
    // Move over to a single array which is big enough to start out with.
    string all_bs_arrays(total_visits, 0);
    
    // Where are we writing to?
    size_t pos = 0;
    
    // Start with a separator for sides 0 and 1.
    // We don't start at run 0 because we can't select(0, BS_SEPARATOR).
    all_bs_arrays[pos++] = BS_SEPARATOR;
    
#ifdef VERBOSE_DEBUG
    cerr << "Baking " << bs_arrays.size() << " sides' arrays..." << endl;
#endif
    
    for(auto& bs_array : bs_arrays) {
        // Stick everything together with a separator at the front of every
        // range.
        all_bs_arrays[pos++] = BS_SEPARATOR;
        for(size_t i = 0; i < bs_array.size(); i++) {
            all_bs_arrays[pos++] = bs_array[i];
        }
        bs_array.clear();
    }
    
#ifdef VERBOSE_DEBUG
    cerr << "B_s: ";
    for(auto entry : all_bs_arrays) { 
        cerr << to_string(entry);
    }
    cerr << endl;
#endif
    
    // Rebuild based on the entire concatenated string.
    construct_im(bs_single_array, all_bs_arrays, 1);
    
    bs_arrays.clear();
#endif
}

void XG::tn_bake() {
    names_str.append("$"); // tail marker to avoid edge case
    construct_im(tn_csa, names_str, 1);
    // build the bv
    bit_vector tn_bv(names_str.size());
    int i = 0;
    for (auto c : names_str) {
        if (c == '$') tn_bv[i] = 1;
        ++i;
    }
    names_str.clear();
    
    // make a compressed version and its supports
    util::assign(tn_cbv, sd_vector<>(tn_bv));
    util::assign(tn_cbv_rank, sd_vector<>::rank_1_type(&tn_cbv));
    util::assign(tn_cbv_select, sd_vector<>::select_1_type(&tn_cbv));
}


void XG::bs_dump(ostream& out) const {
#if GPBWT_MODE == MODE_SDSL
    if(!bs_arrays.empty()) {
        // We still have per-side arrays
        
        for(auto& array : bs_arrays) {
            // For each side in order
            out << "---SEP---" << endl;
            for(auto& entry : array) {
                if(entry == BS_NULL) {
                    // Mark nulls
                    out << "**NULL**" << endl;
                } else {
                    // Output adjusted, actual edge numbers
                    out << entry - 2 << endl;
                }
            }
        }
    } else {
        // We have a single big array
        
        for(size_t i = 0; i < bs_single_array.size(); i++) {
            auto entry = bs_single_array[i];
            if(entry == BS_SEPARATOR) {
                /// Mark separators
                out << "---SEP---" << endl;
            } else if(entry == BS_NULL) {
                // Mark nulls
                out << "**NULL**" << endl;
            } else {
                // Output adjusted, actual edge numbers
                out << entry - 2 << endl;
            }
        }
        
        // Now dump the wavelet tree to a string
        stringstream wt_dump;
        bs_single_array.serialize(wt_dump, nullptr, "");
        
        out << endl << "++++Serialized Bits++++" << endl;
        
        // Turn all the bytes into binary representations
        for(auto& letter : wt_dump.str()) {
            // Make each into a bitset
            bitset<8> bits(letter);
            out << bits << endl;
        }
        
    }
#elif GPBWT_MODE == MODE_DYNAMIC
    auto& bs_single_array = const_cast<rank_select_int_vector&>(this->bs_single_array);
    
    for(size_t i = 0; i < bs_single_array.size(); i++) {
        auto entry = bs_single_array.at(i);
        if(entry == BS_SEPARATOR) {
            /// Mark separators
            out << "---SEP---" << endl;
        } else if(entry == BS_NULL) {
            // Mark nulls
            out << "**NULL**" << endl;
        } else {
            // Output adjusted, actual edge numbers
            out << entry - 2 << endl;
        }
    }
#endif
    
}

size_t XG::count_matches(const thread_t& t) const {
    // This is just a really simple wrapper that does a single extend
    ThreadSearchState state;
    extend_search(state, t);
    return state.count();
}

size_t XG::count_matches(const Path& t) const {
    // We assume the path is a thread and convert.
    thread_t thread;
    for(size_t i = 0; i < t.mapping_size(); i++) {
        // Convert the mapping
        ThreadMapping m = {t.mapping(i).position().node_id(), t.mapping(i).position().is_reverse()};
        
        thread.push_back(m);
    }
    
    // Count matches to the converted thread
    return count_matches(thread);
}

void XG::extend_search(ThreadSearchState& state, const thread_t& t) const {
    
#ifdef VERBOSE_DEBUG
    cerr << "Looking for path: ";
    for(int64_t i = 0; i < t.size(); i++) {
        // For each item in the path
        const ThreadMapping& mapping = t[i];
        int64_t next_id = mapping.node_id;
        bool next_is_reverse = mapping.is_reverse;
        int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;
        cerr << next_side << "; ";
    }
    cerr << endl;
#endif
    
    
    for(int64_t i = 0; i < t.size(); i++) {
        // For each item in the path
        const ThreadMapping& mapping = t[i];
        
        if(state.is_empty()) {
            // Don't bother trying to extend empty things.
            
#ifdef VERBOSE_DEBUG
            cerr << "Nothing selected!" << endl;
#endif
            
            break;
        }
        
        // TODO: make this mapping to side thing a function
        int64_t next_id = mapping.node_id;
        bool next_is_reverse = mapping.is_reverse;
        int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;
        
#ifdef VERBOSE_DEBUG
        cerr << "Extend mapping to " << state.current_side << " range " << state.range_start << " to " << state.range_end << " with " << next_side << endl;
#endif
        
        if(state.current_side == 0) {
            // If the state is a start state, just select the whole node using
            // the node usage count in this orientation. TODO: orientation not
            // really important unless we're going to search during a path
            // addition.
            state.range_start = 0;
            state.range_end = h_civ[node_graph_idx(next_id) * 2 + next_is_reverse];
            
#ifdef VERBOSE_DEBUG
            cerr << "\tFound " << state.range_end << " threads present here." << endl;
            
            int64_t here = node_graph_idx(next_id) * 2 + next_is_reverse;
            cerr << here << endl;
            for(int64_t i = here - 5; i < here + 5; i++) {
                if(i >= 0) {
                    cerr << "\t\t" << (i == here ? "*" : " ") << "h_civ[" << i << "] = " << h_civ[i] << endl;
                }
            }
            
#endif
            
        } else {
            // Else, look at where the path goes to and apply the where_to function to shrink the range down.
            state.range_start = where_to(state.current_side, state.range_start, next_side);
            state.range_end = where_to(state.current_side, state.range_end, next_side);
            
#ifdef VERBOSE_DEBUG
            cerr << "\tFound " << state.range_start << " to " << state.range_end << " threads continuing through." << endl;
#endif
            
        }
        
        // Update the side that the state is on
        state.current_side = next_side;
    }
}

void XG::extend_search(ThreadSearchState& state, const ThreadMapping& t) const {
    // Just make it into a temporary vector
    extend_search(state, thread_t{t});
}

int64_t XG::threads_starting_on_side(int64_t side) const {
    return side_thread_wt.rank(side_thread_wt.size(), side);
}

int64_t XG::thread_starting_at(int64_t side, int64_t offset) const {
    return side_thread_wt.select(offset+1, side)+1;
}

pair<int64_t, int64_t> XG::thread_start(int64_t thread_id, bool is_rev) const {
    int64_t idx = thread_id*2 + is_rev;
    return make_pair(tin_civ[idx], tio_civ[idx]);
}

string XG::thread_name(int64_t thread_id) const {
    // convert to forward thread
    if (thread_id > side_thread_wt.size()/2) {
        thread_id -= side_thread_wt.size()/2;
    }
    return extract(tn_csa,
                   tn_cbv_select(thread_id)+1, // skip delimiter
                   tn_cbv_select(thread_id+1)-1); // to next-1
}

vector<int64_t> XG::threads_named_starting(const string& pattern) const {
    vector<int64_t> results;
    // threads named starting with this, so add the sep character so our occs give us thread ids
    auto occs = locate(tn_csa, "$" + pattern);
    // for each occurrance get the thread id
    for (size_t i = 0; i < occs.size(); ++i) {
        results.push_back(tn_cbv_rank(occs[i])+1);
    }
    return results;
}

XG::ThreadSearchState XG::select_starting(const ThreadMapping& start) const {
    // We need to select just the threads starting at this node with this
    // mapping, rather than those coming in from elsewhere.
    
    // Make a search state with nothing searched
    ThreadSearchState state;
    
    // Say we've searched this ThreadMapping's side.    
    state.current_side = id_rev_to_side(start.node_id, start.is_reverse);
    
    // Threads starting at a node come first, so select from 0 to the number of
    // threads that start there.
    state.range_start = 0;
    state.range_end = ts_civ[state.current_side];
    
    return state;
}

int64_t XG::id_rev_to_side(int64_t id, bool is_rev) const {
    return id_to_rank(id) * 2 + is_rev;
}

pair<int64_t, bool> XG::side_to_id_rev(int64_t side) const {
    return make_pair(side / 2, side % 2);
}

XG::ThreadSearchState XG::select_continuing(const ThreadMapping& start) const {
    // We need to select just the threads coming in from elsewhere, and not
    // those starting here.
    
    // Make a search state with nothing searched
    ThreadSearchState state;
    
    // Say we've searched this ThreadMapping's side.
    state.current_side = id_rev_to_side(start.node_id, start.is_reverse);
    
    // Threads starting at a node come first, so select from past them to the
    // number of threads total on the node.
    state.range_start = ts_civ[state.current_side];
    state.range_end = h_civ[state.current_side];
    
    return state;
}


size_t serialize(XG::rank_select_int_vector& to_serialize, ostream& out,
    sdsl::structure_tree_node* parent, const std::string name) {
#if GPBWT_MODE == MODE_SDSL
    // Just delegate to the SDSL code
    return to_serialize.serialize(out, parent, name);

#elif GPBWT_MODE == MODE_DYNAMIC
    // We need to check to make sure we're actually writing the correct numbers of bytes.
    size_t start = out.tellp();
    
    // We just use the DYNAMIC serialization.  
    size_t written = to_serialize.serialize(out);
    
    // TODO: when https://github.com/nicolaprezza/DYNAMIC/issues/4 is closed,
    // trust the sizes that DYNAMIC reports. For now, second-guess it and just
    // look at how far the stream has actually moved.
    written = (size_t) out.tellp() - start;
    
    // And then do the structure tree stuff
    sdsl::structure_tree_node* child = structure_tree::add_child(parent, name, sdsl::util::class_name(to_serialize));
    sdsl::structure_tree::add_size(child, written);
    
    return written;
#endif
}

void deserialize(XG::rank_select_int_vector& target, istream& in) {
#if GPBWT_MODE == MODE_SDSL
    // We just load using the SDSL deserialization code
    target.load(in);
#elif GPBWT_MODE == MODE_DYNAMIC
    // The DYNAMIC code has the same API
    target.load(in);
#endif
}

bool edges_equivalent(const Edge& e1, const Edge& e2) {
    return ((e1.from() == e2.from() && e1.to() == e2.to() && e1.from_start() == e2.from_start() && e1.to_end() == e2.to_end()) ||
        (e1.from() == e2.to() && e1.to() == e2.from() && e1.from_start() == !e2.to_end() && e1.to_end() == !e2.from_start()));
}

bool relative_orientation(const Edge& e1, const Edge& e2) {
    assert(edges_equivalent(e1, e2));
    
    // Just use the reverse-equivalence check from edges_equivalent.
    return e1.from() == e2.to() && e1.to() == e2.from() && e1.from_start() == !e2.to_end() && e1.to_end() == !e2.from_start();
}

bool arrive_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse) {
    if(e.to() == node_id && (node_is_reverse == e.to_end())) {
        // We can follow the edge forwards and arrive at the correct side of the node
        return false;
    } else if(e.to() == e.from() && e.from_start() != e.to_end()) {
        // Reversing self loop
        return false;
    }
    // Otherwise, since we know the edge goes to the node, we have to take it backward.
    return true;
}

bool depart_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse) {
    if(e.from() == node_id && (node_is_reverse == e.from_start())) {
        // We can follow the edge forwards and arrive at the correct side of the node
        return false;
    } else if(e.to() == e.from() && e.from_start() != e.to_end()) {
        // Reversing self loop
        return false;
    }
    // Otherwise, since we know the edge goes to the node, we have to take it backward.
    return true;
}

Edge make_edge(int64_t from, bool from_start, int64_t to, bool to_end) {
    Edge e;
    e.set_from(from);
    e.set_to(to);
    e.set_from_start(from_start);
    e.set_to_end(to_end);
    
    return e;
}

char reverse_complement(const char& c) {
    switch (c) {
        case 'A': return 'T'; break;
        case 'T': return 'A'; break;
        case 'G': return 'C'; break;
        case 'C': return 'G'; break;
        case 'N': return 'N'; break;
        // Handle the GCSA2 start/stop characters.
        case '#': return '$'; break;
        case '$': return '#'; break;
        default: return 'N';
    }
}

string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        switch (c) {
        case 'A': c = 'T'; break;
        case 'T': c = 'A'; break;
        case 'G': c = 'C'; break;
        case 'C': c = 'G'; break;
        case 'N': c = 'N'; break;
        // Handle the GCSA2 start/stop characters.
        case '#': c = '$'; break;
        case '$': c = '#'; break;
        default: break;
        }
    }
    return rc;
}

void extract_pos(const string& pos_str, int64_t& id, bool& is_rev, size_t& off) {
    // format is id:off for forward, and id:-off for reverse
    // find the colon
    auto s = pos_str.find(":");
    assert(s != string::npos);
    id = stol(pos_str.substr(0, s));
    auto r = pos_str.find("-", s);
    if (r == string::npos) {
        is_rev = false;
        off = stoi(pos_str.substr(s+1, pos_str.size()));
    } else {
        is_rev = true;
        off = stoi(pos_str.substr(r+1, pos_str.size()));
    }
}

void extract_pos_substr(const string& pos_str, int64_t& id, bool& is_rev, size_t& off, size_t& len) {
    // format is id:off:len on forward and id:-off:len on reverse
    auto s = pos_str.find(":");
    assert(s != string::npos);
    id = stol(pos_str.substr(0, s));
    auto r = pos_str.find("-", s);
    if (r == string::npos) {
        is_rev = false;
        // find second colon
        auto t = pos_str.find(":", s+1);
        assert(t != string::npos);
        off = stoi(pos_str.substr(s+1, t-s));
        len = stoi(pos_str.substr(t+1, pos_str.size()));
    } else {
        is_rev = true;
        auto t = pos_str.find(":", r+1);
        assert(t != string::npos);
        off = stoi(pos_str.substr(r+1, t-r+1));
        len = stoi(pos_str.substr(t+1, pos_str.size()));
    }
}

}
