#include "vg.hpp"
#include "stream.hpp"

namespace vg {

using namespace std;


// construct from a stream of protobufs
VG::VG(istream& in, bool showp) {

    // set up uninitialized values
    init();
    show_progress = showp;
    // and if we should show progress
    function<void(uint64_t)> handle_count = [this](uint64_t count) {
        create_progress("loading graph", count);
    };

    // the graph is read in chunks, which are attached to this graph
    uint64_t i = 0;
    function<void(Graph&)> lambda = [this, &i](Graph& g) {
        update_progress(++i);
        // We expect these to not overlap in nodes or edges, so complain if they do.
        extend(g, true);
    };

    stream::for_each(in, lambda, handle_count);

    // store paths in graph
    paths.to_graph(graph);

    destroy_progress();

}

// construct from an arbitrary source of Graph protobuf messages
VG::VG(function<bool(Graph&)>& get_next_graph, bool showp) {
    // set up uninitialized values
    init();
    show_progress = showp;

    // We can't show loading progress since we don't know the total number of
    // subgraphs.

    // Try to load the first graph
    Graph subgraph;
    bool got_subgraph = get_next_graph(subgraph);
    while(got_subgraph) {
        // If there is a valid subgraph, add it to ourselves.
        // We expect these to not overlap in nodes or edges, so complain if they do.
        extend(subgraph, true);
        // Try and load the next subgraph, if it exists.
        got_subgraph = get_next_graph(subgraph);
    }

    // store paths in graph
    paths.to_graph(graph);
}

void VG::serialize_to_ostream(ostream& out, int64_t chunk_size) {

    // save the number of the messages to be serialized into the output file
    int64_t count = graph.node_size() / chunk_size + 1;
    create_progress("saving graph", count);
    // partition the graph into a number of chunks (required by format)
    // constructing subgraphs and writing them to the stream
    function<Graph(uint64_t)> lambda =
        [this, chunk_size](uint64_t i) -> Graph {
        VG g;
        for (int64_t j = i * chunk_size;
             j < (i+1)*chunk_size && j < graph.node_size();
             ++j) {
            Node* node = graph.mutable_node(j);
            // Grab the node and only the edges where it has the lower ID.
            // This prevents duplication of edges in the serialized output.
            nonoverlapping_node_context(node, g);
        }
        // store paths
        g.paths.to_graph(g.graph);

        update_progress(i);
        return g.graph;
    };

    stream::write(out, count, lambda);

    destroy_progress();
}

void VG::serialize_to_file(const string& file_name, int64_t chunk_size) {
    ofstream f(file_name);
    serialize_to_ostream(f);
    f.close();
}

VG::~VG(void) {
    destroy_alignable_graph();
}

VG::VG(void) {
    init();
}

void VG::init(void) {
    gssw_aligner = NULL;
    current_id = 1;
    show_progress = false;
    progress_message = "progress";
    progress = NULL;
}

VG::VG(set<Node*>& nodes, set<Edge*>& edges) {
    init();
    add_nodes(nodes);
    add_edges(edges);
    sort();
}

// check for conflict (duplicate nodes and edges) occurs within add_* functions

void VG::add_nodes(set<Node*>& nodes) {
    for (auto node : nodes) {
        add_node(*node);
    }
}

void VG::add_edges(set<Edge*>& edges) {
    for (auto edge : edges) {
        add_edge(*edge);
    }
}

void VG::add_nodes(vector<Node>& nodes) {
    for (auto& node : nodes) {
        add_node(node);
    }
}

void VG::add_edges(vector<Edge>& edges) {
    for (auto& edge : edges) {
        add_edge(edge);
    }
}

void VG::add_node(Node& node) {
    if (!has_node(node)) {
        Node* new_node = graph.add_node(); // add it to the graph
        *new_node = node; // overwrite it with the value of the given node
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.node_size()-1;
        //cerr << "Added node " << new_node->id() << endl;
    }
}

void VG::add_edge(Edge& edge) {
    if (!has_edge(edge)) {
        Edge* new_edge = graph.add_edge(); // add it to the graph
        *new_edge = edge;
        set_edge(new_edge);
        edge_index[new_edge] = graph.edge_size()-1;
    }
}

int64_t VG::node_count(void) {
    return graph.node_size();
}

int64_t VG::edge_count(void) {
    return graph.edge_size();
}

vector<pair<int64_t, bool>>& VG::edges_start(Node* node) {
    if(node == nullptr) {
        return empty_edge_ends;
    }
    return edges_start(node->id());
}

vector<pair<int64_t, bool>>& VG::edges_start(int64_t id) {
    if(edges_on_start.count(id) == 0) {
        return empty_edge_ends;
    }
    return edges_on_start[id];
}

vector<pair<int64_t, bool>>& VG::edges_end(Node* node) {
    if(node == nullptr) {
        return empty_edge_ends;
    }
    return edges_end(node->id());
}

vector<pair<int64_t, bool>>& VG::edges_end(int64_t id) {
    if(edges_on_end.count(id) == 0) {
        return empty_edge_ends;
    }
    return edges_on_end[id];
}

int VG::start_degree(Node* node) {
    return edges_start(node).size();
}

int VG::end_degree(Node* node) {
    return edges_end(node).size();
}

int VG::left_degree(NodeTraversal node) {
    // If we're backward, the end is on the left. Otherwise, the start is.
    return node.backward ? end_degree(node.node) : start_degree(node.node);
}

int VG::right_degree(NodeTraversal node) {
    // If we're backward, the start is on the right. Otherwise, the end is.
    return node.backward ? start_degree(node.node) : end_degree(node.node);
}

void VG::edges_of_node(Node* node, vector<Edge*>& edges) {
    for(pair<int64_t, bool>& off_start : edges_start(node)) {
        // Go through the edges on this node's start
        Edge* edge = edge_by_sides[NodeSide::pair_from_start_edge(node->id(), off_start)];
        if (!edge) {
            cerr << "error:[VG::edges_of_node] nonexistent edge " << off_start.first << " start <-> "
                 << node->id() << (off_start.second ? " start" : " end") << endl;
            exit(1);
        }
        edges.push_back(edge);
    }

    for(pair<int64_t, bool>& off_end : edges_end(node)) {
        // And on its end
        Edge* edge = edge_by_sides[NodeSide::pair_from_end_edge(node->id(), off_end)];
        if (!edge) {
            cerr << "error:[VG::edges_of_node] nonexistent edge " << off_end.first << " end <-> "
                 << node->id() << (off_end.second ? " end" : " start") << endl;
            exit(1);
        }
        edges.push_back(edge);
    }
}

void VG::edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges) {
    for (set<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        vector<Edge*> ev;
        edges_of_node(*n, ev);
        for (vector<Edge*>::iterator e = ev.begin(); e != ev.end(); ++e) {
            edges.insert(*e);
        }
    }
}

int64_t VG::total_length_of_nodes(void) {
    int64_t length = 0;
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        length += n->sequence().size();
    }
    return length;
}

void VG::build_indexes(void) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        edge_index[e] = i;
        set_edge(e);
    }
}

void VG::clear_indexes(void) {
    node_index.clear();
    node_by_id.clear();
    edge_by_sides.clear();
    edge_index.clear();
    edges_on_start.clear();
    edges_on_end.clear();
}

void VG::clear_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    node_index.clear_no_resize();
    node_by_id.clear_no_resize();
    edge_by_sides.clear_no_resize();
    edge_index.clear_no_resize();
    edges_on_start.clear_no_resize();
    edges_on_end.clear_no_resize();
#else
    clear_indexes();
#endif
}

void VG::resize_indexes(void) {
    node_index.resize(graph.node_size());
    node_by_id.resize(graph.node_size());
    edge_by_sides.resize(graph.edge_size());
    edge_index.resize(graph.edge_size());
    edges_on_start.resize(graph.edge_size());
    edges_on_end.resize(graph.edge_size());
}

void VG::rebuild_indexes(void) {
    //clear_indexes();
    //resize_indexes();
    clear_indexes_no_resize();
    build_indexes();
    paths.rebuild_node_mapping();
}

bool VG::empty(void) {
    return graph.node_size() == graph.edge_size() == 0;
}

bool VG::has_node(Node* node) {
    return node && has_node(node->id());
}

bool VG::has_node(Node& node) {
    return has_node(node.id());
}

bool VG::has_node(int64_t id) {
    return node_by_id.find(id) != node_by_id.end();
}

bool VG::has_edge(Edge* edge) {
    return edge && has_edge(*edge);
}

bool VG::has_edge(Edge& edge) {
    return edge_by_sides.find(NodeSide::pair_from_edge(edge)) != edge_by_sides.end();
}

bool VG::has_edge(const NodeSide& side1, const NodeSide& side2) {
    return edge_by_sides.find(minmax(side1, side2)) != edge_by_sides.end();
}

bool VG::has_edge(const pair<NodeSide, NodeSide>& sides) {
    return has_edge(sides.first, sides.second);
}

// remove duplicated nodes and edges that would occur if we merged the graphs
void VG::remove_duplicated_in(VG& g) {
    vector<Node*> nodes_to_destroy;
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (g.has_node(n)) {
            nodes_to_destroy.push_back(n);
        }
    }
    vector<Edge*> edges_to_destroy;
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (g.has_edge(e)) {
            edges_to_destroy.push_back(e);
        }
    }
    for (vector<Node*>::iterator n = nodes_to_destroy.begin();
         n != nodes_to_destroy.end(); ++n) {
        g.destroy_node(g.get_node((*n)->id()));
    }
    for (vector<Edge*>::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        // Find and destroy the edge that does the same thing in g.
        g.destroy_edge(g.get_edge(NodeSide::pair_from_edge(*e)));
    }
}

void VG::merge_union(VG& g) {
    // remove duplicates, then merge
    remove_duplicated_in(g);
    if (g.graph.node_size() > 0) {
        merge(g.graph);
    }
}

void VG::merge(VG& g) {
    merge(g.graph);
}

// this merges without any validity checks
// this could be rather expensive if the graphs to merge are largely overlapping
void VG::merge(Graph& g) {
    graph.mutable_node()->MergeFrom(g.node());
    graph.mutable_edge()->MergeFrom(g.edge());
    rebuild_indexes();
}

// iterates over nodes and edges, adding them in when they don't already exist
void VG::extend(VG& g, bool warn_on_duplicates) {
    for (int64_t i = 0; i < g.graph.node_size(); ++i) {
        Node* n = g.graph.mutable_node(i);
        if(n->id() == 0) {
            cerr << "[vg] warning: node ID 0 is not allowed. Skipping." << endl;
        } else if (!has_node(n)) {
            add_node(*n);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: node ID " << n->id() << " appears multiple times. Skipping." << endl;
        }
    }
    for (int64_t i = 0; i < g.graph.edge_size(); ++i) {
        Edge* e = g.graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: edge " << e->from() << (e->from_start() ? " start" : " end") << " <-> "
                 << e->to() << (e->to_end() ? " end" : " start") << " appears multiple times. Skipping." << endl;
        }
    }
    paths.append(g.paths);
}

// TODO: unify with above. The only difference is what's done with the paths.
void VG::extend(Graph& graph, bool warn_on_duplicates) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if(n->id() == 0) {
            cerr << "[vg] warning: node ID 0 is not allowed. Skipping." << endl;
        } else if (!has_node(n)) {
            add_node(*n);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: node ID " << n->id() << " appears multiple times. Skipping." << endl;
        }
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        } else if(warn_on_duplicates) {
            cerr << "[vg] warning: edge " << e->from() << (e->from_start() ? " start" : " end") << " <-> "
                 << e->to() << (e->to_end() ? " end" : " start") << " appears multiple times. Skipping." << endl;
        }
    }
    paths.append(graph);
}

// extend this graph by g, connecting the tails of this graph to the heads of the other
// the ids of the second graph are modified for compact representation
void VG::append(VG& g) {

    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();

    // assume we've already compacted the other, or that id compaction doesn't matter
    // just get out of the way
    g.increment_node_ids(max_node_id());

    // get the heads of the other graph, now that we've compacted the ids
    vector<Node*> heads = g.head_nodes();
    // The heads are guaranteed to be forward-oriented.
    vector<int64_t> heads_ids;
    for (Node* n : heads) {
        heads_ids.push_back(n->id());
    }

    // get the current tails of this graph
    vector<Node*> tails = tail_nodes();
    // The tails are also guaranteed to be forward-oriented.
    vector<int64_t> tails_ids;
    for (Node* n : tails) {
        tails_ids.push_back(n->id());
    }

    // add in the other graph
    // note that we don't use merge_union because we are ensured non-overlapping ids
    merge(g);

    /*
    cerr << "this graph size " << node_count() << " nodes " << edge_count() << " edges" << endl;
    cerr << "in append with " << heads.size() << " heads and " << tails.size() << " tails" << endl;
    */

    // now join the tails to heads
    for (int64_t& tail : tails_ids) {
        for (int64_t& head : heads_ids) {
            // Connect the tail to the head with a left to right edge.
            create_edge(tail, head);
        }
    }

    // and join paths that are embedded in the graph, where path names are the same
    paths.append(g.paths);
}

void VG::combine(VG& g) {
    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();
    g.increment_node_ids(max_node_id());
    // now add it into the current graph, without connecting any nodes
    extend(g);
}

int64_t VG::max_node_id(void) {
    int64_t max_id = 0;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (n->id() > max_id) {
            max_id = n->id();
        }
    }
    return max_id;
}

int64_t VG::min_node_id(void) {
    int64_t min_id = max_node_id();
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (n->id() < min_id) {
            min_id = n->id();
        }
    }
    return min_id;
}

void VG::compact_ids(void) {
    hash_map<int64_t, int64_t> new_id;
    int64_t id = 1; // start at 1
    for_each_node([&id, &new_id](Node* n) {
            new_id[n->id()] = id++; });
//#pragma omp parallel for
    for_each_node([&new_id](Node* n) {
            n->set_id(new_id[n->id()]); });
//#pragma omp parallel for
    for_each_edge([&new_id](Edge* e) {
            e->set_from(new_id[e->from()]);
            e->set_to(new_id[e->to()]); });
    paths.for_each_mapping([&new_id](Mapping* m) {
            m->mutable_position()->set_node_id(new_id[m->position().node_id()]);
        });
    rebuild_indexes();
}

void VG::increment_node_ids(int64_t increment) {
    for_each_node_parallel([increment](Node* n) {
            n->set_id(n->id()+increment);
        });
    for_each_edge_parallel([increment](Edge* e) {
            e->set_from(e->from()+increment);
            e->set_to(e->to()+increment);
        });
    rebuild_indexes();
    paths.increment_node_ids(increment);
}

void VG::decrement_node_ids(int64_t decrement) {
    increment_node_ids(-decrement);
}

void VG::swap_node_id(int64_t node_id, int64_t new_id) {
    swap_node_id(node_by_id[node_id], new_id);
}

void VG::swap_node_id(Node* node, int64_t new_id) {

    int edge_n = edge_count();
    int64_t old_id = node->id();
    node->set_id(new_id);
    node_by_id.erase(old_id);

    // we check if the old node exists, and bail out if we're not doing what we expect
    assert(node_by_id.find(new_id) == node_by_id.end());

    // otherwise move to a new id
    node_by_id[new_id] = node;

    // These are sets, so if we try to destroy and recreate the same edge from
    // both ends (i.e. if they both go to this node) we will only do it once.
    set<pair<NodeSide, NodeSide>> edges_to_destroy;
    set<pair<NodeSide, NodeSide>> edges_to_create;

    // Define a function that we will run on every edge this node is involved in
    auto fix_edge = [&](Edge* edge) {
    
        // Destroy that edge
        edges_to_destroy.emplace(NodeSide(edge->from(), !edge->from_start()), NodeSide(edge->to(), edge->to_end()));

        // Make a new edge with our new ID as from or to (or both), depending on which it was before.
        // TODO: Is there a cleaner way to do this?
        if(edge->from() == old_id) {
            if(edge->to() == old_id) {
                edges_to_create.emplace(NodeSide(new_id, !edge->from_start()), NodeSide(new_id, edge->to_end()));
            } else {
                edges_to_create.emplace(NodeSide(new_id, !edge->from_start()), NodeSide(edge->to(), edge->to_end()));
            }
        } else {
            edges_to_create.emplace(NodeSide(edge->from(), !edge->from_start()), NodeSide(new_id, edge->to_end()));
        }
    
    };

    for(pair<int64_t, bool>& other : edges_start(old_id)) {
        // Get the actual Edge
        // We're at a start, so we go to the end of the other node normally, and the start if the other node is backward
        Edge* edge = edge_by_sides[minmax(NodeSide(old_id, false), NodeSide(other.first, !other.second))];
        
        // Plan to fix up its IDs.
        fix_edge(edge);
    }

    for(pair<int64_t, bool>& other : edges_end(old_id)) {
        // Get the actual Edge
        // We're at an end, so we go to the start of the other node normally, and the end if the other node is backward
        Edge* edge = edge_by_sides[minmax(NodeSide(old_id, true), NodeSide(other.first, other.second))];

        // Plan to fix up its IDs.
        fix_edge(edge);
    }

    assert(edges_to_destroy.size() == edges_to_create.size());

    for (auto& e : edges_to_destroy) {
        // Destroy the edge (only one can exist between any two nodes)
        destroy_edge(e.first, e.second);
    }

    for (auto& e : edges_to_create) {
        // Make an edge with the appropriate start and end flags
        create_edge(e.first, e.second);
    }

    assert(edge_n == edge_count());

    // we maintain a valid graph
    // this an expensive check but should work (for testing only)
    //assert(is_valid());

}



// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//

void VG::vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, set<vcflib::VariantAllele> >& altp,
                                int start_pos,
                                int stop_pos,
                                int max_node_size) {

    create_progress("parsing variants", records.size());

    for (int i = 0; i < records.size(); ++i) {
        vcflib::Variant& var = records.at(i);
        // decompose to alts
        bool flat_input_vcf = false; // hack
        map<string, vector<vcflib::VariantAllele> > alternates
            = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());
        for (auto& alleles : alternates) {
            for (auto& allele : alleles.second) {
                 altp[allele.position].insert(allele);
                 if (i % 10000 == 0) {
                     update_progress(altp.size());
                 }
             }
         }
     }
     destroy_progress();
 }

 void VG::slice_alleles(map<long, set<vcflib::VariantAllele> >& altp,
                        int start_pos,
                        int stop_pos,
                        int max_node_size) {

     auto enforce_node_size_limit =
         [this, max_node_size, &altp]
         (int curr_pos, int& last_pos) {
         int last_ref_size = curr_pos - last_pos;
         update_progress(last_pos);
         if (max_node_size && last_ref_size > max_node_size) {
             int div = 2;
             while (last_ref_size/div > max_node_size) {
                 ++div;
             }
             int segment_size = last_ref_size/div;
             int i = 0;
             while (last_pos + i < curr_pos) {
                 altp[last_pos+i];  // empty cut
                 i += segment_size;
                 update_progress(last_pos + i);
             }
         }
     };

     if (max_node_size > 0) {
         create_progress("enforcing node size limit ", (altp.empty()? 0 : altp.rbegin()->first));
         // break apart big nodes
         int last_pos = start_pos;
         for (auto& position : altp) {
             auto& alleles = position.second;
             enforce_node_size_limit(position.first, last_pos);
             for (auto& allele : alleles) {
                 // cut the last reference sequence into bite-sized pieces
                 last_pos = max(position.first + allele.ref.size(), (long unsigned int) last_pos);
            }
        }
        enforce_node_size_limit(stop_pos, last_pos);
        destroy_progress();
    }

}

void VG::dice_nodes(int max_node_size) {
    if (max_node_size) {
        vector<Node*> nodes; nodes.reserve(size());
        for_each_node(
            [this, &nodes](Node* n) {
                nodes.push_back(n);
            });
        auto lambda =
            [this, max_node_size](Node* n) {
            int node_size = n->sequence().size();
            if (node_size > max_node_size) {
                Node* l = NULL;
                Node* r = NULL;
                int div = 2;
                while (node_size/div > max_node_size) {
                    ++div;
                }
                int segment_size = node_size/div;
                int i = 0;
                while (i < node_size) {
                    divide_node(n, i, l, r);
                    n = r;
                    i += segment_size;
                }
            }
        };
        for (int i = 0; i < nodes.size(); ++i) {
            lambda(nodes[i]);
        }
    }
}


void VG::from_alleles(const map<long, set<vcflib::VariantAllele> >& altp,
                      string& seq,
                      string& name) {

    //init();
    this->name = name;

    int tid = omp_get_thread_num();
#ifdef debug
#pragma omp critical (cerr)
    {
        cerr << tid << ": in from_alleles" << endl;
        cerr << tid << ": with " << altp.size() << " vars" << endl;
        cerr << tid << ": and " << seq.size() << "bp" << endl;
        cerr << seq << endl;
    }
#endif


    // maintains the path of the seq in the graph
    map<long, int64_t> seq_node_ids;
    // track the last nodes so that we can connect everything
    // completely when variants occur in succession
    map<long, set<Node*> > nodes_by_end_position;
    map<long, set<Node*> > nodes_by_start_position;


    Node* seq_node = create_node(seq);
    seq_node_ids[0] = seq_node->id();

    for (auto& va : altp) {

        const set<vcflib::VariantAllele>& alleles = va.second;

        // if alleles are empty, we just cut at this point
        if (alleles.empty()) {
            Node* l = NULL; Node* r = NULL;
            divide_path(seq_node_ids, va.first, l, r);
        }

        for (auto allele : alleles) {

            // skip ref-matching alleles; these are not informative
            if (allele.ref == allele.alt) {
                continue;
            }

#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << allele << endl;
#endif

            // 0/1 based conversion happens in offset
            long allele_start_pos = allele.position;
            long allele_end_pos = allele_start_pos + allele.ref.size();
            // for ordering, set insertion start position at +1
            // otherwise insertions at the same position will loop infinitely
            //if (allele_start_pos == allele_end_pos) allele_end_pos++;

            if (allele_start_pos == 0) {
                // ensures that we can handle variation at first position
                // (important when aligning)
                Node* root = create_node("");
                seq_node_ids[-1] = root->id();
                nodes_by_start_position[-1].insert(root);
                nodes_by_end_position[0].insert(root);
            }

            Node* left_seq_node = NULL;
            Node* middle_seq_node = NULL;
            Node* right_seq_node = NULL;

            // make one cut at the ref-path relative start of the allele
            divide_path(seq_node_ids,
                        allele_start_pos,
                        left_seq_node,
                        right_seq_node);

            // if the ref portion of the allele is not empty, then we need to make another cut
            if (!allele.ref.empty()) {
                divide_path(seq_node_ids,
                            allele_end_pos,
                            middle_seq_node,
                            right_seq_node);
            }

            Node* alt_node = NULL;
            // create a new alt node and connect the pieces from before
            if (!allele.alt.empty() && !allele.ref.empty()) {
                //cerr << "both alt and ref have sequence" << endl;

                alt_node = create_node(allele.alt);
                create_edge(left_seq_node, alt_node);
                create_edge(alt_node, right_seq_node);

                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(middle_seq_node);
                //nodes_by_end_position[allele_start_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);
                nodes_by_start_position[allele_start_pos].insert(middle_seq_node);

            } else if (!allele.alt.empty()) { // insertion

                alt_node = create_node(allele.alt);
                create_edge(left_seq_node, alt_node);
                create_edge(alt_node, right_seq_node);
                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);

            } else {// otherwise, we have a deletion

                create_edge(left_seq_node, right_seq_node);
                nodes_by_end_position[allele_end_pos].insert(left_seq_node);
                nodes_by_start_position[allele_start_pos].insert(left_seq_node);

            }

#ifdef debug
#pragma omp critical (cerr)
            {
                if (left_seq_node) cerr << tid << ": left_ref " << left_seq_node->id()
                                        << " " << left_seq_node->sequence() << endl;
                if (middle_seq_node) cerr << tid << ": middle_ref " << middle_seq_node->id()
                                          << " " << middle_seq_node->sequence() << endl;
                if (alt_node) cerr << tid << ": alt_node " << alt_node->id()
                                   << " " << alt_node->sequence() << endl;
                if (right_seq_node) cerr << tid << ": right_ref " << right_seq_node->id()
                                         << " " << right_seq_node->sequence() << endl;
            }
#endif

            if (allele_end_pos == seq.size()) {
                // ensures that we can handle variation at last position (important when aligning)
                Node* end = create_node("");
                seq_node_ids[allele_end_pos] = end->id();
                // for consistency, this should be handled below in the start/end connections
                if (alt_node) {
                    create_edge(alt_node, end);
                }
                if (middle_seq_node) {
                    create_edge(middle_seq_node, end);
                }
            }

            //print_edges();
            /*
            if (!is_valid()) {
                cerr << "graph is invalid after variant " << *a << endl;
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
            */

        }

        map<long, set<Node*> >::iterator ep
            = nodes_by_end_position.find(va.first);
        map<long, set<Node*> >::iterator sp
            = nodes_by_start_position.find(va.first);
        if (ep != nodes_by_end_position.end()
            && sp != nodes_by_start_position.end()) {
            set<Node*>& previous_nodes = ep->second;
            set<Node*>& current_nodes = sp->second;
            for (set<Node*>::iterator n = previous_nodes.begin();
                 n != previous_nodes.end(); ++n) {
                for (set<Node*>::iterator m = current_nodes.begin();
                     m != current_nodes.end(); ++m) {
                    if (node_index.find(*n) != node_index.end()
                        && node_index.find(*m) != node_index.end()
                        && !(previous_nodes.count(*n) && current_nodes.count(*n)
                             && previous_nodes.count(*m) && current_nodes.count(*m))
                        ) {
                        /*
                        cerr << "connecting previous "
                             << (*n)->id() << " @end=" << ep->first << " to current "
                             << (*m)->id() << " @start=" << sp->first << endl;
                        */
                        create_edge(*n, *m);
                    }
                }
            }
        }

        // clean up previous
        while (!nodes_by_end_position.empty() && nodes_by_end_position.begin()->first < va.first) {
            nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
        }

        while (!nodes_by_start_position.empty() && nodes_by_start_position.begin()->first < va.first) {
            nodes_by_start_position.erase(nodes_by_start_position.begin()->first);
        }

    }

    // serialize path
    for (auto& p : seq_node_ids) {
        paths.append_mapping(name, p.second);
    }

    sort();
    compact_ids();

}

void VG::from_gfa(istream& in, bool showp) {
    // c++... split...
    // for line in stdin
    string line;
    auto too_many_fields = [&line]() {
        cerr << "[vg] error: too many fields in line " << endl << line << endl;
        exit(1);
    };

    int64_t id1, id2;
    string seq;
    char side1, side2;
    string cigar;
    string path_name;
    bool is_reverse = false;
    while(std::getline(in, line)) {
        stringstream ss(line);
        string item;
        int field = 0;
        char type = '\0';
        while(std::getline(ss, item, '\t')) {
            switch (field++) {
            case 0:
                type = item[0];
                switch (type) {
                case 'L': break;
                case 'S': break;
                case 'H': break;
                case 'P': break;
                default:
                    cerr << "[vg] error: unrecognized field type " << type << endl;
                    exit(1);
                    break;
                }
                break;
            case 1: id1 = atol(item.c_str()); break;
            case 2: {
                switch (type) {
                case 'S': seq = item; break;
                case 'L': side1 = item[0]; break;
                case 'P': path_name = item; break;
                default: break;
                }
            } break;
            case 3:
                switch (type) {
                case 'L': id2 = atol(item.c_str()); break;
                case 'S': too_many_fields(); break;
                case 'P': is_reverse = (item == "+" ? false : true); break;
                default: break;
                }
                break;
            case 4:
                switch (type) {
                case 'L': side2 = item[0]; break;
                case 'S': too_many_fields(); break;
                case 'P': cigar = item; break;
                default: break;
                }
                break;
            case 5:
                switch (type) {
                case 'L': cigar = item; break;
                case 'S': too_many_fields(); break;
                default: break;
                }
                break;
            default:
                too_many_fields();
                break;
            }
        }

        // now that we've parsed, add to the graph
        if (type == 'S') {
            Node node;
            node.set_sequence(seq);
            node.set_id(id1);
            add_node(node);
        } else if (type == 'L') {
            Edge edge;
            edge.set_from(id1);
            edge.set_to(id2);
            if (side1 == '-') edge.set_from_start(true);
            if (side2 == '-') edge.set_to_end(true);
            add_edge(edge);
        } else if (type == 'P') {
            paths.append_mapping(path_name, id1, is_reverse);
        }
    }
}

void VG::print_edges(void) {
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        int64_t f = e->from();
        int64_t t = e->to();
        cerr << f << "->" << t << " ";
    }
    cerr << endl;
}

void VG::create_progress(const string& message, long count) {
    if (show_progress) {
        progress_message = message;
        create_progress(count);
    }
}

void VG::create_progress(long count) {
    if (show_progress) {
        progress_message.resize(30, ' ');
        progress_count = count;
        last_progress = 0;
        progress = new ProgressBar(progress_count, progress_message.c_str());
        progress->Progressed(0);
    }
}

void VG::update_progress(long i) {
    if (show_progress && progress) {
        if (i <= progress_count
            && (long double) (i - last_progress) / (long double) progress_count >= 0.001
            || i == progress_count) {
#pragma omp critical (progress)
            {
                progress->Progressed(i);
                last_progress = i;
            }
        }
    }
}

void VG::destroy_progress(void) {
    if (show_progress && progress) {
        update_progress(progress_count);
        cerr << endl;
        progress_message = "";
        progress_count = 0;
        delete progress;
        progress = NULL;
    }
}

VG::VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target_region,
       bool target_is_chrom,
       int vars_per_region,
       int max_node_size,
       bool showprog) {

    init();

    omp_set_dynamic(1); // use dynamic scheduling

    show_progress = showprog;

    map<string, VG*> refseq_graph;

    vector<string> targets;
    if (!target_region.empty()) {
        targets.push_back(target_region);
    } else {
        for (vector<string>::iterator r = reference.index->sequenceNames.begin();
             r != reference.index->sequenceNames.end(); ++r) {
            targets.push_back(*r);
        }
    }

    // to scale up, we have to avoid big string memcpys
    // this could be accomplished by some deep surgery on the construction routines
    // however, that could be a silly thing to do,
    // because why break something that's conceptually clear
    // and anyway, we want to break the works into chunks
    //
    // there is a development that could be important
    // our chunk size isn't going to reach into the range where we'll have issues (>several megs)
    // so we'll run this for regions of moderate size, scaling up in the case that we run into a big deletion
    //
    //

    for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {

        //string& seq_name = *t;
        string seq_name;
        string target = *t;
        int start_pos = 0, stop_pos = 0;
        // nasty hack for handling single regions
        if (!target_is_chrom) {
            parse_region(target,
                         seq_name,
                         start_pos, 
                         stop_pos);
            if (stop_pos > 0) {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name, start_pos, stop_pos);
                }
            } else {
                if (variantCallFile.is_open()) {
                    variantCallFile.setRegion(seq_name);
                }
                stop_pos = reference.sequenceLength(seq_name);
            }
        } else {
            // the user said the target is just a sequence name
            // and is unsafe to parse as it may contain ':' or '-'
            // for example "gi|568815592:29791752-29792749"
            if (variantCallFile.is_open()) {
                variantCallFile.setRegion(target);
            }
            stop_pos = reference.sequenceLength(target);
            seq_name = target;
        }
        vcflib::Variant var(variantCallFile);

        vector<vcflib::Variant>* region = NULL;

        // convert from 1-based input to 0-based internal format
        // and handle the case where we are already doing the whole chromosome
        int64_t start = start_pos ? start_pos - 1 : 0;
        int64_t end = start;

        create_progress("loading variants for " + target, stop_pos-start_pos);
        // get records
        vector<vcflib::Variant> records;
        int i = 0;
        while (variantCallFile.is_open() && variantCallFile.getNextVariant(var)) {
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            // only work with DNA sequences
            if (isDNA) {
                var.position -= 1; // convert to 0-based
                records.push_back(var);
            }
            if (++i % 1000 == 0) update_progress(var.position-start_pos);
        }
        destroy_progress();

        map<long,set<vcflib::VariantAllele> > alleles;
        // decompose records int alleles with offsets against our target sequence
        vcf_records_to_alleles(records, alleles, start_pos, stop_pos, max_node_size);
        records.clear(); // clean up

        // enforce a maximum node size
        // by dividing nodes that are > than the max into the smallest number of
        // even pieces that would be smaller than the max
        slice_alleles(alleles, start_pos, stop_pos, max_node_size);

        // store our construction plans
        deque<Plan*> construction;
        // so we can check which graphs we can safely append
        set<VG*> graph_completed;
        // we add and remove from graph_completed, so track count for logging
        int graphs_completed = 0;
        int final_completed = -1; // hm
        // the construction queue
        list<VG*> graphq;
        int graphq_size = 0; // for efficiency
        // ^^^^ (we need to insert/remove things in the middle of the list,
        // but we also need to be able to quickly determine its size)
        // for tracking progress through the chromosome
        map<VG*, unsigned long> graph_end;

        create_progress("planning construction", stop_pos-start_pos);
        // break into chunks
        int chunk_start = start;
        bool invariant_graph = alleles.empty();
        while (invariant_graph || !alleles.empty()) {
            invariant_graph = false;
            auto* new_alleles = new map<long, set<vcflib::VariantAllele> >;
            // our start position is the "offset" we should subtract from the alleles
            // for correct construction
            //chunk_start = (!chunk_start ? 0 : alleles.begin()->first);
            int chunk_end = chunk_start;
            bool clean_end = true;
            for (int i = 0; (i < vars_per_region || !clean_end) && !alleles.empty(); ++i) {
                auto pos = alleles.begin()->first - chunk_start;
                chunk_end = max(chunk_end, (int)alleles.begin()->first);
                auto& pos_alleles = alleles.begin()->second;
                // apply offset when adding to the new alleles
                auto& curr_pos = (*new_alleles)[pos];
                for (auto& allele : pos_alleles) {
                    auto new_allele = allele;
                    int ref_end = new_allele.ref.size() + new_allele.position;
                    // look through the alleles to see if there is a longer chunk
                    if (ref_end > chunk_end) {
                        chunk_end = ref_end;
                    }
                    new_allele.position = pos;
                    curr_pos.insert(new_allele);
                }
                alleles.erase(alleles.begin());
                // TODO here we need to see if we are neighboring another variant
                // and if we are, keep constructing
                if (alleles.begin()->first <= chunk_end) {
                    clean_end = false;
                } else {
                    clean_end = true;
                }
            }
            // record end position, use target end in the case that we are at the end
            if (alleles.empty()) chunk_end = stop_pos;

            // we set the head graph to be this one, so we aren't obligated to copy the result into this object
            // make a construction plan
            Plan* plan = new Plan(graphq.empty() && targets.size() == 1 ? this : new VG,
                                  new_alleles,
                                  reference.getSubSequence(seq_name,
                                                           chunk_start,
                                                           chunk_end - chunk_start),
                                  seq_name);
            chunk_start = chunk_end;
#pragma omp critical (graphq)
            {
                graphq.push_back(plan->graph);
                construction.push_back(plan);
                if (show_progress) graph_end[plan->graph] = chunk_end;
                update_progress(chunk_end);
            }
        }
#ifdef debug
        cerr << omp_get_thread_num() << ": graphq size " << graphq.size() << endl;
#endif
        graphq_size = graphq.size();
        destroy_progress();

        // this system is not entirely general
        // there will be a problem when the regions of overlapping deletions become too large
        // then the inter-dependence of each region will make parallel construction in this way difficult
        // because the chunks will get too large

        // use this function to merge graphs both during and after the construction iteration
        auto merge_first_two_completed_graphs =
            [this, start_pos, &graph_completed, &graphq, &graphq_size, &graph_end, &final_completed](void) {
            // find the first two consecutive graphs which are completed
            VG* first = NULL;
            VG* second = NULL;
//#pragma omp critical (cerr)
//            cerr << omp_get_thread_num() << ": merging" << endl;
#pragma omp critical (graphq)
            {
                auto itp = graphq.begin(); // previous
                auto itn = itp; if (itp != graphq.end()) ++itn; // next
                // scan the graphq to find consecutive entries that are both completed
                while (itp != itn // there is > 1 entry
                       && itn != graphq.end() // we aren't yet at the end
                       && !(graph_completed.count(*itp) // the two we're looking at aren't completed
                            && graph_completed.count(*itn))) {
                    ++itp; ++itn;
                }

                if (itn != graphq.end()) {
                    // we have two consecutive graphs to merge!
                    first = *itp;
                    second = *itn;
                    // unset graph completed for both
                    graph_completed.erase(first);
                    graph_completed.erase(second);
                    graphq.erase(itn);
                    --graphq_size;
                }
            }

            if (first && second) {
                // combine graphs
                first->append(*second);
#pragma omp critical (graphq)
                {
                    if (final_completed != -1) update_progress(final_completed++);
                    graph_completed.insert(first);
                    graph_end.erase(second);
                }
                delete second;
            }
        };

        create_progress("constructing graph", construction.size());

        // (in parallel) construct each component of the graph
#pragma omp parallel for
        for (int i = 0; i < construction.size(); ++i) {

            int tid = omp_get_thread_num();
            Plan* plan = construction.at(i);
#ifdef debug
#pragma omp critical (cerr)
            cerr << tid << ": " << "constructing graph " << plan->graph << " over "
                 << plan->alleles->size() << " variants in " <<plan->seq.size() << "bp "
                 << plan->name << endl;
#endif

            plan->graph->from_alleles(*plan->alleles,
                                      plan->seq,
                                      plan->name);
#pragma omp critical (graphq)
            {
                update_progress(++graphs_completed);
                graph_completed.insert(plan->graph);
#ifdef debug
#pragma omp critical (cerr)
                cerr << tid << ": " << "constructed graph " << plan->graph << endl;
#endif
            }
            // clean up
            delete plan;

            // concatenate chunks of the result graph together
            merge_first_two_completed_graphs();

        }
        destroy_progress();

        // merge remaining graphs
        final_completed = 0;
        create_progress("merging remaining graphs", graphq.size());
#pragma omp parallel
        {
            bool more_to_merge = true;
            while (more_to_merge) {
                merge_first_two_completed_graphs();
                usleep(10);
#pragma omp critical (graphq)
                more_to_merge = graphq_size > 1;
            }
        }
        destroy_progress();

        // parallel end
        // finalize target

        // our target graph should be the only entry in the graphq
        assert(graphq.size() == 1);
        VG* target_graph = graphq.front();

        // store it in our results
        refseq_graph[target] = target_graph;

        create_progress("joining graphs", target_graph->size());
        // clean up "null" nodes that are used for maintaining structure between temporary subgraphs
        target_graph->remove_null_nodes_forwarding_edges();
        destroy_progress();

        // then use topological sorting and re-compression of the id space to make sure that
        create_progress("topologically sorting", target_graph->size());
        target_graph->sort();
        destroy_progress();

        create_progress("compacting ids", target_graph->size());
        // we get identical graphs no matter what the region size is
        target_graph->compact_ids();
        destroy_progress();

    }

    // hack for efficiency when constructing over a single chromosome
    if (refseq_graph.size() == 1) {
        // *this = *refseq_graph[targets.front()];
        // we have already done this because the first graph in the queue is this
    } else {
        // where we have multiple targets
        for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {
            // merge the variants into one graph
            VG& g = *refseq_graph[*t];
            combine(g);
        }
    }
}

void VG::sort(void) {
    if (size() <= 1) return;
    // Topologically sort, which orders and orients all the nodes.
    deque<NodeTraversal> sorted_nodes;
    topological_sort(sorted_nodes);
    deque<NodeTraversal>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.node_size() && n != sorted_nodes.end();
          ++i, ++n) {
        // Put the nodes in the order we got
        swap_nodes(graph.mutable_node(i), (*n).node);
    }
}

size_t VG::size(void) {
    return graph.node_size();
}

size_t VG::length(void) {
    size_t l;
    for_each_node([&l](Node* n) { l+=n->sequence().size(); });
    return l;
}

void VG::swap_nodes(Node* a, Node* b) {
    int aidx = node_index[a];
    int bidx = node_index[b];
    graph.mutable_node()->SwapElements(aidx, bidx);
    node_index[a] = bidx;
    node_index[b] = aidx;
}

Edge* VG::create_edge(NodeTraversal left, NodeTraversal right) {
    // Connect to the start of the left node if it is backward, and the end of the right node if it is backward.
    return create_edge(left.node->id(), right.node->id(), left.backward, right.backward);
}

Edge* VG::create_edge(NodeSide side1, NodeSide side2) {
    // Connect to node 1 (from start if the first side isn't an end) to node 2 (to end if the second side is an end)
    return create_edge(side1.node, side2.node, !side1.is_end, side2.is_end);
}

Edge* VG::create_edge(Node* from, Node* to, bool from_start, bool to_end) {
    return create_edge(from->id(), to->id(), from_start, to_end);
}

Edge* VG::create_edge(int64_t from, int64_t to, bool from_start, bool to_end) {
    //cerr << "creating edge " << from << "->" << to << endl;
    // ensure the edge (or another between the same sides) does not already exist
    Edge* edge = get_edge(NodeSide(from, !from_start), NodeSide(to, to_end));
    if (edge) {
        // The edge we want to make exists.
        return edge;
    }
    // if not, create it
    edge = graph.add_edge();
    edge->set_from(from);
    edge->set_to(to);
    // Only set the backwardness fields if they are true.
    if(from_start) edge->set_from_start(from_start);
    if(to_end) edge->set_to_end(to_end);
    set_edge(edge);
    edge_index[edge] = graph.edge_size()-1;
    //cerr << "created edge " << edge->from() << "->" << edge->to() << endl;
    return edge;
}

Edge* VG::get_edge(const NodeSide& side1, const NodeSide& side2) {
    auto e = edge_by_sides.find(minmax(side1, side2));
    if (e != edge_by_sides.end()) {
        return e->second;
    } else {
        return NULL;
    }
}

Edge* VG::get_edge(const pair<NodeSide, NodeSide>& sides) {
    return get_edge(sides.first, sides.second);
}

Edge* VG::get_edge(const NodeTraversal& left, const NodeTraversal& right) {
    // We went from the right side of left to the left side of right.
    // We used the end of left if if isn't backward, and we used the end of right if it is.
    return get_edge(NodeSide(left.node->id(), !left.backward),
                    NodeSide(right.node->id(), right.backward));
}

void VG::set_edge(Edge* edge) {
    // Note: there must not be an edge between these sides of these nodes already.
    if (!has_edge(edge)) {
        // Note that we might add edges to nonexistent nodes (like in VG::node_context()). That's just fine.

        // Add the edge to the index by node side (edges_on_start, edges_on_end, and edge_by_sides)
        index_edge_by_node_sides(edge);
    }
}

void VG::for_each_edge_parallel(function<void(Edge*)> lambda) {
    create_progress(graph.edge_size());
    int64_t completed = 0;
#pragma omp parallel for shared(completed)
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
        if (progress && completed++ % 1000 == 0) {
#pragma omp critical (progress_bar)
            update_progress(completed);
        }
    }
    destroy_progress();
}

void VG::for_each_edge(function<void(Edge*)> lambda) {
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
    }
}

void VG::destroy_edge(const NodeSide& side1, const NodeSide& side2) {
    destroy_edge(get_edge(side1, side2));
}

void VG::destroy_edge(const pair<NodeSide, NodeSide>& sides) {
    destroy_edge(sides.first, sides.second);
}


void VG::destroy_edge(Edge* edge) {
    //cerr << "destroying edge " << edge->from() << "->" << edge->to() << endl;

    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }

    // first remove the edge from the edge-on-node-side indexes.
    unindex_edge_by_node_sides(edge);

    // get the last edge index (lei) and this edge index (tei)
    int lei = graph.edge_size()-1;
    int tei = edge_index[edge];

    // erase this edge from the index by node IDs.
    // we'll fix up below
    edge_index.erase(edge);

    // Why do we check that lei != tei?
    //
    // It seems, after an inordinate amount of testing and probing,
    // that if we call erase twice on the same entry, we'll end up corrupting the hash_map
    //
    // So, if the element is already at the end of the table,
    // take a little break and just remove the last edge in graph

    // if we need to move the element to the last position in the array...
    if (lei != tei) {

        // get a pointer to the last element
        Edge* last = graph.mutable_edge(lei);

        // erase from our index
        edge_index.erase(last);

        // swap
        graph.mutable_edge()->SwapElements(tei, lei);

        // point to new position
        Edge* nlast = graph.mutable_edge(tei);

        // insert the new edge index position
        edge_index[nlast] = tei;

        // and fix edge indexes for moved edge object
        set_edge(nlast);

    }

    // drop the last position, erasing the node
    graph.mutable_edge()->RemoveLast();

    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }

}

void VG::unindex_edge_by_node_sides(const NodeSide& side1, const NodeSide& side2) {
    unindex_edge_by_node_sides(get_edge(side1, side2));
}

void VG::unindex_edge_by_node_sides(Edge* edge) {
    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }
    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }
    // erase from indexes

    //cerr << "erasing from indexes" << endl;

    //cerr << "Unindexing edge " << edge->from() << "<->" << edge->to() << endl;

    // Remove from the edge by node side pair index
    edge_by_sides.erase(NodeSide::pair_from_edge(edge));

    // Does this edge involve a change of relative orientation?
    bool relative_orientation = edge->from_start() != edge->to_end();

    // Un-index its from node, depending on whether it's attached to the start
    // or end.
    if(edge->from_start()) {
        // The edge is on the start of the from node, so remove it from the
        // start of the from node, with the correct relative orientation for the
        // to node.
        swap_remove(edges_start(edge->from()), make_pair(edge->to(), relative_orientation));
        // removing the sub-indexes if they are now empty
        // we must do this to maintain a valid structure
        if (edges_on_start[edge->from()].empty()) edges_on_start.erase(edge->from());

        //cerr << "Removed " << edge->from() << "-start to " << edge->to() << " orientation " << relative_orientation << endl;
    } else {
        // The edge is on the end of the from node, do remove it form the end of the from node.
        swap_remove(edges_end(edge->from()), make_pair(edge->to(), relative_orientation));
        if (edges_on_end[edge->from()].empty()) edges_on_end.erase(edge->from());

        //cerr << "Removed " << edge->from() << "-end to " << edge->to() << " orientation " << relative_orientation << endl;
    }

    if(edge->from() != edge->to() || edge->from_start() == edge->to_end()) {
        // Same for the to node, if we aren't just on the same node and side as with the from node.
        if(edge->to_end()) {
            swap_remove(edges_end(edge->to()), make_pair(edge->from(), relative_orientation));
            if (edges_on_end[edge->to()].empty()) edges_on_end.erase(edge->to());

            //cerr << "Removed " << edge->to() << "-end to " << edge->from() << " orientation " << relative_orientation << endl;
        } else {
            swap_remove(edges_start(edge->to()), make_pair(edge->from(), relative_orientation));
            if (edges_on_start[edge->to()].empty()) edges_on_start.erase(edge->to());

            //cerr << "Removed " << edge->to() << "-start to " << edge->from() << " orientation "
            //     << relative_orientation << endl;
        }
    }
}

void VG::index_edge_by_node_sides(Edge* edge) {

    // Generate sides, order them, and index the edge by them.
    edge_by_sides[NodeSide::pair_from_edge(edge)] = edge;

    // Index on ends appropriately depending on from_start and to_end.
    bool relative_orientation = edge->from_start() != edge->to_end();

    if(edge->from_start()) {
        edges_on_start[edge->from()].emplace_back(edge->to(), relative_orientation);
    } else {
        edges_on_end[edge->from()].emplace_back(edge->to(), relative_orientation);
    }

    if(edge->from() != edge->to() || edge->from_start() == edge->to_end()) {
        // Only index the other end of the edge if the edge isn't a self-loop on a single side.
        if(edge->to_end()) {
            edges_on_end[edge->to()].emplace_back(edge->from(), relative_orientation);
        } else {
            edges_on_start[edge->to()].emplace_back(edge->from(), relative_orientation);
        }
    }
}

Node* VG::get_node(int64_t id) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(id);
    if (n != node_by_id.end()) {
        return n->second;
    } else {
        throw runtime_error("No node " + to_string(id) + " in graph");
    }
}

Node* VG::create_node(string seq, int64_t id) {
    // create the node
    Node* node = graph.add_node();
    node->set_sequence(seq);
    // ensure we properly update the current_id that's used to generate new ids
    // unless we have a specified id
    if (id == 0) {
        if (current_id == 1) current_id = max_node_id()+1;
        node->set_id(current_id++);
    } else {
        node->set_id(id);
    }
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = graph.node_size()-1;
    //if (!is_valid()) cerr << "graph invalid" << endl;
    return node;
}

void VG::for_each_node_parallel(function<void(Node*)> lambda) {
    create_progress(graph.node_size());
    int64_t completed = 0;
#pragma omp parallel for schedule(dynamic,1) shared(completed)
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
        if (progress && completed++ % 1000 == 0) {
#pragma omp critical (progress_bar)
            update_progress(completed);
        }
    }
    destroy_progress();
}

void VG::for_each_node(function<void(Node*)> lambda) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
    }
}

void VG::for_each_connected_node(Node* node, function<void(Node*)> lambda) {
    // We keep track of nodes to visit.
    set<Node*> to_visit {node};
    // We mark all the nodes we have visited.
    set<Node*> visited;

    while(!to_visit.empty()) {
        // Grab some node to visit
        Node* visiting = *(to_visit.begin());
        to_visit.erase(to_visit.begin());

        // Visit it
        lambda(visiting);
        visited.insert(visiting);

        // Look at all its edges
        vector<Edge*> edges;
        edges_of_node(visiting, edges);
        for(Edge* edge : edges) {
            if(edge->from() != visiting->id() && !visited.count(get_node(edge->from()))) {
                // We found a new node on the from of the edge
                to_visit.insert(get_node(edge->from()));
            } else if(edge->to() != visiting->id() && !visited.count(get_node(edge->to()))) {
                // We found a new node on the to of the edge
                to_visit.insert(get_node(edge->to()));
            }
        }
    }
}

// a graph composed of this node and its edges
void VG::node_context(Node* node, VG& g) {
    // add the node
    g.add_node(*node);
    // and its edges
    vector<pair<int64_t, bool>>& start = edges_start(node->id());
    for (auto& e : start) {
        g.add_edge(*get_edge(NodeSide::pair_from_start_edge(node->id(), e)));
    }
    vector<pair<int64_t, bool>>& end = edges_end(node->id());
    for (auto& e : end) {
        g.add_edge(*get_edge(NodeSide::pair_from_end_edge(node->id(), e)));
    }
    // and its path members
    if (paths.has_node_mapping(node)) {
        auto& node_mappings = paths.get_node_mapping(node);
        for (auto& i : node_mappings) {
            g.paths.append_mapping(i.first, *i.second);
        }
    }
}

// a graph composed of this node and the edges that can be uniquely assigned to it
void VG::nonoverlapping_node_context(Node* node, VG& g) {
    // add the node
    g.add_node(*node);

    auto grab_edge = [&](Edge* e) {
        // What node owns the edge?
        int64_t owner_id = min(e->from(), e->to());
        if(node->id() == owner_id || !has_node(owner_id)) {
            // Either we are the owner, or the owner isn't in the graph to get serialized.
            g.add_edge(*e);
        }
    };

    // Go through all its edges
    vector<pair<int64_t, bool>>& start = edges_start(node->id());
    for (auto& e : start) {
        grab_edge(get_edge(NodeSide::pair_from_start_edge(node->id(), e)));
    }
    vector<pair<int64_t, bool>>& end = edges_end(node->id());
    for (auto& e : end) {
        grab_edge(get_edge(NodeSide::pair_from_end_edge(node->id(), e)));
    }

    // and its path members
    if (paths.has_node_mapping(node)) {
        auto& node_mappings = paths.get_node_mapping(node);
        for (auto& i : node_mappings) {
            g.paths.append_mapping(i.first, *i.second);
        }
    }
}

void VG::destroy_node(int64_t id) {
    destroy_node(get_node(id));
}

void VG::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    //cerr << "destroying node " << node->id() << " degrees " << start_degree(node) << ", " << end_degree(node) << endl;
    // noop on NULL/nonexistent node
    if (!has_node(node)) { return; }
    // remove edges associated with node
    set<pair<NodeSide, NodeSide>> edges_to_destroy;

    for(auto& other_end : edges_start(node)) {
        // Destroy all the edges on its start
        edges_to_destroy.insert(NodeSide::pair_from_start_edge(node->id(), other_end));
    }

    for(auto& other_end : edges_end(node)) {
        // Destroy all the edges on its end
        edges_to_destroy.insert(NodeSide::pair_from_end_edge(node->id(), other_end));
    }

    for (auto& e : edges_to_destroy) {
        //cerr << "Destroying edge " << e.first << ", " << e.second << endl;
        destroy_edge(e.first, e.second);
        //cerr << "Edge destroyed" << endl;
    }

    // assert cleanup
    edges_on_start.erase(node->id());
    edges_on_end.erase(node->id());

    //assert(start_degree(node) == 0);
    //assert(end_degree(node) == 0);

    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = graph.node_size()-1;
    int tni = node_index[node];

    if (lni != tni) {
        // swap this node with the last one
        Node* last = graph.mutable_node(lni);
        graph.mutable_node()->SwapElements(tni, lni);
        Node* nlast = graph.mutable_node(tni);

        // and fix up the indexes
        node_by_id[last->id()] = nlast;
        node_index.erase(last);
        node_index[nlast] = tni;
    }

    // remove this node (which is now the last one) and remove references from the indexes
    node_by_id.erase(node->id());
    node_index.erase(node);
    graph.mutable_node()->RemoveLast();
    //if (!is_valid()) { cerr << "graph is invalid after destroy_node" << endl; exit(1); }
}

void VG::remove_null_nodes(void) {
    vector<int64_t> to_remove;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node->id());
        }
    }
    for (vector<int64_t>::iterator n = to_remove.begin(); n != to_remove.end(); ++n) {
        destroy_node(*n);
    }
}

void VG::remove_null_nodes_forwarding_edges(void) {
    vector<Node*> to_remove;
    int i = 0;
    create_progress(graph.node_size()*2);
    for (i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node);
        }
        update_progress(i);
    }
    for (vector<Node*>::iterator n = to_remove.begin(); n != to_remove.end(); ++n, ++i) {
        remove_node_forwarding_edges(*n);
        update_progress(i);
    }
}

void VG::remove_node_forwarding_edges(Node* node) {
    // Grab all the nodes attached to our start, with true if the edge goes to their start
    vector<pair<int64_t, bool>>& start = edges_start(node);
    // Grab all the nodes attached to our end, with true if the edge goes to their end
    vector<pair<int64_t, bool>>& end = edges_end(node);

    // We instantiate the whole cross product first to avoid working on
    // references to the contents of containers we are modifying. This holds the
    // (node ID, relative orientation) pairs above.
    set<pair<pair<int64_t, bool>, pair<int64_t, bool>>> edges_to_create;

    // Make edges for the cross product of our start and end edges, making sure
    // to maintain relative orientation.
    for(auto& start_pair : start) {
        for(auto& end_pair : end) {
            // We already have the flags for from_start and to_end for the new edge.
            edges_to_create.emplace(start_pair, end_pair);
        }
    }

    for (auto& e : edges_to_create) {
        // make each edge we want to add
        create_edge(e.first.first, e.second.first, e.first.second, e.second.second);
    }

    // remove the node from paths
    if (paths.has_node_mapping(node)) {
        auto& node_mappings = paths.get_node_mapping(node);
        for (auto& p : node_mappings) {
            paths.remove_mapping(p.second);
        }
    }
    // delete the actual node
    destroy_node(node);
}

void VG::remove_orphan_edges(void) {
    set<pair<NodeSide, NodeSide>> edges;
    for_each_edge([this,&edges](Edge* edge) {
            if (!has_node(edge->from())
                || !has_node(edge->to())) {
                edges.insert(NodeSide::pair_from_edge(edge));
            }
        });
    for (auto edge : edges) {
        destroy_edge(edge);
    }
}

void VG::keep_paths(set<string>& path_names, set<string>& kept_names) {
    // Strategy:

    // For each node in our graph
    // If it has node mappings, look at them
    // For each that occurs along a path we want to include
    // Mark that path as kept
    // Mark the node as kept
    // Peek in the path left and right from that mapping
    // Mark the edges traversed as kept
    // Mark all nodes and edges not marked for keeping as for removal
    // Remove all nodes and edges marked for removal
    // Drop all paths from the graph not requested to be kept

    // Previous code here would only keep edges between nodes with adjacent IDs
    // in the set of selected nodes, which is in general incorrect.

    vector<Node*> nodes_to_keep;
    vector<Node*> nodes_to_remove;
    // Holds node sides in smaller-to-larger order, so we can search for edges in it.
    set<pair<NodeSide, NodeSide>> edges_to_keep;

    for_each_node([&](Node* node) {
        // If we don't see anything about this node in our paths, throw it out.
        bool to_keep = false;

        if(paths.has_node_mapping(node)) {
            // This node appears on some paths. Look for appearances on paths we like.
            for(auto& appearance : paths.get_node_mapping(node)) {
                if(path_names.count(appearance.first)) {
                    // We found an appearance of this node on a path we are keeping. It comes with a Mapping*.

                    // Mark the path and node as kept
                    kept_names.insert(appearance.first);
                    to_keep = true;

                    // Walk left along the path and keep the edge we traverse.
                    Mapping* left_neighbor = paths.traverse_left(appearance.second);

                    if(left_neighbor != nullptr) {
                        // We aren't the first thing in the path, we want to keep the edge to our left.
                        // It may not exist, but we can still ask to keep it.

                        // What other node do we go to?
                        int64_t neighbor_id = left_neighbor->position().node_id();

                        if(has_node(neighbor_id)) {
                            // Keep the edge if we actually have the other end.
                            // We know the other end is on this path and will be
                            // kept.

                            // We attach to the end of the previous node if it isn't
                            // backward along the path, and the end of this node if
                            // it is backward along the path.
                            edges_to_keep.insert(minmax(NodeSide(neighbor_id, !left_neighbor->is_reverse()),
                                                        NodeSide(node->id(), appearance.second->is_reverse())));
                        }
                    }

                    // We skip walking right along the path, because if the node
                    // we would find is in the graph, we're going to walk left
                    // from it eventually and find the same edge.
                }
            }
        }

        // Actually decide to keep or delete the node. Keep it if we ever saw it on a path.
        if (to_keep) {
            nodes_to_keep.push_back(node);
        } else {
            nodes_to_remove.push_back(node);
        }

    });

    // Mark any edges we don't keep for destruction
    set<pair<NodeSide, NodeSide>> edges_to_destroy;
    for_each_edge([this, &edges_to_keep, &edges_to_destroy](Edge* edge) {
            auto ep = NodeSide::pair_from_edge(edge);
            if (!edges_to_keep.count(ep)) {
                edges_to_destroy.insert(ep);
            }
        });

    // Destroy all the edges and nodes we don't want
    for (auto edge : edges_to_destroy) {
        destroy_edge(edge.first, edge.second);
    }
    for (auto node : nodes_to_remove) {
        destroy_node(node);
    }

    // Throw out all the paths data for paths we don't want to keep.
    paths.keep_paths(path_names);
}

void VG::keep_path(string& path_name) {
    set<string> s,k; s.insert(path_name);
    keep_paths(s, k);
}

// utilities
void VG::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    //cerr << "dividing node " << node->id() << " at " << pos << endl;

    if (pos < 0 || pos > node->sequence().size()) {
#pragma omp critical (cerr)
        {
            cerr << omp_get_thread_num() << ": cannot divide node " << node->id() << ":" << node->sequence()
                 << " -- position (" << pos << ") is less than 0 or greater than sequence length ("
                 << node->sequence().size() << ")" << endl;
            exit(1);
        }
    }

#ifdef debug
#pragma omp critical (cerr)
    cerr << omp_get_thread_num() << ": in divide_node " << pos << " of " << node->sequence().size() << endl;
#endif


    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    hash_map<int64_t, vector<int64_t> >::const_iterator e;
    // Create edges between the left node (optionally its start) and the right node (optionally its end)
    set<pair<pair<int64_t, bool>, pair<int64_t, bool>>> edges_to_create;

    // replace the connections to the node's start
    for(auto& e : edges_start(node)) {
        // Make an edge to the left node's start from wherever this edge went.
        edges_to_create.emplace(make_pair(e.first, e.second), make_pair(left->id(), false));
    }

    // make our right node
    right = create_node(node->sequence().substr(pos));

    // replace the connections to the node's end
    for(auto& e : edges_end(node)) {
        // Make an edge from the right node's end to wherever this edge went.
        edges_to_create.emplace(make_pair(right->id(), false), make_pair(e.first, e.second));
    }

    // create the edges here as otherwise we will invalidate the iterators
    for (auto& e : edges_to_create) {
        // Swizzle the from_start and to_end bits to the right place.
        create_edge(e.first.first, e.second.first, e.first.second, e.second.second);
    }

    // connect left to right. This edge always goes from end to start.
    create_edge(left, right);

    //cerr << "Dividing into " << left->id() << " and " << right->id() << endl;

    // divide paths
    // note that we can't do this (yet) for non-exact matching paths
    if (paths.has_node_mapping(node)) {
        auto& node_path_mapping = paths.get_node_mapping(node);
        // apply to left and right
        vector<Mapping*> to_divide;
        for (auto& pm : node_path_mapping) {
            string path_name = pm.first;
            Mapping* m = pm.second;
            to_divide.push_back(m);
        }
        for (auto m : to_divide) {
            // we have to divide the mapping
            string path_name = paths.mapping_path_name(m);
            Mapping l, r; divide_invariant_mapping(*m, l, r, pos, left, right);
            // with the mapping divided, insert the pieces where the old one was
            auto mpit = paths.remove_mapping(m);
            // insert right then left (insert puts *before* the iterator)
            mpit = paths.insert_mapping(mpit, path_name, r);
            mpit = paths.insert_mapping(mpit, path_name, l);
        }
    }

    destroy_node(node);

}

// for dividing a path of nodes with an underlying coordinate system
void VG::divide_path(map<long, int64_t>& path, long pos, Node*& left, Node*& right) {

    map<long, int64_t>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = get_node(target->second);

    // nothing to do
    if (node_pos == pos) {
        map<long, int64_t>::iterator n = target; --n;
        left = get_node(n->second);
        right = get_node(target->second);
    } else {
        // divide the target node at our pos
        int diff = pos - node_pos;
        divide_node(old, diff, left, right);
        // left
        path[node_pos] = left->id();
        // right
        path[pos] = right->id();
    }
}

void VG::nodes_prev(NodeTraversal node, vector<NodeTraversal>& nodes) {
    // Get the node IDs that attach to the left of this node, and whether we are
    // attached relatively forward (false) or backward (true)
    vector<pair<int64_t, bool>>& left_nodes = node.backward ? edges_end(node.node) : edges_start(node.node);
    for (auto& prev : left_nodes) {
        // Make a NodeTraversal that is an oriented description of the node attached to our relative left.
        // If we're backward, and it's in the same relative orientation as us, it needs to be backward too.
        nodes.emplace_back(get_node(prev.first), prev.second != node.backward);
    }
}

void VG::nodes_next(NodeTraversal node, vector<NodeTraversal>& nodes) {
    // Get the node IDs that attach to the right of this node, and whether we
    // are attached relatively forward (false) or backward (true)
    vector<pair<int64_t, bool>>& right_nodes = node.backward ? edges_start(node.node) : edges_end(node.node);
    for (auto& next : right_nodes) {
        // Make a NodeTraversal that is an oriented description of the node attached to our relative right.
        // If we're backward, and it's in the same relative orientation as us, it needs to be backward too.
        nodes.emplace_back(get_node(next.first), next.second != node.backward);
    }
}

int VG::node_count_prev(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_prev(n, nodes);
    return nodes.size();
}

int VG::node_count_next(NodeTraversal n) {
    vector<NodeTraversal> nodes;
    nodes_next(n, nodes);
    return nodes.size();
}

void VG::prev_kpaths_from_node(NodeTraversal node, int length, int edge_max, bool edge_bounding,
                               list<NodeTraversal> postfix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes) {
    if (length == 0) return;
    if (edge_bounding && edge_max <= 0) {
        // We recursed in here and hit the max edge depth. Complain to the caller.
        maxed_nodes(node);
        return;
    }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    postfix.push_front(node);
    // Get all the nodes left of this one
    vector<NodeTraversal> prev_nodes;
    nodes_prev(node, prev_nodes);
    if (prev_nodes.empty()) {
        // We can't go any lefter, so we produce this as a path
        list<NodeTraversal> new_path = postfix;
        paths.insert(new_path);
    } // implicit else
    for (NodeTraversal& prev : prev_nodes) {
        if (prev.node->sequence().size() < length) {
            prev_kpaths_from_node(prev,
                                  length - prev.node->sequence().size(),
                                  // Charge 1 against edge_max for every alternative edge we passed up
                                  edge_max - max(left_degree(node)-1, 0),
                                  // but only if we are using edge bounding
                                  edge_bounding,
                                  postfix, paths, maxed_nodes);
        } else {
            // create a path for this node
            list<NodeTraversal> new_path = postfix;
            new_path.push_front(prev);
            paths.insert(new_path);
        }
    }
}

void VG::next_kpaths_from_node(NodeTraversal node, int length, int edge_max, bool edge_bounding,
                               list<NodeTraversal> prefix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes) {
    if (length == 0) return;
    if (edge_bounding && edge_max <= 0) {
        maxed_nodes(node);
        return;
    }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    prefix.push_back(node);
    vector<NodeTraversal> next_nodes;
    nodes_next(node, next_nodes);
    if (next_nodes.empty()) {
        list<NodeTraversal> new_path = prefix;
        paths.insert(new_path);
    } // implicit else
    for (NodeTraversal& next : next_nodes) {
        if (next.node->sequence().size() < length) {
            next_kpaths_from_node(next,
                                  length - next.node->sequence().size(),
                                  // Charge 1 against edge_max for every alternative edge we passed up
                                  edge_max - max(right_degree(node)-1, 0),
                                  // but only if we are using edge bounding
                                  edge_bounding,
                                  prefix, paths, maxed_nodes);
        } else {
            // create a path for this node
            list<NodeTraversal> new_path = prefix;
            new_path.push_back(next);
            paths.insert(new_path);
        }
    }
}

// iterate over the kpaths in the graph, doing something

void VG::for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> prev_maxed,
                        function<void(NodeTraversal)> next_maxed,
                        function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    auto by_node = [k, edge_max, &lambda, &prev_maxed, &next_maxed, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node(by_node);
}

void VG::for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> prev_maxed,
                        function<void(NodeTraversal)> next_maxed,
                        function<void(size_t,Path&)> lambda) {
    auto by_node = [k, edge_max, &lambda, &prev_maxed, &next_maxed, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node(by_node);
}

// parallel versions of above
// this isn't by default because the lambda may have side effects
// that need to be guarded explicitly

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed,
                                 function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    auto by_node = [k, edge_max, &prev_maxed, &next_maxed, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node_parallel(by_node);
}

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> prev_maxed,
                                 function<void(NodeTraversal)> next_maxed,
                                 function<void(size_t,Path&)> lambda) {
    auto by_node = [k, edge_max, &prev_maxed, &next_maxed, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, prev_maxed, next_maxed, lambda);
    };
    for_each_node_parallel(by_node);
}

// per-node kpaths

void VG::for_each_kpath_of_node(Node* n, int k, int edge_max,
                                function<void(NodeTraversal)> prev_maxed,
                                function<void(NodeTraversal)> next_maxed,
                                function<void(size_t,Path&)> lambda) {
    auto apply_to_path = [&lambda, this](list<NodeTraversal>::iterator n, list<NodeTraversal>& p) {
        Path path = create_path(p);

        // Find the index of the node with a scan. We already did O(n) work to make the path.
        // TODO: is there a more efficient way?
        size_t index = 0;
        while(n != p.begin()) {
            --n;
            ++index;
        }

        lambda(index, path);
    };
    for_each_kpath_of_node(n, k, edge_max, prev_maxed, next_maxed, apply_to_path);
}

void VG::for_each_kpath_of_node(Node* node, int k, int edge_max,
                                function<void(NodeTraversal)> prev_maxed,
                                function<void(NodeTraversal)> next_maxed,
                                function<void(list<NodeTraversal>::iterator,list<NodeTraversal>&)> lambda) {
    // get left, then right
    set<list<NodeTraversal> > prev_paths;
    set<list<NodeTraversal> > next_paths;
    list<NodeTraversal> empty_list;
    
    prev_kpaths_from_node(NodeTraversal(node), k, edge_max, (edge_max != 0),
                          empty_list, prev_paths, prev_maxed);
    next_kpaths_from_node(NodeTraversal(node), k, edge_max, (edge_max != 0),
                          empty_list, next_paths, next_maxed);
    // now take the cross and give to the callback
    for (set<list<NodeTraversal> >::iterator p = prev_paths.begin(); p != prev_paths.end(); ++p) {
        for (set<list<NodeTraversal> >::iterator n = next_paths.begin(); n != next_paths.end(); ++n) {
            list<NodeTraversal> path = *p;

            // Find the iterator to this node in the list that will become the
            // path. We know it's the last thing in the prev kpath we made our path from.
            list<NodeTraversal>::iterator this_node = path.end(); this_node--;

            list<NodeTraversal>::const_iterator m = n->begin(); ++m; // skips current node, which is included in *p in the correct orientation
            while (m != n->end()) {
                // Add on all the other nodes from the next kpath.
                path.push_back(*m);
                ++m;
            }

            lambda(this_node, path);
        }
    }
}

void VG::kpaths_of_node(Node* node, set<list<NodeTraversal> >& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    auto collect_path = [&paths](list<NodeTraversal>::iterator n, list<NodeTraversal>& path) {
        paths.insert(path);
    };
    for_each_kpath_of_node(node, length, edge_max, prev_maxed, next_maxed, collect_path);
}

void VG::kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    set<list<NodeTraversal> > unique_paths;
    kpaths_of_node(node, unique_paths, length, edge_max, prev_maxed, next_maxed);
    for (auto& unique_path : unique_paths) {
        Path path = create_path(unique_path);
        paths.push_back(path);
    }
}

// aggregators, when a callback won't work

void VG::kpaths(set<list<NodeTraversal> >& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        kpaths_of_node(node, paths, length, edge_max, prev_maxed, next_maxed);
    }
}

void VG::kpaths(vector<Path>& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    set<list<NodeTraversal> > unique_paths;
    kpaths(unique_paths, length, edge_max, prev_maxed, next_maxed);
    for (auto& unique_path : unique_paths) {
        Path path = create_path(unique_path);
        paths.push_back(path);
    }
}

// path utilities
// these are in this class because attributes of the path (such as its sequence) are a property of the graph

Path VG::create_path(const list<NodeTraversal>& nodes) {
    Path path;
    for (const NodeTraversal& n : nodes) {
        Mapping* mapping = path.add_mapping();
        mapping->mutable_position()->set_node_id(n.node->id());
        // If the node is backward along this path, note it.
        if(n.backward) mapping->set_is_reverse(n.backward);
        // TODO: Is it OK if we say we're at a mapping at offset 0 of this node, backwards? Or should we offset ourselves to the end?
    }
    return path;
}

string VG::path_string(const list<NodeTraversal>& nodes) {
    string seq;
    for (const NodeTraversal& n : nodes) {
        if(n.backward) {
            seq.append(reverse_complement(n.node->sequence()));
        } else {
            seq.append(n.node->sequence());
        }
    }
    return seq;
}

string VG::path_string(Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        Mapping* m = path.mutable_mapping(i);
        Node* n = node_by_id[m->position().node_id()];
        if (m->is_reverse()) {
            seq.append(reverse_complement(n->sequence()));
        } else {
            seq.append(n->sequence());
        }
    }
    return seq;
}

void VG::expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded) {
    for (list<NodeTraversal>::const_iterator n = path.begin(); n != path.end(); ++n) {
        NodeTraversal node = *n;
        int s = node.node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(node);
        }
    }
}

void VG::expand_path(list<NodeTraversal>& path, vector<list<NodeTraversal>::iterator>& expanded) {
    for (list<NodeTraversal>::iterator n = path.begin(); n != path.end(); ++n) {
        int s = (*n).node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(n);
        }
    }
}

void VG::include(const Path& path) {
    // TODO we are failing in two basic ways
    // 1) failing to account for the current node we should be on after edits
    // 2) not maintaining the existing path architecture through edits
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& m = path.mapping(i);
        // TODO is it reversed?
        Node* n = get_node(m.position().node_id());
        int f = m.position().offset();
        int t = 0;
        for (int j = 0; j < m.edit_size(); ++j) {
            const Edit& edit = m.edit(j);
            if (edit.from_length() && edit.from_length() == edit.to_length()) {
                f += edit.from_length();
                t += edit.from_length();
            } else {
                // cut at f, and f + from_length
                Node* l=NULL;
                Node* r=NULL;
                Node* m=NULL;
                Node* c=NULL;
                if (edit.from_length() > 0) {
                    divide_node(n, f, l, m);
                    divide_node(m, edit.from_length(), m, r);
                    n = m; f = 0; t = 0;
                    // TODO there is a quirk (for sanity, efficiency later)
                    // if we "delete" over several nodes, we should join one path for the deletion
                    if (edit.to_length() == 0) {
                        // deletion
                        assert(edit.from_length());
                        create_edge(l, r);
                        n = r;
                    } else {
                        // swap/ SNP
                        c = create_node(edit.sequence());
                        create_edge(l, c);
                        create_edge(c, r);
                        n = r;
                    }
                } else if (!edit.sequence().empty()) { // create a node
                    divide_node(n, f, l, r);
                    n = r; f = 0; t = 0;
                    // insertion
                    assert(edit.to_length());
                    c = create_node(edit.sequence());
                    create_edge(l, c);
                    create_edge(c, r);
                } else {
                    // do nothing for soft clips, where we have to_length and unset from_length
                    f += edit.from_length();
                    t += edit.to_length();
                }
            }
        }
    }
    remove_null_nodes_forwarding_edges();
}

void VG::edit_node(int64_t node_id,
                   const vector<Mapping>& mappings,
                   map<pair<size_t, int64_t>, pair<set<Node*>, set<Node*>>>& cut_trans) {

    //cerr << "editing node " << node_id << endl;

    // assume these are all for the same node
    //int64_t node_id = mappings.front().node_id();
    // convert edits to a series of cut points
    set<int> cut_at;
    // map pairs of cut points to sequences of novel paths joining the cut points
    // which are themselves named
    map<pair<int, int>, map<string, string> > cut_seqs;
    Node* node = get_node(node_id);
    string node_seq = node->sequence();
    for (auto& mapping : mappings) {
        // check that we're really working on one node
        assert(mapping.position().node_id() == node_id);
        int offset = 0;
        const auto& nitr = mapping.metadata().info().find("path_name");
        assert(nitr != mapping.metadata().info().end());
        const string& name = nitr->second.str();
        for (int i = 0; i < mapping.edit_size(); ++i) {
            const Edit& edit = mapping.edit(i);
            //edits[offset].push_back(edit);
            int end = offset + edit.from_length();
            //if (edit_is_match(edit) || edit_is_softclip(edit)) {
            if (edit.sequence().empty() && edit.from_length() == edit.to_length() ||
                (i == 0 || i == mapping.edit_size() - 1) && edit.from_length() == 0) {
                //cerr << "ignoring soft clip" << endl;
                // soft clip, ignore
                // match, ignore
            } else { //if (edit_is_deletion(edit)) {
                //cerr << "o/e " << offset << " " << end << endl;
                cut_seqs[make_pair(offset, end)][name] = edit.sequence();
                cut_at.insert(offset);
                cut_at.insert(end);
            }
            offset = end;
        }
    }
    /*
    for (auto cut : cut_at) {
        cerr << "cut_at: " << cut << endl;
    }
    for (auto cut : cut_seqs) {
        cerr << "cut_seq: " << cut.first.first << "," << cut.first.second << " ";
        for (auto s : cut.second) {
            cerr << " " << s.first << "," << s.second << " ";
        }
        cerr << endl;
    }
    */


    // tricky bit
    // we need to put the node pointers into structures to track the sides of each cut
    // however, as we make the subsequent cut we end up changing the nodes destructively

    Node* n = get_node(node_id);
    Node* l = NULL;
    Node* r = NULL;
    int offset = 0;

    map<int, tuple<set<Node*>, set<Node*>, set<Node*> > > cuts;

    // make the cuts
    set<Node*> emptyset;
    for (auto cut : cut_at) {
        // don't cut if it would be at the end of the sequence
        if (n->sequence().size()) {
            divide_node(n, cut-offset, l, r);
        }
        set<Node*> left{l};
        set<Node*> right{r};
        cuts[cut] = make_tuple(left, emptyset, right);
        // now that we've made the cuts record the mapping between the old node id and the new
        // start and ends

        // adjust last cut (we just cut up the second node in the pair)
        if (offset > 0) {
            get<2>(cuts[offset]) = left;
        }
        offset = cut;
        n = r;
    }

    // add the novel seqs
    // and join deletions
    for (auto& cs : cut_seqs) {
        // get the left and right
        int f = cs.first.first; // we go from after the left node at this cut
        int t = cs.first.second; // to before the second node in this cut
        auto& left = get<0>(cuts[f]);
        auto& center = get<1>(cuts[f]);
        auto& right = get<2>(cuts[t]);
        //cerr << f << "," << t << " " << left.size() << "," << center.size() << "," << right.size() << endl;
        for (auto& p : cs.second) {
            auto& name = get<0>(p);
            auto& seq = get<1>(p);
            if (!seq.empty()) {
                Node* c = create_node(seq);
                // add to path
                paths.append_mapping(name, c->id(), false);
                center.insert(c);
            } else {

                // todo todo!
                // record the deletion
                // ?: does the deletion extend past the node boundary into another node
                // if so, we should record where we end up in the deletion tables
                // so that we can fix up deletions externally after both ends are created
                //cerr << "seq is empty" << endl;

                // if the deletion is local to the node, simply implement it here:
                /*
                if (t < node_seq.size()) {
                    for (auto& lp : left) {
                        for (auto& rp : right) {
                            cerr << "creating edge!" << endl;
                            create_edge(lp, rp);
                        }
                    }
                }
                */
            }
        }
    }

    // now we can record the cut pattern in

    // join in the novel seqs
    // remember: map<pair<int, int>, vector<string> > cut_seqs;
    for (auto& cs : cut_seqs) {
        // get the left and right
        int f = cs.first.first; // we go from after the left node at this cut
        int t = cs.first.second; // to before the second node in this cut
        auto& left = get<0>(cuts[f]);
        auto& center = get<1>(cuts[f]);
        auto& right = get<2>(cuts[t]);
        //cerr << "cut seq " << f << ":" << t << " "
        //<< left.size() << "," << center.size() << "," << right.size() << endl;
        for (auto& seq : cs.second) {
            for (auto& lp : left) {
                for (auto& rp : right) {
                    if (!center.empty()) {
                        for (auto& cp : center) {
                            // make the edge across the left cut
                            create_edge(lp, cp);
                            // record the cut translation (used for stitching deletions)
                            auto& p1 = cut_trans[make_pair(node_id, f)];
                            p1.first.insert(lp); p1.second.insert(cp);
                            //cerr << *p1.first.begin() << " " << *p1.second.begin() << endl;
                            // make the edge across the right cut
                            create_edge(cp, rp);
                            // record the cut translation ...
                            auto& p2 = cut_trans[make_pair(node_id, t)];
                            p2.first.insert(cp); p2.second.insert(rp);
                            //cerr << *p2.first.begin() << " " << *p2.second.begin() << endl;
                        }
                    } else {
                        // simpler case when we don't have a center node
                        // make the edge across the single cut
                        // and record the translation
                        auto& p1 = cut_trans[make_pair(node_id, f)];
                        p1.first.insert(lp);
                        p1.second.insert(rp);
                        auto& p2 = cut_trans[make_pair(node_id, t)];
                        p2.first.insert(lp);
                        p2.second.insert(rp);
                        //cerr << "cut trans " << *p.first.begin() << " " << *p.second.begin() << endl;
                        //cerr << node_id << ":" << f << endl;
                    }
                }
            }
        }
    }

    // now we should have incorporated the edits against a given node
    // deletions spanning multiple nodes must be handled externally using cut_trans
}

void VG::edit(const vector<Path>& paths) {
    map<pair<int64_t, size_t>, pair<int64_t, size_t> > del_f;
    map<pair<int64_t, size_t>, pair<int64_t, size_t> > del_t;
    map<int64_t, vector<Mapping>> mappings; // by node
    bool in_del = false;
    pair<int64_t, size_t> del_start;
    for (auto& path : paths) {
        for (int i = 0; i < path.mapping_size(); ++i) {
            Mapping mapping = path.mapping(i);
            Info info; info.set_str(path.name());
            (*mapping.mutable_metadata()->mutable_info())["path_name"] = info;
            mappings[mapping.position().node_id()].push_back(mapping);

            assert(!mapping_is_total_deletion(mapping)); // not handled
            // ^^ these should be stripped out using simplify_deletions(path)
            int64_t node_id = mapping.position().node_id();

            // find internal deletions
            size_t ioff = 0;
            for (int j = 0; j < mapping.edit_size(); ++j) {
                auto& edit = mapping.edit(j);
                // is it an internal deletion?
                if (edit_is_deletion(edit) && j > 0 && j < mapping.edit_size()-1) {
                    auto dstart = make_pair(node_id, ioff);
                    auto dend =   make_pair(node_id, ioff+edit.from_length());
                    del_f[dstart] = dend;
                    del_t[dend] = dstart;
                }
                ioff += edit.from_length();
            }

            // not just ones that jump nodes

            if (!in_del) {
                if (mapping_starts_in_deletion(mapping) && i > 0) {
                    // this deletion should extend from the end of the last node to this one
                    int64_t last_id = path.mapping(i-1).position().node_id();
                    del_start = make_pair(last_id, get_node(last_id)->sequence().size());
                }
                if (mapping_ends_in_deletion(mapping)) {
                    size_t del_off = get_node(node_id)->sequence().size()
                        - mapping.edit(mapping.edit_size()-1).from_length();
                    del_start = make_pair(node_id, del_off);
                    in_del = true;
                }
            } else {
                pair<int64_t, size_t> del_end;
                if (mapping_starts_in_deletion(mapping)) {
                    size_t end_off = mapping.edit(0).from_length();
                    del_end = make_pair(node_id, end_off);
                    in_del = false;
                } else {
                    del_end = make_pair(node_id, 0);
                    in_del = false;
                }
                // record the deletion
                //cerr << "del_start " << del_start.first << ":" << del_start.second << " -> " << del_end.first << ":" << del_end.second << endl;
                del_f[del_start] = del_end;
                del_t[del_end] = del_start;
            }
        }
    }
    edit(mappings, del_f, del_t);
}

// mappings sorted by node id
void VG::edit(const map<int64_t, vector<Mapping> >& mappings,
              map<pair<int64_t, size_t>, pair<int64_t, size_t> >& del_f,
              map<pair<int64_t, size_t>, pair<int64_t, size_t> >& del_t) {
    // we are adding the edits to the graph
    // doing this correctly requires us to keep track of where the other end of
    // deletions of the starts and ends (or whole) nodes would land
    // in this way we can produce edges that jump between the

    // deletions are not node local
    // so how do we deal with them?
    //  - we will need to process the graph in a partially-ordered way
    //  - or at least the region where they happen
    //  maybe the easiest thing to do is to try to store them in the graph
    // we can try to convert them into a node local event
    // this would be easy using node
    map<int64_t, pair<Node*, Node*> > node_end_map;
    map<pair<size_t, int64_t>, pair<set<Node*>, set<Node*> > > cut_trans;
    for (auto& node : mappings) {
        //cerr << "edit: " << node.first << endl;
        edit_node(node.first, node.second, cut_trans);
    }

    // now resolve the deletions
    for (auto& d : del_f) {
        auto& f = d.first;
        //cerr << "del" << endl;
        //cerr << f.first << ":" << f.second << endl;
        auto& t = d.second;
        //cerr << t.first << ":" << t.second << endl;
        auto& l = cut_trans[f];
        auto& r = cut_trans[t];
        // connect the left side of the from to the to-side of the to
        for (auto& ln : l.first) {
            for (auto& rn : r.second) {
                create_edge(ln, rn);
            }
        }
    }
    remove_null_nodes_forwarding_edges();
}

void VG::node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start) {
    int i = 0;
    for (list<NodeTraversal>::const_iterator n = path.begin(); n != path.end(); ++n) {
        node_start[(*n).node] = i;
        int l = (*n).node->sequence().size();
        i += l;
    }
}

void VG::node_starts_in_path(list<NodeTraversal>& path,
                             map<NodeTraversal*, int>& node_start) {
    int i = 0;
    for (list<NodeTraversal>::iterator n = path.begin(); n != path.end(); ++n) {
        node_start[&(*n)] = i;
        int l = (*n).node->sequence().size();
        i += l;
    }
}

void VG::kpaths_of_node(int64_t node_id, vector<Path>& paths, int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(node_id);
    if (n != node_by_id.end()) {
        Node* node = n->second;
        kpaths_of_node(node, paths, length, edge_max, prev_maxed, next_maxed);
    }
}

// not quite right -- only for exact matching paths
string VG::path_sequence(const Path& path) {
    string sequence;
    for (int i = 0; i < path.mapping_size(); ++i) {
        sequence.append(node_by_id[path.mapping(i).position().node_id()]->sequence());
    }
    return sequence;
}

string VG::random_read(int length, mt19937& rng, int64_t min_id, int64_t max_id, bool either_strand) {
    uniform_int_distribution<int64_t> int64_dist(min_id, max_id);
    int64_t id = int64_dist(rng);
    // We start at the node in its local forward orientation
    NodeTraversal node(get_node(id), false);
    int32_t start_pos = 0;
    if (node.node->sequence().size() > 1) {
        uniform_int_distribution<uint32_t> uint32_dist(0,node.node->sequence().size()-1);
        start_pos = uint32_dist(rng);
    }
    string read = node.node->sequence().substr(start_pos);
    while (read.size() < length) {
        // pick a random downstream node
        vector<NodeTraversal> next_nodes;
        nodes_next(node, next_nodes);
        if (next_nodes.empty()) break;
        uniform_int_distribution<int> next_dist(0, next_nodes.size()-1);
        node = next_nodes.at(next_dist(rng));
        // Put in the node sequence in the correct relative orientation
        read.append(node.backward ? reverse_complement(node.node->sequence()) : node.node->sequence());
    }
    read = read.substr(0, length);
    uniform_int_distribution<int> binary_dist(0, 1);
    if (either_strand && binary_dist(rng) == 1) {
        // We can flip to the other strand (i.e. node's local reverse orientation).
        return reverse_complement(read);
    } else {
        return read;
    }
}

bool VG::is_valid(void) {

    if (node_by_id.size() != graph.node_size()) {
        cerr << "graph invalid: node count is not equal to that found in node by-id index" << endl;
        return false;
    }

    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (node_by_id.find(n->id()) == node_by_id.end()) {
            cerr << "graph invalid: node " << n->id() << " missing from by-id index" << endl;
            return false;
        }
    }

    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        int64_t f = e->from();
        int64_t t = e->to();

        //cerr << "edge " << e << " " << e->from() << "->" << e->to() << endl;

        if (node_by_id.find(f) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " (" << f << "->" << t << ") cannot find node (from) " << f << endl;
            return false;
        }
        if (node_by_id.find(t) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " (" << f << "->" << t << ") cannot find node (to) " << t << endl;
            return false;
        }

        if (!edges_on_start.count(f) && !edges_on_end.count(f)) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in either index for 'from' node " << f << endl;
            return false;
        }

        if (!edges_on_start.count(t) && !edges_on_end.count(t)) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in either index for 'to' node " << t << endl;
            return false;
        }

    }

    for (pair<const int64_t, vector<pair<int64_t, bool>>>& start_and_edges : edges_on_start) {
        for (auto& edge_destination : start_and_edges.second) {
            // We're on the start, so we go to the end if we aren't a reversing edge
            Edge* e = get_edge(NodeSide::pair_from_start_edge(start_and_edges.first, edge_destination));
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if(start_and_edges.first != e->to() && start_and_edges.first != e->from()) {
                // It needs to be attached to the node we looked up
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't have start-indexed node in " << start_and_edges.first << "<->"
                     << edge_destination.first << endl;
                return false;
            }
            if(edge_destination.first != e->to() && edge_destination.first != e->from()) {
                // It also needs to be attached to the node it says it goes to
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't have non-start-indexed node in " << start_and_edges.first << "<->"
                     << edge_destination.first << endl;
                return false;
            }
            if(!((start_and_edges.first == e->to() && !e->to_end()) ||
                 (start_and_edges.first == e->from() && e->from_start()))) {
            
                // The edge needs to actually attach to the start of the node we looked it up for.
                // So at least one of its ends has to be to the start of the correct node.
                // It may also be attached to the end.
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't attach to start of " << start_and_edges.first << endl;
                return false;
            }
            if (!has_node(e->from())) {
                cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
            if (!has_node(e->to())) {
                cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
        }
    }

     for (pair<const int64_t, vector<pair<int64_t, bool>>>& end_and_edges : edges_on_end) {
        for (auto& edge_destination : end_and_edges.second) {
            Edge* e = get_edge(NodeSide::pair_from_end_edge(end_and_edges.first, edge_destination));
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if(end_and_edges.first != e->to() && end_and_edges.first != e->from()) {
                // It needs to be attached to the node we looked up
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't have end-indexed node in " << end_and_edges.first << "<->"
                     << edge_destination.first << endl;
                return false;
            }
            if(edge_destination.first != e->to() && edge_destination.first != e->from()) {
                // It also needs to be attached to the node it says it goes to
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't have non-end-indexed node in " << end_and_edges.first << "<->"
                     << edge_destination.first << endl;
                return false;
            }
            if(!((end_and_edges.first == e->to() && e->to_end()) ||
                 (end_and_edges.first == e->from() && !e->from_start()))) {
            
                // The edge needs to actually attach to the end of the node we looked it up for.
                // So at least one of its ends has to be to the end of the correct node.
                // It may also be attached to the start.
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " doesn't attach to end of " << end_and_edges.first << endl;
                return false;
            }
            if (!has_node(e->from())) {
                cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
            if (!has_node(e->to())) {
                cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
        }
    }

    return true;
}

void VG::to_dot(ostream& out, vector<Alignment> alignments) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    out << "    rankdir=LR;" << endl;
//    out << "    splines=line;" << endl;
    out << "    smoothType=spring;" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        auto node_paths = paths.of_node(n->id());
        if (!paths.empty() && node_paths.empty()) {
            out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\",fontcolor=red,color=red,shape=box];" << endl;
        } else {
            out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\",shape=box];" << endl;
        }
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        auto from_paths = paths.of_node(e->from());
        auto to_paths = paths.of_node(e->to());
        set<string> both_paths;
        std::set_intersection(from_paths.begin(), from_paths.end(),
                              to_paths.begin(), to_paths.end(),
                              std::inserter(both_paths, both_paths.begin()));
        // are both nodes in the same path?
        bool in_path = !(!paths.empty()
                        && (both_paths.empty()
                            || !paths.are_consecutive_nodes_in_path(e->from(), e->to(),
                                                                    *both_paths.begin())));

        // Is the edge in the "wrong" direction for rank constraints?
        bool is_backward = e->from_start() && e->to_end();

        if(is_backward) {
            // Flip the edge around and write it forward.
            Edge* original = e;
            e = new Edge();

            e->set_from(original->to());
            e->set_from_start(!original->to_end());

            e->set_to(original->from());
            e->set_to_end(!original->from_start());
        }

        // display what kind of edge we have using different edge head and tail styles
        // depending on if the edge comes from the start or not
        out << "    " << e->from() << " -> " << e->to();
        out << " [dir=both,";
        if (!in_path) {
            out << "color=red,";
        }
        if (e->from_start()) {
            out << "arrowtail=none,";
            out << "tailport=sw,";
        } else {
            out << "arrowtail=none,";
            out << "tailport=ne,";
        }
        if (e->to_end()) {
            out << "arrowhead=normal,";
            out << "headport=se";
        } else {
            out << "arrowhead=normal,";
            out << "headport=nw";
        }
        out << "];" << endl;

        if(is_backward) {
            // We don't need this duplicate edge
            delete e;
        }
    }
    // add nodes for the alignments and link them to the nodes they match
    int alnid = max_node_id()+1;
    for (auto& aln : alignments) {
        // check direction
        if (!aln.is_reverse()) {
            out << "    " << alnid << " [label=\"+""\",fontcolor=green];" << endl;
            out << "    " << alnid << " -> " << alnid+1 << " [dir=none,color=green];" << endl;
        } else {
            out << "    " << alnid << " [label=\"-""\",fontcolor=purple];" << endl;
            out << "    " << alnid << " -> " << alnid+1 << " [dir=none,color=purple];" << endl;
        }
        alnid++;
        for (int i = 0; i < aln.path().mapping_size(); ++i) {
            const Mapping& m = aln.path().mapping(i);
            //void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
            //string cigar_string(vector<pair<int, char> >& cigar);
            vector<pair<int, char> > cigar;
            mapping_cigar(m, cigar);
            stringstream mapid;
            const string& nodestr = get_node(m.position().node_id())->sequence();
            string mstr = mapping_string(nodestr, m);
            //mapid << alnid << ":" << m.position().node_id() << ":" << cigar_string(cigar);
            mapid << cigar_string(cigar) << ":" << mstr;
            // determine sequence of this portion of the alignment
            // set color based on cigar/mapping relationship
            string nstr;
            if (i == 0) {
                nstr = nodestr.substr(m.position().offset());
            } else if (i == aln.path().mapping_size()-1) {
                nstr = nodestr.substr(0, mapping_from_length(m));
            } else {
                nstr = nodestr;
            }
            string color = "blue";
            if (mstr != nstr) { // some mismatch, indicate with orange color
                color = "orange";
            }
            out << "    " << alnid << " [label=\"" << mapid.str() << "\",fontcolor=" << color << "];" << endl;
            if (i > 0) {
                out << "    " << alnid-1 << " -> " << alnid << "[dir=none,color=" << color << "];" << endl;
            }
            out << "    " << alnid << " -> " << m.position().node_id() << "[dir=none,color=" << color << ", style=invis];" << endl;
            out << "    { rank = same; " << alnid << "; " << m.position().node_id() << "; };" << endl;
            //out << "    " << m.position().node_id() << " -- " << alnid << "[color=" << color << ", style=invis];" << endl;
            alnid++;
        }
        if (!aln.is_reverse()) {
            out << "    " << alnid << " [label=\"-""\",fontcolor=purple];" << endl;
            out << "    " << alnid-1 << " -> " << alnid << " [dir=none,color=purple];" << endl;
        } else {
            out << "    " << alnid << " [label=\"+""\",fontcolor=green];" << endl;
            out << "    " << alnid-1 << " -> " << alnid << " [dir=none,color=green];" << endl;
        }
        alnid++;
    }
    out << "}" << endl;
}

void VG::to_gfa(ostream& out) {
    map<int64_t, vector<string> > sorted_output;
    out << "H" << "\t" << "HVN:Z:1.0" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        stringstream s;
        s << "S" << "\t" << n->id() << "\t" << n->sequence() << "\n";
        auto& node_mapping = paths.get_node_mapping(n->id());
        set<Mapping*> seen;
        for (auto& p : node_mapping) {
            if (seen.count(p.second)) continue;
            else seen.insert(p.second);
            const Mapping& mapping = *p.second;
            string cigar;
            if (mapping.edit_size() > 0) {
                vector<pair<int, char> > cigarv;
                mapping_cigar(mapping, cigarv);
                cigar = cigar_string(cigarv);
            } else {
                // empty mapping edit implies perfect match
                stringstream cigarss;
                cigarss << n->sequence().size() << "M";
                cigar = cigarss.str();
            }
            string orientation = mapping.is_reverse() ? "-" : "+";
            s << "P" << "\t" << n->id() << "\t" << p.first << "\t"
              << orientation << "\t" << cigar << "\n";
        }
        sorted_output[n->id()].push_back(s.str());
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        stringstream s;
        s << "L" << "\t" << e->from() << "\t"
          << (e->from_start() ? "-" : "+") << "\t"
          << e->to() << "\t"
          << (e->to_end() ? "-" : "+") << "\t"
          << "0M" << endl;
        sorted_output[e->from()].push_back(s.str());
    }
    for (auto& chunk : sorted_output) {
        for (auto& line : chunk.second) {
            out << line;
        }
    }
}

void VG::destroy_alignable_graph(void) {
    if (gssw_aligner != NULL) {
        delete gssw_aligner;
    }
}

void VG::connect_node_to_nodes(Node* node, vector<Node*>& nodes, bool from_start) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right, unless instructed otherwise
        create_edge(node, (*n), from_start, false);
    }
}

void VG::connect_nodes_to_node(vector<Node*>& nodes, Node* node, bool to_end) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right, unless instructed otherwise
        create_edge((*n), node, false, to_end);
    }
}

void VG::connect_node_to_nodes(NodeTraversal node, vector<NodeTraversal>& nodes) {
    for (vector<NodeTraversal>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right
        create_edge(node, (*n));
    }
}

void VG::connect_nodes_to_node(vector<NodeTraversal>& nodes, NodeTraversal node) {
    for (vector<NodeTraversal>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        // Connect them left to right
        create_edge((*n), node);
    }
}

// join all subgraphs together to a "null" head node
Node* VG::join_heads(void) {
    // Find the head nodes
    vector<Node*> heads;
    head_nodes(heads);

    // Then create the new node (so it isn't picked up as a head)
    current_id = max_node_id()+1;
    Node* root = create_node("N");

    // Wire it to all the heads and return
    connect_node_to_nodes(root, heads);
    return root;
}

void VG::join_heads(Node* node, bool from_start) {
    vector<Node*> heads;
    head_nodes(heads);
    
    // If the node we have been given shows up as a head, remove it.
    for(auto i = heads.begin(); i != heads.end(); ++i) {
        if(*i == node) {
            heads.erase(i);
            break;
        }
    }
    
    connect_node_to_nodes(node, heads, from_start);
}

void VG::join_tails(Node* node, bool to_end) {
    vector<Node*> tails;
    tail_nodes(tails);
    
    // If the node we have been given shows up as a tail, remove it.
    for(auto i = tails.begin(); i != tails.end(); ++i) {
        if(*i == node) {
            tails.erase(i);
            break;
        }
    }
    
    connect_nodes_to_node(tails, node, to_end);
}

void VG::add_start_and_end_markers(int length, char start_char, char end_char,
                                   Node*& head_node, Node*& tail_node,
                                   int64_t head_id, int64_t tail_id) {

    // first do the head

    if(head_node == nullptr) {
        // We get to create the head node
        string start_string(length, start_char);
        head_node = create_node(start_string, head_id);
    } else {
        // We got a head node
        add_node(*head_node);
    }
    // Attach the new head node to all the existing heads
    join_heads(head_node);

    // then the tail
    if(tail_node == nullptr) {
        string end_string(length, end_char);
        tail_node = create_node(end_string, tail_id);
    } else {
        add_node(*tail_node);
    }
    join_tails(tail_node);
}

void VG::add_single_start_end_marker(int length, char start_char, Node*& head_tail_node, int64_t id) {
    // This set will hold all the nodes we haven't attached yet. But we don't
    // want it to hold the head_tail_node, so we fill it in now.
    set<Node*> unattached;
    for_each_node([&](Node* node) {
        unattached.insert(node);
    });

    if(head_tail_node == nullptr) {
        // We get to create the node. In its forward orientation it's the start node, so we use the start character.
        string start_string(length, start_char);
        head_tail_node = create_node(start_string, id);
    } else {
        // We got a node to use
        add_node(*head_tail_node);
    }

    // We handle the head and tail joining ourselves so we can do connected components.
    vector<Node*> heads;
    head_nodes(heads);
    vector<Node*> tails;
    tail_nodes(tails);

    for(Node* head : heads) {
        if(unattached.count(head)) {
            // This head is unconnected.

            // Mark everything it's attached to as attached
            for_each_connected_node(head, [&](Node* node) {
                unattached.erase(node);
            });
        }

        // Tie it to the start/stop node
        create_edge(head_tail_node, head);
    }

    for(Node* tail : tails) {
        if(unattached.count(tail)) {
            // This tail is unconnected.

            // Mark everything it's attached to as attached
            for_each_connected_node(tail, [&](Node* node) {
                unattached.erase(node);
            });
        }

        // Tie it to the start/stop node, going from the end of the tail to the
        // end of the start/stop node.
        create_edge(tail, head_tail_node, false, true);
    }

    // Find the connected components that aren't attached, if any.
    while(!unattached.empty()) {
        // Grab and attach some unattached node
        Node* to_attach = *(unattached.begin());

        // Mark everything it's attached to as attached
        for_each_connected_node(to_attach, [&](Node* node) {
            unattached.erase(node);
        });

        // Add the edge
        create_edge(head_tail_node, to_attach);
#ifdef debug
        cerr << "Broke into disconnected component at " << to_attach->id() << endl;
#endif
    }
}

Alignment& VG::align(Alignment& alignment) {

    set<int64_t> flipped_nodes;
    orient_nodes_forward(flipped_nodes);

    // to be completely aligned, the graph's head nodes need to be fully-connected to a common root
    Node* root = join_heads();
    // Put the nodes in sort order within the graph
    sort();

    gssw_aligner = new GSSWAligner(graph);
    gssw_aligner->align(alignment);
    delete gssw_aligner;
    gssw_aligner = NULL;

    destroy_node(root);

    flip_nodes(alignment, flipped_nodes, [this](int64_t node_id) {
            // We need to feed in the lengths of nodes, so the offsets in the alignment can be updated.
            return get_node(node_id)->sequence().size();
        });

    return alignment;
}

Alignment VG::align(string& sequence) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment);
}

void VG::for_each_kmer_parallel(int kmer_size,
                                int edge_max,
                                function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                                int stride,
                                bool allow_dups,
                                bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, true, stride, allow_dups, allow_negatives);
}

void VG::for_each_kmer(int kmer_size,
                       int edge_max,
                       function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                       int stride,
                       bool allow_dups,
                       bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, false, stride, allow_dups, allow_negatives);
}

void VG::for_each_kmer_of_node(Node* node,
                               int kmer_size,
                               int edge_max,
                               function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                               int stride,
                               bool allow_dups,
                               bool allow_negatives) {
    _for_each_kmer(kmer_size, edge_max, lambda, false, stride, allow_dups, allow_negatives, node);
}

void VG::_for_each_kmer(int kmer_size,
                        int edge_max,
                        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives,
                        Node* node) {

    // use an LRU cache to clean up duplicates over the last 1mb
    // use one per thread so as to avoid contention
    // If we aren't starting a parallel kmer iteration from here, just fill in 0.
    // TODO: How do we know this is big enough?
    map<int, LRUCache<string, bool>* > lru;
#pragma omp parallel
    {
#pragma omp single
        for (int i = 0; i < (parallel ? omp_get_num_threads() : 1); ++i) {
            lru[i] = new LRUCache<string, bool>(100000);
        }
    }
    // constructs the cache key
    // experiment -- use a struct here
    // We deduplicate kmers based on where they start, where they end (optionally), and where they are viewed from.
    auto make_cache_key = [](string& kmer,
                             int64_t start_node, int start_pos,
                             int64_t view_node, int view_pos,
                             int64_t end_node, int end_pos) -> string {
        string cache_key = kmer;
        cache_key.resize(kmer.size() + sizeof(Node*) + 3*sizeof(int64_t) + 3*sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size(), &start_node, sizeof(int64_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+1*sizeof(int64_t), &start_pos, sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size()+1*sizeof(int64_t)+1*sizeof(int), &view_node, sizeof(int64_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+2*sizeof(int64_t)+1*sizeof(int), &view_pos, sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size()+2*sizeof(int64_t)+2*sizeof(int), &end_node, sizeof(int64_t));
        memcpy((char*)cache_key.c_str()+kmer.size()+3*sizeof(int64_t)+2*sizeof(int), &end_pos, sizeof(int));
        return cache_key;
    };

    auto handle_path = [this,
                        &lambda,
                        kmer_size,
                        stride,
                        allow_dups,
                        allow_negatives,
                        &lru,
                        &make_cache_key,
                        &parallel,
                        &node](list<NodeTraversal>::iterator forward_node, list<NodeTraversal>& forward_path) {
#ifdef debug
        cerr << "Handling path: " << endl;
        for(auto& traversal : forward_path) {
            cerr << "\t" << traversal.node->id() << " " << traversal.backward << endl;
        }
#endif

        // Reserve a place for our reversed path, if we find ourselves needing to flip it around for a negative offset kmer.
        list<NodeTraversal> reversed_path;

        // And one for the reversed version of this NodeTraversal on that path
        list<NodeTraversal>::iterator reversed_node;


        // expand the path into a vector :: 1,1,1,2,2,2,2,3,3 ... etc.
        // this makes it much easier to quickly get all the node matches of each kmer

        // We expand out as iterators in the path list, so we can get the
        // traversal but also distinguish different node instances in the path.
        vector<list<NodeTraversal>::iterator> node_by_path_position;
        expand_path(forward_path, node_by_path_position);

        // Go get the cache for this thread if _for_each_kmer launched threads,
        // and the only one we made (thread 0) if we're running this
        // _for_each_kmer call in a single thread. Remember that _for_each_kmer
        // itself may be called by many MPI threads in parallel.
        auto cache = lru[parallel ? omp_get_thread_num() : 0];
        assert(cache != nullptr);

        map<NodeTraversal*, int> node_start;
        node_starts_in_path(forward_path, node_start);

        // now process the kmers of this sequence
        // by first getting the sequence
        string seq = path_string(forward_path);

        // but bail out if the sequence is shorter than the kmer size
        if (seq.size() < kmer_size) return;

        // and then stepping across the path, finding the kmers, and then implied node overlaps
        for (int i = 0; i <= seq.size() - kmer_size; i+=stride) {

            // get the kmer
            string forward_kmer = seq.substr(i, kmer_size);
            // record when we get a kmer match

            // We'll fill this in if needed.
            string reversed_kmer;

            // execute our callback on each kmer/node/position
            // where node == node
            int j = 0;
            while (j < kmer_size) {
                if (forward_node == node_by_path_position[i+j]) {
                    // At this position along this possible kmer, we're in the
                    // instance of the node we're interested in the kmers of.


                    // Grab the node the kmer started at
                    list<NodeTraversal>::iterator start_node = node_by_path_position[i];
                    // And the one it's going to end at
                    list<NodeTraversal>::iterator end_node = node_by_path_position[i + kmer_size - 1];
                    // Work out how far into its actual starting node this kmer started.
                    size_t start_node_offset = i - node_start[&(*start_node)];

                    if(!allow_negatives && node == nullptr) {
                        // If we do allow negatives, we'll just articulate kmers
                        // from both sides whenever they cross edges. Otherwise,
                        // we only want edge-crossing kmers once, so we should
                        // only announce them to the callback from one of their
                        // ends. We arbitrarily choose the end with the lower
                        // node ID.

                        // We only do this when we aren't getting the kmers of a
                        // specific node.

                        if(forward_node == start_node &&
                           (*start_node).node->id() > (*end_node).node->id()) {
                            // We're on the start, but it's ID is larger than the end's.
                            // Announce the kmer from the end instead.
                            ++j;
                            continue;
                        }

                        if(forward_node == end_node &&
                           (*end_node).node->id() > (*start_node).node->id()) {
                            // We're on the end, but it's ID is larger than the start's.
                            // Announce the kmer from the start instead.
                            ++j;
                            continue;
                        }

                        if((*end_node).node->id() == (*start_node).node->id() &&
                            end_node != start_node &&
                            forward_node == end_node) {

                            // If this kmer starts and ends in different
                            // instances of the same node along the path, only
                            // announce it from the one it starts in. Skip the
                            // one it ends in.
                            ++j;
                            continue;
                        }
                    }

                    // We now know we should announce this kmer from this node.

                    // Work out where this node started along the path
                    int node_position = node_start[&(*forward_node)];
                    // And how far into the node this kmer started
                    int kmer_forward_relative_start = i - node_position;
                    // Negative-offset kmers will be processed, but will be
                    // corrected to the opposite strand if negative offsets are
                    // not allowed.
                    int kmer_reversed_relative_start;
                    // Did we flip?
                    bool reversed = false;
                    if(kmer_forward_relative_start < 0 && !allow_negatives) {
                        // This kmer starts at a negative offset from this node.
                        // We need to announce it with a positive offset.

                        size_t node_length = (*forward_node).node->sequence().size();

                        if(kmer_forward_relative_start + kmer_size > node_length) {
                            // If it doesn't start or end in this node, and is
                            // just passing through, there's no way to
                            // articulate it for this node with a positive
                            // offset, so we just skip to the next character in
                            // the kmer.
                            ++j;
                            continue;
                        }

                        // We know the kmer has its end in this node.

                        // If it ends in this node, we can reverse it and get a
                        // positive offset.

                        if(reversed_kmer.empty()) {
                            // Fill in the reversed kmer the first time we need it.
                            reversed_kmer = reverse_complement(forward_kmer);
                        }

                        if(reversed_path.empty()) {
                            // Only fill in the reversed path the first time we need it.
                            for(NodeTraversal traversal : forward_path) {
                                // Operate on copies here
                                traversal.backward = !traversal.backward;
                                reversed_path.push_front(traversal);
                            }

                            // Fill in the reversed iterator too
                            reversed_node = reversed_path.begin();
                            list<NodeTraversal>::iterator i(forward_node);

                            // If forward_node is the last non-end() item, we
                            // don't want to advance reversed_node at all from
                            // the first item of the reversed path.
                            ++i;

                            while(i != forward_path.end()) {
                                // Walk i towards the end of the forward path,
                                // and reversed_node in from the corresponding
                                // end of the reversed path.
                                ++i;
                                ++reversed_node;
                            }
                        }

                        // Flip the kmer start around to something that will be positive.
                        // We don't need a -1 here.
                        kmer_reversed_relative_start = node_length - (kmer_forward_relative_start + kmer_size);

                        // We flipped.
                        reversed = true;
                    }

                    // Set up some references so we don't need to make local
                    // copies unless we actually did need to reverse the kmer.
                    string& kmer = reversed ? reversed_kmer : forward_kmer;
                    list<NodeTraversal>::iterator& instance = reversed ? reversed_node : forward_node;
                    list<NodeTraversal>& path = reversed ? reversed_path : forward_path;
                    int& kmer_relative_start = reversed ? kmer_reversed_relative_start : kmer_forward_relative_start;

                    // Make sure we aren't disobeying instructions
                    assert(!(kmer_relative_start < 0 && !allow_negatives));

                    // What key in the cache do we say that we processed? We're
                    // going to cache kmers in their forward orientation, even
                    // if they are at negative relative offsets and we aren't
                    // supposed to be using those. Because for_each_kpath only
                    // ever presents a node to this function in its forward
                    // orientation, we'll never have a situation where we should
                    // have done the cache key on the opposite strand.
                    string cache_key;
                    if (allow_dups) {
                        // Duplicate kmers starting at the same place are allowed if the paths go to different places next.
                        // This is deduplicating by node ID so we don't have to worry about instances on the path.
                        // TODO: forward_node won't be passed in as a backward traversal from for_each_kpath, right?

                        // figure out past-the-end-of-the-kmer position and node
                        list<NodeTraversal>::iterator past_end = (i+kmer_size >= node_by_path_position.size())
                            ? path.end()
                            : node_by_path_position[i + kmer_size-1];
                        int node_past_end_position = (past_end == path.end()) ? 0 : i+kmer_size - node_start[&(*past_end)];

#ifdef debug
                        cerr << "Checking for duplicates of " << (*start_node).node->id() << "." << start_node_offset
                             << "-" << forward_kmer  << "-" << (past_end == path.end() ? 0 : (*past_end).node->id()) << "."
                             << node_past_end_position << " viewed from "
                             << (*forward_node).node->id() << " offset " << kmer_forward_relative_start << endl;
#endif

                        cache_key = make_cache_key(forward_kmer, (*start_node).node->id(), start_node_offset,
                                                                 (*forward_node).node->id(), kmer_forward_relative_start,
                                                                 (past_end == path.end() ? 0 : (*past_end).node->id()),
                                                                 node_past_end_position);
                    } else {
                        // Duplicate kmers starting at the same place aren't allowed, no matter where they go after the end.
                        // This is deduplicating by node ID so we don't have to worry about instances on the path.
                        cache_key = make_cache_key(forward_kmer, (*start_node).node->id(), start_node_offset,
                                                                 (*forward_node).node->id(), kmer_forward_relative_start,
                                                                 0, 0);
                    }

                    // See if this kmer is mentioned in the cache already
                    pair<bool, bool> c = cache->retrieve(cache_key);
                    if (!c.second && (*instance).node != NULL) {
                        // TODO: how could we ever get a null node here?
                        // If not, put it in and run on it.
                        cache->put(cache_key, true);

                        lambda(kmer, instance, kmer_relative_start, path, *this);
                    } else {
#ifdef debug
                        cerr << "Skipped " << kmer << " because it was already done" << endl;
#endif
                    }

                }
                ++j;
            }
        }
    };

    auto noop = [](NodeTraversal) { };

    if(node == nullptr) {
        // Look at all the kpaths
        if (parallel) {
            for_each_kpath_parallel(kmer_size, edge_max, noop, noop, handle_path);
        } else {
            for_each_kpath(kmer_size, edge_max, noop, noop, handle_path);
        }
    } else {
        // Look only at kpaths of the specified node
        for_each_kpath_of_node(node, kmer_size, edge_max, noop, noop, handle_path);
    }

}

int VG::path_edge_count(list<NodeTraversal>& path, int32_t offset, int path_length) {
    int edge_count = 0;
    // starting from offset in the first node
    // how many edges do we cross?

    // This is the remaining path length
    int l = path_length;

    // This is the next node we are looking at.
    list<NodeTraversal>::iterator pitr = path.begin();

    // How many bases of the first node can we use?
    int available_in_first_node = (*pitr).node->sequence().size() - offset;

    if(available_in_first_node >= l) {
        // Cross no edges
        return 0;
    }

    l -= available_in_first_node;
    pitr++;
    while (l > 0) {
        // Now we can ignore node orientation
        ++edge_count;
        l -= (*pitr++).node->sequence().size();
    }
    return edge_count;
}

int VG::path_end_node_offset(list<NodeTraversal>& path, int32_t offset, int path_length) {
    // This is the remaining path length
    int l = path_length;

    // This is the next node we are looking at.
    list<NodeTraversal>::iterator pitr = path.begin();

    // How many bases of the first node can we use?
    int available_in_first_node = (*pitr).node->sequence().size() - offset;

    if(available_in_first_node >= l) {
        // Cross no edges
        return available_in_first_node - l;
    }

    l -= available_in_first_node;
    pitr++;
    while (l > 0) {
        l -= (*pitr++).node->sequence().size();
    }
    // Now back out the last node we just took.
    l += (*--pitr).node->sequence().size();

    // Measure form the far end of the last node.
    l = (*pitr).node->sequence().size() - l - 1;

    return l;
}

void VG::kmer_context(string& kmer,
                      int kmer_size,
                      int edge_max,
                      list<NodeTraversal>& path,
                      list<NodeTraversal>::iterator start_node,
                      int32_t start_offset,
                      list<NodeTraversal>::iterator& end_node,
                      int32_t& end_offset,
                      set<char>& prev_chars,
                      set<char>& next_chars,
                      set<pair<pair<int64_t, bool>, int32_t> >& prev_positions,
                      set<pair<pair<int64_t, bool>, int32_t> >& next_positions) {

    // Say we couldn't find an and node. We'll replace this when we do.
    end_node = path.end();

    // Start our iterator at our kmer's starting node instance.
    list<NodeTraversal>::iterator np = start_node;

    if (start_offset == 0) {
        // for each node connected to this one
        // what's its last character?
        // add to prev_chars

        // TODO: do we have to check if we would cross too many edges here too?

        vector<NodeTraversal> prev_nodes;
        nodes_prev(*start_node, prev_nodes);
        for (auto n : prev_nodes) {
            const string& seq = n.node->sequence();
            // We have to find the last chartacter in either orientation.
            prev_chars.insert(n.backward ? reverse_complement(seq[0]) : seq[seq.size()-1]);
            // Also note the previous position (which was always the last character in the orientation we'll be looking at it in)
            prev_positions.insert(make_pair(make_pair(n.node->id(), n.backward), n.node->sequence().size() - 1));
        }
    } else {
        // Grab and point to the previous character in this orientation of this node.
        const string& seq = (*start_node).node->sequence();
        // If we're on the reverse strand, we go one character later in the
        // string than start_offset from its end. Otherwise we go one character
        // earlier than start_offset from its beginning.
        // TODO: Add some methods to get characters at offsets in NodeTraversals
        prev_chars.insert((*start_node).backward ? reverse_complement(seq[seq.size() - start_offset]) : seq[start_offset - 1]);
        prev_positions.insert(make_pair(make_pair((*start_node).node->id(), (*start_node).backward), start_offset - 1));
    }

    // find the kmer end
    int pos = start_offset; // point at start of kmer
    bool first_in_path = true;
    // while we're not through with the path
    while (np != path.end()) {
        NodeTraversal n = *np;
        // Every time we look into a new node, we add in the amount of sequence
        // in that node. So this is the number of bases after the start of our
        // kmer up through the node we're looking at.
        // TODO: rename this to something more descriptive.
        int newpos = pos + n.node->sequence().size();

        // QUESTION:
        // Would the edge_max constraint cause us to drop the next kmer?
        // ANSWER:
        // It will when the count of edges in the implied path would be >edge_max
        // So assemble these paths and answer that question.
        // --- you can assemble any one of them, or simpler,
        // the question to ask is 1) are we losing an edge crossing on the left?
        // 2) are we gaining a new edge crossing on the right?

        if (first_in_path) {
            // Start with newpos equal to the amount of sequence in the node the
            // kmer is on, minus that consumed by the offset to the start of the
            // kmer.
            newpos = n.node->sequence().size() - pos;
            first_in_path = false;
        }
        if (newpos == kmer.size()) {
            // The kmer ended right at the end of a node

            // Save the end node instance and offset for the kmer. Remember that
            // end_offset counts rightward.
            end_node = np;
            end_offset = 0;

            // we might lose the next kmer
            // if the current path crosses the edge max number of edges
            // and the next doesn't lose an edge on the left
            // TODO: implement that

            // Look at all the nodes that might come next.
            vector<NodeTraversal> next_nodes;
            nodes_next(n, next_nodes);
            for (auto m : next_nodes) {
                // How long is this next node?
                size_t node_length = m.node->sequence().size();
                // If the next node is backward, get the rc of its last character. Else get its first.
                next_chars.insert(m.backward ? reverse_complement(m.node->sequence()[node_length - 1]) : m.node->sequence()[0]);
                // We're going to the 0 offset on this node, no matter what orientation that actually is.
                next_positions.insert(make_pair(make_pair(m.node->id(), m.backward), 0));
            }
            break;
        } else if (newpos > kmer.size()) {
            // There is at least one character in this node after the kmer

            // How long is this node?
            size_t node_length = n.node->sequence().size();

            // How far in from the left of the oriented node is the first base not in the kmer?
            // We know this won't be 0, or we'd be in the if branch above.
            int off = node_length - (newpos - kmer.size());

            // Save the end node instance and offset for the kmer. Remember that
            // end_offset counts rightward, and we need to walk in 1 more to
            // actually be on the kmer (which we accomplish by omitting the -1
            // we would usually use when reversing).
            end_node = np;
            end_offset = node_length - off;

            // Fill in the next characters and positions. Remember that off
            // points to the first character after the end of the kmer.
            next_chars.insert(n.backward ?
                reverse_complement(n.node->sequence()[node_length - off - 1]) :
                n.node->sequence()[off]);
            next_positions.insert(make_pair(make_pair(n.node->id(), n.backward), off));
            break;
        } else {
            pos = newpos;
            ++np;
        }
    }

    if(end_node == path.end()) {
        cerr << "Could not find end node for " << kmer << " at " << start_offset << " into "
             << (*start_node).node->id() << " " << (*start_node).backward << endl;
        for(auto t : path) {
            cerr << t.node->id() << " " << t.backward << endl;
        }
        assert(false);
    }
}

void VG::prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node) {
    vector<set<NodeTraversal> > prev_maxed_nodes;
    vector<set<NodeTraversal> > next_maxed_nodes;
#pragma omp parallel
    {
#pragma omp single
        {
            int threads = omp_get_num_threads();
            prev_maxed_nodes.resize(threads);
            next_maxed_nodes.resize(threads);
        }
    }
    auto prev_maxed = [this, &prev_maxed_nodes](NodeTraversal node) {
        int tid = omp_get_thread_num();
        prev_maxed_nodes[tid].insert(node);
    };
    auto next_maxed = [this, &next_maxed_nodes](NodeTraversal node) {
        int tid = omp_get_thread_num();
        next_maxed_nodes[tid].insert(node);
    };
    auto noop = [](list<NodeTraversal>::iterator node, list<NodeTraversal>& path) { };
    for_each_kpath_parallel(path_length, edge_max, prev_maxed, next_maxed, noop);

    // What nodes will we destroy because we got into them with too much complexity?
    set<Node*> to_destroy;

    set<NodeTraversal> prev;
    for (auto& p : prev_maxed_nodes) {
        for (auto node : p) {
            prev.insert(node);
        }
    }
    for (auto node : prev) {
        // This node was prev-maxed, meaning we tried to go left into it and couldn't.
        // We want to connect the end of the head node to everywhere we tried to come left from.

        // We drop any links going into the other side of this node.

        if(node.backward) {
            // Going left into it means coming to its start. True flags on the
            // edges mean connecting to the start of the other nodes.
            for (auto& e : edges_start(node.node)) {
                create_edge(e.first, head_node->id(), e.second, true);
            }
        } else {
            // Going left into it means coming to its end. True flags on the
            // edges mean connecting to the ends of the other nodes.
            for (auto& e : edges_end(node.node)) {
                create_edge(head_node->id(), e.first, false, e.second);
            }
        }

        // Remember to estroy the node. We might come to the same node in two directions.
        to_destroy.insert(node.node);
    }

    set<NodeTraversal> next;
    for (auto& n : next_maxed_nodes) {
        for (auto node : n) {
            next.insert(node);
        }
    }
    for (auto node : next) {
        // This node was next-maxed, meaning we tried to go right into it and couldn't.
        // We want to connect the start of the tail node to everywhere we tried to come right from.

        // We drop any links going into the other side of this node.

        if(node.backward) {
            // Going right into it means coming to its end. True flags on the
            // edges mean connecting to the end of the other nodes.
            for (auto& e : edges_end(node.node)) {
                create_edge(tail_node->id(), e.first, false, e.second);
            }
        } else {
            // Going right into it means coming to its start. True flags on the
            // edges mean connecting to the starts of the other nodes.
            for (auto& e : edges_start(node.node)) {
                create_edge(e.first, head_node->id(), e.second, true);
            }
        }

        // Remember to estroy the node. We might come to the same node in two directions.
        to_destroy.insert(node.node);
    }

    for(Node* n : to_destroy) {
        // Destroy all the nodes we wanted to destroy.
        destroy_node(n);
    }

    for (auto* n : head_nodes()) {
        if (n != head_node) {
            // Fix up multiple heads with a left-to-right edge
            create_edge(head_node, n);
        }
    }
    for (auto* n : tail_nodes()) {
        if (n != tail_node) {
            // Fix up multiple tails with a left-to-right edge
            create_edge(n, tail_node);
        }
    }
}

void VG::prune_short_subgraphs(size_t min_size) {
    list<VG> subgraphs;
    disjoint_subgraphs(subgraphs);
    for (auto& g : subgraphs) {
        vector<Node*> heads;
        g.head_nodes(heads);
        // calculate length
        // if < N
        if (g.total_length_of_nodes() < min_size) {
            //cerr << "removing" << endl;
            g.for_each_node([this](Node* n) {
                    // remove from this graph a node of the same id
                    //cerr << n->id() << endl;
                    this->destroy_node(n->id());
                });
        }
    }
}

/*
// todo
void VG::prune_complex_subgraphs(size_t ) {

}
*/

void VG::collect_subgraph(Node* start_node, set<Node*>& subgraph) {

    // add node to subgraph
    subgraph.insert(start_node);

    set<Node*> checked;
    set<Node*> to_check;
    to_check.insert(start_node);

    while (!to_check.empty()) {
        // for each predecessor of node
        set<Node*> curr_check = to_check;
        to_check.clear();
        for (auto* node : curr_check) {
            if (checked.count(node)) {
                continue;
            } else {
                checked.insert(node);
            }
            vector<NodeTraversal> prev;
            nodes_prev(node, prev);
            for (vector<NodeTraversal>::iterator p = prev.begin(); p != prev.end(); ++p) {
            // if it's not already been examined, collect its neighborhood
                if (!subgraph.count((*p).node)) {
                    subgraph.insert((*p).node);
                    to_check.insert((*p).node);
                }
            }
            // for each successor of node
            vector<NodeTraversal> next;
            nodes_next(node, next);
            for (vector<NodeTraversal>::iterator n = next.begin(); n != next.end(); ++n) {
                if (!subgraph.count((*n).node)) {
                    subgraph.insert((*n).node);
                    to_check.insert((*n).node);
                }
            }
        }
    }
    //cerr << "node " << start_node->id() << " subgraph size " << subgraph.size() << endl;
}

void VG::disjoint_subgraphs(list<VG>& subgraphs) {
    vector<Node*> heads;
    head_nodes(heads);
    map<Node*, set<Node*> > subgraph_by_head;
    map<Node*, set<Node*>* > subgraph_membership;
    // start at the heads, but keep in mind that we need to explore fully
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        if (subgraph_membership.find(*h) == subgraph_membership.end()) {
            set<Node*>& subgraph = subgraph_by_head[*h];
            collect_subgraph(*h, subgraph);
            for (set<Node*>::iterator n = subgraph.begin(); n != subgraph.end(); ++n) {
                subgraph_membership[*n] = &subgraph;
            }
        }
    }
    for (map<Node*, set<Node*> >::iterator g = subgraph_by_head.begin();
         g != subgraph_by_head.end(); ++ g) {
        set<Node*>& nodes = g->second;
        set<Edge*> edges;
        edges_of_nodes(nodes, edges);
        subgraphs.push_back(VG(nodes, edges));
    }
}

void VG::head_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (start_degree(n) == 0) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::head_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    return heads;
}

void VG::tail_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (end_degree(n) == 0) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::tail_nodes(void) {
    vector<Node*> tails;
    tail_nodes(tails);
    return tails;
}

void VG::wrap_with_null_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    Node* head = create_node("");
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        create_edge(head, *h);
    }

    vector<Node*> tails;
    tail_nodes(tails);
    Node* tail = create_node("");
    for (vector<Node*>::iterator t = tails.begin(); t != tails.end(); ++t) {
        create_edge(*t, tail);
    }
}


/**
Order and orient the nodes in the graph using a topological sort.

We use a birirected adaptation of Kahn's topological sort (1962), which can handle components with no heads or tails.

L ← Empty list that will contain the sorted and oriented elements
S ← Set of nodes which have been oriented, but which have not had their downstream edges examined
N ← Set of all nodes that have not yet been put into S

while N is nonempty do
    remove a node from N, orient it arbitrarily, and add it to S
        (In practice, we use "seeds": the heads and any nodes we have seen that had too many incoming edges)
    while S is non-empty do
        remove an oriented node n from S
        add n to tail of L
        for each node m with an edge e from n to m do
            remove edge e from the graph
            if m has no other edges to that side then
                orient m such that the side the edge comes to is first
                remove m from N
                insert m into S
            otherwise
                put an oriented m on the list of arbitrary places to start when S is empty
                    (This helps start at natural entry points to cycles)
    return L (a topologically sorted order and orientation)
*/
void VG::topological_sort(deque<NodeTraversal>& l) {
    //assert(is_valid());

    // using a map instead of a set ensures a stable sort across different systems
    map<int64_t, NodeTraversal> s;

    // We find the head and tails, if there are any
    vector<Node*> heads;
    head_nodes(heads);
    vector<Node*> tails;
    tail_nodes(tails);

    // We'll fill this in with the heads, so we can orient things according to
    // them first, and then arbitrarily. We ignore tails since we only orient
    // right from nodes we pick.
    // Maps from node ID to first orientation we suggested for it.
    map<int64_t, NodeTraversal> seeds;
    for(Node* head : heads) {
        seeds[head->id()] = NodeTraversal(head, false);
    }

    // We will use a std::map copy of node_by_id as our set of unvisited nodes,
    // and remove index entries from it when we visit nodes. It will be rebuilt
    // when we rebuild the indexes later. We know its order will be fixed across
    // systems.
    map<int64_t, Node*> unvisited;
    // Fill it in, we can't use the copy constructor since std::map doesn't speak vg::hash_map
    // TODO: is the vg::hash_map order fixed across systems?
    for_each_node([&](Node* node) {
        unvisited[node->id()] = node;
    });

    // How many nodes have we ordered and oriented?
    int64_t seen = 0;

    while(!unvisited.empty()) {

        // Put something in s. First go through seeds until we can find one
        // that's not already oriented.
        while(s.empty() && !seeds.empty()) {
            // Look at the first seed
            NodeTraversal first_seed = (*seeds.begin()).second;

            if(unvisited.count(first_seed.node->id())) {
                // We have an unvisited seed. Use it
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Starting from seed " << first_seed.node->id() << " orientation " << seeds.front().backward << endl;
#endif

                s[first_seed.node->id()] = first_seed;
                unvisited.erase(first_seed.node->id());
            }
            // Whether we used the seed or not, don't keep it around
            seeds.erase(seeds.begin());
        }

        if(s.empty()) {
            // If we couldn't find a seed, just grab any old node.
            // Since map order is stable across systems, we can take the first node by id and put it locally forward.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Starting from arbitrary node " << unvisited.begin()->first << " locally forward" << endl;
#endif

            s[unvisited.begin()->first] = NodeTraversal(unvisited.begin()->second, false);
            unvisited.erase(unvisited.begin()->first);
        }

        while (!s.empty()) {
            // Grab an oriented node
            NodeTraversal n = s.begin()->second;
            s.erase(n.node->id());
            l.push_back(n);
            ++seen;
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Using oriented node " << n.node->id() << " orientation " << n.backward << endl;
#endif

            // See if it has an edge from its start to the start of some node
            // where both were picked as places to break into cycles. A
            // reversing self loop on a cycle entry point is a special case of
            // this.
            vector<NodeTraversal> prev;
            nodes_prev(n, prev);
            for(NodeTraversal& prev_node : prev) {
                if(!unvisited.count(prev_node.node->id())) {
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\tHas left-side edge to cycle entry point " << prev_node.node->id()
                         << " orientation " << prev_node.backward << endl;
#endif

                    // Unindex it
                    Edge* bad_edge = get_edge(prev_node, n);
                    unindex_edge_by_node_sides(bad_edge);
                }
            }

            // All other connections and self loops are handled by looking off the right side.

            // See what all comes next, minus deleted edges.
            vector<NodeTraversal> next;
            nodes_next(n, next);

            for(NodeTraversal& next_node : next) {

#ifdef debug
#pragma omp critical (cerr)
                cerr << "\tHas edge to " << next_node.node->id() << " orientation " << next_node.backward << endl;
#endif

                // Unindex the edge connecting these nodes in this order and
                // relative orientation, so we can't traverse it with
                // nodes_next/nodes_prev.

                // Grab the edge
                Edge* edge = get_edge(n, next_node);

                // Unindex it
                unindex_edge_by_node_sides(edge);

                if(unvisited.count(next_node.node->id())) {
                    // We haven't already started here as an arbitrary cycle entry point

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tAnd node hasn't been visited yet" << endl;
#endif

                    if(node_count_prev(next_node) == 0) {

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\t\tIs last incoming edge" << endl;
#endif
                        // Keep this orientation and put it here
                        s[next_node.node->id()] = next_node;
                        // Remember that we've visited and oriented this node, so we
                        // don't need to use it as a seed.
                        unvisited.erase(next_node.node->id());

                    } else if(!seeds.count(next_node.node->id())) {
                        // We came to this node in this orientation; when we need a
                        // new node and orientation to start from (i.e. an entry
                        // point to the node's cycle), we might as well pick this
                        // one.
                        // Only take it if we don't already know of an orientation for this node.
                        seeds[next_node.node->id()] = next_node;

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\t\tSuggests seed " << next_node.node->id() << " orientation " << next_node.backward << endl;
#endif
                    }
                } else {
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tAnd node was already visited (to break a cycle)" << endl;
#endif
                }
            }

            // The caller may put us in a progress context with the denominator
            // being the number of nodes in the graph.
            update_progress(seen);
        }
    }

    // There should be no edges left in the index
    if(!edges_on_start.empty() || !edges_on_end.empty()) {
#pragma omp critical (cerr)
        {
            cerr << "Error: edges remainin after topological sort and cycle breaking" << endl;
            std::ofstream out("fail.vg");
            serialize_to_ostream(out);
            out.close();
        }
        exit(1);
    }

    // we have destroyed the graph's edge and node index to ensure its order
    // rebuild the indexes
    rebuild_indexes();
}

void VG::orient_nodes_forward(set<int64_t>& nodes_flipped) {
    // TODO: update paths in the graph when you do this!

    // Clear the flipped nodes set.
    nodes_flipped.clear();

    // First do the topological sort to order and orient
    deque<NodeTraversal> order_and_orientation;
    topological_sort(order_and_orientation);

    // These are the node IDs we've visited so far
    set<int64_t> visited;

    for(auto& traversal : order_and_orientation) {
        // Say we visited this node
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Visiting " << traversal.node->id() << endl;
#endif
        visited.insert(traversal.node->id());

        // Make sure this node is the "from" in all its edges with un-visited nodes.

        if(traversal.backward) {
            // We need to flip this node around.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Flipped node " << traversal.node->id() << endl;
#endif
            // Say we flipped it
            nodes_flipped.insert(traversal.node->id());

            // Flip the sequence
            traversal.node->set_sequence(reverse_complement(traversal.node->sequence()));

        }

        // Get all the edges
        vector<Edge*> node_edges;
        edges_of_node(traversal.node, node_edges);

        for(Edge* edge : node_edges) {
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Found edge " << edge->from() << "->" << edge->to() << endl;
#endif
            // This flag sets if we unindexed the edge. We leave it indexed unless we have to change it.
            bool unindexed = false;

            if(edge->to() == traversal.node->id()) {
                // If the other node in the edge is a node we haven't seen
                // yet in our ordering, make sure this node is the edge
                // from. It's the to now, so flip it.
                if(visited.count(edge->from()) == 0) {
                    // We need to change the edge
                    unindexed = true;
                    // Unindex it
                    unindex_edge_by_node_sides(edge);

                    // Flip the nodes
                    int64_t temp_id = edge->from();
                    edge->set_from(edge->to());
                    edge->set_to(temp_id);

                    // Move the directionality flags, but invert both.
                    bool temp_orientation = !edge->from_start();
                    edge->set_from_start(!edge->to_end());
                    edge->set_to_end(temp_orientation);
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "Reversed edge direction to " << edge->from() << "->" << edge->to() << endl;
#endif
                }
            }

            if(traversal.backward) {
                // We flipped this node around, which means we need to rewire
                // every edge. They already have the right to and from, but
                // their from_start and to_end flags could be wrong.

                // We need to change the edge
                if(!unindexed) {
                    unindexed = true;
                    unindex_edge_by_node_sides(edge);
                }

                // Flip the backwardness flag for the end that this node is on.
                if(edge->to() == traversal.node->id()) {
                    edge->set_to_end(!edge->to_end());
                }
                if(edge->from() == traversal.node->id()) {
                    edge->set_from_start(!edge->from_start());
                }
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Rewired edge " << edge->from() << "->" << edge->to() << endl;
#endif
            }

            if(unindexed) {
                // Reindex the edge by the nodes it connects to. The order of nodes doesn't
                // matter, so it's OK that unindexing doesn't remove edges from
                // the index by pair of nodes (we won't put it in twice).
                index_edge_by_node_sides(edge);
            }

            // It should always work out that the edges are from end to
            // start when we are done, but right now they might not be,
            // because the nodes at the other ends may still need to be
            // flipped themselves.

        }
    }

    // We now know there are no edges from later nodes to earlier nodes. But
    // edges that are reversing, or "from" earlier nodes and "to" later nodes
    // with both end flags set, will cause problems. So remove those.
    // This works out to just clearing any edges with from_start or to_end set.
    // We also need to clear otherwise-normal-looking self loops here.
    vector<Edge*> to_remove;
    for_each_edge([&](Edge* e) {
        if(e->from_start() || e->to_end() || e->from() == e->to()) {
            // gssw won't know how to read this. Get rid of it.
            to_remove.push_back(e);
        }
    });

    for(Edge* edge : to_remove) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Removed cycle edge " << edge->from() << "->" << edge->to() << endl;
#endif
        destroy_edge(edge);
    }

}

} // end namespace
