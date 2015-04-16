#include "vg.hpp"
#include "stream.hpp"

namespace vg {

using namespace std;

// todo embed paths

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
        extend(g);
    };

    stream::for_each(in, lambda, handle_count);

    destroy_progress();

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
            node_context(node, g);
        }
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
    paths._paths = graph.mutable_path();
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
        *new_node = node;
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.node_size()-1;
    }
}

void VG::add_edge(Edge& edge) {
    if (!has_edge(edge)) {
        Edge* new_edge = graph.add_edge(); // add it to the graph
        *new_edge = edge;
        set_edge(edge.from(), edge.to(), new_edge);
        edge_index[new_edge] = graph.edge_size()-1;
    }
}

int64_t VG::node_count(void) {
    return graph.node_size();
}

int64_t VG::edge_count(void) {
    return graph.edge_size();
}

int VG::in_degree(Node* node) {
    auto in = edges_to_from.find(node->id());
    if (in == edges_to_from.end()) {
        return 0;
    } else {
        return in->second.size();
    }
}

int VG::out_degree(Node* node) {
    auto out = edges_from_to.find(node->id());
    if (out == edges_from_to.end()) {
        return 0;
    } else {
        return out->second.size();
    }
}

void VG::edges_of_node(Node* node, vector<Edge*>& edges) {
    hash_map<int64_t, vector<int64_t> >::iterator in = edges_to_from.find(node->id());
    if (in != edges_to_from.end()) {
        vector<int64_t>::iterator e = in->second.begin();
        for (; e != in->second.end(); ++e) {
            Edge* edge = edge_by_id[make_pair(*e, node->id())];
            if (!edge) {
                cerr << "error:[VG::edges_of_node] nonexistent edge " << *e << "->" << node->id() << endl;
                exit(1);
            }
            edges.push_back(edge);
        }
    }
    hash_map<int64_t, vector<int64_t> >::iterator out = edges_from_to.find(node->id());
    if (out != edges_from_to.end()) {
        vector<int64_t>::iterator e = out->second.begin();
        for (; e != out->second.end(); ++e) {
            Edge* edge = edge_by_id[make_pair(node->id(), *e)];
            if (!edge) {
                cerr << "error:[VG::edges_of_node] nonexistent edge " << node->id() << "->" << *e << endl;
                exit(1);
            }
            edges.push_back(edge);
        }
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
        set_edge(e->from(), e->to(), e);
    }
}

void VG::clear_indexes(void) {
    node_index.clear();
    node_by_id.clear();
    edge_by_id.clear();
    edge_index.clear();
    edges_from_to.clear();
    edges_to_from.clear();
}

void VG::clear_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    node_index.clear_no_resize();
    node_by_id.clear_no_resize();
    edge_by_id.clear_no_resize();
    edge_index.clear_no_resize();
    edges_from_to.clear_no_resize();
    edges_to_from.clear_no_resize();
#else
    clear_indexes();
#endif
}

void VG::resize_indexes(void) {
    node_index.resize(graph.node_size());
    node_by_id.resize(graph.node_size());
    edge_by_id.resize(graph.edge_size());
    edge_index.resize(graph.edge_size());
    edges_from_to.resize(graph.edge_size());
    edges_to_from.resize(graph.edge_size());
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
    return edge && has_edge(edge->from(), edge->to());
}

bool VG::has_edge(Edge& edge) {
    return has_edge(edge.from(), edge.to());
}

bool VG::has_edge(int64_t from, int64_t to) {
    return edge_by_id.find(make_pair(from, to)) != edge_by_id.end();
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
        g.destroy_edge(g.get_edge((*e)->from(), (*e)->to()));
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
void VG::extend(VG& g) {
    for (int64_t i = 0; i < g.graph.node_size(); ++i) {
        Node* n = g.graph.mutable_node(i);
        if (!has_node(n)) {
            add_node(*n);
        }
    }
    for (int64_t i = 0; i < g.graph.edge_size(); ++i) {
        Edge* e = g.graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        }
    }
    paths.append(g.paths);
}

void VG::extend(Graph& graph) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (!has_node(n)) {
            add_node(*n);
        }
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
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
    // collect ids as node*'s may change
    vector<int64_t> heads_ids;
    for (vector<Node*>::iterator n = heads.begin(); n != heads.end(); ++n) {
        heads_ids.push_back((*n)->id());
    }

    // get the current tails of this graph
    vector<Node*> tails = tail_nodes();
    vector<int64_t> tails_ids;
    for (vector<Node*>::iterator n = tails.begin(); n != tails.end(); ++n) {
        tails_ids.push_back((*n)->id());
    }

    // add in the other graph
    // note that we don't use merge_union because we are ensured non-overlapping ids
    merge(g);

    /*
    cerr << "this graph size " << node_count() << " nodes " << edge_count() << " edges" << endl;
    cerr << "in append with " << heads.size() << " heads and " << tails.size() << " tails" << endl;
    */

    // now join the tails to heads
    for (vector<int64_t>::iterator t = tails_ids.begin(); t != tails_ids.end(); ++t) {
        for (vector<int64_t>::iterator h = heads_ids.begin(); h != heads_ids.end(); ++h) {
            create_edge(*t, *h);
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
            m->set_node_id(new_id[m->node_id()]);
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

    //cerr << "swapping " << node->id() << " for new id " << new_id << endl;
    int edge_n = edge_count();
    int64_t old_id = node->id();
    node->set_id(new_id);
    node_by_id.erase(old_id);

    // we check if the old node exists, and bail out if we're not doing what we expect
    assert(node_by_id.find(new_id) == node_by_id.end());

    // otherwise move to a new id
    node_by_id[new_id] = node;

    set<pair<int64_t, int64_t> > edges_to_destroy;
    set<pair<int64_t, int64_t> > edges_to_create;

    vector<int64_t>& to = edges_to(old_id);
    for (vector<int64_t>::iterator e = to.begin(); e != to.end(); ++e) {
        edges_to_create.insert(make_pair(*e, new_id));
        edges_to_destroy.insert(make_pair(*e, old_id));
    }

    vector<int64_t>& from = edges_from(old_id);
    for (vector<int64_t>::iterator e = from.begin(); e != from.end(); ++e) {
        edges_to_create.insert(make_pair(new_id, *e));
        edges_to_destroy.insert(make_pair(old_id, *e));
    }

    assert(edges_to_destroy.size() == edges_to_create.size());

    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(e->first, e->second);
    }

    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_create.begin();
         e != edges_to_create.end(); ++e) {
        create_edge(e->first, e->second);
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
                default: break;
                }
                break;
            case 4:
                switch (type) {
                case 'L': side2 = item[0]; break;
                case 'S': too_many_fields(); break;
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
            add_edge(edge);
        } else if (type == 'P') {
            paths.append_mapping(path_name, id1);
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

bool allATGC(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
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
    deque<Node*> sorted_nodes;
    topological_sort(sorted_nodes);
    deque<Node*>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.node_size() && n != sorted_nodes.end();
          ++i, ++n) {
        swap_nodes(graph.mutable_node(i), *n);
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

Edge* VG::create_edge(Node* from, Node* to) {
    return create_edge(from->id(), to->id());
}

Edge* VG::create_edge(int64_t from, int64_t to) {
    //cerr << "creating edge " << from << "->" << to << endl;
    // prevent self-linking (violates DAG/partial ordering property)
    if (to == from) return NULL;
    // ensure the edge does not already exist
    Edge* edge = get_edge(from, to);
    if (edge) return edge;
    // if not, create it
    edge = graph.add_edge();
    edge->set_from(from);
    edge->set_to(to);
    set_edge(from, to, edge);
    edge_index[edge] = graph.edge_size()-1;
    //cerr << "created edge " << edge->from() << "->" << edge->to() << endl;
    return edge;
}

Edge* VG::get_edge(int64_t from, int64_t to) {
    pair_hash_map<pair<int64_t, int64_t>, Edge*>::iterator e = edge_by_id.find(make_pair(from, to));
    if (e != edge_by_id.end()) {
        return e->second;
    } else {
        return NULL;
    }
}

void VG::set_edge(int64_t from, int64_t to, Edge* edge) {
    if (!has_edge(edge)) {
        edge_by_id[make_pair(from, to)] = edge;
        edges_from_to[from].push_back(to);
        edges_to_from[to].push_back(from);
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

vector<int64_t>& VG::edges_from(int64_t id) {
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_from_to.find(id);
    if (e == edges_from_to.end()) {
        return empty_ids;
    } else {
        return e->second;
    }
}

vector<int64_t>& VG::edges_to(int64_t id) {
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_to_from.find(id);
    if (e == edges_to_from.end()) {
        return empty_ids;
    } else {
        return e->second;
    }
}

vector<int64_t>& VG::edges_from(Node* node) {
    return edges_from(node->id());
}

vector<int64_t>& VG::edges_to(Node* node) {
    return edges_to(node->id());
}

void VG::remove_edge_fti(int64_t from, int64_t to) {
    vector<int64_t>& f = edges_from(from);
    swap_remove(f, to);
}

void VG::remove_edge_tfi(int64_t from, int64_t to) {
    vector<int64_t>& t = edges_to(to);
    swap_remove(t, from);
}

void VG::destroy_edge(int64_t from, int64_t to) {
    destroy_edge(get_edge(from, to));
}

void VG::destroy_edge(Edge* edge) {
    //cerr << "destroying edge " << edge->from() << "->" << edge->to() << endl;

    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }
    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }
    // erase from indexes

    //cerr << "erasing from indexes" << endl;

    remove_edge_fti(edge->from(), edge->to());
    remove_edge_tfi(edge->from(), edge->to());

    //assert(edges_from_to[edge->from()].find(edge->to()) == edges_from_to[edge->from()].end());
    //assert(edges_to_from[edge->to()].find(edge->from()) == edges_to_from[edge->to()].end());

    // removing the sub-indexes if they are now empty
    // we must do this to maintain a valid structure
    if (edges_from_to[edge->from()].empty()) edges_from_to.erase(edge->from());
    if (edges_to_from[edge->to()].empty()) edges_to_from.erase(edge->to());

    // erase from edges by moving to end and dropping

    // get the last edge index (lei) and this edge index (tei)
    int lei = graph.edge_size()-1;
    int tei = edge_index[edge];

    // erase this edge from the index
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
        set_edge(nlast->from(), nlast->to(), nlast);

    }

    // drop the last position, erasing the node
    graph.mutable_edge()->RemoveLast();

    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }

}

Node* VG::get_node(int64_t id) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(id);
    if (n != node_by_id.end()) {
        return n->second;
    } else {
        // again... should this throw an error?
        return NULL;
    }
}

Node* VG::create_node(string seq) {
    // create the node
    Node* node = graph.add_node();
    node->set_sequence(seq);
    node->set_id(current_id++);
    //cerr << "creating node " << node->id() << endl;
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

// a graph composed of this node and its edges
void VG::node_context(Node* node, VG& g) {
    // add the node
    g.add_node(*node);
    // and its edges
    vector<int64_t>& to = edges_to(node->id());
    for (vector<int64_t>::iterator e = to.begin(); e != to.end(); ++e) {
        g.add_edge(*get_edge(*e, node->id()));
    }
    vector<int64_t>& from = edges_from(node->id());
    for (vector<int64_t>::iterator e = from.begin(); e != from.end(); ++e) {
        g.add_edge(*get_edge(node->id(), *e));
    }
    // and its path members
    auto& node_mappings = paths.get_node_mapping(node);
    for (auto& i : node_mappings) {
        g.paths.append_mapping(i.first->name(), *i.second);
    }
}

void VG::destroy_node(int64_t id) {
    destroy_node(get_node(id));
}

void VG::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    //cerr << "destroying node " << node->id() << endl;
    // noop on NULL/nonexistent node
    if (!has_node(node)) { return; }
    // remove edges associated with node
    set<pair<int64_t, int64_t> > edges_to_destroy;
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_from_to.find(node->id());
    if (e != edges_from_to.end()) {
        for (vector<int64_t>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(make_pair(node->id(), *f));
        }
    }
    e = edges_to_from.find(node->id());
    if (e != edges_to_from.end()) {
        for (vector<int64_t>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(make_pair(*f, node->id()));
        }
    }
    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(e->first, e->second);
    }
    // assert cleanup
    edges_to_from.erase(node->id());
    edges_from_to.erase(node->id());

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
    vector<int64_t>& to = edges_to(node);
    vector<int64_t>& from = edges_from(node);
    // for edge to
    set<pair<int64_t, int64_t> > edges_to_create;
    for (vector<int64_t>::iterator f = to.begin(); f != to.end(); ++f) {
        for (vector<int64_t>::iterator t = from.begin(); t != from.end(); ++t) {
            // connect
            edges_to_create.insert(make_pair(*f, *t));
        }
    }
    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_create.begin();
         e != edges_to_create.end(); ++e) {
        create_edge(e->first, e->second);
    }
    destroy_node(node);
}

void VG::remove_orphan_edges(void) {
    set<pair<int64_t, int64_t> > edges;
    for_each_edge([this,&edges](Edge* edge) {
            if (!has_node(edge->from())
                || !has_node(edge->to())) {
                edges.insert(make_pair(edge->from(), edge->to()));
            }
        });
    for (auto edge : edges) {
        destroy_edge(edge.first, edge.second);
    }
}

void VG::keep_paths(set<string>& path_names, set<string>& kept_names) {
    // edges have implicit path
    // now... at least ...
    // maybe they shouldn't
    vector<Node*> path;
    vector<Node*> nodes_to_remove;
    for_each_node([this, &kept_names, &path_names, &path, &nodes_to_remove](Node* node) {
            // use set intersection
            bool to_keep = false;
            for (auto& s : paths.of_node(node->id())) {
                if (path_names.count(s)) {
                    kept_names.insert(s);
                    to_keep = true;
                    break;
                }
            }
            if (to_keep) {
                path.push_back(node);
            } else {
                nodes_to_remove.push_back(node);
            }
        });
    set<pair<int64_t, int64_t> > edges_to_keep;
    if (path.size()) {
        Node* prev = path.front();
        for (auto node : path) {
            if (node != prev) {
                edges_to_keep.insert(make_pair(prev->id(), node->id()));
                prev = node;
            }
        }
    }
    set<pair<int64_t, int64_t> > edges_to_destroy;
    for_each_edge([this, &edges_to_keep, &edges_to_destroy](Edge* edge) {
            auto ep = make_pair(edge->from(), edge->to());
            if (!edges_to_keep.count(ep)) {
                edges_to_destroy.insert(ep);
            }
        });
    for (auto edge : edges_to_destroy) {
        destroy_edge(edge.first, edge.second);
    }
    for (auto node : nodes_to_remove) {
        destroy_node(node);
    }
    set<string> names;
    for (auto& s : path_names) {
        names.insert(s);
    }
    paths.keep_paths(names);
}

void VG::keep_path(string& path_name) {
    set<string> s,k; s.insert(path_name);
    keep_paths(s, k);
}

// utilities
void VG::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    //cerr << "dividing node " << node->id() << endl;

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
    set<pair<int64_t, int64_t> > edges_to_create;

    // replace node connections to prev (left)
    e = edges_to_from.find(node->id());
    if (e != edges_to_from.end()) {
        for (vector<int64_t>::const_iterator p = e->second.begin();
             p != e->second.end(); ++p) {
            edges_to_create.insert(make_pair(*p, left->id()));
        }
    }

    // make our right node
    right = create_node(node->sequence().substr(pos,node->sequence().size()-1));

    // replace node connections to next (right)
    e = edges_from_to.find(node->id());
    if (e != edges_from_to.end()) {
        for (vector<int64_t>::const_iterator n = e->second.begin();
             n != e->second.end(); ++n) {
            edges_to_create.insert(make_pair(right->id(), *n));
        }
    }

    // create the edges here as otherwise we will invalidate the iterators
    for (set<pair<int64_t, int64_t> >::iterator c = edges_to_create.begin();
         c != edges_to_create.end(); ++c) {
        create_edge(c->first, c->second);
    }

    // connect left to right
    create_edge(left, right);

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

void VG::nodes_prev(Node* node, vector<Node*>& nodes) {
    vector<int64_t>& from = edges_to(node);
    for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
        nodes.push_back(node_by_id[*f]);
    }
}

void VG::nodes_next(Node* node, vector<Node*>& nodes) {
    vector<int64_t>& to = edges_from(node);
    for (vector<int64_t>::iterator t = to.begin(); t != to.end(); ++t) {
        nodes.push_back(node_by_id[*t]);
    }
}

void VG::prev_kpaths_from_node(Node* node, int length, int edge_max,
                               list<Node*> postfix, set<list<Node*> >& paths) {
    if (length == 0 || edge_max == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    postfix.push_front(node);
    vector<Node*> prev_nodes;
    nodes_prev(node, prev_nodes);
    if (prev_nodes.empty()) {
        list<Node*> new_path = postfix;
        paths.insert(new_path);
    } // implicit else
    for (vector<Node*>::iterator p = prev_nodes.begin(); p != prev_nodes.end(); ++p) {
        if ((*p)->sequence().size() < length) {
            prev_kpaths_from_node(*p, length - (*p)->sequence().size(), edge_max - 1, postfix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = postfix;
            new_path.push_front(*p);
            paths.insert(new_path);
        }
    }
}

void VG::next_kpaths_from_node(Node* node, int length, int edge_max,
                               list<Node*> prefix, set<list<Node*> >& paths) {
    if (length == 0 || edge_max == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    prefix.push_back(node);
    vector<Node*> next_nodes;
    nodes_next(node, next_nodes);
    if (next_nodes.empty()) {
        list<Node*> new_path = prefix;
        paths.insert(new_path);
    } // implicit else
    for (vector<Node*>::iterator n = next_nodes.begin(); n != next_nodes.end(); ++n) {
        if ((*n)->sequence().size() < length) {
            next_kpaths_from_node(*n, length - (*n)->sequence().size(), edge_max - 1, prefix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = prefix;
            new_path.push_back(*n);
            paths.insert(new_path);
        }
    }
}

// iterate over the kpaths in the graph, doing something

void VG::for_each_kpath(int k, int edge_max,
                        function<void(Node*,list<Node*>&)> lambda) {
    auto by_node = [k, edge_max, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, lambda);
    };
    for_each_node(by_node);
}

void VG::for_each_kpath(int k, int edge_max,
                        function<void(Node*,Path&)> lambda) {
    auto by_node = [k, edge_max, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, lambda);
    };
    for_each_node(by_node);
}

// parallel versions of above
// this isn't by default because the lambda may have side effects
// that need to be guarded explicitly

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(Node*,list<Node*>&)> lambda) {
    auto by_node = [k, edge_max, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, lambda);
    };
    for_each_node_parallel(by_node);
}

void VG::for_each_kpath_parallel(int k, int edge_max,
                                 function<void(Node*,Path&)> lambda) {
    auto by_node = [k, edge_max, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, edge_max, lambda);
    };
    for_each_node_parallel(by_node);
}

// per-node kpaths

void VG::for_each_kpath_of_node(Node* n, int k, int edge_max,
                                function<void(Node*,Path&)> lambda) {
    auto apply_to_path = [&lambda, this](Node* n, list<Node*>& p) {
        Path path = create_path(p);
        lambda(n, path);
    };
    for_each_kpath_of_node(n, k, edge_max, apply_to_path);
}

void VG::for_each_kpath_of_node(Node* node, int k, int edge_max,
                                function<void(Node*,list<Node*>&)> lambda) {
    // get left, then right
    set<list<Node*> > prev_paths;
    set<list<Node*> > next_paths;
    list<Node*> empty_list;
    prev_kpaths_from_node(node, k, edge_max, empty_list, prev_paths);
    next_kpaths_from_node(node, k, edge_max, empty_list, next_paths);
    // now take the cross and give to the callback
    for (set<list<Node*> >::iterator p = prev_paths.begin(); p != prev_paths.end(); ++p) {
        for (set<list<Node*> >::iterator n = next_paths.begin(); n != next_paths.end(); ++n) {
            list<Node*> path = *p;
            list<Node*>::const_iterator m = n->begin(); ++m; // skips current node, which is included in *p
            while (m != n->end()) {
                path.push_back(*m);
                ++m;
            }
            lambda(node, path);
        }
    }
}

void VG::kpaths_of_node(Node* node, set<list<Node*> >& paths,
                        int length, int edge_max) {
    auto collect_path = [&paths](Node* n, list<Node*>& path) {
        paths.insert(path);
    };
    for_each_kpath_of_node(node, length, edge_max, collect_path);
}

void VG::kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, int edge_max) {
    set<list<Node*> > unique_paths;
    kpaths_of_node(node, unique_paths, length, edge_max);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

// aggregators, when a callback won't work

void VG::kpaths(set<list<Node*> >& paths, int length, int edge_max) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        kpaths_of_node(node, paths, length, edge_max);
    }
}

void VG::kpaths(vector<Path>& paths, int length, int edge_max) {
    set<list<Node*> > unique_paths;
    kpaths(unique_paths, length, edge_max);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

// path utilities
// these are in this class because attributes of the path (such as its sequence) are a property of the graph

Path VG::create_path(const list<Node*>& nodes) {
    Path path;
    for (list<Node*>::const_iterator n = nodes.begin(); n != nodes.end(); ++n) {
        Mapping* mapping = path.add_mapping();
        mapping->set_node_id((*n)->id());
    }
    return path;
}

string VG::path_string(const list<Node*>& nodes) {
    string seq;
    for (list<Node*>::const_iterator n = nodes.begin(); n != nodes.end(); ++n) {
        seq.append((*n)->sequence());
    }
    return seq;
}

string VG::path_string(Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        Mapping* m = path.mutable_mapping(i);
        Node* n = node_by_id[m->node_id()];
        seq.append(n->sequence());
    }
    return seq;
}

void VG::expand_path(const list<Node*>& path, vector<Node*>& expanded) {
    for (list<Node*>::const_iterator n = path.begin(); n != path.end(); ++n) {
        Node* node = *n;
        int s = node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(node);
        }
    }
}

void VG::include(Path& path) {
    for (int i = 0; i < path.mapping_size(); ++i) {
        Mapping* m = path.mutable_mapping(i);
        Node* n = get_node(m->node_id());
        int f = 0;
        int t = 0;
        for (int j = 0; j < m->edit_size(); ++j) {
            const Edit& edit = m->edit(j);
            if (edit.has_from_length() && !edit.has_to_length()) {
                f += edit.from_length();
                t += edit.from_length();
            } else if (edit.has_from_length() && edit.has_to_length()) {
                // cut at f, and f + from_length
                Node* l=NULL;
                Node* r=NULL;
                Node* m=NULL;
                Node* c=NULL;
                if (edit.from_length()) {
                    divide_node(n, f, l, m);
                    divide_node(m, edit.from_length(), m, r);
                    // there is a quirk (for sanity, efficiency later)
                    // if we "delete" over several nodes, we should join one path for the deletion
                    if (edit.to_length() == 0) {
                        // deletion
                        assert(edit.from_length());
                        create_edge(l, r);
                    } else {
                        // swap/ SNP
                        c = create_node(edit.sequence());
                        create_edge(l, c);
                        create_edge(c, r);
                    }
                } else {
                    divide_node(n, f, l, r);
                    // insertion
                    assert(edit.to_length());
                    c = create_node(edit.sequence());
                    create_edge(l, c);
                    create_edge(c, r);
                }
            } // do nothing for soft clips, where we have to_length and unset from_length
        }
    }
}

void VG::node_starts_in_path(const list<Node*>& path,
                             map<Node*, int>& node_start) {
    int i = 0;
    for (list<Node*>::const_iterator n = path.begin(); n != path.end(); ++n) {
        node_start[*n] = i;
        int l = (*n)->sequence().size();
        i += l;
    }
}

void VG::kpaths_of_node(int64_t node_id, vector<Path>& paths, int length, int edge_max) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(node_id);
    if (n != node_by_id.end()) {
        Node* node = n->second;
        kpaths_of_node(node, paths, length, edge_max);
    }
}

string VG::path_sequence(const Path& path) {
    string sequence;
    for (int i = 0; i < path.mapping_size(); ++i) {
        sequence.append(node_by_id[path.mapping(i).node_id()]->sequence());
    }
    return sequence;
}

string VG::random_read(int length, mt19937& rng, int64_t min_id, int64_t max_id, bool either_strand) {
    uniform_int_distribution<int64_t> int64_dist(min_id, max_id);
    int64_t id = int64_dist(rng);
    Node* node = get_node(id);
    int32_t start_pos = 0;
    if (node->sequence().size() > 1) {
        uniform_int_distribution<uint32_t> uint32_dist(0,node->sequence().size()-1);
        start_pos = uint32_dist(rng);
    }
    string read = node->sequence().substr(start_pos);
    while (read.size() < length) {
        // pick a random downstream node
        vector<Node*> next_nodes;
        nodes_next(node, next_nodes);
        if (next_nodes.empty()) break;
        uniform_int_distribution<int> next_dist(0, next_nodes.size()-1);
        node = next_nodes.at(next_dist(rng));
        read.append(node->sequence());
    }
    read = read.substr(0, length);
    uniform_int_distribution<int> binary_dist(0, 1);
    if (either_strand && binary_dist(rng) == 1) {
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

        if (edges_from_to.find(f) == edges_from_to.end()) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_from_to for node " << f << endl;
            return false;
        }

        if (edges_to_from.find(t) == edges_to_from.end()) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_to_from for node " << t << endl;
            return false;
        }

    }

    //cerr << "there are " << edges_from_to.size() << " items in the edges_from_to index" << endl;
    for (hash_map<int64_t, vector<int64_t> >::iterator f = edges_from_to.begin();
         f != edges_from_to.end(); ++f) {
        vector<int64_t>& to = f->second;
        for (vector<int64_t>::iterator t = to.begin(); t != to.end(); ++t) {
            Edge* e = get_edge(f->first, *t);
            //cerr << "edges_from_to " << e << " " << e->from() << "->" << e->to() << endl;
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if (*t != e->to() || f->first != e->from()) {
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " stored in to_from index under " << f->first << "->" << *t << endl;
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

    //cerr << "there are " << edges_to_from.size() << " items in the edges_to_from index" << endl;
    for (hash_map<int64_t, vector<int64_t> >::iterator t = edges_to_from.begin();
         t != edges_to_from.end(); ++t) {
        vector<int64_t>& from = t->second;
        for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
            Edge* e = get_edge(*f, t->first);
            //cerr << "edges_to_from " << e << " " << e->from() << "->" << e->to() << endl;
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if (t->first != e->to() || *f != e->from()) {
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " stored in to_from index under " << *f << "->" << t->first << endl;
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

    if (head_nodes().empty()) {
        cerr << "graph invalid: no head nodes" << endl;
        return false;
    }

    if (tail_nodes().empty()) {
        cerr << "graph invalid: no tail nodes" << endl;
        return false;
    }

    //cerr << "all is well" << endl;

    return true;
}

void VG::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        auto node_paths = paths.of_node(n->id());
        if (node_paths.empty()) {
            out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\",fontcolor=red];" << endl;          
        } else {
            out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;
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
        if (both_paths.empty()) {
            out << "    " << e->from() << " -> " << e->to() << "[color=red];" << endl;
        } else {
            out << "    " << e->from() << " -> " << e->to() << ";" << endl;
        }
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
        for (auto& name : paths.of_node(n->id())) {
            s << "P" << "\t" << n->id() << "\t" << name << "\n";
        }
        sorted_output[n->id()].push_back(s.str());
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        stringstream s;
        s << "L" << "\t" << e->from() << "\t" << "-" << "\t" << e->to() << "\t" << "+" << "\t" << "0M" << endl;
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

void VG::connect_node_to_nodes(Node* node, vector<Node*>& nodes) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        create_edge(node->id(), (*n)->id());
    }
}

// join all subgraphs together to a "null" head node
Node* VG::join_heads(void) {
    current_id = max_node_id()+1;
    Node* root = create_node("N");
    vector<Node*> heads;
    head_nodes(heads);
    connect_node_to_nodes(root, heads);
    return root;
}

Alignment& VG::align(Alignment& alignment) {

    // to be completely aligned, the graph's head nodes need to be fully-connected to a common root
    Node* root = join_heads();
    sort();

    gssw_aligner = new GSSWAligner(graph);
    gssw_aligner->align(alignment);
    delete gssw_aligner;
    gssw_aligner = NULL;

    destroy_node(root);

    return alignment;
}

Alignment VG::align(string& sequence) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment);
}

void VG::for_each_kmer_parallel(int kmer_size,
                                int edge_max,
                                function<void(string&, Node*, int)> lambda,
                                int stride) {
    _for_each_kmer(kmer_size, edge_max, lambda, true, stride);
}

void VG::for_each_kmer(int kmer_size,
                       int edge_max,
                       function<void(string&, Node*, int)> lambda,
                       int stride) {
    _for_each_kmer(kmer_size, edge_max, lambda, false, stride);
}

void VG::_for_each_kmer(int kmer_size,
                        int edge_max,
                        function<void(string&, Node*, int)> lambda,
                        bool parallel,
                        int stride) {

    // use an LRU cache to clean up duplicates over the last 1mb
    // use one per thread so as to avoid contention
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
    auto make_cache_key = [](string& kmer, Node* node, int start) -> string {
        string cache_key = kmer;
        cache_key.resize(kmer.size() + sizeof(Node*) + sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size(), &node, sizeof(Node*));
        memcpy((char*)cache_key.c_str()+kmer.size()+sizeof(Node*), &start, sizeof(Node*));
        return cache_key;
    };

    auto handle_path = [this,
                        lambda,
                        kmer_size,
                        stride,
                        &lru,
                        &make_cache_key](Node* node, list<Node*>& path) {

        // expand the path into a vector :: 1,1,1,2,2,2,2,3,3 ... etc.
        // this makes it much easier to quickly get all the node matches of each kmer
        vector<Node*> node_by_path_position;
        expand_path(path, node_by_path_position);

        auto cache = lru[omp_get_thread_num()];

        map<Node*, int> node_start;
        node_starts_in_path(path, node_start);
        
        // now process the kmers of this sequence
        // by first getting the sequence
        string seq = path_string(path);

        // and then stepping across the path, finding the kmers, and then implied node overlaps
        for (int i = 0; i <= seq.size() - kmer_size; i+=stride) {

            // get the kmer
            string kmer = seq.substr(i, kmer_size);
            // record when we get a kmer match

            // execute our callback on each kmer/node/position
            // where node == node
            int j = 0;
            while (j < kmer_size) {
                if (node == node_by_path_position[i+j]) {
                    int node_position = node_start[node];
                    int kmer_relative_start = i - node_position;
                    string cache_key = make_cache_key(kmer, node, kmer_relative_start);
                    pair<bool, bool> c = cache->retrieve(cache_key);
                    if (!c.second) {
                        cache->put(cache_key, true);
                        lambda(kmer, node, kmer_relative_start);
                    }
                }
                ++j;
            }
        }
    };

    if (parallel) {
        for_each_kpath_parallel(kmer_size, edge_max, handle_path);
    } else {
        for_each_kpath(kmer_size, edge_max, handle_path);
    }

}

void VG::collect_subgraph(Node* node, set<Node*>& subgraph) {

    // add node to subgraph
    subgraph.insert(node);

    // for each predecessor of node
    vector<Node*> prev;
    nodes_prev(node, prev);
    for (vector<Node*>::iterator p = prev.begin(); p != prev.end(); ++p) {
        // if it's not already been examined, collect its neighborhood
        if (!subgraph.count(*p)) {
            collect_subgraph(*p, subgraph);
        }
    }

    // for each successor of node
    vector<Node*> next;
    nodes_next(node, next);
    for (vector<Node*>::iterator n = next.begin(); n != next.end(); ++n) {
        if (!subgraph.count(*n)) {
            collect_subgraph(*n, subgraph);
        }
    }

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
        if (edges_to_from.find(n->id()) == edges_to_from.end()) {
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
        if (edges_from_to.find(n->id()) == edges_from_to.end()) {
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


    /*
Kahn's topological sort (1962)

L ← Empty list that will contain the sorted elements
S ← Set of all nodes with no incoming edges
while S is non-empty do
    remove a node n from S
    add n to tail of L
    for each node m with an edge e from n to m do
        remove edge e from the graph
        if m has no other incoming edges then
            insert m into S
if graph has edges then
    return error (graph has at least one cycle)
else 
    return L (a topologically sorted order)
    */

void VG::topological_sort(deque<Node*>& l) {
    //assert(is_valid());

    // using a map instead of a set ensures a stable sort across different systems
    map<int64_t, Node*> s;
    vector<Node*> heads;
    head_nodes(heads);
    for (vector<Node*>::iterator n = heads.begin(); n != heads.end(); ++n) {
        s[(*n)->id()] = *n;
    }

    // check that we have heads of the graph
    if (heads.empty() && graph.node_size() > 0) {
        cerr << "error:[VG::topological_sort] No heads of graph given, but graph not empty. "
             << "In-memory indexes of nodes and edges may be out of sync." << endl;
        exit(1);
    }
    int64_t seen = heads.size();

    while (!s.empty()) {
        Node* n = s.begin()->second;
        s.erase(n->id());
        l.push_back(n);
        vector<int64_t>& from = edges_from(n);
        vector<int64_t> to_erase;
        for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
            Node* m = node_by_id[*f];
            ++seen;
            to_erase.push_back(m->id());
            remove_edge_tfi(n->id(), m->id());
            if (edges_to(m->id()).empty()) {
                s[m->id()] = m;
            }
        }
        update_progress(seen);
        for (vector<int64_t>::iterator t = to_erase.begin(); t != to_erase.end(); ++t) {
            remove_edge_fti(n->id(), *t);
        }
    }
    // if we have a cycle, signal an error, as we are not guaranteed an order
    for (hash_map<int64_t, vector<int64_t> >::iterator f = edges_from_to.begin();
         f != edges_from_to.end(); ++f) {
        if (!f->second.empty()) {
            cerr << "error:[VG::topological_sort] graph has a cycle from " << f->first
                 << " to " << f->second.front() << endl
                 << "thread " << omp_get_thread_num() << endl;
#pragma omp critical
            {
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
        }
    }
    for (hash_map<int64_t, vector<int64_t> >::iterator t = edges_to_from.begin();
         t != edges_to_from.end(); ++t) {
        if (!t->second.empty()) {
            cerr << "error:[VG::topological_sort] graph has a cycle to " << t->first
                 << " to " << t->second.front() << endl
                 << "thread " << omp_get_thread_num() << endl;
#pragma omp critical
            {
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
        }
    }
    // we have destroyed the graph's index to ensure its order
    // rebuild the indexes
    rebuild_indexes();
}

} // end namespace
