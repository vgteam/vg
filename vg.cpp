#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;


//VG::VG(void) { };
// construct from protobufs
VG::VG(istream& in) {
    init();
    graph.ParseFromIstream(&in);
    build_indexes();
}

VG::~VG(void) {
    destroy_alignable_graph();
}

VG::VG(void) {
    init();
}

void VG::init(void) {
    gssw_aligner = NULL;
}

VG::VG(set<Node*>& nodes, set<Edge*>& edges) {
    init();
    add_nodes(nodes);
    add_edges(edges);
    topologically_sort_graph();
}

// check for conflict (duplicate nodes and edges) occurs within add_* functions

void VG::add_nodes(set<Node*>& nodes) {
    for (set<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        add_node(**n);
    }
}

void VG::add_edges(set<Edge*>& edges) {
    for (set<Edge*>::iterator e = edges.begin(); e != edges.end(); ++e) {
        add_edge(**e);
    }
}

void VG::add_nodes(vector<Node>& nodes) {
    for (vector<Node>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        add_node(*n);
    }
}

void VG::add_edges(vector<Edge>& edges) {
    for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
        add_edge(*e);
    }
}


void VG::add_node(Node& node) {
    if (!has_node(node)) {
        Node* new_node = graph.add_node(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.node_size()-1;
    }
}

void VG::add_edge(Edge& edge) {
    if (!has_edge(edge)) {
        Edge* new_edge = graph.add_edge(); // add it to the graph
        new_edge->set_from(edge.from());
        new_edge->set_to(edge.to());
        edge_from_to[edge.from()][edge.to()] = new_edge;
        edge_to_from[edge.to()][edge.from()] = new_edge;
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
    map<int64_t, map<int64_t, Edge*> >::iterator in = edge_to_from.find(node->id());
    if (in == edge_to_from.end()) {
        return 0;
    } else {
        return in->second.size();
    }
}

int VG::out_degree(Node* node) {
    map<int64_t, map<int64_t, Edge*> >::iterator out = edge_from_to.find(node->id());
    if (out == edge_from_to.end()) {
        return 0;
    } else {
        return out->second.size();
    }
}

void VG::edges_of_node(Node* node, vector<Edge*>& edges) {
    map<int64_t, map<int64_t, Edge*> >::iterator in = edge_to_from.find(node->id());
    if (in != edge_to_from.end()) {
        map<int64_t, Edge*>::iterator e = in->second.begin();
        for (; e != in->second.end(); ++e) {
            edges.push_back(e->second);
        }
    }
    map<int64_t, map<int64_t, Edge*> >::iterator out = edge_from_to.find(node->id());
    if (out != edge_from_to.end()) {
        map<int64_t, Edge*>::iterator e = out->second.begin();
        for (; e != out->second.end(); ++e) {
            edges.push_back(e->second);
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
        edge_from_to[e->from()][e->to()] = e;
        edge_to_from[e->to()][e->from()] = e;
    }
}

void VG::clear_indexes(void) {
    node_index.clear();
    node_by_id.clear();
    edge_index.clear();
    edge_from_to.clear();
    edge_to_from.clear();
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
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(from);
    return e != edge_from_to.end() && e->second.find(to) != e->second.end();
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

void VG::merge(VG& g) {
    // remove duplicates, then merge
    remove_duplicated_in(g);
    if (g.graph.node_size() > 0) {
        merge(g.graph);
    }
}

// this merges without any validity checks
// this could be rather expensive if the graphs to merge are largely overlapping
void VG::merge(Graph& g) {
    graph.mutable_node()->MergeFrom(g.node());
    graph.mutable_edge()->MergeFrom(g.edge());
    clear_indexes();
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        edge_index[e] = i;
        edge_from_to[e->from()][e->to()] = e;
        edge_to_from[e->to()][e->from()] = e;
    }
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
VG::VG(vcf::VariantCallFile& variantCallFile, FastaReference& reference) {
    init();
    for (vector<string>::iterator r = reference.index->sequenceNames.begin();
         r != reference.index->sequenceNames.end(); ++r) {

        string& seqName = *r;
        map<long, Node*> reference_path;
        //map<long, set<Node*> > nodes; // for maintaining a reference-sorted graph
        string seq = reference.getSequence(seqName);

        Node* ref_node = create_node(seq);
        reference_path[0] = ref_node;

        // track the last nodes so that we can connect everything completely when variants occur in succession
        map<long, set<Node*> > nodes_by_end_position;

        variantCallFile.setRegion(seqName);
        vcf::Variant var(variantCallFile);

        while (variantCallFile.getNextVariant(var)) {

            int current_pos = (long int) var.position - 1;
            // decompose the alt
            bool flat_input_vcf = false; // hack
            map<string, vector<vcf::VariantAllele> > alternates = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());
            for (map<string, vector<vcf::VariantAllele> >::iterator va = alternates.begin(); va !=alternates.end(); ++va) {
                vector<vcf::VariantAllele>& alleles = va->second;

                for (vector<vcf::VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    vcf::VariantAllele& allele = *a;

                    // reference alleles are provided naturally by the reference itself
                    if (allele.ref == allele.alt) {
                        continue;
                    }

                    long allele_start_pos = allele.position - 1;  // 0/1 based conversion... thanks vcflib!
                    long allele_end_pos = allele_start_pos + allele.ref.size();

                    if (allele_start_pos == 0) {
                        Node* root = create_node(""); // ensures that we can handle variation at first position (important when aligning)
                        reference_path[-1] = root;
                    }

                    Node* left_ref_node = NULL;
                    Node* middle_ref_node = NULL;
                    Node* right_ref_node = NULL;

                    // divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {
                    divide_path(reference_path,
                                allele_start_pos,
                                left_ref_node,
                                right_ref_node);

                    //cerr << "nodes: left: " << left_ref_node->id() << " right: " << right_ref_node->id() << endl;

                    // if the ref portion of the allele is not empty, then we need to make another cut
                    if (!allele.ref.empty()) {
                        divide_path(reference_path,
                                    allele_end_pos,
                                    middle_ref_node,
                                    right_ref_node);
                    }

                    Node* alt_node;
                    // create a new alt node and connect the pieces from before
                    if (!allele.alt.empty() && !allele.ref.empty()) {

                        alt_node = create_node(allele.alt);
                        //ref_map.add_node(alt_node, allele_start_pos, );
                        create_edge(left_ref_node, alt_node);
                        create_edge(alt_node, right_ref_node);

                        // XXXXXXXX middle is borked
                        // why do we have to force this edge back in?
                        // ... because it's not in the ref map?? (??)
                        create_edge(left_ref_node, middle_ref_node);

                        nodes_by_end_position[allele_end_pos].insert(alt_node);
                        nodes_by_end_position[allele_end_pos].insert(middle_ref_node);

                    } else if (!allele.alt.empty()) { // insertion

                        alt_node = create_node(allele.alt);
                        create_edge(left_ref_node, alt_node);
                        create_edge(alt_node, right_ref_node);
                        nodes_by_end_position[allele_end_pos].insert(alt_node);
                        nodes_by_end_position[allele_end_pos].insert(left_ref_node);

                    } else {// otherwise, we have a deletion

                        create_edge(left_ref_node, right_ref_node);
                        nodes_by_end_position[allele_end_pos].insert(left_ref_node);

                    }

                    if (allele_end_pos == seq.size()) {
                        // ensures that we can handle variation at first position (important when aligning)
                        Node* end = create_node("");
                        reference_path[allele_end_pos] = end;
                        create_edge(alt_node, end);
                        create_edge(middle_ref_node, end);
                    }

                    // if there are previous nodes, connect them
                    map<long, set<Node*> >::iterator ep = nodes_by_end_position.find(allele_start_pos);
                    if (ep != nodes_by_end_position.end()) {
                        set<Node*>& previous_nodes = ep->second;
                        for (set<Node*>::iterator n = previous_nodes.begin(); n != previous_nodes.end(); ++n) {
                            if (node_index.find(*n) != node_index.end()) {
                                if (middle_ref_node) {
                                    create_edge(*n, middle_ref_node);
                                }
                                create_edge(*n, alt_node);
                            }
                        }
                    }
                    // clean up previous
                    while (nodes_by_end_position.begin()->first < allele_start_pos) {
                        nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
                    }

                    /*
                    if (!is_valid()) {
                        cerr << "graph is invalid after variant" << endl
                             << var << endl;
                        exit(1);
                    }
                    */
                }
            }
        }
    }

    topologically_sort_graph();

}

void VG::topologically_sort_graph(void) {
    list<Node*> sorted_nodes;
    topological_sort(sorted_nodes);
    list<Node*>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.node_size() && n != sorted_nodes.end();
         ++i, ++n) {
        swap_nodes(graph.mutable_node(i), *n);
    }
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
    // prevent self-linking (violates DAG/partial ordering property)
    if (to == from) return NULL;
    // ensure the edge does not already exist
    Edge* edge = edge_from_to[from][to];
    if (edge) return edge;
    // if not, create it
    edge = graph.add_edge();
    edge->set_from(from);
    edge->set_to(to);
    edge_from_to[from][to] = edge;
    edge_to_from[to][from] = edge;
    edge_index[edge] = graph.edge_size()-1;
    return edge;
}

Edge* VG::get_edge(int64_t from, int64_t to) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(from);
    if (e != edge_from_to.end() && e->second.find(to) != e->second.end()) {
        return e->second[to];
    } else {
        // or error?
        return NULL;
    }
}

void VG::destroy_edge(int64_t from, int64_t to) {
    destroy_edge(get_edge(from, to));
}

void VG::destroy_edge(Edge* edge) {
    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }
    //if (!is_valid()) cerr << "graph ain't valid" << endl;
    // erase from indexes
    edge_from_to[edge->from()].erase(edge->to());
    edge_to_from[edge->to()].erase(edge->from());

    // erase from edges by moving to end and dropping
    int lei = graph.edge_size()-1;
    int tei = edge_index[edge];
    Edge* last = graph.mutable_edge(lei);
    edge_index.erase(last);
    edge_index.erase(edge);

    // swap
    graph.mutable_edge()->SwapElements(tei, lei);

    // point to new position
    Edge* nlast = graph.mutable_edge(tei);

    // insert the new edge index position
    edge_index[nlast] = tei;

    // and fix edge indexes for moved edge object
    edge_from_to[nlast->from()][nlast->to()] = nlast;
    edge_to_from[nlast->to()][nlast->from()] = nlast;

    // drop the last position, erasing the node
    graph.mutable_edge()->RemoveLast();

}

Node* VG::get_node(int64_t id) {
    map<int64_t, Node*>::iterator n = node_by_id.find(id);
    if (n != node_by_id.end()) {
        return n->second;
    } else {
        // again... should this throw an error?
        return NULL;
    }
}

// use the VG class to generate ids
Node* VG::create_node(string seq) {
    // create the node
    Node* node = graph.add_node();
    node->set_sequence(seq);
    node->set_id(graph.current_id());
    graph.set_current_id(graph.current_id()+1);
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = graph.node_size()-1;
    return node;
}

void VG::destroy_node(int64_t id) {
    destroy_node(get_node(id));
}

void VG::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    // noop on NULL/nonexistent node
    if (!has_node(node)) { return; }
    // remove edges associated with node
    set<Edge*> edges_to_destroy;
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(f->second);
        }
    }
    e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(f->second);
        }
    }
    for (set<Edge*>::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(*e);
    }
    // assert cleanup
    edge_to_from.erase(node->id());
    edge_from_to.erase(node->id());

    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = graph.node_size()-1;
    int tni = node_index[node];
    Node* last = graph.mutable_node(lni);
    graph.mutable_node()->SwapElements(tni, lni);
    Node* nlast = graph.mutable_node(tni);
    node_by_id[last->id()] = nlast;
    node_index.erase(last);
    node_index[nlast] = tni;
    node_by_id.erase(node->id());
    node_index.erase(node);
    graph.mutable_node()->RemoveLast();
    //if (!is_valid()) cerr << "graph is invalid after destroy_node" << endl;
}

// utilities
void VG::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    //cerr << "divide node " << node->id() << " @" << pos << endl;

    map<int64_t, map<int64_t, Edge*> >::iterator e;

    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    // replace node connections to prev (left)
    e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator p = e->second.begin();
             p != e->second.end(); ++p) {
            create_edge(p->first, left->id());
        }
    }

    // make our right node
    right = create_node(node->sequence().substr(pos,node->sequence().size()-1));

    // replace node connections to next (right)
    e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator n = e->second.begin();
             n != e->second.end(); ++n) {
            create_edge(right->id(), n->first);
        }
    }

    // connect left to right
    create_edge(left, right);

    destroy_node(node);

}

// for dividing a path of nodes with an underlying coordinate system
void VG::divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {

    map<long, Node*>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = target->second;
    
    // nothing to do
    if (node_pos == pos) {

        map<long, Node*>::iterator n = target; --n;
        left = n->second;
        right = target->second;

    } else {

        // divide the target node at our pos
        int diff = pos - node_pos;

        divide_node(old, diff, left, right);

        // left
        path[node_pos] = left;

        // right
        path[pos] = right;
    }
}

void VG::nodes_prev(Node* node, vector<Node*>& nodes) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin(); f != e->second.end(); ++f) {
            nodes.push_back(node_by_id[f->first]);
        }
    }
}

void VG::nodes_next(Node* node, vector<Node*>& nodes) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator t = e->second.begin(); t != e->second.end(); ++t) {
            nodes.push_back(node_by_id[t->first]);
        }
    }
}

void VG::bounded_prev_paths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
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
            bounded_prev_paths_from_node(*p, length - (*p)->sequence().size(), postfix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = postfix;
            new_path.push_front(*p);
            paths.insert(new_path);
        }
    }
}

void VG::bounded_next_paths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
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
            bounded_next_paths_from_node(*n, length - (*n)->sequence().size(), prefix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = prefix;
            new_path.push_back(*n);
            paths.insert(new_path);
        }
    }
}

void VG::bounded_paths(Node* node, set<list<Node*> >& paths, int length) {
    // get left, then right
    set<list<Node*> > prev_paths;
    set<list<Node*> > next_paths;
    list<Node*> empty_list;
    bounded_prev_paths_from_node(node, length, empty_list, prev_paths);
    bounded_next_paths_from_node(node, length, empty_list, next_paths);
    // now take the cross and put into paths
    for (set<list<Node*> >::iterator p = prev_paths.begin(); p != prev_paths.end(); ++p) {
        for (set<list<Node*> >::iterator n = next_paths.begin(); n != next_paths.end(); ++n) {
            list<Node*> path = *p;
            list<Node*>::const_iterator m = n->begin(); ++m; // skips current node, which is included in *p
            while (m != n->end()) {
                path.push_back(*m);
                ++m;
            }
            paths.insert(path);
        }
    }
}

void VG::bounded_paths(Node* node, vector<Path>& paths, int length) {
    set<list<Node*> > unique_paths;
    bounded_paths(node, unique_paths, length);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

void VG::bounded_paths(set<list<Node*> >& paths, int length) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        bounded_paths(node, paths, length);
    }
}

void VG::bounded_paths(vector<Path>& paths, int length) {
    set<list<Node*> > unique_paths;
    bounded_paths(unique_paths, length);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

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

void VG::node_starts_in_path(const list<Node*>& path, map<Node*, int>& node_start) {
    int i = 0;
    for (list<Node*>::const_iterator n = path.begin(); n != path.end(); ++n) {
        node_start[*n] = i;
        int l = (*n)->sequence().size();
        i += l;
    }
}

void VG::bounded_paths(int64_t node_id, vector<Path>& paths, int length) {
    map<int64_t, Node*>::iterator n = node_by_id.find(node_id);
    if (n != node_by_id.end()) {
        Node* node = n->second;
        bounded_paths(node, paths, length);
    }
}

string VG::path_sequence(Path& path) {
    string sequence;
    for (int i = 0; i < path.mapping_size(); ++i) {
        sequence.append(node_by_id[path.mapping(i).node_id()]->sequence());
    }
    return sequence;
}

bool VG::is_valid(void) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        int64_t f = e->from();
        int64_t t = e->to();
        if (node_by_id.find(f) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " cannot find node (from) " << f << endl;
            return false;
        }
        if (node_by_id.find(t) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " cannot find node (to) " << t << endl;
            return false;
        }
        if (edge_from_to.find(f) == edge_from_to.end()) {
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_from_to for node " << f << endl;
            return false;
        }
        if (edge_to_from.find(t) == edge_to_from.end()) {
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_to_from for node " << t << endl;
            return false;
        }
    }
    return true;
}

void VG::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "    " << p->id() << " -> " << n->id() << ";" << endl;
    }
    out << "}" << endl;
}

void VG::to_gfa(ostream& out) {
    out << "H" << "\t" << "HVN:Z:1.0" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        out << "S" << "\t" << n->id() << "\t" << n->sequence() << endl;
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "L" << "\t" << p->id() << "\t" << "-" << "\t" << n->id() << "\t" << "+" << "\t" << "0M" << endl;
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

Node* VG::join_heads(void) {
    Node* root = create_node("N");
    vector<Node*> heads;
    head_nodes(heads);
    connect_node_to_nodes(root, heads);
    return root;
}

Alignment& VG::align(Alignment& alignment) {

    // to be completely aligned, the graph's head nodes need to be fully-connected to a common root
    Node* root = join_heads();

    gssw_aligner = new GSSWAligner(graph);
    gssw_aligner->align(alignment);
    delete gssw_aligner;
    gssw_aligner = NULL;

    // adjust the alignment position to respect the trim of the first node
    Path* path = alignment.mutable_path();
    Node* first_node = node_by_id[path->mutable_mapping(0)->node_id()];
    if (first_node->has_trim()) {
        int32 t = path->target_position();
        path->set_target_position(t - first_node->trim());
    }

    // remove root
    destroy_node(root);
    
    return alignment;
}

Alignment VG::align(string& sequence) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment);
}


void VG::kmers_of(map<string, map<Node*, int> >& kmer_map, int kmer_size) {

    // get the bounded paths of the graph
    // these are paths of up to length kmer_size from a root node to neighbors
    set<list<Node*> > paths;

    // we use half the kmer size + 1 here
    bounded_paths(paths, kmer_size/2 + 1);

    // the kmer_map is the map which we will use to temporarily store the mappings
    // this could be really big before serializing to disk

    // for each path...
    for (set<list<Node*> >::iterator p = paths.begin(); p != paths.end(); ++p) {
        const list<Node*>& path = *p;

        // expand the path into a vector :: 1,1,1,2,2,2,2,3,3 ... etc.
        // this makes it much easier to quickly get all the node matches of each kmer
        vector<Node*> node_by_path_position;
        expand_path(path, node_by_path_position);

        map<Node*, int> node_start;
        node_starts_in_path(path, node_start);
        
        // now process the kmers of this sequence
        // by first getting the sequence
        string seq = path_string(path);

        // and then stepping across the path, finding the kmers, and then implied node overlaps
        for (int i = 0; i < seq.size() - kmer_size; ++i) {

            // get the kmer
            string kmer = seq.substr(i, kmer_size);

            // create our kmer node match set if it doesn't exist
            map<Node*, int>& nodes = kmer_map[kmer];

            // drop the node overlaps of this kmer into our index
            int j = 0;
            while (j < kmer_size) {
                Node* node = node_by_path_position[i+j];
                int node_position = node_start[node];
                int kmer_relative_start = i - node_position;
                nodes[node] = kmer_relative_start;
                j += node->sequence().size();
            }
        }
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
        if (edge_to_from.find(n->id()) == edge_to_from.end()) {
            nodes.push_back(n);
        }
    }
}


    /*
Tarjan's topological sort

L <- Empty list that will contain the sorted nodes
while there are unmarked nodes do
    select an unmarked node n
    visit(n) 
function visit(node n)
    if n has a temporary mark then stop (not a DAG)
    if n is not marked (i.e. has not been visited yet) then
        mark n temporarily
        for each node m with an edge from n to m do
            visit(m)
        mark n permanently
        add n to head of L
    */

void VG::topological_sort(list<Node*>& sorted_nodes) {
    set<Node*> unmarked_nodes;
    set<Node*> temporary_marks;
    for (int i = 0; i < graph.node_size(); ++i) {
        unmarked_nodes.insert(graph.mutable_node(i));
    }
    while (!unmarked_nodes.empty()) {
        Node* node = *(unmarked_nodes.begin());
        visit_node(node,
                   sorted_nodes,
                   unmarked_nodes,
                   temporary_marks);
    }
}

void VG::visit_node(Node* node,
                              list<Node*>& sorted_nodes,
                              set<Node*>& unmarked_nodes,
                              set<Node*>& temporary_marks) {
    if (unmarked_nodes.find(node) != unmarked_nodes.end()) {
        temporary_marks.insert(node);
        map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(node->id());
        if (e != edge_from_to.end()) {
            for (map<int64_t, Edge*>::iterator f = e->second.begin(); f != e->second.end(); ++f) {
                Node* to = node_by_id[f->second->to()];
                visit_node(to,
                           sorted_nodes,
                           unmarked_nodes,
                           temporary_marks);
            }
        }
        unmarked_nodes.erase(node);
        sorted_nodes.push_front(node);
    }
}
