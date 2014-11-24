#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;


//VariantGraph::VariantGraph(void) { };
// construct from protobufs
VariantGraph::VariantGraph(istream& in) {
    init();
    graph.ParseFromIstream(&in);
    // populate by-id node index
    for (int64_t i = 0; i < graph.nodes_size(); ++i) {
        Node* n = graph.mutable_nodes(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
        edge_index[e] = i;
        edge_from_to[e->from()][e->to()] = e;
        edge_to_from[e->to()][e->from()] = e;
    }
}

VariantGraph::~VariantGraph(void) {
    destroy_alignable_graph();
}

VariantGraph::VariantGraph(void) {
    init();
}

void VariantGraph::init(void) {
    gssw_aligner = NULL;
}

// TODO add edges
VariantGraph::VariantGraph(vector<Node>& nodesv) {
    init();
    for (vector<Node>::iterator n = nodesv.begin(); n != nodesv.end(); ++n) {
        Node& node = *n;
        int64_t id = node.id();
        if (graph.current_id() < id) {
            graph.set_current_id(id);
        }
        Node* new_node = graph.add_nodes(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.nodes_size()-1;
    }
}

// these are quite unsafe--- perhaps they should check if there are conflicts?

void VariantGraph::add_nodes(vector<Node>& nodes) {
    for (vector<Node>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        add_node(*n);
    }
}

void VariantGraph::add_edges(vector<Edge>& edges) {
    for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
        add_edge(*e);
    }
}

bool VariantGraph::node_exists(Node& node) {
    return node_by_id.find(node.id()) != node_by_id.end();
}

bool VariantGraph::edge_exists(Edge& edge) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(edge.from());
    return e != edge_from_to.end() && e->second.find(edge.to()) != e->second.end();
}

void VariantGraph::add_node(Node& node) {
    if (!node_exists(node)) {
        Node* new_node = graph.add_nodes(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.nodes_size()-1;
    }
}

void VariantGraph::add_edge(Edge& edge) {
    if (!edge_exists(edge)) {
        Edge* new_edge = graph.add_edges(); // add it to the graph
        new_edge->set_from(edge.from());
        new_edge->set_to(edge.to());
        edge_from_to[edge.from()][edge.to()] = new_edge;
        edge_to_from[edge.to()][edge.from()] = new_edge;
        edge_index[new_edge] = graph.edges_size()-1;
    }
}

void VariantGraph::clear_indexes(void) {
    node_index.clear();
    node_by_id.clear();
    edge_index.clear();
    edge_from_to.clear();
    edge_to_from.clear();
}

void VariantGraph::extend(VariantGraph& g) {
    extend(g.graph);
}

void VariantGraph::extend(Graph& g) {
    graph.mutable_nodes()->MergeFrom(g.nodes());
    graph.mutable_edges()->MergeFrom(g.edges());
    clear_indexes();
    for (int64_t i = 0; i < graph.nodes_size(); ++i) {
        Node* n = graph.mutable_nodes(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
        edge_index[e] = i;
        edge_from_to[e->from()][e->to()] = e;
        edge_to_from[e->to()][e->from()] = e;
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
VariantGraph::VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference) {
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

void VariantGraph::topologically_sort_graph(void) {
    list<Node*> sorted_nodes;
    topological_sort(sorted_nodes);
    list<Node*>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.nodes_size() && n != sorted_nodes.end();
         ++i, ++n) {
        swap_nodes(graph.mutable_nodes(i), *n);
    }
}

void VariantGraph::swap_nodes(Node* a, Node* b) {
    int aidx = node_index[a];
    int bidx = node_index[b];
    graph.mutable_nodes()->SwapElements(aidx, bidx);
    node_index[a] = bidx;
    node_index[b] = aidx;
}

Edge* VariantGraph::create_edge(Node* from, Node* to) {
    return create_edge(from->id(), to->id());
}

Edge* VariantGraph::create_edge(int64_t from, int64_t to) {
    // prevent self-linking (violates DAG/partial ordering property)
    if (to == from) return NULL;
    // ensure the edge does not already exist
    Edge* edge = edge_from_to[from][to];
    if (edge) return edge;
    // if not, create it
    edge = graph.add_edges();
    edge->set_from(from);
    edge->set_to(to);
    edge_from_to[from][to] = edge;
    edge_to_from[to][from] = edge;
    edge_index[edge] = graph.edges_size()-1;
    return edge;
}

void VariantGraph::destroy_edge(Edge* edge) {
    //if (!is_valid()) cerr << "graph ain't valid" << endl;
    // erase from indexes
    edge_from_to[edge->from()].erase(edge->to());
    edge_to_from[edge->to()].erase(edge->from());

    // erase from edges by moving to end and dropping
    int lei = graph.edges_size()-1;
    int tei = edge_index[edge];
    Edge* last = graph.mutable_edges(lei);
    edge_index.erase(last);
    edge_index.erase(edge);

    // swap
    graph.mutable_edges()->SwapElements(tei, lei);

    // point to new position
    Edge* nlast = graph.mutable_edges(tei);

    // insert the new edge index position
    edge_index[nlast] = tei;

    // and fix edge indexes for moved edge object
    edge_from_to[nlast->from()][nlast->to()] = nlast;
    edge_to_from[nlast->to()][nlast->from()] = nlast;

    // drop the last position, erasing the node
    graph.mutable_edges()->RemoveLast();

}

// use the VariantGraph class to generate ids
Node* VariantGraph::create_node(string seq) {
    // create the node
    Node* node = graph.add_nodes();
    node->set_sequence(seq);
    node->set_id(graph.current_id());
    graph.set_current_id(graph.current_id()+1);
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = graph.nodes_size()-1;
    return node;
}

void VariantGraph::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
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
    int lni = graph.nodes_size()-1;
    int tni = node_index[node];
    Node* last = graph.mutable_nodes(lni);
    graph.mutable_nodes()->SwapElements(tni, lni);
    Node* nlast = graph.mutable_nodes(tni);
    node_by_id[last->id()] = nlast;
    node_index.erase(last);
    node_index[nlast] = tni;
    node_by_id.erase(node->id());
    node_index.erase(node);
    graph.mutable_nodes()->RemoveLast();
    //if (!is_valid()) cerr << "graph is invalid after destroy_node" << endl;
}

// utilities
void VariantGraph::divide_node(Node* node, int pos, Node*& left, Node*& right) {

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
void VariantGraph::divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {

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

void VariantGraph::nodes_prev(Node* node, vector<Node*>& nodes) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin(); f != e->second.end(); ++f) {
            nodes.push_back(node_by_id[f->first]);
        }
    }
}

void VariantGraph::nodes_next(Node* node, vector<Node*>& nodes) {
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator t = e->second.begin(); t != e->second.end(); ++t) {
            nodes.push_back(node_by_id[t->first]);
        }
    }
}

void VariantGraph::bounded_prev_paths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    postfix.push_front(node);
    vector<Node*> prev_nodes;
    nodes_prev(node, prev_nodes);
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

void VariantGraph::bounded_next_paths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    prefix.push_back(node);
    vector<Node*> next_nodes;
    nodes_next(node, next_nodes);
    for (vector<Node*>::iterator n = next_nodes.begin(); n != next_nodes.end(); ++n) {
        if ((*n)->sequence().size() < length) {
            bounded_prev_paths_from_node(*n, length - (*n)->sequence().size(), prefix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = prefix;
            new_path.push_back(*n);
            paths.insert(new_path);
        }
    }
}

void VariantGraph::bounded_paths(Node* node, vector<list<Node*> >& paths, int length) {
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
            paths.push_back(path);
        }
    }
}

bool VariantGraph::is_valid(void) {
    for (int i = 0; i < graph.nodes_size(); ++i) {
        Node* n = graph.mutable_nodes(i);
    }
    for (int i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
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

void VariantGraph::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    for (int i = 0; i < graph.nodes_size(); ++i) {
        Node* n = graph.mutable_nodes(i);
        out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;
    }
    for (int i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "    " << p->id() << " -> " << n->id() << ";" << endl;
    }
    out << "}" << endl;
}

void VariantGraph::destroy_alignable_graph(void) {
    if (gssw_aligner != NULL) {
        delete gssw_aligner;
    }
}

Alignment VariantGraph::align(string& sequence, bool reuse_gssw) {

    if (reuse_gssw && !gssw_aligner) {
        gssw_aligner = new GSSWAligner(graph);
    }

    Alignment alignment;
    alignment.set_sequence(sequence);
    gssw_aligner->align(alignment);

    if (!reuse_gssw) {
        delete gssw_aligner;
        gssw_aligner = NULL;
    }

    return alignment;

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

void VariantGraph::topological_sort(list<Node*>& sorted_nodes) {
    set<Node*> unmarked_nodes;
    set<Node*> temporary_marks;
    for (int i = 0; i < graph.nodes_size(); ++i) {
        unmarked_nodes.insert(graph.mutable_nodes(i));
    }
    while (!unmarked_nodes.empty()) {
        Node* node = *(unmarked_nodes.begin());
        visit_node(node,
                   sorted_nodes,
                   unmarked_nodes,
                   temporary_marks);
    }
}

void VariantGraph::visit_node(Node* node,
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
