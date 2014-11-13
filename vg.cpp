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

void VariantGraph::init(void) {
    _gssw_graph = NULL;
}

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
        //Node& new_node = nodes(graph.nodes_size()-1); // get a reference to it
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.nodes_size()-1;
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
    if (_gssw_graph) {
        gssw_graph_destroy(_gssw_graph);
        _gssw_nodes.clear(); // these are freed via gssw_graph_destroy
        free(_gssw_nt_table);
        free(_gssw_score_matrix);
    }
}

gssw_graph* VariantGraph::create_alignable_graph(
    int32_t match,
    int32_t mismatch,
    int32_t gap_open,
    int32_t gap_extension
) {

    _gssw_match = match;
    _gssw_mismatch = mismatch;
    _gssw_gap_open = gap_open;
    _gssw_gap_extension = gap_extension;

    // these are used when setting up the nodes
    // they can be cleaned up via destroy_alignable_graph()
    _gssw_nt_table = gssw_create_nt_table();
	_gssw_score_matrix = gssw_create_score_matrix(_gssw_match, _gssw_mismatch);

    _gssw_graph = gssw_graph_create(graph.nodes_size());

    for (int i = 0; i < graph.nodes_size(); ++i) {
        Node* n = graph.mutable_nodes(i);
        _gssw_nodes[n->id()] = (gssw_node*)gssw_node_create(NULL, n->id(),
                                                           n->sequence().c_str(),
                                                           _gssw_nt_table,
                                                           _gssw_score_matrix);
    }

    for (int i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
        gssw_nodes_add_edge(_gssw_nodes[e->from()], _gssw_nodes[e->to()]);
    }

    list<gssw_node*> sorted_nodes;
    topological_sort(sorted_nodes);

    for (list<gssw_node*>::iterator n = sorted_nodes.begin(); n != sorted_nodes.end(); ++n) {
        gssw_graph_add_node(_gssw_graph, *n);
    }

}

Alignment VariantGraph::align(string& sequence) {
    
    gssw_graph_fill(_gssw_graph, sequence.c_str(),
                    _gssw_nt_table, _gssw_score_matrix,
                    _gssw_gap_open, _gssw_gap_extension, 15, 2);

    //gssw_graph_print_score_matrices(_gssw_graph, sequence.c_str(), sequence.size());
    gssw_graph_mapping* gm = gssw_graph_trace_back (_gssw_graph,
                                                    sequence.c_str(),
                                                    sequence.size(),
                                                    _gssw_match,
                                                    _gssw_mismatch,
                                                    _gssw_gap_open,
                                                    _gssw_gap_extension);


    Alignment alignment;
    alignment.set_sequence(sequence);
    gssw_mapping_to_alignment(gm, alignment);

    //gssw_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);

    return alignment;

}

void VariantGraph::gssw_mapping_to_alignment(gssw_graph_mapping* gm, Alignment& alignment) {

    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    alignment.set_target_position(gm->position);
    Path* path = alignment.mutable_path();

    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;

    for (int i = 0; i < gc->length; ++i, ++nc) {

        Mapping* mapping = path->add_nodes();
        mapping->set_node_id(nc->node->id);
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;

        for (int j=0; j < l; ++j, ++e) {

            Edit* edit = mapping->add_edits();
            edit->set_length(e->length);

            switch (e->type) {
            case 'M':
                edit->set_type(Edit_Type_MATCH);
                break;
            case 'D':
                edit->set_type(Edit_Type_DELETION);
                break;
            case 'I':
                edit->set_type(Edit_Type_INSERTION);
                break;
            case 'S':
                edit->set_type(Edit_Type_SOFTCLIP);
                break;
            default:
                cerr << "error [VariantGraph::gssw_mapping_to_alignment] "
                     << "unsupported cigar op type " << e->type << endl;
                exit(1);
                break;

            }
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

void VariantGraph::topological_sort(list<gssw_node*>& sorted_nodes) {
    set<gssw_node*> unmarked_nodes;
    set<gssw_node*> temporary_marks;
    for (map<int64_t, gssw_node*>::iterator n = _gssw_nodes.begin();
         n != _gssw_nodes.end(); ++n) {
        unmarked_nodes.insert(n->second);
    }
    while (!unmarked_nodes.empty()) {
        gssw_node* node = *(unmarked_nodes.begin());
        visit_node(node,
                   sorted_nodes,
                   unmarked_nodes,
                   temporary_marks);
    }
}

void VariantGraph::visit_node(gssw_node* node,
                list<gssw_node*>& sorted_nodes,
                set<gssw_node*>& unmarked_nodes,
                set<gssw_node*>& temporary_marks) {
    /*
    if (temporary_marks.find(node) != temporary_marks.end()) {
        cerr << "cannot sort graph because it is not a DAG!" << endl;
        exit(1);
    }
    */
    if (unmarked_nodes.find(node) != unmarked_nodes.end()) {
        temporary_marks.insert(node);
        for (int i = 0; i < node->count_next; ++i) {
            visit_node(node->next[i],
                       sorted_nodes,
                       unmarked_nodes,
                       temporary_marks);
        }
        unmarked_nodes.erase(node);
        sorted_nodes.push_front(node);
    }
}

