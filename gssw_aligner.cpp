#include "gssw_aligner.h"

using namespace vg;

GSSWAligner::~GSSWAligner(void) {
    gssw_graph_destroy(graph);
    free(nt_table);
    free(score_matrix);
}

GSSWAligner::GSSWAligner(
    Graph& g,
    int32_t _match,
    int32_t _mismatch,
    int32_t _gap_open,
    int32_t _gap_extension
) {

    match = _match;
    mismatch = _mismatch;
    gap_open = _gap_open;
    gap_extension = _gap_extension;

    // these are used when setting up the nodes
    // they can be cleaned up via destroy_alignable_graph()
    nt_table = gssw_create_nt_table();
	score_matrix = gssw_create_score_matrix(match, mismatch);

    graph = gssw_graph_create(g.nodes_size());

    for (int i = 0; i < g.nodes_size(); ++i) {
        Node* n = g.mutable_nodes(i);
        nodes[n->id()] = (gssw_node*)gssw_node_create(NULL, n->id(),
                                                      n->sequence().c_str(),
                                                      nt_table,
                                                      score_matrix);
    }

    for (int i = 0; i < g.edges_size(); ++i) {
        Edge* e = g.mutable_edges(i);
        gssw_nodes_add_edge(nodes[e->from()], nodes[e->to()]);
    }

    list<gssw_node*> sorted_nodes;
    topological_sort(sorted_nodes);

    for (list<gssw_node*>::iterator n = sorted_nodes.begin(); n != sorted_nodes.end(); ++n) {
        gssw_graph_add_node(graph, *n);
    }

}


void GSSWAligner::align(Alignment& alignment) {

    const string& sequence = alignment.sequence();

    gssw_graph_fill(graph, sequence.c_str(),
                    nt_table, score_matrix,
                    gap_open, gap_extension, 15, 2);

    //gssw_graph_print_score_matrices(_gssw_graph, sequence.c_str(), sequence.size());
    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    sequence.c_str(),
                                                    sequence.size(),
                                                    match,
                                                    mismatch,
                                                    gap_open,
                                                    gap_extension);

    gssw_mapping_to_alignment(gm, alignment);

    //gssw_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);

}

void GSSWAligner::gssw_mapping_to_alignment(gssw_graph_mapping* gm, Alignment& alignment) {

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
                cerr << "error [GSSWAligner::gssw_mapping_to_alignment] "
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

void GSSWAligner::topological_sort(list<gssw_node*>& sorted_nodes) {
    set<gssw_node*> unmarked_nodes;
    set<gssw_node*> temporary_marks;
    for (map<int64_t, gssw_node*>::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
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

void GSSWAligner::visit_node(gssw_node* node,
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
