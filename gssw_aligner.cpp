#include "gssw_aligner.hpp"

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

    graph = gssw_graph_create(g.node_size());

    for (int i = 0; i < g.node_size(); ++i) {
        Node* n = g.mutable_node(i);
        gssw_node* node = (gssw_node*)gssw_node_create(NULL, n->id(),
                                                       n->sequence().c_str(),
                                                       nt_table,
                                                       score_matrix);
        nodes[n->id()] = node;
        gssw_graph_add_node(graph, node);
    }

    for (int i = 0; i < g.edge_size(); ++i) {
        Edge* e = g.mutable_edge(i);
        gssw_nodes_add_edge(nodes[e->from()], nodes[e->to()]);
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
    Path* path = alignment.mutable_path();
    path->set_target_position(gm->position);

    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;

    for (int i = 0; i < gc->length; ++i, ++nc) {

        Mapping* mapping = path->add_mapping();
        mapping->set_node_id(nc->node->id);
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;

        for (int j=0; j < l; ++j, ++e) {

            Edit* edit = mapping->add_edit();
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
