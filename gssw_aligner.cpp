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
        gssw_node* node = (gssw_node*)gssw_node_create(n, n->id(),
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

void GSSWAligner::gssw_mapping_to_alignment(gssw_graph_mapping* gm,
                                            Alignment& alignment) {
    alignment.clear_path();
    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    Path* path = alignment.mutable_path();
    //alignment.set_cigar(graph_cigar(gm));

    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //cerr << "gm->position " << gm->position << endl;
    string& to_seq = *alignment.mutable_sequence();
    //cerr << "-------------" << endl;

    //gssw_graph_print_score_matrices(graph, to_seq.c_str(), to_seq.size(), stderr);


    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        // check that the current alignment has a sensible length

        Node* from_node = (Node*) nc->node->data;
        string& from_seq = *from_node->mutable_sequence();
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(nc->node->id);
        mapping->mutable_position()->set_offset(from_pos);
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        //cerr << from_node->id() << ":" << endl;

        for (int j=0; j < l; ++j, ++e) {
            Edit* edit;
            int32_t length = e->length;
            //cerr << e->length << e->type << endl;
            switch (e->type) {
            case 'M': {
                // do the sequences match?
                // emit a stream of "SNPs" and matches
                int h = from_pos;
                int last_start = from_pos;
                int k = to_pos;
                for ( ; h < from_pos + length; ++h, ++k) {
                    //cerr << h << ":" << k << " " << from_seq[h] << " " << to_seq[k] << endl;
                    if (from_seq[h] != to_seq[k]) {
                        // emit the last "match" region
                        if (h-last_start > 0) {
                            edit = mapping->add_edit();
                            edit->set_from_length(h-last_start);
                            edit->set_to_length(h-last_start);
                        }
                        // set up the SNP
                        edit = mapping->add_edit();
                        edit->set_from_length(1);
                        edit->set_to_length(1);
                        edit->set_sequence(to_seq.substr(k,1));
                        last_start = h+1;
                    }
                }
                // handles the match at the end or the case of no SNP
                if (h-last_start > 0) {
                    edit = mapping->add_edit();
                    edit->set_from_length(h-last_start);
                    edit->set_to_length(h-last_start);
                }
                to_pos += length;
                from_pos += length;
            } break;
            case 'D':
                edit = mapping->add_edit();
                edit->set_from_length(length);
                edit->set_to_length(0);
                from_pos += length;
                break;
            case 'I':
                edit = mapping->add_edit();
                edit->set_from_length(0);
                edit->set_to_length(length);
                edit->set_sequence(to_seq.substr(to_pos, length));
                to_pos += length;
                break;
            case 'S':
                edit = mapping->add_edit();
                edit->set_from_length(0);
                edit->set_to_length(length);
                to_pos += length;
                break;
            default:
                cerr << "error [GSSWAligner::gssw_mapping_to_alignment] "
                     << "unsupported cigar op type " << e->type << endl;
                exit(1);
                break;

            }
        }
        //cerr << "path to_length " << path_to_length(*path) << endl;
    }
}

string GSSWAligner::graph_cigar(gssw_graph_mapping* gm) {
    stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //string& to_seq = *alignment.mutable_sequence();
    s << from_pos << '@';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        Node* from_node = (Node*) nc->node->data;
        s << from_node->id() << ':';
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            s << e->length << e->type;
        }
        if (i + 1 < gc->length) {
            s << ",";
        }
    }
    return s.str();
}
