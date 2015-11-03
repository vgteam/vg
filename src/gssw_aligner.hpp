#ifndef GSSW_ALIGNER_H
#define GSSW_ALIGNER_H

#include <vector>
#include <set>
#include <string>
#include "gssw.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"
#include "path.hpp"

namespace vg {


class GSSWAligner {
public:

    GSSWAligner(
        Graph& g,
        int32_t _match = 2,
        int32_t _mismatch = 2,
        int32_t _gap_open = 3,
        int32_t _gap_extension = 1);

    ~GSSWAligner(void);

    // for construction
    // needed when constructing an alignable graph from the nodes
    void topological_sort(list<gssw_node*>& sorted_nodes);
    void visit_node(gssw_node* node,
                    list<gssw_node*>& sorted_nodes,
                    set<gssw_node*>& unmarked_nodes,
                    set<gssw_node*>& temporary_marks);

    // alignment functions
    void align(Alignment& alignment);
    void gssw_mapping_to_alignment(gssw_graph_mapping* gm, Alignment& alignment);
    string graph_cigar(gssw_graph_mapping* gm);

    // members
    map<int64_t, gssw_node*> nodes;
    gssw_graph* graph;
    int8_t* nt_table;
    int8_t* score_matrix;
    int32_t match;
    int32_t mismatch;
    int32_t gap_open;
    int32_t gap_extension;

};

} // end namespace vg

#endif
