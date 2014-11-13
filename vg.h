#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include "gssw.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"


namespace vg {


class VariantGraph {

public:

    // protobuf-based representation
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    // nodes by id
    map<int64_t, Node*> node_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    map<Node*, int> node_index;

    // edges indexed by nodes they connect
    map<int64_t, map<int64_t, Edge*> > edge_from_to;
    map<int64_t, map<int64_t, Edge*> > edge_to_from;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    map<Edge*, int> edge_index;

    // constructors
    //VariantGraph(void) { };
    // construct from protobufs
    VariantGraph(istream& in);
    VariantGraph(Graph& graph);
    VariantGraph(vector<Node>& nodes);
    ~VariantGraph(void);

    // construct from VCF records
    VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference);

    // use the VariantGraph class to generate ids
    Node* create_node(string seq);
    void destroy_node(Node* node);

    Edge* create_edge(Node* from, Node* to);
    Edge* create_edge(int64_t from, int64_t to);
    void destroy_edge(Edge* edge);

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out);
    bool is_valid(void);

    gssw_graph* create_alignable_graph(
        int32_t match = 2,
        int32_t mismatch = 2,
        int32_t gap_open = 3,
        int32_t gap_extension = 1);
    void destroy_alignable_graph(void);

    Alignment align(string& sequence);
    Alignment& align(Alignment& alignment);
    void gssw_mapping_to_alignment(gssw_graph_mapping* gm, Alignment& alignment);

    // needed when constructing an alignable graph from the nodes
    void topological_sort(list<gssw_node*>& sorted_nodes);

    void visit_node(gssw_node* node,
                    list<gssw_node*>& sorted_nodes,
                    set<gssw_node*>& unmarked_nodes,
                    set<gssw_node*>& temporary_marks);

private:

    void init(void); // setup, ensures that _gssw_graph == NULL on startup
    map<int64_t, gssw_node*> _gssw_nodes;
    gssw_graph* _gssw_graph;
    int8_t* _gssw_nt_table;
    int8_t* _gssw_score_matrix;
    int32_t _gssw_match;
    int32_t _gssw_mismatch;
    int32_t _gssw_gap_open;
    int32_t _gssw_gap_extension;

};


} // end namespace vg

#endif
