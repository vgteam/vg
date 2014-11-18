#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include "gssw.h"
#include "gssw_aligner.h"
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

    void topologically_sort_graph(void);
    void topological_sort(list<Node*>& sorted_nodes);
    void visit_node(Node* node,
                    list<Node*>& sorted_nodes,
                    set<Node*>& unmarked_nodes,
                    set<Node*>& temporary_marks);
    void swap_nodes(Node* a, Node* b);

    Alignment align(string& sequence);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

private:

    void init(void); // setup, ensures that gssw == NULL on startup

};

} // end namespace vg

#endif
