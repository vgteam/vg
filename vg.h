#ifndef VG_H
#define VG_H

#include <vector>
#include <string>
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"


namespace vg {


class VariantGraph : public Graph {

public:

    // nodes by id
    map<int64_t, Node*> node_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    map<Node*, int> node_index;

    // edges indexed by nodes they connect
    map<pair<int64_t, int64_t>, Edge*> edge_by_ids;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    map<Edge*, int> edge_index;

    // constructors
    //VariantGraph(void) { };
    // construct from protobufs
    VariantGraph(ifstream& in);
    VariantGraph(Graph& graph);
    VariantGraph(vector<Node>& nodes);

    // construct from VCF records
    VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference);

    // use the VariantGraph class to generate ids
    Node* create_node(string seq);
    void destroy_node(Node* node);

    Edge* create_edge(int64_t from, int64_t to);
    void destroy_edge(Edge* edge);
    Edge* add_edge(Node* from, Node* to);
    void remove_edge_from_node(Node* node, Edge* edge);

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right);
    void node_replace_prev(Node* node, Node* before, Node* after);
    void node_replace_next(Node* node, Node* before, Node* after);

    //void align(Alignment& alignment);

    //void topological_sort(void); // possibly unnecissary

};


} // end namespace vg

#endif
