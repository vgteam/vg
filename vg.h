#ifndef VG_H
#define VG_H

#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"

namespace vg {


VariantGraph {

public:

    // nodes by id
    map<int64_t, Node*> nodes;

    // edges by id
    //map<int64_t, Edge*> edges;

    // constructors
    VariantGraph(void) { };
    // construct from protobufs
    VariantGraph(Graph& graph);
    VariantGraph(vector<Node>& nodes);

    // construct from VCF records
    VariantGraph(vector<vcf::Variant>& variants, FastaReference& ref);

    // use the VariantGraph class to generate ids
    Node* create_node(string& seq);
    void destroy_node(Node* node);

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_ref_path(map<long, Node*>& ref_path, long pos, Node*& left, Node*& right);
    void replace_prev(Node* node, Node* from, Node* to);
    void replace_next(Node* node, Node* from, Node* to);

    //void align(Alignment& alignment);

    //void topological_sort(void); // possibly unnecissary

};


}

#endif
