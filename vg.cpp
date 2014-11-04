#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;



//VariantGraph::VariantGraph(void) { };
// construct from protobufs
VariantGraph::VariantGraph(Graph& graph) {
    
}

VariantGraph::VariantGraph(vector<Node>& nodes) {

}

// construct from VCF records
VariantGraph::VariantGraph(vector<vcf::Variant>& variants) {

}

// use the VariantGraph class to generate ids
Node* VariantGraph::create_node(string& seq) {

}

void VariantGraph::destroy_node(Node* node) {

}

// utilities
void VariantGraph::divide_node(Node* node, int pos, Node*& left, Node*& right) {

}

void VariantGraph::divide_ref_path(map<long, Node*>& ref_path, long pos, Node*& left, Node*& right) {

}

void VariantGraph::replace_prev(Node* node, Node* from, Node* to) {

}

void VariantGraph::replace_next(Node* node, Node* from, Node* to) {

}

//void align(Alignment& alignment);
