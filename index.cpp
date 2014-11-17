#include "index.h"

namespace vg {

using namespace std;

Index::Index(string& name) {
    options.create_if_missing = true;
    options.error_if_exists = true;
    leveldb::Status status = leveldb::DB::Open(options, name, &db);
    if (!status.ok()) {
        throw indexOpenException;
    }
}

Index::~Index(void) {
    delete db;
}

void Index::put_node(Node& node) {
    string data;
    n.SerializeToString(&data);
    string key = "n" + '\xff' + (string) n.id();
    db->Put(leveldb::WriteOptions(), key, data);
}

void Index::load_graph(VariantGraph& graph) {
    for (int i = 0; i < graph.nodes_size(); ++i) {
        put_node(graph.nodes(i));
        //out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;

    }
    for (int i = 0; i < graph.edges_size(); ++i) {
        Edge* e = graph.mutable_edges(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "    " << p->id() << " -> " << n->id() << ";" << endl;
    }
}

void index_kmers(VariantGraph& graph, int kmer_size = 15) {
// 
}

void index_positions(VariantGraph& graph, map<long, Node*>& node_path, map<long, Edge*>& edge_path) {

}

}
