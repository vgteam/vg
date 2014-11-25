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
    VariantGraph(void);
    VariantGraph(istream& in);
    VariantGraph(Graph& graph);
    VariantGraph(vector<Node>& nodes);
    ~VariantGraph(void);

    void clear_indexes(void);
    void extend(Graph& g);
    void extend(VariantGraph& g);

    bool node_exists(Node& node);
    bool edge_exists(Edge& edge);

    void add_node(Node& node);
    void add_nodes(vector<Node>& nodes);
    void add_edge(Edge& edge);
    void add_edges(vector<Edge>& edges);

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

    Alignment align(string& sequence, bool reuse_gssw = false);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    void bounded_paths(Node* node, vector<Path>& paths, int length);
    void bounded_paths(int64_t node_id, vector<Path>& paths, int length);
    void bounded_prev_paths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths);
    void bounded_next_paths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths);
    void paths(Node* from, Node* to, vector<Path>& paths);
    void paths(int64_t from, int64_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    void nodes_prev(Node* n, vector<Node*>& nodes);
    void nodes_next(Node* n, vector<Node*>& nodes);

    // create paths
    Path create_path(const list<Node*>& nodes);
    Path create_path(const vector<Node*>& nodes);

private:

    void init(void); // setup, ensures that gssw == NULL on startup

};

} // end namespace vg

#endif
