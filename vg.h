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
    VariantGraph(void);
    // construct from protobufs
    VariantGraph(istream& in);
    VariantGraph(Graph& graph);
    // construct from sets of nodes and edges (e.g. subgraph of another graph)
    VariantGraph(set<Node*>& nodes, set<Edge*>& edges);
    // construct from VCF records
    VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference);
    ~VariantGraph(void);
    VariantGraph& operator=(const VariantGraph& other) {
        if (this != &other) {
            // cleanup
            clear_indexes();
            // assign
            graph = other.graph;
            // re-index
            build_indexes();
        }
        return *this;
    }

    void build_indexes(void);
    void clear_indexes(void);
    void extend(Graph& g);
    void extend(VariantGraph& g);

    bool node_exists(Node& node);
    bool edge_exists(Edge& edge);

    void add_node(Node& node);
    void add_nodes(vector<Node>& nodes);
    void add_edge(Edge& edge);
    void add_edges(vector<Edge>& edges);
    void add_nodes(set<Node*>& nodes);
    void add_edges(set<Edge*>& edges);

    int64_t node_count(void);
    int64_t edge_count(void);
    int64_t total_length_of_nodes(void);
    int in_degree(Node* node);
    int out_degree(Node* node);
    void edges_of_node(Node* node, vector<Edge*>& edges);
    void edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges);

    // use the VariantGraph class to generate ids
    Node* create_node(string seq);
    void destroy_node(Node* node);

    // edges
    Edge* create_edge(Node* from, Node* to);
    Edge* create_edge(int64_t from, int64_t to);
    void destroy_edge(Edge* edge);

    // connect node -> nodes
    void connect_node_to_nodes(Node* node, vector<Node*>& nodes);

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

    Alignment& align(Alignment& alignment);
    Alignment align(string& sequence);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    void bounded_paths(vector<Path>& paths, int length);
    void bounded_paths(set<list<Node*> >& paths, int length);
    void bounded_paths(Node* node, set<list<Node*> >& paths, int length);
    void bounded_paths(Node* node, vector<Path>& paths, int length);
    void bounded_paths(int64_t node_id, vector<Path>& paths, int length);
    void bounded_prev_paths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths);
    void bounded_next_paths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(int64_t from, int64_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    // traversal
    void nodes_prev(Node* n, vector<Node*>& nodes);
    void nodes_next(Node* n, vector<Node*>& nodes);

    // paths
    Path create_path(const list<Node*>& nodes);
    Path create_path(const vector<Node*>& nodes);
    string path_string(const list<Node*>& nodes);
    string path_string(Path& path);
    void expand_path(const list<Node*>& path, vector<Node*>& expanded);
    void node_starts_in_path(const list<Node*>& path, map<Node*, int>& node_start);

    // kmers
    void kmers_of(map<string, map<Node*, int> >& kmer_map, int kmer_size);

    // subgraphs
    void disjoint_subgraphs(list<VariantGraph>& subgraphs);
    void head_nodes(vector<Node*>& nodes);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node
    Node* join_heads(void);

private:

    void init(void); // setup, ensures that gssw == NULL on startup

};


} // end namespace vg

#endif
