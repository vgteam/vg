#ifndef VG_HAPLOTYPE_EXTRACTER_HPP_INCLUDED
#define VG_HAPLOTYPE_EXTRACTER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"

using namespace std;
using namespace vg;
using thread_t = vector<xg::XG::ThreadMapping>;

// Walk forward from a node, collecting all haplotypes.  Also do a regular
// subgraph search for all the paths too.  Haplotype thread i will be
// embedded as Paths a path with name thread_i.
// Each path name (including threads) is mapped to a frequency in out_thread_frequencies
void trace_haplotypes_and_paths(xg::XG& index,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph = true);

// Turns an (xg-based) thread_t into a (vg-based) Path
Path path_from_thread_t(thread_t& t);

// Lists all the sub-haplotypes of length extend_distance nodes starting at node
// start_node from the set of haplotypes embedded as thread_t's in xg index.
// Records, for each thread_t t the number of haplotypes of which t is a
// subhaplotype
vector<pair<thread_t,int> > list_haplotypes(xg::XG& index,
            xg::XG::ThreadMapping start_node, int extend_distance);

// writes to subgraph_ostream the subgraph covered by
// the haplotypes in haplotype_list, as well as these haplotypes embedded as
// Paths.  Will output in JSON format if json set to true and Protobuf otherwise.
void output_graph_with_embedded_paths(ostream& subgraph_ostream,
            vector<pair<thread_t,int>>& haplotype_list, xg::XG& index, bool json = true);
// get the graph directly
Graph output_graph_with_embedded_paths(vector<pair<thread_t,int>>& haplotype_list, xg::XG& index);

// writes to annotation_ostream the list of counts of identical subhaplotypes
// using the same ordering as the Paths from output_graph_with_embedded_paths
void output_haplotype_counts(ostream& annotation_ostream,
            vector<pair<thread_t,int>>& haplotype_list, xg::XG& index);

// Adds to a Graph the nodes and edges touched by a thread_t
void thread_to_graph_spanned(thread_t& t, Graph& graph, xg::XG& index);
// Adds to a set of nodes all those touched by thread_t t
void add_thread_nodes_to_set(thread_t& t, set<int64_t>& nodes);
// Adds to a set of edges all those touched by thread_t t
void add_thread_edges_to_set(thread_t& t, set<pair<int,int> >& edges);
// Turns a set of nodes and a set of edges into a Graph
void construct_graph_from_nodes_and_edges(Graph& g, xg::XG& index,
            set<int64_t>& nodes, set<pair<int,int> >& edges);

#endif
