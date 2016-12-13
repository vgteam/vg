#ifndef HAPLOTYPE_EXTRACTER_H
#define HAPLOTYPE_EXTRACTER_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"
#include "haplotypes.hpp"

using namespace std;
using namespace vg;
using thread_t = vector<xg::XG::ThreadMapping>;

void output_graph_with_embedded_paths(string output_path, vector<pair<thread_t,int>>& haplotype_list, xg::XG& index);
void output_weighted_haplotype_list(string output_path, vector<pair<thread_t,int>>& haplotype_list, xg::XG& index, bool likelihoods);
Path path_from_thread_t(thread_t& t);
vector<pair<thread_t,int> > list_haplotypes(xg::XG& index, xg::XG::ThreadMapping start_node, int extend_distance);

void thread_to_graph_spanned(thread_t& t, Graph& graph, xg::XG& index);
void add_thread_nodes_to_set(thread_t& t, set<int64_t>& nodes);
void add_thread_edges_to_set(thread_t& t, set<pair<int,int> >& edges);
void construct_graph_from_nodes_and_edges(Graph& g, xg::XG& index, set<int64_t>& nodes, set<pair<int,int> >& edges);

#endif
