#include <iostream>
#include "vg.hpp"
#include "haplotype_extracter.hpp"
#include "json2pb.h"
#include "xg.hpp"

using namespace std;
using namespace vg;

void trace_haplotypes_and_paths(xg::XG& index,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph) {
  // get our haplotypes
  xg::XG::ThreadMapping n = {start_node, false};
  vector<pair<thread_t,int> > haplotypes = list_haplotypes(index, n, extend_distance);

  if (expand_graph) {
    // get our subgraph and "regular" paths by expanding forward
    *out_graph.add_node() = index.node(start_node);
    index.expand_context(out_graph, extend_distance, true, true, true, false);
  }

  // add a frequency of 1 for each normal path
  for (int i = 0; i < out_graph.path_size(); ++i) {
    out_thread_frequencies[out_graph.path(i).name()] = 1;
  }

  // add our haplotypes to the subgraph, naming ith haplotype "thread_i"
  for (int i = 0; i < haplotypes.size(); ++i) {
    Path p = path_from_thread_t(haplotypes[i].first);
    p.set_name("thread_" + to_string(i));
    *(out_graph.add_path()) = move(p);
    out_thread_frequencies[p.name()] = haplotypes[i].second;
  }
}


void output_haplotype_counts(ostream& annotation_ostream,
            vector<pair<thread_t,int>>& haplotype_list, xg::XG& index) {
  for(int i = 0; i < haplotype_list.size(); i++) {
    annotation_ostream << i << "\t" << haplotype_list[i].second << endl;
  }
}

Graph output_graph_with_embedded_paths(vector<pair<thread_t,int>>& haplotype_list, xg::XG& index) {
  Graph g;
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  for(int i = 0; i < haplotype_list.size(); i++) {
    add_thread_nodes_to_set(haplotype_list[i].first, nodes);
    add_thread_edges_to_set(haplotype_list[i].first, edges);
  }
  construct_graph_from_nodes_and_edges(g, index, nodes, edges);
  for(int i = 0; i < haplotype_list.size(); i++) {
    Path p = path_from_thread_t(haplotype_list[i].first);
    p.set_name(to_string(i));
    *(g.add_path()) = move(p);
  }
  return g;
}
 
void output_graph_with_embedded_paths(ostream& subgraph_ostream,
            vector<pair<thread_t,int>>& haplotype_list, xg::XG& index, bool json) {
  Graph g = output_graph_with_embedded_paths(haplotype_list, index);

  if (json) {
    subgraph_ostream << pb2json(g);
  } else {
    VG subgraph;
    subgraph.extend(g);
    subgraph.serialize_to_ostream(subgraph_ostream);
  }
}

void thread_to_graph_spanned(thread_t& t, Graph& g, xg::XG& index) {
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  nodes.insert(t[0].node_id);
  for(int i = 1; i < t.size(); i++) {
    nodes.insert(t[i].node_id);
    edges.insert(make_pair(xg::make_side(t[i-1].node_id,t[i-1].is_reverse), 
              xg::make_side(t[i].node_id,t[i].is_reverse)));
	}
 	for (auto& n : nodes) {
    *g.add_node() = index.node(n);
 	}
  for (auto& e : edges) {
    Edge edge;
    edge.set_from(xg::side_id(e.first));
    edge.set_from_start(xg::side_is_end(e.first));
    edge.set_to(xg::side_id(e.second));
    edge.set_to_end(xg::side_is_end(e.second));
    *g.add_edge() = edge;
  }
}

void add_thread_nodes_to_set(thread_t& t, set<int64_t>& nodes) {
  for(int i = 0; i < t.size(); i++) {
    nodes.insert(t[i].node_id);
  }
}

void add_thread_edges_to_set(thread_t& t, set<pair<int,int> >& edges) {
  for(int i = 1; i < t.size(); i++) {
    edges.insert(make_pair(xg::make_side(t[i-1].node_id,t[i-1].is_reverse), 
              xg::make_side(t[i].node_id,t[i].is_reverse)));
  }
}

void construct_graph_from_nodes_and_edges(Graph& g, xg::XG& index,
            set<int64_t>& nodes, set<pair<int,int> >& edges) {
  for (auto& n : nodes) {
	   *g.add_node() = index.node(n);
 	}
  for (auto& e : edges) {
    Edge edge;
    edge.set_from(xg::side_id(e.first));
    edge.set_from_start(xg::side_is_end(e.first));
    edge.set_to(xg::side_id(e.second));
    edge.set_to_end(xg::side_is_end(e.second));
    *g.add_edge() = edge;
  }
}

Path path_from_thread_t(thread_t& t) {
	Path toReturn;
	int rank = 1;
	for(int i = 0; i < t.size(); i++) {
		Mapping* mapping = toReturn.add_mapping();

    // Set up the position
    mapping->mutable_position()->set_node_id(t[i].node_id);
    mapping->mutable_position()->set_is_reverse(t[i].is_reverse);

    // Set the rank
    mapping->set_rank(rank++);
  }
  // We're done making the path
  return toReturn;
}

vector<pair<thread_t,int> > list_haplotypes(xg::XG& index,
            xg::XG::ThreadMapping start_node, int extend_distance) {
  vector<pair<thread_t,xg::XG::ThreadSearchState> > search_intermediates;
  vector<pair<thread_t,int> > search_results;
  thread_t first_thread = {start_node};
  xg::XG::ThreadSearchState first_state;
  index.extend_search(first_state,first_thread);
  vector<Edge> edges = start_node.is_reverse ?
            index.edges_on_start(start_node.node_id) :
            index.edges_on_end(start_node.node_id);
  for(int i = 0; i < edges.size(); i++) {
    xg::XG::ThreadMapping next_node;
    next_node.node_id = edges[i].to();
    next_node.is_reverse = edges[i].to_end();
    xg::XG::ThreadSearchState new_state = first_state;
    thread_t t = {next_node};
    index.extend_search(new_state, t);
    thread_t new_thread = first_thread;
    new_thread.push_back(next_node);
    if(!new_state.is_empty()) {
      search_intermediates.push_back(make_pair(new_thread,new_state));
    }
  }
  while(search_intermediates.size() > 0) {
    pair<thread_t,xg::XG::ThreadSearchState> last = search_intermediates.back();
    search_intermediates.pop_back();
    int check_size = search_intermediates.size();
    vector<Edge> edges = last.first.back().is_reverse ?
              index.edges_on_start(last.first.back().node_id) :
              index.edges_on_end(last.first.back().node_id);
    if(edges.size() == 0) {
      search_results.push_back(make_pair(last.first,last.second.count()));
    } else {
      for(int i = 0; i < edges.size(); i++) {
        xg::XG::ThreadMapping next_node;
        next_node.node_id = edges[i].to();
        next_node.is_reverse = edges[i].to_end();
        xg::XG::ThreadSearchState new_state = last.second;
        thread_t next_thread = {next_node};
        index.extend_search(new_state,next_thread);
        thread_t new_thread = last.first;
        new_thread.push_back(next_node);
        if(!new_state.is_empty()) {
          if(new_thread.size() >= extend_distance) {
            search_results.push_back(make_pair(new_thread,new_state.count()));
          } else {
            search_intermediates.push_back(make_pair(new_thread,new_state));
          }
        }
      }
      if(check_size == search_intermediates.size() &&
                last.first.size() < extend_distance - 1) {
        search_results.push_back(make_pair(last.first,last.second.count()));
      }
    }
  }
  return search_results;
}
