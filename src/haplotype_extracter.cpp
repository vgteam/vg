#include <iostream>
#include "vg.hpp"
#include "haplotype_extracter.hpp"
#include "json2pb.h"
#include "xg.hpp"

namespace vg {

using namespace std;

void trace_haplotypes_and_paths(vg::xg::XG& index, const gbwt::GBWT* haplotype_database,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph) {
  // get our haplotypes
  xg::XG::ThreadMapping n = {start_node, false};
  vector<pair<thread_t,int> > haplotypes = haplotype_database ?
    list_haplotypes(index, *haplotype_database, n, extend_distance) :
    list_haplotypes(index, n, extend_distance);

#ifdef debug
  cerr << "Haplotype database " << haplotype_database << " produced " << haplotypes.size() << " haplotypes" << endl;
#endif

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
    Path p = path_from_thread_t(haplotypes[i].first, index);
    p.set_name("thread_" + to_string(i));
    out_thread_frequencies[p.name()] = haplotypes[i].second;
    *(out_graph.add_path()) = move(p);
  }
}


void output_haplotype_counts(ostream& annotation_ostream,
            vector<pair<thread_t,int>>& haplotype_list, vg::xg::XG& index) {
  for(int i = 0; i < haplotype_list.size(); i++) {
    annotation_ostream << i << "\t" << haplotype_list[i].second << endl;
  }
}

Graph output_graph_with_embedded_paths(vector<pair<thread_t,int>>& haplotype_list, vg::xg::XG& index) {
  Graph g;
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  for(int i = 0; i < haplotype_list.size(); i++) {
    add_thread_nodes_to_set(haplotype_list[i].first, nodes);
    add_thread_edges_to_set(haplotype_list[i].first, edges);
  }
  construct_graph_from_nodes_and_edges(g, index, nodes, edges);
  for(int i = 0; i < haplotype_list.size(); i++) {
    Path p = path_from_thread_t(haplotype_list[i].first, index);
    p.set_name(to_string(i));
    *(g.add_path()) = move(p);
  }
  return g;
}
 
void output_graph_with_embedded_paths(ostream& subgraph_ostream,
            vector<pair<thread_t,int>>& haplotype_list, vg::xg::XG& index, bool json) {
  Graph g = output_graph_with_embedded_paths(haplotype_list, index);

  if (json) {
    subgraph_ostream << pb2json(g);
  } else {
    VG subgraph;
    subgraph.extend(g);
    subgraph.serialize_to_ostream(subgraph_ostream);
  }
}

void thread_to_graph_spanned(thread_t& t, Graph& g, vg::xg::XG& index) {
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  nodes.insert(t[0].node_id);
  for(int i = 1; i < t.size(); i++) {
    nodes.insert(t[i].node_id);
    edges.insert(make_pair(make_side(t[i-1].node_id,t[i-1].is_reverse),
              make_side(t[i].node_id,t[i].is_reverse)));
	}
 	for (auto& n : nodes) {
    *g.add_node() = index.node(n);
 	}
  for (auto& e : edges) {
    Edge edge;
    edge.set_from(side_id(e.first));
    edge.set_from_start(side_is_end(e.first));
    edge.set_to(side_id(e.second));
    edge.set_to_end(side_is_end(e.second));
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
    edges.insert(make_pair(make_side(t[i-1].node_id,t[i-1].is_reverse),
              make_side(t[i].node_id,t[i].is_reverse)));
  }
}

void construct_graph_from_nodes_and_edges(Graph& g, vg::xg::XG& index,
            set<int64_t>& nodes, set<pair<int,int> >& edges) {
  for (auto& n : nodes) {
	   *g.add_node() = index.node(n);
 	}
  for (auto& e : edges) {
    Edge edge;
    edge.set_from(side_id(e.first));
    edge.set_from_start(side_is_end(e.first));
    edge.set_to(side_id(e.second));
    edge.set_to_end(side_is_end(e.second));
    *g.add_edge() = edge;
  }
}

Path path_from_thread_t(thread_t& t, vg::xg::XG& index) {
	Path toReturn;
	int rank = 1;
	for(int i = 0; i < t.size(); i++) {
		Mapping* mapping = toReturn.add_mapping();

        // Set up the position
        mapping->mutable_position()->set_node_id(t[i].node_id);
        mapping->mutable_position()->set_is_reverse(t[i].is_reverse);

        // Set up the edits
        Edit* e = mapping->add_edit();
        size_t l = index.node_length(t[i].node_id);
        e->set_from_length(l);
        e->set_to_length(l);

        // Set the rank
        mapping->set_rank(rank++);
    }
    // We're done making the path
    return toReturn;
}

vector<pair<thread_t,int> > list_haplotypes(vg::xg::XG& index,
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

vector<pair<thread_t,int> > list_haplotypes(vg::xg::XG& index, const gbwt::GBWT& haplotype_database,
            xg::XG::ThreadMapping start_node, int extend_distance) {

#ifdef debug
  cerr << "Extracting haplotypes from GBWT" << endl;
#endif

  vector<pair<thread_t,gbwt::SearchState> > search_intermediates;
  vector<pair<thread_t,int> > search_results;
  // We still keep our data as thread_ts full of xg ThreadMappings and convert on the fly.
  thread_t first_thread = {start_node};
  auto first_node = gbwt::Node::encode(start_node.node_id, start_node.is_reverse);
  gbwt::SearchState first_state = haplotype_database.find(first_node);
#ifdef debug
  cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(first_node) << endl;
#endif
  vector<Edge> edges = start_node.is_reverse ?
            index.edges_on_start(start_node.node_id) :
            index.edges_on_end(start_node.node_id);

  // TODO: this is just most of the loop body repeated!
  for(int i = 0; i < edges.size(); i++) {
    xg::XG::ThreadMapping next_node;
    next_node.node_id = edges[i].to();
    next_node.is_reverse = edges[i].to_end();
    auto extend_node = gbwt::Node::encode(next_node.node_id, next_node.is_reverse);
    auto new_state = haplotype_database.extend(first_state, extend_node);
#ifdef debug
    cerr << "Extend state " << first_state << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
    thread_t new_thread = first_thread;
    new_thread.push_back(next_node);
    if(!new_state.empty()) {
#ifdef debug
      cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
      search_intermediates.push_back(make_pair(new_thread,new_state));
    }
  }
  while(search_intermediates.size() > 0) {
    pair<thread_t,gbwt::SearchState> last = search_intermediates.back();
    search_intermediates.pop_back();
    int check_size = search_intermediates.size();
    vector<Edge> edges = last.first.back().is_reverse ?
              index.edges_on_start(last.first.back().node_id) :
              index.edges_on_end(last.first.back().node_id);
    if(edges.size() == 0) {
#ifdef debug
      cerr << "Hit end of graph on state " << last.second << endl;
#endif
      search_results.push_back(make_pair(last.first,last.second.size()));
    } else {
      for(int i = 0; i < edges.size(); i++) {
        xg::XG::ThreadMapping next_node;
        next_node.node_id = edges[i].to();
        next_node.is_reverse = edges[i].to_end();
        auto extend_node = gbwt::Node::encode(next_node.node_id, next_node.is_reverse);
        auto new_state = haplotype_database.extend(last.second, extend_node);
#ifdef debug
        cerr << "Extend state " << last.second << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
        thread_t new_thread = last.first;
        new_thread.push_back(next_node);
        if(!new_state.empty()) {
          if(new_thread.size() >= extend_distance) {
#ifdef debug
            cerr << "\tGot " << new_state.size() << " results at limit; emitting" << endl;
#endif
            search_results.push_back(make_pair(new_thread,new_state.size()));
          } else {
#ifdef debug
            cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
            search_intermediates.push_back(make_pair(new_thread,new_state));
          }
        }
      }
      if(check_size == search_intermediates.size() &&
                last.first.size() < extend_distance - 1) {
        search_results.push_back(make_pair(last.first,last.second.size()));
      }
    }
  }
  return search_results;
}

}
