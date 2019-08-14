#include <iostream>
#include "vg.hpp"
#include "haplotype_extracter.hpp"
#include "json2pb.h"
#include "graph.hpp"

namespace vg {

using namespace std;

void trace_haplotypes_and_paths(const PathHandleGraph& source, const gbwt::GBWT& haplotype_database,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph) {
  // get our haplotypes
  gbwt::node_type n = gbwt::Node::encode(start_node, false);
  vector<pair<thread_t,int> > haplotypes = list_haplotypes(source, haplotype_database, n, extend_distance);

#ifdef debug
  cerr << "Haplotype database " << &haplotype_database << " produced " << haplotypes.size() << " haplotypes" << endl;
#endif

  if (expand_graph) {
      // get our subgraph and "regular" paths by expanding forward
      handle_t handle = source.get_handle(start_node);
      bdsg::HashGraph extractor;
      extractor.create_handle(source.get_sequence(handle), source.get_id(handle));
      // TODO: is expanding only forward really the right behavior here?
      algorithms::expand_context_with_paths(&source, &extractor, extend_distance, true, true, false);
      
      // Convert to Protobuf I guess
      from_path_handle_graph(extractor, out_graph);
  }

  // add a frequency of 1 for each normal path
  for (int i = 0; i < out_graph.path_size(); ++i) {
    out_thread_frequencies[out_graph.path(i).name()] = 1;
  }

  // add our haplotypes to the subgraph, naming ith haplotype "thread_i"
  for (int i = 0; i < haplotypes.size(); ++i) {
    Path p = path_from_thread_t(haplotypes[i].first, source);
    p.set_name("thread_" + to_string(i));
    out_thread_frequencies[p.name()] = haplotypes[i].second;
    *(out_graph.add_path()) = move(p);
  }
}


void output_haplotype_counts(ostream& annotation_ostream,
            vector<pair<thread_t,int>>& haplotype_list) {
  for(int i = 0; i < haplotype_list.size(); i++) {
    annotation_ostream << i << "\t" << haplotype_list[i].second << endl;
  }
}

Graph output_graph_with_embedded_paths(vector<pair<thread_t,int>>& haplotype_list, const HandleGraph& source) {
  Graph g;
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  for(int i = 0; i < haplotype_list.size(); i++) {
    add_thread_nodes_to_set(haplotype_list[i].first, nodes);
    add_thread_edges_to_set(haplotype_list[i].first, edges);
  }
  construct_graph_from_nodes_and_edges(g, source, nodes, edges);
  for(int i = 0; i < haplotype_list.size(); i++) {
    Path p = path_from_thread_t(haplotype_list[i].first, source);
    p.set_name(to_string(i));
    *(g.add_path()) = move(p);
  }
  return g;
}
 
void output_graph_with_embedded_paths(ostream& subgraph_ostream,
            vector<pair<thread_t,int>>& haplotype_list, const HandleGraph& source, bool json) {
  Graph g = output_graph_with_embedded_paths(haplotype_list, source);

  if (json) {
    subgraph_ostream << pb2json(g);
  } else {
    VG subgraph;
    subgraph.extend(g);
    subgraph.serialize_to_ostream(subgraph_ostream);
  }
}

void thread_to_graph_spanned(thread_t& t, Graph& g, const HandleGraph& source) {
  set<int64_t> nodes;
  set<pair<int,int> > edges;
  nodes.insert(gbwt::Node::id(t[0]));
  for(int i = 1; i < t.size(); i++) {
    nodes.insert(gbwt::Node::id(t[i]));
    edges.insert(make_pair(make_side(gbwt::Node::id(t[i-1]),gbwt::Node::is_reverse(t[i-1])),
              make_side(gbwt::Node::id(t[i]),gbwt::Node::is_reverse(t[i]))));
  }
  for (auto& n : nodes) {
    handle_t handle = source.get_handle(n);
    Node* node = g.add_node();
    node->set_sequence(source.get_sequence(handle));
    node->set_id(n);
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
    nodes.insert(gbwt::Node::id(t[i]));
  }
}

void add_thread_edges_to_set(thread_t& t, set<pair<int,int> >& edges) {
  for(int i = 1; i < t.size(); i++) {
    edges.insert(make_pair(make_side(gbwt::Node::id(t[i-1]),gbwt::Node::is_reverse(t[i-1])),
              make_side(gbwt::Node::id(t[i]),gbwt::Node::is_reverse(t[i]))));
  }
}

void construct_graph_from_nodes_and_edges(Graph& g, const HandleGraph& source,
            set<int64_t>& nodes, set<pair<int,int> >& edges) {
  for (auto& n : nodes) {
    handle_t handle = source.get_handle(n);
    Node* node = g.add_node();
    node->set_sequence(source.get_sequence(handle));
    node->set_id(n);
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

Path path_from_thread_t(thread_t& t, const HandleGraph& source) {
	Path toReturn;
	int rank = 1;
	for(int i = 0; i < t.size(); i++) {
		Mapping* mapping = toReturn.add_mapping();

        // Set up the position
        mapping->mutable_position()->set_node_id(gbwt::Node::id(t[i]));
        mapping->mutable_position()->set_is_reverse(gbwt::Node::is_reverse(t[i]));

        // Set up the edits
        Edit* e = mapping->add_edit();
        size_t l = source.get_length(source.get_handle(gbwt::Node::id(t[i])));
        e->set_from_length(l);
        e->set_to_length(l);

        // Set the rank
        mapping->set_rank(rank++);
    }
    // We're done making the path
    return toReturn;
}

vector<pair<thread_t,int> > list_haplotypes(const HandleGraph& source, const gbwt::GBWT& haplotype_database,
            gbwt::node_type start_node, int extend_distance) {

#ifdef debug
    cerr << "Extracting haplotypes from GBWT" << endl;
#endif

    vector<pair<thread_t,gbwt::SearchState> > search_intermediates;
    vector<pair<thread_t,int> > search_results;
    thread_t first_thread = {start_node};
    gbwt::SearchState first_state = haplotype_database.find(start_node);
#ifdef debug
    cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(start_node) << endl;
#endif
    
    source.follow_edges(gbwt_to_handle(source, start_node), false, [&](const handle_t& next) {
        auto extend_node = handle_to_gbwt(source, next);
        auto new_state = haplotype_database.extend(first_state, extend_node);
#ifdef debug
        cerr << "Extend state " << first_state << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
        thread_t new_thread = first_thread;
        new_thread.push_back(extend_node);
        if(!new_state.empty()) {
#ifdef debug
            cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
            search_intermediates.push_back(make_pair(new_thread,new_state));
        }
    });
    
    while(search_intermediates.size() > 0) {
        auto last = search_intermediates.back();
        search_intermediates.pop_back();
        
        int check_size = search_intermediates.size();
        bool any_edges = false;
        source.follow_edges(gbwt_to_handle(source, last.first.back()), false, [&](const handle_t& next) {
            any_edges = true;
            auto extend_node = handle_to_gbwt(source, next);
            auto new_state = haplotype_database.extend(last.second, extend_node);
#ifdef debug
            cerr << "Extend state " << last.second << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
            thread_t new_thread = last.first;
            new_thread.push_back(extend_node);
            if(!new_state.empty()) {
                if(new_thread.size() >= extend_distance) {
#ifdef debug
                    cerr << "\tGot " << new_state.size() << " results at limit; emitting" << endl;
#endif
                    search_results.push_back(make_pair(new_thread,new_state.size()));
                }
                else {
#ifdef debug
                    cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
                    search_intermediates.push_back(make_pair(new_thread,new_state));
                }
            }
        });
        if (!any_edges) {
#ifdef debug
            cerr << "Hit end of graph on state " << last.second << endl;
#endif
            search_results.push_back(make_pair(last.first,last.second.size()));
        }
        else if (check_size == search_intermediates.size() &&
                 last.first.size() < extend_distance - 1) {
            search_results.push_back(make_pair(last.first,last.second.size()));
        }
    }
    
    return search_results;
}

}
