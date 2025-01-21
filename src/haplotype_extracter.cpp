#include <iostream>
#include "vg.hpp"
#include "haplotype_extracter.hpp"
#include "vg/io/json2pb.h"
#include "graph.hpp"

namespace vg {

using namespace std;

static id_t side_id(int64_t side) {
    return abs(side);
}

static bool side_is_end(int64_t side) {
    return side < 0;
}

static int64_t make_side(id_t id, bool is_end) {
    return !is_end ? id : -1 * id;
}


void trace_haplotypes(const PathHandleGraph& source, const gbwt::GBWT& haplotype_database,
                      const handle_t& start_node, function<bool(const vector<gbwt::node_type>&)> stop_fn,
                      MutablePathMutableHandleGraph& out_graph,
                      map<string, int>& out_thread_frequencies) {
  // get our haplotypes
  // TODO: Tell it to keep haplotypes that end before triggering stop_fn!
  vector<pair<thread_t, gbwt::SearchState> > haplotypes = list_haplotypes(source, haplotype_database, start_node, stop_fn);

#ifdef debug
  cerr << "Haplotype database " << &haplotype_database << " produced " << haplotypes.size() << " haplotypes" << endl;
#endif

  // add our haplotypes to the subgraph, naming ith haplotype "thread_i"
  for (int i = 0; i < haplotypes.size(); ++i) {
    std::string path_name = "thread_" + to_string(i);
    path_handle_t path_handle = out_graph.create_path_handle(path_name);
    const thread_t& thread = haplotypes[i].first;
    for (int j = 0; k < thread.size(); j++) {
      // Copy in all the steps
      out_graph.append_step(path_handle, out_graph.get_handle(gbwt::Node::id(thread[j]), gbwt::Node::is_reverse(thread[j])));
    }
  }
}

void trace_paths(const PathHandleGraph& source,
                 const handle_t& start_node, int extend_distance,
                 MutablePathMutableHandleGraph& out_graph,
                 map<string, int>& out_thread_frequencies) {

  // get our subgraph and "regular" paths by expanding forward
  bdsg::HashGraph extractor;
  extractor.create_handle(source.get_sequence(start_node), source.get_id(start_node));
  // TODO: is expanding only forward really the right behavior here?
  algorithms::expand_context_with_paths(&source, &extractor, extend_distance, true, true, false);

  // Copy over nodes and edges
  copy_handle_graph(&extractor, &out_graph);

  for (auto& sense : {PathSense::REFERENCE, PathSense::GENERIC, PathSense::HAPLOTYPE}) {
    extractor.for_each_path_of_sense(sense, [&](const path_handle_t& path_handle) {
        // For every non-haplotype path we pulled, give it a frequency of 1.
        // Including haplotypes, since we didn't trace them to get frequencies.
        out_thread_frequencies[extractor.get_path_name(path_handle)] = 1;

        // Copy it over
        copy_path(&extractor, path_handle, &out_graph);
    });
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

vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > list_haplotypes(const HandleGraph& graph,
                                                                          const gbwt::GBWT& gbwt,
                                                                          handle_t start,
                                                                          function<bool(const vector<gbwt::node_type>&)> stop_fn) {
    
    // Keep track of all the different paths we're extending
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > search_intermediates;
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > search_results;

    // Look up the start node in GBWT and start a thread
    gbwt::node_type start_node = handle_to_gbwt(graph, start);    
    vector<gbwt::node_type> first_thread = {start_node};
    gbwt::SearchState first_state = gbwt.find(start_node);
    
#ifdef debug
    cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(start_node)  << ":"
         << gbwt::Node::is_reverse(start_node) << endl;
#endif

    if (!first_state.empty()) {
        search_intermediates.push_back(make_pair(first_thread, first_state));
    }

    while(!search_intermediates.empty()) {

        // pick up a thread to continue from the queue
        auto last = std::move(search_intermediates.back());
        search_intermediates.pop_back();

        vector<tuple<handle_t, gbwt::node_type, gbwt::SearchState>> next_handle_states;
        graph.follow_edges(gbwt_to_handle(graph, last.first.back()), false, [&](const handle_t& next) {
                // extend the last node of the thread using gbwt
                auto extend_node = handle_to_gbwt(graph, next);
                auto new_state = gbwt.extend(last.second, extend_node);
#ifdef debug
                cerr << "Extend state " << last.second << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
                if (!new_state.empty()) {
                    next_handle_states.push_back(make_tuple(next, extend_node, new_state));
                }                    
            });

        for (auto& nhs : next_handle_states) {
            
            const handle_t& next = get<0>(nhs);
            gbwt::node_type& extend_node = get<1>(nhs);
            gbwt::SearchState& new_state = get<2>(nhs);
                
            vector<gbwt::node_type> new_thread;
            if (&nhs == &next_handle_states.back()) {
                // avoid a copy by re-using the vector for the last thread. this way simple cases
                // like scanning along one path don't blow up to n^2
                new_thread = std::move(last.first);
            } else {
                new_thread = last.first;
            }                        
            new_thread.push_back(extend_node);

            if (stop_fn(new_thread)) {
#ifdef debug
                cerr << "\tGot " << new_state.size() << " results at limit; emitting" << endl;
#endif
                search_results.push_back(make_pair(std::move(new_thread), new_state));
            }
            else {
#ifdef debug
                cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
                search_intermediates.push_back(make_pair(std::move(new_thread), new_state));
            }
        }

        if (next_handle_states.empty()) {
            // TODO: Also handle case where some but not all haplotypes go on, and some end here.
#ifdef debug
            std::cerr <<  "Extend state " << last.second << " to nowhere; haplotype ends!" << std::endl;
#endif
        }

    }
    
    return search_results;
}



}
