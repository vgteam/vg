#include "phase_duplicator.hpp"

namespace vg {

using namespace std;

PhaseDuplicator::PhaseDuplicator(const xg::XG& index) : index(index) {
    // Nothing to do!
}

pair<Graph, vector<Translation>> PhaseDuplicator::duplicate(set<id_t> subgraph, id_t& next_id) const {

    // Allocate space for the result
    pair<Graph, vector<Translation>> to_return;
    // Name the pair members
    Graph& duplicated = to_return.first;
    vector<Translation>& translations = to_return.second;
    
    // TODO: actually do the duplication
    
    // Ship out the result
    return to_return;
}

vector<pair<xg::thread_t, int>> PhaseDuplicator::list_haplotypes(const xg::XG& index, set<id_t> subgraph) const {
    // Adapted from Yohei's "haplotype_extracter"
    
    // This holds extracted sub-haplotypes and the search states for looking right off their ends.
    vector<pair<thread_t, xg::XG::ThreadSearchState>> search_intermediates;
    // This holds the results we will return, as pairs of extracted haplotypes and haplotype counts
    vector<pair<thread_t, int>> search_results;
    
    // Start a search at the first node traversal
    thread_t first_thread = {start_node};
    xg::XG::ThreadSearchState first_state;
    index.extend_search(first_state, first_thread);
    
    // Grab the edges out of that node traversal
    vector<Edge> edges = start_node.is_reverse ? index.edges_on_start(start_node.node_id) : index.edges_on_end(start_node.node_id);
    for (int i = 0; i < edges.size(); i++) {
        // Consider following each edge
        
        // Where would we go?
        xg::XG::ThreadMapping next_node;
        next_node.node_id = edges[i].to();
        next_node.is_reverse = edges[i].to_end();
        
        // How many haplotypes also go there?
        xg::XG::ThreadSearchState new_state = first_state;
        thread_t t = {next_node};
        index.extend_search(new_state, t);
        
        // If any do, remember there as a search state
        thread_t new_thread = first_thread;
        new_thread.push_back(next_node);
        if (!new_state.is_empty()) {
            search_intermediates.push_back(make_pair(new_thread,new_state));
        }
    }
    
    while (search_intermediates.size() > 0) {
        // Until we've finished all our intermediates...
        // Pop one off
        pair<thread_t,xg::XG::ThreadSearchState> last = search_intermediates.back();
        search_intermediates.pop_back();
        
        // Remember how many we had besides thsi one, so we can see if we added
        // any descendants of this one, or if this one is a dead end.
        int check_size = search_intermediates.size();
        
        // Grab the edges out of the last thing in the partial haplotype
        vector<Edge> edges = last.first.back().is_reverse ? index.edges_on_start(last.first.back().node_id) : index.edges_on_end(last.first.back().node_id);
        if (edges.size() == 0) {
            // If there's nowhere to go, we're a dead end, so finish the haplotype here.
            search_results.push_back(make_pair(last.first,last.second.count()));
        } else {
            // There are edges to trace
            for (int i = 0; i < edges.size(); i++) {
                // Follow each of them
                xg::XG::ThreadMapping next_node;
                next_node.node_id = edges[i].to();
                next_node.is_reverse = edges[i].to_end();
                
                // Try searching with where the edge goes
                xg::XG::ThreadSearchState new_state = last.second;
                thread_t next_thread = {next_node};
                index.extend_search(new_state,next_thread);
                
                // If there's anything there...
                thread_t new_thread = last.first;
                new_thread.push_back(next_node);
                if (!new_state.is_empty()) {
                    if (new_thread.size() >= extend_distance) {
                        // We found something but it's too far away. Stop here.
                        search_results.push_back(make_pair(new_thread,new_state.count()));
                    } else {
                        // We found something and it's close enough to explore. Go explore it.
                        search_intermediates.push_back(make_pair(new_thread,new_state));
                    }
                }
            }
            if (check_size == search_intermediates.size() && last.first.size() < extend_distance - 1) {
                // If we didn't hit a dead end or our extension limit but no
                // haplotypes actually went anywhere, then we need to count what
                // we have so far as a haplotype.
                search_results.push_back(make_pair(last.first,last.second.count()));
            }
        }
    }
    return search_results;
}

}
