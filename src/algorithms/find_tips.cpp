#include "find_tips.hpp"

namespace vg {
namespace algorithms {

using namespace std;

vector<handle_t> head_nodes(const HandleGraph* g) {
    vector<handle_t> to_return;
    g->for_each_handle([&](const handle_t& found) {
        // For each (locally forward) node
        
        bool no_left_edges = true;
        g->follow_edges(found, true, [&](const handle_t& ignored) {
            // We found a left edge!
            no_left_edges = false;
            // We only need one
            return false;
        });
        
        if (no_left_edges) {
            to_return.push_back(found);
        }
    });
    
    return to_return;
    
}

vector<handle_t> tail_nodes(const HandleGraph* g) {
    vector<handle_t> to_return;
    g->for_each_handle([&](const handle_t& found) {
        // For each (locally forward) node
        
        bool no_right_edges = true;
        g->follow_edges(found, false, [&](const handle_t& ignored) {
            // We found a right edge!
            no_right_edges = false;
            // We only need one
            return false;
        });
        
        if (no_right_edges) {
            to_return.push_back(found);
        }
    });
    
    return to_return;
    
}

vector<handle_t> find_tips(const HandleGraph* g) {
    // Start with the heads
    vector<handle_t> tips = head_nodes(g);
    vector<handle_t> tails = tail_nodes(g);
    tips.reserve(tips.size() + tails.size());
    for (auto tip : tails) {
        // And add all the tails backward
        tips.push_back(g->flip(tip));
    }
    return tips;
}

}
}
