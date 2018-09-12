#include "strongly_connected_components.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    // recursion-free version of Tarjan's strongly connected components algorithm
    // https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    // Generalized to bidirected graphs as described (confusingly) in
    // "Decomposition of a bidirected graph into strongly connected components and
    // its signed poset structure", by Kazutoshi Ando, Satoru Fujishige, and Toshio
    // Nemoto. http://www.sciencedirect.com/science/article/pii/0166218X95000683
    
    // The best way to think about that paper is that the edges are vectors in a
    // vector space with number of dimensions equal to the number of nodes in the
    // graph, and an edge attaching to the end a node is the positive unit vector in
    // its dimension, and an edge attaching to the start of node is the negative
    // unit vector in its dimension.
    
    // The basic idea is that you just consider the orientations as different nodes,
    // and the edges as existing between both pairs of orientations they connect,
    // and do connected components on that graph. Since we don't care about
    // "consistent" or "inconsistent" strongly connected components, we just put a
    // node in a component if either orientation is in it. But bear in mind that
    // both orientations of a node might not actually be in the same strongly
    // connected component in a bidirected graph, so now the components may overlap.
    vector<unordered_set<id_t>> strongly_connected_components(const HandleGraph* handle_graph) {
        
#ifdef debug
        cerr << "Computing strongly connected components" << endl;
#endif
        
        // What node visit step are we on?
        int64_t index = 0;
        // What's the search root from which a node was reached?
        unordered_map<handle_t, handle_t> roots;
        // At what index step was each node discovered?
        unordered_map<handle_t, int64_t> discover_idx;
        // We need our own copy of the DFS stack
        vector<handle_t> stack;
        // And our own set of nodes already on the stack
        unordered_set<handle_t> on_stack;
        // What components did we find? Because of the way strongly connected
        // components generalizes, both orientations of a node always end up in the
        // same component.
        vector<unordered_set<id_t>> components;
        
        dfs(*handle_graph,
        [&](const handle_t& trav) {
            // When a NodeTraversal is first visited
#ifdef debug
            cerr << "First visit to " << handle_graph->get_id(trav) << " orientation " << handle_graph->get_is_reverse(trav) << endl;
#endif
            // It is its own root
            roots[trav] = trav;
            // We discovered it at this step
            discover_idx[trav] = index++;
            // And it's on the stack
            stack.push_back(trav);
            on_stack.insert(trav);
        },
        [&](const handle_t& trav) {
            // When a NodeTraversal is done being recursed into
#ifdef debug
            cerr << "Finishing " << handle_graph->get_id(trav) << " orientation " << handle_graph->get_is_reverse(trav) << endl;
#endif
            // Go through all the NodeTraversals reachable reading onwards from this traversal.
            handle_graph->follow_edges(trav, false, [&](const handle_t& next) {
#ifdef debug
                cerr << "\tCould next reach " << handle_graph->get_id(next) << " orientation " << handle_graph->get_is_reverse(next) << endl;
#endif
                if (on_stack.count(next)) {
                    // If any of those NodeTraversals are on the stack already
#ifdef debug
                    cerr << "\t\tIt is already on the stack, so maybe we want its root" << endl;
#endif
                    auto& node_root = roots[trav];
                    auto& next_root = roots[next];
#ifdef debug
                    cerr << "\t\t\tWe have root " << handle_graph->get_id(node_root) << " orientation "
                        << handle_graph->get_is_reverse(node_root)
                        << " discovered at time " << discover_idx[node_root] << endl;
                    cerr << "\t\t\tThey have root " << handle_graph->get_id(next_root) << " orientation "
                        << handle_graph->get_is_reverse(next_root)
                        << " discovered at time " << discover_idx[next_root] << endl;
#endif
                    // Adopt the root of the NodeTraversal that was discovered first.
                    roots[trav] = discover_idx[node_root] < discover_idx[next_root] ?
                    node_root : next_root;
#ifdef debug
                    cerr << "\t\t\tWinning root: " << handle_graph->get_id(roots[trav]) << " orientation "
                        << handle_graph->get_is_reverse(roots[trav]) << endl;
#endif
                }
                return true;
            });
            
            if (roots[trav] == trav) {
                // If we didn't find a better root
#ifdef debug
                cerr << "\tWe are our own best root, so glom up everything under us" << endl;
#endif
                handle_t other;
                components.emplace_back();
                auto& component = components.back();
                do
                {
                    // Grab everything that was put on the DFS stack below us
                    // and put it in our component.
                    other = stack.back();
                    stack.pop_back();
                    on_stack.erase(other);
                    component.insert(handle_graph->get_id(other));
#ifdef debug
                    cerr << "\t\tSnarf up node " << handle_graph->get_id(other) << " from handle in orientation "
                        << handle_graph->get_is_reverse(other) << endl;
#endif
                } while (other != trav);
            }
        },
        vector<handle_t>(), unordered_set<handle_t>());
        
        return components;
    }

}
}
