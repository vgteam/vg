///
///  \file integrated_snarl_finder.cpp
///
///

#include "integrated_snarl_finder.hpp"

#include <bdsg/overlays/overlay_helper.hpp>
#include <structures/union_find.hpp>

namespace vg {

using namespace std;

class IntegratedSnarlFinder::MergedAdjacencyGraph {
protected:
    /// Hold onto the backing HandleGraph
    const HandleGraph* graph;
    
    /// Keep a vectorizable overlay over it to let us map between handles
    /// and union-find indices via handle ranking. The handles are all at index
    /// (rank - 1) * 2 + is_reverse.
    ///
    /// We rely on handles in the vectorizable overlay and handles in the
    /// backing graph being identical.
    ///
    /// TODO: make not-mutable when https://github.com/vgteam/libbdsg/issues/63
    /// is fixed.
    mutable bdsg::VectorizableOverlayHelper overlay_helper;
    
    /// Keep a union-find over the ranks of the merged oriented handles that
    /// make up each component. Runs with include_children=true so we can find
    /// all the members of each group.
    ///
    /// Needs to be mutable because union-find find operations do internal tree
    /// massaging and aren't const.
    /// TODO: this makes read operations not thread safe!
    mutable structures::UnionFind union_find;
    
    /// Get the rank corresponding to the given handle, in the union-find.
    /// Our ranks are 0-based.
    size_t uf_rank(handle_t into) const;
    
    /// Get the handle with the given rank in union-find space.
    /// Our ranks are 0-based.
    handle_t uf_handle(size_t rank) const;
    
public:
    /// Make a MergedAdjacencyGraph representing the graph of adjacency components of the given HandleGraph.
    MergedAdjacencyGraph(const HandleGraph* graph);
    
    /// Copy a MergedAdjacencyGraph by re-doing all the merges. Uses its own internal vectorization.
    MergedAdjacencyGraph(const MergedAdjacencyGraph& other);
    
    /// Given handles reading into two components, a and b, merge them into a single component.
    void merge(handle_t into_a, handle_t into_b);
    
    /// Find the handle heading the component that the given handle is in.
    handle_t find(handle_t into) const; 
    
    /// For each item other than the head in each component, calls the iteratee
    /// with the head and the other item. Does not call the iteratee for
    /// single-item components.
    void for_each_membership(const function<void(handle_t, handle_t)>& iteratee) const;
};

size_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_rank(handle_t into) const {
    // We need to 0-base the backing rank, space it out, and make the low bit orientation
    return (overlay_helper.get()->id_to_rank(graph->get_id(into)) - 1) * 2 + (size_t) graph->get_is_reverse(into);
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_handle(size_t rank) const {
    // We need to take the high bits and than make it 1-based, and get the orientation from the low bit
    return graph->get_handle(overlay_helper.get()->rank_to_id(rank / 2 + 1), rank % 2);
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const HandleGraph* graph) : graph(graph),
    overlay_helper(), union_find(graph->get_node_count() * 2, true) {
    
    // Make sure we have our vectorizable version of the graph.
    // TODO: remove const_cast when https://github.com/vgteam/libbdsg/issues/64 is fixed.
    overlay_helper.apply(const_cast<HandleGraph*>(graph));
    
    // TODO: we want the adjacency components that are just single edges
    // between two handles (i.e. trivial snarls) to be implicit, so we don't
    // have to do O(n) work for so much of the graph. But to do that we need a
    // union-find that lets us declare it over a potentially large space
    // without filling it all in.
    
    // So we do this the easy way and compute all the merges for all adjacency
    // components, including tiny/numerous ones, right now.
    
    // If we ever change this, we should also make MergedAdjacencyGraph
    // stackable, to save a copy when we want to do further merges but keep the
    // old state.
    
    graph->for_each_edge([&](const handlegraph::edge_t& e) {
        // Get the inward-facing version of the second handle
        auto into_b = graph->flip(e.second);
        
        // Merge to create initial adjacency components
        merge(e.first, into_b);
    });
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const MergedAdjacencyGraph& other) : MergedAdjacencyGraph(other.graph) {
    other.for_each_membership([&](handle_t head, handle_t member) {
        // For anything in a component, other than its head, do the merge with the head.
        merge(head, member);
    });
}

void IntegratedSnarlFinder::MergedAdjacencyGraph::merge(handle_t into_a, handle_t into_b) {
    // Get ranks and merge
    union_find.union_groups(uf_rank(into_a), uf_rank(into_b));
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::find(handle_t into) const {
    // Get rank, find head, and get handle
    return uf_handle(union_find.find_group(uf_rank(into)));
}
    
void IntegratedSnarlFinder::MergedAdjacencyGraph::for_each_membership(const function<void(handle_t, handle_t)>& iteratee) const {
    // We do this weird iteration because it's vaguely efficient in the union-find we use.
    vector<vector<size_t>> uf_components = union_find.all_groups();
    
    for (auto& component : uf_components) {
        // For each component
        for (size_t i = 1; i < component.size(); i++) {
            // For everything other than the head, announce with the head.
            iteratee(uf_handle(component[0]), uf_handle(component[i]));
        }
    }
}

    




////////////////////////////////////////////////////////////////////////////////////////////





IntegratedSnarlFinder::IntegratedSnarlFinder(const PathHandleGraph& graph) : graph(&graph) {
    // Nothing to do!
}

void IntegratedSnarlFinder::for_each_snarl_including_trivial(const function<void(handle_t, handle_t)>& iteratee) const {
    // Do the actual snarl finding work and then feed the iteratee our snarls.
    
    // We need a union-find over the adjacency components of the graph, in which we will build the cactus graph.
    MergedAdjacencyGraph cactus(graph);
    
    // It magically gets the adjacency components itself.
    
    // TODO: Now we need to do the 3 edge connected component merging, using Tsin's algorithm.
    
    // Then we need to copy the base Cactus graph so we can make the bridge forest
    MergedAdjacencyGraph forest(cactus);
    
    // TODO: Then we need to actually condense simple cycles in order to make
    // the bridge forest.
    //
    // TODO: We should find the lengths of the simple cycles in bases at the
    // same time.
    //
    // TODO: We also want to find the longest simple cycle for each connected
    // component of the bridge forest.
    
    // TODO: Then find the lengths of the bridge edge spanning trees. We want
    // to be able to know the lengths no matter where we root, so I guess for
    // each node we need to store the longest length in each direction.
    
    // TODO: For each bridge forest tree, find the longest bridge edge path and
    // the longest simple cycle, and root with whichever is longer.
    
    // TODO: Recursively root things down the bridge edge trees off what we've
    // already rooted.
    
    // TODO: When we encounter a snarl, feed it to the iteratee. 
}

SnarlManager IntegratedSnarlFinder::find_snarls() {
    // Start with an empty SnarlManager
    SnarlManager snarl_manager;
    
    for_each_snarl_including_trivial([&](handle_t start_inward, handle_t end_outward) {
        // For every snarl, including the trivial ones
        
        // Make a Protobuf version of it
        vg::Snarl proto_snarl;
        
        // Convert boundary handles to Visits
        proto_snarl.mutable_start()->set_node_id(graph->get_id(start_inward));
        proto_snarl.mutable_start()->set_backward(graph->get_is_reverse(start_inward));
        proto_snarl.mutable_end()->set_node_id(graph->get_id(end_outward));
        proto_snarl.mutable_end()->set_backward(graph->get_is_reverse(end_outward));
        
        // TODO: Determine snarl metadata somehow. Do a postorder traversal? Or make the manager somehow do it?
    
        // Add the Protobuf version of the snarl to the SNarl Manager
        snarl_manager.add_snarl(proto_snarl);
    });
    
    // Let the snarl manager compute all its indexes
    snarl_manager.finish();
    
    // Give it back
    return snarl_manager;
}

}
