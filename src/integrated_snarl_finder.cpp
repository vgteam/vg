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
    VectorizableOverlayHelper overlay_helper;
    
    /// Keep a union-find over the ranks of the merged oriented handles that
    /// make up each component. Runs with include_children=true so we can find
    /// all the members of each group.
    structures::UnionFind union_find;
    
    /// If not null, we represent a set of further merges on top of this
    /// backing MergedAdjacencyGraph. The backing graph MAY NOT CHANGE as long
    /// as we are alive, and we DO NOT OWN IT.
    ///
    /// Our union-find will be over the ranked heads of components in the base
    /// MergedAdjacencyGraph, and our overlay helper will be unused.
    const MergedAdjacencyGraph* base;
    
    /// If we are an overlay ourselves, we define a dense rank space over the
    /// heads of the base's components.
    /// TODO: When we get a good sparse union-find, get rid of this.
    vector<handle_t> base_components;
    /// And we map from handle to index in base_components.
    unordered_map<handle_t, size_t> base_ranks;
    // TODO: is keeping all this really better than just copying the union-find?
    
    /// Get the rank corresponding to the given handle, in the union-find.
    /// Our ranks are 0-based.
    size_t uf_rank(handle_t into) const;
    
    /// Get the handle with the given rank in union-find space.
    /// Our ranks are 0-based.
    handle_t uf_handle(size_t rank) const;
    
public:
    /// Make a MergedAdjacencyGraph representing the graph of adjacency components of the given HandleGraph.
    MergedAdjacencyGraph(const HandleGraph* graph);
    
    /// Make a MergedAdjacencyGraph representing further merges of the
    /// components of the given immutable MergedAdjacencyGraph.
    MergedAdjacencyGraph(const MergedAdjacencyGraph* base);
    
    /// Given handles reading into two components, a and b, merge them into a single component.
    void merge(handle_t into_a, handle_t into_b);
    
    /// Find the handle heading the component that the given handle is in.
    handle_t find(handle_t into) const; 
    
    /// Count the number of components in the MergedAdjacencyGraph
    size_t size() const;
};

size_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_rank(handle_t into) const {
    if (base == nullptr) {
        // We have a vectorizable overlay
        return (overlay_helper.get()->id_to_rank(graph->get_id(into)) - 1) * 2 + (size_t) graph->get_is_reverse(into);
    } else {
        // We're really an overlay ourselves
        
    }
}

handle_t IntegratedSnarlFinder::MergedAdjacencyGraph::uf_handle(size_t rank) const {
    if (base == nullptr) {
        return graph->get_handle(overlay_helper.get()->rank_to_id(rank / 2 + 1), rank % 2);
    } else {
    }
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const HandleGraph* graph) : graph(graph),
    overlay_helper(), union_find(overlay_helper.apply(graph)->get_node_count() * 2, true), base(nullptr), base_components(), base_ranks()  {
    
    // TODO: we want the adjacency components that are just single edges
    // between two handles (i.e. trivial snarls) to be implicit, so we don't
    // have to do O(n) work for so much of the graph. But to do that we need a
    // union-find that lets us declare it over a potentially large space
    // without filling it all in.
    
    // So we do this the easy way and compute all the merges for all adjacency
    // components, including tiny/numerous ones, right now.
    
    graph->for_each_edge([&](const handlegraph::edge_t& e) {
        // Get the inward-facing version of the second handle
        auto into_b = graph->flip(e.second);
        
        // Merge to create initial adjacency components
        merge(e.first, into_b)
    });
}

IntegratedSnarlFinder::MergedAdjacencyGraph::MergedAdjacencyGraph(const MergedAdjacencyGraph* base) : graph(base->graph),
    overlay_helper(), union_find(base->size(), true), base(base)  {
    
    // Nothing to do!
}

IntegratedSnarlFinder::MergedAdjacencyGraph::merge(
    
IntegratedSnarlFinder::IntegratedSnarlFinder(const PathHandleGraph& graph) : graph(&graph) {
    // Nothing to do!
}

void IntegratedSnarlFinder::for_each_snarl_including_trivial(const function<void(handle_t, handle_t)>& iteratee) const {
    // Do the actual snarl finding work and then feed the iteratee our snarls.
    
    // We need our base graph to be vectorizable, so we can have ranks for all the handles, so we can do a union-find over them.
    PathVectorizableOverlayHelper overlay_helper;
    // Get or create a vectorizability facet of the base graph, which we can use to go between ranks and handles.
    VectorizableHandleGraph vectorized_graph = overlay_helper.apply(graph);
    
    
    // We need a union-find over the whole graph.
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
