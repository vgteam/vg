#include "genotypekit.hpp"

namespace vg {

using namespace std;


CactusSiteFinder::CactusSiteFinder(VG& graph, const string& hint_path_name): graph(graph), hint_path_name(hint_path_name) {
    // Make sure the graph is sorted.
    // cactus needs the nodes to be sorted in order to find a source and sink.
    graph.sort();
}

void CactusSiteFinder::for_each_site_parallel(const function<void(const NestedSite)>& lambda) {

    // Set up our output vector
    vector<NestedSite> to_return;
    
    // get endpoints using node ranks
    pair<NodeSide, NodeSide> source_sink = graph.paths.has_path(hint_path_name) ? 
        get_cactus_source_sink(graph, hint_path_name)
        : get_cactus_source_sink(graph);
        
    // Don't keep going if we can't find sources/sinks
    assert(graph.has_node(source_sink.first.node));
    assert(graph.has_node(source_sink.second.node));

    // Get the bubble tree in Cactus format
    BubbleTree bubble_tree = cactusbubble_tree(graph, source_sink);

    // Convert to NestedSites
    
    // We use this to hold the NestedSites that are children until their parents
    // are ready to be converted.
    map<BubbleTree::Node*, NestedSite> converted_children;

    bubble_tree.for_each_postorder([&](BubbleTree::Node* node) {
        // Process children before parents so we can embed them in the parent.
        
        Bubble& bubble = node->v;
        if (node != bubble_tree.root) {
            // If we aren't the root node of the tree, we need to be a NestedSite
            
            // We're going to fill in this NestedSite.
            NestedSite& to_fill = converted_children[node];
            
            // Set up the start and end
            NodeTraversal start(graph.get_node(bubble.start.node), !bubble.start.is_end);
            NodeTraversal end(graph.get_node(bubble.end.node), bubble.end.is_end);
            // Make sure to preserve original endpoint
            // ordering, because swapping them without flipping their
            // orientation flags will make an inside-out site.
            to_fill.start = start;
            to_fill.end = end;
            
            for(id_t node_id : bubble.contents) {
                // Convert all the directly contained nodes to pointers
                to_fill.nodes.insert(graph.get_node(node_id));
            }
            
            for(BubbleTree::Node* child_node : node->children) {
                // Attach all the children by moving them out of our map.
                assert(converted_children.count(child_node));
                to_fill.children.emplace_back(std::move(converted_children[child_node]));
                converted_children.erase(child_node);
                
                // Fill in child borders with the NodeTraversals leading into the children.
                to_fill.child_border_index[to_fill.children.back().start] = to_fill.children.size() - 1;
                to_fill.child_border_index[to_fill.children.back().end.reverse()] = to_fill.children.size() - 1;
                
            }
            
            // Now do all the edges
            
            for(Node* internal_node : to_fill.nodes) {
                // Collect edges on internal nodes
                if(internal_node == to_fill.start.node || internal_node == to_fill.end.node) {
                    // Look only at internal nodes (not start or end of site)
                    continue;
                }
                
                // Since these aren't the start or end nodes of either this site
                // or any child site, all the edges on them must be part of this
                // site.
                auto node_edges = graph.edges_of(internal_node);
                // Dump them all in
                to_fill.edges.insert(node_edges.begin(), node_edges.end());
            }

            for(auto& child : to_fill.children) {
                // Now do edges between sites linked by just an edge. For each child site...
                
                // Get the outer side of the start traversal
                NodeSide start_outer_side = NodeSide(child.start.node->id(), child.start.backward);
                
                // Get the connected sides
                auto start_connected_sides = graph.sides_of(start_outer_side);
                
                for(auto& connected_side : start_connected_sides) {
                    // Pull in all the edges between the outer side of the
                    // contained site and everything else. They must be within
                    // this parent site.
                    to_fill.edges.insert(graph.get_edge(start_outer_side, connected_side));
                }
                
                // Repeat for the end of the child site
                
                // The outer side of the end of the site is a right side if the
                // child side doesn't end with a backwards node.
                NodeSide end_outer_side = NodeSide(child.end.node->id(), !child.end.backward);
                auto end_connected_sides = graph.sides_of(end_outer_side);
                for(auto& connected_side : end_connected_sides) {
                    to_fill.edges.insert(graph.get_edge(end_outer_side, connected_side));
                }
                
            }
            
            // Finally do edges on the inside sides of this site's start and
            // end. Those are the only ones not yet covered.
            NodeSide start_inner_side = NodeSide(to_fill.start.node->id(), !to_fill.start.backward);
            for(auto& connected_side : graph.sides_of(start_inner_side)) {
                // Wire up the inner side of the start node to all its connected sides
                to_fill.edges.insert(graph.get_edge(start_inner_side, connected_side));
            }
            
            NodeSide end_inner_side = NodeSide(to_fill.end.node->id(), to_fill.end.backward);
            for(auto& connected_side : graph.sides_of(end_inner_side)) {
                // Wire up the inner side of the end node to all its connected sides
                to_fill.edges.insert(graph.get_edge(end_inner_side, connected_side));
            }
            
            // Now we're done with all the edges.
            
        } 
    });
    
    // Now emit all the top-level sites
    
    for(auto& child_and_site : converted_children) {
    
        // OpenMP doesn't like to let us use the reference, even though
        // we know it will survive the tasks. We grab a pointer to make
        // it happy.
        auto* lambda_pointer = &lambda;
        
        // Ditto for the reference into converted_children
        auto* site = &child_and_site.second;
        
        // Operate on the site. Make sure to move into task stack storage.
        #pragma omp task
        {
            (*lambda_pointer)(std::move(*site));
        }
        
    }
        
    // Don't return until all the lambda tasks are done.
    #pragma omp taskwait    
    
}

double FixedGenotypePriorCalculator::calculate_log_prior(const Genotype& genotype) {
    // Are all the alleles the same?
    bool all_same = true;
    // What is the common allele number (-1 for unset)
    int allele_value = -1;
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele in the genotype
        if(allele_value == -1) {
            // On the first one, grab it. Everyone else will have to match.
            allele_value = genotype.allele(i);
        }
        
        if(allele_value != genotype.allele(i)) {
            // There are two distinct allele values in here
            all_same = false;
            break;
        }
    }
    
    // Return the appropriate prior depending on whether the alleles are all the
    // same (homozygous) or not (heterozygous).
    return all_same ? homozygous_prior_ln : heterozygous_prior_ln;
}

}
