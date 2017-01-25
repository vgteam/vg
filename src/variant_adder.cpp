#include "variant_adder.hpp"
#include "path_index.hpp"

namespace vg {

using namespace std;

VariantAdder::VariantAdder(VG& graph) : graph(graph) {
    // Nothing to do!
}

void VariantAdder::add_variants(vcflib::VariantCallFile* vcf) {
    // We collect all the edits we need to make, and then make them all at
    // once.
    vector<Path> edits_to_make;
    
    // We need indexes of all the paths that variants happen on.
    map<string, PathIndex> indexes;
    
    // We need a function to grab the index for a path
    auto get_path_index = [&](const string& path_name) -> PathIndex& {
        if (!indexes.count(path_name)) {
            // Not already made. Generate it.
            indexes.emplace(piecewise_construct,
                forward_as_tuple(path_name), // Make the key
                forward_as_tuple(graph, path_name, true)); // Make the PathIndex
        }
        return indexes.at(path_name);
    };
    
    // Make a buffer
    WindowedVcfBuffer buffer(vcf, variant_range);
    
    while(buffer.next()) {
        // For each variant in its context of nonoverlapping variants
        vcflib::Variant* variant;
        vector<vcflib::Variant*> before;
        vector<vcflib::Variant*> after;
        tie(before, variant, after) = buffer.get_nonoverlapping();
    
        // Where is it?
        auto& variant_path_name = variant->sequenceName;
        auto& variant_path_offset = variant->position; // Already made 0-based by the buffer
        
        auto variant_ref_length = variant->ref.size();
        
        if (!graph.paths.has_path(variant_path_name)) {
            // Explode if the path is not in the vg
            
            cerr << "error:[vg add] could not find path" << variant_path_name << " in graph" << endl;
            throw runtime_error("Missing path " + variant_path_name);
        }

        // Grab the path index
        auto& index = get_path_index(variant_path_name);
    
        // Extract its left and right context from the appropriate path in the graph
        // On the left we want either 100 bases or all the bases before the first ref base.
        size_t left_context_length = max(min((int64_t)100, (int64_t) variant_path_offset - (int64_t) 1), (int64_t) 0);
        // On the right we want either 100 bases or all the bases after the last ref base.
        size_t right_context_length = min(index.sequence.size() - variant_path_offset - variant_ref_length, (size_t) 100);
    
        string left_context = index.sequence.substr(variant_path_offset - left_context_length, left_context_length);
        string right_context = index.sequence.substr(variant_path_offset + variant_ref_length, right_context_length);
        
        // Find the node that the variant falls on
        NodeSide center = index.at_position(variant_path_offset);
        
        // Extract its graph context for realignment
        VG context;
        graph.nonoverlapping_node_context_without_paths(graph.get_node(center.node), context);
        // TODO: how many nodes should this be?
        // TODO: write/copy the search from xg so we can do this by node length.
        graph.expand_context(context, 10, false);
        
        for (auto& alt : variant->alt) {
            // For each non-ref alt
        
            // Paste it in with the left and right context
            stringstream alt_stream;
            alt_stream << left_context << alt << right_context;
            
            // Align it to the subgraph.
            
            // TODO: we want to give a full-length bonus to keep the
            // ends attached when a variant is near the end of the
            // reference. But we can't without turning on pinned mode,
            // which is incorrect because we don't have a read that we
            // know runs up to the end of the graph.
            Alignment aln = context.align(alt_stream.str(), 0, false, false, 30);
            
            // Add the path to our collection of paths to add
            edits_to_make.push_back(aln.path());
        }
    }
    
    // Then at the end of the VCF, edit the graph
    graph.edit(edits_to_make);
}

}
