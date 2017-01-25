#include "variant_adder.hpp"

namespace vg {

using namespace std;

VariantAdder::VariantAdder(VG& graph) : graph(graph) {
    // Nothing to do!
}


    
// We need a function to grab the index for a path
PathIndex& VariantAdder::get_path_index(const string& path_name) {
    if (!indexes.count(path_name)) {
        // Not already made. Generate it.
        indexes.emplace(piecewise_construct,
            forward_as_tuple(path_name), // Make the key
            forward_as_tuple(graph, path_name, true)); // Make the PathIndex
    }
    return indexes.at(path_name);
}

void VariantAdder::add_variants(vcflib::VariantCallFile* vcf) {
    // We collect all the edits we need to make, and then make them all at
    // once.
    vector<Path> edits_to_make;
    
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
        
        if (!graph.paths.has_path(variant_path_name)) {
            // Explode if the path is not in the vg
            
            cerr << "error:[vg add] could not find path" << variant_path_name << " in graph" << endl;
            throw runtime_error("Missing path " + variant_path_name);
        }
        
        // Grab the path index
        auto& index = get_path_index(variant_path_name);
        
        // Find the node that the variant falls on
        NodeSide center = index.at_position(variant_path_offset);
        
        // Extract its graph context for realignment
        VG context;
        graph.nonoverlapping_node_context_without_paths(graph.get_node(center.node), context);
        // TODO: how many nodes should this be?
        // TODO: write/copy the search from xg so we can do this by node length.
        graph.expand_context(context, 10, false);
        
        // Make the list of all the local variants in one vector
        vector<vcflib::Variant*> local_variants{before};
        local_variants.push_back(variant);
        copy(after.begin(), after.end(), back_inserter(local_variants));
        
#ifdef debug
        cerr << "Local variants: ";
        for (auto* v : local_variants) {
            cerr << v->sequenceName << ":" << v->position << " ";
        }
        cerr << endl;
#endif
        
        // Where does the group start?
        size_t group_start = local_variants.front()->position;
        // And where does it end (exclusive)?
        size_t group_end = local_variants.back()->position + local_variants.back()->ref.size();
        
        // Get the leading and trailing ref sequence on either side of this group of variants (to pin the outside variants down).

        // On the left we want either flank_range bases or all the bases before the first base in the group.
        size_t left_context_length = max(min((int64_t)flank_range, (int64_t) group_start - (int64_t) 1), (int64_t) 0);
        // On the right we want either flank_range bases or all the bases after the last base in the group.
        size_t right_context_length = min(index.sequence.size() - group_end, (size_t) flank_range);
    
        string left_context = index.sequence.substr(group_start - left_context_length, left_context_length);
        string right_context = index.sequence.substr(group_end, right_context_length);
        
        // Get the unique haplotypes
        auto haplotypes = get_unique_haplotypes(local_variants);
        
#ifdef debug
        cerr << "Have " << haplotypes.size() << " haplotypes for variant " << *variant << endl;
#endif
        
        for (auto& haplotype : haplotypes) {
            // For each haplotype
            
#ifdef debug
            cerr << "Haplotype ";
            for (auto& allele_number : haplotype) {
                cerr << allele_number << " ";
            }
            cerr << endl;
#endif
            
            // Make its combined string
            stringstream to_align;
            to_align << left_context << haplotype_to_string(haplotype, local_variants) << right_context;
            
#ifdef debug
            cerr << "Align " << to_align.str() << endl;
#endif
            
            // Align it
            // TODO: we want to give a full-length bonus to keep the
            // ends attached when a variant is near the end of the
            // reference. But we can't without turning on pinned mode,
            // which is incorrect because we don't have a read that we
            // know runs up to the end of the graph.
            Alignment aln = context.align(to_align.str(), 0, false, false, 30);
            
            // Queue it up to edit the graph
            edits_to_make.push_back(aln.path());
        }
    }
    
    // Then at the end of the VCF, edit the graph
    graph.edit(edits_to_make);
}

set<vector<int>> VariantAdder::get_unique_haplotypes(const vector<vcflib::Variant*>& variants) const {
    set<vector<int>> haplotypes;
    
    if (variants.empty()) {
        // Nothing's there
        return haplotypes;
    }
    
    for (auto& sample_name : variants.front()->sampleNames) {
        // For every sample
        
        // Make its haplotype(s) on the region. We have a map from haplotype
        // number to actual vector. We'll tack stuff on the ends when they are
        // used, then throw out any that aren't full-length.
        map<size_t, vector<int>> sample_haplotypes;
        
        
        for (auto* variant : variants) {
            // Get the genotype for each sample
            auto genotype = variant->getGenotype(sample_name);
            
            // Fake it being phased
            replace(genotype.begin(), genotype.end(), '/', '|');
            
            auto alts = vcflib::decomposePhasedGenotype(genotype);
            
            for (size_t phase = 0; phase < alts.size(); phase++) {
                // Stick each allele number at the end of its appropriate phase
                sample_haplotypes[phase].push_back(alts[phase]);
            }
        }
        
        for (auto& kv : sample_haplotypes) {
            auto& haplotype = kv.second;
            // For every haplotype in this sample
            if (haplotype.size() != variants.size()) {
                // If it's not the full length, it means some variants don't
                // have it. Skip.
                continue;
            }
            
            // Otherwise, add it to the set of observed haplotypes
            haplotypes.insert(haplotype);
        }
    }
    
    // After processing all the samples, return the unique haplotypes
    return haplotypes;
    
}

string VariantAdder::haplotype_to_string(const vector<int>& haplotype, const vector<vcflib::Variant*>& variants) {
    // We'll fill this in with variants and separating sequences.
    stringstream result;
    
    // These lists need to be in 1 to 1 correspondence
    assert(haplotype.size() == variants.size());
    
    if (variants.empty()) {
        // No variants means no string representation.
        return "";
    }
    
    // Do the first variant
    result << variants.front()->alleles.at(haplotype.at(0));
    
    for (size_t i = 1; i < variants.size(); i++) {
        // For each subsequent variant
        auto* variant = variants.at(i);
        auto* last_variant = variants.at(i - 1);
        
        // Do the intervening sequence.
        // Where does that sequence start?
        size_t sep_start = last_variant->position + last_variant->ref.size();
        // And how long does it run?
        size_t sep_length = variant->position - sep_start;
        
        // Pull out the separator sequence and tack it on.
        result << get_path_index(variant->sequenceName).sequence.substr(sep_start, sep_length);

        // Then put the appropriate allele of this variant
        result << variant->alleles.at(haplotype.at(i));
    }
    
    return result.str();
}

}





















