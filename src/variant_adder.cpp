#include "variant_adder.hpp"

namespace vg {

using namespace std;

VariantAdder::VariantAdder(VG& graph) : graph(graph), sync(graph) {
    // Nothing to do!
}

void VariantAdder::add_variants(vcflib::VariantCallFile* vcf) {
    
    // Make a buffer
    WindowedVcfBuffer buffer(vcf, variant_range);
    
    // Count how many variants we have done
    size_t variants_processed = 0;
    
    while(buffer.next()) {
        // For each variant in its context of nonoverlapping variants
        vcflib::Variant* variant;
        vector<vcflib::Variant*> before;
        vector<vcflib::Variant*> after;
        tie(before, variant, after) = buffer.get_nonoverlapping();
    
        // Where is it?
        auto variant_path_name = vcf_to_fasta(variant->sequenceName);
        auto& variant_path_offset = variant->position; // Already made 0-based by the buffer
        
        // Make the list of all the local variants in one vector
        vector<vcflib::Variant*> local_variants{before};
        local_variants.push_back(variant);
        copy(after.begin(), after.end(), back_inserter(local_variants));
        
#ifdef debug
        cerr << "Local variants: ";
        for (auto* v : local_variants) {
            cerr << vcf_to_fasta(v->sequenceName) << ":" << v->position << " ";
        }
        cerr << endl;
#endif
        
        // Where does the group start?
        size_t group_start = local_variants.front()->position;
        // And where does it end (exclusive)?
        size_t group_end = local_variants.back()->position + local_variants.back()->ref.size();
        
        
        // Grab the sequence of the path, which won't change
        const string& path_sequence = sync.get_path_sequence(variant_path_name);
        
        // Get the leading and trailing ref sequence on either side of this group of variants (to pin the outside variants down).

        // On the left we want either flank_range bases or all the bases before the first base in the group.
        size_t left_context_length = max(min((int64_t)flank_range, (int64_t) group_start), (int64_t) 0);
        // On the right we want either flank_range bases or all the bases after the last base in the group.
        size_t right_context_length = min(path_sequence.size() - group_end, (size_t) flank_range);
    
        string left_context = path_sequence.substr(group_start - left_context_length, left_context_length);
        string right_context = path_sequence.substr(group_end, right_context_length);
        
        // Get the unique haplotypes
        auto haplotypes = get_unique_haplotypes(local_variants, &buffer);
        
        // Track the total bp of haplotypes
        size_t total_haplotype_bases = 0;
        
        // Track the total graph size for the alignments
        size_t total_graph_bases = 0;
        
#ifdef debug
        cerr << "Have " << haplotypes.size() << " haplotypes for variant " << *variant << endl;
#endif
        
        for (auto& haplotype : haplotypes) {
            // For each haplotype
            
            // TODO: since we lock repeatedly, neighboring variants will come in
            // in undefined order and our result is nondeterministic.
            
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

            // Count all the bases
            total_haplotype_bases += to_align.str().size();
            
            // Make a request to lock the subgraph
            GraphSynchronizer::Lock lock(sync, variant_path_name, variant_path_offset, (variant_range + flank_range * 2));
            
#ifdef debug
            cerr << "Waiting for lock on " << variant_path_name << ":" << variant_path_offset << endl;
#endif
            
            // Block until we get it
            lock_guard<GraphSynchronizer::Lock> guard(lock);
            
#ifdef debug
            cerr << "Got lock on " << variant_path_name << ":" << variant_path_offset << endl;
#endif            
                
#ifdef debug
            cerr << "Got " << lock.get_subgraph().length() << " bp in " << lock.get_subgraph().size() << " nodes" << endl;
#endif
                
            // Record the size of graph we're aligning to in bases
            total_graph_bases += lock.get_subgraph().length();
            
            // Do the alignment
            Alignment aln = lock.get_subgraph().align(to_align.str(), 0, false, false, 30);

#ifdef debug            
            cerr << "Alignment: " << pb2json(aln) << endl;
#endif
            
            // Make this path's edits to the original graph and get the
            // translations.
            auto translations = lock.apply_edit(aln.path());
                
#ifdef debug
            cerr << "Translations: " << endl;
            for (auto& t : translations) {
                cerr << "\t" << pb2json(t) << endl;
            }
#endif
                
        }
        
        if (variants_processed++ % 1000 == 0) {
            cerr << "Variant " << variants_processed << ": " << haplotypes.size() << " haplotypes at "
                << variant->sequenceName << ":" << variant->position << ": "
                << (total_haplotype_bases / haplotypes.size()) << " bp vs. "
                << (total_graph_bases / haplotypes.size()) << " bp haplotypes vs. graphs average" << endl;
        }
    }
    
    
}

set<vector<int>> VariantAdder::get_unique_haplotypes(const vector<vcflib::Variant*>& variants, WindowedVcfBuffer* cache) const {
    set<vector<int>> haplotypes;
    
    if (variants.empty()) {
        // Nothing's there
        return haplotypes;
    }
    
    for (size_t sample_index = 0; sample_index < variants.front()->sampleNames.size(); sample_index++) {
        // For every sample
        auto& sample_name = variants.front()->sampleNames[sample_index];
        
        // Make its haplotype(s) on the region. We have a map from haplotype
        // number to actual vector. We'll tack stuff on the ends when they are
        // used, then throw out any that aren't full-length.
        map<size_t, vector<int>> sample_haplotypes;
        
        
        for (auto* variant : variants) {
            // Get the genotype for each sample
            const vector<int>* genotype;
            
            if (cache != nullptr) {
                // Use the cache provided by the buffer
                genotype = &cache->get_parsed_genotypes(variant).at(sample_index);
            } else {
                // Parse from the variant ourselves
                auto genotype_string = variant->getGenotype(sample_name);
            
                // Fake it being phased
                replace(genotype_string.begin(), genotype_string.end(), '/', '|');
                
                genotype = new vector<int>(vcflib::decomposePhasedGenotype(genotype_string));
            }
            
#ifdef debug
            cerr << "Genotype of " << sample_name << " at " << variant->position << ": ";
            for (auto& alt : *genotype) {
                cerr << alt << " ";
            }
            cerr << endl;
#endif
            
            for (size_t phase = 0; phase < genotype->size(); phase++) {
                // For each phase in the genotype
                
                // Get the allele number and ignore missing data
                int allele_index = (*genotype)[phase];
                if (allele_index == vcflib::NULL_ALLELE) {
                    allele_index = 0;
                }
                
                // Stick each allele number at the end of its appropriate phase
                sample_haplotypes[phase].push_back(allele_index);
            }
            
            if (cache == nullptr) {
                // We're responsible for this vector
                delete genotype;
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
        result << sync.get_path_sequence(vcf_to_fasta(variant->sequenceName)).substr(sep_start, sep_length);

        // Then put the appropriate allele of this variant.
        result << variant->alleles.at(haplotype.at(i));
    }
    
    return result.str();
}

}





















