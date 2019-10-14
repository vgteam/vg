#include "mcmc_caller.hpp"  
#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"

namespace vg{

    /**
     * MCMCCaller : Inherits from VCFOutputCaller    
     */         
    MCMCCaller::MCMCCaller(const PathPositionHandleGraph& graph,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets,
                           const vector<size_t>& ref_path_lengths,
                           ostream& out_stream) :
    graph(graph),sample_name(sample_name), ref_paths(ref_paths), ref_path_offsets(ref_path_offsets),
    ref_path_lengths(ref_path_lengths),out_stream(out_stream) {
        

    }

    
    MCMCCaller:: ~MCMCCaller(){

    }


    void process_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const string& ref_path_name) const {
            // convert traversal to string
    function<string(const SnarlTraversal&)> trav_string = [&](const SnarlTraversal& trav) {
        string seq;
        for (int i = 0; i < trav.visit_size(); ++i) {
            seq += graph.get_sequence(graph.get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
        }
        return seq;
    };

    vcflib::Variant out_variant;

    // when calling alt/alt, the reference traversal doesn't end up in called_traversals.
    // this should get changed, but in the meantime we add it back here (as we need it for
    // the VCF output)
    // udpate: the reference traversal will be there when re-genotyping, but we can leave this logic
    // in case we want to ever add an option to toggle this.
    vector<SnarlTraversal> site_traversals;
    vector<int> site_genotype;
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == 0) {
            site_traversals.push_back(called_traversals[i]);
            break;
        }
    }
    if (site_traversals.empty()) {
        site_traversals.push_back(get_reference_traversal(snarl, trav_finder));
    }
    out_variant.ref = trav_string(site_traversals[0]);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;    
    allele_to_gt[out_variant.ref] = 0;    
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == 0) {
            site_genotype.push_back(0);
        } else {
            string allele_string = trav_string(called_traversals[i]);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_traversals.push_back(called_traversals[i]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    for (auto& allele_gt : allele_to_gt) {
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;
    }

    // fill out the rest of the variant
    out_variant.sequenceName = ref_path_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get_ref_position(snarl, ref_path_name).first + ref_offsets.find(ref_path_name)->second + 1; 
    out_variant.id = std::to_string(snarl.start().node_id()) + "_" + std::to_string(snarl.end().node_id());
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    
    stringstream vcf_gt;
    for (int i = 0; i < site_genotype.size(); ++i) {
        vcf_gt << site_genotype[i];
        if (i != site_genotype.size() - 1) {
            vcf_gt << "/";
        }
    }
    genotype_vector.push_back(vcf_gt.str());

    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true);
    flatten_common_allele_ends(out_variant, false);
        

    }

    


}  





