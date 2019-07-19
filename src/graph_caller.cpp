#include "graph_caller.hpp"

namespace vg {

GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager,
                         ostream& out_stream) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager), out_stream(out_stream) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::call_top_level_snarls(bool recurse_on_fail) {

    // Used to recurse on children of parents that can't be called
    vector<const Snarl*> snarl_queue;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {
        bool was_called = call_snarl(*snarl);
        if (!was_called && recurse_on_fail) {
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
            snarl_queue.insert(snarl_queue.end(), children.begin(), children.end());
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_snarl_parallel(process_snarl);

    // Then recurse on any children the snarl caller failed to handle
    while (!snarl_queue.empty()) {
        vector<const Snarl*> cur_queue;
        std::swap(snarl_queue, cur_queue);
#pragma omp parallel for
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
  
}

void GraphCaller::vcf_header(const PathHandleGraph& graph, const vector<path_handle_t>& contigs) const {
    out_stream << "##fileformat=VCFv4.2" << endl;
    for (auto contig : contigs) {
        size_t length = 0;
        for (handle_t handle : graph.scan_path(contig)) {
            length += graph.get_length(handle);
        }
        out_stream << "##contig=<ID=" << graph.get_path_name(contig) << ",length=" << length << ">" << endl;
    }
}

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const vector<string>& ref_paths,
                           ostream& out_stream) :
    GraphCaller(snarl_caller, snarl_manager, out_stream),
    graph(graph),
    regions(regions),
    input_vcf(variant_file),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths) {
    
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // get our traversals out of the finder
    vector<pair<SnarlTraversal, vector<int>>> alleles;
    vector<vcflib::Variant*> variants;
    std::tie(alleles, variants) = traversal_finder.find_allele_traversals(snarl);

    if (!alleles.empty()) {

        // hmm, maybe find a way not to copy?
        vector<SnarlTraversal> travs;
        travs.reserve(alleles.size());
        for (const auto& ta : alleles) {
            travs.push_back(ta.first);
        }

        // find the reference traversal
        // todo: is it the reference always first?
        int ref_trav_idx = -1;
        for (int i = 0; i < alleles.size() && ref_trav_idx < 0; ++i) {
            if (std::all_of(alleles[i].second.begin(), alleles[i].second.end(), [](int x) {return x == 0;})) {
                ref_trav_idx = i;
            }
        }

        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype = snarl_caller.genotype(snarl, travs, ref_trav_idx, 2);

        assert(trav_genotype.empty() || trav_genotype.size() == 2);

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            string vcf_genotype;;
            if (trav_genotype.empty()) {
                vcf_genotype = ".";
            } else {
                // map our traversal genotype to a vcf variant genotype
                // using the information out of the traversal finder
                for (int j = 0; j < trav_genotype.size(); ++j) {
                    int trav_allele = trav_genotype[j];
                    int vcf_allele = alleles[trav_allele].second[i];
                    vcf_genotype += std::to_string(vcf_allele);
                    if (j < trav_genotype.size() - 1) {
                        vcf_genotype += "/";
                    }
                }

                // testing
                out_stream << *variants[i] << " ==> " << vcf_genotype << endl;
            }
        }
        return true;
    }

    return false;

}



}

