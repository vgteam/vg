#include "vcf_genotyper.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"

//#define debug

namespace vg {

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<int>& ref_path_ploidies,
                           FastaReference* ref_fasta,
                           FastaReference* ins_fasta,
                           AlignmentEmitter* aln_emitter,
                           bool traversals_only,
                           bool gaf_output,
                           size_t trav_padding) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    input_vcf(variant_file),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths, ref_fasta, ins_fasta, snarl_caller.get_skip_allele_fn()),
    traversals_only(traversals_only),
    gaf_output(gaf_output) {

    scan_contig_lengths();

    assert(ref_paths.size() == ref_path_ploidies.size());
    for (int i = 0; i < ref_paths.size(); ++i) {
        path_to_ploidy[ref_paths[i]] = ref_path_ploidies[i];
    }
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // could be that our graph is a subgraph of the graph the snarls were computed from
    // so bypass snarls we can't process
    if (!graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        return false;
    }

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

        // find a path range corresponding to our snarl by way of the VCF variants.
        tuple<string, size_t, size_t> ref_positions = get_ref_positions(variants);

        // just print the traversals if requested
        if (traversals_only) {
            assert(gaf_output);
            // todo: can't get ref position here without pathposition graph
            emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, "", -1);
            return true;
        }
        
        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype;
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, path_to_ploidy[get<0>(ref_positions)],
                                                                        get<0>(ref_positions),  make_pair(get<1>(ref_positions), get<2>(ref_positions)));

        assert(trav_genotype.size() <= 2);

        if (gaf_output) {
            // todo: can't get ref position here without pathposition graph
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, "", -1);
            return true;
        }

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            vector<int> vcf_alleles;
            set<int> used_vcf_alleles;
            string vcf_genotype;
            vector<SnarlTraversal> vcf_traversals(variants[i]->alleles.size());            
            if (trav_genotype.empty()) {
                vcf_genotype = "./.";
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
                    if (!used_vcf_alleles.count(vcf_allele)) {                    
                        vcf_alleles.push_back(vcf_allele);
                        used_vcf_alleles.insert(vcf_allele);
                        vcf_traversals[vcf_allele] = travs[trav_allele];
                    }
                }
                // add traversals that correspond to vcf genotypes that are not
                // present in the traversal_genotypes
                for (int j = 0; j < travs.size(); ++j) {
                    int vcf_allele = alleles[j].second[i];
                    if (!used_vcf_alleles.count(vcf_allele)) {
                        vcf_traversals[vcf_allele] = travs[j];
                        used_vcf_alleles.insert(vcf_allele);
                    }
                }
            }
            // create an output variant from the input one
            vcflib::Variant out_variant;
            out_variant.sequenceName = variants[i]->sequenceName;
            out_variant.position = variants[i]->position;
            out_variant.id = variants[i]->id;
            out_variant.ref = variants[i]->ref;
            out_variant.alt = variants[i]->alt;
            out_variant.alleles = variants[i]->alleles;
            out_variant.filter = "PASS";
            out_variant.updateAlleleIndexes();

            // add the genotype
            out_variant.format.push_back("GT");
            auto& genotype_vector = out_variant.samples[sample_name]["GT"];
            genotype_vector.push_back(vcf_genotype);

            // add some info
            snarl_caller.update_vcf_info(snarl, vcf_traversals, vcf_alleles, trav_call_info, sample_name, out_variant);

            // print the variant
            add_variant(out_variant);
        }
        return true;
    }
    
    return false;

}

string VCFGenotyper::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    assert(contig_length_overrides.empty()); // using this override makes no sense

    // get the contig length overrides from the VCF
    vector<size_t> vcf_contig_lengths;
    auto length_map = scan_contig_lengths();
    for (int i = 0; i < ref_paths.size(); ++i) {
        vcf_contig_lengths.push_back(length_map[ref_paths[i]]);
    }
    
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, vcf_contig_lengths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

tuple<string, size_t, size_t> VCFGenotyper::get_ref_positions(const vector<vcflib::Variant*>& variants) const {
    // if there is more than one path in our snarl (unlikely for most graphs we'll vcf-genoetype)
    // then we return the one with the biggest interval
    map<string, pair<size_t, size_t>> path_offsets;
    for (const vcflib::Variant* var : variants) {
        if (path_offsets.count(var->sequenceName)) {
            pair<size_t, size_t>& record = path_offsets[var->sequenceName];
            record.first = std::min((size_t)var->position, record.first);
            record.second = std::max((size_t)var->position + var->ref.length(), record.second);
        } else {
            path_offsets[var->sequenceName] = make_pair(var->position, var->position + var->ref.length());
        }
    }

    string ref_path;
    size_t ref_range_size = 0;
    pair<size_t, size_t> ref_range;
    for (auto& path_offset : path_offsets) {
        size_t len = path_offset.second.second - path_offset.second.first;
        if (len > ref_range_size) {
            ref_range_size = len;
            ref_path = path_offset.first;
            ref_range = path_offset.second;
        }
    }

    return make_tuple(ref_path, ref_range.first, ref_range.second);
}

unordered_map<string, size_t> VCFGenotyper::scan_contig_lengths() const {

    unordered_map<string, size_t> ref_lengths;
    
    // copied from dumpContigsFromHeader.cpp in vcflib
    vector<string> headerLines = split(input_vcf.header, "\n");
    for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++) {
        if((*it).substr(0,8) == "##contig"){
            string contigInfo = (*it).substr(10, (*it).length() -11);
            vector<string> info = split(contigInfo, ",");
            string id;
            int64_t length = -1;
            for(vector<string>::iterator sub = info.begin(); sub != info.end(); sub++) {
                vector<string> subfield = split((*sub), "=");
                if(subfield[0] == "ID"){
                    id = subfield[1];
                }
                if(subfield[0] == "length"){
                    length = parse<int>(subfield[1]);
                }
            }
            if (!id.empty() && length >= 0) {
                ref_lengths[id] = length;
            }
        }
    }

    return ref_lengths;
}
}

