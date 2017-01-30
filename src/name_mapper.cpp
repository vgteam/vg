/**
 * \file name_mapper.cpp: contains implementations for the NameMapper to convert
 * between VCF and graph/FASTA contig names.
 */


#include "name_mapper.hpp"


namespace vg {

using namespace std;

void NameMapper::add_name_mapping(const string& vcf_name, const string& fasta_name) {
    // Fill in both one-way maps.
    // TODO: C++ doesn't have a 2-way map right?
    vcf_to_fasta_renames[vcf_name] = fasta_name;
    fasta_to_vcf_renames[fasta_name] = vcf_name;
#ifdef debug
    cerr << "Added rename of " << vcf_name << " to " << fasta_name << endl;
#endif
}

string NameMapper::vcf_to_fasta(const string& vcf_name) const {
    return vcf_to_fasta_renames.count(vcf_name) ? vcf_to_fasta_renames.at(vcf_name) : vcf_name;
}

string NameMapper::fasta_to_vcf(const string& fasta_name) const {
    return fasta_to_vcf_renames.count(fasta_name) ? fasta_to_vcf_renames.at(fasta_name) : fasta_name;
}

}

