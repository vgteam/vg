#ifndef VG_NAME_MAPPER_HPP
#define VG_NAME_MAPPER_HPP

/** \file
 * name_mapper.hpp: defines a class that maps back and forth from VCF-space
 * contig names to graph or FASTA-space contig names.
 */

#include <map>
#include <string>

namespace vg {

using namespace std;

/**
 * Class to do name mapping (or to mix in and provide name mapping functionality
 * to other classes.
 */
class NameMapper {

public:

    /**
     * Add a name mapping between a VCF contig name and a FASTA sequence name or
     * graph path name. Both must be unique.
     */
    void add_name_mapping(const string& vcf_name, const string& fasta_name);
    
    /**
     * Convert the given VCF contig name to a FASTA sequence or graph path name,
     * through the rename mappings.
     */
    string vcf_to_fasta(const string& vcf_name) const;
    
    /**
     * Convert the given FASTA sequence name or graph path name to a VCF contig
     * name, through the rename mappings.
     */
    string fasta_to_vcf(const string& fasta_name) const;
    
protected:
    
    /// This map maps from VCF sequence names to FASTA sequence names. If a
    /// VCF sequence name doesn't appear in here, it gets passed through
    /// unchanged.
    map<string, string> vcf_to_fasta_renames;
    
    /// This is the reverse map from FASTA sequence name to VCF sequence name.
    map<string, string> fasta_to_vcf_renames;

};

}

#endif
