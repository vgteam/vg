#ifndef VG_MSA_CONVERTER_HPP_INCLUDED
#define VG_MSA_CONVERTER_HPP_INCLUDED

/** \file
 * msa_converter.hpp: contains a class that can construct VGs from clustal MSAs
 */

#include <vector>
#include <set>
#include <map>
#include <cstdlib>
#include <functional>
#include <regex>

#include <vg/vg.pb.h>

#include "vg.hpp"

namespace vg {

using namespace std;

    class MSAConverter : public Progressive  {
    public:
        
        MSAConverter();
        ~MSAConverter();
        
        void load_alignments(istream& in, string format = "fasta");
        
        VG make_graph(bool keep_paths = true, size_t max_node_length = numeric_limits<size_t>::max());
        
    private:
        
        vector<unordered_map<string, string>> alignments;
        
    };

}

#endif
