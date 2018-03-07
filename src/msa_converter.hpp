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

#include "vg.pb.h"

#include "vg.hpp"

namespace vg {

using namespace std;

    class MSAConverter {
    public:
        // TODO: currently only MAF format, should write other sockets too
        MSAConverter(istream& in, size_t max_node_length = numeric_limits<size_t>::max());
        ~MSAConverter();
        
        VG make_graph(bool keep_paths = true);
        
    private:
        
        unordered_map<string, string> alignments;
        size_t max_node_length;
        
    };

}

#endif
