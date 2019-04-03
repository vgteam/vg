#ifndef VG_GFF_READER_HPP_INCLUDED
#define VG_GFF_READER_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <limits>
#include <functional>

namespace vg {
    
    using namespace std;
    
    /**
     * A package of the information contained in a GFF3 record. The null "." entries in a
     * a GFF are parsed into empty strings or the default values of the numerical fields
     * as given below.
     */
    struct GFFRecord {
    public:
        GFFRecord() = default;
        ~GFFRecord() = default;
        
        string sequence_id;
        string source;
        string type;
        // 0-based indexing, unlike the actual GFF standard
        int64_t start = -1;
        // 0-based, inclusive
        int64_t end = -1;
        double score = numeric_limits<double>::quiet_NaN();
        bool strand_is_rev = false;
        int32_t phase = -1;
        string attributes;
        
        map<string, string> parse_attributes();
    };
    
    /**
     * A class that can parse and iterate over a GFF3 file.
     */
    class GFFReader {
    public:
        GFFReader(istream& in);
        ~GFFReader() = default;
        
        void for_each_gff_record(function<void(const GFFRecord&)>& lambda);
        
    private:
        istream& in;
    };

}

#endif
