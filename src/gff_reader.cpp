#include "gff_reader.hpp"

namespace vg {
    
    map<string, string> GFFRecord::parse_attributes() {
        
        map<string, string> parsed_attributes;
        stringstream attr_stream(attributes);
        
        string buffer;
        while (attr_stream.good()) {
            getline(attr_stream, buffer, ';');
            
            stringstream split_stream(buffer);
            
            string attr_type;
            string attr_value;
            
            getline(split_stream, attr_type, '=');
            getline(split_stream, attr_value, '\0');
            
            parsed_attributes[attr_type] = attr_value;
            
            buffer.clear();
        }
        
        return parsed_attributes;
    }
    
    GFFReader::GFFReader(istream& in) : in(in) {
        
    }
    
    void GFFReader::for_each_gff_record(function<void(const GFFRecord&)>& lambda) {
        
        while (in.good()) {
            // skip header lines
            if (in.peek() == '#') {
                in.ignore(numeric_limits<streamsize>::max(), '\n');
                continue;
            }
            
            GFFRecord record;
            
            string buffer;
            char* ignored;
            
            // parse sequence ID
            getline(in, buffer, '\t');
            if (buffer.empty()) {
                continue;
            }
            else if (buffer != ".") {
                record.sequence_id = std::move(buffer);
            }
            buffer.clear();
            
            // parse data source
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.source = std::move(buffer);
            }
            buffer.clear();
            
            // parse type of annotation
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.type = std::move(buffer);
            }
            buffer.clear();
            
            // parse start coordinate
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.start = strtol(buffer.c_str(), &ignored, 10) - 1;
            }
            buffer.clear();
            
            // parse end coordinate
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.end = strtol(buffer.c_str(), &ignored, 10) - 1;
            }
            buffer.clear();
            
            // parse score
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.score = strtod(buffer.c_str(), &ignored);
            }
            buffer.clear();
            
            // parse strand
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.strand_is_rev = (buffer == "-");
            }
            buffer.clear();
            
            // parse phase
            getline(in, buffer, '\t');
            if (buffer != ".") {
                record.phase = stoi(buffer);
            }
            buffer.clear();
            
            // parse annotations (but leave as an unparsed string)
            getline(in, buffer, '\n');
            if (buffer != ".") {
                record.attributes = std::move(buffer);
            }
            
            // execute the iteratee
            lambda(record);
        }
    }
        
}
