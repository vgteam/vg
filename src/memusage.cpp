#include "memusage.hpp"

#include <istream>
#include <fstream>
#include <sstream>

#include <sys/time.h>
#include <sys/resource.h>

namespace vg {

using namespace std;

string get_proc_status_value(const string& name) {

    ifstream status_file("/proc/self/status");
    
    string line;
    while (status_file.good()) {
        // Grab each line
        getline(status_file, line);
        
        // Find the first colon
        size_t colon_pos = line.find(':');
        if (colon_pos == string::npos) {
            // No colon found. Try the next line.
            continue;
        }
        
        if (line.substr(0, colon_pos) != name) {
            // This isn't what we care about
            continue;
        }
        
        if (line.size() == colon_pos + 1) {
            // There's nothing after the colon
            continue;
        }
        
        // Find the first non-whitespace after the colon
        size_t value_pos = line.find_first_not_of(" \t", colon_pos + 1);
        if (value_pos == string::npos) {
            // No value
            continue;
        }
        
        // Get the value text and return it
        return line.substr(value_pos);
    }
    
    return "";
    
}

size_t get_max_rss_kb() {
    // This isn't in /proc, we have to get it ourselves
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    return usage.ru_maxrss;
}


size_t get_max_vmem_kb() {
    string value = get_proc_status_value("VmPeak");
    
    if (value == "") {
        return 0;
    }
    
    stringstream sstream(value);
    
    size_t result = 0;
    
    sstream >> result;
    
    return result;
    
}

size_t get_current_vmem_kb() {
    string value = get_proc_status_value("VmSize");
    
    if (value == "") {
        return 0;
    }
    
    stringstream sstream(value);
    
    size_t result = 0;
    
    sstream >> result;
    
    return result;
}


}
