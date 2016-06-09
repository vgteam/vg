#include "utility.hpp"

namespace vg {

static const char complement[256] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 16
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 24
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 32
                                     'N', 'N', 'N', '$', '#', 'N', 'N', 'N', // 40 GCSA stop/start characters
                                     'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', // 48
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 56
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 64
                                     'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', // 72
                                     'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', // 80
                                     'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S', // 88
                                     'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N', // 96
                                     'N', 't', 'v', 'g', 'h', 'N', 'N', 'c', // 104
                                     'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N', // 112
                                     'N', 'q', 'y', 'w', 'a', 'a', 'b', 's', // 120
                                     'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N', // 128
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 136
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 144
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 152
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 160
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 168
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 176
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 184
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 192
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 200
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 208
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 216
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 224
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 232
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 240
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 248
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};// 256
    
char reverse_complement(const char& c) {
    return complement[c];
}

string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        c = complement[c];
    }
    return rc;
}

int get_thread_count(void) {
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    return thread_count;
}

std::vector<std::string> &split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
    char* tok;
    char cchars [s.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, s.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }
    return elems;
}
std::vector<std::string> split_delims(const std::string &s, const std::string& delims) {
    std::vector<std::string> elems;
    return split_delims(s, delims, elems);
}

const std::string sha1sum(const std::string& data) {
    SHA1 checksum;
    checksum.update(data);
    return checksum.final();
}

const std::string sha1head(const std::string& data, size_t head) {
    return sha1sum(data).substr(0, head);
}

string wrap_text(const string& str, size_t width) {
    stringstream w;
    size_t j = 0;
    for (auto c : str) {
        if (j++ > 50) {
            if (c == ' ') {
                w << "\n";
                j = 0;
            } else {
                w << c;
            }
        } else {
            w << c;
        }
    }
    return w.str();
}

bool is_number(const std::string& s) {
    return !s.empty()
        && std::find_if(s.begin(), s.end(),
                        [](char c) { return !std::isdigit(c); }) == s.end();
}

bool allATGC(const string& s) {
    for (string::const_iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

string nonATGCNtoN(const string& s) {
    auto n = s;
    for (string::iterator c = n.begin(); c != n.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C' && b != 'N') {
            *c = 'N';
        }
    }
    return n;
}

void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar) {
    for (const auto& edit : mapping.edit()) {
        if (edit.from_length() && edit.from_length() == edit.to_length()) {
// *matches* from_length == to_length, or from_length > 0 and offset unset
            // match state
            cigar.push_back(make_pair(edit.from_length(), 'M'));
        } else {
            // mismatch/sub state
// *snps* from_length == to_length; sequence = alt
            if (edit.from_length() == edit.to_length()) {
                cigar.push_back(make_pair(edit.from_length(), 'M'));
            } else if (edit.from_length() == 0 && edit.sequence().empty()) {
// *skip* from_length == 0, to_length > 0; implies "soft clip" or sequence skip
                cigar.push_back(make_pair(edit.to_length(), 'S'));
            } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
                int32_t del = edit.from_length() - edit.to_length();
                int32_t eq = edit.to_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(del, 'D'));
            } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
                int32_t ins = edit.to_length() - edit.from_length();
                int32_t eq = edit.from_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(ins, 'I'));
            }
        }
    }
}

string mapping_string(const string& source, const Mapping& mapping) {
    string result;
    int p = mapping.position().offset();
    for (const auto& edit : mapping.edit()) {
        // mismatch/sub state
// *matches* from_length == to_length, or from_length > 0 and offset unset
// *snps* from_length == to_length; sequence = alt
        // mismatch/sub state
        if (edit.from_length() == edit.to_length()) {
            if (!edit.sequence().empty()) {
                result += edit.sequence();
            } else {
                result += source.substr(p, edit.from_length());
            }
            p += edit.from_length();
        } else if (edit.from_length() == 0 && edit.sequence().empty()) {
// *skip* from_length == 0, to_length > 0; implies "soft clip" or sequence skip
            //cigar.push_back(make_pair(edit.to_length(), 'S'));
        } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
            result += edit.sequence();
            p += edit.from_length();
        } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
            result += edit.sequence();
            p += edit.from_length();
        }
    }
    return result;
}

string cigar_string(vector<pair<int, char> >& cigar) {
    vector<pair<int, char> > cigar_comp;
    pair<int, char> cur = make_pair(0, '\0');
    for (auto& e : cigar) {
        if (cur == make_pair(0, '\0')) {
            cur = e;
        } else {
            if (cur.second == e.second) {
                cur.first += e.first;
            } else {
                cigar_comp.push_back(cur);
                cur = e;
            }
        }
    }
    cigar_comp.push_back(cur);
    stringstream cigarss;
    for (auto& e : cigar_comp) {
        cigarss << e.first << e.second;
    }
    return cigarss.str();
}

string tmpfilename(const string& base) {
    string tmpname = base + "XXXXXXXX";
    // hack to use mkstemp to get us a safe temporary file name
    int fd = mkstemp(&tmpname[0]);
    if(fd != -1) {
        // we don't leave it open; we are assumed to open it again externally
        close(fd);
    } else {
        cerr << "[vg utility.cpp]: couldn't create temp file on base "
             << base << " : " << tmpname << endl;
        exit(1);
    }
    return tmpname;
}

string get_or_make_variant_id(vcflib::Variant variant) {

     if(!variant.id.empty() && variant.id != ".") {
        // We assume all the actually filled in ID fields in a VCF are unique.
        return variant.id;
    } else {
        // Synthesize a name for the variant
        
        // Let's just hash
        SHA1 hasher;
        
        // Turn the variant back into a string line and hash it.
        // Note that this keeps the modified 0-based position.
        std::stringstream variant_stringer;
        variant_stringer << variant;
        hasher.update(variant_stringer.str());
        
        // Name the variant with the hex hash. Will be unique unless two
        // identical variant lines are in the file.
        return hasher.final();
        
    }
}

double median(std::vector<int> &v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    int vn = v[n];
    if (v.size()%2 == 1) {
        return vn;
    } else {
        std::nth_element(v.begin(), v.begin()+n-1, v.end());
        return 0.5*(vn+v[n-1]);
    }
}

}
