#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <cstring>
#include <algorithm>
#include <unistd.h>
#include "vg.pb.h"
#include "sha1.hpp"
#include "Variant.h"

namespace vg {

using namespace std;

char reverse_complement(const char& c);
string reverse_complement(const string& seq);
int get_thread_count(void);
string wrap_text(const string& str, size_t width);

// split a string on any character found in the string of delimiters (delims)
std::vector<std::string>& split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string> split_delims(const std::string &s, const std::string& delims);

const std::string sha1sum(const std::string& data);
const std::string sha1head(const std::string& data, size_t head);

bool allATGC(string& s);
void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);

template<typename T, typename V>
set<T> map_keys_to_set(const map<T, V>& m) {
    set<T> r;
    for (auto p : m) r.insert(p.first);
    return r;
}

// pairwise maximum
template<typename T>
vector<T> pmax(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c;
    assert(a.size() == b.size());
    c.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(c),
                   [](T a, T b) { return std::max<T>(a, b); });
    return c;
}

// maximum of all vectors
template<typename T>
vector<T> vpmax(const std::vector<std::vector<T>>& vv) {
    std::vector<T> c;
    if (vv.empty()) return c;
    c = vv.front();
    typename std::vector<std::vector<T> >::const_iterator v = vv.begin();
    ++v; // skip the first element
    for ( ; v != vv.end(); ++v) {
        c = pmax(c, *v);
    }
    return c;
}

string tmpfilename(const string& base);

// Code to detect if a variant lacks an ID and give it a unique but repeatable
// one.
string get_or_make_variant_id(vcflib::Variant variant);

}

#endif
