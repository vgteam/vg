#include "entropy.hpp"

namespace vg {

using namespace std;

double entropy(const string& st) {
    return entropy(st.c_str(), st.length());
}

double entropy(const char* st, size_t len) {
    map<char, int> freqs;
    for (size_t i = 0; i < len; ++i) {
        ++freqs[st[i]];
    }
    double ent = 0;
    double ln2 = log(2);
    for (auto& f : freqs) {
        double freq = (double)f.second/len;
        ent += freq * log(freq)/ln2;
    }
    ent = -ent;
    return ent;
}

}
