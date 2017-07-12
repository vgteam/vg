#ifndef VG_ENTROPY_HPP_INCLUDED
#define VG_ENTROPY_HPP_INCLUDED

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <map>

namespace vg {

using namespace std;

double entropy(const string& st);
double entropy(const char* st, size_t len);

}

#endif
