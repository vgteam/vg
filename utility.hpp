#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <omp.h>

namespace vg {

using namespace std;

string reverse_complement(const string& seq);
int get_thread_count(void);

}

#endif
