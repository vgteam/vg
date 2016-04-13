#ifndef VG_SIMULATOR_H
#define VG_SIMULATOR_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "json2pb.h"

namespace vg {

using namespace std;

class Sampler {

public:

    xg::XG* xgidx;
    mt19937 rng;
    Sampler(xg::XG* x, int seed = 0) : xgidx(x) {
        if (!seed) {
            seed = time(NULL);
        }
        rng.seed(seed);
    }

    pos_t position(void);
    string sequence(size_t length);
    Alignment alignment(size_t length);
    char pos_char(pos_t pos);
    map<pos_t, char> next_pos_chars(pos_t pos);

    Alignment mutate(const Alignment& aln,
                     double base_error,
                     double indel_error);

    Edit mutate_edit(const Edit& edit,
                     const Mapping& mapping,
                     double base_error,
                     double indel_error,
                     const string& bases,
                     uniform_real_distribution<double>& rprob,
                     uniform_int_distribution<int>& rbase);

};

}

#endif
