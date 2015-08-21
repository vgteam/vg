#ifndef VG_COLORS_H
#define VG_COLORS_H

#include <vector>
#include <random>

namespace vg {

using namespace std;

class Colors {
    mt19937 rng;
public:
    const vector<string> colors = { "/dark28/1", "/dark28/2", "/dark28/3", "/dark28/4", "/dark28/5", "/dark28/6", "/dark28/7", "/dark28/8" };
    Colors(void) { };
    Colors(int seed_val) {
        rng.seed(seed_val);
    };
    ~Colors(void) { };
    string random(void) {
        uniform_int_distribution<int> dist(0, 7);
        return colors[dist(rng)];
    }
};

}

#endif

