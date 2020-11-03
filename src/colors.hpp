#ifndef VG_COLORS_HPP_INCLUDED
#define VG_COLORS_HPP_INCLUDED

/**
 * \file colors.hpp
 * Includes facilities for generating random GraphViz colors seeded by strings.
 * Used for colorizing paths in visualizations.
 */

#include <vector>
#include <random>

namespace vg {

using namespace std;

/**
 * Generates random GraphViz colors seeded by strings.
 * Used for colorizing paths in visualizations.
 */
class Colors {
    mt19937 rng;
public:
    const vector<string> colors = { "/dark28/1", "/dark28/2", "/dark28/3", "/dark28/4", "/dark28/5", "/dark28/6", "/dark28/7", "/dark28/8" };
    Colors(void) { };
    Colors(int seed_val) {
        rng.seed(seed_val);
    };
    ~Colors(void) { };
    string hashed(const string& str) {
        std::hash<std::string> hash_fn;
        std::size_t str_hash = hash_fn(str);
        size_t i = str_hash % colors.size();
        return colors[i];
    }
    string random(void) {
        uniform_int_distribution<int> dist(0, 7);
        return colors[dist(rng)];
    }
};

}

#endif

