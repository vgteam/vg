/**
 * \file pad_band.cpp
 * 
 * Defines implementation for band padding functions for banded global alignment.
 */

#include "pad_band.hpp"

#include <cmath>

namespace vg {
namespace algorithms {


std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_random_walk_internal(double band_padding_multiplier, size_t band_padding_memo_size, size_t max_padding,
                                                                                          const std::function<size_t(const Alignment&, const HandleGraph&)> size_function) {
    // We're goign to capture this vector by value into the closure
    std::vector<size_t> band_padding_memo;
    
    // Fill it in to initialize
    band_padding_memo.resize(band_padding_memo_size);
    for (size_t i = 0; i < band_padding_memo.size(); i++) {
        band_padding_memo[i] = std::min<size_t>(max_padding, size_t(band_padding_multiplier * sqrt(i)) + 1);
    }
    
    function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding =
        [band_padding_multiplier, band_padding_memo, size_function, max_padding](const Alignment& seq, const HandleGraph& graph) {
            size_t size = size_function(seq, graph);
            return size < band_padding_memo.size() ? band_padding_memo.at(size)
                                                   : std::min<size_t>(max_padding, size_t(band_padding_multiplier * sqrt(size)) + 1);
    };
    
    // And return the closure which now owns the memo table.
    return choose_band_padding;
}

std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_random_walk(double band_padding_multiplier, size_t band_padding_memo_size, size_t max_padding) {
    return pad_band_random_walk_internal(band_padding_multiplier, band_padding_memo_size, max_padding,
                                        [](const Alignment& seq, const HandleGraph&) {
        return seq.sequence().size();
    });
}

std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_min_random_walk(double band_padding_multiplier, size_t band_padding_memo_size, size_t max_padding) {
    return pad_band_random_walk_internal(band_padding_multiplier, band_padding_memo_size, max_padding,
                                         [](const Alignment& seq, const HandleGraph& graph) {
        return std::min(seq.sequence().size(), graph.get_total_length());
    });
}

std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_constant(size_t band_padding) {
    // don't dynamically choose band padding, shim constant value into a function type
    function<size_t(const Alignment&,const HandleGraph&)> constant_padding = [band_padding](const Alignment& seq, const HandleGraph& graph) {
        return band_padding;
    };

    return constant_padding;
}

}
}
