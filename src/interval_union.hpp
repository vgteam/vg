/** \file
 * interval_union.hpp: defines the union of a collection of intervals that can
 * be constructed incrementally
 */
#ifndef VG_INTERVAL_UNION_HPP_INCLUDED
#define VG_INTERVAL_UNION_HPP_INCLUDED

#include <map>
#include <vector>
#include <utility>

namespace vg {

using namespace std;

/**
 * The union of a collection of intervals
 */
class IntervalUnion
{
public:
    IntervalUnion() = default;
    ~IntervalUnion() = default;

    // add a new interval to the union
    void add(size_t begin, size_t end);
    inline void add(const pair<size_t, size_t>& interval);

    // return the union to an empy starting state
    void clear();

    // query an interval's length of overlap with the union
    size_t overlap(size_t begin, size_t end) const;
    inline size_t overlap(const pair<size_t, size_t>& interval) const;

    // get the total size of the intervals in the union
    size_t total_size() const;

    // return the number of disjoint intervals that comprise the union
    size_t component_size() const;

    // get the intervals that comprise the union
    vector<pair<size_t, size_t>> get_union() const;

private:
    
    map<size_t, size_t> irvl_union;

    size_t union_size = 0;
     
};

/*
 * Inline implementations
 */

void IntervalUnion::add(const pair<size_t, size_t>& interval) {
    add(interval.first, interval.second);
}

size_t IntervalUnion::overlap(const pair<size_t, size_t>& interval) const {
    return overlap(interval.first, interval.second);
}

}

#endif // VG_INTERVAL_UNION_HPP_INCLUDED
