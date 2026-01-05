/**
 * \file interval_union.cpp: contains the implementation of IntervalUnion
 */


#include "interval_union.hpp"

#include <cassert>
#include <iostream>

namespace vg {

using namespace std;

void IntervalUnion::add(size_t begin, size_t end) {

    assert(begin <= end);
    if (begin == end) {
        // empty interval, can be ignore
        return;
    }

    auto after = irvl_union.upper_bound(begin);
    auto merge_from = irvl_union.end();
    if (after != irvl_union.begin()) {
        // get the nearest interval to the left, possibly at the same index
        auto before = after;
        --before;
        if (before->second >= begin) {
            // merge this interval into the earlier interval
            if (end > before->second) {
                union_size += end - before->second;
                before->second = end;
                merge_from = before;
            }
        }
        else {
            // the start of this interval is outside the one to the left
            union_size += (end - begin);
            merge_from = irvl_union.emplace_hint(after, begin, end);
        }
    }
    else {
        // there is no interval to the left
        union_size += (end - begin);
        merge_from = irvl_union.emplace_hint(after, begin, end);
    }

    // TODO: it might also be possible to do a lazy-update version of this data structure
    // to get amortized O(log n) add

    if (merge_from != irvl_union.end()) {
        // try to merge this interval with the ones to the right
        auto next = merge_from;
        ++next;
        while (next != irvl_union.end() && next->first <= merge_from->second) {
            union_size -= (min(merge_from->second, next->second) - next->first);
            merge_from->second = max(merge_from->second, next->second);
            irvl_union.erase(next);
            next = merge_from;
            ++next;
        }
    }
}

void IntervalUnion::clear()  {
    union_size = 0;
    irvl_union.clear();
}

size_t IntervalUnion::overlap(size_t begin, size_t end) const {

    // start to the left of the query interval
    auto iter = irvl_union.upper_bound(begin);
    if (iter != irvl_union.begin()) {
        --iter;
    }

    // TODO: there's a worst-case O(log n) solution using a sum-augmented range tree
    // to get the contribution from strictly contained intervals, but rebalancing it
    // would be a pain

    // add up the length of the intervals that this overlaps
    size_t total_overlap = 0;
    while (iter != irvl_union.end() && iter->first < end) {
        if (iter->second > begin && end > iter->first) {
            total_overlap += (min(end, iter->second) - max(begin, iter->first));
        }
        ++iter;
    }

    return total_overlap;
}

size_t IntervalUnion::total_size() const {
    return union_size;
}

size_t IntervalUnion::component_size() const {
    return irvl_union.size();
}

vector<pair<size_t, size_t>> IntervalUnion::get_union() const {
    vector<pair<size_t, size_t>> intervals;
    intervals.reserve(irvl_union.size());
    for (const auto& interval : irvl_union) {
        intervals.push_back(interval);
    }
    return intervals;
}

}

