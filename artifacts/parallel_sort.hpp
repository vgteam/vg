#ifndef VG_PARALLEL_SORT_HPP_INCLUDED
#define VG_PARALLEL_SORT_HPP_INCLUDED

/**
 * \file parallel_sort.hpp
 * A sort that runs on several threads where the standard library has a parallel sort.
 */

#include <algorithm>

// GCC's standard library (libstdc++) has a parallel sort; others, such as clang's libc++ on macOS,
// do not.
#if defined(__GLIBCXX__) && defined(_OPENMP)
#include <parallel/algorithm>
#endif

namespace vg {

/// Sort [begin, end) by `less`, on several threads with libstdc++ and OpenMP and on one thread
/// otherwise. Neither sort is stable, so the two can put elements that `less` does not order in
/// different orders: `less` should order every two elements that are used differently.
template<typename Iterator, typename Less>
void parallel_sort(Iterator begin, Iterator end, const Less& less) {
#if defined(__GLIBCXX__) && defined(_OPENMP)
    __gnu_parallel::sort(begin, end, less);
#else
    std::sort(begin, end, less);
#endif
}

}

#endif
