#ifndef VG_ALGORITHMS_FIND_TRANSLATION_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_TRANSLATION_HPP_INCLUDED

/**
 * \file find_translation.hpp
 *
 * Defines an algorithm for finding the NamedNodeBackTranslation associated with a handle graph, if any.
 */

#include "../handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/**
 * Find the NamedNodeBackTranslation defining e.g. GFA segment space for the given handle graph, if any exists.
 * Works on GFAs that have been loaded and finds their path space information.
 */
const NamedNodeBackTranslation* find_translation(const HandleGraph* graph);

}
}

#endif
