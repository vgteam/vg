#ifndef VG_ALGORITHMS_BACK_TRANSLATE_HPP_INCLUDED
#define VG_ALGORITHMS_BACK_TRANSLATE_HPP_INCLUDED

/**
 * \file back_translate.hpp
 *
 * Defines an algorithm for translating an Alignment from node ID space to named segment space. 
 */

#include "../handle.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace algorithms {
using namespace std;

/**
 * Translate the given Path in place from node ID space to named segment space.
 */
void back_translate_in_place(const NamedNodeBackTranslation* translation, Path& path);

/**
 * Translate the given Snarl in place from node ID space to named segment space.
 */
void back_translate_in_place(const NamedNodeBackTranslation* translation, Snarl& snarl);

/**
 * Translate the given SnarlTraversal in place from node ID space to named
 * segment space. Visits that end up being to snarls where both boundaries are
 * from the same orientation of the same segment will be removed. Multiple
 * visits in a row to the same orientation of the same segment will be elided.
 *
 * TODO: Work out a way to preserve traversals around cycles while not
 * repeating ourselves for visits to diffrent pieces of the same chopped
 * segment.
 */
void back_translate_in_place(const NamedNodeBackTranslation* translation, SnarlTraversal& traversal);

/**
 * Translate the given Visit in place from node ID space to named segment space.
 */
void back_translate_in_place(const NamedNodeBackTranslation* translation, Visit& visit);

}
}

#endif
