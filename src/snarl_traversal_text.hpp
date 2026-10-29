#ifndef VG_SNARL_TRAVERSAL_TEXT_HPP_INCLUDED
#define VG_SNARL_TRAVERSAL_TEXT_HPP_INCLUDED

/** \file
 * The text of a snarl and of a traversal through one, as the VCF callers write them.
 */

#include <string>

#include "handle.hpp"
#include "snarls.hpp"
#include "../artifacts/vcf_record.hpp"

namespace vg {

using namespace std;

/// Return the handle of a visit to a node: the node, oriented as the visit reads it.
handle_t visit_handle(const HandleGraph& graph, const Visit& visit);

/// Return the name of a snarl, as `site_name` gives it from the snarl's boundary visits: ">3435<12222",
/// or "(>3435<12222)" with `in_brackets`. Nodes are named by `translation` if it is not null.
string snarl_name(const Snarl& snarl, const NodeTranslation* translation, bool in_brackets = false);

/// Return the text of a traversal: the sequence of each node it visits, in the orientation of the
/// visit, and the name of each child snarl it visits, in brackets (see `snarl_name`).
string trav_string(const HandleGraph& graph, const SnarlTraversal& trav,
                   const NodeTranslation* translation);

/// Return the traversal of `travs` that each of `numbers` gives, or an empty traversal for a
/// negative number, which stands for the "*" allele.
vector<SnarlTraversal> traversals_of(const vector<SnarlTraversal>& travs,
                                     const vector<int>& numbers);

}

#endif
