/**
 * \file vg_gaf.hpp
 *
 * GAM <==> GAF conversion
 */

#pragma once

#include <vg/vg.pb.h>
#include "gafkluge.hpp"
#include "handle.hpp"
#include "alignment.hpp"
#include "json2pb.h"

namespace vg {

/// Alignemnt -> GAF
gafkluge::GafRecord aln2gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar = true, bool base_quals = true);

/// GAF -> Alignment
Alignment gaf2aln(const HandleGraph& graph, const gafkluge::GafRecord& gaf);

/// Stream a GAF from a file.  Can be bgzipped
int gaf_for_each(const string& filename, function<void(Alignment&)> lambda,
                 const HandleGraph& graph);

/// Stream a GAF from a file... in parallel.  Can be bgzipped
int gaf_for_each_parallel(const string& filename, function<void(Alignment&)> lambda,
                          const HandleGraph& graph);

/// Stream a GAF from a file ... in pairs... in parallel.  Can be bgzipped
int gaf_for_each_itnerleaved_pair_parallel(const string& filename, function<void(Alignment&, Alignment&)> lambda2,
                                           const HandleGraph& graph);

}
