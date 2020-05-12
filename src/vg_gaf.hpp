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
gafkluge::GafRecord aln2gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar = true);

/// GAF -> Alignment
Alignment gaf2aln(const HandleGraph& graph, const gafkluge::GafRecord& gaf);

}
