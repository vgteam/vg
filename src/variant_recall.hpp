#ifndef VG_VARIANT_RECALL_HPP_INCLUDED
#define VG_VARIANT_RECALL_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <vector>
#include <list>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "deconstructor.hpp"
#include "hash_map.hpp"
#include "types.hpp"
#include "genotypekit.hpp"

namespace vg {

using namespace std;


// Genotype known variants from a VCF file.
void variant_recall(VG* graph,
                    vcflib::VariantCallFile* vars,
                    FastaReference* ref_genome,
                    vector<FastaReference*> insertions,
                    string gamfile);
// Genotype new SVs from a GAM
void genotype_svs(VG* graph, 
                  string gamfile, string refpath);

}

#endif
