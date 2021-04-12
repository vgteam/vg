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
// Genotype new SVs from a GAM, and mark them as T/N for somatic calling
// void agument_SVs(VG* graph,
//                  string gamfile,
//                  string aug_prefix,
//                  bool do_split,
//                  bool do_paired);
// void augment_svs(VG* graph, 
//                    vector<string> gamfiles,
//                    vector<string> aug_prefixes,
//                    bool do_split_remapping,
//                    boool do_paired_ends);
// void somatic_call(VG* augmented, string refpath, vector<pair<string, string> > som_samples);

// void somatic_call(VG* augmented, vector<pair<string, string> > som_samples);

}

#endif
