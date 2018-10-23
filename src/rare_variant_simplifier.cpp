#include "rare_variant_simplifier.hpp"

namespace vg {

using namespace std;

RareVariantSimplifier::RareVariantSimplifier(VG& graph, VcfBuffer& variant_source) : Progressive(), graph(graph), variant_source(variant_source) {
    // Nothing to do!
}

void RareVariantSimplifier::simplify() {
    // This holds the IDs of all the nodes we want to keep around
    unordered_set<id_t> to_keep;

    // Mark the entire reference path(s) as to-keep
    

    // For each variant

    // If it is sufficiently common, mark all its alt path nodes as to-keep

    // Otherwise delete all its alt paths and als its ref path

    // After going through all the variants, delete all nodes that aren't to-keep
}

}
