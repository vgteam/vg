#ifndef SNARLMAP_HPP
#define SNARLMAP_HPP

#include "algorithms/min_cut_graph.hpp"
#include "snarls.hpp"
#include <structures/union_find.hpp>
#include "sparse_union_find.hpp"
#include <unordered_map>
#include "phased_genome.hpp"
#include <vg/vg.pb.h>

namespace vg{
    using namespace std;
    using vg::algorithms::Graph;

class SnarlMap{
       
    
    
public:
    SnarlMap(SnarlManager& snarls, const vector<multipath_alignment_t>& reads, unique_ptr<PhasedGenome>& phased_genome);

    unordred_map<pair<size_t, size_t>, size_t> make_snarl_map();

    Graph make_snarl_graph(unordred_map<pair<size_t, size_t>, size_t> map);

private:
SnarlManager& snarls;
const vector<multipath_alignment_t>& reads;
unique_ptr<PhasedGenome>& phased_genome;
// maps snarl pair to edge weight that overlaps both snarls
// snarls numbers are paired in ascending order 
unordred_map<pair<size_t, size_t>, size_t> map;

    
        
};


}