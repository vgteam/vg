#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?
#include <algorithm>
#include <string>
#include "../snarls.hpp"
namespace vg {

/**
 * evaluates the performance of the snarl normalizer.
 * 42794 - another snarl with the issue where the path jumps over a node of interest? Is this my fault or..? 
 * 43899
 * 71109
 * stops for several minuts at 2049699 b/c large size snarl now at 1:40
 * Stats of the new graph:number of top_level snarls in graph: 3418
 * number of total snarls in graph: 7443
 * nodes	4218503
 * edges	4227433
 * length	134506805
 */

void evaluate_normalized_snarls(ifstream &snarl_stream) {
    cerr << "evaluate_normalized_snarls" << endl;
    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);

    // Use this code to count number of snarls in graph.
    int top_count = 0;
    for (const Snarl* snarl : snarl_manager->top_level_snarls()){
        top_count++;
    }
    cerr << "number of top_level snarls in graph: " << top_count << endl;
        int general_count = 0;
    snarl_manager->for_each_snarl_preorder([&](const vg::Snarl * ignored){
        general_count++;
    });
    cerr << "number of total snarls in graph: " << general_count << endl;

}
}