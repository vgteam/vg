#include "deconstructor.hpp"

using namespace std;


namespace vg {
    Deconstructor::Deconstructor(){

    }
    Deconstructor::~Deconstructor(){
    }

    void Deconstructor::deconstruct(string refpath, vg::VG* graph){
        // Load graph
        // Find snarls
        // Spit header
        // For each snarl
        //    enumerate refs / alts
        //    convert to VCF
        //    spit out
    }

    void Deconstructor::deconstruct(vector<string> ref_path, vg::VG* graph){
        // for each ref path
        // Load graph
        // Find snarls
        // Spit header
        for (auto p : ref_path){
            SnarlFinder* snarl_finder = new CactusUltrabubbleFinder(*graph, p, true);
            SnarlManager snarl_manager = snarl_finder->find_snarls();
            vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
            SimpleConsistencyCalculator scc;
            TraversalFinder* trav_finder = new PathBasedTraversalFinder(*graph);
        }


        // For each snarl
        //    enumerate refs / alts
        //    convert to VCF
        //    spit out
    }
}

