#include "deconstructor.hpp"

using namespace std;


namespace vg {
    Deconstructor::Deconstructor(){

    }
    Deconstructor::~Deconstructor(){
    }

    string Deconstructor::vcf_header(string refpath, vg::VG* graph){
        stringstream ret;

        return ret.str();
    }

    // string Deconstructor::get_ref_allele(vector<SnarlTraversal> alleles, string refpath, vg::VG* graph){

    // }

    void Deconstructor::deconstruct(string refpath, vg::VG* graph){
        PathIndex pp(*graph, refpath);
        // Load graph
        // Find snarls
        // Spit header
        // For each snarl
        //    enumerate refs / alts
        //    convert to VCF
        //    spit out
        SnarlFinder* snarl_finder = new CactusUltrabubbleFinder(*graph, refpath, true);
        SnarlManager snarl_manager = snarl_finder->find_snarls();
        vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
        SimpleConsistencyCalculator scc;
        TraversalFinder* trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);
        for (const Snarl* snarl : snarl_roots ){
            vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
            // Pack snarl traversals into allele fields
            // get ref position
            // produce variant
        }
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
            TraversalFinder* trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);
            for (const Snarl* snarl : snarl_roots ){
                vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
                
            }
        }


        // For each snarl
        //    enumerate refs / alts
        //    convert to VCF
        //    spit out
    }
}

