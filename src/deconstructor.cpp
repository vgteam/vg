#include "deconstructor.hpp"

using namespace std;


namespace vg {
    Deconstructor::Deconstructor(){

    }
    Deconstructor::~Deconstructor(){
    }

    /**
    * Takes in a vector of snarltraversals
    * returns their sequences as a vector<string>
    * returns a boolean hasRef
    * if a reference path is present, hasRef is set to true and the first
    * string in the vector is the reference allele
    * otherwise, hasRef is set to false and all strings are alt alleles.
    */
    pair<bool, vector<string> > Deconstructor::get_alleles(vector<SnarlTraversal> travs, string refpath, vg::VG* graph){
        vector<string> ret;
        bool hasRef = false;

        // Check if we have a PathIndex for this path
        bool path_indexed = pindexes.find(refpath) != pindexes.end();
        for (auto t : travs){
            stringstream t_allele;

            // Get ref path index
            if (path_indexed){
                PathIndex* pind = pindexes[refpath];
                // Check nodes of traversals Visits
                // if they're all on the ref path,
                // then this Snarltraversal is the ref allele.
                bool is_ref = true;
                for (auto v : t.visits()){
                    if (!pind->path_contains_node(v.node_id())){
                        is_ref = false;
                    }
                    t_allele << graph->get_node(v.node_id())->sequence();
                }
                if (is_ref){
                    ret.insert(ret.begin(), t_allele.str());
                    hasRef = true;
                }
                else{
                    ret.push_back(t_allele.str());
                }

            }

            else{
                // All alleles are alt alleles
                // Just make our strings and push them back.
                for (auto v : t.visits()){
                    t_allele << graph->get_node(v.node_id())->sequence();
                }
                ret.push_back(t_allele.str());
            }
            
        }
        return make_pair(hasRef, ret);

    }

    void Deconstructor::deconstruct(string refpath, vg::VG* graph){
        
     

        // Create path index for the contig if we don't have one.
        if (pindexes.find(refpath) == pindexes.end()){
            pindexes[refpath] = new PathIndex(*graph, refpath, false);
        }

        // Spit header
        // Set contig to refpath
        // Set program field
        // Set the version and data
        // Set info field, if needed
        // Make the header line
        // open a VCF file

        if (!headered){
            vcflib::VariantCallFile outvcf;
            stringstream stream;
            stream << "##fileformat=VCFv4.2" << endl;
            stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << "Sample" << endl;
            
            string hstr = stream.str();
            assert(outvcf.openForOutput(hstr));
            cout << outvcf.header << endl;
            this->headered = true;
        }

        // Find snarls
        // Snarls are variant sites ("bubbles")
        SnarlFinder* snarl_finder = new CactusUltrabubbleFinder(*graph, refpath, true);
        SnarlManager snarl_manager = snarl_finder->find_snarls();
        vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
        SimpleConsistencyCalculator scc;
        TraversalFinder* trav_finder = new ExhaustiveTraversalFinder(*graph, snarl_manager);
        for (const Snarl* snarl: snarl_roots){
            vcflib::Variant v;
            // SnarlTraversals are the (possible) alleles of our variant site.
            vector<SnarlTraversal> travs = trav_finder->find_traversals(*snarl);
            // write variant's sequenceName (VCF contig)
            v.sequenceName = refpath;
            // Set position based on the lowest position in the snarl.
            pair<size_t, bool> pos_orientation_start = pindexes[refpath]->by_id[snarl->start().node_id()];
            pair<size_t, bool> pos_orientation_end = pindexes[refpath]->by_id[snarl->end().node_id()];
            
            v.position = (pos_orientation_start.first < pos_orientation_end.first) ? pos_orientation_start.first : pos_orientation_end.first;
            std::pair<bool, vector<string> > t_alleles = get_alleles(travs, refpath, graph);
            if (t_alleles.first){
                v.alleles.insert(v.alleles.begin(), t_alleles.second[0]);
                v.ref = t_alleles.second[0];
                for (int i = 1; i < t_alleles.second.size(); i++){
                    v.alleles.push_back(t_alleles.second[i]);
                    v.alt.push_back(t_alleles.second[i]);
                }
            }
            else{
                cerr << "NO REFERENCE ALLELE FOUND" << endl;
                v.alleles.insert(v.alleles.begin(), ".");
                for (int i = 0; i < t_alleles.second.size(); i++){
                    v.alleles.push_back(t_alleles.second[i]);
                    v.alt.push_back(t_alleles.second[i]);
                }
            }
            v.updateAlleleIndexes();
            cerr << v << endl;

        }
        

    }

    /**
    * Convenience wrapper function for deconstruction of multiple paths.
    */
    void Deconstructor::deconstruct(vector<string> ref_paths, vg::VG* graph){

        for (auto path : ref_paths){
            deconstruct(path, graph);
        }

    }
}

