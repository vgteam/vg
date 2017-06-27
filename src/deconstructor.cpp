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
        vector<SnarlTraversal> ordered_traversals;
        bool hasRef = false;

        bool normalize_indels = false;

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
                
                // Get the middle of the traversal that doesn't include the
                // boundary nodes
                auto iter = t.visits().begin();
                iter++;
                auto end = t.visits().end();
                end--;
                for (; iter != end; iter++){
                    auto v = *iter;
                    if (!pind->path_contains_node(v.node_id())){
                        is_ref = false;
                    }
                    if (v.node_id() == 0){
                        continue;
                    }
                    t_allele << graph->get_node(v.node_id())->sequence();
                }

                string t_str = t_allele.str();
                if (t_str == ""){
                    normalize_indels = true;
                }
                if (is_ref){
                    ret.insert(ret.begin(), t_str);
                    ordered_traversals.insert(ordered_traversals.begin(), t);
                    hasRef = true;
                }
                else{
                    ret.push_back(t_str);
                    ordered_traversals.push_back(t);
                }
                
            }

            else{
                // All alleles are alt alleles
                // Just make our strings and push them back.
                
                // Get the middle of the traversal that doesn't include the
                // boundary nodes
                auto iter = t.visits().begin();
                iter++;
                auto end = t.visits().end();
                end--;
                for (; iter != end; iter++){
                    auto v = *iter;
                    t_allele << graph->get_node(v.node_id())->sequence();
                }
                ret.push_back(t_allele.str());
                ordered_traversals.push_back(t);
            }
        }
            // If we haev indels to normalize, loop over our alleles
            // normalize each string to VCF-friendly format (i.e. clip one ref base
            // on the left side and put it in the ref field and the alt field).
            if (normalize_indels){
                for (int i = 0; i < ret.size(); ++i){
                // Get the reference base to the left of the variant.
                // If our empty allele is the reference (and we have a reference),
                // put our new-found ref base in the 0th index of alleles vector.
                // Then, prepend that base to each allele in our alleles vector.
                    SnarlTraversal t = ordered_traversals[i];
                    id_t start_id = t.visits(0).node_id();
                    id_t end_id = t.visits(t.visits_size() - 1).node_id();
                    pair<size_t, bool> pos_orientation_start = pindexes[refpath]->by_id[start_id];
                    pair<size_t, bool> pos_orientation_end = pindexes[refpath]->by_id[end_id];
                    bool use_start = pos_orientation_start.first < pos_orientation_end.first;
                    bool rev = use_start ? pos_orientation_start.second : pos_orientation_end.second;
                    string pre_node_seq = use_start ? graph->get_node(start_id)->sequence() :
                                            graph->get_node(end_id)->sequence();
                    string pre_variant_base = rev ? string(1, pre_node_seq[0]) : string(1, pre_node_seq[pre_node_seq.length() - 1]);
                    ret[i].insert(0, pre_variant_base);
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
            bool use_start = pos_orientation_start.first < pos_orientation_end.first;
            size_t node_pos = (use_start ? pos_orientation_start.first : pos_orientation_end.first);
            v.position = node_pos +(use_start ? graph->get_node(snarl->start().node_id())->sequence().length() : graph->get_node(snarl->end().node_id())->sequence().length());
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

