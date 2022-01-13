#include <cstdint>
#include "variant_recall.hpp"
#include "cactus_snarl_finder.hpp"
#include "traversal_finder.hpp"
#include "augment.hpp"
#include "translator.hpp"
#include "xg.hpp"
#include "utility.hpp"
#include "filter.hpp"
#include "statistics.hpp"
#include "path_index.hpp"

using namespace xg;

namespace vg {

using namespace std;

/**
 // Takes a graph and two GAMs, one tumor and one normal
 // Locates existing variation supported by the tumor and annotate it with a path
 // Then overlay the normal sample
 // Use a depthmap of snarltraversal transforms, one for tumor, one for normal
 // which we can use to count the normal and tumor alleles
 void somatic_genotyper(VG* graph, string tumorgam, string normalgam);

 // Do smart augment, maintaining a depth map for tumor/normal perfect matches
 // and then editing in all of the SV reads (after normalization) with a T/N_ prefix
 // Then, get our Snarls
 // count reads supporting each and genotype
 void somatic_caller_genotyper(VG* graph, string tumorgam, string normalgam);
**/

//void smart_augment(VG* graph, string gamfile);

void genotype_svs(VG* graph,
                             string gamfile,
                             string refpath){
    // Open up our GAM file
    ifstream gamstream;
    gamstream.open(gamfile);
    if (!gamstream.good()){
        cerr << "GAM file is no good" << endl;
        exit(2);
    }
    Filter filt;

    //DepthMap dm(graph);
    vector<pair<Alignment, Alignment> > sv_reads;
    // Pull out all of our boring reads and just load them in a depth map
    // Collect our SV-supporting reads and load them into a local map for PE manipulation
    // This includes:
    // Softclipped (later split) reads
    // One end anchored (one mate is mapped)
    // Unmapped reads
    // Reads with internal mismatches (cleanly anchored reads)
    // everted pairs and split-flips
    // We trust that the relevant flags are set by FILTER
    vector<Path> direct_ins;
    set<NodeSide> spare_nodesides;
    std::function<void(Alignment&, Alignment&)> readfunc = [&](Alignment& a, Alignment& b){
        bool toss_into_sv_map = false;


        if (filt.mark_sv_alignments(a,b)){
            sv_reads.push_back(make_pair(a, b));
        }
                                
        else if (filt.mark_smallVariant_alignments(a, b)){
            direct_ins.push_back(a.path());
            direct_ins.push_back(b.path());
        }

    };
    vg::io::for_each_interleaved_pair_parallel(gamstream, readfunc);
    vector<Translation> transls;
    if (refpath != ""){
        augment(graph, direct_ins, "GAM", &transls);

        XG xg_index;
        xg_index.from_path_handle_graph(*graph); // Index the graph so deconstruct can get path positions
        Deconstructor decon;
        CactusSnarlFinder finder(xg_index);
        SnarlManager snarl_manager = finder.find_snarls();
        decon.deconstruct({refpath}, &xg_index, &snarl_manager, false, 1, false, 10000, false, false, false);
    }
    direct_ins.clear();



    // Now the weird bit
    // Transform our SV reads into clean calls
    // Merge calls
    // Then call and genotype them
    // Inversions
    // Insertions
    // Deletions
    // Duplications
    vector<BREAKPOINT> insert_bps;
    vector<BREAKPOINT> del_bps;
    vector<BREAKPOINT> inversion_bps;

    // Detect INV, DEL, INS using split-reads
    std::function<void(Alignment& first, Alignment& second)> splitreadfunc = [&](Alignment& first, Alignment& second){
                
    };

    std::function<void(Alignment& first)> se_splitreadfunc = [&](Alignment& first){
                
    };
    // Detect INV, DEL, INS using paired-ends
    // Includes local assembly using fermilite
    std::function<void(Alignment& first, Alignment& second)> pairedendfunc = [&](Alignment& first, Alignment& second){
        // grab unmapped reads
        // grab OEAs
        // desperately try to assemble them after converting to bseqs
        // Remap
        // prepare to edit
    };
}

    
/**
 * run with : vg genotype -L -V v.vcf -I i.fa -R ref.fa 
 */
void variant_recall(VG* graph,
                               vcflib::VariantCallFile* vars,
                               FastaReference* ref_genome,
                               vector<FastaReference*> insertions,
                               string gamfile){
    // Store variant->name index
    map<string, vcflib::Variant> hash_to_var;
    set<int64_t> variant_nodes;

    bool use_snarls = false;

    std::function<double(int, int)> scale_read_counts = [&](int allele_len, int read_len){
        double scale = (double) allele_len / (double) read_len;
        if (scale < 1.0){
            return 1.0;
        }
        else{
            return scale;
        }
    };

    // Store a list of node IDs each variant covers
    map<string, unordered_set<int64_t> > allele_name_to_node_id;
        
    unordered_map<int64_t, string> node_to_variant;

        

    map<string, set<string> > allele_name_to_alignment_name;
        

    // A dumb depth map, for each node in the graph, of substantial size
    unordered_map<int64_t, int32_t> node_id_to_depth;
    node_id_to_depth.reserve(10000);

    //Cache graph paths
    unordered_map<string, list<mapping_t> > gpaths( (graph->paths)._paths.begin(), (graph->paths)._paths.end()  );

    // To allow non-flat alleles, we want to use SnarlTraversals rather than
    // paths.
    unordered_map<string, set<int64_t> > snarl_name_to_node_set;
    map<string, set<string> > traversal_name_to_alignment_names;
    unordered_map<string, SnarlTraversal> name_to_traversal;
    map<string, Snarl> name_to_snarl;
    map<string, vector<SnarlTraversal> > snarl_name_to_traversals;

    if (use_snarls){
        SnarlFinder* snarl_finder = new CactusSnarlFinder(*graph);
        SnarlManager snarl_manager = snarl_finder->find_snarls();
        vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();
        SimpleConsistencyCalculator scc;
        TraversalFinder* trav_finder = new PathBasedTraversalFinder(*graph, snarl_manager);
        for (const Snarl* snarl : snarl_roots ){
            vector<SnarlTraversal> travs =  trav_finder->find_traversals(*snarl);
            snarl_name_to_traversals[ snarl->name() ] = travs;
            name_to_snarl[snarl->name()] = *snarl;
            for (auto x : travs){
                name_to_traversal[x.name()] = x;
                for (int i = 0; i < x.visit_size(); i++){
                    snarl_name_to_node_set[snarl->name()].insert( x.visit(i).node_id());
                }    
            }
        }
        cerr << "Snarls processed, and we have this many: " << snarl_name_to_node_set.size() << endl;
    }
        

    // For each variant in VCF:
    vcflib::Variant var;
    while(vars->getNextVariant(var)){
        // DON'T Adjust the position offset, canonicalize any structural variants,
        // and get the sha1 hash of the variant and store that in a map for later.
        // var.position -= 1;
        var.canonicalize(*ref_genome, insertions, true, 0);
        string var_id = make_variant_id(var);

        hash_to_var[ var_id ] = var;


        // If we're just using the GAM, build a map Node -> variant and
        // a map alt_path_id (i.e. "allele name") -> NodeID. We'll use these to count
        // mappings to nodes later.
        for (int alt_ind = 0; alt_ind <= var.alt.size(); alt_ind++){
            string alt_id = "_alt_" + var_id + "_" + std::to_string(alt_ind);
            list<mapping_t> x_path = gpaths[ alt_id ];
            for (mapping_t x_m : x_path){
                allele_name_to_node_id[ alt_id ].insert(x_m.node_id());
                node_to_variant[x_m.node_id()] = alt_id;
                variant_nodes.insert(x_m.node_id());
            }
        }

    }

    vcflib::VariantCallFile outvcf;
    stringstream stream;
    stream << "##fileformat=VCFv4.2" << endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    stream << "##FORMAT=<ID=GP,Number=1,Type=String,Description=\"Genotype Probability\">" << endl;
    stream << "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << endl;
    stream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">" << endl;
    stream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV Type\">" << endl;
    stream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"SV End\">" << endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << endl;
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << endl;
    stream << "##FORMAT=<ID=AL,Number=.,Type=Float,Description=\"Allelic likelihoods for the ref and alt alleles in the order listed\">" << endl;

    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << "Sample" << endl;

    string hstr = stream.str();
    assert(outvcf.openForOutput(hstr));
    cout << outvcf.header << endl;

    std::function<bool(const Mapping& m)> sufficient_matches = [&](const Mapping& m){
        int matches = 0;
        int tot_len = 0;
        for (int i = 0; i < m.edit_size(); ++i){
            Edit e = m.edit(i);
            if (e.to_length() == e.from_length() && e.sequence().empty()){
                matches += e.to_length();
            }
            tot_len += e.to_length();
        }
        return ( (double) matches / (double) tot_len) > 0.85;
    };

    std::function<void(Alignment&)> incr = [&](const Alignment& a){
        bool anchored = false;
        bool contained = false;
        int64_t node_for_var = 0;
        for (int i = 0; i < a.path().mapping_size(); i++){
            int64_t node_id = a.path().mapping(i).position().node_id();
            if (variant_nodes.count(node_id) && a.mapping_quality() > 20 && sufficient_matches(a.path().mapping(i))){
#pragma omp atomic write
                contained = true;
                node_for_var = node_id;

            }
            else if (!variant_nodes.count(node_id)){
#pragma omp atomic write
                anchored = true;
            }
        }

        if (contained & anchored){

#pragma omp critical
            {
                // Get our variant allele hash
                string all_str = node_to_variant[node_for_var];

                // Add alignment to the variant's set of alignments
                allele_name_to_alignment_name[all_str].insert(a.name());

            }
                
        }
    };

    std::function<void(const Alignment&)> index_incr = [&allele_name_to_alignment_name, &sufficient_matches, &node_to_variant, &variant_nodes](const Alignment& a){
        bool anchored = false;
        bool contained = false;
        int64_t node_for_var = 0;
        for (int i = 0; i < a.path().mapping_size(); i++){
            int64_t node_id = a.path().mapping(i).position().node_id();
            if (variant_nodes.count(node_id) && a.mapping_quality() > 20 && sufficient_matches(a.path().mapping(i))){
#pragma omp atomic write
                contained = true;
#pragma omp atomic write
                node_for_var = node_id;
            }
            else if (!variant_nodes.count(node_id)){
#pragma omp atomic write
                anchored = true;
            }
        }

        if (contained & anchored){
#pragma omp critical
            {
                // Get our variant allele hash
                string all_str = node_to_variant[node_for_var];

                // Add alignment to the variant's set of alignments
                allele_name_to_alignment_name[all_str].insert(a.name());

            }
                
        }
    };

    SimpleConsistencyCalculator scc;
    std::function<void(Alignment&)> count_traversal_supports = [&](const Alignment& a){
        vector<int64_t> node_list;
        for (int i = 0; i < a.path().mapping_size(); i++){
            int64_t node_id = a.path().mapping(i).position().node_id();
            if (sufficient_matches(a.path().mapping(i))){
                node_list.push_back(node_id);
            }
        }
        std::sort(node_list.begin(), node_list.end());
        for (auto n_to_s : name_to_snarl){
            std::vector<int64_t>::iterator intersection_begin;
            std::vector<int64_t>::iterator intersection_end = set_intersection(node_list.begin(), node_list.end(),
                                                                               snarl_name_to_node_set[n_to_s.first].begin(),
                                                                               snarl_name_to_node_set[n_to_s.first].end(),
                                                                               intersection_begin);
            vector<int64_t> inter (intersection_begin, intersection_end);
            if (inter.size() != 0){
                vector<bool> consistencies = scc.calculate_consistency(n_to_s.second,
                                                                       snarl_name_to_traversals[n_to_s.first], a);
                for (int c = 0; c < consistencies.size(); c++){
                        
                    if (consistencies[c]){
                        SnarlTraversal st = snarl_name_to_traversals[n_to_s.first][c];
                        traversal_name_to_alignment_names[st.name()].insert(a.name()); 
                    }
                }
                    
            } 
                
            
        }
            
            
    };
    // open our gam, count our reads, close our gam.
    if (use_snarls){
        ifstream gamstream(gamfile);
        if (gamstream.good()){
            vg::io::for_each(gamstream, count_traversal_supports);
        }
        gamstream.close();
    }
    else{
        ifstream gamstream(gamfile);
        if (gamstream.good()){
            vg::io::for_each(gamstream, incr);
        }
        else{
            cerr << "GAM stream is bad " << gamfile << endl;
            exit(9);
        }
        gamstream.close();
    }


    std::function<long double(int64_t)> fac = [](int64_t t){
        long double result = 1.0;
        for (int64_t i = 1; i <= t; ++i){
            result *= (long double) i; 
        }
        return result;
    };


    std::function<long double(int64_t)> big_fac = [](int64_t t){
        long double result = 1.0;
        int64_t to_calc = t;
        if (t > 15){
            // Screw calculating anything and return a big number.
            to_calc = 15;
        }
        for (int64_t i = 1; i <= to_calc; ++i){
            result *= (long double) i; 
        }
        return result;
    };



    // Creates a three-element vector, ref, het, alt
    std::function<long double(int64_t, int64_t, double)> get_binoms = [fac, big_fac](int64_t ref_count, int64_t alt_count, long double geno_prior){
        long double top_term = ( big_fac(alt_count + ref_count) );
        long double bottom_term = (big_fac(ref_count) * big_fac(alt_count));

        long double dat_binom = top_term / bottom_term;
        return dat_binom * pow(geno_prior, alt_count) * pow((1.0 - geno_prior), ref_count);
    };

    std::function<pair<long double, int>(int64_t, int64_t, double)> do_math = [get_binoms, fac](int64_t ref_count, int64_t alt_count, double allele_prior = 0.333){

        // Genotype priors: HOMOZYOUS REF : HET : HOMOZYGOUS ALT
        vector<double> geno_priors {0.1, 0.4, 0.8};
        vector<long double> data_probs (geno_priors.size());
        vector<long double> geno_probs (geno_priors.size());

#ifdef DEBUG
//cerr << "Ref count " << ref_count << endl;
        //           cerr << "Alt count " << alt_count << endl;
#endif
        // calculate genotype data probs
        for (int i = 0; i < geno_priors.size(); i++){
            long double x = get_binoms(ref_count, alt_count, geno_priors[i]) * allele_prior; 
#ifdef DEBUG
//                cerr << "Prob: " << x << endl;
#endif
            data_probs[i] = x;
        }
        long double sum_data_probs = std::accumulate(data_probs.begin(), data_probs.end(), 0.0);


        long double big_prob = 0.0;
        // 0 = ref, 1 = het, 2 = alt
        int geno_index = -1;


        vector<long double> prob_cache;
        for (int i = 0; i < geno_priors.size(); ++i){
            double tmp_binom = data_probs[i];
            double prob = tmp_binom / sum_data_probs;
            prob_cache.push_back(prob);
            if (prob > big_prob ){
                big_prob = prob;
                geno_index = i;
            }
            else if (prob == big_prob){
                cerr << "EQUAL PROBABILITIES OF GENOTYPES." << endl;
                geno_index = -1;
                break;
            }
        }

#ifdef DEBUG
        //          cerr << "Big prob: " << big_prob << endl
        //             <<  "geno_index " << geno_index << endl;
#endif

        return std::make_pair(big_prob, geno_index);
    };

    string sampleName = "Sample";
    for (auto it : hash_to_var){
        //cerr << it.second.position << " ";
        it.second.setVariantCallFile(outvcf);
        it.second.format.push_back("GT");
        auto& genotype_vector = it.second.samples[sampleName]["GT"];
        vector<int64_t> read_counts(it.second.alt.size() + 1, 0);
        for (int i = 0; i <= it.second.alt.size(); ++i){
            int64_t readsum = 0;
            string alt_id = "_alt_" + it.first + "_" + std::to_string(i);
            if (!use_snarls){
                readsum = allele_name_to_alignment_name[ alt_id ].size();
            }
            else{
                readsum = traversal_name_to_alignment_names[ alt_id ].size();
            }
            read_counts[i] = readsum;
            it.second.info["AD"].push_back(std::to_string(readsum));
        }

        pair<double, int> prob_and_geno_index = do_math(  read_counts[0], read_counts[1], 0.333);
        if (prob_and_geno_index.second == 0){
            genotype_vector.push_back("0/0");
        }
        else if (prob_and_geno_index.second == 1){
            genotype_vector.push_back("0/1");
        }
        else if(prob_and_geno_index.second == 2){
            genotype_vector.push_back("1/1");
        }
        else{
            genotype_vector.push_back("./.");
        }

        it.second.info["GP"].push_back(std::to_string(prob_and_geno_index.first));


        cout << it.second << endl;

    }

}

// void genotype(void variant_recall(VG* graph,
//         vcflib::VariantCallFile* vars,
//         FastaReference* ref_genome,
//         vector<FastaReference*> insertions,
//         string gamfile, bool isIndex){

//         }


}
