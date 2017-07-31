#include <cstdint>
#include "genotyper.hpp"


namespace vg {

    using namespace std;

    // genotyper:
    // use graph and reads to:
    // - Augment graph
    // - Find superbubbles or cactus branches to determine sites
    // - Generate proposals for paths through each site (from reads?)
    // - Compute affinities of each read for each proposed path through a site
    // - Compute diploid genotypes for each site
    // - Output as vcf or as native format

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

    void Genotyper::genotype_svs(VG* graph,
                        string gamfile,
                        string refpath){
            // Open up our GAM file
            ifstream gamstream;
            gamstream.open(gamfile);
            if (!gamstream.good()){
                cerr << "GAM file is no good" << endl;
                exit(2);
            }
            SRPE srrp;

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


                if (srrp.ff.mark_sv_alignments(a,b)){
                    sv_reads.push_back(make_pair(a, b));
                }
                                
                else if (srrp.ff.mark_smallVariant_alignments(a, b)){
                    direct_ins.push_back(a.path());
                    direct_ins.push_back(b.path());
                }

            };
            stream::for_each_interleaved_pair_parallel(gamstream, readfunc);
            vector<Translation> transls;
            if (refpath != ""){
                transls = graph->edit(direct_ins); // TODO could maybe use edit_fast??
               
                Deconstructor decon;
                decon.deconstruct(refpath, graph);
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
    void Genotyper::variant_recall(VG* graph,
            vcflib::VariantCallFile* vars,
            FastaReference* ref_genome,
            vector<FastaReference*> insertions,
            string gamfile, bool isIndex){
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
        unordered_map<string, list<Mapping> > gpaths( (graph->paths)._paths.begin(), (graph->paths)._paths.end()  );

        // To allow non-flat alleles, we want to use SnarlTraversals rather than
        // paths.
        unordered_map<string, set<int64_t> > snarl_name_to_node_set;
        map<string, set<string> > traversal_name_to_alignment_names;
        unordered_map<string, SnarlTraversal> name_to_traversal;
        map<string, Snarl> name_to_snarl;
        map<string, vector<SnarlTraversal> > snarl_name_to_traversals;

        if (use_snarls){
            SnarlFinder* snarl_finder = new CactusUltrabubbleFinder(*graph, "", true);
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
                    for (int i = 0; i < x.visits_size(); i++){
                        snarl_name_to_node_set[snarl->name()].insert( x.visits(i).node_id());
                    }    
                }
            }
            cerr << "Snarls processed, and we have this many: " << snarl_name_to_node_set.size() << endl;
        }
        

        // For each variant in VCF:
        vcflib::Variant var;
        while(vars->getNextVariant(var)){
            // Adjust the position offset, canonicalize any structural variants,
            // and get the sha1 hash of the variant and store that in a map for later.
            var.position -= 1;
            var.canonicalize_sv(*ref_genome, insertions, -1);
            string var_id = make_variant_id(var);

            hash_to_var[ var_id ] = var;

            if (!isIndex){
                // If we're just using the GAM, build a map Node -> variant and
                // a map alt_path_id (i.e. "allele name") -> NodeID. We'll use these to count
                // mappings to nodes later.
                for (int alt_ind = 0; alt_ind <= var.alt.size(); alt_ind++){
                    string alt_id = "_alt_" + var_id + "_" + std::to_string(alt_ind);
                    list<Mapping> x_path = gpaths[ alt_id ];
                    for (Mapping x_m : x_path){
                        allele_name_to_node_id[ alt_id ].insert(x_m.position().node_id());
                        node_to_variant[x_m.position().node_id()] = alt_id;
                        variant_nodes.insert(x_m.position().node_id());
                    }
                }

            }
            else{
                // We're using a GAM index, so we have to just get all alignments mapping
                // to a set of nodes.
                for (int alt_ind = 0; alt_ind <= var.alt.size(); alt_ind++){
                    string alt_id = "_alt_" + var_id + "_" + std::to_string(alt_ind);
                    list<Mapping> x_path = gpaths[ alt_id ];
                    for (Mapping x_m : x_path){
                        variant_nodes.insert(x_m.position().node_id());
                    }
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

        std::function<bool(const Mapping& m)> perfect_matches = [&](const Mapping& m){

            for (int i = 0; i < m.edit_size(); i++){
                Edit e = m.edit(i);
                if (e.to_length() != e.from_length() || !e.sequence().empty()){
                    return false;
                }
                return true;
            }
        };


        std::function<void(Alignment& a)> incr = [&](const Alignment& a){
            for (int i = 0; i < a.path().mapping_size(); i++){
                int64_t node_id = a.path().mapping(i).position().node_id();
                if (variant_nodes.count(node_id) && a.mapping_quality() > 20 && sufficient_matches(a.path().mapping(i))){
#pragma omp critical
                    {
                        // Get our variant allele hash
                        string all_str = node_to_variant[node_id];
                        // Add alignment to the variant's set of alignments
                        allele_name_to_alignment_name[all_str].insert(a.name());
                    }
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
        if (!isIndex){
            ifstream gamstream(gamfile);
            if (gamstream.good()){
                stream::for_each(gamstream, incr);
            }
            gamstream.close();
        }
        else if (use_snarls && !isIndex){
            ifstream gamstream(gamfile);
            if (gamstream.good()){
                stream::for_each(gamstream, count_traversal_supports);
            }
            gamstream.close();
        }
        else{
            Index gamindex;
            gamindex.open_read_only(gamfile);
            //gamindex.for_alignment_to_nodes(variant_nodes, incr);
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

    // void Genotyper::genotype(void Genotyper::variant_recall(VG* graph,
    //         vcflib::VariantCallFile* vars,
    //         FastaReference* ref_genome,
    //         vector<FastaReference*> insertions,
    //         string gamfile, bool isIndex){

    //         }


    void Genotyper::run(VG& graph,
            vector<Alignment>& alignments,
            ostream& out,
            string ref_path_name,
            string contig_name,
            string sample_name,
            string augmented_file_name,
            bool use_cactus,
            bool subset_graph,
            bool show_progress,
            bool output_vcf,
            bool output_json,
            int length_override,
            int variant_offset) {

        if(ref_path_name.empty()) {
            // Guess the ref path name
            if(graph.paths.size() == 1) {
                // Autodetect the reference path name as the name of the only path
                ref_path_name = (*graph.paths._paths.begin()).first;
            } else {
                ref_path_name = "ref";
            }
        }

        if(output_vcf && show_progress) {
#pragma omp critical (cerr)
            cerr << "Calling against path " << ref_path_name << endl;
        }

        if(sample_name.empty()) {
            // Set a default sample name
            sample_name = "SAMPLE";
        }

        // Make sure they have unique names.
        set<string> names_seen;
        // We warn about duplicate names, but only once.
        bool duplicate_names_warned = false;
        for(size_t i = 0; i < alignments.size(); i++) {
            if(alignments[i].name().empty()) {
                // Generate a name
                alignments[i].set_name("_unnamed_alignment_" + to_string(i));
            }
            if(names_seen.count(alignments[i].name())) {
                // This name is duplicated
                if(!duplicate_names_warned) {
                    // Warn, but only once
                    cerr << "Warning: duplicate alignment names present! Example: " << alignments[i].name() << endl;
                    duplicate_names_warned = true;
                }

                // Generate a new name
                // TODO: we assume this is unique
                alignments[i].set_name("_renamed_alignment_" + to_string(i));
                assert(!names_seen.count(alignments[i].name()));
            }
            names_seen.insert(alignments[i].name());
        }
        names_seen.clear();

        // Suck out paths
        vector<Path> paths;
        for(auto& alignment : alignments) {
            // Copy over each path, naming it after its alignment
            // and trimming so that it begins and ends with a match to avoid
            // creating a bunch of stubs.
            Path path = trim_hanging_ends(alignment.path());
            path.set_name(alignment.name());
            paths.push_back(path);
        }

        // Run them through vg::edit() to add them to the graph. Save the translations.
        vector<Translation> augmentation_translations = graph.edit(paths);
        translator.load(augmentation_translations);

        if(show_progress) {
#pragma omp critical (cerr)
            cerr << "Augmented graph; got " << augmentation_translations.size() << " translations" << endl;
        }

        // Make sure that we actually have an index for traversing along paths.
        graph.paths.rebuild_mapping_aux();

        if(!augmented_file_name.empty()) {
            ofstream augmented_stream(augmented_file_name);
            graph.serialize_to_ostream(augmented_stream);
            augmented_stream.close();
        }

        // store the reads that are embedded in the augmented graph, by their unique names
        map<string, Alignment*> reads_by_name;
        for(auto& alignment : alignments) {
            reads_by_name[alignment.name()] = &alignment;
            // Make sure to replace the alignment's path with the path it has in the augmented graph
            list<Mapping>& mappings = graph.paths.get_path(alignment.name());
            alignment.mutable_path()->clear_mapping();
            for(auto& mapping : mappings) {
                // Copy over all the transformed mappings
                *alignment.mutable_path()->add_mapping() = mapping;
            }
        }
#pragma omp critical (cerr)
        cerr << "Converted " << alignments.size() << " alignments to embedded paths" << endl;


        // We need to decide if we want to work on the full graph or just on the subgraph that has any support.

        // Set up our genotypekit members.
        SnarlFinder* snarl_finder;

        // Find all the sites in either the main graph or the subset
        vector<Genotyper::Site> sites;

        if(subset_graph) {
            // We'll collect the supported subset of the original graph
            set<Node*> supported_nodes;
            set<Edge*> supported_edges;

            for(auto& name_and_read : reads_by_name) {
                // Go through all the paths followed by reads
                auto& path = name_and_read.second->path();
                for(size_t i = 0; i < path.mapping_size(); i++) {
                    // Look at all the nodes we visit along this read
                    id_t node_id = path.mapping(i).position().node_id();
                    // Make sure they are all supported
                    supported_nodes.insert(graph.get_node(node_id));

                    if(i > 0) {
                        // We also need the edge from the last mapping to this one.
                        // Make the two sides we connected moving from the last mapping to this one.
                        NodeSide last(path.mapping(i - 1).position().node_id(), !path.mapping(i - 1).position().is_reverse());
                        NodeSide here(node_id, path.mapping(i).position().is_reverse());

                        Edge* edge = graph.get_edge(last, here);

                        if(edge == nullptr) {
                            cerr << "Error! Edge " << last << " to " << here
                                << " from path " << name_and_read.first << " is missing!" << endl;
                            exit(1);
                        }

                        // We know the graph will have the edge
                        supported_edges.insert(edge);
                    }

                }
            }

            // We also want to support all nodes and edges used by the reference path.
            // TODO: once Cactus can root without hints, we can discard this
            if(graph.paths.has_path(ref_path_name)) {
                // We actually have a reference path, so get it for traversing.
                list<Mapping>& ref_mappings = graph.paths.get_path(ref_path_name);
                // We need to remember the previous mapping for finding edges
                list<Mapping>::iterator last_mapping = ref_mappings.end();
                for(list<Mapping>::iterator mapping = ref_mappings.begin(); mapping != ref_mappings.end(); ++mapping) {
                    // For each mapping along the reference path

                    // What node is it on?
                    id_t node_id = mapping->position().node_id();
                    // Make sure it is supported
                    supported_nodes.insert(graph.get_node(node_id));

                    if(last_mapping != ref_mappings.end()) {
                        // We're coming from another mapping and need to support the edge

                        NodeSide last(last_mapping->position().node_id(), !last_mapping->position().is_reverse());
                        NodeSide here(node_id, mapping->position().is_reverse());

                        Edge* edge = graph.get_edge(last, here);

                        if(edge == nullptr) {
                            cerr << "Error! Edge " << last << " to " << here
                                << " from path " << ref_path_name << " is missing!" << endl;
                            exit(1);
                        }

                        // We know the graph will have the edge
                        supported_edges.insert(edge);

                    }

                    // Save the iterator so we can get the next edge
                    last_mapping = mapping;

                }
            }

            // Make the subset graph of only supported nodes and edges (which will
            // internally contain copies of all of them).
            VG subset;
            subset.add_nodes(supported_nodes);
            subset.add_edges(supported_edges);

            if(graph.paths.has_path(ref_path_name)) {
                // Copy over the reference path
                subset.paths.extend(graph.paths.path(ref_path_name));
            }


            if(show_progress) {
#pragma omp critical (cerr)
                cerr << "Looking at subset of " << subset.size() << " nodes" << endl;
            }

            // Unfold/unroll, find the superbubbles, and translate back. Note that
            // we can only use Cactus with the ref path if it survived the
            // subsetting.

            sites = use_cactus ? (
                    subset.paths.has_path(ref_path_name) ?
                    find_sites_with_cactus(subset, ref_path_name)
                    : find_sites_with_cactus(subset)
                    )
                : find_sites_with_supbub(subset);

            for(auto& site : sites) {
                // Translate all the NodeTraversals back to node pointers in the
                // non-subset graph
                site.start.node = graph.get_node(site.start.node->id());
                site.end.node = graph.get_node(site.end.node->id());
            }

        } else {
            // Don't need to mess around with creating a subset.

            if(show_progress) {
#pragma omp critical (cerr)
                cerr << "Looking at graph of " << graph.size() << " nodes" << endl;
            }

            // Unfold/unroll, find the superbubbles, and translate back.
            graph.sort();
            sites = use_cactus ? find_sites_with_cactus(graph, ref_path_name)
            : find_sites_with_supbub(graph);



//             bool filter_trivials = true;
//             snarl_finder =  new CactusUltrabubbleFinder(graph, "", filter_trivials);

//             // Get our snarls.
//             SnarlManager snarl_manager = snarl_finder->find_snarls();
//             vector<const Snarl*> snarl_roots = snarl_manager.top_level_snarls();

//             // Enumerate traversals through our snarls
//             // TODO : handle reversing sites in sites
            // TraversalFinder* trav_finder = new ExhaustiveTraversalFinder(graph, snarl_manager);
            // for (auto x : snarl_roots){
            //     vector<SnarlTraversal> site_traversals = trav_finder->find_traversals(*x);
            //     Site current;
            //     current.start.node = graph.get_node(x->start().node_id());
            //     current.end.node = graph.get_node(x->end().node_id());
            //     set<id_t> contents;
            //     for (auto tr : site_traversals){
            //         for (auto v : tr.visits()){
            //             contents.insert(v.node_id());
            //         }
            //     }
            //     current.contents = contents;
            //     sites.push_back(current);
            // }


         }

        if(show_progress) {
#pragma omp critical (cerr)
            {
                cerr << "Found " << sites.size() << " superbubbles" << endl;
                // for (auto x : sites){
                //     cout << x.start.node->id() << " ";
                //     for (auto y : x.contents){
                //         cout << y << " ";
                //     }
                //     cout << x.end.node->id() << endl;
                // }

                // exit(1);
            }

    }

        // We're going to count up all the affinities we compute
        size_t total_affinities = 0;

        // We need a buffer for output
        vector<vector<Locus>> buffer;
        int thread_count = get_thread_count();
        buffer.resize(thread_count);

        // If we're doing VCF output we need a VCF header
        vcflib::VariantCallFile* vcf = nullptr;
        // And a reference index tracking the primary path
        PathIndex* reference_index = nullptr;
        if(output_vcf) {
            // Build a reference index on our reference path
            // Make sure to capture the sequence
            reference_index = new PathIndex(graph, ref_path_name, true);

            // Start up a VCF
            vcf = start_vcf(cout, *reference_index, sample_name, contig_name, length_override);
        }

        // We want to do this in parallel, but we can't #pragma omp parallel for over a std::map
#pragma omp parallel shared(total_affinities)
        {
#pragma omp single nowait
            {
                for(auto it = sites.begin(); it != sites.end(); it++) {
                    // For each site in parallel

#pragma omp task firstprivate(it) shared(total_affinities)
                    {

                        auto& site = *it;

                        // Report the site to our statistics code
                        report_site(site, reference_index);

                        int tid = omp_get_thread_num();

                        // Get all the paths through the site supported by enough reads, or by real named paths
                        vector<list<NodeTraversal>> paths = get_paths_through_site(graph, site, reads_by_name);

                        if(paths.size() == 0) {
                            // TODO: this compensates for inside-out sites from
                            // Cactus. Make Cactus give us sites that can actually
                            // be traversed through.

                            // Flip the site around and try again
                            swap(site.start, site.end);
                            vector<list<NodeTraversal>> reverse_paths = get_paths_through_site(graph, site, reads_by_name);
                            if(reverse_paths.size() != 0) {
                                // We actually got some paths. Use them
                                swap(paths, reverse_paths);
#pragma omp critical (cerr)
                                cerr << "Warning! Corrected inside-out site " << site.end << " - " << site.start << endl;
                            } else {
                                // Put original start and end back for complaining
                                swap(site.start, site.end);
                            }
                        }

                        if(reference_index != nullptr &&
                                reference_index->by_id.count(site.start.node->id()) && 
                                reference_index->by_id.count(site.end.node->id())) {
                            // This site is on the reference (and we are indexing a reference because we are going to vcf)

                            // Where do the start and end nodes fall in the reference?
                            auto start_ref_appearance = reference_index->by_id.at(site.start.node->id());
                            auto end_ref_appearance = reference_index->by_id.at(site.end.node->id());

                            // Are the ends running with the reference (false) or against it (true)
                            auto start_rel_orientation = (site.start.backward != start_ref_appearance.second);
                            auto end_rel_orientation = (site.end.backward != end_ref_appearance.second);

                            if(show_progress) {
                                // Determine where the site starts and ends along the reference path
#pragma omp critical (cerr)
                                cerr << "Site " << site.start << " - " << site.end << " runs reference " <<
                                    start_ref_appearance.first << " to " <<
                                    end_ref_appearance.first << endl;

                                if(!start_rel_orientation && !end_rel_orientation &&
                                        end_ref_appearance.first < start_ref_appearance.first) {
                                    // The site runs backward in the reference (but somewhat sensibly).
#pragma omp critical (cerr)
                                    cerr << "Warning! Site runs backwards!" << endl;
                                }

                            }

                        }

                        // Even if it looks like there's only one path, it might not
                        // be the reference path, because the reference path might
                        // not have passed the min recurrence filter. So we can't
                        // skip things yet.
                        if(paths.empty()) {
                            // Don't do anything for superbubbles with no routes through
                            if(show_progress) {
#pragma omp critical (cerr)
                                cerr << "Site " << site.start << " - " << site.end << " has " << paths.size() <<
                                    " alleles: skipped for having no alleles" << endl;
                            }
                        } else {

                            if(show_progress) {
#pragma omp critical (cerr)
                                cerr << "Site " << site.start << " - " << site.end << " has " << paths.size() << " alleles" << endl;
                                for(auto& path : paths) {
                                    // Announce each allele in turn
#pragma omp critical (cerr)
                                    cerr << "\t" << traversals_to_string(path) << endl;
                                }
                            }

                            // Compute the lengths of all the alleles
                            set<size_t> allele_lengths;
                            for(auto& path : paths) {
                                allele_lengths.insert(traversals_to_string(path).size());
                            }

                            // Get the affinities for all the paths
                            map<Alignment*, vector<Genotyper::Affinity>> affinities;

                            if(allele_lengths.size() > 1 && realign_indels) {
                                // This is an indel, because we can change lengths. Use the slow route to do idnel realignment.
                                affinities = get_affinities(graph, reads_by_name, site, paths);
                            } else {
                                // Just use string comparison. Don't re-align when
                                // length can't change, or when indle realignment is
                                // off.
                                affinities = get_affinities_fast(graph, reads_by_name, site, paths);
                            }

                            if(show_progress) {
                                // Sum up all the affinity counts by consistency flags
                                map<string, size_t> consistency_combo_counts;

                                // And average raw scores by alleles they are
                                // consistent with, for things consistent with just
                                // one allele.
                                vector<double> score_totals(paths.size());
                                vector<size_t> score_counts(paths.size());

                                for(auto& alignment_and_affinities : affinities) {
                                    // For every alignment, make a string describing which alleles it is consistent with.
                                    string consistency;

                                    // How many alleles are we consstent with?
                                    size_t consistent_allele_count = 0;
                                    // And which one is it, if it's only one?
                                    int chosen = -1;

                                    for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
                                        auto& affinity = alignment_and_affinities.second.at(i);
                                        if(affinity.consistent) {
                                            // Consistent alleles get marked with a 1
                                            consistency.push_back('1');

                                            // Say we're consistent with an allele
                                            chosen = i;
                                            consistent_allele_count++;

                                        } else {
                                            // Inconsistent ones get marked with a 0
                                            consistency.push_back('0');
                                        }
                                    }

                                    if(consistent_allele_count == 1) {
                                        // Add in the non-normalized score for the average
                                        score_totals.at(chosen) += alignment_and_affinities.second.at(chosen).score;
                                        score_counts.at(chosen)++;
                                    }

#ifdef debug
#pragma omp critical (cerr)
                                    cerr << consistency << ": " << alignment_and_affinities.first->sequence() << endl;
#endif


                                    // Increment the count for that pattern
                                    consistency_combo_counts[consistency]++;
                                }

#pragma omp critical (cerr)
                                {
                                    cerr << "Support patterns:" << endl;
                                    for(auto& combo_and_count : consistency_combo_counts) {
                                        // Spit out all the counts for all the combos
                                        cerr << "\t" << combo_and_count.first << ": " << combo_and_count.second << endl;
                                    }

                                    cerr << "Average scores for unique support:" << endl;
                                    for(size_t i = 0; i < score_totals.size(); i++) {
                                        // Spit out average scores of uniquely supporting reads for each allele that has them.
                                        if(score_counts.at(i) > 0) {
                                            cerr << "\t" << traversals_to_string(paths.at(i)) << ": "
                                                << score_totals.at(i) / score_counts.at(i) << endl;
                                        } else {
                                            cerr << "\t" << traversals_to_string(paths.at(i)) << ": --" << endl;
                                        }
                                    }

                                }
                            }

                            for(auto& alignment_and_affinities : affinities) {
#pragma omp critical (total_affinities)
                                total_affinities += alignment_and_affinities.second.size();
                            }

                            // Get a genotyped locus in the original frame
                            Locus genotyped = genotype_site(graph, site, paths, affinities);
                            if (output_vcf) {
                                // Get 0 or more variants from the superbubble
                                vector<vcflib::Variant> variants =
                                    locus_to_variant(graph, site, *reference_index, *vcf, genotyped, sample_name);
                                for(auto& variant : variants) {
                                    // Fix up all the variants
                                    if(!contig_name.empty()) {
                                        // Override path name
                                        variant.sequenceName = contig_name;
                                    } else {
                                        // Keep path name
                                        variant.sequenceName = ref_path_name;
                                    }
                                    variant.position += variant_offset;

#pragma omp critical(cout)
                                    cout << variant << endl;
                                }
                            } else {
                                // project into original graph
                                genotyped = translator.translate(genotyped);
                                // record a consistent name based on the start and end position of the first allele
                                stringstream name;
                                if (genotyped.allele_size() && genotyped.allele(0).mapping_size()) {
                                    name << make_pos_t(genotyped.allele(0).mapping(0).position())
                                        << "_"
                                        << make_pos_t(genotyped
                                                .allele(0)
                                                .mapping(genotyped.allele(0).mapping_size()-1)
                                                .position());
                                }
                                genotyped.set_name(name.str());
                                if (output_json) {
                                    // Dump in JSON
#pragma omp critical (cout)
                                    cout << pb2json(genotyped) << endl;
                                } else {
                                    // Write out in Protobuf
                                    buffer[tid].push_back(genotyped);
                                    stream::write_buffered(cout, buffer[tid], 100);
                                }
                            }
                        }
                    }
                }
            }
        }           

        if(!output_json && !output_vcf) {
            // Flush the protobuf output buffers
            for(int i = 0; i < buffer.size(); i++) {
                stream::write_buffered(cout, buffer[i], 0);
            }
        } 


        if(show_progress) {
#pragma omp critical (cerr)
            cerr << "Computed " << total_affinities << " affinities" << endl;
        }

        // Dump statistics before the sites go away, so the pointers won't be dangling
        print_statistics(cerr);

        if(output_vcf) {
            delete vcf;
            delete reference_index;
        }

    }


    pair<pair<int64_t, int64_t>, bool> Genotyper::get_site_reference_bounds(const Site& site, const PathIndex& index) {
        // Grab the start and end node IDs.
        auto first_id = site.start.node->id();
        auto last_id = site.end.node->id();

        if(!index.by_id.count(first_id) || !index.by_id.count(last_id)) {
            // Site isn;t actually on the reference path so return a sentinel.
            return make_pair(make_pair(-1, -1), false);
        }

        // The position we have stored for this start node is the first
        // position along the reference at which it occurs. Our bubble
        // goes forward in the reference, so we must come out of the
        // opposite end of the node from the one we have stored.
        auto referenceIntervalStart = index.by_id.at(first_id).first + site.start.node->sequence().size();

        // The position we have stored for the end node is the first
        // position it occurs in the reference, and we know we go into
        // it in a reference-concordant direction, so we must have our
        // past-the-end position right there.
        auto referenceIntervalPastEnd = index.by_id.at(last_id).first;

        // Is this bubble articulated backwards relative to the reference?
        bool site_is_reverse = false;

        if(referenceIntervalStart > referenceIntervalPastEnd) {
            // Everything we know about the site is backwards relative to the reference. Flip it around frontways.
            site_is_reverse = true;
            swap(first_id, last_id);
            // Recalculate reference positions Use the end node, which we've now
            // made first_id, to get the length offset to the start of the actual
            // internal variable bit.
            referenceIntervalStart = index.by_id.at(first_id).first + site.end.node->sequence().size();
            referenceIntervalPastEnd = index.by_id.at(last_id).first;
        }

        return make_pair(make_pair(referenceIntervalStart, referenceIntervalPastEnd), site_is_reverse);
    }

    /**
     * Turn the given path (which must be a thread) into an allele. Drops the first
     * and last mappings and looks up the sequences for the nodes of the others.
     */
    string allele_to_string(VG& graph, const Path& allele) {
        stringstream stream;

        for(size_t i = 1; i < allele.mapping_size() - 1; i++) {
            // Get the sequence for each node
            string node_string = graph.get_node(allele.mapping(i).position().node_id())->sequence();

            if(allele.mapping(i).position().is_reverse()) {
                // Flip it
                node_string = reverse_complement(node_string);
            }
            // Add it to the stream
            stream << node_string;
        }

        return stream.str();
    }

    int Genotyper::alignment_qual_score(VG& graph, const Site& site, const Alignment& alignment) {
        if(alignment.quality().empty()) {
            // Special case: qualities not given. Assume something vaguely sane so
            // we can genotype without quality.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "No base qualities. Assuming default quality of " << default_sequence_quality << endl;
#endif
            return default_sequence_quality;
        }

        // Go through all the qualities in the site
        // TODO: can we do this in place?
        string relevant_qualities = get_qualities_in_site(graph, site, alignment);

        if(relevant_qualities.empty()) {
            // No qualities available internal to the site for this read. Must be a
            // pure-deletion allele.
#ifdef debug
#pragma omp critical (cerr)
            cerr << "No internal qualities. Assuming default quality of " << default_sequence_quality << endl;
#endif
            // TODO: look at bases on either side of the deletion.
            return default_sequence_quality;
        }

        double total = 0;
        for(auto& quality : relevant_qualities) {
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Quality: " << (int)quality << endl;
#endif
            total += quality;
        }
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Count: " << relevant_qualities.size() << endl;
#endif
        // Make the total now actually be an average
        total /= relevant_qualities.size();
        return round(total);
    }

    vector<Genotyper::Site> Genotyper::find_sites_with_supbub(VG& graph) {

        // Set up our output vector
        vector<Site> to_return;

        // Unfold the graph
        // Copy the graph and unfold the copy. We need to hold this translation
        // table from new node ID to old node and relative orientation.
        unordered_map<vg::id_t, pair<vg::id_t, bool>> unfold_translation;
        auto transformed = graph.unfold(unfold_max_length, unfold_translation);

        // Fix up any doubly reversed edges
        transformed.flip_doubly_reversed_edges();

        // Now dagify the graph. We need to hold this translation table from new
        // node ID to old node and relative orientation.
        unordered_map<vg::id_t, pair<vg::id_t, bool>> dag_translation;
        transformed = transformed.dagify(dagify_steps, dag_translation);

        // Compose the complete translation
        unordered_map<vg::id_t, pair<vg::id_t, bool>> overall_translation = transformed.overlay_node_translations(dag_translation, unfold_translation);
        dag_translation.clear();
        unfold_translation.clear();

        // Find the superbubbles in the DAG
        map<pair<id_t, id_t>, vector<id_t>> superbubbles = vg::superbubbles(transformed);

        for(auto& superbubble : superbubbles) {

            // Translate the superbubble coordinates into NodeTraversals.
            // This is coming from a DAG so we only need to look at the translation for orientation.
            auto& start_translation = overall_translation[superbubble.first.first];
            NodeTraversal start(graph.get_node(start_translation.first), start_translation.second);

            auto& end_translation = overall_translation[superbubble.first.second];
            NodeTraversal end(graph.get_node(end_translation.first), end_translation.second);

            // Make a Site and tell it where to start and end
            Site site;
            site.start = start;
            site.end = end;

            for(auto id : superbubble.second) {
                // Translate each ID and put it in the set
                Node* node = graph.get_node(overall_translation[id].first);
            }

            // Save the site
            to_return.emplace_back(std::move(site));
        }

        // Give back the collection of sites
        return to_return;    

    }

    vector<Genotyper::Site> Genotyper::find_sites_with_cactus(VG& graph, const string& ref_path_name) {

        // Set up our output vector
        vector<Site> to_return;

        // cactus needs the nodes to be sorted in order to find a source and sink
        graph.sort();

        // todo: use deomposition instead of converting tree into flat structure
        BubbleTree* bubble_tree = ultrabubble_tree(graph);

        bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {
                Bubble& bubble = node->v;
                // cut root to be consistent with superbubbles()
                if (node != bubble_tree->root) {
                set<id_t> nodes{bubble.contents.begin(), bubble.contents.end()};
                NodeTraversal start(graph.get_node(bubble.start.node), !bubble.start.is_end);
                NodeTraversal end(graph.get_node(bubble.end.node), bubble.end.is_end);
                // Fill in a Site. Make sure to preserve original endpoint
                // ordering, because swapping them without flipping their
                // orientation flags will make an inside-out site.
                Site site;
                site.start = start;
                site.end = end;
                swap(site.contents, nodes);
                // Save the site
                to_return.emplace_back(std::move(site));
                }
                });

        delete bubble_tree;

        return to_return;
    }

    vector<list<NodeTraversal>> Genotyper::get_paths_through_site(VG& graph, const Site& site,
            const map<string, Alignment*>& reads_by_name) {
        // We're going to emit traversals supported by any paths in the graph.

        // Put all our subpaths in here to deduplicate them by sequence they spell
        // out. And to count occurrences. Note that the occurrence count will be
        // boosted to min_recurrence if a non-read path in the graph supports a
        // certain traversal string, so we don't end up dropping unsupported ref
        // alleles.
        map<string, pair<list<NodeTraversal>, int>> results;

#ifdef debug
#pragma omp critical (cerr)
        cerr << "Looking for paths between " << site.start << " and " << site.end << endl;
#endif

        if(graph.paths.has_node_mapping(site.start.node) && graph.paths.has_node_mapping(site.end.node)) {
            // If we have some paths that visit both ends (in some orientation)

            // Get all the mappings to the end node, by path name
            auto& endmappings_by_name = graph.paths.get_node_mapping(site.end.node);

            for(auto& name_and_mappings : graph.paths.get_node_mapping(site.start.node)) {
                // Go through the paths that visit the start node

                // Grab their names
                auto& name = name_and_mappings.first;

                if(!endmappings_by_name.count(name_and_mappings.first)) {
                    //cerr << "no endmappings match" << endl;
                    // No path by this name has any mappings to the end node. Skip
                    // it early.
                    continue;
                }

                for(auto* mapping : name_and_mappings.second) {
                    // Start at each mapping in the appropriate orientation

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "Trying mapping of read/path " << name_and_mappings.first << endl;
#endif

                    // How many times have we gone to the next mapping looking for a
                    // mapping to the end node in the right orientation?
                    size_t traversal_count = 0;

                    // Do we want to go left (true) or right (false) from this
                    // mapping? If start is a forward traversal and we found a
                    // forward mapping, we go right. If either is backward we go
                    // left, and if both are backward we go right again.
                    bool traversal_direction = mapping->position().is_reverse() != site.start.backward;

                    // What orientation would we want to find the end node in? If
                    // we're traveling backward, we expect to find it in the
                    // opposite direction to the one we were given.
                    bool expected_end_orientation = site.end.backward != traversal_direction;

                    // We're going to fill in this list with traversals.
                    list<NodeTraversal> path_traversed;

                    // And we're going to fill this with the sequence
                    stringstream allele_stream;

                    while(mapping != nullptr && traversal_count < max_path_search_steps) {
                        // Traverse along until we hit the end traversal or take too
                        // many steps

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\tTraversing " << pb2json(*mapping) << endl;
#endif

                        // Say we visit this node along the path, in this orientation
                        path_traversed.push_back(NodeTraversal(graph.get_node(mapping->position().node_id()),
                                    mapping->position().is_reverse() != traversal_direction));

                        // Stick the sequence of the node (appropriately oriented) in the stream for the allele sequence
                        string seq = graph.get_node(mapping->position().node_id())->sequence();
                        allele_stream << (path_traversed.back().backward ? reverse_complement(seq) : seq);

                        if(mapping->position().node_id() == site.end.node->id() && mapping->position().is_reverse() == expected_end_orientation) {
                            // We have stumbled upon the end node in the orientation we wanted it in.
                            if(results.count(allele_stream.str())) {
                                // It is already there! Increment the observation count.
#ifdef debug
#pragma omp critical (cerr)
                                cerr << "\tFinished; got known sequence " << allele_stream.str() << endl;
#endif

                                if(reads_by_name.count(name)) {
                                    // We are a read. Just increment count
                                    results[allele_stream.str()].second++;
                                } else {
                                    // We are a named path (like "ref")
                                    if(results[allele_stream.str()].second < min_recurrence) {
                                        // Ensure that this allele doesn't get
                                        // eliminated, since ref or some other named
                                        // path supports it.
                                        results[allele_stream.str()].second = min_recurrence;
                                    } else {
                                        results[allele_stream.str()].second++;
                                    }
                                }
                            } else {
                                // Add it in. Give it a count of 1 if we are a read,
                                // and a count of min_recurrence (so it doesn't get
                                // filtered later) if we are a named non-read path
                                // (like "ref").
                                results[allele_stream.str()] = make_pair(path_traversed,
                                        reads_by_name.count(name) ? 1 : min_recurrence);
#ifdef debug
#pragma omp critical (cerr)
                                cerr << "\tFinished; got novel sequence " << allele_stream.str() << endl;
#endif
                            }

                            if(reads_by_name.count(name)) {
                                // We want to log stats on reads that read all the
                                // way through sites. But since we may be called
                                // multiple times we need to send the unique read
                                // name too.
                                report_site_traversal(site, name);
                            }

                            // Then try the next embedded path
                            break;
                        }

                        // Otherwise just move to the right (or left)
                        if(traversal_direction) {
                            // We're going backwards
                            mapping = graph.paths.traverse_left(mapping);
                        } else {
                            // We're going forwards
                            mapping = graph.paths.traverse_right(mapping);
                        }
                        // Tick the counter so we don't go really far on long paths.
                        traversal_count++;

                    }


                }
            }

        }

        // Now collect the unique results
        vector<list<NodeTraversal>> to_return;

        for(auto& result : results) {
            // Break out each result
            const string& seq = result.first;
            auto& traversals = result.second.first;
            auto& count = result.second.second;

            if(count < min_recurrence) {
                // We don't have enough initial hits for this sequence to justify
                // trying to re-align the rest of the reads. Skip it. Note that the
                // reference path (and other named paths) will stuff in at least
                // min_recurrence to make sure we don't throw out their alleles.
                continue;
            }

            // Send out each list of traversals
            to_return.emplace_back(std::move(traversals));
        }

        return to_return;
    }

    template<typename T> inline void set_intersection(const unordered_set<T>& set_1, const unordered_set<T>& set_2,
            unordered_set<T>* out_intersection ) {
        bool set_1_smaller = set_1.size() < set_2.size();
        const unordered_set<T>& smaller_set = set_1_smaller ? set_1 : set_2;
        const unordered_set<T>& larger_set = set_1_smaller ? set_2 : set_1;

        *out_intersection = unordered_set<T>();
        unordered_set<T>& intersection = *out_intersection;
        for (T item : smaller_set) {
            if (larger_set.count(item)) {
                intersection.insert(item);
            }
        }
    }


    // TODO properly handle cycles inside superbubble by including multiplicity of an edge in a path
    void Genotyper::edge_allele_labels(const VG& graph,
            const Site& site,
            const vector<list<NodeTraversal>>& superbubble_paths,
            unordered_map<pair<NodeTraversal, NodeTraversal>,
            unordered_set<size_t>,
            hash_oriented_edge>* out_edge_allele_sets)
    {
        // edges are indicated by the pair of node ids they connect and what orientation
        // they start and end in (true indicates forward)
        *out_edge_allele_sets = unordered_map<pair<NodeTraversal, NodeTraversal>, unordered_set<size_t>, hash_oriented_edge>();
        unordered_map<pair<NodeTraversal, NodeTraversal>, unordered_set<size_t>, hash_oriented_edge>& edge_allele_sets = *out_edge_allele_sets;

        for (size_t i = 0; i < superbubble_paths.size(); i++) {
            list<NodeTraversal> path = superbubble_paths[i];
            // start at second node so we can look at edge leading into it
            auto iter = path.begin();
            iter++;
            // can stop before last node because only interested in case where allele is ambiguous
            auto last_node = path.end();
            last_node--;
            for (; iter != last_node; iter++) {

                auto prev_iter = iter;
                prev_iter--;

                NodeTraversal node = *iter;
                NodeTraversal prev_node = *prev_iter;

                pair<NodeTraversal, NodeTraversal> edge = make_pair(prev_node, node);

                if (!edge_allele_sets.count(edge)) {
                    edge_allele_sets.emplace(edge, unordered_set<size_t>());
                }

                // label the edge with this allele (indicated by its position in the vector)
                edge_allele_sets.at(edge).insert(i);
            }

        }
    }

    // find the log conditional probability of each ambiguous allele set given each true allele
    void Genotyper::allele_ambiguity_log_probs(const VG& graph,
            const Site& site,
            const vector<list<NodeTraversal>>& superbubble_paths,
            const unordered_map<pair<NodeTraversal, NodeTraversal>,
            unordered_set<size_t>,
            hash_oriented_edge>& edge_allele_sets,
            vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>* out_allele_ambiguity_probs)
    {
        *out_allele_ambiguity_probs = vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>();
        vector<unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>>& ambiguous_allele_probs = *out_allele_ambiguity_probs;
        ambiguous_allele_probs.reserve(superbubble_paths.size());

        for (size_t i = 0; i < superbubble_paths.size(); i++) {
            ambiguous_allele_probs[i] = unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>();
            unordered_map<vector<size_t>, double, hash_ambiguous_allele_set>& allele_probs = ambiguous_allele_probs[i];
            list<NodeTraversal> path = superbubble_paths[i];

            // consider both prefixes and suffixes that partially cross the superbubble
            for (bool forward : {true, false}) {
                // the set of alleles that this prefix of the current allele is consistent with
                unordered_set<size_t> prefix_allele_set;

                // find the first and last node to iterate through in this orientation
                // stop one node before last node because only interested in case where allele is ambiguous
                list<NodeTraversal>::iterator iter;
                list<NodeTraversal>::iterator final;
                if (forward) {
                    // extend prefix forward through the superbubble
                    iter = path.begin();
                    final = path.end();
                    final--;
                }
                else {
                    // extend suffix backward through the superbubble
                    iter = path.end();
                    iter--;
                    final = path.begin();
                }
                // iterate forwards or backwards along edges of path
                while (true) {
                    // get the two nodes of the edge in the order they were entered into the allele label map
                    NodeTraversal node;
                    NodeTraversal next_node;
                    auto next_iter = iter;
                    if (forward) {
                        next_iter++;
                        node = *iter;
                        next_node = *next_iter;
                    }
                    else {
                        next_iter--;
                        node = *next_iter;
                        next_node = *iter;
                    }
                    if (next_iter == final) {
                        break;
                    }

                    pair<NodeTraversal, NodeTraversal> edge = make_pair(node, next_node);
                    const unordered_set<size_t>& edge_allele_set = edge_allele_sets.at(edge);

                    if (prefix_allele_set.empty()) {
                        // first edge in path, consistent with all alleles edge is labeled with
                        prefix_allele_set = edge_allele_set;
                    }
                    else {
                        // take set intersection of prefix alleles and the edge's allele labels
                        unordered_set<size_t> new_prefix_allele_set;
                        set_intersection<size_t>(prefix_allele_set, edge_allele_set, &new_prefix_allele_set);
                        prefix_allele_set = new_prefix_allele_set;
                    }

                    // convert unordered set into a sorted vector for consistent hash-key behavior
                    vector<size_t> allele_set_key;
                    allele_set_key.reserve(prefix_allele_set.size());
                    for (size_t allele : prefix_allele_set) {
                        allele_set_key.emplace_back(allele);
                    }
                    sort(allele_set_key.begin(), allele_set_key.end());

                    // add the length of the sequence to the probability of that allele set (will normalize later)
                    if (allele_probs.count(allele_set_key)) {
                        allele_probs.at(allele_set_key) += node.node->sequence().length();
                    }
                    else {
                        allele_probs.at(allele_set_key) = node.node->sequence().length();
                    }

                    // iterate forward or backward through bubble
                    if (forward) {
                        iter++;
                    }
                    else {
                        iter--;
                    }
                }
            }

            // normalize lengths to probabilities (must sum to 1)
            size_t total_ambiguous_length = 0;
            for (auto& allele_length : allele_probs) {
                total_ambiguous_length += allele_length.second;
            }
            for (auto& allele_length : allele_probs) {
                allele_length.second = log(allele_length.second / total_ambiguous_length);
            }
        }
    }




    map<Alignment*, vector<Genotyper::Affinity>>
        Genotyper::get_affinities(VG& graph,
                const map<string, Alignment*>& reads_by_name,
                const Site& site,
                const vector<list<NodeTraversal>>& superbubble_paths) {

            // Grab our thread ID, which determines which aligner we get.
            int tid = omp_get_thread_num();

            // We're going to build this up gradually, appending to all the vectors.
            map<Alignment*, vector<Affinity>> to_return;

            // What reads are relevant to this superbubble?
            set<string> relevant_read_names;

#ifdef debug
#pragma omp critical (cerr)
            cerr << "Superbubble contains " << site.contents.size() << " nodes" << endl;
#endif

            for(auto id : site.contents) {
                // For every node in the superbubble, what paths visit it?
                if(graph.paths.has_node_mapping(id)) {
                    auto& mappings_by_name = graph.paths.get_node_mapping(id);
                    for(auto& name_and_mappings : mappings_by_name) {
                        // For each path visiting the node
                        if(reads_by_name.count(name_and_mappings.first)) {
                            // This path is a read, so add the name to the set if it's not
                            // there already
                            relevant_read_names.insert(name_and_mappings.first);
                        }    
                    }
                }
            }

            // What IDs are visited by these reads?
            unordered_set<id_t> relevant_ids;

            for(auto& name : relevant_read_names) {
                // Get the mappings for each read
                auto& mappings = graph.paths.get_path(name);
                for(auto& mapping : mappings) {
                    // Add in all the nodes that are visited
                    relevant_ids.insert(mapping.position().node_id());
                }
            }

            for(auto id : site.contents) {
                // Throw out all the IDs that are also used in the superbubble itself
                relevant_ids.erase(id);
            }

#ifdef debug
#pragma omp critical (cerr)
            cerr << relevant_read_names.size() << " reads visit an additional " << relevant_ids.size() << " nodes" << endl;
#endif

            // Make a vg graph with all the nodes used by the reads relevant to the
            // superbubble, but outside the superbubble itself.
            VG surrounding;
            for(auto id : relevant_ids) {
                // For all the IDs in the surrounding material

                if(min_recurrence != 0 && 
                        (!graph.paths.has_node_mapping(id) || 
                         graph.paths.get_node_mapping(id).size() < min_recurrence)) {
                    // Skip nodes in the graph that have too little support. In practice
                    // this means we'll restrict ourselves to supported, known nodes.
                    // TODO: somehow do the same for edges.
                    continue;
                }

                // Add each node and its edges to the new graph. Ignore dangling edges.
                // We'll keep edges dangling to the superbubble anchor nodes.
                surrounding.add_node(*graph.get_node(id));
                surrounding.add_edges(graph.edges_of(graph.get_node(id)));
            }

            for(auto& path : superbubble_paths) {
                // Now for each superbubble path, make a copy of that graph with it in
                VG allele_graph(surrounding);

                for(auto it = path.begin(); it != path.end(); ++it) {
                    // Add in every node on the path to the new allele graph
                    allele_graph.add_node(*(*it).node);

                    // Add in just the edge to the previous node on the path
                    if(it != path.begin()) {
                        // There is something previous on the path.
                        auto prev = it;
                        --prev;
                        // Make an edge
                        Edge path_edge;
                        // And hook it to the correct side of the last node
                        path_edge.set_from((*prev).node->id());
                        path_edge.set_from_start((*prev).backward);
                        // And the correct side of the next node
                        path_edge.set_to((*it).node->id());
                        path_edge.set_to_end((*it).backward);

                        assert(graph.has_edge(path_edge));

                        // And add it in
                        allele_graph.add_edge(path_edge);
                    }
                }

                // Get rid of dangling edges
                allele_graph.remove_orphan_edges();

#ifdef debug_verbose
#pragma omp critical (cerr)
                cerr << "Align to " << pb2json(allele_graph.graph) << endl;
#endif

                // Grab the sequence of the path we are trying the reads against, so we
                // can check for identity across the site and not just globally for the
                // read.
                auto path_seq = traversals_to_string(path);

                for(auto& name : relevant_read_names) {
                    // For every read that touched the superbubble, grab its original
                    // Alignment pointer.
                    Alignment* read = reads_by_name.at(name);

                    // Look to make sure it touches more than one node actually in the
                    // superbubble, or a non-start, non-end node. If it just touches the
                    // start or just touches the end, it can't be informative.
                    set<id_t> touched_set;
                    // Will this read be informative?
                    bool informative = false;            
                    for(size_t i = 0; i < read->path().mapping_size(); i++) {
                        // Look at every node the read touches
                        id_t touched = read->path().mapping(i).position().node_id();
                        if(site.contents.count(touched)) {
                            // If it's in the superbubble, keep it
                            touched_set.insert(touched);
                        }
                    }

                    if(touched_set.size() >= 2) {
                        // We touch both the start and end, or an internal node.
                        informative = true;
                    } else {
                        // Throw out the start and end nodes, if we touched them.
                        touched_set.erase(site.start.node->id());
                        touched_set.erase(site.end.node->id());
                        if(!touched_set.empty()) {
                            // We touch an internal node
                            informative = true;
                        }
                    }

                    if(!informative) {
                        // We only touch one of the start and end nodes, and can say nothing about the superbubble. Try the next read.
                        // TODO: mark these as ambiguous/consistent with everything (but strand?)
                        continue;
                    }

                    // If we get here, we know this read is informative as to the internal status of this superbubble.
                    Alignment aligned_fwd;
                    Alignment aligned_rev;
                    // We need a way to get graph node sizes to reverse these alignments
                    auto get_node_size = [&](id_t id) {
                        return graph.get_node(id)->sequence().size();
                    };
                    if(read->sequence().size() == read->quality().size()) {
                        // Re-align a copy to this graph (using quality-adjusted alignment).
                        // TODO: actually use quality-adjusted alignment
                        aligned_fwd = allele_graph.align(*read);
                        aligned_rev = allele_graph.align(reverse_complement_alignment(*read, get_node_size));
                    } else {
                        // If we don't have the right number of quality scores, use un-adjusted alignment instead.
                        aligned_fwd = allele_graph.align(*read);
                        aligned_rev = allele_graph.align(reverse_complement_alignment(*read, get_node_size));
                    }
                    // Pick the best alignment, and emit in original orientation
                    Alignment aligned = (aligned_rev.score() > aligned_fwd.score()) ? reverse_complement_alignment(aligned_rev, get_node_size) : aligned_fwd;

#ifdef debug
#pragma omp critical (cerr)
                    cerr << path_seq << " vs " << aligned.sequence() << ": " << aligned.score() << endl;

#endif

#ifdef debug_verbose
#pragma omp critical (cerr)
                    cerr << "\t" << pb2json(aligned) << endl;
#endif

                    // Compute the score per base. TODO: is this at all comparable
                    // between quality-adjusted and non-quality-adjusted reads?
                    double score_per_base = (double)aligned.score() / aligned.sequence().size();

                    // Save the score (normed per read base) and orientation
                    // We'll normalize the affinities later to enforce the max of 1.0.
                    Affinity affinity(score_per_base, aligned_rev.score() > aligned_fwd.score());

                    // Compute the unnormalized likelihood of the read given the allele graph.
                    if(read->sequence().size() == read->quality().size()) {
                        // Use the quality-adjusted default scoring system
                        affinity.likelihood_ln = quality_aligner.score_to_unnormalized_likelihood_ln(aligned.score());
                    } else {
                        // We will have aligned without quality adjustment, so interpret
                        // score in terms of the normal scoring parameters.
                        affinity.likelihood_ln = normal_aligner.score_to_unnormalized_likelihood_ln(aligned.score());
                    }

                    // Get the NodeTraversals for the winning alignment through the site.
                    auto read_traversal = get_traversal_of_site(graph, site, aligned.path());

                    if(affinity.is_reverse) {
                        // We really traversed this site backward. Flip it around.
                        read_traversal.reverse();
                        for(auto& item : read_traversal) {
                            // Flip around every traversal as well as reversing their order.
                            item = item.reverse();
                        }

                    }

                    // Decide we're consistent if the alignment's string across the site
                    // matches the string for the allele, anchored at the appropriate
                    // ends.

                    // Get the string this read spells out in its best alignment to this allele
                    auto seq = traversals_to_string(read_traversal);

                    // Now decide if the read's seq supports this path.
                    if(read_traversal.front() == site.start && read_traversal.back() == site.end) {
                        // Anchored at both ends.
                        // Need an exact match. Record if we have one or not.
                        affinity.consistent = (seq == path_seq);
                    } else if(read_traversal.front() == site.start) {
                        // Anchored at start only.
                        // seq needs to be a prefix of path_seq
                        auto difference = std::mismatch(seq.begin(), seq.end(), path_seq.begin());
                        // If the first difference is the past-the-end of the prefix, then it's a prefix
                        affinity.consistent = (difference.first == seq.end());
                    } else if(read_traversal.back() == site.end) {
                        // Anchored at end only.
                        // seq needs to be a suffix of path_seq
                        auto difference = std::mismatch(seq.rbegin(), seq.rend(), path_seq.rbegin());
                        // If the first difference is the past-the-rend of the suffix, then it's a suffix
                        affinity.consistent = (difference.first == seq.rend());
                    } else {
                        // This read doesn't touch either end. This might happen if the
                        // site is very large. Just assume it's consistent and let
                        // scoring work it out.
#pragma omp critical (cerr)
                        cerr << "Warning: realigned read " << aligned.sequence() << " doesn't touch either end of its site!" << endl;
                        affinity.consistent = true;
                    }

                    if(score_per_base < min_score_per_base) {
                        // Say we can't really be consistent with this if we have such a
                        // terrible score.
                        affinity.consistent = false;
                    }

                    // Grab the identity and save it for this read and superbubble path
                    to_return[read].push_back(affinity);

                }
            }

            for(auto& name : relevant_read_names) {
                // For every read that touched the superbubble, mark it consistent only
                // with its best-score alleles that don't mismatch in the allele.

                // So basically make everything that isn't normalized affinity 1.0
                // inconsistent if it wasn't already.

                // Grab its original Alignment pointer.
                Alignment* read = reads_by_name.at(name);

                // Which is the best affinity we can get while being consistent?
                double best_consistent_affinity = 0;
                // And which is the best affinity we can get overall?
                double best_affinity = 0;
                for(auto& affinity : to_return[read]) {
                    if(affinity.consistent && affinity.affinity > best_consistent_affinity) {
                        // Look for the max affinity found on anything already consistent
                        best_consistent_affinity = affinity.affinity;
                    }

                    if(affinity.affinity > best_affinity) {
                        // And the max affinity overall, for normalizing to correct
                        // affinity range of 0 to 1.
                        best_affinity = affinity.affinity;
                    }
                }

                for(auto& affinity : to_return[read]) {
                    if(affinity.affinity < best_consistent_affinity) {
                        // Mark all the Affinities that are worse than the best
                        // consistent one as not actually being consistent.
                        affinity.consistent = false;
                    }

                    if(best_affinity == 0) {
                        affinity.affinity = 1.0;
                    } else {
                        // Normalize to the best affinity overall being 1.
                        affinity.affinity /= best_affinity;
                    }
                }

            }

            // After scoring all the reads against all the versions of the superbubble,
            // return the affinities
            return to_return;
        }


    list<NodeTraversal> Genotyper::get_traversal_of_site(VG& graph, const Site& site, const Path& path) {

        // We'll fill this in
        list<NodeTraversal> to_return;

        for(size_t i = 0; i < path.mapping_size(); i++) {
            // Make a NodeTraversal version of the Mapping
            NodeTraversal traversal(graph.get_node(path.mapping(i).position().node_id()), path.mapping(i).position().is_reverse());

            if(site.contents.count(traversal.node->id())) {
                // We're inside the bubble. This is super simple when we have the contents!
                to_return.push_back(traversal);
            }

        }
        return to_return;
    }

    string Genotyper::traversals_to_string(const list<NodeTraversal>& path) {
        stringstream seq;
        for(auto& traversal : path) {
            // Stick in each sequence in order, with orientation
            seq << (traversal.backward ? reverse_complement(traversal.node->sequence()) : traversal.node->sequence());
        }
        return seq.str();
    }

    map<Alignment*, vector<Genotyper::Affinity> >
        Genotyper::get_affinities_fast(VG& graph,
                const map<string, Alignment*>& reads_by_name,
                const Site& site,
                const vector<list<NodeTraversal> >& superbubble_paths,
                bool allow_internal_alignments) {

            // We're going to build this up gradually, appending to all the vectors.
            map<Alignment*, vector<Affinity>> to_return;

            // What reads are relevant to this superbubble?
            set<string> relevant_read_names;

#ifdef debug
#pragma omp critical (cerr)
            cerr << "Superbubble contains " << site.contents.size() << " nodes" << endl;
#endif

            // Convert all the Paths used for alleles back to their strings.
            vector<string> allele_strings;
            for(auto& path : superbubble_paths) {
                // Convert all the Paths used for alleles back to their strings.
                allele_strings.push_back(traversals_to_string(path));
            }

            for(auto id : site.contents) {
                // For every node in the superbubble, what paths visit it?
                if(graph.paths.has_node_mapping(id)) {
                    auto& mappings_by_name = graph.paths.get_node_mapping(id);
                    for(auto& name_and_mappings : mappings_by_name) {
                        // For each path visiting the node
                        if(reads_by_name.count(name_and_mappings.first)) {
                            // This path is a read, so add the name to the set if it's not
                            // there already
                            relevant_read_names.insert(name_and_mappings.first);
                        }    
                    }
                }
            }

            for(auto name : relevant_read_names) {
                // For each relevant read, work out a string for the superbubble and whether
                // it's anchored on each end.

                // Make an Affinity to fill in
                Affinity base_affinity;

                // Get the NodeTraversals for this read through this site.
                auto read_traversal = get_traversal_of_site(graph, site, reads_by_name.at(name)->path());

                if(read_traversal.front() == site.end.reverse() || read_traversal.back() == site.start.reverse()) {
                    // We really traversed this site backward. Flip it around.
                    read_traversal.reverse();
                    for(auto& item : read_traversal) {
                        // Flip around every traversal as well as reversing their order.
                        item = item.reverse();
                    }

                    // We're on the reverse strand
                    base_affinity.is_reverse = true;
                }

                if(read_traversal.size() == 1 && (read_traversal.front() == site.start || read_traversal.back() == site.end)) {
                    // This read only touches the head or tail of the site, and so
                    // cannot possibly be informative.
                    cerr << "Non-informative site being removed " << endl
                        << read_traversal.front() << " to " << read_traversal.back() << endl;  
                    continue;
                }

                size_t total_supported = 0;

                // Get the string it spells out
                auto seq = traversals_to_string(read_traversal);

#ifdef debug
#pragma omp critical (cerr)
                cerr << "Consistency of " << reads_by_name.at(name)->sequence() << endl;
#endif

                // Now decide if the read's seq supports each path.
                for(auto& path_seq : allele_strings) {
                    // We'll make an affinity for this allele
                    Affinity affinity = base_affinity;
                    if(read_traversal.front() == site.start && read_traversal.back() == site.end) {
                        // Anchored at both ends.
                        // Need an exact match. Record if we have one or not.
                        affinity.consistent = (seq == path_seq);
                    } else if(read_traversal.front() == site.start) {
                        // Anchored at start only.
                        // seq needs to be a prefix of path_seq
                        auto difference = std::mismatch(seq.begin(), seq.end(), path_seq.begin());
                        // If the first difference is the past-the-end of the prefix, then it's a prefix
                        affinity.consistent = (difference.first == seq.end());
                    } else if(read_traversal.back() == site.end) {
                        // Anchored at end only.
                        // seq needs to be a suffix of path_seq
                        auto difference = std::mismatch(seq.rbegin(), seq.rend(), path_seq.rbegin());
                        // If the first difference is the past-the-rend of the suffix, then it's a suffix
                        affinity.consistent = (difference.first == seq.rend());
                    } else {
                        // This read doesn't touch either end.
#pragma omp critical (cerr)
                        cerr << "Warning: read doesn't touch either end of its site!" << endl;
                        if (allow_internal_alignments){

                        }
                    }

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t" << path_seq << " vs observed " << (read_traversal.front() == site.start) << " " << seq << " " << (read_traversal.back() == site.end) << ": " << affinity.consistent << endl;
#endif

                    // Fake a weight
                    affinity.affinity = (double)affinity.consistent;
                    to_return[reads_by_name.at(name)].push_back(affinity);

                    // Add in to the total if it supports this
                    total_supported += affinity.consistent;
                }

                if(total_supported == 0 && min_recurrence <= 1) {
                    // This is weird. Doesn't match anything and we had no excuse to remove alleles.
#pragma omp critical (cerr)
                    cerr << "Warning! Bubble sequence " << seq << " supports nothing!" << endl;
                }


            }


            // After scoring all the reads against all the versions of the superbubble,
            // return the affinities
            return to_return;
        }

    double Genotyper::get_genotype_log_likelihood(VG& graph, const Site& site, const vector<int>& genotype, const vector<pair<Alignment*, vector<Affinity>>>& alignment_consistency) {
        // For each genotype, calculate P(observed reads | genotype) as P(all reads
        // that don't support an allele from the genotype are mismapped or
        // miscalled) * P(all reads that do support alleles from the genotype ended
        // up apportioned across the alleles as they are)

        // This works out to the product over all reads that don't support either
        // alleles of 1 - ((1 - MAPQ) * (1 - P(bases called wrong)), times the
        // likelihood of the observed (censored by multi-support) read counts coming
        // from the alleles they support, and the strands they are observed on.

        // TODO: handle contamination like Freebayes

        // This is the log probability that all reads that don't support either allele in this genotype are wrong.
        double all_non_supporting_wrong = prob_to_logprob(1);

        // This is the log probability that all the reads that do support alleles in this genotype were drawn from the genotype they support.
        double all_supporting_drawn = prob_to_logprob(1);

#ifdef debug
        if(genotype.size() == 2) {
#pragma omp critical (cerr)
            cerr << "Calculating P(a" << genotype[0] << "/a" << genotype[1] << ")" << endl;
        }
#endif

        // We want to keep count of how many reads support each allele in each
        // orientation. TODO: we might be sort of double-counting reads that support
        // multiple alleles, because we'll be treating them as having independent
        // orientations per allele, when really they're very likely to be oriented
        // the same way.

        // Maps from int allele number (as appears in Genotype) to total reads
        // forward and reverse giving support.
        map<int, pair<int, int>> strand_count_by_allele_and_orientation;

        for(auto& read_and_consistency : alignment_consistency) {
            // For each read, work out if it supports a genotype we have or not.

            // Split out the alignment from its consistency flags
            auto& read = *read_and_consistency.first;
            auto& consistency = read_and_consistency.second;

            // How many of the alleles in our genotype is it consistent with?
            int consistent_alleles = 0;
            // We only count each allele in the genotype once.
            set<int> alleles_seen;
            for(int allele : genotype) {
                if(alleles_seen.count(allele)) {
                    // Counted up consistency with this allele already.
                    continue;
                }
                alleles_seen.insert(allele);

                // For each unique allele in the genotype...

                if(consistency.size() > allele && consistency.at(allele).consistent) {
                    // We're consistent with this allele
                    consistent_alleles++;
                    // And in this orientation
                    if(consistency.at(allele).is_reverse) {
                        // Consistent with reverse
                        strand_count_by_allele_and_orientation[allele].second++;
                    } else {
                        // Consistent with forward
                        strand_count_by_allele_and_orientation[allele].first++;
                    }
                }
            }

            auto read_qual = alignment_qual_score(graph, site, read);

#ifdef debug
#pragma omp critical (cerr)
            cerr << "Read (qual score " << read_qual << ") consistent with " << consistent_alleles
                << " genotype alleles observed." << endl;
#endif

            if(consistent_alleles == 0) {
                // This read is inconsistent with all the alleles in the genotype,
                // so, given the genotype, the read must be sequenced or mapped
                // wrong.

                double logprob_wrong;
                if(use_mapq) {
                    // Compute P(mapped wrong or called wrong) = P(not (mapped right and called right)) = P(not (not mapped wrong and not called wrong))
                    logprob_wrong = logprob_invert(logprob_invert(phred_to_logprob(read.mapping_quality())) +
                            logprob_invert(phred_to_logprob(read_qual)));
                } else {
                    // Compute P(called wrong).
                    logprob_wrong = phred_to_logprob(read_qual);
                }

#ifdef debug
#pragma omp critical (cerr)
                cerr << "P(wrong) = " << logprob_to_prob(logprob_wrong) << endl;
#endif
                all_non_supporting_wrong += logprob_wrong;
            } else {
                // This read is consistent with some of the alleles in the genotype,
                // so we must have drawn one of those alleles when sequencing.

                // We account for this in a framework where we consider the reads
                // indistinguishable, so ignore it for now.
            }

        }

        // Multiply in in the probability that the supporting reads all came from
        // the strands they are on.
        double strands_as_specified = prob_to_logprob(1);
        // Each strand is equally likely
        vector<double> probs_by_orientation = {0.5, 0.5};
        for(auto& kv : strand_count_by_allele_and_orientation) {
            // For the forward and reverse strand counts for all the alleles
            auto& forward_count = kv.second.first;
            auto& reverse_count = kv.second.second;

            // Convert to a vector to satisfy the multinomial PMF function.
            vector<int> obs = {forward_count, reverse_count};

            // Get the log probability under multinomial.

            // TODO: we're counting oriented reads supporting multiple alleles
            // multiple times. Maybe we should look at orientation overall somehow?
            // Or treat orientations of alleles as the thing we do our master
            // censored multinomial likelihood over?
            double logprob = multinomial_sampling_prob_ln(probs_by_orientation, obs);

#ifdef debug
            cerr << "Allele "  << kv.first << " supported by " << forward_count << " forward, "
                << reverse_count << " reverse (P=" << logprob_to_prob(logprob) << ")" << endl;
#endif
            strands_as_specified += logprob;


        }

        // Multiply in probability that the reads came from alleles they support,
        // treating reads as indistinguishable and using a multinomial/binomial
        // model.
        double alleles_as_specified = prob_to_logprob(1);
        if(genotype.size() == 2 && genotype.at(0) != genotype.at(1)) {
            // For diploid heterozygotes, we can easily handle multi-support as
            // censorship. We know that at least the reads that only support allele
            // 0 are from allele 0, and that at most the reads that only support
            // allele 0 and those that support both alleles all are from allele 0.
            // We end up summing over a normalized choose (since the success
            // probability is known and fixed at 1/2), as described in Frey and
            // Marrero (2008) "A surprising MLE for Interval-Censored Binomial
            // Data". Diploid homozygotes don't need any of this logic, and keep the
            // default probability of 1 above for the support distribution across
            // alleles.

            // Work out how many reads support allele 0 only.
            int first_only_reads = 0;
            // And how many support both
            int ambiguous_reads = 0;
            // And how many total reads there are (# of trials).
            int total_reads = 0;

            for(auto& read_and_consistency : alignment_consistency) {
                auto& consistency = read_and_consistency.second;

                if(consistency.size() <= max(genotype.at(0), genotype.at(1))) {
                    // No consistency info calculated for this useless uninformative
                    // read.
                    continue;
                }

                if(consistency.at(genotype.at(0)).consistent) {
                    // Read is consistent with first allele
                    if(consistency.at(genotype.at(1)).consistent) {
                        // And also second, so it's ambiguous
                        ambiguous_reads++;
#ifdef debug
                        cerr << "Ambiguous read: " << read_and_consistency.first->sequence() << endl;
                        for(int i = 0; i < consistency.size(); i++) {
                            cerr << "\t" << i << ": " << consistency[i].consistent << endl;
                        }
#endif
                    } else {
                        // And only first allele
                        first_only_reads++;
                    }
                    total_reads++;
                } else if(consistency.at(genotype.at(1)).consistent) {
                    // It's consistent with only the second allele.
                    total_reads++;
                }
                // Don't count reads inconsistent with the genotype in this analysis.
            }

            // Now do the likelihood. We know each atom will be weighted by the same
            // factor (for assigning all the reads at 50% probability) so we can
            // pull it out of the sum.
            double log_atom_weight = prob_to_logprob(0.5) * total_reads;

            // We calculate the probability of each atom, then sum, and then weight.
            vector<double> unweighted_atom_logprobs;

            for(int i = first_only_reads; i <= first_only_reads + ambiguous_reads; i++) {
                // For each possible actual number of reads from the first allele,
                // add in the probability of that atom.
                auto unweighted_atom_logprob = choose_ln(total_reads, i);
                unweighted_atom_logprobs.push_back(unweighted_atom_logprob);

#ifdef debug
#pragma omp critical (cerr)
                cerr << "P(" << i << " from first allele) = " << logprob_to_prob(unweighted_atom_logprob)
                    << " * " << logprob_to_prob(log_atom_weight) << " = " 
                    << logprob_to_prob(unweighted_atom_logprob + log_atom_weight) << endl;
#endif

            }

            // Weight all the atoms with the shared weight, and then multiply in (to
            // probability of 1) the probability of observing this range of possible
            // totals of reads from each of the two alleles.
            alleles_as_specified = (logprob_sum(unweighted_atom_logprobs) + log_atom_weight);

#ifdef debug
#pragma omp critical (cerr)
            cerr << "P(" << first_only_reads << " to "  << (first_only_reads + ambiguous_reads) << " from first allele) = "
                << logprob_to_prob(logprob_sum(unweighted_atom_logprobs)) << " * " << logprob_to_prob(log_atom_weight) 
                << " = " << logprob_to_prob(alleles_as_specified) << endl;
#endif

        } else if(genotype.size() != 2) {
#pragma omp critical (cerr)
            cerr << "Warning: not accounting for allele assignment likelihood in non-diploid genotype!" << endl;
        }

        // Now we've looked at all the reads, so AND everything together
        double total_logprob = all_non_supporting_wrong + all_supporting_drawn + strands_as_specified + alleles_as_specified;

#ifdef debug
        if(genotype.size() == 2) {
#pragma omp critical (cerr)
            cerr << "logP(a" << genotype[0] << "/a" << genotype[1] << ") = " << all_non_supporting_wrong << " + "
                << all_supporting_drawn << " + " << strands_as_specified << " + " << alleles_as_specified << " = "
                << total_logprob << endl;
        }
#endif

        return total_logprob;


    }

    double Genotyper::get_genotype_log_prior(const vector<int>& genotype) {
        assert(genotype.size() == 2);


        // Priors are boring: certain amount for het, inverse of that for everyone else
        if(genotype[0] != genotype[1]) {
            // This is a het!
            return het_prior_logprob;
        } else {
            // This is a homozygote. Much more common.
            return logprob_invert(het_prior_logprob);
        }
    }

    string Genotyper::get_qualities_in_site(VG& graph, const Site& site, const Alignment& alignment) {
        // We'll fill this in.
        stringstream to_return;

        // Are we currently in the site?
        bool in_site = false;
        // What NodeTraversal do we need to see to leave?
        NodeTraversal expected;

        // Where are we in the quality string?
        size_t quality_pos = 0;

        for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
            // For every mapping in the path in order
            auto& mapping = alignment.path().mapping(i);

            // What NodeTraversal is this?
            NodeTraversal traversal(graph.get_node(mapping.position().node_id()), mapping.position().is_reverse());

            if(!in_site) {
                // If we aren't in the site, we may be entering
                if(traversal == site.start) {
                    // We entered through the start
                    in_site = true;
                    // We'll leave at the end
                    expected = site.end;
                } else if(traversal == site.end.reverse()) {
                    // We entered through the end
                    in_site = true;
                    // We'll leave when we hit the start in reverse
                    expected = site.start.reverse();
                }
            }

            for(size_t j = 0; j < mapping.edit_size(); j++) {
                // For every edit
                auto& edit = mapping.edit(j);

                if(in_site && mapping.position().node_id() != site.start.node->id()
                        && mapping.position().node_id() != site.end.node->id()) {
                    // We're in the site but not on the start or end nodes.
                    // TODO: qualities for a deletion/insertion?
                    for(size_t k = 0; k < edit.to_length(); k++) {
                        // Take each quality value from the edit and add it to our collection to return
                        if(quality_pos >= alignment.quality().size()) {
                            // If we've run out of quality values, give back no
                            // qualities, because base qualities aren't really being
                            // used.
                            return "";
                        }
                        to_return << (char)alignment.quality().at(quality_pos);
                        quality_pos++;
                    }
                } else {
                    // Skip this edit's qualities 
                    quality_pos += edit.to_length();
                }
            }

            if(in_site && traversal == expected) {
                // This was the node we were supposed to leave the site at.
                in_site = false;
            }
        }

        return to_return.str();

    }

    Locus Genotyper::genotype_site(VG& graph,
            const Site& site,
            const vector<list<NodeTraversal>>& superbubble_paths,
            const map<Alignment*, vector<Affinity>>& affinities) {

        // Freebayes way (improved with multi-support)

        // We're going to populate this locus
        Locus to_return;

        for(auto& path : superbubble_paths) {
            // Convert each allele to a Path and put it in the locus
            *to_return.add_allele() = path_from_node_traversals(path);
        }

#ifdef debug
#pragma omp critical (cerr)
        cerr << "Looking between " << site.start << " and " << site.end << endl;
#endif

        // We'll fill this in with the alignments for this site and their consistency-with-alleles flags.
        vector<pair<Alignment*, vector<Affinity>>> alignment_consistency;

        // We fill this in with totals of reads supporting alleles
        vector<int> reads_consistent_with_allele(superbubble_paths.size(), 0);
        // And this with the same thing split out by forward and reverse strand
        vector<pair<int, int>> strand_support_for_allele(superbubble_paths.size(), make_pair(0, 0));

        // We'll store affinities by read name and allele here, for printing later.
        map<string, vector<Affinity>> debug_affinities;

        // We track overall forward and reverse strand reads, of reads that
        // support any allele.
        size_t overall_forward_reads = 0;
        size_t overall_reverse_reads = 0;

        for(auto& alignment_and_affinities : affinities) {
            // We need to clip down to just the important quality values        
            Alignment& alignment = *alignment_and_affinities.first;

            // Hide all the affinities where we can pull them later
            debug_affinities[alignment.name()] = alignment_and_affinities.second;

            // Affinities already know whether they are consistent with an allele. Don't second-guess them.
            // Fine even with alignment; no read we embedded should ever have non-perfect identity.

            // We'll set these if the read supports anything in a forward or reverse
            // orientation.
            bool is_forward = false;
            bool is_reverse = false;

            for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
                // Count up reads consistent with each allele
                if(alignment_and_affinities.second.at(i).consistent) {
                    // This read is consistent with this allele
                    reads_consistent_with_allele[i]++;
                    if(alignment_and_affinities.second.at(i).is_reverse) {
                        // It is on the reverse strand
                        strand_support_for_allele[i].second++;
                        is_reverse = true;
                    } else {
                        // It is on the forward strand
                        strand_support_for_allele[i].first++;
                        is_forward = true;
                    }
                }
            }

            if(is_forward) {
                if(is_reverse) {
                    // This is weird
#pragma omp critical (cerr)
                    cerr << "Warning! Read supports alleles as both forward and reverse!" << endl;
                    // Just call it forward
                }
                // This read supports an allele forward, so call it a forward read for the site
                overall_forward_reads++;
            } else if(is_reverse) {
                // This read supports an allele reverse, so call it a reverse read for the site
                overall_reverse_reads++;
            } else if(min_recurrence <= 1) {
                // Reads generally ought to support at least one allele, unless we
                // have weird softclips or they were for elided non-recurrent
                // alleles.
#pragma omp critical (cerr)
                cerr << "Warning! Read supports no alleles!" << endl;
            }

            // Save the alignment and its affinities, which we use to get GLs.
            alignment_consistency.push_back(alignment_and_affinities);

        }

#ifdef debug
        for(int i = 0; i < superbubble_paths.size(); i++) {
            // Build a useful name for the allele
            stringstream allele_name;
            for(auto& traversal : superbubble_paths[i]) {
                allele_name << traversal.node->id() << ",";
            }
#pragma omp critical (cerr)
            {
                cerr << "a" << i << "(" << allele_name.str() << "): " << reads_consistent_with_allele[i] << "/" << affinities.size() << " reads consistent" << endl;
                for(auto& read_and_consistency : alignment_consistency) {
                    if(read_and_consistency.second.size() > i && 
                            read_and_consistency.second[i].consistent &&
                            read_and_consistency.first->sequence().size() < 30) {
                        // Dump all the short consistent reads
                        cerr << "\t" << read_and_consistency.first->sequence() << " " << debug_affinities[read_and_consistency.first->name()][i].affinity << endl;
                    }
                }
            }
        }
#endif

        // We'll go through all the genotypes, fill in their probabilities, put them
        // in here, and then sort them to find the best.
        vector<Genotype> genotypes_sorted;

        for(int allele1 = 0; allele1 < superbubble_paths.size(); allele1++) {
            // For each first allele in the genotype
            for(int allele2 = 0; allele2 <= allele1; allele2++) {
                // For each second allele so we get all order-independent combinations

                // Make the combo
                vector<int> genotype_vector = {allele1, allele2};

                // Compute the log probability of the data given the genotype
                double log_likelihood = get_genotype_log_likelihood(graph, site, genotype_vector, alignment_consistency);

                // Compute the prior
                double log_prior = get_genotype_log_prior(genotype_vector);

                // Apply Bayes Rule
                double log_posterior_unnormalized = log_likelihood + log_prior;

#ifdef debug
#pragma omp critical (cerr)
                {
                    cerr << "P(obs | a" << allele1 << "/a" << allele2 << ") = " << logprob_to_prob(log_likelihood) <<
                        " (" << log_likelihood << ")" << endl;
                    cerr << "P(a" << allele1 << "/a" << allele2 << ") = " << logprob_to_prob(log_prior) <<
                        " (" << log_prior << ")" << endl;
                    cerr << "P(a" << allele1 << "/a" << allele2 << " | obs) * P(obs) = " <<
                        logprob_to_prob(log_posterior_unnormalized) << " (" << log_posterior_unnormalized << ")" << endl;
                }
#endif

                // Fill in the actual Genotype object
                Genotype genotype;
                genotype.set_log_likelihood(log_likelihood);
                genotype.set_log_prior(log_prior);
                genotype.set_log_posterior(log_posterior_unnormalized);

                for(auto allele_id : genotype_vector) {
                    // Copy over all the indexes of alleles in the genotype
                    genotype.add_allele(allele_id);
                }

                // Put it in to sort
                genotypes_sorted.push_back(genotype);
            }
        }

        // Sort the genotypes in order of descending log posterior.
        sort(genotypes_sorted.begin(), genotypes_sorted.end(), [](const Genotype& a, const Genotype& b) {
                return a.log_posterior() > b.log_posterior();
                });

        for(size_t i = 0; i < superbubble_paths.size(); i++) {
            // For each allele, make a support
            Support* support = to_return.add_support();
            // Set forward and reverse depth
            support->set_forward(strand_support_for_allele[i].first);
            support->set_reverse(strand_support_for_allele[i].second);
        }

        for(auto& genotype : genotypes_sorted) {
            // Add a genotype to the Locus for every one we looked at, in order by descending posterior
            *to_return.add_genotype() = genotype;
        }

        // Set up total support for overall depth
        Support* overall_support = to_return.mutable_overall_support();
        overall_support->set_forward(overall_forward_reads);
        overall_support->set_reverse(overall_reverse_reads);

        // Now we've populated the genotype so return it.
        return to_return;
    }

    void Genotyper::write_vcf_header(std::ostream& stream, const std::string& sample_name, const std::string& contig_name, size_t contig_size) {
        stream << "##fileformat=VCFv4.2" << std::endl;
        stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << std::endl;
        stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
        stream << "##INFO=<ID=XSBB,Number=1,Type=Integer,Description=\"Superbubble Bases\">" << std::endl;
        stream << "##INFO=<ID=XSBN,Number=1,Type=Integer,Description=\"Superbubble Nodes\">" << std::endl;
        stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
        stream << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
        stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
        stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << std::endl;
        // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
        stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << std::endl;
        if(!contig_name.empty()) {
            // Announce the contig as well.
            stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << std::endl;
        }
        stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;
    }

    vcflib::VariantCallFile* Genotyper::start_vcf(std::ostream& stream, const PathIndex& index, const string& sample_name, const string& contig_name, size_t contig_size) {
        // Generate a vcf header. We can't make Variant records without a
        // VariantCallFile, because the variants need to know which of their
        // available info fields or whatever are defined in the file's header, so
        // they know what to output.
        // Handle length override if specified.
        std::stringstream headerStream;
        write_vcf_header(headerStream, sample_name, contig_name, contig_size > 0 ? contig_size : index.sequence.size());

        // Load the headers into a new VCF file object
        vcflib::VariantCallFile* vcf = new vcflib::VariantCallFile();
        std::string headerString = headerStream.str();
        assert(vcf->openForOutput(headerString));

        // Spit out the header
        stream << headerStream.str();

        // Give back the created VCF
        return vcf;
    }

    vector<vcflib::Variant>
        Genotyper::locus_to_variant(VG& graph,
                const Site& site,
                const PathIndex& index,
                vcflib::VariantCallFile& vcf,
                const Locus& locus,
                const string& sample_name) {

            // Make a vector to fill in
            vector<vcflib::Variant> to_return;

            // Make a new variant
            vcflib::Variant variant;
            // Attach it to the VCF
            variant.setVariantCallFile(vcf);
            // Fake the quality
            variant.quality = 0;

            // Make sure we have stuff
            if(locus.allele_size() == 0) {
                throw runtime_error("Can't turn an empty genotype into VCF");
            }
            if(locus.allele(0).mapping_size() == 0) {
                throw runtime_error("Can't turn an empty allele into VCF");
            }

            // Get the superbubble    
            auto first_id = site.start.node->id();
            auto last_id = site.end.node->id();

            if(!index.by_id.count(first_id) || !index.by_id.count(last_id)) {
                // We need to be anchored to the primary path to make a variant
#pragma omp critical (cerr)
                cerr << "Warning: Superbubble endpoints not on reference!" << endl;
                // If not return no variant
                return to_return;
            }

            // Compute the reference region occupied by the site, accounting for
            // orientation.
            auto bounds = get_site_reference_bounds(site, index);

            // Where does this bubble start and end in the reference?
            auto referenceIntervalStart = bounds.first.first;
            auto referenceIntervalPastEnd = bounds.first.second;

            // Is this bubble articulated backwards relative to the reference?
            bool site_is_reverse = bounds.second;
            if(site_is_reverse) {
                // Make sure our first and last IDs are actually accurate.
                swap(first_id, last_id);
            }

            // Get the string for the reference allele
            string ref_string = index.sequence.substr(referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);

            // And make strings for all the locus's alleles
            vector<string> allele_strings;

            for(size_t i = 0; i < locus.allele_size(); i++) {
                // Get the string for each allele
                string allele = allele_to_string(graph, locus.allele(i));
                if(site_is_reverse) {
                    // Flip the alleles to match the reference orientation if necessary.
                    allele = reverse_complement(allele);
                }
                allele_strings.push_back(allele);
            }

            // See if any alleles are empty
            bool empty_alleles = ref_string.empty();
            for(auto& allele : allele_strings) {
                if(allele == "") {
                    empty_alleles = true;
                }
            }

            // Fix them up
            if(empty_alleles) {
                // Grab the character before our site
                string prefix = index.sequence.substr(referenceIntervalStart - 1, 1);
                for(auto& allele : allele_strings) {
                    // Prepend it to every allele
                    allele = prefix + allele;
                }
                // Also prepend it to the reference string
                ref_string = prefix + ref_string;

                // Budge the variant over
                referenceIntervalStart--;
            }

            // Trim fixed characters off the right
            while(ref_string.size() > 1) {
                // Grab the character at the end of the ref sequence
                char fixed = ref_string.back();

                // Set this to false if not all the alt strings also have this character
                // at the end, free to be chopped off.
                bool all_have_character = true;
                for(auto& allele_string : allele_strings) {
                    if(allele_string.size() <= 1 || allele_string.back() != fixed) {
                        // String is too short to trim or ends in the wrong character
                        all_have_character = false;
                        break;
                    }
                }

                if(all_have_character) {
                    // Trim it off
                    ref_string.pop_back();
                    for(auto& allele_string : allele_strings) {
                        allele_string.pop_back();
                    }

                    // Record we budged the end of the interval left.
                    referenceIntervalPastEnd--;
                } else {
                    // Done trimming
                    break;
                }
            }

            // Trim fixed characters off the left
            while(ref_string.size() > 1) {
                // Grab the character at the start of the ref sequence
                char fixed = ref_string.front();

                // Set this to false if not all the alt strings also have this character
                // at the start, free to be chopped off.
                bool all_have_character = true;
                for(auto& allele_string : allele_strings) {
                    if(allele_string.size() <= 1 || allele_string.front() != fixed) {
                        // String is too short to trim or starts with the wrong character
                        all_have_character = false;
                        break;
                    }
                }

                if(all_have_character) {
                    // Trim it off
                    // TODO: this is going to be O(n^2)
                    ref_string.erase(0, 1);
                    for(auto& allele_string : allele_strings) {
                        allele_string.erase(0, 1);
                    }

                    // Record that we budged the reference sequence start right.
                    referenceIntervalStart++;
                } else {
                    // Done trimming
                    break;
                }
            }

            // Make the ref allele
            create_ref_allele(variant, ref_string);

            // Make a vector of supports by assigned VCF alt number
            vector<Support> support_by_alt;

            // This maps from locus allele index to VCF record alt number
            vector<int> allele_to_alt;

            // This records the max alt number used. We might have more alts than
            // alleles if the reference was never considered as an allele.
            int max_alt_number = 0;

            for(size_t i = 0; i < locus.allele_size(); i++) {
                // For each allele

                // Add it/find its number if it already exists (i.e. is the ref)
                int alt_number = add_alt_allele(variant, allele_strings[i]);

                // See if it's a new max
                max_alt_number = max(max_alt_number, alt_number);

                // Remember what VCF alt number it got
                allele_to_alt.push_back(alt_number);

                if(i < locus.support_size()) {
                    if(support_by_alt.size() <= alt_number) {
                        // Make sure we have a slot to put the support in
                        support_by_alt.resize(alt_number + 1);
                    }
                    // Put it there
                    support_by_alt[alt_number] = locus.support(i);
                }
            }

            // Get the best genotype
            assert(locus.genotype_size() > 0);
            const Genotype& best_genotype = locus.genotype(0);
            // And its support
            assert(locus.support_size() > 0);
            const Support& best_support = locus.support(0);
            // TODO: right now we only handle diploids
            assert(best_genotype.allele_size() == 2);

            // Compose the ML genotype
            variant.format.push_back("GT");
            auto& genotype_out = variant.samples[sample_name]["GT"];
            // Translate each allele to a VCF alt number, and put them in a string with the right separator
            genotype_out.push_back(to_string(allele_to_alt[best_genotype.allele(0)]) + (best_genotype.is_phased() ? "|" : "/")  +
                    to_string(allele_to_alt[best_genotype.allele(1)]));

            // Make sure that the called alleles have sufficient support on each strand
            for(size_t i = 0; i < best_genotype.allele_size(); i++) {
                // Check each allele marked present
                if(locus.support(best_genotype.allele(i)).forward() < min_consistent_per_strand ||
                        locus.support(best_genotype.allele(i)).reverse() < min_consistent_per_strand) {
                    // If there's not enough support for that allele in an orientation, skip the site. 

#pragma omp critical (cerr)
                    cerr << "Warning: dropping locus from VCF due to insufficient per-strand support "
                        << locus.support(best_genotype.allele(i)).forward() << ", " 
                        << locus.support(best_genotype.allele(i)).reverse() << endl;

                    return to_return;
                }
            }

            // Put the total depth overall (not double-counting)
            string depth_string = std::to_string(locus.overall_support().forward() + locus.overall_support().reverse());
            variant.format.push_back("DP");
            variant.samples[sample_name]["DP"].push_back(depth_string);
            variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth

            // Also the site statistics
            // Superbubble bases
            size_t superbubble_bases = 0;
            for(auto& node_id : site.contents) {
                superbubble_bases += graph.get_node(node_id)->sequence().size();
            }
            variant.info["XSBB"].push_back(to_string(superbubble_bases));
            // Superbubble nodes
            variant.info["XSBN"].push_back(to_string(site.contents.size()));

            variant.format.push_back("GQ");
            if(locus.genotype_size() > 1) {
                // Compute a quality from the difference between the best and second-
                // best genotype posteriors. Really this should be:

                // P(genotype is wrong) = sum(P(genotype) over other genotypes) / sum(P(genotype) over all genotypes)

                // When best genotype is much more probable than second best, which is
                // much more probable than all the rest, this approximation woks well.
                variant.samples[sample_name]["GQ"].push_back(to_string(
                            logprob_to_phred(locus.genotype(1).log_posterior() - best_genotype.log_posterior())));
            } else {
                // This is very unlikely to be wrong. It can only be wrong if all the
                // reads in support of ref missed the haplotype on which an alt is.
                // TODO: this isn't exactly right; we should somehow account here for
                // reads we threw out for being on non-recurrent alleles...
                int total_reads = best_support.forward() + best_support.reverse();
                // Compute the likelihood that we missed everything, multiply it by a
                // prior of 5% for just being completely wrong, and treat that as the
                // posterior for the second best genotype.
                double all_missed_logprob = prob_to_logprob(0.5) * total_reads + prob_to_logprob(0.05);
                variant.samples[sample_name]["GQ"].push_back(to_string(logprob_to_phred(all_missed_logprob - best_genotype.log_posterior())));
            }

            // Compose the allele-specific depth
            variant.format.push_back("AD");
            for(auto& support : support_by_alt) {
                // Add the forward and reverse support together and use that for AD for the allele.
                variant.samples[sample_name]["AD"].push_back(to_string(support.forward() + support.reverse()));
            }

            // Do support for ref and alt alleles by strand
            variant.format.push_back("SB");
            for(auto& support : support_by_alt) {
                // Add the forward and reverse support in sequence, for ref and all the alts.
                // TODO: make this really only have the alt that's called.
                variant.samples[sample_name]["SB"].push_back(to_string(site_is_reverse ? support.reverse() : support.forward()));
                variant.samples[sample_name]["SB"].push_back(to_string(site_is_reverse ? support.forward() : support.reverse()));
            }


            // Work out genotype log likelihoods
            // Make a vector to shuffle them into VCF order. Fill it with inf for impossible genotypes.
            vector<double> log_likelihoods((max_alt_number * (max_alt_number + 1)) / 2 + max_alt_number + 1,
                    numeric_limits<double>::infinity());
            for(size_t i = 0; i < locus.genotype_size(); i++) {
                // For every genotype, calculate its VCF-order index based on the VCF allele numbers

                // TODO: we can only do diploids
                assert(locus.genotype(i).allele_size() == 2);

                // We first need the low and high alt numbers
                size_t low_alt = allele_to_alt.at(locus.genotype(i).allele(0));
                size_t high_alt = allele_to_alt.at(locus.genotype(i).allele(1));
                if(low_alt > high_alt) {
                    // Flip them to be the right way around
                    swap(low_alt, high_alt);
                }

                // Compute the position as (sort of) specified in the VCF spec
                size_t index = (high_alt * (high_alt + 1)) / 2 + low_alt;
                // Store the log likelihood
                log_likelihoods.at(index) = locus.genotype(i).log_likelihood();
#ifdef debug
#pragma omp critical (cerr)
                cerr << high_alt << "/" << low_alt << ": " << index << " = " << pb2json(locus.genotype(i)) << endl;
#endif
            }

            variant.format.push_back("PL");
            for(auto& log_likelihood : log_likelihoods) {
                // Add all the likelihood strings, normalizing against the best
                // TODO: the best may not actually be the chosen genotype, because we genotype on posteriors.
                variant.samples[sample_name]["PL"].push_back(to_string(logprob_to_phred(log_likelihood - best_genotype.log_likelihood())));
            }

            // Set the variant position (now that we have budged it left if necessary)
            variant.position = referenceIntervalStart + 1;

            // Return the variant, since we managed to make it
            to_return.push_back(variant);
            return to_return;

        }

    void Genotyper::report_site(const Site& site, const PathIndex* index) {
        if(site.contents.size() == 2) {
            // Skip degenerate sites
            return;
        }

        // Remember that we have the site
#pragma omp critical (all_sites)
        all_sites.insert(&site);

        if(index != nullptr) {
            // We have an index of the reference.

            // Figure out its reference length and log it.
            auto bounds = get_site_reference_bounds(site, *index);

            if(bounds.first.first == -1) {
                // It's not really on that path
                return;
            }

            int64_t length = bounds.first.second - bounds.first.first;

#pragma omp critical (site_reference_length_histogram)
            site_reference_length_histogram[length]++;

        }

    }

    void Genotyper::report_site_traversal(const Site& site, const string& name) {
        if(site.contents.size() == 2) {
            // Skip degenerate sites
            return;
        }

        // Mark this read as traversing this site
#pragma omp critical (site_traversals)
        site_traversals[&site].insert(name);
    }

    void Genotyper::print_statistics(ostream& out) {
        // Dump our stats to the given ostream.

        out << "Statistics:" << endl;
        out << "Number of Non-Degenerate Sites: " << all_sites.size() << endl;

        // How many sites were actually traversed by reads?
        size_t sites_traversed = 0;
        for(const Site* site : all_sites) {
            // For every site
            if(site_traversals.count(site) && site_traversals.at(site).size() > 0) {
                // If it has a set of read names and the set is nonempty, it was traversed
                sites_traversed++;
            }
        }
        out << "Sites traversed by reads: " << sites_traversed << endl;

        // How many sites are on the reference? Only those that have defined lengths
        size_t sites_on_reference = 0;
        for(auto& length_and_count : site_reference_length_histogram) {
            sites_on_reference += length_and_count.second;
        }
        out << "Sites on reference: " << sites_on_reference << endl;

        // What's the length distribution?
        out << "Site length distribution: " << endl;
        for(auto& length_and_count : site_reference_length_histogram) {
            // Dump length and count as a TSV bit.
            out << length_and_count.first << "\t" << length_and_count.second << endl;
        }

    }

}


