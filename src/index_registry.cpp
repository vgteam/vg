// index_registry.cpp: index registry system implementation

#include "index_registry.hpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <cctype>
#include <cstdio>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/odgi.hpp>
#include <xg.hpp>
#include <gbwt/variants.h>
#include <vg/io/vpkg.hpp>
#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>

#include "vg.hpp"
#include "vg_set.hpp"
#include "handle.hpp"
#include "utility.hpp"
#include "constructor.hpp"
#include "hash_map.hpp"
#include "haplotype_indexer.hpp"
#include "phase_unfolder.hpp"
#include "gbwt_helper.hpp"
#include "kmer.hpp"
#include "transcriptome.hpp"
#include "integrated_snarl_finder.hpp"
#include "min_distance.hpp"
#include "gfa.hpp"

#include "io/save_handle_graph.hpp"

#include "algorithms/gfa_to_handle.hpp"
#include "algorithms/prune.hpp"

//#define debug_index_registry
//#define debug_index_registry_setup
//#define debug_index_registry_recipes
//#define debug_index_registry_path_state

namespace std {
    
/// Convert IndexNames to strings, without defining it for all things sharing
/// the same underlying type.
static string to_string(const vg::IndexGroup& name) {
    stringstream ss;
    for (auto it = name.begin(); it != name.end(); ++it) {
        if (it != name.begin()) {
            ss << " + ";
        }
        ss << *it;
    }
    return ss.str();
}
    
}

namespace vg {


IndexingParameters::MutableGraphImplementation IndexingParameters::mut_graph_impl = HashGraph;
int IndexingParameters::max_node_size = 32;
int IndexingParameters::pruning_max_node_degree = 128;
int IndexingParameters::pruning_walk_length = 24;
int IndexingParameters::pruning_max_edge_count = 3;
int IndexingParameters::pruning_min_component_size = 33;
int IndexingParameters::gcsa_initial_kmer_length = gcsa::Key::MAX_LENGTH;
int IndexingParameters::gcsa_doubling_steps = gcsa::ConstructionParameters::DOUBLING_STEPS;
int IndexingParameters::gbwt_insert_batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
int IndexingParameters::gbwt_sampling_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
bool IndexingParameters::bidirectional_haplo_tx_gbwt = false;
string IndexingParameters::gff_feature_name = "exon";
string IndexingParameters::gff_transcript_tag = "transcript_id";
bool IndexingParameters::verbose = false;

// return file size in bytes
int64_t get_file_size(const string& filename) {
    // get the file size
    ifstream infile(filename);
    infile.seekg(0, ios::end);
    return infile.tellg();
}

// quickly guess the number of variants in a VCF file based on the filesize
double approx_num_vars(const string& vcf_filename) {
    
    // read the file to get the number of samples
    htsFile* vcf_file = hts_open(vcf_filename.c_str(),"rb");
    bcf_hdr_t* header = bcf_hdr_read(vcf_file);
    size_t num_samples = bcf_hdr_nsamples(header);
    bcf_hdr_destroy(header);
    hts_close(vcf_file);
    
    int64_t file_size = get_file_size(vcf_filename);
    
    bool is_gzipped = false;
    if (vcf_filename.size() > 2 && vcf_filename.substr(vcf_filename.size() - 3, 3) == ".gz") {
        is_gzipped = true;
    }
    
    // TODO: bcf coefficient
    // a shitty regression that Jordan fit on the human autosomes, gives a very rough
    // estimate of the number of variants contained in a VCF
    double coef = is_gzipped ? 12.41 : 0.2448;
    return (coef * file_size) / num_samples;
}

// the ratio to the HashGraph memory usage, as estimated by a couple of graphs
// that Jordan had laying around when he wrote this
// in case it's useful later: XG is ~0.231
double format_multiplier() {
    switch (IndexingParameters::mut_graph_impl) {
        case IndexingParameters::HashGraph:
            return 1.0;
        case IndexingParameters::ODGI:
            return 0.615;
        case IndexingParameters::PackedGraph:
            return 0.187;
        case IndexingParameters::VG:
            return 2.91;
        default:
            cerr << "error:[IndexRegistry] unrecognized mutable graph implementation format" << endl;
            exit(1);
            return 0.0;
    }
}

// approximate the memory of a graph that would be constructed with all of these
// inputs
int64_t approx_graph_memory(const vector<string>& fasta_filenames, const vector<string>& vcf_filenames) {

    // compute the size of the reference and the approximate number of
    // variants in the VCF
    int64_t ref_size = 0;
    double num_vars = 0.0;
    for (const auto& fasta : fasta_filenames) {
        ref_size += get_file_size(fasta);
    }
    for (const auto& vcf : vcf_filenames) {
        num_vars += approx_num_vars(vcf);
    }
        
    // estimates made by regressing the memory usage of a linear reference on the size
    // of the FASTA and then regressing the difference in memory usage between the linear
    // reference and 1000GP graph on the number of variants in the VCF, all using human
    // chromosomes with outliers removed
    double linear_memory = 30.4483 * ref_size;
    double var_memory = 2242.90 * num_vars;
    double hash_graph_memory_usage = linear_memory + var_memory;
    return hash_graph_memory_usage * format_multiplier();
}

int64_t approx_graph_memory(const string& fasta_filename, const string& vcf_filename) {
    return approx_graph_memory(vector<string>(1, fasta_filename), vector<string>(1, vcf_filename));
}

// estimate the amount of memory of a graph constructed from all of these inputs
int64_t approx_graph_memory(const vector<string>& gfa_filenames) {
    
    int64_t total_size = 0;
    for (const auto& gfa : gfa_filenames) {
        total_size += get_file_size(gfa);
    }
    
    // factor estimated by regression on 1000GP graphs of human chromosomes
    int64_t hash_graph_memory_usage = 13.17 * total_size;
    return hash_graph_memory_usage * format_multiplier();
}

int64_t approx_graph_memory(const string& gfa_filename) {
    return approx_graph_memory(vector<string>(1, gfa_filename));
}

// returns true if the GTF/GFF has any non-header lines
bool transcript_file_nonempty(const string& transcripts) {
    ifstream strm(transcripts);
    string line;
    while (strm.good()) {
        getline(strm, line);
        if (!line.empty() && line[0] != '#') {
            return true;
        }
        line.clear();
    }
    return false;
}

IndexRegistry VGIndexes::get_vg_index_registry() {
    
    IndexRegistry registry;
    
    /*********************
     * Register all of the VG indexes and input files
     ***********************/
    
    // TODO: we need separate suffixes for co-created indexes
    
    /// Data files
    registry.register_index("Reference FASTA", "fasta");
    registry.register_index("VCF", "vcf");
    registry.register_index("VCF w/ Phasing", "phased.vcf");
    registry.register_index("Insertion Sequence FASTA", "insertions.fasta");
    registry.register_index("Reference GFA", "gfa");
    registry.register_index("GTF/GFF", "gff");
    
    /// True indexes
    registry.register_index("VG", "vg");
    registry.register_index("VG w/ Variant Paths", "varpaths.vg");
    registry.register_index("Pruned VG", "pruned.vg");
    registry.register_index("Spliced VG", "spliced.vg");
    registry.register_index("Spliced VG w/ Variant Paths", "spliced.varpaths.vg");
    registry.register_index("Spliced VG w/ Transcript Paths", "spliced.txpaths.vg");
    registry.register_index("Pruned Spliced VG", "spliced.pruned.vg");
    
    registry.register_index("XG", "xg");
    registry.register_index("Spliced XG", "spliced.xg");
    
    registry.register_index("Unjoined Transcript Origin Table", "unjoined.txorigin.tsv");
    registry.register_index("Transcript Origin Table", "txorigin.tsv");
    
    registry.register_index("MaxNodeID", "maxid.txt");
    registry.register_index("Spliced MaxNodeID", "spliced.maxid.txt");
    registry.register_index("Unfolded NodeMapping", "mapping");
    registry.register_index("Haplotype-Pruned VG", "haplopruned.vg");
    registry.register_index("Unfolded Spliced NodeMapping", "spliced.mapping");
    registry.register_index("Haplotype-Pruned Spliced VG", "spliced.haplopruned.vg");
    registry.register_index("GCSA", "gcsa");
    registry.register_index("LCP", "gcsa.lcp");
    registry.register_index("Spliced GCSA", "spliced.gcsa");
    registry.register_index("Spliced LCP", "spliced.gcsa.lcp");
    
    registry.register_index("GBWT", "gbwt");
    registry.register_index("Spliced GBWT", "spliced.gbwt");
    registry.register_index("Haplotype-Transcript GBWT", "haplotx.gbwt");
    
    registry.register_index("Snarls", "snarls");
    registry.register_index("Spliced Snarls", "spliced.snarls");
    
    registry.register_index("Distance Index", "dist");
    registry.register_index("Spliced Distance Index", "spliced.dist");
    
    
    /*********************
     * A few handy lambda functions
     ***********************/
    
    auto init_in = [](ifstream& in, const string& name) {
        in.open(name);
        if (!in) {
            cerr << "error:[IndexRegistry] could not open input file " << name << endl;
            exit(1);
        }
    };
    auto init_out = [](ofstream& out, const string& name) {
        out.open(name);
        if (!out) {
            cerr << "error:[IndexRegistry] could not write output to " << name << endl;
            exit(1);
        }
    };
    auto init_in_out = [](fstream& strm, const string& name) {
        strm.open(name);
        if (!strm) {
            cerr << "error:[IndexRegistry] could not open " << name << endl;
            exit(1);
        }
    };
    
    auto init_mutable_graph = [&]() -> unique_ptr<MutablePathDeletableHandleGraph> {
        unique_ptr<MutablePathDeletableHandleGraph> graph;
        switch (IndexingParameters::mut_graph_impl) {
            case IndexingParameters::HashGraph:
                graph = make_unique<bdsg::HashGraph>();
                break;
            case IndexingParameters::ODGI:
                graph = make_unique<bdsg::ODGI>();
                break;
            case IndexingParameters::PackedGraph:
                graph = make_unique<bdsg::PackedGraph>();
                break;
            case IndexingParameters::VG:
                graph = make_unique<VG>();
                break;
            default:
                cerr << "error:[IndexRegistry] unrecognized mutable graph implementation format" << endl;
                exit(1);
                break;
        }
        return graph;
    };
    
    
    
    
    /*********************
     * Register all recipes
     ***********************/
    
    ////////////////////////////////////
    // VCF Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering VCF recipes" << endl;
#endif
    
    // alias a phased VCF as an unphased one
    registry.register_recipe({"VCF"}, {"VCF w/ Phasing"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        alias_graph.register_alias(*constructing.begin(), inputs[0]);
        return vector<vector<string>>(1, inputs.front()->get_filenames());
    });
    
    ////////////////////////////////////
    // VG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering VG recipes" << endl;
#endif
    
    // meta-recipe for removing variant paths
    auto strip_variant_paths = [&](const vector<const IndexFile*>& inputs,
                                   const IndexingPlan* plan,
                                   const IndexGroup& constructing) {
        
        assert(inputs.size() == 1);
        assert(constructing.size() == 1);
        
        auto chunk_filenames = inputs.at(0)->get_filenames();
        auto output_index = *constructing.begin();
        
        vector<vector<string>> all_outputs(constructing.size());
        vector<string>& output_names = all_outputs[0];
        for (int i = 0; i < chunk_filenames.size(); ++i) {
            // test streams for I/O
            ifstream infile;
            init_in(infile, chunk_filenames[i]);
            
            string output_name = plan->output_filepath(output_index, i, chunk_filenames.size());
            
            ofstream outfile;
            init_out(outfile, output_name);
            
            // FIXME: this crashes as a MutablePathHandleGraph for some reason...
            unique_ptr<MutablePathMutableHandleGraph> graph
                = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(infile);
            
            // gather handles to the alt allele paths
            vector<path_handle_t> alt_paths;
            graph->for_each_path_handle([&](const path_handle_t& path) {
                auto name = graph->get_path_name(path);
                if (!name.empty() && name.substr(0, 5) == "_alt_") {
                    alt_paths.push_back(path);
                }
            });
            
            // delete them
            for (auto path : alt_paths) {
                graph->destroy_path(path);
            }
            
            // and save the graph
            vg::io::save_handle_graph(graph.get(), outfile);
            
            output_names.push_back(output_name);
        }
        
        // return the filename(s)
        return all_outputs;
    };
    
    // strip alt allele paths from a graph that has them
    registry.register_recipe({"VG"}, {"VG w/ Variant Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Stripping allele paths from VG." << endl;
        }
        
        return strip_variant_paths(inputs, plan, constructing);
    });
        
    // meta-recipe for creating a VG from a GFA
    auto construct_from_gfa = [&](const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexGroup& constructing,
                                  nid_t* max_node_id_out) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing VG graph from GFA input." << endl;
        }
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        
        assert(inputs.size() == 1);
        auto output_index = *constructing.begin();
        auto input_filenames = inputs.at(0)->get_filenames();
        if (input_filenames.size() > 1) {
            cerr << "error:[IndexRegistry] Graph construction does not support multiple GFAs at this time." << endl;
            exit(1);
        }
        auto input_filename = input_filenames.front();
        
        string output_name = plan->output_filepath(output_index);
        ofstream outfile;
        init_out(outfile, output_name);
        auto graph = init_mutable_graph();
        
        // make the graph from GFA
        try {
            algorithms::gfa_to_path_handle_graph(input_filename, graph.get(), true,
                                                 IndexingParameters::mut_graph_impl == IndexingParameters::ODGI);
        }
        catch (algorithms::GFAFormatError& e) {
            cerr << "error:[IndexRegistry] Input GFA is not usuable in VG." << endl;
            cerr << e.what() << endl;
            exit(1);
        }
        
        if (max_node_id_out) {
            *max_node_id_out = graph->max_node_id();
        }
        
        // save the file
        vg::io::save_handle_graph(graph.get(), outfile);
        
        // return the filename
        all_outputs[0].push_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"VG"}, {"Reference GFA"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_from_gfa(inputs, plan, constructing, nullptr);
    });
    
    // A meta-recipe to make VG and spliced VG files using the Constructor
    // Expects inputs to be ordered: FASTA, VCF[, GTF/GFF][, Insertion FASTA]
    auto construct_with_constructor = [&](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          const IndexGroup& constructing,
                                          bool alt_paths,
                                          bool has_transcripts) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing";
            if (has_transcripts) {
                cerr << " spliced";
            }
            cerr << " VG graph from FASTA and VCF input." << endl;
        }
        
        
        assert(constructing.size() == 2);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_max_id = *constructing.begin();
        auto output_graph = *constructing.rbegin();
        auto& max_id_names = all_outputs[0];
        auto& graph_names = all_outputs[1];
        
        bool has_ins_fasta = false;
        if (!has_transcripts) {
            assert(inputs.size() == 2 || inputs.size() == 3);
            has_ins_fasta = (inputs.size() == 3);
        }
        else {
            assert(inputs.size() == 3 || inputs.size() == 4);
            has_ins_fasta = (inputs.size() == 4);
        }
        
        
        // unpack the inputs
        vector<string> ref_filenames, vcf_filenames, insertions, transcripts;
        {
            size_t i = 0;
            if (has_transcripts) {
                transcripts = inputs[i++]->get_filenames();
            }
            if (has_ins_fasta) {
                insertions = inputs[i++]->get_filenames();
            }
            ref_filenames = inputs[i++]->get_filenames();
            vcf_filenames = inputs[i++]->get_filenames();
        }
        
        if (ref_filenames.size() != 1 && vcf_filenames.size() != 1 &&
            ref_filenames.size() != vcf_filenames.size()) {
            cerr << "[IndexRegistry]: When constructing graph from multiple FASTAs and multiple VCFs, the FASTAs and VCFs must be matched 1-to-1, but input contains " <<  inputs[0]->get_filenames().size() << " FASTA files and " << inputs[1]->get_filenames().size() << " VCF files." << endl;
            exit(1);
        }
        if (has_transcripts) {
            if ((transcripts.size() != 1 && vcf_filenames.size() != 1 &&
                 transcripts.size() != vcf_filenames.size()) ||
                (transcripts.size() != 1 && ref_filenames.size() != 1 &&
                 transcripts.size() != ref_filenames.size())) {
                cerr << "[IndexRegistry]: When constructing graph from multiple GTF/GFFs and multiple FASTAs or VCFs, the GTF/GFFs and the FASTAs/VCFs must be matched 1-to-1, but input contains " <<  transcripts.size() << " GTF/GFF files, " <<  ref_filenames.size() << " FASTA files, and " << vcf_filenames.size() << " VCF files." << endl;
                exit(1);
            }
        }
        
        // TODO: allow contig renaming through Constructor::add_name_mapping
        // TODO: actually do the chunking based on available memory and in parallel if possible
        
        nid_t max_node_id = 0;
        size_t i = 0, j = 0, k = 0;
        while (i < ref_filenames.size() && j < vcf_filenames.size()) {
            
            // init and configure the constructor
            Constructor constructor;
            constructor.do_svs = true;
            constructor.alt_paths = alt_paths;
            constructor.max_node_size = IndexingParameters::max_node_size;
            constructor.show_progress = IndexingParameters::verbose;
            if (ref_filenames.size() != 1 && vcf_filenames.size() == 1) {
                // we have multiple FASTA but only 1 VCF, so we'll limit the
                // constructor to the contigs of this FASTA for this run
                FastaReference ref;
                ref.open(ref_filenames[i]);
                for (const string& seqname : ref.index->sequenceNames) {
                    constructor.allowed_vcf_names.insert(seqname);
                }
            }
            
            string output_name = plan->output_filepath(output_graph, max(i, j),
                                                       max(ref_filenames.size(), vcf_filenames.size()));
            ofstream outfile;
            init_out(outfile, output_name);
            
            auto graph = init_mutable_graph();
            
            vector<string> fasta(1, ref_filenames[i]);
            vector<string> vcf(1, vcf_filenames[j]);
            
            // do the construction
            constructor.construct_graph(fasta, vcf, insertions, graph.get());
                        
            if (!transcripts.empty()) {
                
                ifstream infile_tx;
                init_in(infile_tx, transcripts[k]);
                
                // are we broadcasting the transcripts from one chunk to many?
                bool broadcasting_txs = transcripts.size() != max(ref_filenames.size(),
                                                                  vcf_filenames.size());
                
                vector<string> path_names;
                if (broadcasting_txs) {
                    // get the path names in case we need to report them later for debug output
                    graph->for_each_path_handle([&](const path_handle_t& path) {
                        path_names.push_back(graph->get_path_name(path));
                    });
                }
                
                // give away ownership of the graph to the Transcriptome
                Transcriptome transcriptome(move(graph));
                transcriptome.error_on_missing_path = !broadcasting_txs;
                transcriptome.use_reference_paths = true;
                transcriptome.feature_type = IndexingParameters::gff_feature_name;
                transcriptome.transcript_tag = IndexingParameters::gff_transcript_tag;
                
                // add the splice edges
                auto dummy = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
                size_t transcripts_added = transcriptome.add_transcript_splice_junctions(infile_tx, dummy);
                                
                if (broadcasting_txs && !path_names.empty() && transcripts_added == 0
                    && transcript_file_nonempty(transcripts[k])) {
                    cerr << "warning:[IndexRegistry] no matching paths from transcript file " << transcripts[k] << " were found in graph chunk containing the following paths:" << endl;
                    for (const string& path_name : path_names) {
                        cerr << "\t" << path_name << endl;
                    }
                }
                
                max_node_id = max(max_node_id, transcriptome.splice_graph().max_node_id());
                
                // save the file
                transcriptome.write_splice_graph(&outfile);
            }
            else {
                max_node_id = max(max_node_id, graph->max_node_id());
                
                // save the file
                vg::io::save_handle_graph(graph.get(), outfile);
            }
            
            graph_names.push_back(output_name);
            
            // awkward incrementation that allows 1x1, 1xN, Nx1, and NxN matching
            // of VCFs, FASTAs, and GTF/GFFs
            if (ref_filenames.size() == 1 && vcf_filenames.size() == 1) {
                break;
            }
            if (ref_filenames.size() > 1) {
                ++i;
            }
            if (vcf_filenames.size() > 1) {
                ++j;
            }
            if (transcripts.size() > 1) {
                ++k;
            }
        }
        
        if (ref_filenames.size() == 1 && vcf_filenames.size() != 1) {
            
            // we will have added components for all of the paths that are also
            // VCF sequences, but we may have FASTA sequences that don't
            // correspond to any VCF sequence (e.g. unlocalized contigs of
            // decoys)
            
            // TODO: how can we keep track of all of the contigs seen in the Constructor
            // so that we know which ones to include in this component?
            
//            string output_name = (plan->prefix(constructing)
//                                  + "." + to_string(inputs[1]->get_filenames().size())
//                                  + "." + plan->suffix(constructing));
            
            // FIXME: punting on this for now
            // FIXME: also need to add a dummy VCF, but how to add it to const index?
            // FIXME: also need to locate these contigs in the GTF/GFF
            // maybe i should add it as a second index, along with this graph...
        }
        
        if (graph_names.size() > 1) {
            // join the ID spaces
            VGset graph_set(graph_names);
            max_node_id = graph_set.merge_id_space();
        }
        
        // save the max node id as a simple text file
        auto max_id_name = plan->output_filepath(output_max_id);
        ofstream max_id_outfile;
        init_out(max_id_outfile, max_id_name);
        max_id_outfile << max_node_id;
        
        max_id_names.push_back(max_id_name);
        
        // return the filename(s)
        return all_outputs;
    };
    
    // the specific instantiations of the meta-recipe above
    registry.register_recipe({"MaxNodeID", "VG"}, {"Insertion Sequence FASTA", "Reference FASTA", "VCF"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, false);
    });
    registry.register_recipe({"MaxNodeID", "VG"}, {"Reference FASTA", "VCF"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, false);
    });
    registry.register_recipe({"MaxNodeID", "VG w/ Variant Paths"}, {"Insertion Sequence FASTA", "Reference FASTA", "VCF w/ Phasing"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, false);
    });
    registry.register_recipe({"MaxNodeID", "VG w/ Variant Paths"}, {"Reference FASTA", "VCF w/ Phasing"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, false);
    });
    
#ifdef debug_index_registry_setup
    cerr << "registering Spliced VG recipes" << endl;
#endif
    
    ////////////////////////////////////
    // Spliced VG Recipes
    ////////////////////////////////////
    
    registry.register_recipe({"Spliced VG"}, {"Spliced VG w/ Variant Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Stripping allele paths from spliced VG." << endl;
        }
        
        return strip_variant_paths(inputs, plan, constructing);
    });
    
    // TODO: spliced vg from GFA input
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG"},
                             {"GTF/GFF", "Insertion Sequence FASTA", "Reference FASTA", "VCF"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG"},
                             {"GTF/GFF", "Reference FASTA", "VCF"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, true);
    });
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Variant Paths"},
                             {"GTF/GFF", "Insertion Sequence FASTA", "Reference FASTA", "VCF w/ Phasing"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Variant Paths"},
                             {"GTF/GFF", "Reference FASTA", "VCF w/ Phasing"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, true);
    });
    
    
    ////////////////////////////////////
    // XG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering XG recipes" << endl;
#endif
    
    registry.register_recipe({"XG"}, {"Reference GFA"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing XG graph from GFA input." << endl;
        }
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_index = *constructing.begin();
        auto gfa_names = inputs.front()->get_filenames();
        if (gfa_names.size() > 1) {
            cerr << "error:[IndexRegistry] Graph construction does not support multiple GFAs at this time." << endl;
            exit(1);
        }
        
        string output_name = plan->output_filepath(output_index);
        ofstream outfile;
        init_out(outfile, output_name);
        
        xg::XG xg_index;
        xg_index.from_gfa(gfa_names.front());
        
        vg::io::save_handle_graph(&xg_index, outfile);
        
        // return the filename
        all_outputs[0].emplace_back(output_name);
        return all_outputs;
    });
    
    auto make_xg_from_graph = [&](const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexGroup& constructing) {
        
        assert(inputs.size() == 1);
        assert(constructing.size() == 1);
        auto output_index = *constructing.begin();
        vector<vector<string>> all_outputs(constructing.size());
        
        auto graph_filenames = inputs.at(0)->get_filenames();
        
        string output_name = plan->output_filepath(output_index);
        ofstream outfile;
        init_out(outfile, output_name);
        
        xg::XG xg_index;
        
        if (graph_filenames.size() == 1) {
            // we do the one-graph conversion directly, which is more efficient than the
            // VGset option
            
            // test streams for I/O
            ifstream infile;
            init_in(infile, graph_filenames.front());
            
            unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            
            xg_index.from_path_handle_graph(*graph);
        }
        else {
            // the inefficient 3-pass, multi-graph construction algorithm
            
            // make a mutable copy of the graph names
            vector<string> graph_files;
            for (const string& graph_file : graph_filenames) {
                // test for I/O while we're at it
                ifstream infile;
                init_in(infile, graph_file);
                
                graph_files.push_back(graph_file);
            }
            
            VGset graph_set(graph_files);
            graph_set.to_xg(xg_index);
        }
        
        vg::io::save_handle_graph(&xg_index, outfile);
        
        // return the filename
        all_outputs[0].emplace_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"XG"}, {"VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing XG graph from VG graph." << endl;
        }
        return make_xg_from_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced XG"}, {"Spliced VG w/ Transcript Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing spliced XG graph from spliced VG graph." << endl;
        }
        return make_xg_from_graph(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // MaxNodeID Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering MaxNodeID recipes" << endl;
#endif
    
    // meta-recipe to write max node id down to a file
    auto write_max_node_id = [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Determining node ID interval." << endl;
        }
        
        // TODO: this is pretty unoptimized in that we have to load the whole graph just
        // to read the max node id
        
        assert(constructing.size() == 1);
        assert(inputs.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_index = *constructing.begin();
        auto graph_files = inputs.at(0)->get_filenames();
        
        // test I/O
        for (const string& graph_file : graph_files) {
            ifstream infile;
            init_in(infile, graph_file);
        }
        string output_name = plan->output_filepath(output_index);
        ofstream outfile;
        init_out(outfile, output_name);
        
        VGset graph_set(graph_files);
        nid_t max_node_id = graph_set.max_node_id();
        
        outfile << max_node_id;
        
        all_outputs[0].push_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"MaxNodeID"}, {"VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return write_max_node_id(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced MaxNodeID"}, {"Spliced VG w/ Transcript Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return write_max_node_id(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // GBWT Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering GBWT recipes" << endl;
#endif
    
    // merge multiple GBWTs if there are multiple, otherwise leave in place
    auto merge_gbwts = [&](const vector<string>& gbwt_names,
                           const IndexingPlan* plan,
                           const IndexName& constructing_name) {
        if (gbwt_names.size() > 1) {
            if (IndexingParameters::verbose) {
                cerr << "[IndexRegistry]: Merging contig GBWTs." << endl;
            }
            // we also need to merge the GBWTs
            
            string merged_gbwt_name = plan->output_filepath(constructing_name);
            ofstream outfile;
            init_out(outfile, merged_gbwt_name);
            
            vector<gbwt::GBWT> gbwt_indexes(gbwt_names.size());
            for (size_t i = 0; i < gbwt_names.size(); ++i) {
                load_gbwt(gbwt_names[i], gbwt_indexes[i], IndexingParameters::verbose);
            }
            gbwt::GBWT merged(gbwt_indexes);
            merged.serialize(outfile);
            return merged_gbwt_name;
        }
        else {
            // note: we don't need to register an alias here because it all happens
            // internally to one index's recipe
            return gbwt_names.front();
        }
    };
    
    // TODO: add Jouni's parallel chunked workflow here
    // meta-recipe to make GBWTs
    auto make_gbwt = [&](const vector<const IndexFile*>& inputs,
                         const IndexingPlan* plan,
                         const IndexGroup& constructing) {
        
        assert(inputs.size() == 2);
        
        vector<string> graph_filenames, vcf_filenames;
        // the two recipes that use this meta-recipe have the alphabetical order
        // of their inputs swapped relative to each other
        if (*constructing.begin() == "GBWT") {
            graph_filenames = inputs[1]->get_filenames();
            vcf_filenames = inputs[0]->get_filenames();
        }
        else {
            graph_filenames = inputs[0]->get_filenames();
            vcf_filenames = inputs[1]->get_filenames();
        }
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        
        auto output_index = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        if ((graph_filenames.size() != 1 && graph_filenames.size() != vcf_filenames.size()) ||
            (vcf_filenames.size() != 1 && graph_filenames.size() != vcf_filenames.size())) {
            cerr << "[IndexRegistry]: When constructing GBWT from multiple graphs and multiple VCFs, the graphs and VCFs must be matched 1-to-1, but input contains " <<  graph_filenames.size() << " graphs and " << vcf_filenames.size() << " VCF files." << endl;
            exit(1);
        }
        
        if (IndexingParameters::verbose) {
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        
        vector<string> gbwt_names;
        
        // construct a GBWT from the i-th VCF
        auto gbwt_job = [&](size_t i, const PathHandleGraph& graph) {
            string gbwt_name;
            if (vcf_filenames.size() != 1) {
                // multiple components, so make a temp file that we will merge later
                gbwt_name = temp_file::create();
            }
            else {
                // one component, so we will actually save the output
                gbwt_name = plan->output_filepath(output_index);
            }
            
            // make this critical so we don't end up with a race on the verbosity
            unique_ptr<HaplotypeIndexer> haplotype_indexer;
#pragma omp critical
            {
                haplotype_indexer = unique_ptr<HaplotypeIndexer>(new HaplotypeIndexer());
                // HaplotypeIndexer resets this in its constructor
                if (IndexingParameters::verbose) {
                    gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
                }
                else {
                    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                }
            }
            haplotype_indexer->show_progress = IndexingParameters::verbose;
            
            vector<string> parse_files = haplotype_indexer->parse_vcf(vcf_filenames[i],
                                                                      graph);
            
            unique_ptr<gbwt::DynamicGBWT> gbwt_index = haplotype_indexer->build_gbwt(parse_files);
            
            vg::io::VPKG::save(*gbwt_index, gbwt_name);
            
            gbwt_names.push_back(gbwt_name);
        };
        
        if (graph_filenames.size() == 1) {
            // we only have one graph, so we can save time by loading it only one time
            // test streams for I/O
            ifstream infile;
            init_in(infile, graph_filenames.front());
            
            unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            
            for (size_t i = 0; i < vcf_filenames.size(); ++i) {
                gbwt_job(i, *graph);
            }
        }
        else {
            // FIXME: it doesn't seem possible to do many component graphs with 1 VCF
            for (size_t i = 0; i < graph_filenames.size(); ++i) {
                ifstream infile;
                init_in(infile, graph_filenames[i]);
                unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
                gbwt_job(i, *graph);
            }
        }
        
        // merge GBWTs if necessary
        output_names.push_back(merge_gbwts(gbwt_names, plan, output_index));
        // return filename
        return all_outputs;
    };
    
    registry.register_recipe({"GBWT"}, {"VCF w/ Phasing", "VG w/ Variant Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing GBWT from VG graph and phased VCF input." << endl;
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        return make_gbwt(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GBWT"}, {"Spliced VG w/ Variant Paths", "VCF w/ Phasing"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing GBWT from spliced VG graph and phased VCF input." << endl;
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        return make_gbwt(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Haplotype-Transcript GBWT", "Spliced VG w/ Transcript Paths", "Unjoined Transcript Origin Table", },
                             {"GTF/GFF", "Spliced GBWT", "Spliced VG", },
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing haplotype-transcript GBWT and finishing spliced VG." << endl;
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        assert(constructing.size() == 3);
        vector<vector<string>> all_outputs(constructing.size());
        IndexName output_haplo_tx, output_tx_table, output_tx_graph;
        {
            int i = 0;
            for (auto output_index : constructing) {
                if (i == 0) {
                    output_haplo_tx = output_index;
                }
                else if (i == 1) {
                    output_tx_graph = output_index;
                }
                else {
                    output_tx_table = output_index;
                }
                ++i;
            }
        }
        auto& haplo_tx_gbwt_names = all_outputs[0];
        auto& tx_graph_names = all_outputs[1];
        auto& tx_table_names = all_outputs[2];
        
        assert(inputs.size() == 3);
        auto tx_filenames = inputs[0]->get_filenames();
        auto gbwt_filenames = inputs[1]->get_filenames();
        auto graph_filenames = inputs[2]->get_filenames();
        
        assert(gbwt_filenames.size() == 1);
        auto gbwt_filename = gbwt_filenames.front();
        
        unique_ptr<gbwt::GBWT> haplotype_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_filename);
        
        // TODO: i can't find where in the building code you actually ensure this...
        assert(haplotype_index->bidirectional());
        
        // TODO: repetitive with original spliced VG construction...
        
        vector<string> gbwt_chunk_names;
        for (size_t i = 0, j = 0; i < graph_filenames.size(); ++i) {
            string tx_graph_name = plan->output_filepath(output_tx_graph, i, graph_filenames.size());
            ofstream tx_graph_outfile;
            init_out(tx_graph_outfile, tx_graph_name);
            
            string gbwt_name;
            if (graph_filenames.size() != 1) {
                // multiple components, so make a temp file that we will merge later
                gbwt_name = temp_file::create();
            }
            else {
                // one component, so we will actually save the output
                gbwt_name = plan->output_filepath(output_haplo_tx, i, graph_filenames.size());
            }
            
            string info_table_name = plan->output_filepath(output_tx_table, i, graph_filenames.size());
            ofstream info_outfile;
            init_out(info_outfile, info_table_name);
            
            ifstream infile_graph, infile_tx;
            init_in(infile_graph, graph_filenames[i]);
            init_in(infile_tx, tx_filenames[j]);
            
            // are we using 1 transcript file for multiple graphs?
            bool broadcasting_txs = (graph_filenames.size() != tx_filenames.size());
                        
            unique_ptr<MutablePathDeletableHandleGraph> graph
                = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(infile_graph);
            
            vector<string> path_names;
            if (broadcasting_txs) {
                // get the path names in case we need to report them later for debug output
                graph->for_each_path_handle([&](const path_handle_t& path) {
                    path_names.push_back(graph->get_path_name(path));
                });
            }
            
            Transcriptome transcriptome(move(graph));
            transcriptome.error_on_missing_path = !broadcasting_txs;
            transcriptome.feature_type = IndexingParameters::gff_feature_name;
            transcriptome.transcript_tag = IndexingParameters::gff_transcript_tag;
            
            // load up the transcripts and add edges on the reference path
            size_t transcripts_added = transcriptome.add_transcript_splice_junctions(infile_tx, haplotype_index);
            
            if (broadcasting_txs && !path_names.empty() && transcripts_added == 0
                && transcript_file_nonempty(tx_filenames[j])) {
                cerr << "warning:[IndexRegistry] no matching paths from transcript file " << tx_filenames[j] << " were found in graph chunk containing the following paths:" << endl;
                for (const string& path_name : path_names) {
                    cerr << "\t" << path_name << endl;
                }
            }
            
            // go back to the beginning of the transcripts
            infile_tx.clear();
            infile_tx.seekg(0);
            
            // add egdes on other haplotypes
            size_t num_transcripts_projected = transcriptome.add_transcripts(infile_tx, *haplotype_index);
            
            // init the haplotype transcript GBWT
            size_t node_width = gbwt::bit_length(gbwt::Node::encode(transcriptome.splice_graph().max_node_id(), true));
            gbwt::GBWTBuilder gbwt_builder(node_width,
                                           IndexingParameters::gbwt_insert_batch_size,
                                           IndexingParameters::gbwt_sampling_interval);
            // actually build it
            transcriptome.add_transcripts_to_gbwt(&gbwt_builder,
                                                  false, // don't output transcript paths as a FASTA
                                                  IndexingParameters::bidirectional_haplo_tx_gbwt);
            
            // save the haplotype transcript GBWT
            gbwt_builder.finish();
            vg::io::VPKG::save(gbwt_builder.index, gbwt_name);
            
            // write transcript origin info table
            transcriptome.write_info(&info_outfile, *haplotype_index, false);
            
            // save the graph with the transcript paths added
            transcriptome.write_splice_graph(&tx_graph_outfile);
            
            tx_graph_names.push_back(tx_graph_name);
            tx_table_names.push_back(info_table_name);
            gbwt_chunk_names.push_back(gbwt_name);
            
            if (tx_filenames.size() > 1) {
                ++j;
            }
        }
        
        haplo_tx_gbwt_names.push_back(merge_gbwts(gbwt_chunk_names, plan, output_haplo_tx));
        
        return all_outputs;
    });
    
    ////////////////////////////////////
    // Info Table Recipes
    ////////////////////////////////////
    
    registry.register_recipe({"Transcript Origin Table"}, {"Unjoined Transcript Origin Table"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Joining transcript origin table." << endl;
        }
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_index = *constructing.begin();
        
        assert(inputs.size() == 1);
        auto input_table_names = inputs[0]->get_filenames();
        
        if (input_table_names.size() == 1) {
            alias_graph.register_alias(output_index, inputs[0]);
            all_outputs[0] = input_table_names;
            return all_outputs;
        }
        string output_name = plan->output_filepath(output_index);
        
        ofstream outfile;
        init_out(outfile, output_name);
        // join the tables into one
        for (size_t i = 0; i < inputs[0]->get_filenames().size(); ++i) {
            ifstream infile(inputs[0]->get_filenames()[i]);
            string line;
            if (i != 0) {
                // skip the header
                infile.ignore(numeric_limits<streamsize>::max(), '\n');
            }
            else {
                getline(infile, line);
                outfile << line;
            }
            while (infile.good()) {
                getline(infile, line);
                outfile << endl << line;
            }
        }
        all_outputs[0].push_back(output_name);
        return all_outputs;
    });
    
    ////////////////////////////////////
    // Pruned VG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering pruning recipes" << endl;
#endif
    
    // meta-recipe for pruning with/without GBWT
    auto prune_graph = [&](const vector<const IndexFile*>& inputs,
                           const IndexingPlan* plan,
                           const IndexGroup& constructing) {
        
        // we only want to focus on two specific recipes
        assert(inputs.size() == 2 || inputs.size() == 3);
        bool using_haplotypes = inputs.size() == 3;
        
        vector<string> graph_names;
        string gbwt_name, max_node_name;
        {
            size_t i = 0;
            if (using_haplotypes) {
                auto gbwt_names = inputs[i++]->get_filenames();
                assert(gbwt_names.size() == 1);
                gbwt_name = gbwt_names.front();
            }
            auto max_node_id_names = inputs[i++]->get_filenames();
            assert(max_node_id_names.size() == 1);
            max_node_name = max_node_id_names.front();
            graph_names = inputs[i++]->get_filenames();
        }
        
        if (using_haplotypes) {
            assert(constructing.size() == 2);
        }
        else {
            assert(constructing.size() == 1);
        }
        
        vector<vector<string>> all_outputs(constructing.size());
        auto& pruned_graph_names = all_outputs[0];
        auto output_pruned_graph = *constructing.begin();
        
        // test streams for I/O
        ifstream infile_gbwt, infile_max_id;
        init_in(infile_max_id, max_node_name);
        unique_ptr<gbwt::GBWT> gbwt_index;
        if (using_haplotypes) {
            init_in(infile_gbwt, gbwt_name);
            gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(infile_gbwt);
        }
        
        // read the max node ID (across all chunks)
        nid_t max_node_id;
        infile_max_id >> max_node_id;
        
        string mapping_name;
        if (using_haplotypes) {
            
            auto output_mapping = *constructing.rbegin();
            
            gcsa::NodeMapping mapping(max_node_id + 1);
            
            mapping_name = plan->output_filepath(output_mapping);
            
            ofstream mapping_file;
            init_out(mapping_file, mapping_name);
            mapping.serialize(mapping_file);
            
            all_outputs[1].push_back(mapping_name);
        }
        
        // TODO: is it possible to do this in parallel? i'm worried about races
        // on the node mapping, and it seems like it definitely wouldn't work
        // for haplotype unfolding, which needs to modify the NodeMapping
        for (size_t i = 0; i < graph_names.size(); ++i) {
            
            ifstream infile_vg;
            init_in(infile_vg, graph_names[i]);
            
            string vg_output_name = plan->output_filepath(output_pruned_graph, i, graph_names.size());
            
            ofstream outfile_vg;
            init_out(outfile_vg, vg_output_name);
            
            unique_ptr<MutablePathDeletableHandleGraph> graph
                = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(infile_vg);
            
            // destroy all paths, which might be made inconsistent
            vector<path_handle_t> paths;
            paths.reserve(graph->get_path_count());
            graph->for_each_path_handle([&](const path_handle_t& path) {
                paths.push_back(path);
            });
            for (auto path : paths) {
                graph->destroy_path(path);
            }
            
            // prune the graph based on topology
            size_t removed_high_degree, removed_complex, removed_subgraph;
            if (IndexingParameters::pruning_max_node_degree != 0) {
                removed_high_degree = algorithms::remove_high_degree_nodes(*graph, IndexingParameters::pruning_max_node_degree);
            }
            removed_complex = algorithms::prune_complex_with_head_tail(*graph, IndexingParameters::pruning_walk_length,
                                                               IndexingParameters::pruning_max_edge_count);
            removed_subgraph = algorithms::prune_short_subgraphs(*graph, IndexingParameters::pruning_min_component_size);
            
            if ((removed_high_degree != 0 || removed_complex != 0 || removed_subgraph != 0)
                && (!paths.empty() || using_haplotypes)) {
                // we've removed from this graph but there are paths/threads we could use
                // to restore the graph
                
                // TODO: in a single component graph, it would be more efficient to load
                // an XG rather than the mutable graph for this step
                ifstream infile_unpruned_vg;
                init_in(infile_unpruned_vg, graph_names[i]);
                
                unique_ptr<PathHandleGraph> unpruned_graph
                    = vg::io::VPKG::load_one<PathHandleGraph>(infile_unpruned_vg);
                
                if (!using_haplotypes) {
                    // we can bring back edges on embedded paths
                    
                    // Make an empty GBWT index to pass along
                    gbwt::GBWT empty_gbwt;
                    PhaseUnfolder unfolder(*unpruned_graph, empty_gbwt, max_node_id + 1);
                    unfolder.restore_paths(*graph, IndexingParameters::verbose);
                }
                else {
                    // we can expand out complex regions using haplotypes as well as paths
                    // TODO: it's a bit inelegant that i keep overwriting the mapping...
                    
                    PhaseUnfolder unfolder(*unpruned_graph, *gbwt_index, max_node_id + 1);
                    unfolder.read_mapping(mapping_name);
                    unfolder.unfold(*graph, IndexingParameters::verbose);
                    unfolder.write_mapping(mapping_name);
                }
            }
            
            vg::io::save_handle_graph(graph.get(), outfile_vg);
            
            pruned_graph_names.push_back(vg_output_name);
        }
        
        return all_outputs;
    };
    
    registry.register_recipe({"Pruned VG"}, {"MaxNodeID", "VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Haplotype-Pruned VG", "Unfolded NodeMapping"}, {"GBWT", "MaxNodeID", "VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing with GBWT unfolding." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Pruned Spliced VG"}, {"Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of spliced VG to prepare for GCSA indexing." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    // TODO: would it be better to use the Haplotype-Transcript GBWT, or maybe to join them?
    // the splice edges will be covered by the transcript paths, so it won't be too bad
    registry.register_recipe({"Haplotype-Pruned Spliced VG", "Unfolded Spliced NodeMapping"},
                             {"Spliced GBWT", "Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of spliced VG to prepare for GCSA indexing with GBWT unfolding." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // GCSA + LCP Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering GCSA recipes" << endl;
#endif
    
    // meta-recipe for GCSA indexing with or without unfolded input
    auto construct_gcsa = [&](const vector<const IndexFile*>& inputs,
                              const IndexingPlan* plan,
                              const IndexGroup& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing GCSA/LCP indexes." << endl;
        }
        
        assert(inputs.size() == 1 || inputs.size() == 2);
        bool unfolded = inputs.size() == 2;
        auto graph_filenames = inputs[0]->get_filenames();
        string mapping_filename;
        if (unfolded) {
            auto mapping_filenames = inputs[1]->get_filenames();
            assert(mapping_filenames.size() == 1);
            mapping_filename = mapping_filenames.front();
        }
        
        assert(constructing.size() == 2);
        vector<vector<string>> all_outputs(constructing.size());
        auto& gcsa_names = all_outputs[0];
        auto& lcp_names = all_outputs[1];
        
        auto output_gcsa = *constructing.begin();
        auto output_lcp = *constructing.rbegin();
        
        // test streams for I/O
        ifstream infile_mapping;
        if (unfolded) {
            init_in(infile_mapping, mapping_filename);
        }
        string gcsa_output_name = plan->output_filepath(output_gcsa);
        string lcp_output_name = plan->output_filepath(output_lcp);
        ofstream outfile_gcsa;
        ofstream outfile_lcp;
        init_out(outfile_gcsa, gcsa_output_name);
        init_out(outfile_lcp, lcp_output_name);
        
        // configure GCSA to use scratch in the general temp directory
        gcsa::TempFile::setDirectory(temp_file::get_dir());
        if (IndexingParameters::verbose) {
            gcsa::Verbosity::set(gcsa::Verbosity::BASIC);
        }
        else {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }
        auto params = gcsa::ConstructionParameters();
        params.doubling_steps = IndexingParameters::gcsa_doubling_steps;
                
#ifdef debug_index_registry_recipes
        cerr << "enumerating k-mers for input pruned graphs:" << endl;
        for (auto& name : graph_filenames) {
            cerr << "\t" << name << endl;
        }
#endif
        
        VGset graph_set(graph_filenames);
        size_t kmer_bytes = params.getLimitBytes();
        vector<string> dbg_names = graph_set.write_gcsa_kmers_binary(IndexingParameters::gcsa_initial_kmer_length,
                                                                     kmer_bytes);
        
#ifdef debug_index_registry_recipes
        cerr << "making GCSA2" << endl;
#endif
        
        // construct the indexes (giving empty mapping name is sufficient to make
        // indexing skip the unfolded code path)
        gcsa::InputGraph input_graph(dbg_names, true, gcsa::Alphabet(),
                                     mapping_filename);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        
        // clean up the k-mer files
        for (auto dbg_name : dbg_names) {
            temp_file::remove(dbg_name);
        }
        
#ifdef debug_index_registry_recipes
        cerr << "saving GCSA/LCP pair" << endl;
#endif
        
        vg::io::VPKG::save(gcsa_index, gcsa_output_name);
        vg::io::VPKG::save(lcp_array, lcp_output_name);
        
        gcsa_names.push_back(gcsa_output_name);
        lcp_names.push_back(lcp_output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"GCSA", "LCP"}, {"Haplotype-Pruned VG", "Unfolded NodeMapping"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"GCSA", "LCP"}, {"Pruned VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GCSA", "Spliced LCP"}, {"Haplotype-Pruned Spliced VG", "Unfolded Spliced NodeMapping"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GCSA", "Spliced LCP"}, {"Pruned Spliced VG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // Snarls Recipes
    ////////////////////////////////////
    
    // meta-recipe to find snarls
    auto find_snarls = [&](const vector<const IndexFile*>& inputs,
                           const IndexingPlan* plan,
                           const IndexGroup& constructing) {
        
        assert(inputs.size() == 1);
        auto xg_filenames = inputs[0]->get_filenames();
        assert(xg_filenames.size() == 1);
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_snarls = *constructing.begin();
        auto& snarl_names = all_outputs[0];
        
        ifstream infile;
        init_in(infile, xg_filename);
        
        string output_name = plan->output_filepath(output_snarls);
        ofstream outfile;
        init_out(outfile, output_name);
        
        // load graph
        unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
        
        // find snarls
        unique_ptr<SnarlFinder> snarl_finder = unique_ptr<SnarlFinder>(new IntegratedSnarlFinder(*graph));
        SnarlManager snarl_manager = snarl_finder->find_snarls_parallel();
        
        // traverse snarl tree and write them out
        vector<Snarl> buffer;
        for (auto root : snarl_manager.top_level_snarls()) {
            vector<const Snarl*> stack(1, root);
            while (!stack.empty()) {
                const Snarl* snarl = stack.back();
                stack.pop_back();
                
                buffer.push_back(*snarl);
                vg::io::write_buffered(outfile, buffer, 1024);
                
                for (const Snarl* child_snarl : snarl_manager.children_of(snarl)) {
                    stack.push_back(child_snarl);
                }
            }
        }
        // flush
        vg::io::write_buffered(outfile, buffer, 0);
        
        snarl_names.push_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"Snarls"}, {"XG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Finding snarls in graph." << endl;
        }
        return find_snarls(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced Snarls"}, {"Spliced XG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Finding snarls in spliced graph." << endl;
        }
        return find_snarls(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // Distance Index Recipes
    ////////////////////////////////////
    
    // meta-recipe to make distance index
    auto make_distance_index = [&](const vector<const IndexFile*>& inputs,
                                   const IndexingPlan* plan,
                                   const IndexGroup& constructing) {
        
        assert(inputs.size() == 2);
        auto snarls_filenames = inputs[0]->get_filenames();
        auto xg_filenames = inputs[1]->get_filenames();
        assert(xg_filenames.size() == 1);
        assert(snarls_filenames.size() == 1);
        auto snarls_filename = snarls_filenames.front();
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto dist_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        ifstream infile_graph, infile_snarls;
        init_in(infile_graph, xg_filename);
        init_in(infile_snarls, snarls_filename);
        string output_name = plan->output_filepath(dist_output);
        ofstream outfile;
        init_out(outfile, output_name);
        
        unique_ptr<HandleGraph> graph = vg::io::VPKG::load_one<HandleGraph>(infile_graph);
        unique_ptr<SnarlManager> snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(infile_snarls));
        
        MinimumDistanceIndex distance_index(graph.get(), snarl_manager.get());
        
        vg::io::VPKG::save(distance_index, output_name);
        
        output_names.push_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"Distance Index"}, {"Snarls", "XG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Making distance index." << endl;
        }
        return make_distance_index(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced Distance Index"}, {"Spliced Snarls", "Spliced XG"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Making distance index for a spliced graph." << endl;
        }
        return make_distance_index(inputs, plan, constructing);
    });
    
    return registry;
}

vector<IndexName> VGIndexes::get_default_map_indexes() {
    vector<IndexName> indexes{
        "XG",
        "GCSA",
        "LCP"
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_mpmap_indexes() {
    vector<IndexName> indexes{
        "Spliced XG",
        "Spliced Distance Index",
        "Spliced GCSA",
        "Spliced LCP",
        "Haplotype-Transcript GBWT",
        "Transcript Origin Table"
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_giraffe_indexes() {
    vector<IndexName> indexes{
        "GBWT",
        "GBWTGraph",
        "Distance",
        "Minimizer"
    };
    return indexes;
}

bool IndexingPlan::is_intermediate(const IndexName& identifier) const {
    if (registry->get_index(identifier)->was_provided_directly()) {
        // It's not an intermediate if it is input
        return false;
    }
    // Or if it is directly requested
    return !targets.count(identifier);
}
    
string IndexingPlan::output_filepath(const IndexName& identifier) const {
    return output_filepath(identifier, 0, 1);
}

string IndexingPlan::output_filepath(const IndexName& identifier, size_t chunk, size_t num_chunks) const {
    
    string filepath;
    if (registry->keep_intermediates || !is_intermediate(identifier)) {
        // we're saving this file, put it at the output prefix
        filepath = registry->output_prefix;
    } else {
        // we're not saving this file, make it temporary
        filepath = registry->get_work_dir() + "/" + sha1sum(identifier);
    }
    if (num_chunks > 1) {
        // we add digits to make the suffix unique for this chunk (the setup disallows suffixes
        // that start with digits)
        filepath += "." + to_string(chunk);
    }
    filepath += "." + registry->get_index(identifier)->get_suffix();
    return filepath;
}
 
const vector<RecipeName>& IndexingPlan::get_steps() const {
    return steps;
}

IndexRegistry::~IndexRegistry() {
    if (!work_dir.empty()) {
        // Clean up our work directory with its temporary indexes.
        temp_file::remove(work_dir);
        work_dir.clear();
    }
}

IndexRegistry::IndexRegistry(IndexRegistry&& other) :
    index_registry(std::move(other.index_registry)),
    recipe_registry(std::move(other.recipe_registry)),
    registered_suffixes(std::move(other.registered_suffixes)),
    work_dir(std::move(other.work_dir)),
    output_prefix(std::move(other.output_prefix)),
    keep_intermediates(std::move(other.keep_intermediates)) {
    
    // Make sure other doesn't delete our work dir when it goes away
    other.work_dir.clear();
}

IndexRegistry& IndexRegistry::operator=(IndexRegistry&& other) {
    index_registry = std::move(other.index_registry);
    recipe_registry = std::move(other.recipe_registry);
    registered_suffixes = std::move(other.registered_suffixes);
    work_dir = std::move(other.work_dir);
    output_prefix = std::move(other.output_prefix);
    keep_intermediates = std::move(other.keep_intermediates);
    
    // Make sure other doesn't delete our work dir when it goes away
    other.work_dir.clear();
    
    return *this;
}

void IndexRegistry::set_prefix(const string& prefix) {
    this->output_prefix = prefix;
}

void IndexRegistry::set_intermediate_file_keeping(bool keep_intermediates) {
    this->keep_intermediates = keep_intermediates;
}

void IndexRegistry::make_indexes(const vector<IndexName>& identifiers) {
    
    // figure out the best plan to make the objectives from the inputs
    IndexGroup identifier_group(identifiers.begin(), identifiers.end());
    auto plan = make_plan(identifier_group);
    
    // to store the results of indexes we create
    map<IndexGroup, vector<vector<string>>> indexing_results;
    
    // to keep track of which indexes are aliases of others
    AliasGraph alias_graph;
    
    // execute the plan
    for (const auto& step : plan.steps) {
        indexing_results[step.first] = execute_recipe(step, &plan, alias_graph);
        assert(indexing_results[step.first].size() == step.first.size());
        auto it = step.first.begin();
        for (const auto& results : indexing_results[step.first]) {
            auto index = get_index(*it);
            if (!index->is_finished()) {
                // the index wasn't already provided directly
                index->assign_constructed(results);
            }
            ++it;
        }
    }
#ifdef debug_index_registry
    cerr << "finished executing recipes, resolving aliases" << endl;
#endif
    
    auto aliases = alias_graph.non_intermediate_aliases(&plan, keep_intermediates);
    for (auto& alias_record : aliases) {
        IndexName aliasee;
        vector<IndexName> aliasors;
        tie(aliasee, aliasors) = alias_record;
        
#ifdef debug_index_registry
        cerr << "index " << aliasee << " is aliased by:" << endl;
        for (const auto& aliasor : aliasors) {
            cerr << "\t" << aliasor << endl;
        }
#endif
        
        // if the index it itself non-intermediate, it will be in the list of aliases.
        // otherwise, we can alias one index by moving instead of copying
        auto f = find(aliasors.begin(), aliasors.end(), aliasee);
        bool can_move = f == aliasors.end();
        if (!can_move) {
            // just remove the "alias" so we don't need to deal with it
            std::swap(*f, aliasors.back());
            aliasors.pop_back();
        }
        
        const auto& aliasee_filenames = get_index(aliasee)->get_filenames();
        
        // copy aliases for any that we need to
        for (size_t i = can_move; i < aliasors.size(); ++i) {
            for (size_t j = 0; j < aliasee_filenames.size(); ++j) {
                
                auto copy_filename = plan.output_filepath(aliasors[i], j, aliasee_filenames.size());
                ifstream aliasee(aliasee_filenames[j], std::ios::binary);
                ofstream aliasor(copy_filename, std::ios::binary);
                
                aliasor << aliasee.rdbuf();
            }
        }
        // if we can move the aliasee (i.e. it is intermediate), then make
        // one index by moving instead of copying
        if (can_move) {
            for (size_t j = 0; j < aliasee_filenames.size(); ++j) {
                auto move_filename = plan.output_filepath(aliasors[0], j, aliasee_filenames.size());
                rename(aliasee_filenames[j].c_str(), move_filename.c_str());
            }
        }
    }
    
    // prepare the index registry to go again, if necessary
    reset();
}

void IndexRegistry::register_index(const IndexName& identifier, const string& suffix) {
    // Add this index to the registry
    if (identifier.empty()) {
        cerr << "error:[IndexRegistry] indexes must have a non-empty identifier" << endl;
        exit(1);
    }
    if (suffix.empty()) {
        cerr << "error:[IndexRegistry] indexes must have a non-empty suffix" << endl;
        exit(1);
    }
    if (isdigit(suffix.front())) {
        // this ensures that we can add numbers to the suffix to create a unique suffix
        // for chunked workflows
        cerr << "error:[IndexRegistry] suffixes cannot start with a digit" << endl;
        exit(1);
    }
    if (index_registry.count(identifier)) {
        cerr << "error:[IndexRegistry] index registry contains a duplicated identifier: " << identifier << endl;
        exit(1);
    }
    if (registered_suffixes.count(suffix)) {
        cerr << "error:[IndexRegistry] index registry contains a duplicated suffix: " << suffix << endl;
        exit(1);
    }
    index_registry[identifier] = unique_ptr<IndexFile>(new IndexFile(identifier, suffix));
    registered_suffixes.insert(suffix);
}


void IndexRegistry::provide(const IndexName& identifier, const string& filename) {
    if (!index_registry.count(identifier)) {
        cerr << "error:[IndexRegistry] cannot provide unregistered index: " << identifier << endl;
        exit(1);
    }
    provide(identifier, vector<string>(1, filename));
}

void IndexRegistry::provide(const IndexName& identifier, const vector<string>& filenames) {
    get_index(identifier)->provide(filenames);
}

vector<IndexName> IndexRegistry::completed_indexes() const {
    vector<IndexName> indexes;
    for (const auto& index : index_registry) {
        if (index.second->is_finished()) {
            indexes.push_back(index.first);
        }
    }
    return indexes;
}

RecipeName IndexRegistry::register_recipe(const vector<IndexName>& identifiers,
                                          const vector<IndexName>& input_identifiers,
                                          const RecipeFunc& exec) {
    
    for (const IndexName& identifier : identifiers) {
        if (!index_registry.count(identifier)) {
            cerr << "error:[IndexRegistry] cannot register recipe for unregistered index " << identifier << endl;
            exit(1);
        }
    }
    
    // test that the input identifiers are in alphabetical order
    // this is an easy-to-troubleshoot check that lets us use IndexGroup's internally and
    // still provide the vector<IndexFile*> in the same order as the input identifiers to
    // the RecipeFunc and in the recipe declaration
    IndexGroup input_group(input_identifiers.begin(), input_identifiers.end());
    IndexGroup output_group(identifiers.begin(), identifiers.end());
    {
        if (input_group.size() != input_identifiers.size()) {
            cerr << "error:[IndexRegistry] recipe has duplicate inputs" << endl;
            exit(1);
        }
        size_t i = 0;
        for (const auto& sorted_identifier : input_group) {
            if (sorted_identifier != input_identifiers[i]) {
                cerr << "error:[IndexRegistry] recipe has inputs that are not provided in alphabetical order" << endl;
                exit(1);
            }
            ++i;
        }
        if (output_group.size() != identifiers.size()) {
            cerr << "error:[IndexRegistry] recipe has duplicate outputs" << endl;
            exit(1);
        }
        i = 0;
        for (const auto& sorted_identifier : output_group) {
            if (sorted_identifier != identifiers[i]) {
                cerr << "error:[IndexRegistry] recipe has outputs that are not provided in alphabetical order" << endl;
                exit(1);
            }
            ++i;
        }
    }
    
    vector<const IndexFile*> inputs;
    for (const auto& input_identifier : input_identifiers) {
        if (!index_registry.count(input_identifier)) {
            cerr << "error:[IndexRegistry] cannot register recipe from unregistered index " << input_identifier << endl;
            exit(1);
        }
        inputs.push_back(get_index(input_identifier));
    }
#ifdef debug_index_registry_setup
    cerr << "registering recipe for " << to_string(output_group) << endl;
    cerr << "inputs:" << endl;
    for (const auto& input : inputs) {
        cerr << "\t" << input->get_identifier() << endl;
    }
#endif
    
    recipe_registry[output_group].emplace_back(inputs, exec);
    RecipeName name(output_group, recipe_registry[output_group].size() - 1);
    
    if (output_group.size() > 1) {
        // add unboxing recipes at the same priority level as the full recipe
        auto it = output_group.begin();
        for (size_t i = 0; i < identifiers.size(); ++i) {
            register_recipe({*it}, identifiers,
                            [=](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
                return vector<vector<string>>(1, inputs[i]->get_filenames());
            });
            ++it;
        }
    }
    
    return name;
}

//void IndexRegistry::register_joint_recipe(const vector<IndexName>& identifiers,
//                                          const vector<IndexName>& input_identifiers,
//                                          const JointRecipeFunc& exec) {
//    // We're going to generate a bunch of single-index recipes where the first
//    // one to run calls the joint recipe, and other ones to run just return
//    // their slice of the joint recipe's return value.
//
//    // We need all the joint recipe names, one for each identifier we generate
//    vector<RecipeName> names;
//
//    // We need a place to hold the return values we can carry around by value.
//    shared_ptr<vector<vector<string>>> results(std::make_shared<vector<vector<string>>>());
//
//    for (size_t i = 0; i < identifiers.size(); i++) {
//        IndexName being_generated = identifiers[i];
//
//        // Create a recipe that invokes the joint recipe.
//        RecipeFunc stub = [i, results, being_generated, &exec](const vector<const IndexFile*>& inputs, const IndexingPlan* plan, const IndexName& constructing) -> vector<string> {
//            if (results->empty()) {
//                // Invoke the actual logic, passing along the plan, and fill in results
//                *results = exec(inputs, plan);
//                // TODO: handle parallel invocations?
//            }
//
//            // Get our slice of the result file list.
//            return results->at(i);
//        };
//
//        names.push_back(register_recipe(being_generated, input_identifiers, stub));
//    }
//
//    // Remember that these are a joint recipe.
//    simplifications.emplace_back(input_identifiers, names);
//}

IndexFile* IndexRegistry::get_index(const IndexName& identifier) {
    return index_registry.at(identifier).get();
}

const IndexFile* IndexRegistry::get_index(const IndexName& identifier) const {
    return index_registry.at(identifier).get();
}

bool IndexRegistry::all_finished(const vector<const IndexFile*>& inputs) const {
    IndexGroup group;
    for (auto input : inputs) {
        group.insert(input->get_identifier());
    }
    return all_finished(group);
}

bool IndexRegistry::all_finished(const IndexGroup& indexes) const {
    bool finished = true;
    for (const auto& index_name : indexes) {
        if (!get_index(index_name)->is_finished()) {
            finished = false;
            break;
        }
    }
    return finished;
}

void IndexRegistry::reset() {
    for (pair<const IndexName, unique_ptr<IndexFile>>& index : index_registry) {
        index.second->reset();
    }
}

string IndexRegistry::get_work_dir() {
    if (work_dir.empty()) {
        // Ensure the directory exists
        work_dir = temp_file::create_directory();
    }
    return work_dir;
}

vector<IndexGroup> IndexRegistry::dependency_order() const {
    
#ifdef debug_index_registry
    cerr << "finding topological order in dependency graph" << endl;
#endif
    
    // assign each index file an index in a vector (arbitrarily) and build the dependency graph
    map<IndexGroup, size_t> graph_idx;
    vector<IndexGroup> graph_label;
    vector<vector<size_t>> dependency_graph;
    // add nodes for the index groups
    for (const auto& recipe_record : recipe_registry) {
        if (!graph_idx.count(recipe_record.first)) {
            graph_idx[recipe_record.first] = graph_label.size();
            graph_label.push_back(recipe_record.first);
            dependency_graph.emplace_back();
        }
        for (auto output : recipe_record.first) {
            IndexGroup singleton_output{output};
            if (!graph_idx.count(singleton_output)) {
                graph_idx[singleton_output] = graph_label.size();
                graph_label.push_back(singleton_output);
                dependency_graph.emplace_back();
            }
        }
        for (const auto& recipe : recipe_record.second) {
            for (auto input : recipe.input_group()) {
                IndexGroup singleton_input{input};
                if (!graph_idx.count(singleton_input)) {
                    graph_idx[singleton_input] = graph_label.size();
                    graph_label.push_back(singleton_input);
                    dependency_graph.emplace_back();
                }
            }
        }
    }
    
    // add nodes for the recipes and recipe edges
    size_t recipe_node_start = dependency_graph.size();
    size_t recipe_num = 0;
    for (const auto& recipe_record : recipe_registry) {
        for (const auto& recipe : recipe_record.second) {
            IndexName recipe_label = "Recipe " + to_string(recipe_num++);
            //cerr << "adding edges for " << recipe_label << endl;
            //cerr << "\tinputs " << to_string(recipe.input_group()) << endl;
            //cerr << "\toutputs " << to_string(recipe_record.first) << endl;
            graph_label.push_back({recipe_label});
            dependency_graph.emplace_back();
            if (recipe_record.first.size() == 1 &&
                recipe.input_group().count(*recipe_record.first.begin())) {
                // this is an unboxing recipe, only link to the collective input, not individual ingredients
                dependency_graph[graph_idx.at(recipe.input_group())].push_back(dependency_graph.size() - 1);
                //cerr << "\tedge " << to_string(recipe.input_group()) << " -> " << recipe_label << endl;
            }
            else {
                for (auto index_name : recipe.input_group()) {
                    IndexGroup singleton_input{index_name};
                    dependency_graph[graph_idx.at(singleton_input)].push_back(dependency_graph.size() - 1);
                    //cerr << "\tedge " << index_name << " -> " << recipe_label << endl;
                }
            }
            //cerr << "\tedge " << recipe_label << " -> " << to_string(recipe_record.first) << endl;
            dependency_graph.back().push_back(graph_idx.at(recipe_record.first));
        }
    }
    
    // deduplicate any edges
    for (auto& adj : dependency_graph) {
        sort(adj.begin(), adj.end());
        adj.resize(unique(adj.begin(), adj.end()) - adj.begin());
    }
    
#ifdef debug_index_registry
    cerr << "dependency graph:" << endl;
    for (size_t i = 0; i < dependency_graph.size(); ++i) {
        cerr << to_string(graph_label[i]) << endl;
        for (auto j : dependency_graph[i]) {
            cerr << "\t" << to_string(graph_label[j]) << endl;
        }
    }
#endif
    
    // kahn's algorithm to determine a topological order
    vector<size_t> in_degree(dependency_graph.size(), 0);
    for (auto& adj : dependency_graph) {
        for (size_t i : adj) {
            ++in_degree[i];
        }
    }
    
    vector<size_t> stack;
    for (size_t i = 0; i < dependency_graph.size(); ++i) {
        if (in_degree[i] == 0) {
            stack.push_back(i);
        }
    }
    
    vector<size_t> order;
    while (!stack.empty()) {
        size_t i = stack.back();
        stack.pop_back();
        order.push_back(i);
        for (size_t j : dependency_graph[i]) {
            --in_degree[j];
            if (in_degree[j] == 0) {
                stack.push_back(j);
            }
        }
    }
    
    if (order.size() != dependency_graph.size()) {
        cerr << "error:[IndexFile] index dependency graph is not a DAG" << endl;
        
#ifdef debug_index_registry
        // do DFS to find the cycle
        bool found_cycle = false;
        for (size_t i = 0; i < dependency_graph.size() && !found_cycle; ++i) {
            
            vector<bool> traversed(dependency_graph.size(), false);
            vector<bool> stacked(dependency_graph.size(), false);
            // records of (node, next edge to take)
            vector<pair<size_t, size_t>> stack;
            stack.emplace_back(i, 0);
            stacked[i] = true;
            while (!stack.empty()) {
                if (stack.back().second == dependency_graph[stack.back().first].size()) {
                    traversed[stack.back().first] = false;
                    stack.pop_back();
                    continue;
                }
                traversed[stack.back().first] = true;
                size_t next = dependency_graph[stack.back().first][stack.back().second++];
                if (traversed[next]) {
                    size_t j = stack.size() - 1;
                    cerr << "found cycle:" << endl;
                    cerr << "\t" << to_string(graph_label[next]) << endl;
                    while (stack[j].first != next) {
                        cerr << "\t" << to_string(graph_label[stack[j].first]) << endl;
                        --j;
                    }
                    cerr << "\t" << to_string(graph_label[stack[j].first]) << endl;
                    found_cycle = true;
                    break;
                }
                if (!stacked[next]) {
                    stack.emplace_back(next, 0);
                    stacked[next] = true;
                }
            }
        }
#endif
        exit(1);
    }
    
    // convert to return format
    vector<IndexGroup> ordered_identifiers;
    for (size_t i = 0; i < order.size(); ++i) {
        if (order[i] < recipe_node_start) {
            ordered_identifiers.push_back(graph_label[order[i]]);
        }
    }
    
#ifdef debug_index_registry
    cerr << "final order:" << endl;
    for (const auto& identifier : ordered_identifiers) {
        cerr << "\t" << to_string(identifier) << endl;
    }
#endif
    
    return ordered_identifiers;
}

IndexingPlan IndexRegistry::make_plan(const IndexGroup& end_products) const {
    
#ifdef debug_index_registry
    cerr << "generating plan for indexes:" << endl;
    for (const auto& product : end_products) {
        cerr << "\t" << product << endl;
    }
#endif
    
    // get the dependency ordering of the indexes
    vector<IndexGroup> identifier_order = dependency_order();
    map<IndexGroup, size_t> dep_order_of_identifier;
    for (size_t i = 0; i < identifier_order.size(); ++i) {
        dep_order_of_identifier[identifier_order[i]] = i;
    }
    
    // TODO: I'm sure there's a more elegant implementation of this algorithm
    set<RecipeName> plan_elements;
    for (const auto& product : end_products) {
#ifdef debug_index_registry
        cerr << "making a plan for end product " << product << endl;
#endif
        // make a singleton group for the recipe graph
        IndexGroup product_group{product};
        
        // records of (identifier, lowest level requester, ordinal index of recipe selected)
        vector<tuple<size_t, size_t, size_t>> plan_path;
        
        // map dependency priority to lowest level priority that requested this and
        // the number of requesters
        map<size_t, pair<size_t, size_t>, greater<size_t>> queue;
        
        // update the queue to request the inputs of a recipe from the final index on the plan path
        auto request_from_back = [&]() {
            auto make_request = [&](const IndexGroup& inputs) {
                size_t dep_order = dep_order_of_identifier[inputs];
                auto f = queue.find(dep_order);
                if (f == queue.end()) {
                    // no lower-level index has requested this one yet
                    queue[dep_order] = pair<size_t, size_t>(get<0>(plan_path.back()), 1);
                }
                else {
                    // record that one more index is requesting this one
                    f->second.second++;
                }
                
            };
            
            // the index at the back of the plan path is making the request
            auto& requester = identifier_order[get<0>(plan_path.back())];
            // get its next recipe
            auto inputs = recipe_registry.at(requester).at(get<2>(plan_path.back())).input_group();
            
#ifdef debug_index_registry
            cerr << "index(es) " << to_string(requester) << " can be made by a recipe requiring " << to_string(inputs) << endl;
#endif
            
            if (requester.size() == 1 && inputs.count(*requester.begin())) {
                // this is an unboxing recipe, request the whole previous group
                make_request(inputs);
            }
            else {
                // this is not an unboxing recipe, request all of the recipe inputs separately
                for (auto& input_index : inputs) {
                    IndexGroup singleton_input{input_index};
                    make_request(singleton_input);
                }
            }
        };
        
        // update the queue to remove requests to the inputs of a recipe from the final index on the plan path
        auto unrequest_from_back = [&]() {
            auto make_unrequest = [&](const IndexGroup& inputs) {
                size_t input_dep_order = dep_order_of_identifier[inputs];
                auto q = queue.find(input_dep_order);
                if (q != queue.end()) {
                    // there is now one fewer index requesting this index as input
                    --q->second.second;
                    if (q->second.second == 0) {
                        // this is the only index that's requesting this queued index,
                        // so we can remove it from the queue
                        queue.erase(q);
                    }
                }
            };
            
            auto& requester = identifier_order[get<0>(plan_path.back())];
            
#ifdef debug_index_registry
            cerr << "retracting request from " << to_string(requester) << endl;
#endif
            
            if (!all_finished(requester) && recipe_registry.count(requester)) {
                // this index was using a recipe, we need to update its dependencies
                // that are currently in the queue
                
                if (get<2>(plan_path.back()) < recipe_registry.at(requester).size()) {
                    // there's a recipe left
                    auto inputs = recipe_registry.at(requester).at(get<2>(plan_path.back())).input_group();
                    
#ifdef debug_index_registry
                    cerr << to_string(requester) << " made requests from recipe requiring " << to_string(inputs) << endl;
#endif
                    
                    if (requester.size() == 1 && inputs.count(*requester.begin())) {
                        // this is an unboxing recipe, unrequest the whole previous group
                        make_unrequest(inputs);
                    }
                    else {
                        // this is not an unboxing recipe, unrequest all of the recipe inputs separately
                        for (auto& input_index : inputs) {
                            IndexGroup singleton_input{input_index};
                            make_unrequest(singleton_input);
                        }
                    }
                }
            }
        };
        
        // init the queue
        queue[dep_order_of_identifier[product_group]] = pair<size_t, size_t>(identifier_order.size(), 1);
        
        while (!queue.empty()) {
#ifdef debug_index_registry_path_state
            cerr << "new iteration, path:" << endl;
            for (auto pe : plan_path) {
                cerr << "\t" << to_string(identifier_order[get<0>(pe)]) << ", requester: " << (get<1>(pe) == identifier_order.size() ? string("PLAN TARGET") : to_string(identifier_order[get<1>(pe)])) << ", recipe " << get<2>(pe) << endl;
            }
            cerr << "state of queue:" << endl;
            for (auto q : queue) {
                cerr << "\t" << to_string(identifier_order[q.first]) << ", requester: " << (q.second.first == identifier_order.size() ? string("PLAN TARGET") : to_string(identifier_order[q.second.first])) << ", num requesters " << q.second.second << endl;
            }
#endif
            
            // get the latest file in the dependency order
            // that we have left to build
            auto it = queue.begin();
            plan_path.emplace_back(it->first, it->second.first, 0);
            
#ifdef debug_index_registry
            cerr << "dequeue " << to_string(identifier_order[it->first]) << " requested from " << (it->second.first == identifier_order.size() ? string("PLAN TARGET") : to_string(identifier_order[it->second.first])) << " and " << (it->second.second - 1) << " other indexes" << endl;
#endif
            queue.erase(it);
            auto index_group = identifier_order[get<0>(plan_path.back())];
            
            // TODO: am i completely sure that no clobbering will happen if only some of the
            // inputs are provided?
            if (all_finished(index_group)) {
                // this index has been provided, we don't need to use a recipe
#ifdef debug_index_registry
                cerr << "index has been provided as input" << endl;
#endif
                continue;
            }
            else if (recipe_registry.count(index_group)) {
                // there are recipes to make this index, add the requests for the first one
                request_from_back();
            }
            else {
                // we've reached a file that needs to be provided but we don't have it,
                // so now we backtrack until hitting something that has remaining
                // lower priority recipes
#ifdef debug_index_registry
                cerr << "index " << to_string(index_group) << " cannot be made from existing inputs, need to back prune" << endl;
#endif
                
                // prune to requester and advance to its next recipe, as many times as necessary until
                // requester has remaining un-tried recipes
                while (!plan_path.empty() &&
                       (recipe_registry.count(identifier_order[get<0>(plan_path.back())]) ?
                        get<2>(plan_path.back()) >= recipe_registry.at(identifier_order[get<0>(plan_path.back())]).size() : true)) {
                    // there are no remaining recipes to build the last index in the plan
                    
                    // remove items off the plan path until we get to the index that requested
                    // this one
                    size_t requester = get<1>(plan_path.back());
#ifdef debug_index_registry
                    cerr << "pruning path to previous requester: " << (requester == identifier_order.size() ? "PLAN TARGET" : to_string(identifier_order[requester])) << endl;
#endif
                    while (!plan_path.empty() && get<0>(plan_path.back()) != requester) {
                        unrequest_from_back();
                        plan_path.pop_back();
                    }
                    
                    if (!plan_path.empty()) {
                        // the requester should now use its next highest priority recipe
                        unrequest_from_back();
                        // advance to the next recipe
                        ++get<2>(plan_path.back());
                    }
                }
                
                if (!plan_path.empty()) {
                    request_from_back();
                }
            }
            
        }
        
#ifdef debug_index_registry
        cerr << "final plan path for index " << product << ":" << endl;
        for (auto path_elem : plan_path) {
            cerr << "\t" << to_string(identifier_order[get<0>(path_elem)]) << ", from " << (get<1>(path_elem) == identifier_order.size() ? "PLAN START" : to_string(identifier_order[get<1>(path_elem)])) << ", recipe " << get<2>(path_elem) << endl;
        }
#endif
        
        if (plan_path.empty()) {
            // we don't have enough of the inputs to create this index
            throw InsufficientInputException(product, *this);
        }
        
        // record the elements of this plan
        for (size_t i = 0; i < plan_path.size(); ++i) {
            plan_elements.emplace(identifier_order[get<0>(plan_path[i])], get<2>(plan_path[i]));
        }
    }
   
    // Now fill in the plan struct that the recipes need to know how to run.
    IndexingPlan plan;
    
    // Copy over the end products
    std::copy(end_products.begin(), end_products.end(), std::inserter(plan.targets, plan.targets.begin()));
   
    // convert the aggregated plan elements into a forward ordered plan
    std::copy(plan_elements.begin(), plan_elements.end(), std::back_inserter(plan.steps));
    sort(plan.steps.begin(), plan.steps.end(), [&](const RecipeName& a, const RecipeName& b) {
        return dep_order_of_identifier[a.first] < dep_order_of_identifier[b.first];
    });
#ifdef debug_index_registry
    cerr << "full plan including provided files:" << endl;
    for (auto plan_elem : plan.steps) {
        cerr << "\t" << to_string(plan_elem.first) << " " << plan_elem.second << endl;
    }
#endif
    
//    // Now simplify the plan by using joint recipes if possible.
//
//    // First we need to know all the indexes being created, and when
//    map<IndexName, size_t> make_at_step;
//    for (size_t i = 0; i < plan.steps.size(); i++) {
//        make_at_step.emplace(plan.steps[i].first, i);
//    }
//
//    // We also need to know what steps we've already simplified, or are inputs.
//    // We don't want to apply overlapping simplifications.
//    vector<bool> fixed_step(plan.steps.size(), false);
//
//    for (size_t i = 0; i < plan.steps.size(); i++) {
//        auto& recipe = plan.steps.at(i);
//        if (get_index(recipe.first)->is_finished()) {
//            // This is already provided and ineligible for simplification.
//            fixed_step[i] = true;
//        }
//    }
//
//    for (auto& simplification : simplifications) {
//        // For each set of output indexes from a simplification
//
//#ifdef debug_index_registry
//        cerr << "Consider simplification to jointly make:" << endl;
//        for (auto& recipe: simplification.second) {
//            cerr << "\t" << to_string(recipe.first) << endl;
//        }
//#endif
//
//        // Determine if we are making all the products of the simplification,
//        // and those products have not been involved in prior simplifications
//        bool making_all_products_unsimplified = true;
//        // And if so, the first step at which we are making any
//        size_t first_step = numeric_limits<size_t>::max();
//        for (auto& product_recipe : simplification.second) {
//            const IndexName& product_name = product_recipe.first;
//
//            auto found = make_at_step.find(product_name);
//            if (found == make_at_step.end()) {
//                // We aren't making this product
//
//#ifdef debug_index_registry
//                cerr << "We are not making " << to_string(product_name) << endl;
//#endif
//
//                making_all_products_unsimplified = false;
//                break;
//            }
//
//            if (fixed_step[found->second]) {
//                // We are making this product but we already simplified it or took it as input
//
//#ifdef debug_index_registry
//                cerr << "We cannot further simplify making " << to_string(product_name) << endl;
//#endif
//
//                making_all_products_unsimplified = false;
//                break;
//            }
//
//#ifdef debug_index_registry
//            cerr << "We are making " << to_string(product_name) << " at step " << found->second << endl;
//#endif
//
//            first_step = min(first_step, found->second);
//        }
//
//        if (!making_all_products_unsimplified) {
//            // This simplification can't be used becuase it makes extra
//            // products, or products that are already simplified.
//
//#ifdef debug_index_registry
//            cerr << "We are not making all the products for this simplification, or some products cannot be further simplified" << endl;
//#endif
//
//            continue;
//        }
//
//#ifdef debug_index_registry
//        cerr << "To simplify, all inputs will need to be available before step " << first_step << endl;
//#endif
//
//        // See what we have available before the first step
//        set<IndexName> available_in_time;
//        for (size_t i = 0; i < first_step; i++) {
//            available_in_time.insert(plan.steps[i].first);
//        }
//
//        // See if it's all the inputs the simplification needs
//        bool all_available = true;
//        for (auto& needed : simplification.first) {
//            if (!available_in_time.count(needed)) {
//#ifdef debug_index_registry
//                cerr << "We are not making " << to_string(needed) << " in time or at all." << endl;
//#endif
//                all_available = false;
//                break;
//            }
//        }
//
//        if (!all_available) {
//            // This simplification can't be used because not all its inputs are available in time.
//
//#ifdef debug_index_registry
//            cerr << "Not all inputs will be available in time." << endl;
//#endif
//
//            continue;
//        }
//
//#ifdef debug_index_registry
//        cerr << "All inputs will be available in time. Apply simplification!" << endl;
//#endif
//
//        for (auto& recipe : simplification.second) {
//            // Replace each relevant step with the corresponding joint step for that index.
//            size_t step_to_simplify = make_at_step.at(recipe.first);
//            plan.steps.at(step_to_simplify) = recipe;
//            fixed_step[step_to_simplify] = true;
//        }
//    }
//
//#ifdef debug_index_registry
//    cerr << "plan after simplification:" << endl;
//    for (auto plan_elem : plan.steps) {
//        cerr << "\t" << to_string(plan_elem.first) << " " << plan_elem.second << endl;
//    }
//#endif

    // Now remove the input data from the plan
    plan.steps.resize(remove_if(plan.steps.begin(), plan.steps.end(), [&](const RecipeName& recipe_choice) {
        return all_finished(recipe_choice.first);
    }) - plan.steps.begin());
    
    // The plan has methods that can come back and modify the registry.
    // We're not going to call any of them, but we have to hand off a non-const
    // pointer to ourselves so the plan can modify us later.
    plan.registry = const_cast<IndexRegistry*>(this);
    
    return plan;
}

vector<vector<string>> IndexRegistry::execute_recipe(const RecipeName& recipe_name, const IndexingPlan* plan,
                                                     AliasGraph& alias_graph) {
    const auto& recipes = recipe_registry.at(recipe_name.first);
    assert(recipe_name.second < recipes.size());
    const auto& index_recipe = recipes.at(recipe_name.second);
    if (recipe_name.first.size() > 1 || !index_recipe.input_group().count(*recipe_name.first.begin())) {
        // we're not in an unboxing recipe (in which case not all of the indexes might have been
        // unboxed yet, in which case they appear unfinished)
        for (auto input : index_recipe.inputs) {
            assert(input->is_finished());
        }
    }
    return index_recipe.execute(plan, alias_graph, recipe_name.first);;
}

string IndexRegistry::to_dot() const {
    return to_dot(vector<IndexName>());
}

string IndexRegistry::to_dot(const vector<IndexName>& targets) const {
    
    
    stringstream strm;
    strm << "digraph recipegraph {" << endl;
    
    set<IndexGroup> plan_targets;
    for (const auto& target : targets) {
        plan_targets.insert({target});
    }
    set<RecipeName> plan_elements;
    set<IndexGroup> plan_indexes;
    if (!targets.empty()) {
        IndexingPlan plan;
        try {
            IndexGroup target_group(targets.begin(), targets.end());
            plan = make_plan(target_group);
        }
        catch (InsufficientInputException ex) {
#ifdef debug_index_registry
            cerr << ex.what() << endl;
#endif
            strm << "labelloc=\"t\";" << endl;
            strm << "label=\"Insufficient input to create targets\";" << endl;
        }
        for (const auto& plan_elem : plan.steps) {
            plan_elements.insert(plan_elem);
            plan_indexes.insert(plan_elem.first);
        }
    }
    
    // gather all singletons and products of recipes, which will be the index nodes
    set<IndexGroup> all_indexes;
    for (const auto& index_record : index_registry) {
        all_indexes.insert({index_record.first});
    }
    for (const auto& recipe_record : recipe_registry) {
        all_indexes.insert(recipe_record.first);
    }
    
    map<IndexGroup, string> index_to_dot_id;
    size_t index_idx = 0;
    for (const auto& index_group : all_indexes) {
        index_to_dot_id[index_group] = "I" + to_string(index_idx);
        ++index_idx;
        strm << index_to_dot_id[index_group] << "[label=\"" << to_string(index_group) << "\" shape=box";
        if (all_finished(index_group)) {
            strm << " style=\"filled,bold\" fillcolor=lightgray";
        }
        else if (plan_targets.count(index_group)) {
            strm << " style=\"filled,bold\" fillcolor=lightblue";
        }
        else if (plan_indexes.count(index_group)) {
            strm << " style=bold";
        }
        strm << "];" << endl;
    }
    string unselected_col = targets.empty() ? "black" : "gray33";
    size_t recipe_idx = 0;
    for (const auto& recipe_record : recipe_registry) {
        const auto& recipes = recipe_record.second;
        for (size_t priority_idx = 0; priority_idx < recipes.size(); ++priority_idx, ++recipe_idx) {
            const auto& recipe = recipes[priority_idx];
            string recipe_dot_id = "R" + to_string(recipe_idx);
            bool recipe_in_plan = plan_elements.count(RecipeName(recipe_record.first, priority_idx));
            if (recipe_in_plan) {
                strm << recipe_dot_id << "[label=\"" << priority_idx << "\" shape=circle style=bold];" << endl;
                strm << recipe_dot_id << " -> " << index_to_dot_id[recipe_record.first] << "[style=bold];" << endl;
            }
            else {
                strm << recipe_dot_id << "[label=\"" << priority_idx << "\" shape=circle];" << endl;
                strm << recipe_dot_id << " -> " << index_to_dot_id[recipe_record.first] << " [color=" << unselected_col << "];" << endl;
            }
            auto input_group = recipe.input_group();
            if (recipe_record.first.size() == 1 && input_group.count(*recipe_record.first.begin())) {
                // unboxing recipe, link directly to group
                if (recipe_in_plan) {
                    strm << index_to_dot_id[input_group] << " -> " << recipe_dot_id << "[style=bold];" << endl;
                }
                else {
                    strm << index_to_dot_id[input_group] << " -> " << recipe_dot_id << " [color=" << unselected_col << "];" << endl;
                }
            }
            else {
                // not an unboxing recipe, link to singletons
                for (const auto& input : input_group) {
                    if (recipe_in_plan) {
                        strm << index_to_dot_id[IndexGroup{input}] << " -> " << recipe_dot_id << "[style=bold];" << endl;
                    }
                    else {
                        strm << index_to_dot_id[IndexGroup{input}] << " -> " << recipe_dot_id << " [color=" << unselected_col << "];" << endl;
                    }
                }
            }
            
        }
    }
    strm << "}" << endl;
    return strm.str();
}

IndexFile::IndexFile(const IndexName& identifier, const string& suffix) : identifier(identifier), suffix(suffix) {
    // nothing more to do
}

bool IndexFile::is_finished() const {
    return !filenames.empty();
}

const IndexName& IndexFile::get_identifier() const {
    return identifier;
}

const string& IndexFile::get_suffix() const {
    return suffix;
}

const vector<string>& IndexFile::get_filenames() const {
    return filenames;
}

void IndexFile::provide(const vector<string>& filenames) {
    // append all filenames
    // TODO: would it be better to sometimes error check that the file isn't a duplicate?
    for (const string& filename : filenames) {
        this->filenames.emplace_back(filename);
    }
    provided_directly = true;
}

void IndexFile::assign_constructed(const vector<string>& filenames) {
    this->filenames = filenames;
    provided_directly = false;
}

bool IndexFile::was_provided_directly() const {
    return provided_directly;
}

void IndexFile::reset() {
    filenames.clear();
    provided_directly = false;
}

IndexRecipe::IndexRecipe(const vector<const IndexFile*>& inputs,
                         const RecipeFunc& exec) :
    exec(exec), inputs(inputs)
{
    // nothing more to do
}

vector<vector<string>> IndexRecipe::execute(const IndexingPlan* plan, AliasGraph& alias_graph,
                                            const IndexGroup& constructing) const {
    return exec(inputs, plan, alias_graph, constructing);
}

IndexGroup IndexRecipe::input_group() const {
    IndexGroup group;
    for (auto input : inputs) {
        group.insert(input->get_identifier());
    }
    return group;
}

void AliasGraph::register_alias(const IndexName& aliasor, const IndexFile* aliasee) {
    assert(aliasee->get_identifier() != aliasor);
    graph[aliasee->get_identifier()].emplace_back(aliasor);
}

vector<pair<IndexName, vector<IndexName>>> AliasGraph::non_intermediate_aliases(const IndexingPlan* plan,
                                                                                bool keep_all) const {
    
#ifdef debug_index_registry
    cerr << "finding non intermediate aliases in alias graph" << endl;
    for (const auto& adj : graph) {
        cerr << adj.first << ":" << endl;
        for (const auto& dest : adj.second) {
            cerr << "\t" << dest << endl;
        }
    }
#endif
    
    vector<pair<IndexName, vector<IndexName>>> aliases;
    
    // find the heads in the graph (the origins of aliasing chains)
    unordered_set<IndexName> heads;
    for (const auto& adj : graph) {
        heads.insert(adj.first);
    }
    for (const auto& adj : graph) {
        for (const auto& dest : adj.second) {
            if (heads.count(dest)) {
                heads.erase(dest);
            }
        }
    }
    
    for (const auto& head : heads) {
        
#ifdef debug_index_registry
        cerr << "starting a DFS from head index " << head << endl;
#endif
        
        // do DFS out from this head to identify aliasors
        vector<IndexName> non_inmdt_aliasors;
        vector<IndexName> stack(1, head);
        unordered_set<IndexName> stacked{head};
        while (!stack.empty()) {
            auto here = stack.back();
            stack.pop_back();
            if (!plan->is_intermediate(here) || keep_all) {
                non_inmdt_aliasors.push_back(here);
            }
            if (graph.count(here)) {
                for (const auto& dest : graph.at(here)) {
                    if (!stacked.count(dest)) {
                        stack.push_back(dest);
                        stacked.insert(dest);
                    }
                }
            }
        }
        
        if (!non_inmdt_aliasors.empty()) {
            aliases.emplace_back(head, move(non_inmdt_aliasors));
        }
    }
#ifdef debug_index_registry
    cerr << "identified aliases" << endl;
    for (const auto& alias_record : aliases) {
        cerr << alias_record.first << ":" << endl;
        for (auto aliasor : alias_record.second) {
            cerr << "\t" << aliasor << endl;
        }
    }
#endif
    return aliases;
}

InsufficientInputException::InsufficientInputException(const IndexName& target,
                                                       const IndexRegistry& registry) :
    runtime_error("Insufficient input to create " + target), target(target), inputs(registry.completed_indexes())
{
    // nothing else to do
}

const char* InsufficientInputException::what() const throw () {
    stringstream ss;
    ss << "Inputs" << endl;
    for (const auto& input : inputs) {
        ss << "\t" << input << endl;
    }
    ss << "are insufficient to create target index " << target << endl;
    return ss.str().c_str();
}

}

