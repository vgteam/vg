// index_registry.cpp: index registry system implementation

#include "index_registry.hpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <map>

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
#include "source_sink_overlay.hpp"

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
static string to_string(const vg::IndexName& name) {
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

int64_t get_file_size(const string& filename) {
    // get the file size
    ifstream infile(filename);
    infile.seekg(0, ios::end);
    return infile.tellg();
}

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

double format_multiplier() {
    
    // the ratio to the HashGraph memory usage, as estimated by a couple of graphs
    // that Jordan had laying around when he wrote this
    // in case it's useful later: XG is ~0.231
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


IndexingParameters::MutableGraphImplementation IndexingParameters::mut_graph_impl = HashGraph;
int IndexingParameters::max_node_size = 32;
int IndexingParameters::pruning_max_node_degree = 128;
int IndexingParameters::pruning_walk_length = 24;
int IndexingParameters::pruning_max_edge_count = 3;
int IndexingParameters::pruning_min_component_size = 33;
int IndexingParameters::gcsa_initial_kmer_length = gcsa::Key::MAX_LENGTH;
int IndexingParameters::gcsa_doubling_steps = gcsa::ConstructionParameters::DOUBLING_STEPS;
bool IndexingParameters::verbose = false;

IndexRegistry VGIndexes::get_vg_index_registry() {
    
    IndexRegistry registry;
    
    /*********************
     * Register all of the VG indexes and input files
     ***********************/
    
    /// Data files
    registry.register_index({"Reference FASTA"}, "fasta");
    registry.register_index({"VCF"}, "vcf");
    registry.register_index({"Phased VCF"}, "phased.vcf");
    registry.register_index({"Insertion Sequence FASTA"}, "insertions.fasta");
    registry.register_index({"Reference GFA"}, "gfa");
    
    /// True indexes
    registry.register_index({"VG + Variant Paths"}, "varpaths.vg");
    registry.register_index({"VG"}, "vg");
    registry.register_index({"XG"}, "xg");
    registry.register_index({"GBWT"}, "gbwt");
    registry.register_index({"Pruned VG"}, "pruned.vg");
    // TODO: the suffixes are pretty awkward without unboxing, which makes
    // it so that we have to register every combination as a separate index
    registry.register_index({"MaxNodeID"}, "maxid.txt");
    registry.register_index({"UnfoldedNodeMapping"}, "unfolded.mapping");
    registry.register_index({"Haplotype-Pruned VG"}, "haplo.vg");
    registry.register_index({"Haplotype-Pruned VG", "UnfoldedNodeMapping"}, "haplopruned.vg");
    registry.register_index({"GCSA", "LCP"}, "gcsa");
    
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
    registry.register_recipe({"VCF"}, {{"Phased VCF"}},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                const IndexName& constructing) {
        return inputs.front()->get_filenames();
    });
    
    ////////////////////////////////////
    // VG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering VG recipes" << endl;
#endif
    
    // strip alt allele paths from a graph that has them
    registry.register_recipe({"VG"}, {{"VG + Variant Paths"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Stripping allele paths from VG." << endl;
        }
        
        vector<string> output_names;
        for (int i = 0; i < inputs.at(0)->get_filenames().size(); ++i) {
            // test streams for I/O
            ifstream infile;
            init_in(infile, inputs.at(0)->get_filenames()[i]);
            
            string output_name = plan->prefix(constructing);
            if (inputs.front()->get_filenames().size() != 1) {
                output_name += "." + to_string(i);
            }
            output_name += "." + plan->suffix(constructing);
            output_names.push_back(output_name);
            
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
        }
        
        // return the filename(s)
        return output_names;
    });
        
    // meta-recipe for creating a VG from a GFA
    auto construct_from_gfa = [&](const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing,
                                  nid_t* max_node_id_out) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing VG graph from GFA input." << endl;
        }
        string output_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
        ofstream outfile;
        init_out(outfile, output_name);
        auto graph = init_mutable_graph();
        
        // make the graph from GFA
        try {
            algorithms::gfa_to_path_handle_graph(inputs.at(0)->get_filenames().front(), graph.get(), true,
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
        return vector<string>(1, output_name);
    };
    
    registry.register_recipe({"VG"}, {{"Reference GFA"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return construct_from_gfa(inputs, plan, constructing, nullptr);
    });
    
    // A meta-recipe to make VG files using the Constructor
    // Expects inputs to be ordered: FASTA, VCF[, Insertion FASTA]
    auto construct_with_constructor = [&](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          const IndexName& constructing,
                                          bool alt_paths,
                                          nid_t* max_node_id_out) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing VG graph from FASTA and VCF input." << endl;
        }
        assert(inputs.size() == 2 || inputs.size() == 3);
        if ((inputs[0]->get_filenames().size() != 1 && inputs[0]->get_filenames().size() != inputs[1]->get_filenames().size()) ||
            (inputs[1]->get_filenames().size() != 1 && inputs[0]->get_filenames().size() != inputs[1]->get_filenames().size())) {
            cerr << "[IndexRegistry]: When constructing graph from multiple FASTAs and multiple VCFs, the FASTAs and VCFs must be matched 1-to-1, but input contains " <<  inputs[0]->get_filenames().size() << " FASTA files and " << inputs[1]->get_filenames().size() << " VCF files." << endl;
            exit(1);
        }
        
        vector<string> insertions;
        if (inputs.size() == 3) {
            insertions = inputs.back()->get_filenames();
        }
        
        // TODO: allow contig renaming through Constructor::add_name_mapping
        // TODO: actually do the chunking based on available memory and in parallel if possible
        
        nid_t max_node_id = 0;
        vector<string> output_names;
        size_t i = 0, j = 0;
        while (i < inputs[0]->get_filenames().size() && j < inputs[1]->get_filenames().size()) {
            
            // init and configure the constructor
            Constructor constructor;
            constructor.do_svs = true;
            constructor.alt_paths = alt_paths;
            constructor.max_node_size = IndexingParameters::max_node_size;
            constructor.show_progress = IndexingParameters::verbose;
            if (inputs[0]->get_filenames().size() != 1 && inputs[1]->get_filenames().size() == 1) {
                // we have multiple FASTA but only 1 VCF, so we'll limit the
                // constructor to the contigs of this FASTA for this run
                FastaReference ref;
                ref.open(inputs[0]->get_filenames()[i]);
                for (const string& seqname : ref.index->sequenceNames) {
                    constructor.allowed_vcf_names.insert(seqname);
                }
            }
            
            string output_name = plan->prefix(constructing);
            if (inputs[0]->get_filenames().size() != 1
                || inputs[1]->get_filenames().size() != 1) {
                output_name += "." + to_string(max(i, j));
            }
            output_name += "." + plan->suffix(constructing);
            ofstream outfile;
            init_out(outfile, output_name);
            
            auto graph = init_mutable_graph();
            
            vector<string> fasta(1, inputs[0]->get_filenames()[i]);
            vector<string> vcf(1, inputs[1]->get_filenames()[j]);
            
            // do the construction
            constructor.construct_graph(fasta, vcf,
                                        insertions,
                                        graph.get());
            
            // save the file
            vg::io::save_handle_graph(graph.get(), outfile);
            
            output_names.push_back(output_name);
            
            max_node_id = max(max_node_id, graph->max_node_id());
            
            // awkward incrementation that allows 1x1, 1xN, Nx1, and NxN matching
            // of VCFs and FASTAs
            if (inputs[0]->get_filenames().size() == 1 && inputs[1]->get_filenames().size() == 1) {
                break;
            }
            if (inputs[0]->get_filenames().size() > 1) {
                ++i;
            }
            if (inputs[1]->get_filenames().size() > 1) {
                ++j;
            }
        }
        
        if (inputs[0]->get_filenames().size() == 1 && inputs[1]->get_filenames().size() != 1) {
            
            // we will have added components for all of the paths that are also
            // VCF sequences, but we may have FASTA sequences that don't
            // correspond to any VCF sequence (e.g. unlocalized contigs of
            // decoys)
            
            // TODO: how can we keep track of all of the contigs seen in the Constructor
            // so that we know which ones to include in this component?
            
            string output_name = (plan->prefix(constructing)
                                  + "." + to_string(inputs[1]->get_filenames().size())
                                  + "." + plan->suffix(constructing));
            
            // FIXME: punting on this for now
            // FIXME: also need to add a dummy VCF, but how to add it to const index?
            // maybe i should add it as a second index, along with this graph...
        }
        
        if (output_names.size() > 1) {
            // join the ID spaces
            VGset graph_set(output_names);
            max_node_id = graph_set.merge_id_space();
        }
        
        if (max_node_id_out) {
            // TODO: maybe use this for initializing a NodeMapping
            *max_node_id_out = max_node_id;
        }
        
        // return the filename(s)
        return output_names;
    };
    
    // the specific instantiations of the meta-recipe above
    registry.register_recipe({"VG"}, {{"Reference FASTA"}, {"VCF"}, {"Insertion Sequence FASTA"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, nullptr);
    });
    registry.register_recipe({"VG"}, {{"Reference FASTA"}, {"VCF"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, nullptr);
    });
    registry.register_recipe({"VG + Variant Paths"}, {{"Reference FASTA"}, {"Phased VCF"}, {"Insertion Sequence FASTA"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, nullptr);
    });
    registry.register_recipe({"VG + Variant Paths"}, {{"Reference FASTA"}, {"Phased VCF"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, nullptr);
    });
    
    
    ////////////////////////////////////
    // XG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering XG recipes" << endl;
#endif
    
    registry.register_recipe({"XG"}, {{"Reference GFA"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing XG graph from GFA input." << endl;
        }
        string output_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
        ofstream outfile;
        init_out(outfile, output_name);
        
        xg::XG xg_index;
        xg_index.from_gfa(inputs.front()->get_filenames().front());
        
        vg::io::save_handle_graph(&xg_index, outfile);
        
        // return the filename
        return vector<string>(1, output_name);
    });
    
    registry.register_recipe({"XG"}, {{"VG"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing XG graph from VG graph." << endl;
        }
        
        string output_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
        ofstream outfile;
        init_out(outfile, output_name);
        
        xg::XG xg_index;
        
        if (inputs.at(0)->get_filenames().size() == 1) {
            // we do the one-graph conversion directly, which is more efficient than the
            // VGset option
            
            // test streams for I/O
            ifstream infile;
            init_in(infile, inputs.at(0)->get_filenames().front());
            
            unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            
            xg_index.from_path_handle_graph(*graph);
        }
        else {
            // the inefficient 3-pass, multi-graph construction algorithm
            
            // make a mutable copy of the graph names
            vector<string> graph_files;
            for (const string& graph_file : inputs.at(0)->get_filenames()) {
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
        return vector<string>(1, output_name);
    });
    
    ////////////////////////////////////
    // MaxNodeID Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering MaxNodeID recipes" << endl;
#endif
    
    registry.register_recipe({"MaxNodeID"}, {{"VG"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Determining node ID interval." << endl;
        }
        // TODO: this is pretty unoptimized in that we have to load the whole graph just
        // to read the max node id
        
        // test I/O and make non-const copy of file names
        vector<string> graph_files;
        for (const string& graph_file : inputs.at(0)->get_filenames()) {
            ifstream infile;
            init_in(infile, graph_file);
            graph_files.push_back(graph_file);
        }
        string output_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
        ofstream outfile;
        init_out(outfile, output_name);
        
        VGset graph_set(graph_files);
        nid_t max_node_id = graph_set.max_node_id();
        
        outfile << max_node_id;
        
        return vector<string>(1, output_name);
    });
    
    ////////////////////////////////////
    // GBWT Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering GBWT recipes" << endl;
#endif
    
    // TODO: add Jouni's chunked workflow here
    registry.register_recipe({"GBWT"}, {{"VG + Variant Paths"}, {"Phased VCF"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing GBWT from VG graph and phased VCF input." << endl;
        }
        if ((inputs[0]->get_filenames().size() != 1 && inputs[0]->get_filenames().size() != inputs[1]->get_filenames().size()) ||
            (inputs[1]->get_filenames().size() != 1 && inputs[0]->get_filenames().size() != inputs[1]->get_filenames().size())) {
            cerr << "[IndexRegistry]: When constructing GBWT from multiple graphs and multiple VCFs, the graphs and VCFs must be matched 1-to-1, but input contains " <<  inputs[0]->get_filenames().size() << " graphs and " << inputs[1]->get_filenames().size() << " VCF files." << endl;
            exit(1);
        }
        
        vector<string> gbwt_names;
        
        auto gbwt_job = [&](size_t i, const PathHandleGraph& graph) {
            string gbwt_name;
            if (inputs[1]->get_filenames().size() != 1) {
                // multiple components, so make a temp file that we will merge later
                gbwt_name = temp_file::create();
            }
            else {
                // one component, so we will actually save the output
                gbwt_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
            }
            
            HaplotypeIndexer haplotype_indexer;
            haplotype_indexer.show_progress = IndexingParameters::verbose;
            
            vector<string> parse_files = haplotype_indexer.parse_vcf(inputs[1]->get_filenames()[i],
                                                                     graph);
            // we don't want to do this here, since we're reusing the graph
            //graph.reset(); // Save memory by deleting the graph.
            unique_ptr<gbwt::DynamicGBWT> gbwt_index = haplotype_indexer.build_gbwt(parse_files);
            
            vg::io::VPKG::save(*gbwt_index, gbwt_name);
            
            gbwt_names.push_back(gbwt_name);
        };
        
        if (inputs[0]->get_filenames().size() == 1) {
            
            // we only have one graph, so we can save time by loading it only one time
            
            // test streams for I/O
            ifstream infile;
            init_in(infile, inputs.at(0)->get_filenames().front());
            
            unique_ptr<PathHandleGraph> graph
                = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            
            for (size_t i = 0; i < inputs[1]->get_filenames().size(); ++i) {
                gbwt_job(i, *graph);
            }
        }
        else {
            // FIXME: it doesn't seem possible to do many component graphs with 1 VCF
            for (size_t i = 0; i < inputs[0]->get_filenames().size(); ++i) {
                ifstream infile;
                init_in(infile, inputs.at(0)->get_filenames()[i]);
                unique_ptr<PathHandleGraph> graph
                    = vg::io::VPKG::load_one<PathHandleGraph>(infile);
                gbwt_job(i, *graph);
            }
        }
        
        if (gbwt_names.size() > 1) {
            if (IndexingParameters::verbose) {
                cerr << "[IndexRegistry]: Merging contig GBWTs." << endl;
            }
            // we also need to merge the GBWTs
            
            string merged_gbwt_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
            ofstream outfile;
            init_out(outfile, merged_gbwt_name);
            
            vector<gbwt::GBWT> gbwt_indexes(gbwt_names.size());
            for (size_t i = 0; i < gbwt_names.size(); ++i) {
                load_gbwt(gbwt_names[i], gbwt_indexes[i], IndexingParameters::verbose);
            }
            gbwt::GBWT merged(gbwt_indexes);
            merged.serialize(outfile);
            
            return vector<string>(1, merged_gbwt_name);
        }
        else {
            return gbwt_names;
        }
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
                           const IndexName& constructing) {
        
        // we only want to focus on two specific recipes
        assert(inputs.size() == 2 || inputs.size() == 3);
        bool using_haplotypes = inputs.size() == 3;
        
        // test streams for I/O
        ifstream infile_gbwt, infile_max_id;
        init_in(infile_max_id, inputs.at(1)->get_filenames().front());
        if (using_haplotypes) {
            init_in(infile_gbwt, inputs.at(2)->get_filenames().front());
        }
        
        // read the max node ID (across all chunks)
        nid_t max_node_id;
        infile_max_id >> max_node_id;
        
        string mapping_output_name;
        if (using_haplotypes) {
            
            gcsa::NodeMapping mapping(max_node_id + 1);
            
            mapping_output_name = plan->prefix(constructing) + "." + plan->suffix(constructing) + ".mapping";
            
            ofstream mapping_file;
            init_out(mapping_file, mapping_output_name);
            mapping.serialize(mapping_file);
        }
        
        // TODO: is it possible to do this in parallel? i'm worried about races
        // on the node mapping, and it seems like it definitely wouldn't work
        // for haplotype unfolding, which needs to modify the NodeMapping
        vector<string> output_names;
        for (size_t i = 0; i < inputs.at(0)->get_filenames().size(); ++i) {
            
            ifstream infile_vg;
            init_in(infile_vg, inputs.at(0)->get_filenames()[i]);
            
            string vg_output_name = plan->prefix(constructing);
            if (inputs.at(0)->get_filenames().size() != 1) {
                vg_output_name += "." + to_string(i);
            }
            vg_output_name += "." + plan->suffix(constructing);
            
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
                init_in(infile_unpruned_vg, inputs.at(0)->get_filenames()[i]);
                
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
                    // we can expand out complex regions using haplotypes
                    
                    // copy the node mapping so that we have a copy that we can alter without
                    // potentially messing up the input for other recipes
                    unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(infile_gbwt);
                    PhaseUnfolder unfolder(*unpruned_graph, *gbwt_index, max_node_id + 1);
                    unfolder.read_mapping(mapping_output_name);
                    unfolder.unfold(*graph, IndexingParameters::verbose);
                    unfolder.write_mapping(mapping_output_name);
                }
            }
            
            vg::io::save_handle_graph(graph.get(), outfile_vg);
            
            output_names.push_back(vg_output_name);
        }
                
        if (using_haplotypes) {
            output_names.push_back(mapping_output_name);
        }
        return output_names;
    };
    
    registry.register_recipe({"Pruned VG"}, {{"VG"}, {"MaxNodeID"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Haplotype-Pruned VG", "UnfoldedNodeMapping"}, {{"VG"}, {"MaxNodeID"}, {"GBWT"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing with GBWT unfolding." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // Unboxing
    ////////////////////////////////////
    
    // TODO: in order to have automatic unboxing we would need to have the
    // filenames in separate vectors rather than counting on knowing the
    // layout within the vector
    registry.register_recipe({"Haplotype-Pruned VG"}, {{"Haplotype-Pruned VG", "UnfoldedNodeMapping"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        const auto& filenames = inputs.front()->get_filenames();
        return vector<string>(filenames.begin(), filenames.end() - 1);
    });
    
    registry.register_recipe({"UnfoldedNodeMapping"}, {{"Haplotype-Pruned VG", "UnfoldedNodeMapping"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        return vector<string>(1, inputs.front()->get_filenames().back());
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
                              const IndexName& constructing) {
        
        if (IndexingParameters::verbose) {
            cerr << "[IndexRegistry]: Constructing GCSA/LCP indexes." << endl;
        }
        
        bool unfolded = inputs.size() == 2;
        
        // test streams for I/O
        ifstream infile_mapping;
        string mapping_filename;
        if (unfolded) {
            mapping_filename = inputs.at(1)->get_filenames().back();
            init_in(infile_mapping, mapping_filename);
        }
        string gcsa_output_name = plan->prefix(constructing) + "." + plan->suffix(constructing);
        string lcp_output_name = gcsa_output_name + ".lcp";
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
        
        vector<string> graph_filenames = inputs.at(0)->get_filenames();
        
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
        
        vg::io::VPKG::save(gcsa_index, gcsa_output_name);
        vg::io::VPKG::save(lcp_array, lcp_output_name);
        
        return vector<string>{gcsa_output_name, lcp_output_name};
    };
    
    registry.register_recipe({"GCSA", "LCP"}, {{"Haplotype-Pruned VG"}, {"UnfoldedNodeMapping"}},
                            [&](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                const IndexName& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"GCSA", "LCP"}, {{"Pruned VG"}},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 const IndexName& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    
    return registry;
}

vector<IndexName> VGIndexes::get_default_map_indexes() {
    vector<IndexName> indexes{
        {"XG"},
        {"GCSA", "LCP"}
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_mpmap_indexes() {
    vector<IndexName> indexes{
        {"Spliced XG"},
        {"Spliced Distance"},
        {"Spliced GCSA", "Spliced LCP"},
        {"Haplotype-Transcript GBWT"}
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_giraffe_indexes() {
    vector<IndexName> indexes{
        {"GBWT"},
        {"GBWTGraph"},
        {"Distance"},
        {"Minimizer"}
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
    
string IndexingPlan::prefix(const IndexName& identifier) const {
    // TODO: recipes that are simply aliasing a more general file will sometimes
    // ignore the prefix, for this system to work those can never be final outputs...
    string index_prefix;
    if (registry->keep_intermediates || !is_intermediate(identifier)) {
        // we're saving this file, put it at the output prefix
        index_prefix = registry->output_prefix;
    } else {
        // we're not saving this file, make it temporary
        index_prefix = registry->get_work_dir() + "/" + sha1sum(to_string(identifier));
    }
    return index_prefix;
}
    
string IndexingPlan::suffix(const IndexName& identifier) const {
    return registry->get_index(identifier)->get_suffix();
}
    

IndexRegistry::~IndexRegistry() {
    if (!work_dir.empty()) {
        // Clean up our work directory with its temporary indexes.
        temp_file::remove(work_dir);
        work_dir.clear();
    }
}

IndexRegistry::IndexRegistry(IndexRegistry&& other) :
    registry(std::move(other.registry)),
    registered_suffixes(std::move(other.registered_suffixes)),
    work_dir(std::move(other.work_dir)),
    output_prefix(std::move(other.output_prefix)),
    keep_intermediates(std::move(other.keep_intermediates)) {
    
    // Make sure other doesn't delete our work dir when it goes away
    other.work_dir.clear();
}

IndexRegistry& IndexRegistry::operator=(IndexRegistry&& other) {
    registry = std::move(other.registry);
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
    auto plan = make_plan(identifiers);
    
    // execute the plan
    for (auto& step : plan.steps) {
        
        auto index = get_index(step.first);
        
        index->execute_recipe(step.second, &plan);
    }
    
    // clean up intermediate files
    if (!keep_intermediates) {
        // collect the names of everything we want to keep (i.e. considered non-
        // intermediate by at least one index)
        unordered_set<string> filenames_to_keep;
        for (const auto& registered_index : registry) {
            if (!plan.is_intermediate(registered_index.second->get_identifier())) {
                for (auto& filename : registered_index.second->get_filenames()) {
                    filenames_to_keep.insert(filename);
                }
            }
        }
        
        // delete everything else
        for (const auto& registered_index : registry) {
            for (auto& filename : registered_index.second->get_filenames()) {
                if (!filenames_to_keep.count(filename)) {
                    std::remove(filename.c_str());
                }
            }
        }
    }
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
    if (registry.count(identifier)) {
        cerr << "error:[IndexRegistry] index registry contains a duplicated identifier: " << to_string(identifier) << endl;
        exit(1);
    }
    if (registered_suffixes.count(suffix)) {
        cerr << "error:[IndexRegistry] index registry contains a duplicated suffix: " << suffix << endl;
        exit(1);
    }
    registry[identifier] = unique_ptr<IndexFile>(new IndexFile(identifier, suffix));
    registered_suffixes.insert(suffix);
}


void IndexRegistry::provide(const IndexName& identifier, const string& filename) {
    provide(identifier, vector<string>(1, filename));
}

void IndexRegistry::provide(const IndexName& identifier, const vector<string>& filenames) {
    get_index(identifier)->provide(filenames);
}

vector<IndexName> IndexRegistry::completed_indexes() const {
    vector<IndexName> indexes;
    for (const auto& index : registry) {
        if (index.second->is_finished()) {
            indexes.push_back(index.first);
        }
    }
    return indexes;
}

RecipeName IndexRegistry::register_recipe(const IndexName& identifier,
                                          const vector<IndexName>& input_identifiers,
                                          const RecipeFunc& exec) {
    if (!registry.count(identifier)) {
        cerr << "error:[IndexRegistry] cannot register recipe for unregistered index " << to_string(identifier) << endl;
        exit(1);
    }
    vector<const IndexFile*> inputs;
    for (const auto& input_identifier : input_identifiers) {
        if (!registry.count(input_identifier)) {
            cerr << "error:[IndexRegistry] cannot register recipe from unregistered index " << to_string(input_identifier) << endl;
            exit(1);
        }
        inputs.push_back(get_index(input_identifier));
    }
    return get_index(identifier)->add_recipe(inputs, exec);
}

void IndexRegistry::register_joint_recipe(const vector<IndexName>& identifiers,
                                          const vector<IndexName>& input_identifiers,
                                          const JointRecipeFunc& exec) {
    // We're going to generate a bunch of single-index recipes where the first
    // one to run calls the joint recipe, and other ones to run just return
    // their slice of the joint recipe's return value.
    
    // We need all the joint recipe names, one for each identifier we generate
    vector<RecipeName> names;
    
    // We need a place to hold the return values we can carry around by value.
    shared_ptr<vector<vector<string>>> results(std::make_shared<vector<vector<string>>>());
    
    for (size_t i = 0; i < identifiers.size(); i++) {
        IndexName being_generated = identifiers[i];
        
        // Create a recipe that invokes the joint recipe.
        RecipeFunc stub = [i, results, being_generated, &exec](const vector<const IndexFile*>& inputs, const IndexingPlan* plan, const IndexName& constructing) -> vector<string> {
            if (results->empty()) {
                // Invoke the actual logic, passing along the plan, and fill in results
                *results = exec(inputs, plan);
                // TODO: handle parallel invocations?
            }
            
            // Get our slice of the result file list.
            return results->at(i);
        };
        
        names.push_back(register_recipe(being_generated, input_identifiers, stub));
    }
    
    // Remember that these are a joint recipe.
    simplifications.emplace_back(input_identifiers, names);
}

IndexFile* IndexRegistry::get_index(const IndexName& identifier) {
    return registry.at(identifier).get();
}

const IndexFile* IndexRegistry::get_index(const IndexName& identifier) const {
    return registry.at(identifier).get();
}

string IndexRegistry::get_work_dir() {
    if (work_dir.empty()) {
        // Ensure the directory exists
        work_dir = temp_file::create_directory();
    }
    return work_dir;
}

vector<IndexName> IndexRegistry::dependency_order() const {
    
#ifdef debug_index_registry
    cerr << "finding topological order in dependency graph" << endl;
#endif
    
    // assign each index file an index in a vector (arbitrarily)
    map<IndexName, size_t> graph_idx;
    vector<IndexName> graph_label;
    for (const auto& idx_file : registry) {
        graph_idx[idx_file.first] = graph_label.size();
        graph_label.push_back(idx_file.first);
    }
    
    // build the dependency graph
    vector<vector<size_t>> dependency_graph(graph_label.size());
    for (size_t i = 0; i < dependency_graph.size(); ++i) {
        auto index = get_index(graph_label[i]);
        for (const auto& recipe : index->get_recipes()) {
            for (auto input : recipe.inputs) {
                dependency_graph[graph_idx[input->get_identifier()]].push_back(i);
            }
        }
    }
    
    // deduplicate any edges
    for (auto& adj : dependency_graph) {
        sort(adj.begin(), adj.end());
        adj.resize(unique(adj.begin(), adj.end()) - adj.begin());
    }
    
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
        exit(1);
    }
    
    // convert to return format
    vector<IndexName> ordered_identifiers(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        ordered_identifiers[i] = graph_label[order[i]];
    }
    
#ifdef debug_index_registry
    for (const auto& identifier : ordered_identifiers) {
        cerr << "\t" << to_string(identifier) << endl;
    }
#endif
    
    return ordered_identifiers;
}

IndexingPlan IndexRegistry::make_plan(const vector<IndexName>& end_products) const {
    
#ifdef debug_index_registry
    cerr << "generating plan for indexes:" << endl;
    for (const auto& product : end_products) {
        cerr << "\t" << to_string(product) << endl;
    }
#endif
    
    // get the dependency ordering of the indexes
    vector<IndexName> identifier_order = dependency_order();
    map<IndexName, size_t> dep_order_of_identifier;
    for (size_t i = 0; i < identifier_order.size(); ++i) {
        dep_order_of_identifier[identifier_order[i]] = i;
    }
    
    // TODO: I'm sure there's a more elegant implementation of this algorithm
    set<RecipeName> plan_elements;
    for (const auto& product : end_products) {
#ifdef debug_index_registry
        cerr << "making a plan for end product " << to_string(product) << endl;
#endif
        
        // records of (identifier, lowest level requester, ordinal index of recipe selected)
        vector<tuple<size_t, size_t, size_t>> plan_path;
        
        // map dependency priority to lowest level priority that requested this and
        // the number of requesters
        map<size_t, pair<size_t, size_t>, greater<size_t>> queue;
        queue[dep_order_of_identifier[product]] = pair<size_t, size_t>(identifier_order.size(), 1);
        
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
            
            if (get_index(identifier_order[get<0>(plan_path.back())])->is_finished()) {
                // this index has been provided, we don't need to use a recipe
#ifdef debug_index_registry
                cerr << "file has been provided as input" << endl;
#endif
                continue;
            }
            else if (!get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().empty()) {
                
                // this index can be created by a recipe
                const auto& recipe = get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().front();
#ifdef debug_index_registry
                cerr << "index can be made by a recipe requiring";
                for (auto input : recipe.inputs) {
                    cerr << " " << to_string(input->get_identifier());
                }
                cerr << endl;
#endif
                
                for (auto input : recipe.inputs) {
                    size_t dep_order = dep_order_of_identifier[input->get_identifier()];
                    auto f = queue.find(dep_order);
                    if (f == queue.end()) {
                        // no lower-level index has requested this one yet
                        queue[dep_order] = pair<size_t, size_t>(get<0>(plan_path.back()), 1);
                    }
                    else {
                        // record that one more index is requesting this one
                        f->second.second++;
                    }
                    
                }
            }
            else {
                // we've reached a file that needs to be provided but we don't have it,
                // so now we backtrack until hitting something that has remaining
                // lower priority recipes
#ifdef debug_index_registry
                cerr << "file cannot be made from existing inputs" << endl;
#endif
                while (!plan_path.empty() &&
                       get<2>(plan_path.back()) == get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().size()) {
                    // there are no remaining recipes to build the last index in the plan
                    
                    // remove items off the plan path until we get to the index that requested
                    // this one
                    size_t requester = get<1>(plan_path.back());
#ifdef debug_index_registry
                    cerr << "pruning path to previous requester: " << (requester == identifier_order.size() ? "PLAN TARGET" : to_string(identifier_order[requester])) << endl;
#endif
                    while (!plan_path.empty() && get<0>(plan_path.back()) != requester) {
                        
                        auto index = get_index(identifier_order[get<0>(plan_path.back())]);
                        if (!index->is_finished() && get<2>(plan_path.back()) < index->get_recipes().size()) {
                            // this index was using a recipe, we need to update its dependencies
                            // that are currently in the queue
                            const auto& recipe = index->get_recipes().at(get<2>(plan_path.back()));
                            for (auto input : recipe.inputs) {
                                size_t input_dep_order = dep_order_of_identifier[input->get_identifier()];
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
                            }
                        }
                        
                        plan_path.pop_back();
                    }
                    
                    if (!plan_path.empty()) {
                        // the requester should now use its next highest priority recipe
                        auto index = get_index(identifier_order[get<0>(plan_path.back())]);
                        if (!index->is_finished() && get<2>(plan_path.back()) < index->get_recipes().size()) {
                            // this index was using a recipe, we need to update its dependencies
                            // that are currently in the queue
                            const auto& recipe = index->get_recipes().at(get<2>(plan_path.back()));
                            for (auto input : recipe.inputs) {
                                size_t input_dep_order = dep_order_of_identifier[input->get_identifier()];
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
                            }
                        }
                        ++get<2>(plan_path.back());
                    }
                }
                
                if (!plan_path.empty()) {
                    const auto& recipe = get_index(identifier_order[get<0>(plan_path.back())])->get_recipes()[get<2>(plan_path.back())];
                    
#ifdef debug_index_registry
                    cerr << "advancing to recipe " << get<2>(plan_path.back()) << " for index " << to_string(identifier_order[get<0>(plan_path.back())]) << ", which requires";
                    for (auto input : recipe.inputs) {
                        cerr << " " << to_string(input->get_identifier());
                    }
                    cerr << endl;
#endif
                    for (auto input : recipe.inputs) {
                        size_t dep_order = dep_order_of_identifier[input->get_identifier()];
                        auto f = queue.find(dep_order);
                        if (f == queue.end()) {
                            // no lower-level index has requested this one yet
                            queue[dep_order] = pair<size_t, size_t>(get<0>(plan_path.back()), 1);
                        }
                        else {
                            // record that one more index is requesting this one
                            f->second.second++;
                        }
                        
                    }
                }
            }
            
        }
        
#ifdef debug_index_registry
        cerr << "final plan path for index " << to_string(product) << ":" << endl;
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
    
    // Now simplify the plan by using joint recipes if possible.
    
    // First we need to know all the indexes being created, and when
    map<IndexName, size_t> make_at_step;
    for (size_t i = 0; i < plan.steps.size(); i++) {
        make_at_step.emplace(plan.steps[i].first, i);
    }
    
    // We also need to know what steps we've already simplified, or are inputs.
    // We don't want to apply overlapping simplifications.
    vector<bool> fixed_step(plan.steps.size(), false);
    
    for (size_t i = 0; i < plan.steps.size(); i++) {
        auto& recipe = plan.steps.at(i);
        if (get_index(recipe.first)->is_finished()) {
            // This is already provided and ineligible for simplification.
            fixed_step[i] = true;
        }
    }
    
    for (auto& simplification : simplifications) {
        // For each set of output indexes from a simplification
        
#ifdef debug_index_registry
        cerr << "Consider simplification to jointly make:" << endl;
        for (auto& recipe: simplification.second) {
            cerr << "\t" << to_string(recipe.first) << endl;
        }
#endif
        
        // Determine if we are making all the products of the simplification,
        // and those products have not been involved in prior simplifications
        bool making_all_products_unsimplified = true;
        // And if so, the first step at which we are making any
        size_t first_step = numeric_limits<size_t>::max();
        for (auto& product_recipe : simplification.second) {
            const IndexName& product_name = product_recipe.first;
            
            auto found = make_at_step.find(product_name);
            if (found == make_at_step.end()) {
                // We aren't making this product
                
#ifdef debug_index_registry
                cerr << "We are not making " << to_string(product_name) << endl;
#endif
                
                making_all_products_unsimplified = false;
                break;
            }
            
            if (fixed_step[found->second]) {
                // We are making this product but we already simplified it or took it as input
                
#ifdef debug_index_registry
                cerr << "We cannot further simplify making " << to_string(product_name) << endl;
#endif
                
                making_all_products_unsimplified = false;
                break;
            }
            
#ifdef debug_index_registry
            cerr << "We are making " << to_string(product_name) << " at step " << found->second << endl;
#endif
            
            first_step = min(first_step, found->second);
        }
        
        if (!making_all_products_unsimplified) {
            // This simplification can't be used becuase it makes extra
            // products, or products that are already simplified.
            
#ifdef debug_index_registry
            cerr << "We are not making all the products for this simplification, or some products cannot be further simplified" << endl;
#endif
            
            continue;
        }
        
#ifdef debug_index_registry
        cerr << "To simplify, all inputs will need to be available before step " << first_step << endl;
#endif
        
        // See what we have available before the first step
        set<IndexName> available_in_time;
        for (size_t i = 0; i < first_step; i++) {
            available_in_time.insert(plan.steps[i].first);
        }
        
        // See if it's all the inputs the simplification needs
        bool all_available = true;
        for (auto& needed : simplification.first) {
            if (!available_in_time.count(needed)) {
#ifdef debug_index_registry
                cerr << "We are not making " << to_string(needed) << " in time or at all." << endl;
#endif
                all_available = false;
                break;
            }
        }
        
        if (!all_available) {
            // This simplification can't be used because not all its inputs are available in time.
            
#ifdef debug_index_registry
            cerr << "Not all inputs will be available in time." << endl;
#endif
            
            continue;
        }
        
#ifdef debug_index_registry
        cerr << "All inputs will be available in time. Apply simplification!" << endl;
#endif
        
        for (auto& recipe : simplification.second) {
            // Replace each relevant step with the corresponding joint step for that index.
            size_t step_to_simplify = make_at_step.at(recipe.first);
            plan.steps.at(step_to_simplify) = recipe;
            fixed_step[step_to_simplify] = true;
        }
    }
    
#ifdef debug_index_registry
    cerr << "plan after simplification:" << endl;
    for (auto plan_elem : plan.steps) {
        cerr << "\t" << to_string(plan_elem.first) << " " << plan_elem.second << endl;
    }
#endif

    // Now remove the input data from the plan
    plan.steps.resize(remove_if(plan.steps.begin(), plan.steps.end(), [&](const RecipeName& recipe_choice) {
        return get_index(recipe_choice.first)->is_finished();
    }) - plan.steps.begin());
    
    // The plan has methods that can come back and modify the registry.
    // We're not going to call any of them, but we have to hand off a non-const
    // pointer to ourselves so the plan can modify us later.
    plan.registry = const_cast<IndexRegistry*>(this);
    
    return plan;
}

string IndexRegistry::to_dot() const {
    return to_dot(vector<IndexName>());
}

string IndexRegistry::to_dot(const vector<IndexName>& targets) const {
    
    
    stringstream strm;
    strm << "digraph recipegraph {" << endl;
    
    set<IndexName> plan_targets(targets.begin(), targets.end());
    set<RecipeName> plan_elements;
    set<IndexName> plan_indexes;
    if (!targets.empty()) {
        IndexingPlan plan;
        try {
            plan = make_plan(targets);
        }
        catch (InsufficientInputException ex) {
            strm << "labelloc=\"t\";" << endl;
            strm << "label=\"Insufficient input to create targets\";" << endl;
        }
        for (const auto& plan_elem : plan.steps) {
            plan_elements.insert(plan_elem);
            plan_indexes.insert(plan_elem.first);
        }
    }
    
    map<IndexName, string> index_to_dot_id;
    size_t index_idx = 0;
    for (const auto& index_file : registry) {
        index_to_dot_id[index_file.first] = "I" + to_string(index_idx);
        ++index_idx;
        strm << index_to_dot_id[index_file.first] << "[label=\"" << to_string(index_file.first) << "\" shape=box";
        if (index_file.second->is_finished()) {
            strm << " style=\"filled,bold\" fillcolor=lightgray";
        }
        else if (plan_targets.count(index_file.first)) {
            strm << " style=\"filled,bold\" fillcolor=lightblue";
        }
        else if (plan_indexes.count(index_file.first)) {
            strm << " style=bold";
        }
        strm << "];" << endl;
    }
    string unselected_col = targets.empty() ? "black" : "gray33";
    size_t recipe_idx = 0;
    for (const auto& index_file : registry) {
        const auto& recipes = index_file.second->get_recipes();
        for (size_t priority_idx = 0; priority_idx < recipes.size(); ++priority_idx, ++recipe_idx) {
            const auto& recipe = recipes[priority_idx];
            string recipe_dot_id = "R" + to_string(recipe_idx);
            if (plan_elements.count(make_pair(index_file.first, priority_idx))) {
                strm << recipe_dot_id << "[label=\"" << priority_idx << "\" shape=circle style=bold];" << endl;
                strm << recipe_dot_id << " -> " << index_to_dot_id[index_file.first] << "[style=bold];" << endl;
            }
            else {
                strm << recipe_dot_id << "[label=\"" << priority_idx << "\" shape=circle];" << endl;
                strm << recipe_dot_id << " -> " << index_to_dot_id[index_file.first] << " [color=" << unselected_col << "];" << endl;
            }
            for (const auto& input : recipe.inputs) {
                if (plan_elements.count(make_pair(index_file.first, priority_idx))) {
                    strm << index_to_dot_id[input->get_identifier()] << " -> " << recipe_dot_id << "[style=bold];" << endl;
                }
                else {
                    strm << index_to_dot_id[input->get_identifier()] << " -> " << recipe_dot_id << " [color=" << unselected_col << "];" << endl;
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

const vector<IndexRecipe>& IndexFile::get_recipes() const {
    return recipes;
}

void IndexFile::provide(const vector<string>& filenames) {
    // append all filenames
    // TODO: would it be better to sometimes error check that the file isn't a duplicate?
    for (const string& filename : filenames) {
        this->filenames.emplace_back(filename);
    }
    provided_directly = true;
}

bool IndexFile::was_provided_directly() const {
    return provided_directly;
}

void IndexFile::execute_recipe(size_t recipe_priority, const IndexingPlan* plan) {
    assert(recipe_priority < recipes.size());
    auto& recipe = recipes[recipe_priority];
    for (auto input : recipe.inputs) {
        assert(input->is_finished());
    }
    filenames = recipe.execute(plan, get_identifier());
}

RecipeName IndexFile::add_recipe(const vector<const IndexFile*>& inputs,
                                 const RecipeFunc& exec) {
    recipes.emplace_back(inputs, exec);
    return RecipeName(identifier, recipes.size() - 1);
}

IndexRecipe::IndexRecipe(const vector<const IndexFile*>& inputs,
                         const RecipeFunc& exec) :
    exec(exec), inputs(inputs)
{
    // nothing more to do
}

vector<string> IndexRecipe::execute(const IndexingPlan* plan, const IndexName& constructing) {
    return exec(inputs, plan, constructing);
}

InsufficientInputException::InsufficientInputException(const IndexName& target,
                                                       const IndexRegistry& registry) :
    runtime_error("Insufficient input to create " + to_string(target)), target(target), inputs(registry.completed_indexes())
{
    // nothing else to do
}

const char* InsufficientInputException::what() const throw () {
    stringstream ss;
    ss << "Inputs" << endl;
    for (const auto& input : inputs) {
        ss << "\t" << to_string(input) << endl;
    }
    ss << "are insufficient to create target index " << to_string(target) << endl;
    return ss.str().c_str();
}

}

