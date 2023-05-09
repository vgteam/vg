// index_registry.cpp: index registry system implementation

#include "index_registry.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <thread>
#include <mutex>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <omp.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <xg.hpp>
#include <gbwt/variants.h>
#include <gbwtgraph/index.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/path_cover.h>
#include <gbwtgraph/gfa.h>
#include <vg/io/vpkg.hpp>
#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>
#include <Fasta.h>
#include <htslib/tbx.h>

#include "vg.hpp"
#include "vg_set.hpp"
#include "handle.hpp"
#include "utility.hpp"
#include "constructor.hpp"
#include "hash_map.hpp"
#include "haplotype_indexer.hpp"
#include "phase_unfolder.hpp"
#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "gcsa_helper.hpp"
#include "flat_file_back_translation.hpp"
#include "kmer.hpp"
#include "transcriptome.hpp"
#include "integrated_snarl_finder.hpp"
#include "snarl_distance_index.hpp"
#include "gfa.hpp"
#include "job_schedule.hpp"
#include "path.hpp"

#include "io/save_handle_graph.hpp"

#include "algorithms/gfa_to_handle.hpp"
#include "algorithms/prune.hpp"
#include "algorithms/component.hpp"
#include "algorithms/find_translation.hpp"

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


IndexingParameters::MutableGraphImplementation IndexingParameters::mut_graph_impl = PackedGraph;
int IndexingParameters::max_node_size = 32;
int IndexingParameters::pruning_max_node_degree = 128;
int IndexingParameters::pruning_walk_length = 24;
int IndexingParameters::pruning_max_edge_count = 3;
int IndexingParameters::pruning_min_component_size = 33;
double IndexingParameters::pruning_walk_length_increase_factor = 1.5;
double IndexingParameters::pruning_max_node_degree_decrease_factor = 0.75;
int IndexingParameters::gcsa_initial_kmer_length = gcsa::Key::MAX_LENGTH;
int IndexingParameters::gcsa_doubling_steps = gcsa::ConstructionParameters::DOUBLING_STEPS;
int64_t IndexingParameters::gcsa_size_limit = 2ll * 1024ll * 1024ll * 1024ll * 1024ll;
int64_t IndexingParameters::gbwt_insert_batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
int IndexingParameters::gbwt_insert_batch_size_increase_factor = 10;
int IndexingParameters::gbwt_sampling_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
bool IndexingParameters::bidirectional_haplo_tx_gbwt = false;
string IndexingParameters::gff_feature_name = "exon";
string IndexingParameters::gff_transcript_tag = "transcript_id";
bool IndexingParameters::use_bounded_syncmers = false;
int IndexingParameters::minimizer_k = 29;
int IndexingParameters::minimizer_w = 11;
int IndexingParameters::minimizer_s = 18;
int IndexingParameters::path_cover_depth = gbwtgraph::PATH_COVER_DEFAULT_N;
int IndexingParameters::giraffe_gbwt_downsample = gbwtgraph::LOCAL_HAPLOTYPES_DEFAULT_N;
int IndexingParameters::downsample_threshold = 3;
int IndexingParameters::downsample_context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
double IndexingParameters::max_memory_proportion = 0.75;
double IndexingParameters::thread_chunk_inflation_factor = 2.0;
IndexingParameters::Verbosity IndexingParameters::verbosity = IndexingParameters::Basic;

void copy_file(const string& from_fp, const string& to_fp) {
    ifstream from_file(from_fp, std::ios::binary);
    ofstream to_file(to_fp, std::ios::binary);
    if (!from_file) {
        cerr << "error:[IndexRegistry] Couldn't open input file " << from_fp << endl;
        exit(1);
    }
    if (!to_file) {
        cerr << "error:[IndexRegistry] Couldn't open output file " << to_fp << endl;
        exit(1);
    }
    to_file << from_file.rdbuf();
}

// return file size in bytes
int64_t get_file_size(const string& filename) {
    // get the file size
    ifstream infile(filename);
    infile.seekg(0, ios::end);
    return infile.tellg();
}

bool is_gzipped(const string& filename) {
    if (filename.size() > 2 && filename.substr(filename.size() - 3, 3) == ".gz") {
        return true;
    }
    return false;
}

int64_t get_num_samples(const string& vcf_filename) {
    htsFile* vcf_file = hts_open(vcf_filename.c_str(),"rb");
    if (!vcf_file) {
        cerr << "error:[IndexRegistry]: Failed to open VCF file: " << vcf_filename << endl;
        exit(1);
    }
    bcf_hdr_t* header = bcf_hdr_read(vcf_file);
    int64_t num_samples = bcf_hdr_nsamples(header);
    bcf_hdr_destroy(header);
    hts_close(vcf_file);
    return num_samples;
}

// quickly guess the number of variants in a VCF file based on the filesize
double approx_num_vars(const string& vcf_filename) {
    
    int64_t num_samples = get_num_samples(vcf_filename);
    int64_t file_size = get_file_size(vcf_filename);
    
    // TODO: bcf coefficient
    // a shitty regression that Jordan fit on the human autosomes, gives a very rough
    // estimate of the number of variants contained in a VCF
    if (is_gzipped(vcf_filename)) {
        // square root got a pretty good fit, for whatever reason
        return 0.255293 * file_size / sqrt(std::max(num_samples, (int64_t) 1));
    }
    else {
        return 0.192505 * file_size / std::max(num_samples, (int64_t) 1);
    }
}

// the ratio to the HashGraph memory usage, as estimated by a couple of graphs
// that Jordan had laying around when he wrote this
// in case it's useful later: XG is ~0.231
double format_multiplier() {
    switch (IndexingParameters::mut_graph_impl) {
        case IndexingParameters::HashGraph:
            return 1.0;
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

// approximate the memory of a graph that would be constructed with FASTAs and VCFs
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

vector<int64_t> each_approx_graph_memory(const vector<string>& fasta_filenames,
                                         const vector<string>& vcf_filenames) {
    
    auto n = max(fasta_filenames.size(), vcf_filenames.size());
    assert(fasta_filenames.size() == 1 || fasta_filenames.size() == n);
    assert(vcf_filenames.empty() || vcf_filenames.size() == 1 || vcf_filenames.size() == n);
    
    double total_ref_size = 0;
    vector<double> ref_sizes(fasta_filenames.size());
    for (int64_t i = 0; i < ref_sizes.size(); ++i) {
        double ref_size = get_file_size(fasta_filenames[i]);
        ref_sizes[i] = ref_size;
        total_ref_size += ref_size;
    }
    double total_var_count = 0;
    vector<int64_t> var_counts(vcf_filenames.size());
    for (int64_t i = 0; i < vcf_filenames.size(); ++i) {
        int64_t var_count = approx_num_vars(vcf_filenames[i]);
        var_counts[i] = var_count;
        total_var_count += var_count;
    }
    
    vector<int64_t> approx_memories(n);
    for (int64_t i = 0; i < n; ++i) {
        
        double ref_size, var_count;
        if (vcf_filenames.empty()) {
            var_count = 0;
        }
        else if (vcf_filenames.size() == 1) {
            var_count = total_var_count * (ref_sizes[i] / total_ref_size);
        }
        else {
            var_count = var_counts[i];
        }
        if (fasta_filenames.size() == 1 && total_var_count != 0.0) {
            ref_size = total_ref_size * (var_counts[i] / total_var_count);
        }
        else if (fasta_filenames.size() == 1) {
            ref_size = total_ref_size;
        }
        else {
            ref_size = ref_sizes[i];
        }
        
        // TODO: repetitive with previous function, magic constants
        double linear_memory = 30.4483 * ref_size;
        double var_memory = 2242.90 * var_count;
        double hash_graph_memory_usage = linear_memory + var_memory;
        approx_memories[i] = hash_graph_memory_usage * format_multiplier();
    }
    return approx_memories;
}

int64_t approx_graph_memory(const string& fasta_filename, const string& vcf_filename) {
    return approx_graph_memory(vector<string>(1, fasta_filename), vector<string>(1, vcf_filename));
}

// estimate the amount of memory of a GFA constructed graph
int64_t approx_graph_memory(const string& gfa_filename) {
    int64_t hash_graph_memory_usage = 13.17 * get_file_size(gfa_filename);
    return hash_graph_memory_usage * format_multiplier();}

int64_t approx_gbwt_memory(const string& vcf_filename) {
    return 21.9724 * log(std::max(get_num_samples(vcf_filename), (int64_t) 1)) * approx_num_vars(vcf_filename);
}

int64_t approx_graph_load_memory(const string& graph_filename) {
    // TODO: separate regressions for different graph types
    // this one was done on hash graphs, which probably have a larger expansion
    int64_t hash_graph_memory_usage = 12.52059 * get_file_size(graph_filename);
    return hash_graph_memory_usage * format_multiplier();
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

// return all of the contigs with variants in a VCF by iterating through
// the whole damn thing (SQ lines are not required, unfortunately)
vector<string> vcf_contigs(const string& filename) {
    
    htsFile* vcf = hts_open(filename.c_str(),"rb");
    if (vcf == nullptr) {
        cerr << "error:[IndexRegistry] Could not open VCF" << filename << endl;
    }
    
    bcf_hdr_t* header = bcf_hdr_read(vcf);
    unordered_set<string> contigs;
    bcf1_t* bcf_record = bcf_init();
    while (bcf_read(vcf, header, bcf_record) >= 0) {
        
        const char* chrom = bcf_hdr_id2name(header, bcf_record->rid);
        
        contigs.emplace(chrom);
    }
    bcf_destroy(bcf_record);
    vector<string> return_val(contigs.begin(), contigs.end());
    
    
    bcf_hdr_destroy(header);
    hts_close(vcf);
    
    sort(return_val.begin(), return_val.end());
    return return_val;
}

/*********************
 * Indexing helper functions
 ***********************/

// These can't be local lambdas in our indexer setup function because they
// would go away when the setup function returns.

static void init_in(ifstream& in, const string& name) {
    in.open(name);
    if (!in) {
        cerr << "error:[IndexRegistry] could not open input file '" << name << "'" << endl;
        exit(1);
    }
}

static void init_out(ofstream& out, const string& name) {
    out.open(name);
    if (!out) {
        cerr << "error:[IndexRegistry] could not write output to '" << name << "'" << endl;
        exit(1);
    }
}

static void init_in_out(fstream& strm, const string& name) {
    strm.open(name);
    if (!strm) {
        cerr << "error:[IndexRegistry] could not open '" << name << "'" << endl;
        exit(1);
    }
}

static auto init_mutable_graph() -> unique_ptr<MutablePathDeletableHandleGraph> {
    unique_ptr<MutablePathDeletableHandleGraph> graph;
    switch (IndexingParameters::mut_graph_impl) {
        case IndexingParameters::HashGraph:
            graph = make_unique<bdsg::HashGraph>();
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
}

// execute a function in another process and return true if successful
// REMEMBER TO SAVE ANY INDEXES CONSTRUCTED TO DISK WHILE STILL INSIDE THE LAMBDA!!
bool execute_in_fork(const function<void(void)>& exec) {
    
    // we have to clear out the pool of waiting OMP threads (if any) so that they won't
    // be copied with the fork and create deadlocks/races
    omp_pause_resource_all(omp_pause_soft);
    
    pid_t pid = fork();
    
    if (pid == -1) {
        cerr << "error:[IndexRegistry] failed to fork process" << endl;
        exit(1);
    }
    else if (pid == 0) {
        // this is the child process that will actually make the indexes
        
        // we want the pre-existing temp files to live beyond when this process exits
        temp_file::forget();
        
        exec();
                
        // end the child process successfully
        exit(0);
    } else {
        // This is the parent
        if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
            cerr << "[IndexRegistry]: Forked into child process with PID " << pid << "." << endl;
        }
    }
    
    // allow the child to finish
    int child_stat;
    waitpid(pid, &child_stat, 0); // 0 waits until the process fully exits
    
    // pass through signal-based exits
    if (WIFSIGNALED(child_stat)) {
        cerr << "error:[IndexRegistry] Child process " << pid << " signaled with status " << child_stat << " representing signal " << WTERMSIG(child_stat) << endl;
        if (raise(WTERMSIG(child_stat)) == 0) {
            // TODO: on Mac, raise isn't guaranteed to not return before the handler if it succeeds.
            // Also the signal might not be one that necessarily kills us.
            exit(1);
        } else {
            // We couldn't send ourselves the signal.
            exit(1);
        }
    }
    
    assert(WIFEXITED(child_stat));
    
    if (WEXITSTATUS(child_stat) != 0) {
        cerr << "warning:[IndexRegistry] Child process " << pid << " failed with status " << child_stat << " representing exit code " << WEXITSTATUS(child_stat) << endl;
        return false;
    }
    
    return true;
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
    registry.register_index("Reference GFA w/ Haplotypes", "haplo.gfa");
    registry.register_index("GTF/GFF", "gff");
    registry.register_index("Haplotype GTF/GFF", "haplo.gff");
    
    /// Chunked inputs
    registry.register_index("Chunked Reference FASTA", "chunked.fasta");
    registry.register_index("Chunked VCF", "chunked.vcf.gz");
    registry.register_index("Chunked VCF w/ Phasing", "phased.chunked.vcf.gz");
    registry.register_index("Chunked GTF/GFF", "chunked.gff");
    
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
    registry.register_index("NamedNodeBackTranslation", "segments.tsv");
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
    registry.register_index("Giraffe GBWT", "giraffe.gbwt");
    
    registry.register_index("Spliced Snarls", "spliced.snarls");
    
    registry.register_index("Giraffe Distance Index", "dist");
    registry.register_index("Spliced Distance Index", "spliced.dist");
    
    registry.register_index("GBWTGraph", "gg");
    registry.register_index("GBZ", "gbz");
    registry.register_index("Giraffe GBZ", "giraffe.gbz");
    
    registry.register_index("Minimizers", "min");
    
    /*********************
     * Register all recipes
     ***********************/
     
    // Note that recipes MAY NOT CAPTURE ANYTHING BY REFERENCE ([&]) from this scope!
    // This scope is on the stack and will go away by the time the recipes actually run!
    
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
    registry.register_recipe({"Chunked VCF"}, {"Chunked VCF w/ Phasing"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        alias_graph.register_alias(*constructing.begin(), inputs[0]);
        return vector<vector<string>>(1, inputs.front()->get_filenames());
    });
    
    ////////////////////////////////////
    // GFA Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering GFA recipes" << endl;
#endif
    
    // alias a phased GFA as an unphased one
    registry.register_recipe({"Reference GFA"}, {"Reference GFA w/ Haplotypes"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        alias_graph.register_alias(*constructing.begin(), inputs[0]);
        return vector<vector<string>>(1, inputs.front()->get_filenames());
    });
    
    ////////////////////////////////////
    // Chunking Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering chunking recipes" << endl;
#endif
    
    // meta recipe for with/out phasing and with/out transcripts
    auto chunk_contigs = [](const vector<const IndexFile*>& inputs,
                            const IndexingPlan* plan,
                            AliasGraph& alias_graph,
                            const IndexGroup& constructing,
                            bool has_gff,
                            bool phased_vcf) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Chunking inputs for parallelism." << endl;
        }
                        
        // boilerplate
        assert(inputs.size() == 1 || inputs.size() == 2 || inputs.size() == 3);
        assert(constructing.size() == inputs.size());
        vector<string> fasta_filenames, vcf_filenames, tx_filenames;
        bool has_vcf = inputs.size() == 3 || (inputs.size() == 2 && !has_gff);
        {
            int i = 0;
            if (has_gff) {
                tx_filenames = inputs[i++]->get_filenames();
            }
            fasta_filenames = inputs[i++]->get_filenames();
            if (has_vcf) {
                vcf_filenames = inputs[i++]->get_filenames();
            }
        }
        vector<vector<string>> all_outputs(constructing.size());
        string output_fasta, output_vcf, output_tx;
        {
            auto it = constructing.begin();
            if (has_gff) {
                output_tx = *it;
                ++it;
            }
            output_fasta = *it;
            ++it;
            if (has_vcf) {
                output_vcf = *it;
            }
        }
        auto& output_fasta_names = all_outputs[has_gff ? 1 : 0];
        
#ifdef debug_index_registry_recipes
        cerr << "chunking with vcf? " << has_vcf << ", with gff? " << has_gff << endl;
#endif
        
        // let's do this first, since it can detect data problems
        
        // i really hate to do two whole passes over the VCFs, but it's hard to see how to not
        // do this to be able to distinguish when we need to wait for a contig to become available
        // and when it simply doesn't have any variants in the VCFs (especially since it seems
        // to be not uncommon for the header to have these sequences included, such as alt scaffolds)
        vector<vector<string>> vcf_contigs_with_variants(vcf_filenames.size());
        vector<set<string>> vcf_samples(vcf_filenames.size());
        {
            // do the bigger jobs first to reduce makespan
            
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < vcf_filenames.size(); ++i) {
                                
                tbx_t* tabix_index = nullptr;
                for (string tabix_name : {vcf_filenames[i] + ".tbi", vcf_filenames[i] + ".csi"}) {
                    struct stat stat_tbi, stat_vcf;
                    if (stat(tabix_name.c_str(), &stat_tbi) != 0) {
                        // the tabix doesn't exist
                        continue;
                    }
                    stat(vcf_filenames[i].c_str(), &stat_vcf);
                    if (stat_vcf.st_mtime > stat_tbi.st_mtime) {
                        cerr << "warning:[IndexRegistry] Tabix index " + tabix_name + " is older than VCF " + vcf_filenames[i] + " and will not be used. Consider recreating this tabix index to speed up index creation.\n";
                        continue;
                    }
                    
                    tabix_index = tbx_index_load(tabix_name.c_str());
                    if (tabix_index == nullptr) {
                        cerr << "error:[IndexRegistry] failed to load tabix index " << tabix_index << endl;
                        exit(1);
                    }
                }
                
                // get the sample set
                htsFile* vcf = bcf_open(vcf_filenames[i].c_str(), "r");
                bcf_hdr_t* header = bcf_hdr_read(vcf);
                for (int j = 0; j < bcf_hdr_nsamples(header); ++j) {
                    vcf_samples[i].insert(header->samples[j]);
                }
                
                if (tabix_index != nullptr) {
                    // we have a tabix index, so we can make contigs query more efficiently
                    int num_seq_names;
                    const char** seq_names = tbx_seqnames(tabix_index, &num_seq_names);
                    for (int j = 0; j < num_seq_names; ++j) {
                        vcf_contigs_with_variants[i].push_back(seq_names[j]);
                    }
                    free(seq_names);
                    bcf_hdr_destroy(header);
                    int close_err_code = hts_close(vcf);
                    if (close_err_code != 0) {
                        cerr << "error:[IndexRegistry] encountered error closing VCF " << vcf_filenames[i] << endl;
                        exit(1);
                    }
                    continue;
                }
                
                // no tabix index, so we have to do the full scan
                
                bcf1_t* vcf_rec = bcf_init();
                string curr_contig;
                int err_code = bcf_read(vcf, header, vcf_rec);
                while (err_code == 0) {
                    const char* chrom = bcf_hdr_id2name(header, vcf_rec->rid);
                    if (!curr_contig.empty()) {
                        if (curr_contig < chrom) {
                            curr_contig = chrom;
                            vcf_contigs_with_variants[i].push_back(chrom);
                        }
                        else if (curr_contig > chrom) {
                            cerr << "error:[IndexRegistry] Contigs in VCF must be in ASCII-lexicographic order. Encountered contig '" << chrom << "' after contig '" << curr_contig << "' in VCF file" << vcf_filenames[i] << "." << endl;
                            exit(1);
                        }
                    }
                    else {
                        curr_contig = chrom;
                        vcf_contigs_with_variants[i].push_back(chrom);
                    }
                    
                    err_code = bcf_read(vcf, header, vcf_rec);
                }
                if (err_code != -1) {
                    cerr << "error:[IndexRegistry] failed to read from VCF " << vcf_filenames[i] << endl;
                    exit(1);
                }
                // we'll be moving on to a different file, so we won't demand that these
                // be in order anymore
                bcf_destroy(vcf_rec);
                bcf_hdr_destroy(header);
                err_code = hts_close(vcf);
                if (err_code != 0) {
                    cerr << "error:[IndexRegistry] encountered error closing VCF " << vcf_filenames[i] << endl;
                    exit(1);
                }
            }
        }
        
#ifdef debug_index_registry_recipes
        cerr << "contigs that have variants in the VCFs:" << endl;
        for (auto& vcf_contigs : vcf_contigs_with_variants) {
            for (auto& contig : vcf_contigs) {
                cerr << "\t" << contig << endl;
            }
        }
#endif
        
        // consolidate this for easy look up later
        unordered_set<string> contigs_with_variants;
        for (auto& vcf_contigs : vcf_contigs_with_variants) {
            for (auto& contig : vcf_contigs) {
                contigs_with_variants.insert(contig);
            }
        }
        
        unordered_map<string, int64_t> seq_files;
        unordered_map<string, int64_t> seq_lengths;
        // records of (length, name)
        priority_queue<pair<int64_t, string>> seq_queue;
                
        for (int64_t i = 0; i < fasta_filenames.size(); ++i) {
            FastaReference ref;
            ref.open(fasta_filenames[i]);
            for (const auto& idx_entry : *ref.index) {
                seq_files[idx_entry.first] = i;
                seq_lengths[idx_entry.first] = idx_entry.second.length;
                seq_queue.emplace(idx_entry.second.length, idx_entry.first);
            }
        }
        
        // we'll partition sequences that have the same samples (chunking the the VCFs
        // ultimately requires that we do this)
        map<set<string>, vector<string>> sample_set_contigs;
        for (int i = 0; i < vcf_samples.size(); ++i) {
            auto& contigs = sample_set_contigs[vcf_samples[i]];
            for (auto& contig : vcf_contigs_with_variants[i]) {
                contigs.push_back(contig);
            }
        }
        // move these lists of contigs into more convenient data structures
        vector<vector<string>> contig_groups;
        unordered_map<string, int64_t> contig_to_group;
        for (auto it = sample_set_contigs.begin(); it != sample_set_contigs.end(); ++it) {
            for (const auto& contig : it->second) {
                if (contig_to_group.count(contig)) {
                    cerr << "error:[IndexRegistry] Contig " << contig << " is found in multiple VCFs with different samples" << endl;
                    exit(1);
                }
                contig_to_group[contig] = contig_groups.size();
            }
            contig_groups.emplace_back(move(it->second));
        }
        
#ifdef debug_index_registry_recipes
        cerr << "contigs by sample group" << endl;
        for (int i = 0; i < contig_groups.size(); ++i) {
            cerr << "group " << i << endl;
            for (auto contig : contig_groups[i]) {
                cerr << "\t" << contig << endl;
            }
        }
#endif
                
        // we'll greedily assign contigs to the smallest bucket (2-opt bin packing algorithm)
        // modified to ensure that we can bucket contigs with the same sample sets together
        
        // one of the threads gets used up to do scheduling after this initial chunking
        int num_threads = get_thread_count();
        // we'll let it go a bit larger so we can take advantage of dynamic scheduling
        int max_num_buckets = max<int>(contig_groups.size(),
                                       ceil(IndexingParameters::thread_chunk_inflation_factor * num_threads));
        int num_buckets = 0;
        int groups_without_bucket = contig_groups.size();
        size_t num_sample_groups = max<size_t>(contig_groups.size(), 1);
        // buckets of contigs, grouped by sample groups
        vector<vector<vector<string>>> sample_group_buckets(num_sample_groups);
        // records of (total length, bucket index), grouped by sample gorups
        vector<priority_queue<pair<int64_t, int64_t>, vector<pair<int64_t, int64_t>>,
                              greater<pair<int64_t, int64_t>>>> bucket_queues(num_sample_groups);
        while (!seq_queue.empty()) {
            int64_t length;
            string seq_name;
            tie(length, seq_name) = seq_queue.top();
            seq_queue.pop();
            
            int64_t group = 0;
            if (contig_to_group.count(seq_name)) {
                // this contig has variant samples associated with it, so we need to
                // group it in with them
                group = contig_to_group[seq_name];
            }
            else {
                // this contig has no variants associated, we can put it in whichever
                // bucket is smallest
                for (int64_t i = 0; i < bucket_queues.size(); ++i) {
                    int64_t min_bucket_length = numeric_limits<int64_t>::max();
                    if (bucket_queues[i].empty()) {
                        min_bucket_length = 0;
                        group = i;
                    }
                    else if (bucket_queues[i].top().first < min_bucket_length) {
                        min_bucket_length = bucket_queues[i].top().first;
                        group = i;
                    }
                    
                }
            }
            
            // always make sure there's enough room in our budget of buckets
            // to make one for each sample group
            if (bucket_queues[group].empty() || num_buckets < max_num_buckets - groups_without_bucket) {
                // make a new bucket
                if (bucket_queues[group].empty()) {
                    groups_without_bucket--;
                }
                auto& group_buckets = sample_group_buckets[group];
                bucket_queues[group].emplace(seq_lengths[seq_name], group_buckets.size());
                group_buckets.emplace_back(1, seq_name);
                num_buckets++;
            }
            else {
                // add to the smallest bucket
                int64_t total_length, b_idx;
                tie(total_length, b_idx) = bucket_queues[group].top();
                bucket_queues[group].pop();
                sample_group_buckets[group][b_idx].emplace_back(seq_name);
                bucket_queues[group].emplace(total_length + length, b_idx);
            }
        }
        
        // merge the list of sample group buckets and collate with their
        // FASTA of origin
        vector<vector<pair<string, int64_t>>> buckets;
        for (auto& group_buckets : sample_group_buckets) {
            for (auto& bucket : group_buckets) {
                buckets.emplace_back();
                auto& new_bucket = buckets.back();
                for (auto& contig : bucket) {
                    new_bucket.emplace_back(move(contig), seq_files[contig]);
                }
            }
        }
        
        // sort the buckets in descending order by sequence length so that the
        // biggest jobs get dynamically scheduled first
        sort(buckets.begin(), buckets.end(),
             [&](const vector<pair<string, int64_t>>& a, const vector<pair<string, int64_t>>& b) {
            size_t len_a = 0, len_b = 0;
            for (const auto& contig : a) {
                len_a += seq_lengths.at(contig.first);
            }
            for (const auto& contig : b) {
                len_b += seq_lengths.at(contig.first);
            }
            return len_a > len_b;
        });
                
        // to look bucket index up from a contig
        unordered_map<string, int64_t> contig_to_idx;
        for (int64_t i = 0; i < buckets.size(); ++i) {
            for (auto& bucket_item : buckets[i])  {
                contig_to_idx[bucket_item.first] = i;
            }
        }
        
        // sort contigs of each bucket lexicographically so that they occur in the order
        // we'll discover them in the VCF
        for (auto& bucket : buckets) {
            sort(bucket.begin(), bucket.end());
        }
        
#ifdef debug_index_registry_recipes
        cerr << "assigned contigs into buckets:" << endl;
        for (int i = 0; i < buckets.size(); ++i) {
            cerr << "bucket " << i << endl;
            for (auto& contig : buckets[i]) {
                cerr << "\t" << contig.first << ": " << seq_lengths[contig.first] << endl;
            }
        }
#endif
        
        
        if (buckets.size() == fasta_filenames.size()
            && (!has_vcf || buckets.size() == vcf_filenames.size())
            && (!has_gff || buckets.size() == tx_filenames.size())) {
            // it looks like we might have just recapitulated the original chunking, let's check to make sure
            
            // does each bucket come from exactly one FASTA file?
            bool all_buckets_match = true;
            for (int64_t i = 0; i < fasta_filenames.size() && all_buckets_match; ++i) {
                FastaReference ref;
                ref.open(fasta_filenames[i]);
                int64_t bucket_idx = -1;
                for (const auto& idx_entry : *ref.index) {
                    if (bucket_idx == -1) {
                        bucket_idx = contig_to_idx.at(idx_entry.second.name);
                    }
                    else if (contig_to_idx.at(idx_entry.second.name) != bucket_idx) {
                        all_buckets_match = false;
                        break;
                    }
                }
            }
            
            // TODO: shouldn't I also check the correspondence on the input GTFs/VCFs?
            
            if (all_buckets_match) {
                // there's no need for chunking, just alias them
#ifdef debug_index_registry_recipes
                cerr << "chunking matches input files, no need to re-chunk" << endl;
#endif
                
                if (has_gff) {
                    all_outputs[0] = tx_filenames;
                    alias_graph.register_alias(output_tx, inputs[0]);
                }
                
                output_fasta_names = fasta_filenames;
                alias_graph.register_alias(output_fasta, inputs[has_gff ? 1 : 0]);
                
                if (has_vcf) {
                    all_outputs[all_outputs.size() - 1] = vcf_filenames;
                    alias_graph.register_alias(output_vcf, inputs[inputs.size() - 1]);
                }
                return all_outputs;
            }
        }
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Chunking FASTA(s)." << endl;
        }
        
        output_fasta_names.resize(buckets.size());
        if (has_vcf) {
            all_outputs[all_outputs.size() - 1].resize(buckets.size());
        }
        
        // make FASTA sequences for each bucket
        // the threading here gets to be pretty simple because the fai allows random access
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t i = 0; i < buckets.size(); ++i) {
            
            auto chunk_fasta_name = plan->output_filepath(output_fasta, i, buckets.size());
            auto chunk_fai_name = chunk_fasta_name + ".fai";
            output_fasta_names[i] = chunk_fasta_name;
                        
            ofstream outfile_fasta, outfile_fai;
            init_out(outfile_fasta, chunk_fasta_name);
            init_out(outfile_fai, chunk_fai_name);
            for (auto& assigned_seq : buckets[i]) {
                string contig = assigned_seq.first;
                int64_t ref_idx = assigned_seq.second;
                FastaReference ref;
                ref.open(fasta_filenames[ref_idx]);
                
                auto entry = ref.index->entry(contig);
                int64_t length = entry.length;
                int64_t line_length = entry.line_blen; // the base length
                
                // copy over the FASTA sequence
                outfile_fasta << '>' << contig << '\n';
                int64_t seq_start = outfile_fasta.tellp();
                int64_t j = 0;
                while (j < length) {
                    int64_t end = min<int64_t>(j + line_length, length);
                    outfile_fasta << ref.getSubSequence(contig, j, end - j) << '\n';
                    j = end;
                }
                
                // add an FAI entry
                outfile_fai << contig << '\t' <<  length << '\t' << seq_start << '\t' << line_length << '\t' << line_length + 1 << endl;
            }
        }
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Chunking VCF(s)." << endl;
        }
        
        // open all of the input VCF files
        vector<tuple<htsFile*, bcf_hdr_t*, bcf1_t*>> input_vcf_files(vcf_filenames.size());
        vector<atomic<bool>> input_checked_out_or_finished(vcf_filenames.size());
        for (int64_t i = 0; i < input_vcf_files.size(); ++i) {
            htsFile* vcf = bcf_open(vcf_filenames[i].c_str(), "r");
            if (!vcf) {
                cerr << "error:[IndexRegistry] failed to open VCF " << vcf_filenames[i] << endl;
                exit(1);
            }
            bcf_hdr_t* header = bcf_hdr_read(vcf);
            bcf1_t* vcf_rec = bcf_init();
            int err_code = bcf_read(vcf, header, vcf_rec);
            if (err_code == -1) {
                // this vcf is empty, actually
                input_checked_out_or_finished[i].store(true);
            }
            else if (err_code < 0) {
                cerr << "error:[IndexRegistry] failed to read VCF " << vcf_filenames[i] << endl;
                exit(1);
            }
            input_vcf_files[i] = make_tuple(vcf, header, vcf_rec);
        }
        
        unordered_map<string, int64_t> contig_to_vcf_idx;
        for (int64_t i = 0; i < vcf_contigs_with_variants.size(); ++i) {
            for (const auto& contig : vcf_contigs_with_variants[i]) {
                contig_to_vcf_idx[contig] = i;
            }
        }
        
        if (has_vcf) {
            
            auto& output_vcf_names = all_outputs.back();
            
            // see if we can identify any chunked VCFs that are identical to our input VCFs
            
            // records of (input vcf index, bucket index)
            vector<pair<int64_t, int64_t>> copiable_vcfs;
            for (int64_t i = 0; i < buckets.size(); ++i) {
                int64_t prev_vcf_idx = -1;
                int64_t count = 0;
                for (auto& contig : buckets[i]) {
                    if (contig_to_vcf_idx.count(contig.first)) {
                        int64_t vcf_idx = contig_to_vcf_idx[contig.first];
                        if (prev_vcf_idx == -1 || prev_vcf_idx == vcf_idx) {
                            prev_vcf_idx = vcf_idx;
                            count++;
                        }
                        else {
                            // we've seen a second input VCF, mark a sentinel and stop looking
                            count = -1;
                            break;
                        }
                    }
                }
                if (prev_vcf_idx >= 0 && count == vcf_contigs_with_variants[prev_vcf_idx].size()) {
                    // we saw all and only contigs from one VCF, we can just copy it
                    copiable_vcfs.emplace_back(prev_vcf_idx, i);
                }
            }
            
#ifdef debug_index_registry_recipes
            cerr << "identified " << copiable_vcfs.size() << " copiable VCFs:" << endl;
            for (const auto& copiable_vcf : copiable_vcfs) {
                cerr << "\tinput " << copiable_vcf.first << " " << vcf_filenames[copiable_vcf.first] << " -> bucket " << copiable_vcf.second << endl;
            }
#endif
            output_vcf_names.resize(buckets.size());
            
            // check if we can do a sort-of-aliasing for VCFs, since they're the most time-
            // consuming part of the chunking
            if (copiable_vcfs.size() == vcf_filenames.size()) {
                // all of the input VCFs could be copied to 1 bucket, so we'll just alias
                // them and make dummies for the rest
                for (auto vcf_copy : copiable_vcfs) {
                    int64_t input_idx, output_idx;
                    tie(input_idx, output_idx) = vcf_copy;
                    output_vcf_names[output_idx] = vcf_filenames[input_idx];
                }
                for (int64_t i = 0; i < output_vcf_names.size(); ++i) {
                    if (output_vcf_names[i].empty()) {
                        // this bucket didn't receive a VCF chunk, let's make a dummy VCF for it
                        auto output_vcf_name = plan->output_filepath(output_vcf, i, buckets.size());
                        htsFile* vcf = bcf_open(output_vcf_name.c_str(), "wz");
                        bcf_hdr_t* header = bcf_hdr_init("w");
                        // this is to satisfy HaplotypeIndexer, which doesn't like sample-less VCFs
                        if (phased_vcf) {
                            int sample_add_code = bcf_hdr_add_sample(header, "dummy");
                            if (sample_add_code != 0) {
                                cerr << "error:[IndexRegistry] error initializing VCF header" << endl;
                                exit(1);
                            }
                        }
                        int hdr_write_err_code = bcf_hdr_write(vcf, header);
                        if (hdr_write_err_code != 0) {
                            cerr << "error:[IndexRegistry] error writing VCF header to " << output_vcf_name << endl;
                            exit(1);
                        }
                        bcf_hdr_destroy(header);
                        int close_err_code = hts_close(vcf);
                        if (close_err_code != 0) {
                            cerr << "error:[IndexRegistry] encountered error closing VCF " << output_vcf_name << endl;
                            exit(1);
                        }
                        output_vcf_names[i] = output_vcf_name;
                    }
                }
                // register that this is an alias
                alias_graph.register_alias(output_vcf, inputs[inputs.size() - 1]);
#ifdef debug_index_registry_recipes
                cerr << "pseudo-aliased VCFs with filenames:" << endl;
                for (const auto& filename : output_vcf_names) {
                    cerr << "\t" << filename << endl;
                }
#endif
            }
            else {
                
                // trackers for whether we can write to a bucket's vcf
                vector<atomic<bool>> bucket_checked_out_or_finished(buckets.size());
                for (int64_t i = 0; i < buckets.size(); ++i) {
                    bucket_checked_out_or_finished[i].store(false);
                }
                
                // if we can copy over a vcf, we don't want to check it out for reading/writing
                for (auto copiable_vcf : copiable_vcfs) {
                    input_checked_out_or_finished[copiable_vcf.first].store(true);
                    bucket_checked_out_or_finished[copiable_vcf.second].store(true);
                }
                
#ifdef debug_index_registry_recipes
                cerr << "initializing chunked VCFs for output" << endl;
#endif
                
                // the output files
                vector<pair<htsFile*, bcf_hdr_t*>> bucket_vcfs(buckets.size());
                for (int64_t i = 0; i < buckets.size(); ++i) {
                    auto output_vcf_name = plan->output_filepath(output_vcf, i, buckets.size());
                    output_vcf_names[i] = output_vcf_name;
                    
                    if (bucket_checked_out_or_finished[i].load()) {
                        // we can copy to make this file, so we don't need to initialize
                        // a file
                        continue;
                    }
                    
                    // open to write in bgzipped format
                    htsFile* vcf_out = bcf_open(output_vcf_name.c_str(), "wz");
                    bcf_hdr_t* header_out = bcf_hdr_init("w");
                    
                    // identify which input VCFs we'll be pulling from
                    unordered_set<int64_t> vcf_indexes;
                    for (const auto& contig : buckets[i]) {
                        if (contig_to_vcf_idx.count(contig.first)) {
                            vcf_indexes.insert(contig_to_vcf_idx[contig.first]);
                        }
                    }
#ifdef debug_index_registry_recipes
                    cerr << "bucket " << i << " will add samples from input VCFs:" << endl;
                    for (auto j : vcf_indexes) {
                        cerr << "\t" << j << endl;
                    }
#endif
                    // merge will all the input headers
                    unordered_set<string> samples_added;
                    for (auto vcf_idx : vcf_indexes) {
                        
                        auto input_vcf_file = input_vcf_files[vcf_idx];
                        bcf_hdr_t* header_in = get<1>(input_vcf_file);
                        header_out = bcf_hdr_merge(header_out, header_in);
                        if (header_out == nullptr) {
                            cerr << "error:[IndexRegistry] error merging VCF header" << endl;
                            exit(1);
                        }
                        
                        // add the samples from every header
                        for (int64_t j = 0; j < bcf_hdr_nsamples(header_in); ++j) {
                            const char* sample = header_in->samples[j];
                            if (!samples_added.count(sample)) {
                                // TODO: the header has its own dictionary, so this shouldn't be necessary,
                                // but the khash_t isn't documented very well
                                samples_added.insert(sample);
                                // the sample hasn't been added yet
                                int sample_err_code = bcf_hdr_add_sample(header_out, header_in->samples[j]);
                                // returns a -1 if the sample is already included, which we expect
                                if (sample_err_code != 0) {
                                    cerr << "error:[IndexRegistry] error adding samples to VCF header" << endl;
                                    exit(1);
                                }
                            }
                        }
                    }
                    
                    // documentation in htslib/vcf.h says that this has to be called after adding samples
                    int sync_err_code = bcf_hdr_sync(header_out);
                    if (sync_err_code != 0) {
                        cerr << "error:[IndexRegistry] error syncing VCF header" << endl;
                        exit(1);
                    }
                    if (phased_vcf && bcf_hdr_nsamples(header_out) == 0) {
                        cerr << "warning:[IndexRegistry] VCF inputs from file(s)";
                        for (auto vcf_idx : vcf_indexes) {
                            cerr << " " << vcf_filenames[vcf_idx];
                        }
                        cerr << " have been identified as phased but contain no samples. Are these valid inputs?" << endl;
                        
                        // let's add a dummy so that HaplotypeIndexer doesn't get mad later
                        int sample_add_code = bcf_hdr_add_sample(header_out, "dummy");
                        if (sample_add_code != 0) {
                            cerr << "error:[IndexRegistry] error initializing VCF header" << endl;
                            exit(1);
                        }
                        // and re-sync, not sure if necessary, but it will be cheap regardless
                        sync_err_code = bcf_hdr_sync(header_out);
                        if (sync_err_code != 0) {
                            cerr << "error:[IndexRegistry] error syncing VCF header" << endl;
                            exit(1);
                        }
                    }
                    int hdr_write_err_code = bcf_hdr_write(vcf_out, header_out);
                    if (hdr_write_err_code != 0) {
                        cerr << "error:[IndexRegistry] error writing VCF header to " << output_vcf_name << endl;
                        exit(1);
                    }
                    
                    // remember these so that we can check them out later
                    bucket_vcfs[i] = make_pair(vcf_out, header_out);
                }
                
                // the parallel iteration in here is pretty complicated because contigs from
                // the input VCFs are being shuffled among the output bucket VCFs, and contigs
                // need to be both read and written in lexicographic order. the mutexes here
                // let the threads shift between reading and writing from different pairs of VCFs.
                // hopefully this high-contention process won't cause too many problems since
                // copying each contig takes up a relatively large amount of time
                
                // a mutex to lock the process of checking whether the next contig the thread
                // needs is exposed
                mutex input_vcf_mutex;
                // a mutex to lock the process of switching to a new bucket
                mutex output_vcf_mutex;
                // to keep track of which contig in the bucket we're looking for next
                vector<size_t> contig_idx(buckets.size(), 0);
                // how many buckets we've finished so far
                atomic<int64_t> buckets_finished(0);
                vector<thread> workers;
                for (int64_t i = 0; i < num_threads; ++i) {
                    // Worker must not capture i; it will be out of scope!
                    workers.emplace_back([&]() {
                        int64_t bucket_idx = -1;
                        while (buckets_finished.load() < buckets.size()) {
                            // check if any of the input VCFs need to be moved past a contig that isn't
                            // in our reference
                            input_vcf_mutex.lock();
                            int64_t contig_skip_idx = -1;
                            for (int64_t j = 0; j < input_vcf_files.size(); ++j) {
                                if (input_checked_out_or_finished[j].load()) {
                                    continue;
                                }
                                
                                const char* chrom = bcf_hdr_id2name(get<1>(input_vcf_files[j]),
                                                                    get<2>(input_vcf_files[j])->rid);
                                // check this index over the FASTA sequence lengths for the chromosome
                                if (!seq_lengths.count(chrom)) {
                                    contig_skip_idx = j;
                                    input_checked_out_or_finished[j].store(true);
                                }
                            }
                            input_vcf_mutex.unlock();
                            
                            if (contig_skip_idx != -1) {
                                // we found a contig in the VCF that isn't present in the FASTA, we'll have to skip it
                                
                                auto& input_vcf_file = input_vcf_files[contig_skip_idx];
                                string skip_contig = bcf_hdr_id2name(get<1>(input_vcf_file),
                                                                     get<2>(input_vcf_file)->rid);
                                cerr << "warning:[IndexRegistry] Skipping contig " + skip_contig + ", which is found in VCF(s) but not reference.\n";
                                
                                
                                // keep reading until end of file or a different contig
                                int read_err_code = 0;
                                while (read_err_code >= 0) {
                                    string contig = bcf_hdr_id2name(get<1>(input_vcf_file),
                                                                    get<2>(input_vcf_file)->rid);
                                    if (contig != skip_contig) {
                                        break;
                                    }
                                    
                                    read_err_code = bcf_read(get<0>(input_vcf_file), get<1>(input_vcf_file), get<2>(input_vcf_file));
                                }
                                
                                // check the input back out unless we've finished it
                                if (read_err_code >= 0) {
                                    input_checked_out_or_finished[contig_skip_idx].store(false);
                                }
                                continue;
                            }
                            
                            // select an output VCF corresponding to a bucket
                            int64_t copy_from_idx = -1, copy_to_idx = -1;
                            bool found_bucket = false;
                            output_vcf_mutex.lock();
                            if (!copiable_vcfs.empty()) {
                                // there are copiable VCFs remaining, do these first
                                tie(copy_from_idx, copy_to_idx) = copiable_vcfs.back();
                                copiable_vcfs.pop_back();
                            }
                            else {
                                // start iteration at 1 so we always advance to a new bucket if possible
                                for (int64_t j = 1; j <= buckets.size(); ++j) {
                                    int64_t next_bucket_idx = (bucket_idx + j) % buckets.size();
                                    if (!bucket_checked_out_or_finished[next_bucket_idx].load()) {
                                        bucket_checked_out_or_finished[next_bucket_idx].store(true);
                                        bucket_idx = next_bucket_idx;
                                        found_bucket = true;
                                        break;
                                    }
                                }
                            }
                            output_vcf_mutex.unlock();
                            
                            if (copy_from_idx >= 0) {
#ifdef debug_index_registry_recipes
                                cerr << "direct copying " + vcf_filenames[copy_from_idx] + " to " + output_vcf_names[copy_to_idx] + "\n";
#endif
                                // we can copy an entire file on this iteration instead of parsing
                                copy_file(vcf_filenames[copy_from_idx], output_vcf_names[copy_to_idx]);
                                if (file_exists(vcf_filenames[copy_from_idx] + ".tbi")) {
                                    // there's also a tabix, grab that as well
                                    copy_file(vcf_filenames[copy_from_idx] + ".tbi", output_vcf_names[copy_to_idx] + ".tbi");
                                }
                                // this bucket is now totally finished
                                buckets_finished.fetch_add(1);
                                continue;
                            }
                            
                            if (!found_bucket) {
                                // it's now possible for all buckets to be checked out simultaneously
                                // by other threads, so there's no more need to have this thread running
#ifdef debug_index_registry_recipes
                                cerr << "thread exiting\n";
#endif
                                return;
                            }
                            
                            auto& ctg_idx = contig_idx[bucket_idx];
                            
                            if (!contigs_with_variants.count(buckets[bucket_idx][ctg_idx].first)) {
                                // this contig doesn't have variants in any of the VCFs, so we skip it
                                ++ctg_idx;
                                if (ctg_idx == buckets[bucket_idx].size()) {
                                    buckets_finished.fetch_add(1);
                                }
                                else {
                                    bucket_checked_out_or_finished[bucket_idx].store(false);
                                }
                                continue;
                            }
                            
                            htsFile* vcf_out = bucket_vcfs[bucket_idx].first;
                            bcf_hdr_t* header_out = bucket_vcfs[bucket_idx].second;
                            
                            // check if any of the VCFs' next contig is the next one we want for
                            // this bucket (and lock other threads out from checking simultaneously)
                            int64_t input_idx = -1;
                            input_vcf_mutex.lock();
                            for (int64_t j = 0; j < input_vcf_files.size(); ++j) {
                                if (input_checked_out_or_finished[j].load()) {
                                    continue;
                                }
                                
                                // what is the next contig in this VCF?
                                const char* chrom = bcf_hdr_id2name(get<1>(input_vcf_files[j]),
                                                                    get<2>(input_vcf_files[j])->rid);
                                if (buckets[bucket_idx][ctg_idx].first == chrom) {
                                    input_idx = j;
                                    input_checked_out_or_finished[j].store(true);
                                    break;
                                }
                            }
                            input_vcf_mutex.unlock();
                            
                            if (input_idx < 0) {
                                // other threads need to get through earlier contigs until this bucket's next
                                // contig is exposed
                                bucket_checked_out_or_finished[bucket_idx].store(false);
                                continue;
                            }
                            
                            // we've checked out one of the input vcfs, now we can read from it
                            auto& input_vcf_file = input_vcf_files[input_idx];
                            
                            int read_err_code = 0;
                            while (read_err_code >= 0) {
                                
                                const char* chrom = bcf_hdr_id2name(get<1>(input_vcf_file), get<2>(input_vcf_file)->rid);
                                if (buckets[bucket_idx][ctg_idx].first != chrom) {
                                    break;
                                }
                                
    // FIXME: i'm not sure how important it is to handle these malformed VCFs it is
    //                            // read the "END" info field to see if we need to repair it (this seems to be a problem
    //                            // in the grch38 liftover variants from 1kg)
    //                            int32_t* end_dst = NULL;
    //                            int num_end;
    //                            int end_err_code = bcf_get_info_int32(get<1>(input_vcf_file), get<2>(input_vcf_file), "END",
    //                                                                  &end_dst, &num_end);
    //                            if (end_err_code >= 0) {
    //                                // there is an END tag to read
    //                                int64_t end = *end_dst;
    //                                // note: we can query alleles without bcf_unpack, because it will have already
    //                                // unpacked up to info fields
    //                                // calculate it the way the spec says to
    //                                int64_t calc_end = get<2>(input_vcf_file)->pos + strlen(get<2>(input_vcf_file)->d.allele[0]) - 1;
    //                                if (end != calc_end) {
    //                                    string msg = "warning:[IndexRegistry] fixing \"END\" of variant " + buckets[bucket_idx][ctg_idx].first + " " + to_string(get<2>(input_vcf_file)->pos) + " from " + to_string(end) + " to " + to_string(calc_end) + "\n";
    //#pragma omp critical
    //                                    cerr << msg;
    //
    //                                    int update_err_code = bcf_update_info_int32(get<1>(input_vcf_file), get<2>(input_vcf_file), "END",
    //                                                                                &calc_end, 1);
    //                                    if (update_err_code < 0) {
    //                                        cerr << "error:[IndexRegistry] failed to update \"END\"" << endl;
    //                                        exit(1);
    //                                    }
    //                                }
    //                                free(end_dst);
    //                            }
                                
                                bcf_translate(header_out, get<1>(input_vcf_file), get<2>(input_vcf_file));
                                
                                int write_err_code = bcf_write(vcf_out, header_out, get<2>(input_vcf_file));
                                if (write_err_code != 0) {
                                    cerr << "error:[IndexRegistry] error writing VCF line to " << output_vcf_names[bucket_idx] << endl;
                                    exit(1);
                                }
                                
                                read_err_code = bcf_read(get<0>(input_vcf_file), get<1>(input_vcf_file), get<2>(input_vcf_file));
                            }
                            
                            if (read_err_code >= 0) {
                                // there's still more to read, it's just on different contigs
                                input_checked_out_or_finished[input_idx].store(false);
                            }
                            else if (read_err_code != -1) {
                                // we encountered a real error
                                cerr << "error:[IndexRegistry] error reading VCF file " << vcf_filenames[input_idx] << endl;
                                exit(1);
                            }
                            
                            // we finished this contig
                            ++ctg_idx;
                            if (ctg_idx == buckets[bucket_idx].size()) {
                                buckets_finished.fetch_add(1);
                            }
                            else {
                                bucket_checked_out_or_finished[bucket_idx].store(false);
                            }
                        }
                    });
                }
                
                // barrier sync
                for (auto& worker : workers) {
                    worker.join();
                }
                
                // close out files
                for (int64_t i = 0; i < input_vcf_files.size(); ++i) {
                    auto vcf_file = input_vcf_files[i];
                    bcf_destroy(get<2>(vcf_file));
                    bcf_hdr_destroy(get<1>(vcf_file));
                    int err_code = hts_close(get<0>(vcf_file));
                    if (err_code != 0) {
                        cerr << "error:[IndexRegistry] encountered error closing VCF " << vcf_filenames[i] << endl;
                        exit(1);
                    }
                }
                for (int64_t i = 0; i < bucket_vcfs.size(); ++i) {
                    if (!bucket_vcfs[i].second) {
                        // we didn't open this VCF (probably because we just copied it)
                        continue;
                    }
                    bcf_hdr_destroy(bucket_vcfs[i].second);
                    int close_err_code = hts_close(bucket_vcfs[i].first);
                    if (close_err_code != 0) {
                        cerr << "error:[IndexRegistry] encountered error closing VCF " << output_vcf_names[i] << endl;
                        exit(1);
                    }
                }
            }
            
            // TODO: move this into the same work queue as the rest of the VCF chunking?
            // tabix index
#pragma omp parallel for schedule(dynamic, 1)
            for (int64_t i = 0; i < buckets.size(); ++i) {
                // tabix-index the bgzipped VCF we just wrote
                
                if (file_exists(output_vcf_names[i] + ".tbi")) {
                    // the tabix already exists
                    continue;
                }
                
                // parameters inferred from tabix main's sourcecode
                int min_shift = 0;
                tbx_conf_t conf = tbx_conf_vcf;
                int tabix_err_code = tbx_index_build(output_vcf_names[i].c_str(), min_shift, &conf);
                if (tabix_err_code == -2) {
                    cerr << "error:[IndexRegistry] output VCF is not bgzipped: " << output_vcf_names[i] << endl;
                    exit(1);
                }
                else if (tabix_err_code != 0) {
                    cerr << "warning:[IndexRegistry] could not tabix index VCF " + output_vcf_names[i] + "\n";
                }
            }
        }
        
        if (has_gff) {
            
            if (IndexingParameters::verbosity != IndexingParameters::None) {
                cerr << "[IndexRegistry]: Chunking GTF/GFF(s)." << endl;
            }
            
            auto& output_gff_names = all_outputs[0];
            vector<ofstream> tx_files_out(buckets.size());
            for (int64_t i = 0; i < buckets.size(); ++i) {
                output_gff_names.emplace_back(plan->output_filepath(output_tx, i, buckets.size()));
                init_out(tx_files_out[i], output_gff_names.back());
            }
            
            // mutexes to lock the process of checking out for writing/reading
            mutex gff_out_mutex;
            
            // we'll thread by input files in this case
            vector<thread> tx_workers;
            vector<atomic<bool>> chunk_gff_checked_out(buckets.size());
            for (int64_t i = 0; i < chunk_gff_checked_out.size(); ++i) {
                chunk_gff_checked_out[i].store(false);
            }
            
            atomic<int64_t> input_gffs_read(0);
            for (int64_t i = 0; i < num_threads; ++i) {
                // Worker must not capture i; it will be out of scope! 
                tx_workers.emplace_back([&]() {
                    while (input_gffs_read.load() < tx_filenames.size()) {
                        
                        int64_t idx = input_gffs_read.fetch_add(1);
                        
                        if (idx >= tx_filenames.size()) {
                            break;
                        }

                        ifstream infile_tx;
                        init_in(infile_tx, tx_filenames[idx]);
                        
                        ofstream tx_chunk_out;
                        
                        int64_t prev_chunk_idx = -1;
                        while (infile_tx.good()) {
                            
                            string line;
                            getline(infile_tx, line);
                            
                            stringstream line_strm(line);
                            string chrom;
                            getline(line_strm, chrom, '\t');
                            if (chrom.empty() || chrom.front() == '#') {
                                // skip header
                                continue;
                            }
                            
                            auto it = contig_to_idx.find(chrom);
                            if (it == contig_to_idx.end()) {
                                cerr << "error:[IndexRegistry] contig " << chrom << " from GTF/GFF " << tx_filenames[idx] << " is not found in reference" << endl;
                                exit(1);
                            }
                            int64_t chunk_idx = it->second;
                            if (chunk_idx != prev_chunk_idx) {
                                // we're transitioning between chunks, so we need to check the chunk
                                // out for writing
                                
                                // release the old chunk
                                if (prev_chunk_idx >= 0) {
                                    tx_chunk_out.close();
                                    tx_chunk_out.clear();
                                    chunk_gff_checked_out[prev_chunk_idx].store(false);
                                }
                                
                                // keep trying to check the new chunk until succeeding
                                bool success = false;
                                while (!success) {
                                    // only one thread can try to check out at a time
                                    gff_out_mutex.lock();
                                    if (!chunk_gff_checked_out[chunk_idx].load()) {
                                        // the chunk is free to be written to
                                        chunk_gff_checked_out[chunk_idx].store(true);
                                        success = true;
                                    }
                                    gff_out_mutex.unlock();
                                    if (!success) {
                                        // wait for a couple seconds to check again if we can write
                                        // to the file
                                        this_thread::sleep_for(chrono::seconds(1));
                                    }
                                    else {
                                        // open for writing, starting from the end
                                        tx_chunk_out.open(output_gff_names[chunk_idx], ios_base::ate);
                                        if (!tx_chunk_out) {
                                            cerr << "error:[IndexRegistry] could not open " << output_gff_names[chunk_idx] << " for appending" << endl;
                                            exit(1);
                                        }
                                    }
                                }
                            }
                            
                            // copy the line to the chunk
                            tx_chunk_out << line << '\n';
                            
                            prev_chunk_idx = chunk_idx;
                        }
                        
                        // release the last chunk we were writing to
                        if (prev_chunk_idx >= 0) {
                            tx_chunk_out.close();
                            chunk_gff_checked_out[prev_chunk_idx].store(false);
                        }
                    }
                });
            }
            
            // barrier sync
            for (auto& worker : tx_workers) {
                worker.join();
            }
        }
        
        
        return all_outputs;
    };
    
    // call the meta recipe
    registry.register_recipe({"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF w/ Phasing"}, {"GTF/GFF", "Reference FASTA", "VCF w/ Phasing"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, true, true);
    });
    registry.register_recipe({"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF"}, {"GTF/GFF", "Reference FASTA", "VCF"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, true, false);
    });
    registry.register_recipe({"Chunked Reference FASTA", "Chunked VCF w/ Phasing"}, {"Reference FASTA", "VCF w/ Phasing"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, false, true);
    });
    registry.register_recipe({"Chunked Reference FASTA", "Chunked VCF"}, {"Reference FASTA", "VCF"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, false, false);
    });
    registry.register_recipe({"Chunked GTF/GFF", "Chunked Reference FASTA"}, {"GTF/GFF", "Reference FASTA"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, true, false);
    });
    registry.register_recipe({"Chunked Reference FASTA"}, {"Reference FASTA"},
                             [=](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return chunk_contigs(inputs, plan, alias_graph, constructing, false, false);
    });
    
    
    ////////////////////////////////////
    // VG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering VG recipes" << endl;
#endif
    
    // meta-recipe for removing variant paths
    auto strip_variant_paths = [](const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexGroup& constructing) {
        
        assert(inputs.size() == 1);
        assert(constructing.size() == 1);
        
        auto chunk_filenames = inputs.at(0)->get_filenames();
        auto output_index = *constructing.begin();
        vector<vector<string>> all_outputs(constructing.size());
        
        auto& output_names = all_outputs.front();
        output_names.resize(chunk_filenames.size());
        auto strip_chunk = [&](int64_t i) {
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
                if (Paths::is_alt(graph->get_path_name(path))) {
                    alt_paths.push_back(path);
                }
            });
            
            // delete them
            for (auto path : alt_paths) {
                graph->destroy_path(path);
            }
            
            // and save the graph
            vg::io::save_handle_graph(graph.get(), outfile);
            
            output_names[i] = output_name;
        };
        
        // approximate the time and memory use for each chunk
        vector<pair<int64_t, int64_t>> approx_job_requirements;
        for (auto& chunk_filename : chunk_filenames) {
            approx_job_requirements.emplace_back(get_file_size(chunk_filename),
                                                 approx_graph_load_memory(chunk_filename));
        }
        
        JobSchedule schedule(approx_job_requirements, strip_chunk);
        schedule.execute(plan->target_memory_usage());
        
        // return the filename(s)
        return all_outputs;
    };
    
    // strip alt allele paths from a graph that has them
    registry.register_recipe({"VG"}, {"VG w/ Variant Paths"},
                             [strip_variant_paths](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Stripping allele paths from VG." << endl;
        }
        
        return strip_variant_paths(inputs, plan, constructing);
    });
        
    // meta-recipe for creating a VG and its segment space from a GFA
    auto construct_from_gfa = [&](const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing VG graph from GFA input." << endl;
        }
        
        assert(constructing.size() == 3);
        vector<vector<string>> all_outputs(constructing.size());
        
        assert(inputs.size() == 1);
        assert(constructing.size() == 3);
        IndexName output_max_id = "MaxNodeID";
        assert(constructing.count(output_max_id));
        IndexName output_translation = "NamedNodeBackTranslation";
        assert(constructing.count(output_translation));
        IndexName output_index = "VG";
        assert(constructing.count(output_index));
        auto input_filenames = inputs.at(0)->get_filenames();
        if (input_filenames.size() > 1) {
            cerr << "error:[IndexRegistry] Graph construction does not support multiple GFAs at this time." << endl;
            exit(1);
        }
        auto input_filename = input_filenames.front();
        
        string output_name = plan->output_filepath(output_index);
        string translation_name = plan->output_filepath(output_translation);
        string max_id_name = plan->output_filepath(output_max_id);
        // The graph and max ID are streams we start here.
        ofstream outfile, max_id_outfile;
        init_out(outfile, output_name);
        init_out(max_id_outfile, max_id_name);
        auto graph = init_mutable_graph();
        
        // make the graph from GFA, and save segment info to the translation file if there is nontrivial segment info.
        try {
            algorithms::gfa_to_path_handle_graph(input_filename, graph.get(), numeric_limits<int64_t>::max(), translation_name);
        }
        catch (algorithms::GFAFormatError& e) {
            cerr << "error:[IndexRegistry] Input GFA is not usable in VG." << endl;
            cerr << e.what() << endl;
            exit(1);
        }
        
        // Now we need to append some splits to the output file.
        ofstream translation_outfile;
        translation_outfile.open(translation_name, std::ios_base::app);
        if (!translation_outfile) {
            cerr << "error:[IndexRegistry] could not append output to " << translation_name << endl;
            exit(1);
        }
        
        handlealgs::chop(*graph, IndexingParameters::max_node_size, [&](nid_t old_id, size_t offset, size_t rev_offset, handle_t new_node) {
#pragma omp critical (translation_outfile)
            {
                // Write each cut to a line in the translation file, after the segment names are defined.
                translation_outfile << "K\t" << old_id << "\t" << offset << "\t" << rev_offset << "\t" << graph->get_id(new_node) << std::endl;
            }
        });
        
        // save the graph
        vg::io::save_handle_graph(graph.get(), outfile);
        // and the max id
        max_id_outfile << graph->max_node_id();
        
        // return the filenames
        all_outputs[0].push_back(max_id_name);
        all_outputs[1].push_back(translation_name);
        all_outputs[2].push_back(output_name);
        return all_outputs;
    };
    
    
    
    // A meta-recipe to make VG and spliced VG files using the Constructor
    // Expects inputs to be ordered: FASTA, VCF[, GTF/GFF][, Insertion FASTA]
    auto construct_with_constructor = [](const vector<const IndexFile*>& inputs,
                                         const IndexingPlan* plan,
                                         const IndexGroup& constructing,
                                         bool alt_paths,
                                         bool has_transcripts,
                                         bool has_variants) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing";
            if (has_transcripts) {
                cerr << " spliced";
            }
            cerr << " VG graph from FASTA";
            if (has_variants) {
                cerr << " and VCF";
            }
            cerr << " input." << endl;
        }
        
        assert(constructing.size() == 2);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_max_id = *constructing.begin();
        auto output_graph = *constructing.rbegin();
        auto& max_id_names = all_outputs[0];
        auto& graph_names = all_outputs[1];
        
        bool has_ins_fasta = false;
        if (1 + int(has_transcripts) + int(has_variants) != inputs.size()) {
            assert(2 + int(has_transcripts) + int(has_variants) == inputs.size());
            has_ins_fasta = true;
        }
                
        // unpack the inputs
        vector<string> ref_filenames, vcf_filenames, insertions, transcripts;
        {
            size_t i = 0;
            if (has_transcripts) {
                transcripts = inputs[i++]->get_filenames();
            }
            ref_filenames = inputs[i++]->get_filenames();
            if (has_variants) {
                vcf_filenames = inputs[i++]->get_filenames();
            }
            if (has_ins_fasta) {
                insertions = inputs[i++]->get_filenames();
            }
        }
        
        if (has_ins_fasta) {
            if (insertions.size() > 1) {
                cerr << "error:[IndexRegistry] can only provide one FASTA for insertion sequences" << endl;
                exit(1);
            }
            
            // make sure this FASTA has an fai index before we get into all the parallel stuff
            FastaReference ins_ref;
            ins_ref.open(insertions.front());
        }
                
        if (has_variants && ref_filenames.size() != 1 && vcf_filenames.size() != 1 &&
            ref_filenames.size() != vcf_filenames.size()) {
            cerr << "[IndexRegistry]: When constructing graph from multiple FASTAs and multiple VCFs, the FASTAs and VCFs must be matched 1-to-1, but input contains " <<  ref_filenames.size() << " FASTA files and " << vcf_filenames.size() << " VCF files." << endl;
            exit(1);
        }
        if (has_transcripts && transcripts.size() != 1 && ref_filenames.size() != 1 &&
            transcripts.size() != ref_filenames.size()) {
            cerr << "[IndexRegistry]: When constructing graph from multiple GTF/GFFs and multiple FASTAs, the GTF/GFFs and the FASTAs must be matched 1-to-1, but input contains " <<  transcripts.size() << " GTF/GFF files and " <<  ref_filenames.size() << " FASTA files." << endl;
            exit(1);
        }
        if (has_transcripts && has_variants && transcripts.size() != 1 && vcf_filenames.size() != 1 &&
            transcripts.size() != vcf_filenames.size()) {
            cerr << "[IndexRegistry]: When constructing graph from multiple GTF/GFFs and multiple VCFs, the GTF/GFFs and the VCFs must be matched 1-to-1, but input contains " <<  transcripts.size() << " GTF/GFF files and " <<  vcf_filenames.size() << " VCF files." << endl;
            exit(1);
        }
        
        // are we broadcasting the transcripts from one chunk to many?
        bool broadcasting_txs = transcripts.size() != max(ref_filenames.size(),
                                                          vcf_filenames.size());
                
        // TODO: this estimate should include splice edges too
        vector<pair<int64_t, int64_t>> approx_job_requirements;
        {
            size_t i = 0;
            for (auto approx_mem : each_approx_graph_memory(ref_filenames, vcf_filenames)) {
                int64_t approx_time;
                if (!vcf_filenames.empty() && vcf_filenames.size() != 1) {
                    approx_time = get_file_size(vcf_filenames[i]);
                }
                else {
                    approx_time = get_file_size(ref_filenames[i]);
                }
                approx_job_requirements.emplace_back(approx_time, approx_mem);
                ++i;
            }
        }
        
#ifdef debug_index_registry_recipes
        cerr << "approximate chunk requirements:" << endl;
        for (size_t i = 0; i < approx_job_requirements.size(); ++i) {
            auto requirement = approx_job_requirements[i];
            cerr << "\tchunk " << i << " -- time: " << requirement.first << ", memory: " << requirement.second << endl;
        }
#endif
        graph_names.resize(max(ref_filenames.size(), vcf_filenames.size()));
        vector<pair<nid_t, nid_t>> node_id_ranges(graph_names.size());
        auto make_graph = [&](int64_t idx) {
#ifdef debug_index_registry_recipes
            cerr << "making graph chunk " << idx << endl;
#endif
            
            auto ref_filename = ref_filenames.size() == 1 ? ref_filenames[0] : ref_filenames[idx];
            
#ifdef debug_index_registry_recipes
            cerr << "constructing graph with Constructor for ref " << ref_filename << endl;
#endif
            
            // init and configure the constructor
            Constructor constructor;
            constructor.do_svs = true;
            constructor.alt_paths = alt_paths;
            constructor.max_node_size = IndexingParameters::max_node_size;
            constructor.show_progress = IndexingParameters::verbosity >= IndexingParameters::Debug;
            
            if (ref_filenames.size() != 1 && vcf_filenames.size() == 1) {
                // we have multiple FASTA but only 1 VCF, so we'll limit the
                // constructor to the contigs of this FASTA for this run
                FastaReference ref;
                ref.open(ref_filename);
                for (const string& seqname : ref.index->sequenceNames) {
                    constructor.allowed_vcf_names.insert(seqname);
                }
            }
            else if (!vcf_filenames.empty() && vcf_filenames.size() != 1 && ref_filenames.size() == 1) {
                // we have multiple VCFs but only 1 FASTA, so we'll limit the
                // constructor to the contigs of this VCF for this run
                
                // unfortunately there doesn't seem to be a good way to do this without
                // iterating over the entire file:
                for (const auto& contig : vcf_contigs(vcf_filenames[idx])) {
                    constructor.allowed_vcf_names.insert(contig);
                }
            }
            
            string output_name = plan->output_filepath(output_graph, idx,
                                                       max(ref_filenames.size(), vcf_filenames.size()));
            ofstream outfile;
            init_out(outfile, output_name);
            
            auto graph = init_mutable_graph();
            
            vector<string> fasta(1, ref_filename);
            vector<string> vcf;
            if (!vcf_filenames.empty()) {
                vcf.emplace_back(vcf_filenames.size() == 1 ? vcf_filenames[0] : vcf_filenames[idx]);
            }
            
            // do the construction
            constructor.construct_graph(fasta, vcf, insertions, graph.get());
            
#ifdef debug_index_registry_recipes
            cerr << "resulting graph has " << graph->get_node_count() << " nodes" << endl;
#endif
                               
            
            if (!transcripts.empty()) {
                
                auto transcript_filename = transcripts[transcripts.size() == 1 ? 0 : idx];
                
#ifdef debug_index_registry_recipes
                cerr << "adding transcripts from " << transcript_filename << endl;
#endif
                
                ifstream infile_tx;
                init_in(infile_tx, transcript_filename);
                
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
                transcriptome.feature_type = IndexingParameters::gff_feature_name;
                transcriptome.transcript_tag = IndexingParameters::gff_transcript_tag;
                
                // add the splice edges
                auto dummy = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
                size_t transcripts_added = transcriptome.add_reference_transcripts(vector<istream *>({&infile_tx}), dummy, false, false);
                
                if (broadcasting_txs && !path_names.empty() && transcripts_added == 0
                    && transcript_file_nonempty(transcripts[idx])) {
                    cerr << "warning:[IndexRegistry] no matching paths from transcript file " << transcript_filename << " were found in graph chunk containing the following paths:" << endl;
                    for (const string& path_name : path_names) {
                        cerr << "\t" << path_name << endl;
                    }
                }
                
                node_id_ranges[idx] = make_pair(transcriptome.graph().min_node_id(),
                                                transcriptome.graph().max_node_id());
                
                // save the file
                transcriptome.write_graph(&outfile);
            }
            else {
                
                node_id_ranges[idx] = make_pair(graph->min_node_id(), graph->max_node_id());
                
                // save the file
                vg::io::save_handle_graph(graph.get(), outfile);
            }
            
            graph_names[idx] = output_name;
        };
        
        // TODO: allow contig renaming through Constructor::add_name_mapping
        
        // construct the jobs in parallel, trying to use multithreading while also
        // restraining memory usage
        JobSchedule schedule(approx_job_requirements, make_graph);
        schedule.execute(plan->target_memory_usage());
        
        // merge the ID spaces if we need to
        vector<nid_t> id_increment(1, 1 - node_id_ranges[0].first);
        if (graph_names.size() > 1 || id_increment.front() != 0) {
            // the increments we'll need to make each ID range non-overlapping
            for (int i = 1; i < node_id_ranges.size(); ++i) {
                id_increment.push_back(node_id_ranges[i - 1].second + id_increment[i - 1] - node_id_ranges[i].first + 1);
            }
            
            vector<pair<int64_t, int64_t>> approx_job_requirements;
            for (int i = 0; i < node_id_ranges.size(); ++i) {
                approx_job_requirements.emplace_back(node_id_ranges[i].second - node_id_ranges[i].first,
                                                     approx_graph_load_memory(graph_names[i]));
            }
            
#ifdef debug_index_registry_recipes
            cerr << "computed node ID increments for chunks:" << endl;
            for (int i = 0; i < id_increment.size(); ++i) {
                cerr << "\t[" << node_id_ranges[i].first << ", " << node_id_ranges[i].second  << "] + " << id_increment[i] << endl;
            }
#endif
            
            // do the incrementation in parallel
            auto increment_node_ids = [&](int64_t idx) {
                
                // load the graph
                ifstream infile;
                init_in(infile, graph_names[idx]);
                unique_ptr<MutablePathMutableHandleGraph> graph
                    = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(infile);
                
                // adjust the IDs
                graph->increment_node_ids(id_increment[idx]);
                
                // save back to the same file
                ofstream outfile;
                init_out(outfile, graph_names[idx]);
                vg::io::save_handle_graph(graph.get(), outfile);
            };
            
            JobSchedule schedule(approx_job_requirements, increment_node_ids);
            schedule.execute(plan->target_memory_usage());
        }
        
        // save the max node id as a simple text file
        auto max_id_name = plan->output_filepath(output_max_id);
        ofstream max_id_outfile;
        init_out(max_id_outfile, max_id_name);
        nid_t max_node_id = node_id_ranges.back().second + id_increment.back();
        max_id_outfile << max_node_id;
        
        max_id_names.push_back(max_id_name);
        
        // return the filename(s)
        return all_outputs;
    };
    
    // the specific instantiations of the meta-recipe above
    registry.register_recipe({"MaxNodeID", "VG w/ Variant Paths"}, {"Chunked Reference FASTA", "Chunked VCF w/ Phasing", "Insertion Sequence FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, false, true);
    });
    registry.register_recipe({"MaxNodeID", "VG w/ Variant Paths"}, {"Chunked Reference FASTA", "Chunked VCF w/ Phasing"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, false, true);
    });
    registry.register_recipe({"MaxNodeID", "NamedNodeBackTranslation", "VG"}, {"Reference GFA"},
                             [&](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_from_gfa(inputs, plan, constructing);
    });
    registry.register_recipe({"MaxNodeID", "VG"}, {"Chunked Reference FASTA", "Chunked VCF", "Insertion Sequence FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, false, true);
    });
    registry.register_recipe({"MaxNodeID", "VG"}, {"Chunked Reference FASTA", "Chunked VCF"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                                          const IndexingPlan* plan,
                                                          AliasGraph& alias_graph,
                                                          const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, false, true);
    });
    registry.register_recipe({"MaxNodeID", "VG"}, {"Chunked Reference FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, false, false);
    });
    
#ifdef debug_index_registry_setup
    cerr << "registering Spliced VG recipes" << endl;
#endif
    
    ////////////////////////////////////
    // Spliced VG Recipes
    ////////////////////////////////////
    
    registry.register_recipe({"Spliced VG"}, {"Spliced VG w/ Variant Paths"},
                             [strip_variant_paths](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Stripping allele paths from spliced VG." << endl;
        }
        
        return strip_variant_paths(inputs, plan, constructing);
    });
    
    // TODO: spliced vg from GFA input
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Variant Paths"},
                             {"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF w/ Phasing", "Insertion Sequence FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                                          const IndexingPlan* plan,
                                                          AliasGraph& alias_graph,
                                                          const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, true, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Variant Paths"},
                             {"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF w/ Phasing"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                                          const IndexingPlan* plan,
                                                          AliasGraph& alias_graph,
                                                          const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, true, true, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG"},
                             {"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF", "Insertion Sequence FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, true, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG"},
                             {"Chunked GTF/GFF", "Chunked Reference FASTA", "Chunked VCF"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, true, true);
    });
    
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG"},
                             {"Chunked GTF/GFF", "Chunked Reference FASTA"},
                             [construct_with_constructor](const vector<const IndexFile*>& inputs,
                                                          const IndexingPlan* plan,
                                                          AliasGraph& alias_graph,
                                                          const IndexGroup& constructing) {
        return construct_with_constructor(inputs, plan, constructing, false, true, false);
    });
    
    
    ////////////////////////////////////
    // XG Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering XG recipes" << endl;
#endif
    
    // TODO: currently disabling this to ensure, but I'd prefer to make a separate
    // semantic XG for a node-chopped variety and handle the pipeline differences
    // with simplifications
    
//    registry.register_recipe({"XG"}, {"Reference GFA"},
//                             [](const vector<const IndexFile*>& inputs,
//                                const IndexingPlan* plan,
//                                AliasGraph& alias_graph,
//                                const IndexGroup& constructing) {
//        if (IndexingParameters::verbosity != IndexingParameters::None) {
//            cerr << "[IndexRegistry]: Constructing XG graph from GFA input." << endl;
//        }
//        assert(constructing.size() == 1);
//        vector<vector<string>> all_outputs(constructing.size());
//        auto output_index = *constructing.begin();
//        auto gfa_names = inputs.front()->get_filenames();
//        if (gfa_names.size() > 1) {
//            cerr << "error:[IndexRegistry] Graph construction does not support multiple GFAs at this time." << endl;
//            exit(1);
//        }
//
//        string output_name = plan->output_filepath(output_index);
//        ofstream outfile;
//        init_out(outfile, output_name);
//
//        xg::XG xg_index;
//        xg_index.from_gfa(gfa_names.front());
//
//        vg::io::save_handle_graph(&xg_index, outfile);
//
//        // return the filename
//        all_outputs[0].emplace_back(output_name);
//        return all_outputs;
//    });
    
    auto make_xg_from_graph = [](const vector<const IndexFile*>& inputs,
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
                             [make_xg_from_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing XG graph from VG graph." << endl;
        }
        return make_xg_from_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced XG"}, {"Spliced VG w/ Transcript Paths"},
                             [make_xg_from_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing spliced XG graph from spliced VG graph." << endl;
        }
        return make_xg_from_graph(inputs, plan, constructing);
    });
    
    ////////////////////////////////////
    // MaxNodeID Recipes
    ////////////////////////////////////
    
    ////////////////////////////////////
    // GBWT Recipes
    ////////////////////////////////////
    
#ifdef debug_index_registry_setup
    cerr << "registering GBWT recipes" << endl;
#endif
    
    // merge multiple GBWTs if there are multiple, otherwise leave in place
    auto merge_gbwts = [](const vector<string>& gbwt_names,
                           const IndexingPlan* plan,
                           const IndexName& constructing_name) {
        if (gbwt_names.size() > 1) {
            if (IndexingParameters::verbosity != IndexingParameters::None) {
                cerr << "[IndexRegistry]: Merging contig GBWTs." << endl;
            }
            // we also need to merge the GBWTs
            
            string merged_gbwt_name = plan->output_filepath(constructing_name);
            ofstream outfile;
            init_out(outfile, merged_gbwt_name);
            
            vector<gbwt::GBWT> gbwt_indexes(gbwt_names.size());
            for (size_t i = 0; i < gbwt_names.size(); ++i) {
                load_gbwt(gbwt_indexes[i], gbwt_names[i], IndexingParameters::verbosity >= IndexingParameters::Debug);
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
    
    // meta-recipe to make GBWTs
    auto make_gbwt = [merge_gbwts](const vector<const IndexFile*>& inputs,
                        bool include_named_paths,
                        const IndexingPlan* plan,
                        const IndexGroup& constructing) {
        
        assert(inputs.size() == 2);
        
        auto vcf_filenames = inputs[0]->get_filenames();
        auto graph_filenames = inputs[1]->get_filenames();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        
        auto output_index = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        if ((graph_filenames.size() != 1 && graph_filenames.size() != vcf_filenames.size()) ||
            (vcf_filenames.size() != 1 && graph_filenames.size() != vcf_filenames.size())) {
            cerr << "[IndexRegistry]: When constructing GBWT from multiple graphs and multiple VCFs, the graphs and VCFs must be matched 1-to-1, but input contains " <<  graph_filenames.size() << " graphs and " << vcf_filenames.size() << " VCF files." << endl;
            exit(1);
        }
        if (vcf_filenames.size() == 1 && graph_filenames.size() != 1) {
            // FIXME: it should at least try to join the graph chunks together
            cerr << "[IndexRegistry]: GBWT construction currently does not support broadcasting 1 VCF to multiple graph chunks." << endl;
            exit(1);
        }
        
        if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        
        int64_t target_memory_usage = plan->target_memory_usage();
        vector<pair<int64_t, int64_t>> approx_job_requirements;
        
        vector<string> gbwt_names(vcf_filenames.size());
        unique_ptr<PathHandleGraph> broadcast_graph;
        if (graph_filenames.size() == 1) {
            // we only have one graph, so we can save time by loading it only one time
            // test streams for I/O
            ifstream infile;
            init_in(infile, graph_filenames.front());
            
            // we don't want to double-count the graph's contribution to memory in separate jobs, so we
            // subtract it once from the target memory use
            target_memory_usage = max<int64_t>(0, target_memory_usage - approx_graph_load_memory(graph_filenames.front()));
            
            // estimate the time and memory requirements
            for (auto vcf_filename : vcf_filenames) {
                approx_job_requirements.emplace_back(get_file_size(vcf_filename), approx_gbwt_memory(vcf_filename));
            }
            
            // load the graph
            broadcast_graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            
        }
        else {
            // estimate the time and memory requirements
            for (int64_t i = 0; i < vcf_filenames.size(); ++i) {
                approx_job_requirements.emplace_back(get_file_size(vcf_filenames[i]),
                                                     approx_gbwt_memory(vcf_filenames[i]) + approx_graph_load_memory(graph_filenames[i]));
            }
            
        }
        
        // Prepare a single shared haplotype indexer, since everything on it is thread safe.
        // Make this critical so we don't end up with a race on the verbosity
        unique_ptr<HaplotypeIndexer> haplotype_indexer;
#pragma omp critical
        {
            haplotype_indexer = unique_ptr<HaplotypeIndexer>(new HaplotypeIndexer());
            // HaplotypeIndexer resets this in its constructor
            if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
                gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
            }
            else {
                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
            }
        }
        haplotype_indexer->show_progress = IndexingParameters::verbosity >= IndexingParameters::Debug;
        // from the toil-vg best practices
        haplotype_indexer->force_phasing = true;
        haplotype_indexer->discard_overlaps = true;
        
        // If we're using a single graph, we're going to need to do each VCF's
        // named paths in its job, and then come back and do the rest. So we
        // need to know which paths still need to be done (and whether there are any).
        // So we first make a set of all the paths that need doing, and then we
        // clear them out when they're done, and if any are left we know we
        // need a job to do them.
        unordered_set<path_handle_t> broadcast_graph_paths_to_do;
        if (include_named_paths && broadcast_graph) {
            broadcast_graph->for_each_path_handle([&](const path_handle_t& path_handle) {
                // Look at all the paths in advance
                if (broadcast_graph->is_empty(path_handle) || Paths::is_alt(broadcast_graph->get_path_name(path_handle))) {
                    // Skip empty paths and alt allele paths
                    return;
                }
                // Keep the rest.
                broadcast_graph_paths_to_do.insert(path_handle);
            });
        }
        
        // construct a GBWT from the i-th VCF
        auto gbwt_job = [&](size_t i) {
            string gbwt_name;
            if (vcf_filenames.size() != 1) {
                // multiple components, so make a temp file that we will merge later
                gbwt_name = temp_file::create();
            }
            else {
                // one component, so we will actually save the output
                gbwt_name = plan->output_filepath(output_index);
            }
            
            // load the contig graph if necessary
            unique_ptr<PathHandleGraph> contig_graph;
            if (graph_filenames.size() != 1) {
                ifstream infile;
                init_in(infile, graph_filenames[i]);
                contig_graph = vg::io::VPKG::load_one<PathHandleGraph>(infile);
            }
            
            auto graph = graph_filenames.size() == 1 ? broadcast_graph.get() : contig_graph.get();
            
            // Parse the VCFs for this job
            vector<string> parse_files = haplotype_indexer->parse_vcf(vcf_filenames[i],
                                                                      *graph);
            
            // Build the GBWT from the parse files and the graph.
            // For fast merging later, we need to ensure that all threads on a single contig end up in the same initial GBWT.
            // So, if we have just one graph, only threads visited by the VCF can go in.
            // Then at the end, if there are non-alt paths left over, we add another job to make a GBWT just of those paths.
            // Otherwise, if we have one graph per job, all threads from the graph can go in.
            unique_ptr<gbwt::DynamicGBWT> gbwt_index = haplotype_indexer->build_gbwt(parse_files, 
                                                                                     "GBWT" + std::to_string(i),
                                                                                     include_named_paths ? graph : nullptr,
                                                                                     nullptr,
                                                                                     include_named_paths && (bool)broadcast_graph);
            
            save_gbwt(*gbwt_index, gbwt_name, IndexingParameters::verbosity == IndexingParameters::Debug);
            
            gbwt_names[i] = gbwt_name;
            
            if (include_named_paths && broadcast_graph) {
                // We have to check off the paths we embeded in this job.
                for (size_t contig_number = 0; contig_number < gbwt_index->metadata.contigs(); contig_number++) {
                    // Go through all contig names in the metadata
                    string contig_name = gbwt_index->metadata.contig(contig_number);
                    
                    if (graph->has_path(contig_name)) {
                        // And get the graph path
                        path_handle_t contig_path = graph->get_path_handle(contig_name);
                        #pragma omp critical (broadcast_graph_paths_done)
                        {
                            // Check it off in a thread-safe way.
                            // TODO: Will this be too much locking and unlocking when we do transcripts?
                            broadcast_graph_paths_to_do.erase(contig_path);
                        }
                    }
                }
            }
        };
        
        {
            // Do all the GBWT jobs
            JobSchedule schedule(approx_job_requirements, gbwt_job);
            schedule.execute(target_memory_usage);
        }
        
        if (include_named_paths && broadcast_graph && !broadcast_graph_paths_to_do.empty()) {
            // We're Back for One Last Job.
            // We need to embed these remaining paths that weren't VCF contigs.
            
            // There's no VCF to load, so our memory estimate is just 0. The graph is loaded already.
            // TODO: can we improve this?
            approx_job_requirements.clear();
            approx_job_requirements.emplace_back(0, 0);
            
            // This job has exclusive use of our data structures.
            auto one_last_job = [&](size_t ignored) {
                // Make a temp file
                string gbwt_name = temp_file::create();
                
                // Make a GBWT of the remaining graph paths.
                unique_ptr<gbwt::DynamicGBWT> gbwt_index = haplotype_indexer->build_gbwt({}, 
                                                                                         "Leftovers",
                                                                                         broadcast_graph.get(),
                                                                                         &broadcast_graph_paths_to_do);
                
                // And save it in the temp file
                save_gbwt(*gbwt_index, gbwt_name, IndexingParameters::verbosity == IndexingParameters::Debug);
                
                // And add it as one final GBWT.
                gbwt_names.push_back(gbwt_name);
            };
            
            JobSchedule schedule(approx_job_requirements, one_last_job);
            schedule.execute(target_memory_usage);
        }
        
        // merge GBWTs if necessary
        output_names.push_back(merge_gbwts(gbwt_names, plan, output_index));
        // return filename
        return all_outputs;
    };
    
    registry.register_recipe({"GBWT"}, {"Chunked VCF w/ Phasing", "VG w/ Variant Paths"},
                             [make_gbwt](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing GBWT from VG graph and phased VCF input." << endl;
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        return make_gbwt(inputs, true, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GBWT"}, {"Chunked VCF w/ Phasing", "Spliced VG w/ Variant Paths"},
                             [make_gbwt](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing GBWT from spliced VG graph and phased VCF input." << endl;
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        // TODO: If we include named paths here, we then generate
        // haplotype-specific transcripts following them when making the
        // "Haplotype-Transcript GBWT". It's not clear that that's correct, and
        // the "Spliced GBWT" never feeds into a GBZ, so we leave them out for
        // now. 
        return make_gbwt(inputs, false, plan, constructing);
    });
    
    // Giraffe will prefer to use a downsampled haplotype GBWT if possible
    registry.register_recipe({"Giraffe GBWT"}, {"GBWT", "XG"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Downsampling full GBWT." << endl;
        }
        
        assert(inputs.size() == 2);
        auto gbwt_filenames = inputs[0]->get_filenames();
        auto xg_filenames = inputs[1]->get_filenames();
        assert(gbwt_filenames.size() == 1);
        assert(xg_filenames.size() == 1);
        auto gbwt_filename = gbwt_filenames.front();
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto& output_names = all_outputs[0];
        auto sampled_gbwt_output = *constructing.begin();
        
        ifstream infile_xg, infile_gbwt;
        init_in(infile_xg, xg_filename);
        init_in(infile_gbwt, gbwt_filename);

        auto output_name = plan->output_filepath(sampled_gbwt_output);
        ofstream outfile_sampled_gbwt;
        init_out(outfile_sampled_gbwt, output_name);
        
        auto xg_index = vg::io::VPKG::load_one<xg::XG>(infile_xg);
        auto gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(infile_gbwt);

        // Downsample only if it would reduce the number of haplotypes sufficiently.
        size_t threshold = IndexingParameters::giraffe_gbwt_downsample * IndexingParameters::downsample_threshold;
        bool downsample = (gbwt_index->hasMetadata() && gbwt_index->metadata.haplotypes() >= threshold);


        bool success;
        if (downsample) {
            // Downsample the haplotypes and generate a path cover of components without haplotypes.
            
            // We need to drop paths that are alt allele paths and not pass them
            // through from a graph that has them to the synthesized GBWT.
            std::function<bool(const path_handle_t&)> path_filter = [&xg_index](const path_handle_t& path) {
                return !Paths::is_alt(xg_index->get_path_name(path));
            };
            
            // clang wants this one cast to function first for some reason?
            function<void(void)> exec = [&]() {
                gbwt::GBWT cover = gbwtgraph::local_haplotypes(*xg_index, *gbwt_index,
                                                               IndexingParameters::giraffe_gbwt_downsample,
                                                               IndexingParameters::downsample_context_length,
                                                               IndexingParameters::gbwt_insert_batch_size,
                                                               IndexingParameters::gbwt_sampling_interval,
                                                               true, // Also include named paths from the graph
                                                               &path_filter,
                                                               IndexingParameters::verbosity >= IndexingParameters::Debug);
                save_gbwt(cover, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
            };
            success = execute_in_fork(exec);
        }
        else {
            // Augment the GBWT with a path cover of components without haplotypes.
            if (IndexingParameters::verbosity != IndexingParameters::None) {
                cerr << "[IndexRegistry]: Not enough haplotypes; augmenting the full GBWT instead." << endl;
            }
            
            success = execute_in_fork([&]() {
                gbwt::DynamicGBWT dynamic_index(*gbwt_index);
                gbwt_index.reset();
                gbwtgraph::augment_gbwt(*xg_index, dynamic_index,
                                        IndexingParameters::path_cover_depth,
                                        IndexingParameters::downsample_context_length,
                                        IndexingParameters::gbwt_insert_batch_size,
                                        IndexingParameters::gbwt_sampling_interval,
                                        IndexingParameters::verbosity >= IndexingParameters::Debug);
                gbwt::GBWT cover = gbwt::GBWT(dynamic_index);
                save_gbwt(cover, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
            });
        }
        
        if (!success) {
            IndexingParameters::gbwt_insert_batch_size *= IndexingParameters::gbwt_insert_batch_size_increase_factor;
            throw RewindPlanException("[IndexRegistry]: Exceeded GBWT insert buffer size, expanding and reattempting.", {"Giraffe GBWT"});
        }
        
        output_names.push_back(output_name);
        return all_outputs;
    });
    
    // do a greedy haplotype cover if we don't have haplotypes
    registry.register_recipe({"Giraffe GBWT"}, {"XG"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing a greedy path cover GBWT" << endl;
        }
        
        assert(inputs.size() == 1);
        auto xg_filenames = inputs[0]->get_filenames();
        assert(xg_filenames.size() == 1);
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto& output_names = all_outputs[0];
        auto path_cover_output = *constructing.begin();
        
        ifstream infile_xg;
        init_in(infile_xg, xg_filename);
        
        auto output_name = plan->output_filepath(path_cover_output);
        ofstream outfile_path_cover_gbwt;
        init_out(outfile_path_cover_gbwt, output_name);
        
        auto xg_index = vg::io::VPKG::load_one<xg::XG>(infile_xg);
        
        auto comp_sizes = algorithms::component_sizes(*xg_index);
        size_t max_comp_size = *max_element(comp_sizes.begin(), comp_sizes.end());
        
        // We need to drop paths that are alt allele paths and not pass them
        // through from a graph that has them to the synthesized GBWT.
        std::function<bool(const path_handle_t&)> path_filter = [&xg_index](const path_handle_t& path) {
            return !Paths::is_alt(xg_index->get_path_name(path));
        };
        
        // make a GBWT from a greedy path cover
        bool success = execute_in_fork([&]() {
            gbwt::GBWT cover = gbwtgraph::path_cover_gbwt(*xg_index,
                                                          IndexingParameters::path_cover_depth,
                                                          IndexingParameters::downsample_context_length,
                                                          std::max<gbwt::size_type>(IndexingParameters::gbwt_insert_batch_size, 20 * max_comp_size), // buffer size recommendation from Jouni
                                                          IndexingParameters::gbwt_sampling_interval,
                                                          true, // Also include named paths from the graph
                                                          &path_filter,
                                                          IndexingParameters::verbosity >= IndexingParameters::Debug);
            
            save_gbwt(cover, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
        });
        
        if (!success) {
            IndexingParameters::gbwt_insert_batch_size *= IndexingParameters::gbwt_insert_batch_size_increase_factor;
            throw RewindPlanException("[IndexRegistry]: Exceeded GBWT insert buffer size, expanding and reattempting.", {"Giraffe GBWT"});
        }
        
        output_names.push_back(output_name);
        return all_outputs;
    });
    
    // meta-recipe to either add transcripts paths or also make HST collections
    auto do_vg_rna = [merge_gbwts](const vector<const IndexFile*>& inputs,
                                   const IndexingPlan* plan,
                                   AliasGraph& alias_graph,
                                   const IndexGroup& constructing) {
    
        assert(constructing.size() == 3 || constructing.size() == 1);
        bool making_hsts = constructing.size() == 3;
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            if (making_hsts) {
                cerr << "[IndexRegistry]: Constructing haplotype-transcript GBWT and finishing spliced VG." << endl;
            }
            else {
                cerr << "[IndexRegistry]: Finishing spliced VG." << endl;
            }
            if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
                gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
            }
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        vector<vector<string>> all_outputs(constructing.size());
        IndexName output_haplo_tx, output_tx_table, output_tx_graph;
        {
            int i = 0;
            for (auto output_index : constructing) {
                if (i == 0 && making_hsts) {
                    output_haplo_tx = output_index;
                }
                else if (i == 0 && !making_hsts) {
                    output_tx_graph = output_index;
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
        //auto& haplo_tx_gbwt_names = all_outputs[0];
        auto& tx_graph_names = all_outputs[making_hsts ? 1 : 0];
        //auto& tx_table_names = all_outputs[2];
        
        vector<string> tx_filenames, gbwt_filenames, graph_filenames;
        string gbwt_filename;
        unique_ptr<gbwt::GBWT> haplotype_index;
        vector<string> gbwt_chunk_names;
        if (making_hsts) {
            tx_filenames = inputs[0]->get_filenames();
            auto gbwt_filenames = inputs[1]->get_filenames();
            graph_filenames = inputs[2]->get_filenames();
            
            assert(gbwt_filenames.size() == 1);
            gbwt_filename = gbwt_filenames.front();
            
            haplotype_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_filename);
            
            // TODO: i can't find where in the building code you actually ensure this...
            assert(haplotype_index->bidirectional());
            
            // the HST tables
            all_outputs[2].resize(graph_filenames.size());
            
            gbwt_chunk_names.resize(graph_filenames.size());
        }
        else {
            tx_filenames = inputs[0]->get_filenames();
            graph_filenames = inputs[1]->get_filenames();
        }
        
        tx_graph_names.resize(graph_filenames.size());
        
        auto haplo_tx_job = [&](int64_t i) {
            
            string tx_graph_name = plan->output_filepath(output_tx_graph, i, graph_filenames.size());
            ofstream tx_graph_outfile;
            init_out(tx_graph_outfile, tx_graph_name);
            
            string gbwt_name, info_table_name;
            if (making_hsts) {
                if (graph_filenames.size() != 1) {
                    // multiple components, so make a temp file that we will merge later
                    gbwt_name = temp_file::create();
                }
                else {
                    // one component, so we will actually save the output
                    gbwt_name = plan->output_filepath(output_haplo_tx, i, graph_filenames.size());
                }
            }
            
            int64_t j = tx_filenames.size() > 1 ? i : 0;
            
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
            size_t transcripts_added = transcriptome.add_reference_transcripts(vector<istream *>({&infile_tx}), haplotype_index, false, true);
            
            if (broadcasting_txs && !path_names.empty() && transcripts_added == 0
                && transcript_file_nonempty(tx_filenames[j])) {
                cerr << "warning:[IndexRegistry] no matching paths from transcript file " << tx_filenames[j] << " were found in graph chunk containing the following paths:" << endl;
                for (const string& path_name : path_names) {
                    cerr << "\t" << path_name << endl;
                }
            }
            
            if (making_hsts) {
                
                // go back to the beginning of the transcripts
                infile_tx.clear();
                infile_tx.seekg(0);
                
                // add edges on other haplotypes
                size_t num_transcripts_projected = transcriptome.add_haplotype_transcripts(vector<istream *>({&infile_tx}), *haplotype_index, false);
                
                // init the haplotype transcript GBWT
                size_t node_width = gbwt::bit_length(gbwt::Node::encode(transcriptome.graph().max_node_id(), true));
                bool success = execute_in_fork([&]() {
                    gbwt::GBWTBuilder gbwt_builder(node_width,
                                                   IndexingParameters::gbwt_insert_batch_size,
                                                   IndexingParameters::gbwt_sampling_interval);
                    // actually build it
                    transcriptome.add_transcripts_to_gbwt(&gbwt_builder, IndexingParameters::bidirectional_haplo_tx_gbwt, false);
                    
                    // save the haplotype transcript GBWT
                    gbwt_builder.finish();
                    save_gbwt(gbwt_builder.index, gbwt_name, IndexingParameters::verbosity == IndexingParameters::Debug);
                });
                if (!success) {
                    IndexingParameters::gbwt_insert_batch_size *= IndexingParameters::gbwt_insert_batch_size_increase_factor;
                    throw RewindPlanException("[IndexRegistry]: Exceeded GBWT insert buffer size, expanding and reattempting.",
                                              {"Haplotype-Transcript GBWT"});
                }
                
                // write transcript origin info table
                info_table_name = plan->output_filepath(output_tx_table, i, graph_filenames.size());
                ofstream info_outfile;
                init_out(info_outfile, info_table_name);
                transcriptome.write_transcript_info(&info_outfile, *haplotype_index, false);
            }
            
            // save the graph with the transcript paths added
            transcriptome.write_graph(&tx_graph_outfile);
            
            tx_graph_names[i] = tx_graph_name;
            
            if (making_hsts) {
                gbwt_chunk_names[i] = gbwt_name;
                all_outputs[2][i] = info_table_name;
            }
        };
        
        // we'll hold the gbwt in memory, so take it out of our memory budget
        int64_t target_memory_usage = plan->target_memory_usage();
        if (making_hsts) {
            target_memory_usage = max<int64_t>(0, target_memory_usage - get_file_size(gbwt_filename));
        }
        
        vector<pair<int64_t, int64_t>> approx_job_requirements;
        for (int64_t i = 0; i < graph_filenames.size(); ++i) {
            // FIXME: this should also include the approximate memory of the haplotype transcript
            approx_job_requirements.emplace_back(get_file_size(graph_filenames[i]),
                                                 approx_graph_load_memory(graph_filenames[i]));
        }
        
        JobSchedule schedule(approx_job_requirements, haplo_tx_job);
        schedule.execute(target_memory_usage);
        
        if (making_hsts) {
            // merge the GBWT chunks
            all_outputs[0].push_back(merge_gbwts(gbwt_chunk_names, plan, output_haplo_tx));
        }
        
        return all_outputs;
    };
    
    // TODO: somewhat repetitive with non-GBZ pipeline, but also some notable differences...
    auto gbz_vg_rna = [](const vector<const IndexFile*>& inputs,
                         const IndexingPlan* plan,
                         AliasGraph& alias_graph,
                         const IndexGroup& constructing) {
        
        assert(constructing.size() == 4 || constructing.size() == 2);
        bool making_hsts = constructing.size() == 4;
        assert(inputs.size() == 2 || inputs.size() == 3);
        bool projecting_transcripts = (inputs.size() == 3);
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            if (making_hsts) {
                cerr << "[IndexRegistry]: Constructing haplotype-transcript GBWT and spliced graph from GBZ-format graph." << endl;
            }
            else {
                cerr << "[IndexRegistry]: Adding splice junctions to GBZ-format graph." << endl;
            }
            if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
                gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
            }
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
                
        vector<vector<string>> all_outputs(constructing.size());
        IndexName output_haplo_tx, output_tx_table, output_tx_graph, output_max_id;
        if (making_hsts) {
            int i = 0;
            for (auto output_index : constructing) {
                if (i == 0) {
                    output_haplo_tx = output_index;
                }
                else if (i == 1) {
                    output_max_id = output_index;
                }
                else if (i == 2) {
                    output_tx_graph = output_index;
                }
                else {
                    output_tx_table = output_index;
                }
                ++i;
            }
        }
        else {
            output_max_id = *constructing.begin();
            output_tx_graph = *constructing.rbegin();
        }

        //auto& haplo_tx_gbwt_names = all_outputs[0];
        auto& max_id_names = all_outputs[making_hsts ? 1 : 0];
        auto& tx_graph_names = all_outputs[making_hsts ? 2 : 1];
        //auto& tx_table_names = all_outputs[2];
        
        auto gbz_filenames = inputs[0]->get_filenames();
        auto tx_filenames = inputs[1]->get_filenames();
        vector<string> haplo_tx_filenames;
        if (projecting_transcripts) {
            haplo_tx_filenames = inputs[2]->get_filenames();
        }
        assert(gbz_filenames.size() == 1);
        auto gbz_filename = gbz_filenames.front();
        
        vector<ifstream> infiles_tx, infiles_haplo_tx;
        for (auto& tx_filename : tx_filenames) {
            infiles_tx.emplace_back();
            init_in(infiles_tx.back(), tx_filename);
        }
        for (auto& haplo_tx_filename : haplo_tx_filenames) {
            infiles_haplo_tx.emplace_back();
            init_in(infiles_haplo_tx.back(), haplo_tx_filename);
        }
        
        string max_id_name = plan->output_filepath(output_max_id);
        ofstream max_id_outfile;
        init_out(max_id_outfile, max_id_name);
        
        string tx_graph_name = plan->output_filepath(output_tx_graph);
        ofstream tx_graph_outfile;
        init_out(tx_graph_outfile, tx_graph_name);
        
        // load, convert, and discard the GBZ
        unique_ptr<gbwt::GBWT> haplotype_index;
        auto tx_graph = init_mutable_graph();
        {
            unique_ptr<gbwtgraph::GBZ> gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(gbz_filename);
            // copy topology
            handlealgs::copy_handle_graph(&(gbz->graph), tx_graph.get());
            // copy ref paths
            gbz->graph.for_each_path_matching({PathSense::GENERIC, PathSense::REFERENCE}, {}, {}, [&](const path_handle_t& path) {
                handlegraph::algorithms::copy_path(&(gbz->graph), path, tx_graph.get());
            });
            // copy the gbwt
            haplotype_index = make_unique<gbwt::GBWT>(gbz->index);
        }
                
        // hand over the graph
        Transcriptome transcriptome(move(tx_graph));
        transcriptome.error_on_missing_path = true;
        transcriptome.feature_type = IndexingParameters::gff_feature_name;
        transcriptome.transcript_tag = IndexingParameters::gff_transcript_tag;
        
        // gather the GTF file pointers
        vector<istream *> tx_file_ptrs;
        for (auto& tx_file : infiles_tx) {
            tx_file_ptrs.push_back(&tx_file);
        }
        for (auto& tx_file : infiles_haplo_tx) {
            tx_file_ptrs.push_back(&tx_file);
        }
        
        if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
            gbwt::Verbosity::set(gbwt::Verbosity::BASIC);
        }
        else {
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        }
        
        // add the splice edges and ref transcript paths
        size_t transcripts_added = transcriptome.add_reference_transcripts(tx_file_ptrs, haplotype_index,
                                                                           !projecting_transcripts, // use haplotypes as refs for GTF?
                                                                           projecting_transcripts); // update the GBWT for the new node IDs?
        
        if (making_hsts) {
            auto& gbwt_names = all_outputs.front();
            auto& info_table_names = all_outputs.back();
            
            string gbwt_name = plan->output_filepath(output_haplo_tx);
            
            // add edges on other haplotypes
            if (projecting_transcripts) {
                
                // go back to the beginning of the transcripts
                for (auto& tx_file : infiles_tx) {
                    tx_file.clear();
                    tx_file.seekg(0);
                }
                for (auto& tx_file : infiles_haplo_tx) {
                    tx_file.clear();
                    tx_file.seekg(0);
                }
                
                transcriptome.add_haplotype_transcripts(tx_file_ptrs, *haplotype_index, false);
            }
            
            // init the haplotype transcript GBWT
            size_t node_width = gbwt::bit_length(gbwt::Node::encode(transcriptome.graph().max_node_id(), true));
            bool success = execute_in_fork([&]() {
                gbwt::GBWTBuilder gbwt_builder(node_width,
                                               IndexingParameters::gbwt_insert_batch_size,
                                               IndexingParameters::gbwt_sampling_interval);
                // actually build it
                transcriptome.add_transcripts_to_gbwt(&gbwt_builder, IndexingParameters::bidirectional_haplo_tx_gbwt, false);
                
                // save the haplotype transcript GBWT
                gbwt_builder.finish();
                save_gbwt(gbwt_builder.index, gbwt_name, IndexingParameters::verbosity == IndexingParameters::Debug);
            });
            if (!success) {
                IndexingParameters::gbwt_insert_batch_size *= IndexingParameters::gbwt_insert_batch_size_increase_factor;
                throw RewindPlanException("[IndexRegistry]: Exceeded GBWT insert buffer size, expanding and reattempting.",
                                          {"Haplotype-Transcript GBWT"});
            }
            
            // write transcript origin info table
            string info_table_name = plan->output_filepath(output_tx_table);
            ofstream info_outfile;
            init_out(info_outfile, info_table_name);
            transcriptome.write_transcript_info(&info_outfile, *haplotype_index, false);
            
            gbwt_names.push_back(gbwt_name);
            info_table_names.push_back(info_table_name);
        }
        
        
        // save the graph with the transcript paths added
        transcriptome.write_graph(&tx_graph_outfile);
        tx_graph_names.push_back(tx_graph_name);
        
        // write the max ID as well
        max_id_outfile << transcriptome.graph().max_node_id();
        max_id_names.push_back(max_id_name);
        
        return all_outputs;
    };
    
    
    auto vg_rna_gbz_graph_only =
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"}, {"GBZ", "GTF/GFF", "Haplotype GTF/GFF"},
                             [gbz_vg_rna](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          AliasGraph& alias_graph,
                                          const IndexGroup& constructing) {
        return gbz_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    auto vg_rna_gbz_liftover_graph_only =
    registry.register_recipe({"Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"}, {"GBZ", "GTF/GFF"},
                             [gbz_vg_rna](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          AliasGraph& alias_graph,
                                          const IndexGroup& constructing) {
        return gbz_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    auto vg_rna_gbz_full =
    registry.register_recipe({"Haplotype-Transcript GBWT", "Spliced MaxNodeID", "Spliced VG w/ Transcript Paths", "Unjoined Transcript Origin Table"},
                             {"GBZ", "GTF/GFF", "Haplotype GTF/GFF"},
                             [gbz_vg_rna](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          AliasGraph& alias_graph,
                                          const IndexGroup& constructing) {
        return gbz_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    auto vg_rna_gbz_liftover_full =
    registry.register_recipe({"Haplotype-Transcript GBWT", "Spliced MaxNodeID", "Spliced VG w/ Transcript Paths", "Unjoined Transcript Origin Table"},
                             {"GBZ", "GTF/GFF"},
                             [gbz_vg_rna](const vector<const IndexFile*>& inputs,
                                          const IndexingPlan* plan,
                                          AliasGraph& alias_graph,
                                          const IndexGroup& constructing) {
        return gbz_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    auto vg_rna_graph_only =
    registry.register_recipe({"Spliced VG w/ Transcript Paths"},
                             {"Chunked GTF/GFF", "Spliced VG"},
                             [do_vg_rna](const vector<const IndexFile*>& inputs,
                                         const IndexingPlan* plan,
                                         AliasGraph& alias_graph,
                                         const IndexGroup& constructing) {
   
        return do_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    auto vg_rna_full =
    registry.register_recipe({"Haplotype-Transcript GBWT", "Spliced VG w/ Transcript Paths", "Unjoined Transcript Origin Table"},
                             {"Chunked GTF/GFF", "Spliced GBWT", "Spliced VG"},
                             [do_vg_rna](const vector<const IndexFile*>& inputs,
                                         const IndexingPlan* plan,
                                         AliasGraph& alias_graph,
                                         const IndexGroup& constructing) {
        
        return do_vg_rna(inputs, plan, alias_graph, constructing);
    });
    
    // if both the full and graph-only are required, only do the full
    registry.register_generalization(vg_rna_full, vg_rna_graph_only);
    registry.register_generalization(vg_rna_gbz_full, vg_rna_gbz_graph_only);
    registry.register_generalization(vg_rna_gbz_liftover_full, vg_rna_gbz_liftover_graph_only);
    
    ////////////////////////////////////
    // Info Table Recipes
    ////////////////////////////////////
    
    registry.register_recipe({"Transcript Origin Table"}, {"Unjoined Transcript Origin Table"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
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
                getline(infile, line);
            }
            while (infile.good()) {
                getline(infile, line);
                if (!line.empty()) {
                    outfile << line << '\n';
                }
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
    auto prune_graph = [](const vector<const IndexFile*>& inputs,
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
        
        pruned_graph_names.resize(graph_names.size());
        
        mutex unfold_lock;
        
        auto prune_job = [&](int64_t i) {
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
                    unfolder.restore_paths(*graph, IndexingParameters::verbosity >= IndexingParameters::Debug);
                }
                else {
                    // we can expand out complex regions using haplotypes as well as paths
                    
                    // TODO: can't do this fully in parallel because each chunk needs to modify
                    // the same mapping
                    // TODO: it's a bit inelegant that i keep overwriting the mapping...
                    PhaseUnfolder unfolder(*unpruned_graph, *gbwt_index, max_node_id + 1);
                    unfold_lock.lock();
                    unfolder.read_mapping(mapping_name);
                    unfolder.unfold(*graph, IndexingParameters::verbosity >= IndexingParameters::Debug);
                    unfolder.write_mapping(mapping_name);
                    unfold_lock.unlock();
                }
            }
            
            vg::io::save_handle_graph(graph.get(), outfile_vg);
            
            pruned_graph_names[i] = vg_output_name;
        };
        
        int64_t target_memory_usage = plan->target_memory_usage();
        
        if (using_haplotypes) {
            // we only need to load the GBWT once, so we take it out of the shared budget
            target_memory_usage -= get_file_size(gbwt_name);
        }
        vector<pair<int64_t, int64_t>> approx_job_requirements;
        for (int64_t i = 0; i < graph_names.size(); ++i) {
            // for paths, double the memory because we'll probably need to re-load the graph to restore paths
            approx_job_requirements.emplace_back(get_file_size(graph_names[i]),
                                                 (using_haplotypes ? 1 : 2) * approx_graph_load_memory(graph_names[i]));
        }
        
        JobSchedule schedule(approx_job_requirements, prune_job);
        schedule.execute(target_memory_usage);
        
        return all_outputs;
    };
    
    registry.register_recipe({"Pruned VG"}, {"MaxNodeID", "VG"},
                             [prune_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Haplotype-Pruned VG", "Unfolded NodeMapping"}, {"GBWT", "MaxNodeID", "VG"},
                             [prune_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Pruning complex regions of VG to prepare for GCSA indexing with GBWT unfolding." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Pruned Spliced VG"}, {"Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"},
                             [prune_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Pruning complex regions of spliced VG to prepare for GCSA indexing." << endl;
        }
        // call the meta-recipe
        return prune_graph(inputs, plan, constructing);
    });
    
    // TODO: would it be better to use the Haplotype-Transcript GBWT, or maybe to join them?
    // the splice edges will be covered by the transcript paths, so it won't be too bad
    registry.register_recipe({"Haplotype-Pruned Spliced VG", "Unfolded Spliced NodeMapping"},
                             {"Spliced GBWT", "Spliced MaxNodeID", "Spliced VG w/ Transcript Paths"},
                             [prune_graph](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
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
    auto construct_gcsa = [](const vector<const IndexFile*>& inputs,
                             const IndexingPlan* plan,
                             const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
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
        
        if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
            gcsa::Verbosity::set(gcsa::Verbosity::BASIC);
        }
        else {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }
        auto params = gcsa::ConstructionParameters();
        params.setSteps(IndexingParameters::gcsa_doubling_steps);
        params.setLimitBytes(IndexingParameters::gcsa_size_limit);
                
#ifdef debug_index_registry_recipes
        cerr << "enumerating k-mers for input pruned graphs:" << endl;
        for (auto& name : graph_filenames) {
            cerr << "\t" << name << endl;
        }
#endif
        // if indexing fails, we'll rewind to whichever of these we used
        IndexGroup pruned_graphs{"Pruned VG", "Pruned Spliced VG", "Haplotype-Pruned VG", "Haplotype-Pruned Spliced VG"};

        VGset graph_set(graph_filenames);
        size_t kmer_bytes = params.getLimitBytes();
        vector<string> dbg_names;
        try {
            dbg_names = graph_set.write_gcsa_kmers_binary(IndexingParameters::gcsa_initial_kmer_length, kmer_bytes);
        }
        catch (SizeLimitExceededException& ex) {
            // update pruning params
            IndexingParameters::pruning_walk_length *= IndexingParameters::pruning_walk_length_increase_factor;
            IndexingParameters::pruning_max_node_degree *= IndexingParameters::pruning_max_node_degree_decrease_factor;
            string msg = "[IndexRegistry]: Exceeded disk use limit while generating k-mers. "
                         "Rewinding to pruning step with more aggressive pruning to simplify the graph.";
            throw RewindPlanException(msg, pruned_graphs);
        }
        
        bool success = execute_in_fork([&]() {
#ifdef debug_index_registry_recipes
            cerr << "making GCSA2 at " << gcsa_output_name << " and " << lcp_output_name << " after writing de Bruijn graph files to:" << endl;
            for (auto dbg_name : dbg_names) {
                cerr << "\t" << dbg_name << endl;
            }
#endif
            
            // construct the indexes (giving empty mapping name is sufficient to make
            // indexing skip the unfolded code path)
            gcsa::InputGraph input_graph(dbg_names, true, gcsa::Alphabet(),
                                         mapping_filename);
            gcsa::GCSA gcsa_index(input_graph, params);
            gcsa::LCPArray lcp_array(input_graph, params);
            
#ifdef debug_index_registry_recipes
            cerr << "saving GCSA/LCP pair" << endl;
#endif
            
            save_gcsa(gcsa_index, gcsa_output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
            save_lcp(lcp_array, lcp_output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
        });
        
        // clean up the k-mer files
        for (auto dbg_name : dbg_names) {
            temp_file::remove(dbg_name);
        }
        
        if (!success) {
            // the indexing was not successful, presumably because of exponential disk explosion
            
            // update pruning params
            IndexingParameters::pruning_walk_length *= IndexingParameters::pruning_walk_length_increase_factor;
            IndexingParameters::pruning_max_node_degree *= IndexingParameters::pruning_max_node_degree_decrease_factor;
            string msg = "[IndexRegistry]: Exceeded disk use limit while performing k-mer doubling steps. "
                         "Rewinding to pruning step with more aggressive pruning to simplify the graph.";
            throw RewindPlanException(msg, pruned_graphs);
        }
        
        gcsa_names.push_back(gcsa_output_name);
        lcp_names.push_back(lcp_output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"GCSA", "LCP"}, {"Haplotype-Pruned VG", "Unfolded NodeMapping"},
                             [construct_gcsa](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"GCSA", "LCP"}, {"Pruned VG"},
                             [construct_gcsa](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GCSA", "Spliced LCP"}, {"Haplotype-Pruned Spliced VG", "Unfolded Spliced NodeMapping"},
                             [construct_gcsa](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        // execute meta recipe
        return construct_gcsa(inputs, plan, constructing);
    });
    
    registry.register_recipe({"Spliced GCSA", "Spliced LCP"}, {"Pruned Spliced VG"},
                             [construct_gcsa](const vector<const IndexFile*>& inputs,
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
    auto find_snarls = [](const HandleGraph& graph,
                          const IndexingPlan* plan,
                          const IndexGroup& constructing) {
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto output_snarls = *constructing.begin();
        auto& snarl_names = all_outputs[0];
        
        string output_name = plan->output_filepath(output_snarls);
        ofstream outfile;
        init_out(outfile, output_name);
                
        // find snarls
        unique_ptr<SnarlFinder> snarl_finder = unique_ptr<SnarlFinder>(new IntegratedSnarlFinder(graph));
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

    
    // TODO: disabling so that we can distinguish giraffe graphs that may have
    // different node IDs from the GFA import
//    registry.register_recipe({"Snarls"}, {"XG"},
//                             [find_snarls](const vector<const IndexFile*>& inputs,
//                                 const IndexingPlan* plan,
//                                 AliasGraph& alias_graph,
//                                 const IndexGroup& constructing) {
//        if (IndexingParameters::verbosity != IndexingParameters::None) {
//            cerr << "[IndexRegistry]: Finding snarls in graph." << endl;
//        }
//        return find_snarls(inputs, plan, constructing);
//    });
    
    registry.register_recipe({"Spliced Snarls"}, {"Spliced XG"},
                             [find_snarls](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Finding snarls in spliced graph." << endl;
        }
        
        assert(inputs.size() == 1);
        auto graph_filenames = inputs[0]->get_filenames();
        assert(graph_filenames.size() == 1);
        auto graph_filename = graph_filenames.front();
        
        ifstream infile;
        init_in(infile, graph_filename);
        unique_ptr<HandleGraph> graph = vg::io::VPKG::load_one<HandleGraph>(infile);
        
        return find_snarls(*graph, plan, constructing);
    });
    
    ////////////////////////////////////
    // Distance Index Recipes
    ////////////////////////////////////
    
    
    // meta-recipe to make distance index
    auto make_distance_index = [](const HandleGraph& graph,
                                  const IndexingPlan* plan,
                                  const IndexGroup& constructing) {
        
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto dist_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        string output_name = plan->output_filepath(dist_output);
        ofstream outfile;
        init_out(outfile, output_name);
        
        SnarlDistanceIndex distance_index;
        IntegratedSnarlFinder snarl_finder(graph);
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        distance_index.serialize(output_name);
        
        output_names.push_back(output_name);
        return all_outputs;
    };
    
    registry.register_recipe({"Giraffe Distance Index"}, {"Giraffe GBZ"},
                             [make_distance_index](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing distance index for Giraffe." << endl;
        }
        
        assert(inputs.size() == 1);
        auto& gbz_filenames = inputs[0]->get_filenames();
        assert(gbz_filenames.size() == 1);
        auto gbz_filename = gbz_filenames.front();
        
        ifstream infile_gbz;
        init_in(infile_gbz, gbz_filename);
        unique_ptr<gbwtgraph::GBZ> gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(infile_gbz);
        
        return make_distance_index(gbz->graph, plan, constructing);
    });
    
    registry.register_recipe({"Spliced Distance Index"}, {"Spliced XG"},
                             [make_distance_index](const vector<const IndexFile*>& inputs,
                                 const IndexingPlan* plan,
                                 AliasGraph& alias_graph,
                                 const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing distance index for a spliced graph." << endl;
        }
        
        assert(inputs.size() == 1);
        auto& graph_filenames = inputs[0]->get_filenames();
        assert(graph_filenames.size() == 1);
        auto graph_filename = graph_filenames.front();
        
        ifstream infile_graph;
        init_in(infile_graph, graph_filename);
        
        unique_ptr<HandleGraph> graph = vg::io::VPKG::load_one<HandleGraph>(infile_graph);
        
        return make_distance_index(*graph, plan, constructing);
    });
    
    ////////////////////////////////////
    // GBZ Recipes
    ////////////////////////////////////
    
    registry.register_recipe({"Giraffe GBZ"}, {"GBZ"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        alias_graph.register_alias(*constructing.begin(), inputs[0]);
        return vector<vector<string>>(1, inputs.front()->get_filenames());
    });
    
    registry.register_recipe({"GBZ"}, {"Reference GFA w/ Haplotypes"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing a GBZ from GFA input." << endl;
        }
        
        assert(inputs.size() == 1);
        if (inputs[0]->get_filenames().size() != 1) {
            cerr << "error:[IndexRegistry] Graph construction does not support multiple GFAs at this time." << endl;
            exit(1);
        }
        auto gfa_filename = inputs[0]->get_filenames().front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto gbz_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        string output_name = plan->output_filepath(gbz_output);
        
        gbwtgraph::GFAParsingParameters params = get_best_gbwtgraph_gfa_parsing_parameters();
        // TODO: there's supposedly a heuristic to set batch size that could perform better than this global param,
        // but it would be kind of a pain to update it like we do the global param
        params.batch_size = IndexingParameters::gbwt_insert_batch_size;
        params.sample_interval = IndexingParameters::gbwt_sampling_interval;
        params.max_node_length = IndexingParameters::max_node_size;
        params.show_progress = IndexingParameters::verbosity == IndexingParameters::Debug;
        
        bool success = execute_in_fork([&]() {
            
            // jointly generate the GBWT and record sequences
            unique_ptr<gbwt::GBWT> gbwt_index;
            unique_ptr<gbwtgraph::SequenceSource> seq_source;
            tie(gbwt_index, seq_source) = gbwtgraph::gfa_to_gbwt(gfa_filename, params);
            
            // convert sequences into gbwt graph
            gbwtgraph::GBWTGraph gbwt_graph(*gbwt_index, *seq_source);
            
            // save together as a GBZ
            save_gbz(*gbwt_index, gbwt_graph, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
        });
        if (!success) {
            IndexingParameters::gbwt_insert_batch_size *= IndexingParameters::gbwt_insert_batch_size_increase_factor;
            throw RewindPlanException("[IndexRegistry]: Exceeded GBWT insert buffer size, expanding and reattempting.",
                                      {"Giraffe GBZ"});
        }
        
        output_names.push_back(output_name);
        return all_outputs;
    });

    registry.register_recipe({"Giraffe GBZ"}, {"GBWTGraph", "Giraffe GBWT"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Combining Giraffe GBWT and GBWTGraph into GBZ." << endl;
        }

        assert(inputs.size() == 2);
        auto gbwt_filenames = inputs[1]->get_filenames();
        auto gg_filenames = inputs[0]->get_filenames();
        assert(gbwt_filenames.size() == 1);
        assert(gg_filenames.size() == 1);
        auto gbwt_filename = gbwt_filenames.front();
        auto gg_filename = gg_filenames.front();

        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto gbz_output = *constructing.begin();
        auto& output_names = all_outputs[0];

        gbwtgraph::GBZ gbz;
        load_gbz(gbz, gbwt_filename, gg_filename, IndexingParameters::verbosity == IndexingParameters::Debug);

        string output_name = plan->output_filepath(gbz_output);
        save_gbz(gbz, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);

        output_names.push_back(output_name);
        return all_outputs;
    });
    
    // Thses used to be a GBWTGraph recipe, but we don't want to produce GBWTGraphs anymore.
    
    registry.register_recipe({"Giraffe GBZ"}, {"Giraffe GBWT", "NamedNodeBackTranslation", "XG"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing GBZ using NamedNodeBackTranslation." << endl;
        }
        
        assert(inputs.size() == 3);
        auto gbwt_filenames = inputs[0]->get_filenames();
        auto translation_filenames = inputs[1]->get_filenames();
        auto xg_filenames = inputs[2]->get_filenames();
        assert(gbwt_filenames.size() == 1);
        assert(translation_filenames.size() == 1);
        assert(xg_filenames.size() == 1);
        auto gbwt_filename = gbwt_filenames.front();
        auto translation_filename = translation_filenames.front(); 
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto gbz_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        ifstream infile_xg;
        init_in(infile_xg, xg_filename);
        auto xg_index = vg::io::VPKG::load_one<xg::XG>(infile_xg);
        
        ifstream infile_translation;
        init_in(infile_translation, translation_filename);
        // There's only one implementation we can use here at the moment, so
        // don't bother with the normal loader/saver system.
        FlatFileBackTranslation translation(infile_translation);

        gbwtgraph::GBZ gbz;
        load_gbwt(gbz.index, gbwt_filename, IndexingParameters::verbosity == IndexingParameters::Debug);
        // TODO: could add simplification to replace XG index with a gbwt::SequenceSource here
        gbz.graph = gbwtgraph::GBWTGraph(gbz.index, *xg_index, &translation);

        string output_name = plan->output_filepath(gbz_output);
        save_gbz(gbz, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);

        output_names.push_back(output_name);
        return all_outputs;
    });

    registry.register_recipe({"Giraffe GBZ"}, {"Giraffe GBWT", "XG"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing GBZ." << endl;
        }
        
        assert(inputs.size() == 2);
        auto gbwt_filenames = inputs[0]->get_filenames();
        auto xg_filenames = inputs[1]->get_filenames();
        assert(gbwt_filenames.size() == 1);
        assert(xg_filenames.size() == 1);
        auto gbwt_filename = gbwt_filenames.front();
        auto xg_filename = xg_filenames.front();
        
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto gbz_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        
        ifstream infile_xg;
        init_in(infile_xg, xg_filename);
        auto xg_index = vg::io::VPKG::load_one<xg::XG>(infile_xg);

        gbwtgraph::GBZ gbz;
        load_gbwt(gbz.index, gbwt_filename, IndexingParameters::verbosity == IndexingParameters::Debug);
        // TODO: could add simplification to replace XG index with a gbwt::SequenceSource here
        gbz.graph = gbwtgraph::GBWTGraph(gbz.index, *xg_index, algorithms::find_translation(xg_index.get()));

        string output_name = plan->output_filepath(gbz_output);
        save_gbz(gbz, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);

        output_names.push_back(output_name);
        return all_outputs;
    });

    ////////////////////////////////////
    // Minimizers Recipes
    ////////////////////////////////////

    // FIXME We may not always want to store the minimizer index. Rebuilding the index may be
    // faster than loading it from a network drive.
    registry.register_recipe({"Minimizers"}, {"Giraffe Distance Index", "Giraffe GBZ"},
                             [](const vector<const IndexFile*>& inputs,
                                const IndexingPlan* plan,
                                AliasGraph& alias_graph,
                                const IndexGroup& constructing) {
        if (IndexingParameters::verbosity != IndexingParameters::None) {
            cerr << "[IndexRegistry]: Constructing minimizer index." << endl;
        }
        
        // TODO: should the distance index input be a joint simplification to avoid serializing it?
        
        assert(inputs.size() == 2);
        auto dist_filenames = inputs[0]->get_filenames();
        auto gbz_filenames = inputs[1]->get_filenames();
        assert(dist_filenames.size() == 1);
        assert(gbz_filenames.size() == 1);
        auto dist_filename = dist_filenames.front();
        auto gbz_filename = gbz_filenames.front();
                
        assert(constructing.size() == 1);
        vector<vector<string>> all_outputs(constructing.size());
        auto minimizer_output = *constructing.begin();
        auto& output_names = all_outputs[0];
        

        ifstream infile_gbz;
        init_in(infile_gbz, gbz_filename);
        auto gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(infile_gbz);
        
        ifstream infile_dist;
        init_in(infile_dist, dist_filename);
        auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(dist_filename);
        gbwtgraph::DefaultMinimizerIndex minimizers(IndexingParameters::minimizer_k,
                                                    IndexingParameters::use_bounded_syncmers ?
                                                        IndexingParameters::minimizer_s :
                                                        IndexingParameters::minimizer_w,
                                                    IndexingParameters::use_bounded_syncmers);
                
        gbwtgraph::index_haplotypes(gbz->graph, minimizers, [&](const pos_t& pos) -> gbwtgraph::payload_type {
            return MIPayload::encode(get_minimizer_distances(*distance_index, pos));
        });
        
        string output_name = plan->output_filepath(minimizer_output);
        save_minimizer(minimizers, output_name, IndexingParameters::verbosity == IndexingParameters::Debug);
        
        output_names.push_back(output_name);
        return all_outputs;
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
        "Spliced LCP"
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_rpvg_indexes() {
    vector<IndexName> indexes{
        "Spliced XG",
        "Haplotype-Transcript GBWT",
        "Transcript Origin Table"
    };
    return indexes;
}

vector<IndexName> VGIndexes::get_default_giraffe_indexes() {
    vector<IndexName> indexes{
        "Giraffe Distance Index",
        "Giraffe GBZ",
        "Minimizers"
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

int64_t IndexingPlan::target_memory_usage() const {
    return IndexingParameters::max_memory_proportion * registry->get_target_memory_usage();
}
    
string IndexingPlan::output_filepath(const IndexName& identifier) const {
    return output_filepath(identifier, 0, 1);
}

string IndexingPlan::output_filepath(const IndexName& identifier, size_t chunk, size_t num_chunks) const {
    
    string filepath;
    if (registry->keep_intermediates ||
        (!is_intermediate(identifier) && !registry->get_index(identifier)->was_provided_directly())) {
        // we're saving this file, put it at the output prefix
        filepath = registry->output_prefix;
    }
    else {
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

set<RecipeName> IndexingPlan::dependents(const IndexName& identifier) const {
    
    set<RecipeName> dependent_steps;
    
    // seed the successors with the query
    IndexGroup successor_indexes{identifier};
    
    for (const auto& step : steps) {
                
        // TODO: should this behavior change if some of the inputs were provided directly?
        
        // collect inputs and outputs
        const auto& outputs = step.first;
        IndexGroup involved = registry->get_recipe(step).input_group();
        involved.insert(outputs.begin(), outputs.end());
        
        for (const auto& index : involved) {
            if (successor_indexes.count(index)) {
                // this is a step when a successor was either created or used
                dependent_steps.insert(step);
                // outputs are also successors
                successor_indexes.insert(outputs.begin(), outputs.end());
                break;
            }
        }
    }
    return dependent_steps;
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

string IndexRegistry::get_prefix() const {
    return this->output_prefix;
}

void IndexRegistry::set_intermediate_file_keeping(bool keep_intermediates) {
    this->keep_intermediates = keep_intermediates;
}

void IndexRegistry::make_indexes(const vector<IndexName>& identifiers) {
    
    // figure out the best plan to make the objectives from the inputs
    IndexGroup identifier_group(identifiers.begin(), identifiers.end());
    auto plan = make_plan(identifier_group);
    
    // to keep track of which indexes are aliases of others
    AliasGraph alias_graph;
    
    list<RecipeName> steps_remaining(plan.get_steps().begin(), plan.get_steps().end());
    list<RecipeName> steps_completed;
    
    // execute the plan
    while (!steps_remaining.empty()) {
        // get the next step
        auto step = move(steps_remaining.front());
        steps_remaining.pop_front();
        steps_completed.push_back(step);
        
        // do the recipe
        try {
            auto recipe_results = execute_recipe(step, &plan, alias_graph);
            
            // the recipe executed successfully
            assert(recipe_results.size() == step.first.size());
            
            // record the results
            auto it = step.first.begin();
            for (const auto& results : recipe_results) {
                auto index = get_index(*it);
                // don't overwrite directly-provided inputs
                if (!index->was_provided_directly()) {
                    // and assign the new (or first) ones
                    index->assign_constructed(results);
                }
                ++it;
            }
        }
        catch (RewindPlanException& ex) {
            
            // the recipe failed, but we can rewind and retry following the recipe with
            // modified parameters (which should have been set by the exception-throwing code)
            if (IndexingParameters::verbosity != IndexingParameters::None) {
                cerr << ex.what() << endl;
            }
            // gather the recipes we're going to need to re-attempt
            const auto& rewinding_indexes = ex.get_indexes();
            set<RecipeName> dependent_recipes;
            for (const auto& index_name : rewinding_indexes) {
                assert(index_registry.count(index_name));
                for (const auto& recipe : plan.dependents(index_name)) {
                    dependent_recipes.insert(recipe);
                }
            }
            
            // move rewound steps back onto the queue
            vector<list<RecipeName>::iterator> to_move;
            for (auto it = steps_completed.rbegin(); it != steps_completed.rend(); ++it) {
                if (dependent_recipes.count(*it)) {
                    to_move.push_back(--it.base());
                }
            }
            for (auto& it : to_move) {
                steps_remaining.emplace_front(*it);
                steps_completed.erase(it);
            }
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
        
        // if the index is itself non-intermediate, it will be in the list of aliases.
        // otherwise, we can alias one index by moving instead of copying
        auto f = find(aliasors.begin(), aliasors.end(), aliasee);
        bool can_move = f == aliasors.end() && !get_index(aliasee)->was_provided_directly();
        if (!can_move) {
            // just remove the "alias" so we don't need to deal with it
            std::swap(*f, aliasors.back());
            aliasors.pop_back();
        }
        
        const auto& aliasee_filenames = get_index(aliasee)->get_filenames();
        
        // copy aliases for any that we need to (start past index 0 if we can move it)
        for (size_t i = can_move; i < aliasors.size(); ++i) {
            for (size_t j = 0; j < aliasee_filenames.size(); ++j) {
                
                auto copy_filename = plan.output_filepath(aliasors[i], j, aliasee_filenames.size());
                copy_file(aliasee_filenames[j], copy_filename);
            }
        }
        // if we can move the aliasee (i.e. it is intermediate), then make
        // one index by moving instead of copying
        if (can_move) {
            for (size_t j = 0; j < aliasee_filenames.size(); ++j) {
                auto move_filename = plan.output_filepath(aliasors[0], j, aliasee_filenames.size());
                int code = rename(aliasee_filenames[j].c_str(), move_filename.c_str());
                if (code) {
                    // moving failed (maybe because the files on separate drives?) fall back on copying
                    copy_file(aliasee_filenames[j], move_filename);
                }
            }
        }
    }
    
    // Keep all the indexes around. If you want to re-use the object for a
    // different set of indexes, you will need to call reset() yourself.
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
    provide(identifier, vector<string>(1, filename));
}

void IndexRegistry::provide(const IndexName& identifier, const vector<string>& filenames) {
    if (IndexingParameters::verbosity >= IndexingParameters::Debug) {
        cerr << "[IndexRegistry]: Provided: " << identifier << endl;
    }
    if (!index_registry.count(identifier)) {
        cerr << "error:[IndexRegistry] cannot provide unregistered index: " << identifier << endl;
        exit(1);
    }
    get_index(identifier)->provide(filenames);
}

bool IndexRegistry::available(const IndexName& identifier) const {
    if (!index_registry.count(identifier)) {
        // Index is not registered
        return false;
    }
    const IndexFile* index = get_index(identifier);
    if (!index->is_finished()) {
        // Index is not made
        return false;
    }
    return true;
}

vector<string> IndexRegistry::require(const IndexName& identifier) const {
    if (!index_registry.count(identifier)) {
        cerr << "error:[IndexRegistry] cannot require unregistered index: " << identifier << endl;
        exit(1);
    }
    const IndexFile* index = get_index(identifier);
    if (!index->is_finished()) {
        cerr << "error:[IndexRegistry] do not have and did not make index: " << identifier << endl;
        exit(1);
    }
    return index->get_filenames();
}

void IndexRegistry::set_target_memory_usage(int64_t bytes) {
    target_memory_usage = bytes;
}

int64_t IndexRegistry::get_target_memory_usage() const {
    return target_memory_usage;
}

// from https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
int64_t IndexRegistry::get_system_memory() {
    int64_t pages = sysconf(_SC_PHYS_PAGES);
    int64_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
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
    // this is an easy-to-troubleshoot check that lets us use IndexGroup's (which are ordered set)
    // internally and still provide the vector<IndexFile*> in the same order as the input identifiers to
    // the RecipeFunc and in the recipe declaration.
    // i.e. this helps ensure that the order of the indexes that you code in the recipe declaration
    // is the order that they will continue to be given throughout the registry
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
    
    bool first_group_entry = !recipe_registry.count(output_group);
    recipe_registry[output_group].emplace_back(inputs, exec);
    RecipeName name(output_group, recipe_registry[output_group].size() - 1);
        
    if (output_group.size() > 1 && first_group_entry) {
        // add unboxing recipes at the same priority level as the full recipe
        auto it = output_group.begin();
        for (size_t i = 0; i < identifiers.size(); ++i) {
#ifdef debug_index_registry_setup
            cerr << "registering unboxing recipe from " << to_string(output_group) << " to " << *it << endl;
#endif
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

void IndexRegistry::register_generalization(const RecipeName& generalizer, const RecipeName& generalizee) {
    for (const auto& index_name : generalizee.first) {
        if (!generalizer.first.count(index_name)) {
            cerr << "error:[IndexRegistry] registered a generalization that does not contain generalizee's output " << index_name << endl;
            exit(1);
        }
    }
    const auto& generalizer_recipe = recipe_registry.at(generalizer.first).at(generalizer.second);
    const auto& generalizee_recipe = recipe_registry.at(generalizee.first).at(generalizee.second);
    for (const auto& index_name : generalizee_recipe.input_group()) {
        if (!generalizer_recipe.input_group().count(index_name)) {
            cerr << "error:[IndexRegistry] registered a generalization that does not contain generalizee's input " << index_name << endl;
            exit(1);
        }
    }
    
    generalizations[generalizee] = generalizer;
}

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

bool IndexRegistry::vcf_is_phased(const string& filepath) {
    
    if (IndexingParameters::verbosity >= IndexingParameters::Basic) {
        cerr << "[IndexRegistry]: Checking for phasing in VCF(s)." << endl;
    }
    
    
    htsFile* file = hts_open(filepath.c_str(), "rb");
    if (!file) {
        cerr << "error:[IndexRegistry]: Failed to open VCF file: " << filepath << endl;
        exit(1);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(file);
    int phase_set_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "PS");
    // note: it seems that this is not necessary for expressing phasing after all
    int num_samples = bcf_hdr_nsamples(hdr);
    
    // iterate over records
    bcf1_t* line = bcf_init();
    bool found_phased = false;
    // check about 30k non-haploid variants before concluding that the VCF isn't phased
    // TODO: will there be contig ordering biases that make this a bad assumption?
    constexpr int nonhap_vars_to_check = 1 << 15;
    int nonhap_iter = 0;
    while (bcf_read(file, hdr, line) >= 0 && nonhap_iter < nonhap_vars_to_check && !found_phased) {
        if (phase_set_id >= 0) {
            if (phase_set_id == BCF_HT_INT) {
                // phase sets are integers
                int num_phase_set_arr = 0;
                int32_t* phase_sets = NULL;
                int num_phase_sets = bcf_get_format_int32(hdr, line, "PS", &phase_sets, &num_phase_set_arr);
                for (int i = 0; i < num_phase_sets && !found_phased; ++i) {
                    found_phased = phase_sets[i] != 0;
                }
                free(phase_sets);
            }
            else if (phase_set_id == BCF_HT_STR) {
                // phase sets are strings
                int num_phase_set_arr = 0;
                char** phase_sets = NULL;
                int num_phase_sets = bcf_get_format_string(hdr, line, "PS", &phase_sets, &num_phase_set_arr);
                for (int i = 0; i < num_phase_sets && !found_phased; ++i) {
                    found_phased = strcmp(phase_sets[i], ".") != 0;
                }
                if (phase_sets) {
                    // all phase sets are concatenated in one malloc's char*, pointed to by the first pointer
                    free(phase_sets[0]);
                }
                // free the array of pointers
                free(phase_sets);
            }
        }
        
        // init a genotype array
        int32_t* genotypes = nullptr;
        int arr_size = 0;
        // and query it
        int num_genotypes = bcf_get_genotypes(hdr, line, &genotypes, &arr_size);
        if (num_genotypes >= 0) {
            // we got genotypes, check to see if they're phased.
            // We know we can't have genotypes if there are 0 samples.
            int ploidy = num_genotypes / num_samples;
            if (ploidy > 1) {
                for (int i = 0; i < num_genotypes && !found_phased; i += ploidy) {
                    for (int j = 0; j < ploidy && !found_phased; ++j) {
                        if (genotypes[i + j] == bcf_int32_vector_end) {
                            // sample has lower ploidy
                            break;
                        }
                        if (bcf_gt_is_missing(genotypes[i + j])) {
                            continue;
                        }
                        if (bcf_gt_is_phased(genotypes[i + j])) {
                            // the VCF expresses phasing, we can
                            found_phased = true;;
                        }
                    }
                }
                ++nonhap_iter;
            }
        }
        
        free(genotypes);
    }
    if (nonhap_iter == 0 && num_samples > 0) {
        // We looked at some samples and none of them had any non-haploid genotypes.
        // Assume the entire VCF is haploid, which are trivially phased
        found_phased = true;
    }
    // clean up
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(file);
    return found_phased;
}

bool IndexRegistry::gfa_has_haplotypes(const string& filepath) {
    if (IndexingParameters::verbosity >= IndexingParameters::Basic) {
        cerr << "[IndexRegistry]: Checking for haplotype lines in GFA." << endl;
    }
    ifstream strm(filepath);
    if (!strm) {
        cerr << "error:[IndexRegistry] Could not open GFA file " << filepath << endl;
        exit(1);
    }
    while (strm.good()) {
        char line_type = strm.get();
        if (line_type == 'W') {
            return true;
        }
        strm.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    return false;
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
        
        // records of (identifier, requesters, ordinal index of recipe selected)
        vector<tuple<size_t, set<size_t>, size_t>> plan_path;
        
        // map dependency priority to requesters
        map<size_t, set<size_t>, greater<size_t>> queue;
        
        auto num_recipes = [&](const IndexGroup& indexes) {
            int64_t num = 0;
            if (recipe_registry.count(indexes)) {
                num = recipe_registry.at(indexes).size();
            }
            return num;
        };
        
        // update the queue to request the inputs of a recipe from the final index on the plan path
        auto request_from_back = [&]() {
            
            // the index at the back of the plan path is making the request
            auto& requester = identifier_order[get<0>(plan_path.back())];
            // get its next recipe
            auto inputs = recipe_registry.at(requester).at(get<2>(plan_path.back())).input_group();
            
#ifdef debug_index_registry
            cerr << "index(es) " << to_string(requester) << " can be made by a recipe requiring " << to_string(inputs) << endl;
#endif
            
            if (requester.size() == 1 && inputs.count(*requester.begin())) {
                // this is an unboxing recipe, request the whole previous group
                queue[dep_order_of_identifier[inputs]].insert(get<0>(plan_path.back()));
            }
            else {
                // this is not an unboxing recipe, request all of the recipe inputs separately
                for (auto& input_index : inputs) {
                    IndexGroup singleton_input{input_index};
                    queue[dep_order_of_identifier[singleton_input]].insert(get<0>(plan_path.back()));
                }
            }
        };
        
        // place the final step in the plan path back in the queue
        auto requeue_back = [&]() {
#ifdef debug_index_registry
            cerr << "requeueing " << to_string(identifier_order[get<0>(plan_path.back())]) << ", requested by:" << endl;
            for (auto d : get<1>(plan_path.back())) {
                cerr << "\t" << to_string(identifier_order[d]) << endl;
            }
#endif
            // TODO: is this check necessary?
            if (!get<1>(plan_path.back()).empty()) {
                queue[get<0>(plan_path.back())] = get<1>(plan_path.back());
            }
            plan_path.pop_back();
        };
        
        // update the queue to remove requests to the inputs of a recipe from the final index on the plan path
        auto unrequest_from_back = [&]() {
            
            auto make_unrequest = [&](const IndexGroup& inputs) {
                auto it = queue.find(dep_order_of_identifier[inputs]);
                it->second.erase(get<0>(plan_path.back()));
                if (it->second.empty()) {
#ifdef debug_index_registry
                    cerr << "\t\tremoved final request to " << to_string(identifier_order[it->first]) << ", dequeuing" << endl;
#endif
                    queue.erase(it);
                }
            };
            
            auto& requester = identifier_order[get<0>(plan_path.back())];
            
            if (!all_finished(requester) && recipe_registry.count(requester)) {
                // this index was using a recipe, we need to update its dependencies
                // that are currently in the queue
                
#ifdef debug_index_registry
                cerr << "retracting requests from " << to_string(requester) << ", recipe " << get<2>(plan_path.back()) << endl;
#endif
                auto inputs = recipe_registry.at(requester).at(get<2>(plan_path.back())).input_group();
                
#ifdef debug_index_registry
                cerr << "\tmade requests from recipe requiring " << to_string(inputs) << endl;
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
#ifdef debug_index_registry
            else {
                cerr << "no need to retract requests from " << to_string(requester) << endl;
            }
#endif
        };
        
        // init the queue
        queue[dep_order_of_identifier[product_group]] = set<size_t>();
        
        while (!queue.empty()) {
#ifdef debug_index_registry_path_state
            cerr << "new iteration, path:" << endl;
            for (auto pe : plan_path) {
                cerr << "\t" << to_string(identifier_order[get<0>(pe)]) << ", requesters:";
                if (get<1>(pe).empty())  {
                    cerr << " PLAN TARGET";
                }
                else {
                    for (auto d : get<1>(pe)) {
                        cerr << " " << to_string(identifier_order[d]);
                    }
                }
                cerr << ", recipe " << get<2>(pe) << endl;
            }
            cerr << "state of queue:" << endl;
            for (auto q : queue) {
                cerr << "\t" << to_string(identifier_order[q.first]) << ", requesters:";
                if (q.second.empty())  {
                    cerr << " PLAN TARGET";
                }
                else {
                    for (auto d : q.second) {
                        cerr << " " << to_string(identifier_order[d]);
                    }
                }
                cerr << endl;
            }
#endif
            
            // get the latest file in the dependency order that we have left to build
            auto it = queue.begin();
            plan_path.emplace_back(it->first, it->second, 0);
            
#ifdef debug_index_registry
            cerr << "dequeue " << to_string(identifier_order[it->first]) << " requested from:" << endl;
            if (it->second.empty()) {
                cerr << "\tPLAN TARGET" << endl;
            }
            else {
                for (auto requester : it->second) {
                    cerr << "\t" << to_string(identifier_order[requester]) << endl;
                }
            }
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
                cerr << "index " << to_string(index_group) << " cannot be made from existing inputs, need to backtrack" << endl;
#endif
                
                // prune to requester and advance to its next recipe, as many times as necessary until
                // requester has remaining un-tried recipes
                // note: if we're backtracking from a data file it might not have recipes
                while (get<2>(plan_path.back()) >= num_recipes(identifier_order[get<0>(plan_path.back())])) {
                    // there are no remaining recipes to build the last index in the plan
                    
                    if (get<1>(plan_path.back()).empty()) {
                        // this is the product of the plan path, and we're out of recipes for it
                        throw InsufficientInputException(product, *this);
                    }
                    
                    // remove items off the plan path until we get to the first index that requested
                    // this one
                    size_t requester = *get<1>(plan_path.back()).rbegin();
                    
#ifdef debug_index_registry
                    cerr << "no remaining recipes for " << to_string(identifier_order[get<0>(plan_path.back())]) << ", pruning to earliest requester: " << to_string(identifier_order[requester]) << endl;
#endif
                    
                    requeue_back(); // nothing to unrequest from the first one, which is past its last recipe
                    while (get<0>(plan_path.back()) != requester) {
                        unrequest_from_back();
                        requeue_back();
                    }
                    
                    // advance to the next recipe
                    unrequest_from_back();
                    ++get<2>(plan_path.back());
                    
#ifdef debug_index_registry
                    cerr << "advance to recipe " << get<2>(plan_path.back()) << " for " << to_string(identifier_order[get<0>(plan_path.back())]) << endl;
#endif
                }
                
                // we pulled back far enough that we found an index with a lower-priority recipe
                // remaining
                request_from_back();
            }
            
        }
        
#ifdef debug_index_registry
        cerr << "final plan path for index " << product << ":" << endl;
        for (auto path_elem : plan_path) {
            cerr << "\t" << to_string(identifier_order[get<0>(path_elem)]) << ", recipe " << get<2>(path_elem) << ", from:" << endl;
            for (auto d : get<1>(path_elem)) {
                cerr << "\t\t" << to_string(identifier_order[d]) << endl;
            }
        }
#endif
        
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
    cerr << "plan before applying generalizations:" << endl;
    for (auto plan_elem : plan.steps) {
        cerr << "\t" << to_string(plan_elem.first) << " " << plan_elem.second << endl;
    }
#endif
    
    // remove generalizees if we used their generalizers
    set<RecipeName> plan_set(plan.steps.begin(), plan.steps.end());
    plan.steps.resize(remove_if(plan.steps.begin(), plan.steps.end(), [&](const RecipeName& recipe) {
        return generalizations.count(recipe) && plan_set.count(generalizations.at(recipe));
    }) - plan.steps.begin());
    
#ifdef debug_index_registry
    cerr << "full plan including provided files:" << endl;
    for (auto plan_elem : plan.steps) {
        cerr << "\t" << to_string(plan_elem.first) << " " << plan_elem.second << endl;
    }
#endif

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

const IndexRecipe& IndexRegistry::get_recipe(const RecipeName& recipe_name) const {
    const auto& recipes = recipe_registry.at(recipe_name.first);
    assert(recipe_name.second < recipes.size());
    return recipes.at(recipe_name.second);
}

vector<vector<string>> IndexRegistry::execute_recipe(const RecipeName& recipe_name, const IndexingPlan* plan,
                                                     AliasGraph& alias_graph) {
    const auto& index_recipe = get_recipe(recipe_name);
    if (recipe_name.first.size() > 1 || !index_recipe.input_group().count(*recipe_name.first.begin())) {
        // we're not in an unboxing recipe (in which case not all of the indexes might have been
        // unboxed yet, in which case they appear unfinished)
        for (auto input : index_recipe.inputs) {
            assert(input->is_finished());
        }
    }
#ifdef debug_index_registry_recipes
    cerr << "executing recipe " << recipe_name.second << " for " << to_string(recipe_name.first) << endl;
#endif
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
    map<RecipeName, string> recipe_to_dot_id;
    for (const auto& recipe_record : recipe_registry) {
        const auto& recipes = recipe_record.second;
        for (size_t priority_idx = 0; priority_idx < recipes.size(); ++priority_idx, ++recipe_idx) {
            const auto& recipe = recipes[priority_idx];
            string recipe_dot_id = "R" + to_string(recipe_idx);
            recipe_to_dot_id[RecipeName(recipe_record.first, priority_idx)] = recipe_dot_id;
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
    for (const auto& generalization_record : generalizations) {
        strm << recipe_to_dot_id.at(generalization_record.first) << " -> " << recipe_to_dot_id.at(generalization_record.second) << " [style=dashed color=" << unselected_col << "];" << endl;
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
                                                       const IndexRegistry& registry) noexcept :
    runtime_error("Insufficient input to create " + target), target(target), inputs(registry.completed_indexes())
{
    // nothing else to do
    stringstream ss;
    ss << "Inputs" << endl;
    for (const auto& input : inputs) {
        ss << "\t" << input << endl;
    }
    ss << "are insufficient to create target index " << target << endl;
    msg = ss.str();
}

const char* InsufficientInputException::what() const noexcept {
    return msg.c_str();
}


RewindPlanException::RewindPlanException(const string& msg, const IndexGroup& rewind_to) noexcept : msg(msg), indexes(rewind_to) {
    // nothing else to do
}

const char* RewindPlanException::what() const noexcept {
    return msg.c_str();
}

const IndexGroup& RewindPlanException::get_indexes() const noexcept {
    return indexes;
}

}

