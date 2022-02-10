/** \file sim_main.cpp
 *
 * Defines the "vg sim" subcommand, which generates potential reads from a graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>
#include <algorithm>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../aligner.hpp"
#include "../gbwt_helper.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "../sampler.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

// Gets the transcript IDs and TPM values from an RSEM output .tsv file
vector<pair<string, double>> parse_rsem_expression_file(istream& rsem_in) {
    vector<pair<string, double>> return_val;
    string line;
    // skip the header line
    getline(rsem_in, line);
    line.clear();
    while (getline(rsem_in, line)) {
        vector<string> tokens;
        stringstream strm(line);
        string token;
        while (getline(strm, token, '\t')) {
            tokens.push_back(move(token));
            token.clear();
        }
        if (tokens.size() != 8) {
            cerr << "[vg sim] error: Cannot parse transcription file. Expected 8-column TSV file as produced by RSEM, got " << tokens.size() << " columns." << endl;
            exit(1);
        }
        return_val.emplace_back(tokens[0], parse<double>(tokens[5]));
        line.clear();
    }
    return return_val;
}

// Gets the trancript path name, the original transcript name, and the haplotype count from the vg rna -i file
vector<tuple<string, string, size_t>> parse_haplotype_transcript_file(istream& haplo_tx_in) {
    vector<tuple<string, string, size_t>> return_val;
    string line;
    // skip the header line
    getline(haplo_tx_in, line);
    line.clear();
    while (getline(haplo_tx_in, line)) {
        vector<string> tokens;
        stringstream strm(line);
        string token;
        while (getline(strm, token, '\t')) {
            tokens.push_back(move(token));
            token.clear();
        }
        if (tokens.size() != 5) {
            cerr << "[vg sim] error: Cannot parse haplotype transcript file. Expected 5-column TSV file as produced by vg rna -i, got " << tokens.size() << " columns." << endl;
            exit(1);
        }
        // contributing haplotypes are separeted by commas
        size_t haplo_count = 1 + std::count(tokens[4].begin(), tokens[4].end(), ',');
        return_val.emplace_back(tokens[0], tokens[2], haplo_count);
        line.clear();
    }
    return return_val;
}

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options]" << endl
         << "Samples sequences from the xg-indexed graph." << endl
         << endl
         << "basic options:" << endl
         << "    -x, --xg-name FILE          use the graph in FILE (required)" << endl
         << "    -n, --num-reads N           simulate N reads or read pairs" << endl
         << "    -l, --read-length N         simulate reads of length N" << endl
         << "    -r, --progress              show progress information" << endl
         << "output options:" << endl
         << "    -a, --align-out             write alignments in GAM-format" << endl
         << "    -J, --json-out              write alignments in json" << endl
         << "simulation parameters:" << endl
         << "    -F, --fastq FILE            match the error profile of NGS reads in FILE, repeat for paired reads (ignores -l,-f)" << endl
         << "    -I, --interleaved           reads in FASTQ (-F) are interleaved read pairs" << endl
         << "    -s, --random-seed N         use this specific seed for the PRNG" << endl
         << "    -e, --sub-rate FLOAT        base substitution rate (default 0.0)" << endl
         << "    -i, --indel-rate FLOAT      indel rate (default 0.0)" << endl
         << "    -d, --indel-err-prop FLOAT  proportion of trained errors from -F that are indels (default 0.0)" << endl
         << "    -S, --scale-err FLOAT       scale trained error probabilities from -F by this much (default 1.0)" << endl
         << "    -f, --forward-only          don't simulate from the reverse strand" << endl
         << "    -p, --frag-len N            make paired end reads with given fragment length N" << endl
         << "    -v, --frag-std-dev FLOAT    use this standard deviation for fragment length estimation" << endl
         << "    -N, --allow-Ns              allow reads to be sampled from the graph with Ns in them" << endl
         << "    -t, --threads               number of compute threads (only when using FASTQ with -F) [1]" << endl
         << "simulate from paths:" << endl
         << "    -P, --path PATH             simulate from this path (may repeat; cannot also give -T)" << endl
         << "    -A, --any-path              simulate from any path (overrides -P)" << endl
         << "    -m, --sample-name NAME      simulate from this sample (may repeat; requires -g)" << endl
         << "    -R, --ploidy-regex RULES    use the given comma-separated list of colon-delimited REGEX:PLOIDY rules to assign" << endl
         << "                                ploidies to contigs not visited by the selected samples, or to all contigs simulated" << endl
         << "                                from if no samples are used. Unmatched contigs get ploidy 2." << endl
         << "    -g, --gbwt-name FILE        use samples from this GBWT index" << endl
         << "    -T, --tx-expr-file FILE     simulate from an expression profile formatted as RSEM output (cannot also give -P)" << endl
         << "    -H, --haplo-tx-file FILE    transcript origin info table from vg rna -i (required for -T on haplotype transcripts)" << endl
         << "    -u, --unsheared             sample from unsheared fragments" << endl
         << "    -E, --path-pos-file FILE    output a TSV with sampled position on path of each read (requires -F)" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    string xg_name;
    int num_reads = 1;
    int read_length = 100;
    bool progress = false;
    int threads = 1;

    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;
    bool forward_only = false;
    bool align_out = false;
    bool json_out = false;
    int fragment_length = 0;
    double fragment_std_dev = 0;
    bool reads_may_contain_Ns = false;
    bool strip_bonuses = false;
    bool interleaved = false;
    bool unsheared_fragments = false;
    double indel_prop = 0.0;
    double error_scale_factor = 1.0;
    string fastq_name;
    string fastq_2_name;
    string path_pos_filename;

    // What path should we sample from? Empty string = the whole graph.
    vector<string> path_names;
    bool any_path = false;

    // Sample from GBWT threads.
    std::vector<std::string> sample_names;
    std::string gbwt_name;
    
    // When sampling from paths or GBWT threads, what ploidy should we assign to each path?
    // Represented as a list of regexes (to match the whole path name) and ploidies.
    // The first rule to match wins.
    // When using GBWT threads, only applies to contigs with no threads in any sample.
    // Each thread that does exist is ploidy 1.
    std::vector<std::pair<std::regex, double>> ploidy_rules;

    // Alternatively, which transcripts with how much expression?
    string rsem_file_name;
    vector<pair<string, double>> transcript_expressions;
    // If we made haplotype trancripts, we'll need a translation layer onto the
    // expression profile
    string haplotype_transcript_file_name;
    vector<tuple<string, string, size_t>> haplotype_transcripts;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"progress", no_argument, 0, 'r'},
            {"fastq", required_argument, 0, 'F'},
            {"interleaved", no_argument, 0, 'I'},
            {"path", required_argument, 0, 'P'},
            {"any-path", no_argument, 0, 'A'},
            {"sample-name", required_argument, 0, 'm'},
            {"ploidy-regex", required_argument, 0, 'R'},
            {"gbwt-name", required_argument, 0, 'g'},
            {"tx-expr-file", required_argument, 0, 'T'},
            {"haplo-tx-file", required_argument, 0, 'H'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {"json-out", no_argument, 0, 'J'},
            {"allow-Ns", no_argument, 0, 'N'},
            {"unsheared", no_argument, 0, 'u'},
            {"sub-rate", required_argument, 0, 'e'},
            {"indel-rate", required_argument, 0, 'i'},
            {"indel-err-prop", required_argument, 0, 'd'},
            {"scale-err", required_argument, 0, 'S'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {"path-usage", required_argument, 0, 'E'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hrl:n:s:e:i:fax:Jp:v:Nud:F:P:Am:R:g:T:H:S:It:E:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;
            
        case 'r':
            progress = true;
            break;

        case 'F':
            if (fastq_name.empty()) {
                fastq_name = optarg;
            }
            else if (fastq_2_name.empty()) {
                fastq_2_name = optarg;
            }
            else {
                cerr << "error: cannot provide more than 2 FASTQs to train simulator" << endl;
                exit(1);
            }
            break;
            
        case 'I':
            interleaved = true;
            break;
            
        case 'P':
            path_names.push_back(optarg);
            break;

        case 'A':
            any_path = true;
            break;

        case 'm':
            sample_names.push_back(optarg);
            break;
            
        case 'R':
            for (auto& rule : split_delims(optarg, ",")) {
                // For each comma-separated rule
                auto parts = split_delims(rule, ":");
                if (parts.size() != 2) {
                    cerr << "error: ploidy rules must be REGEX:PLOIDY" << endl;
                    exit(1);
                }
                try {
                    // Parse the regex
                    std::regex match(parts[0]);
                    double weight = parse<double>(parts[1]);
                    // Save the rule
                    ploidy_rules.emplace_back(match, weight);
                } catch (const std::regex_error& e) {
                    // This is not a good regex
                    cerr << "error: unacceptable regular expression \"" << parts[0] << "\": " << e.what() << endl;
                    exit(1);
                }
            }
            break;

        case 'g':
            gbwt_name = optarg;
            break;

        case 'T':
            rsem_file_name = optarg;
            break;
                
        case 'H':
            haplotype_transcript_file_name = optarg;
            break;

        case 'l':
            read_length = parse<int>(optarg);
            break;

        case 'n':
            num_reads = parse<int>(optarg);
            break;

        case 's':
            seed_val = parse<int>(optarg);
            if (seed_val == 0) {
                // Don't let the user specify seed 0 as we will confuse it with no deterministic seed.
                cerr << "error[vg sim]: seed 0 cannot be used. Omit the seed option if you want nondeterministic results." << endl;
                exit(1);
            }
            break;

        case 'e':
            base_error = parse<double>(optarg);
            break;

        case 'i':
            indel_error = parse<double>(optarg);
            break;
            
        case 'd':
            indel_prop = parse<double>(optarg);
            break;
            
        case 'S':
            error_scale_factor = parse<double>(optarg);
            break;

        case 'f':
            forward_only = true;
            break;

        case 'a':
            align_out = true;
            break;

        case 'J':
            json_out = true;
            align_out = true;
            break;

        case 'N':
            reads_may_contain_Ns = true;
            break;
                
        case 'u':
            unsheared_fragments = true;
            break;

        case 'p':
            fragment_length = parse<int>(optarg);
            break;

        case 'v':
            fragment_std_dev = parse<double>(optarg);
            break;
                
        case 't':
            threads = parse<int>(optarg);
            break;
                
        case 'E':
            path_pos_filename = optarg;
            break;
            
        case 'h':
        case '?':
            help_sim(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }
    
    omp_set_num_threads(threads);
    
    // We'll fill this in with ploidies for each path in path_names
    std::vector<double> path_ploidies;
    // When we need to consult the ploidy rules about a contig nemr we call this function.
    auto consult_ploidy_rules = [&](const std::string& name) {
        for (auto& rule : ploidy_rules) {
            if (std::regex_match(name, rule.first)) {
                // This rule should apply to this contig
                return rule.second;
            }
        }
        // Unmatched contigs get ploidy 2.
        // 1 makes no sense in the context of a genomic reference.
        // 0 makes no sense for --all-paths which consults the rules for all the names.
        return 2.0;
    };

    if (xg_name.empty()) {
        cerr << "[vg sim] error: we need a graph to sample reads from" << endl;
        return 1;
    }
    if (!gbwt_name.empty() && sample_names.empty() && rsem_file_name.empty()) {
        cerr << "[vg sim] error: --gbwt-name requires --sample-name or --tx-expr-file" << endl;
        return 1;
    }
    if (!sample_names.empty() && gbwt_name.empty()) {
        cerr << "[vg sim] error: --sample-name must be used with --gbwt-name" << endl;
        return 1;
    }
    if (!gbwt_name.empty() && !rsem_file_name.empty() && !haplotype_transcript_file_name.empty()) {
        cerr << "[vg sim] error: using --gbwt-name requires that HSTs be included --tx-expr-file, combination with --haplo-tx-file is not implemented" << endl;
        return 1;
    }

    if (!rsem_file_name.empty()) {
        if (progress) {
            std::cerr << "Reading transcription profile from " << rsem_file_name << std::endl;
        }
        ifstream rsem_in(rsem_file_name);
        if (!rsem_in) {
            cerr << "[vg sim] error: could not open transcription profile file " << rsem_file_name << endl;
            return 1;
        }
        transcript_expressions = parse_rsem_expression_file(rsem_in);
    }
    
    if (!haplotype_transcript_file_name.empty()) {
        if (progress) {
            std::cerr << "Reading haplotype transcript file " << haplotype_transcript_file_name << std::endl;
        }
        ifstream haplo_tx_in(haplotype_transcript_file_name);
        if (!haplo_tx_in) {
            cerr << "[vg sim] error: could not open haplotype transcript file " << haplotype_transcript_file_name << endl;
            return 1;
        }
        haplotype_transcripts = parse_haplotype_transcript_file(haplo_tx_in);
    }

    if (progress) {
        std::cerr << "Loading graph " << xg_name << std::endl;
    }
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    
    if (!path_pos_filename.empty() && fastq_name.empty()) {
        cerr << "[vg sim] error: path usage table is not available unless using trained simulation (-F)" << endl;
        exit(1);
    }
    
    // Deal with path names. Do this before we create paths to represent threads.
    if (any_path) {
        if (progress) {
            std::cerr << "Selecting all " << path_handle_graph->get_path_count() << " paths" << std::endl;
        }
        if (path_handle_graph->get_path_count() == 0) {
            cerr << "[vg sim] error: the graph does not contain paths" << endl;
            return 1;
        }
        path_names.clear();
        path_handle_graph->for_each_path_handle([&](const path_handle_t& handle) {
            // For each path in the graph
            auto name = path_handle_graph->get_path_name(handle);
            // Simulate from it
            path_names.push_back(name);
            // At ploidy defined by the rules (default 2)
            path_ploidies.push_back(consult_ploidy_rules(name));
        });
    } else if (!path_names.empty()) {
        if (progress) {
            std::cerr << "Checking " << path_names.size() << " selected paths" << std::endl;
        }
        for (auto& path_name : path_names) {
            if (path_handle_graph->has_path(path_name) == false) {
                cerr << "[vg sim] error: path \""<< path_name << "\" not found in index" << endl;
                return 1;
            }
            // Synthesize ploidies for explicitly specified paths
            path_ploidies.push_back(consult_ploidy_rules(path_name));
        }
    }

    // Deal with GBWT threads
    if (!gbwt_name.empty()) {
        
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);
        if (!(gbwt_index->hasMetadata()) || !(gbwt_index->metadata.hasSampleNames()) || !(gbwt_index->metadata.hasPathNames())) {
            std::cerr << "[vg sim] error: GBWT index does not contain sufficient metadata" << std::endl;
            return 1;
        }
        
        // we will add these threads to the graph as named paths and index them for easy look up
        hash_map<gbwt::size_type, size_t> sample_id_to_idx;
        
        // this is used to keep track of named paths that we want to simulate from, even if they
        // don't have sample information
        std::unordered_set<std::string> unvisited_paths;
        
        if (!sample_names.empty()) {
            // we're consulting the provided sample names to determine which threads to include
            
            // We need to track the contigs that have not had any threads in any sample
            if (!ploidy_rules.empty()) {
                // We actually want to visit them, so we have to find them
                if (progress) {
                    std::cerr << "Inventorying contigs" << std::endl;
                }
                path_handle_graph->for_each_path_handle([&](const path_handle_t& handle) {
                    // For each path in the graph
                    auto name = path_handle_graph->get_path_name(handle);
                    if (!Paths::is_alt(name)) {
                        // TODO: We assume that if it isn't an alt path it represents a contig!
                        // TODO: We may need to change this when working with graphs with multiple sets of primary paths, or other extra paths.
                        unvisited_paths.insert(name);
                    }
                });
            }
            
            if (progress) {
                std::cerr << "Checking " << sample_names.size() << " samples" << std::endl;
            }
            
            for (std::string& sample_name : sample_names) {
                gbwt::size_type id = gbwt_index->metadata.sample(sample_name);
                if (id >= gbwt_index->metadata.samples()) {
                    std::cerr << "[vg sim] error: sample \"" << sample_name << "\" not found in the GBWT index" << std::endl;
                    return 1;
                }
                sample_id_to_idx[id] = sample_id_to_idx.size();
            }
        }
        else {
            // we are consulting the transcript expression table to decide which threads to include
            for (const auto& transcript_expression  : transcript_expressions) {
                gbwt::size_type id = gbwt_index->metadata.sample(transcript_expression.first);
                if (id >= gbwt_index->metadata.samples()) {
                    std::cerr << "[vg sim] error: haplotype-specific transcript \"" << transcript_expression.first << "\" not found in the GBWT index" << std::endl;
                    return 1;
                }
                auto idx = sample_id_to_idx.size();
                sample_id_to_idx[id] = idx;
            }
        }
        
        MutablePathMutableHandleGraph* mutable_graph = dynamic_cast<MutablePathMutableHandleGraph*>(path_handle_graph.get());
        if (mutable_graph == nullptr) {
            if (progress) {
                std::cerr << "Converting the graph into HashGraph" << std::endl;
            }
            mutable_graph = new bdsg::HashGraph();
            handlealgs::copy_path_handle_graph(path_handle_graph.get(), mutable_graph);
            path_handle_graph.reset(mutable_graph);
        }
        if (progress) {
            std::cerr << "Inserting " << sample_id_to_idx.size() << " GBWT threads into the graph" << std::endl;
        }
        
        size_t inserted = 0;
        for (gbwt::size_type i = 0; i < gbwt_index->metadata.paths(); i++) {
            auto& path = gbwt_index->metadata.path(i);
            auto it = sample_id_to_idx.find(path.sample);
            if (it != sample_id_to_idx.end()) {
                std::string path_name = insert_gbwt_path(*mutable_graph, *gbwt_index, i);
                if (!path_name.empty()) {
                    if (!sample_names.empty()) {
                        // assign this haplotype a single ploidy
                        
                        // We managed to make a path for this thread
                        path_names.push_back(path_name);
                        // It should have ploidy 1
                        path_ploidies.push_back(1.0);
                        
                        if (!unvisited_paths.empty()) {
                            // Remember that the contig this path is on is visited
                            auto contig_name = gbwt_index->metadata.contig(path.contig);
                            unvisited_paths.erase(contig_name);
                        }
                    }
                    else {
                        // update the transcript name so we can assign it expression
                        // later down
                        transcript_expressions[it->second].first = path_name;
                    }
                    // Remember we inserted a path
                    inserted++;
                }
            }
        }
        if (progress) {
            std::cerr << "Inserted " << inserted << " paths" << std::endl;
        }
        if (!unvisited_paths.empty()) {
            // There are unvisited contigs we want to sample from too
            for (auto& name : unvisited_paths) {
                // Sample from each
                path_names.push_back(name);
                // With the rule-determined ploidy
                path_ploidies.push_back(consult_ploidy_rules(name));
            }
            if (progress) {
                std::cerr << "Also sampling from " << unvisited_paths.size() << " paths representing unvisited contigs" << std::endl;
            }
        }
    }
    
    
    if (haplotype_transcript_file_name.empty()) {
        if (progress && !transcript_expressions.empty()) {
            std::cerr << "Checking " << transcript_expressions.size() << " transcripts" << std::endl;
        }
        for (auto& transcript_expression : transcript_expressions) {
            if (!path_handle_graph->has_path(transcript_expression.first)) {
                cerr << "[vg sim] error: transcript path for \""<< transcript_expression.first << "\" not found in index" << endl;
                cerr << "if you embedded haplotype-specific transcripts in the graph, you may need the haplotype transcript file from vg rna -i" << endl;
                return 1;
            }
        }
    }
    else {
        if (progress) {
            std::cerr << "Checking " << haplotype_transcripts.size() << " haplotype transcripts" << std::endl;
        }
        for (auto& haplotype_transcript : haplotype_transcripts) {
            if (!path_handle_graph->has_path(get<0>(haplotype_transcript))) {
                cerr << "[vg sim] error: transcript path for \""<< get<0>(haplotype_transcript) << "\" not found in index" << endl;
                return 1;
            }
        }
    }
    
    if (progress) {
        std::cerr << "Creating path position overlay" << std::endl;
    }
    bdsg::PathPositionVectorizableOverlayHelper overlay_helper;
    PathPositionHandleGraph* xgidx = dynamic_cast<PathPositionHandleGraph*>(overlay_helper.apply(path_handle_graph.get()));
    
    unique_ptr<vg::io::ProtobufEmitter<Alignment>> aln_emitter;
    if (align_out && !json_out) {
        // Make an emitter to emit Alignments
        aln_emitter = unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
    }

    if (progress) {
        std::cerr << "Simulating " << (fragment_length > 0 ? "read pairs" : "reads") << std::endl;
        std::cerr << "--num-reads " << num_reads << std::endl;
        std::cerr << "--read-length " << read_length << std::endl;
        if (align_out) {
            std::cerr << "--align-out" << std::endl;
        }
        if (json_out) {
            std::cerr << "--json-out" << std::endl;
        }
        if (!fastq_name.empty()) {
            std::cerr << "--fastq " << fastq_name << std::endl;
            if (!fastq_2_name.empty()) {
                std::cerr << "--fastq " << fastq_2_name << std::endl;
            }
            if (interleaved) {
                std::cerr << "--interleaved" << std::endl;
            }
        } else {
            if (base_error > 0.0) {
                std::cerr << "--sub-rate " << base_error << std::endl;
            }
        }
        if (indel_error > 0.0) {
            std::cerr << "--indel-rate " << indel_error << std::endl;
        }
        if (!fastq_name.empty()) {
            if (indel_prop > 0.0) {
                std::cerr << "--indel-err-prop " << indel_prop << std::endl;
            }
            if (error_scale_factor != 1.0) {
                std::cerr << "--scale-err " << error_scale_factor << std::endl;
            }
        }
        if (forward_only) {
            std::cerr << "--forward-only" << std::endl;
        }
        if (fragment_length > 0) {
            std::cerr << "--frag-len " << fragment_length << std::endl;
            if (fragment_std_dev > 0.0) {
                std::cerr << "--frag-std-dev " << fragment_std_dev << std::endl;
            }
        }
        if (reads_may_contain_Ns) {
            std::cerr << "--allow-Ns" << std::endl;
        }
    }

    if (fastq_name.empty()) {
        // Use the fixed error rate sampler
        
        if (unsheared_fragments) {
            cerr << "Unsheared fragment option only available when simulating from FASTQ-trained errors" << endl;
        }
        
        // Make a sample to sample reads with
        Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns, path_names, path_ploidies, transcript_expressions, haplotype_transcripts);
        
        // initialize an aligner
        Aligner rescorer(default_score_matrix, default_gap_open, default_gap_extension,
                         default_full_length_bonus, vg::default_gc_content);

        // We define a function to score a using the aligner
        auto rescore = [&] (Alignment& aln) {
            // Score using exact distance.
            aln.set_score(rescorer.score_contiguous_alignment(aln, strip_bonuses));
        };
        
        size_t max_iter = 1000;
        int nonce = 1;
        for (int i = 0; i < num_reads; ++i) {
            // For each read we are going to generate
            
            if (fragment_length) {
                // fragment_lenght is nonzero so make it two paired reads
                auto alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (alns.front().sequence().size() < read_length
                        || alns.back().sequence().size() < read_length) {
                        // If our read was too short, try again
                        alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(alns.front());
                    rescore(alns.back());
                    
                    if (json_out) {
                        cout << pb2json(alns.front()) << endl;
                        cout << pb2json(alns.back()) << endl;
                    } else {
                        aln_emitter->write_copy(alns.front());
                        aln_emitter->write_copy(alns.back());
                    }
                } else {
                    cout << alns.front().sequence() << "\t" << alns.back().sequence() << endl;
                }
            } else {
                // Do single-end reads
                auto aln = sampler.alignment_with_error(read_length, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (aln.sequence().size() < read_length) {
                        // If our read is too short, try again
                        auto aln_prime = sampler.alignment_with_error(read_length, base_error, indel_error);
                        if (aln_prime.sequence().size() > aln.sequence().size()) {
                            // But only keep the new try if it is longer
                            aln = aln_prime;
                        }
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(aln);
                    
                    if (json_out) {
                        cout << pb2json(aln) << endl;
                    } else {
                        aln_emitter->write_copy(aln);
                    }
                } else {
                    cout << aln.sequence() << endl;
                }
            }
        }
        
    }
    else {
        // Use the trained error rate
        
        Aligner aligner(default_score_matrix, default_gap_open, default_gap_extension,
                         default_full_length_bonus, vg::default_gc_content);
        
        NGSSimulator sampler(*xgidx,
                             fastq_name,
                             fastq_2_name,
                             interleaved,
                             path_names,
                             path_ploidies,
                             transcript_expressions,
                             haplotype_transcripts,
                             base_error,
                             indel_error,
                             indel_prop,
                             fragment_length ? fragment_length : std::numeric_limits<double>::max(), // suppresses warnings about fragment length
                             fragment_std_dev ? fragment_std_dev : 0.000001, // eliminates errors from having 0 as stddev without substantial difference
                             error_scale_factor,
                             !reads_may_contain_Ns,
                             unsheared_fragments,
                             seed_val);
        
        if (!path_pos_filename.empty()) {
            sampler.connect_to_position_file(path_pos_filename);
        }
        
        unique_ptr<AlignmentEmitter> alignment_emitter = get_non_hts_alignment_emitter("-", json_out ? "JSON" : "GAM",
                                                                                       map<string, int64_t>(), get_thread_count());
        
        // static scheduling could produce some degradation in speed, but i think it should make
        // the output deterministic (except for ordering)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < num_reads; i++) {
            if (fragment_length) {
                pair<Alignment, Alignment> read_pair = sampler.sample_read_pair();
                read_pair.first.set_score(aligner.score_contiguous_alignment(read_pair.first, strip_bonuses));
                read_pair.second.set_score(aligner.score_contiguous_alignment(read_pair.second, strip_bonuses));
                
                if (align_out) {
                    alignment_emitter->emit_pair(std::move(read_pair.first), std::move(read_pair.second));
                }
                else {
                    string line = read_pair.first.sequence() + '\t' + read_pair.second.sequence() + '\n';
#pragma omp critical
                    cout << line;
                }
            }
            else {
                Alignment read = sampler.sample_read();
                read.set_score(aligner.score_contiguous_alignment(read, strip_bonuses));
                
                if (align_out) {
                    alignment_emitter->emit_single(std::move(read));
                }
                else {
                    string line = read.sequence() + '\n';
#pragma omp critical
                    cout << line;
                }
            }
        }
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_sim("sim", "simulate reads from a graph", TOOLKIT, main_sim);

