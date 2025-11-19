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
vector<pair<string, double>> parse_rsem_expression_file(const Logger& logger, istream& rsem_in) {
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
            tokens.push_back(std::move(token));
            token.clear();
        }
        if (tokens.size() != 8) {
            logger.error() << "Cannot parse transcription file. "
                           << "Expected 8-column TSV file as produced by RSEM, got "
                           << tokens.size() << " columns." << std::endl;
        }
        return_val.emplace_back(tokens[0], parse<double>(tokens[5]));
        line.clear();
    }
    return return_val;
}

// Gets the trancript path name, the original transcript name, and the haplotype count from the vg rna -i file
vector<tuple<string, string, size_t>> parse_haplotype_transcript_file(const Logger& logger, istream& haplo_tx_in) {
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
            tokens.push_back(std::move(token));
            token.clear();
        }
        if (tokens.size() != 5) {
            logger.error() << "Cannot parse haplotype transcript file. "
                           << "Expected 5-column TSV file as produced by vg rna -i, got "
                           << tokens.size() << " columns." << std::endl;
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
         << "Samples sequences from an XG-indexed graph." << endl
         << endl
         << "basic options:" << endl
         << "  -h, --help                  print this help message to stderr and exit" << endl
         << "  -x, --xg-name FILE          use the graph in FILE (required)" << endl
         << "  -n, --num-reads N           simulate N reads or read pairs" << endl
         << "  -l, --read-length N         simulate reads of length N" << endl
         << "  -r, --progress              show progress information" << endl
         << "output options:" << endl
         << "  -a, --align-out             write alignments in GAM-format" << endl
         << "  -q, --fastq-out             write reads in FASTQ format" << endl
         << "  -J, --json-out              write alignments in JSON-format GAM (implies -a)" << endl
         << "      --multi-position        annotate with multiple reference positions" << endl
         << "simulation parameters:" << endl
         << "  -F, --fastq FILE            match the error profile of NGS reads in FILE," << endl
         << "                              repeat for paired reads (ignores -l,-f)" << endl
         << "  -I, --interleaved           reads in FASTQ (-F) are interleaved read pairs" << endl
         << "  -s, --random-seed N         use this specific seed for the PRNG" << endl
         << "  -e, --sub-rate FLOAT        base substitution rate [0.0]" << endl
         << "  -i, --indel-rate FLOAT      indel rate [0.0]" << endl
         << "  -d, --indel-err-prop FLOAT  proportion of trained errors from -F" << endl
         << "                              that are indels [0.01]" << endl
         << "  -S, --scale-err FLOAT       scale trained error probs from -F by FLOAT [1.0]" << endl
         << "  -f, --forward-only          don't simulate from the reverse strand" << endl
         << "  -p, --frag-len N            make paired end reads with fragment length N" << endl
         << "  -v, --frag-std-dev FLOAT    use this standard deviation" << endl
         << "                              for fragment length estimation" << endl
         << "  -N, --allow-Ns              allow reads to be sampled with Ns in them" << endl
         << "      --max-tries N           attempt sampling operations up to N times [100]" << endl
         << "  -t, --threads N             number of compute threads (only when using -F) [1]" << endl
         << "simulate from paths:" << endl
         << "  -P, --path NAME             simulate from this path" << endl
         << "                              (may repeat; cannot also give -T)" << endl
         << "  -A, --any-path              simulate from any path (overrides -P)" << endl
         << "  -m, --sample-name NAME      simulate from this sample (may repeat)" << endl
         << "  -R, --ploidy-regex RULES    use this comma-separated list of colon-delimited" << endl
         << "                              REGEX:PLOIDY rules to assign ploidies to contigs" << endl
         << "                              not visited by the selected samples, or to all" << endl
         << "                              contigs simulated from if no samples are used." << endl
         << "                              Unmatched contigs get ploidy 2" << endl
         << "  -g, --gbwt-name FILE        use samples from this GBWT index" << endl
         << "  -T, --tx-expr-file FILE     simulate from an expression profile formatted as" << endl
         << "                              RSEM output (cannot also give -P)" << endl
         << "  -H, --haplo-tx-file FILE    transcript origin info table from vg rna -i" << endl
         << "                              (required for -T on haplotype transcripts)" << endl
         << "  -u, --unsheared             sample from unsheared fragments" << endl
         << "  -E, --path-pos-file FILE    output a TSV with sampled position on path" << endl
         << "                              of each read (requires -F)" << endl;
}

int main_sim(int argc, char** argv) {
    Logger logger("vg sim");

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    constexpr int OPT_MULTI_POSITION = 1000;
    constexpr int OPT_MAX_TRIES = 1001;

    string xg_name;
    int num_reads = 1;
    int read_length = 100;
    bool progress = false;

    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;
    bool forward_only = false;
    bool align_out = false;
    bool json_out = false;
    bool fastq_out = false;
    bool multi_position_annotations = false;
    int fragment_length = 0;
    double fragment_std_dev = 0;
    bool reads_may_contain_Ns = false;
    size_t max_tries = 100;
    bool strip_bonuses = false;
    bool interleaved = false;
    bool unsheared_fragments = false;
    double indel_prop = 0.01;
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
            {"fastq-out", no_argument, 0, 'q'},
            {"json-out", no_argument, 0, 'J'},
            {"multi-position", no_argument, 0, OPT_MULTI_POSITION},
            {"allow-Ns", no_argument, 0, 'N'},
            {"max-tries", required_argument, 0, OPT_MAX_TRIES},
            {"unsheared", no_argument, 0, 'u'},
            {"sub-rate", required_argument, 0, 'e'},
            {"indel-rate", required_argument, 0, 'i'},
            {"indel-err-prop", required_argument, 0, 'd'},
            {"scale-err", required_argument, 0, 'S'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {"path-pos-file", required_argument, 0, 'E'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?rl:n:s:e:i:fax:qJp:v:Nud:F:P:Am:R:g:T:H:S:It:E:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = require_exists(logger, optarg);
            break;
            
        case 'r':
            progress = true;
            break;

        case 'F':
            assign_fastq_files(logger, optarg, fastq_name, fastq_2_name);
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
                    logger.error() << "ploidy rules must be REGEX:PLOIDY" << std::endl;
                }
                try {
                    // Parse the regex
                    std::regex match(parts[0]);
                    double weight = parse<double>(parts[1]);
                    // Save the rule
                    ploidy_rules.emplace_back(match, weight);
                } catch (const std::regex_error& e) {
                    // This is not a good regex
                    logger.error() << "unacceptable regular expression\n\""
                                   << parts[0] << "\":\n" << e.what() << std::endl;
                }
            }
            break;

        case 'g':
            gbwt_name = require_exists(logger, optarg);
            break;

        case 'T':
            rsem_file_name = require_exists(logger, optarg);
            break;
                
        case 'H':
            haplotype_transcript_file_name = require_exists(logger, optarg);
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
                logger.error() << "seed 0 cannot be used. Omit the seed option "
                               << "if you want nondeterministic results." << std::endl;
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

        case 'q':
            fastq_out = true;
            break;

        case 'J':
            json_out = true;
            align_out = true;
            break;
            
        case OPT_MULTI_POSITION:
            multi_position_annotations = true;
            break;

        case 'N':
            reads_may_contain_Ns = true;
            break;
            
        case OPT_MAX_TRIES:
            max_tries = parse<size_t>(optarg);
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
            set_thread_count(logger, optarg);
            break;
                
        case 'E':
            path_pos_filename = require_exists(logger, optarg);
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

    if (align_out && fastq_out) {
        logger.error() << "only one output format (-a/-J or -q) can be selected." << std::endl;
    }

    if (xg_name.empty()) {
        logger.error() << "we need a graph to sample reads from" << std::endl;
    }
    if (!gbwt_name.empty() && sample_names.empty() && rsem_file_name.empty()) {
        logger.error() << "--gbwt-name requires --sample-name or --tx-expr-file" << std::endl;
    }
    if (!gbwt_name.empty() && !rsem_file_name.empty() && !haplotype_transcript_file_name.empty()) {
        // TODO: This message doesn't really make sense.
        logger.error() << "using --gbwt-name requires that HSTs be included --tx-expr-file; "
                       << "combination with --haplo-tx-file is not implemented" << std::endl;
    }

    if (!rsem_file_name.empty()) {
        if (progress) {
            logger.info() << "Reading transcription profile from "
                          << rsem_file_name << std::endl;
        }
        ifstream rsem_in(rsem_file_name);
        transcript_expressions = parse_rsem_expression_file(logger, rsem_in);
    }
    
    if (!haplotype_transcript_file_name.empty()) {
        if (progress) {
            logger.info() << "Reading haplotype transcript file "
                          << haplotype_transcript_file_name << std::endl;
        }
        ifstream haplo_tx_in(haplotype_transcript_file_name);
        haplotype_transcripts = parse_haplotype_transcript_file(logger, haplo_tx_in);
    }

    if (progress) {
        logger.info() << "Loading graph " << xg_name << std::endl;
    }
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    
    if (!path_pos_filename.empty() && fastq_name.empty()) {
        logger.error() << "path usage table is not available unless using trained simulation (-F)" << std::endl;
    }
    
    if (fastq_name.empty() && unsheared_fragments) {
        logger.error() << "unsheared fragment option only available "
                       << "when simulating from FASTQ-trained errors" << std::endl;
    }
    
    // Deal with path names. Do this before we create paths to represent threads.
    if (any_path) {
        if (progress) {
            logger.info() << "Selecting all " << path_handle_graph->get_path_count()
                          << " paths" << std::endl;
        }
        if (path_handle_graph->get_path_count() == 0) {
            logger.error() << "the graph does not contain paths" << std::endl;
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
            logger.info() << "Checking " << path_names.size() << " selected paths" << std::endl;
        }
        for (auto& path_name : path_names) {
            if (path_handle_graph->has_path(path_name) == false) {
                logger.error() << "path \"" << path_name << "\" not found in index" << std::endl;
            }
            // Synthesize ploidies for explicitly specified paths
            path_ploidies.push_back(consult_ploidy_rules(path_name));
        }
    }
    
    // We may add some paths to our graph. If so, we need to ignore them when
    // annotating with path positions, because they will be useless.
    unordered_set<string> inserted_path_names;
   
    // We need to track the contigs that have generic sense paths in the graph,
    // but don't have any paths with that name as the locus in any selected
    // sample. If ploidy rules are used, we will simulate from them as well, at
    // those ploidies.
    std::unordered_set<std::string> unvisited_contigs;
    if (!sample_names.empty() && !ploidy_rules.empty()) {
        // We need to track the contigs that have not had any threads in any sample.
        
        // We actually want to visit them, so we have to find them
        if (progress) {
            logger.info() << "Inventorying contigs" << std::endl;
        }
        // We assume the generic paths in the graph are contigs ("chr1", "chr2", etc.)
        path_handle_graph->for_each_path_of_sense(PathSense::GENERIC, [&](const path_handle_t& handle) {
            // For each path in the graph
            auto name = path_handle_graph->get_path_name(handle);
            if (!Paths::is_alt(name)) {
                unvisited_contigs.insert(name);
            }
        });
    }

    // Deal with GBWT threads
    if (!gbwt_name.empty()) {
        if (progress) {
            logger.info() << "Loading GBWT index " << gbwt_name << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);
        if (!(gbwt_index->hasMetadata()) || !(gbwt_index->metadata.hasSampleNames()) 
                                              || !(gbwt_index->metadata.hasPathNames())) {
            logger.error() << "GBWT index does not contain sufficient metadata" << std::endl;
        }
        
        // we will add these threads to the graph as named paths and index them for easy look up
        hash_map<gbwt::size_type, size_t> sample_id_to_idx;
        
        if (!sample_names.empty()) {
            // we're consulting the provided sample names to determine which threads to include
            if (progress) {
                logger.info() << "Checking " << sample_names.size() << " samples" << std::endl;
            }
            for (std::string& sample_name : sample_names) {
                gbwt::size_type id = gbwt_index->metadata.sample(sample_name);
                if (id >= gbwt_index->metadata.samples()) {
                    logger.error() << "sample \"" << sample_name 
                                   << "\" not found in the GBWT index" << std::endl;
                }
                auto idx = sample_id_to_idx.size();
                sample_id_to_idx[id] = idx;
            }
        }
        else {
            // we are consulting the transcript expression table to decide which threads to include
            for (const auto& transcript_expression  : transcript_expressions) {
                gbwt::size_type id = gbwt_index->metadata.sample(transcript_expression.first);
                if (id >= gbwt_index->metadata.samples()) {
                    logger.error() << "haplotype-specific transcript \"" << transcript_expression.first
                                   << "\" not found in the GBWT index" << std::endl;
                }
                auto idx = sample_id_to_idx.size();
                sample_id_to_idx[id] = idx;
            }
        }
        
        MutablePathMutableHandleGraph* mutable_graph \
            = dynamic_cast<MutablePathMutableHandleGraph*>(path_handle_graph.get());
        if (mutable_graph == nullptr) {
            if (progress) {
                logger.info() << ": Converting the graph into HashGraph" << std::endl;
            }
            mutable_graph = new bdsg::HashGraph();
            handlealgs::copy_path_handle_graph(path_handle_graph.get(), mutable_graph);
            path_handle_graph.reset(mutable_graph);
        }
        if (progress) {
            logger.info() << "Inserting " << sample_id_to_idx.size()
                          << " GBWT threads into the graph" << std::endl;
        }
        
        for (gbwt::size_type i = 0; i < gbwt_index->metadata.paths(); i++) {
            auto& path = gbwt_index->metadata.path(i);
            auto it = sample_id_to_idx.find(path.sample);
            if (it != sample_id_to_idx.end()) {
                std::string path_name = insert_gbwt_path(*mutable_graph, *gbwt_index, i);
                if (!path_name.empty()) {
                    // path was successfully added
                    if (!sample_names.empty()) {
                        // assign this haplotype a ploidy of 1
                        
                        // We managed to make a path for this thread
                        path_names.push_back(path_name);
                        // It should have ploidy 1
                        path_ploidies.push_back(1.0);
                        
                        if (!unvisited_contigs.empty()) {
                            // Remember that the contig this path is on is visited
                            auto contig_name = gbwt_index->metadata.contig(path.contig);
                            unvisited_contigs.erase(contig_name);
                        }
                    }
                    else {
                        // update the transcript name so we can assign it expression
                        // later down
                        transcript_expressions[it->second].first = path_name;
                    }
                    // Remember we inserted a path
                    inserted_path_names.insert(path_name);
                }
            }
        }
        if (progress) {
            logger.info() << "Inserted " << inserted_path_names.size()
                               << " paths" << std::endl;
        }
    } else {
        // We're not using a separate GBWT, so when asked to simulate from a
        // sample, pull that sample from the base graph.
        if (!sample_names.empty()) {
            // we're consulting the provided sample names to determine which threads to include
            
            // Make a query set of sample names
            std::unordered_set<std::string> sample_name_set(sample_names.begin(), sample_names.end());

            if (sample_name_set.size() != sample_names.size()) {
                // Do a quick check for duplicates since we've bothered to make the set.
                logger.error() << "Of the " << sample_names.size()
                               << " samples, there are only " << sample_name_set.size()
                               << " distinct values. Remove the duplicate entries." << std::endl;
            }

            if (progress) {
                logger.info() << "Finding matching paths for "
                              << sample_name_set.size() << " samples" << std::endl;
            }

            // Also keep a set of sample names actually seen
            std::unordered_set<std::string> seen_sample_names;
            // And a count of sample-based paths we picked
            size_t sample_path_count = 0;
            path_handle_graph->for_each_path_matching(nullptr, &sample_name_set, nullptr, [&](const path_handle_t path) {
                // Remember we saw this sample
                seen_sample_names.insert(path_handle_graph->get_sample_name(path));
                
                if (!unvisited_contigs.empty()) {
                    // Remember we visited this path's contig, if it has a full-contig generic path lying around.
                    unvisited_contigs.erase(path_handle_graph->get_locus_name(path));
                }

                // Remember to simulate from this path by name with ploidy 1.
                path_names.push_back(path_handle_graph->get_path_name(path));
                path_ploidies.push_back(1.0);
                ++sample_path_count;
            });

            if (seen_sample_names.size() != sample_name_set.size()) {
                // TODO: Use std::set_difference in C++17
                auto error_msg = logger.error();
                error_msg << "Some samples requested are not in the graph:";
                for (auto& s : sample_name_set) {
                    if (!seen_sample_names.count(s)) {
                        error_msg << " " << s;
                    }
                }
            }

            if (progress) {
                logger.info() << "Using " << sample_path_count
                              << " sample paths" << std::endl;
            }
        }
    }

    if (!unvisited_contigs.empty()) {
        // There are unvisited contigs we want to sample from too
        for (auto& name : unvisited_contigs) {
            // Sample from each
            path_names.push_back(name);
            // With the rule-determined ploidy
            path_ploidies.push_back(consult_ploidy_rules(name));
        }
        if (progress) {
            logger.info() << "Also sampling from " << unvisited_contigs.size()
                          << " paths representing unvisited contigs" << std::endl;
        }
    }
    
    if (haplotype_transcript_file_name.empty()) {
        if (!transcript_expressions.empty()) {
            if (progress) {
                logger.info() << "Checking " << transcript_expressions.size()
                              << " transcripts" << std::endl;
            }
            for (auto& transcript_expression : transcript_expressions) {
                if (!path_handle_graph->has_path(transcript_expression.first)) {
                    logger.error() << "transcript path \"" << transcript_expression.first 
                                   << "\" not found in index. If you embedded haplotype-specific transcripts "
                                   << "in the graph, you may need the "
                                   << "haplotype transcript file from vg rna -i" << std::endl;
                }
            }
        }
    }
    else {
        if (progress) {
            logger.info() << "Checking " << haplotype_transcripts.size()
                          << " haplotype transcripts" << std::endl;
        }
        for (auto& haplotype_transcript : haplotype_transcripts) {
            if (!path_handle_graph->has_path(get<0>(haplotype_transcript))) {
                logger.error() << "transcript path for \"" << get<0>(haplotype_transcript) 
                               << "\" not found in index." << std::endl;
            }
        }
    }
    
    if (progress) {
        logger.info() << "Creating path position overlay" << std::endl;
    }
    
    bdsg::ReferencePathVectorizableOverlayHelper overlay_helper;
    // Work out the list of path names we definitely need position queries on,
    // even if they aren't reference sense.
    // Just put all the target path names in there, since it can't hurt.
    std::unordered_set<std::string> extra_path_names(path_names.begin(), path_names.end());
    // Apply the overlay to ensure position lookups are efficient.
    PathPositionHandleGraph* xgidx = dynamic_cast<PathPositionHandleGraph*>(
        overlay_helper.apply(
            path_handle_graph.get(),
            extra_path_names
        )
    );
    
    // We want to store the inserted paths as a set of handles, which are
    // easier to hash than strings for lookup.
    unordered_set<path_handle_t> inserted_path_handles;
    if (!inserted_path_names.empty()) {
        if (progress) {
            logger.info() << "Finding inserted paths" << std::endl;
        }
        for (auto& name : inserted_path_names) {
            inserted_path_handles.insert(xgidx->get_path_handle(name));
        }
    }
    
    unique_ptr<AlignmentEmitter> alignment_emitter;
    if (align_out) {
        // We're writing in an alignment format
        alignment_emitter = get_non_hts_alignment_emitter("-", json_out ? "JSON" : "GAM",
                                                          map<string, int64_t>(), get_thread_count());
    }
    // Otherwise we're just dumping sequence strings; leave it null.
    
    if (progress) {
        auto opt_info = logger.info();
        opt_info << "Simulating " << (fragment_length > 0 ? "read pairs" : "reads") << std::endl
                 << "--num-reads " << num_reads << std::endl
                 << "--read-length " << read_length << std::endl;
        if (align_out) {
            opt_info << "--align-out" << std::endl;
        }
        if (json_out) {
            opt_info << "--json-out" << std::endl;
        }
        if (fastq_out) {
            opt_info << "--fastq-out" << std::endl;
        }
        if (!fastq_name.empty()) {
            opt_info << "--fastq " << fastq_name << std::endl;
            if (!fastq_2_name.empty()) {
                opt_info << "--fastq " << fastq_2_name << std::endl;
            }
            if (interleaved) {
                opt_info << "--interleaved" << std::endl;
            }
        } else {
            if (base_error > 0.0) {
                opt_info << "--sub-rate " << base_error << std::endl;
            }
        }
        if (indel_error > 0.0) {
            opt_info << "--indel-rate " << indel_error << std::endl;
        }
        if (!fastq_name.empty()) {
            if (indel_prop > 0.0) {
                opt_info << "--indel-err-prop " << indel_prop << std::endl;
            }
            if (error_scale_factor != 1.0) {
                opt_info << "--scale-err " << error_scale_factor << std::endl;
            }
        }
        if (forward_only) {
            opt_info << "--forward-only" << std::endl;
        }
        if (fragment_length > 0) {
            opt_info << "--frag-len " << fragment_length << std::endl;
            if (fragment_std_dev > 0.0) {
                opt_info << "--frag-std-dev " << fragment_std_dev << std::endl;
            }
        }
        if (reads_may_contain_Ns) {
            opt_info << "--allow-Ns" << std::endl;
        }
        if (max_tries != 100) {
            opt_info << "--max-tries" << max_tries << std::endl;
        }
    }
    
    unique_ptr<AbstractReadSampler> sampler;
    if (fastq_name.empty()) {
        // Use the fixed error rate sampler
        
        if (unsheared_fragments) {
            logger.warn() << "Unsheared fragments option is only available "
                          << "when simulating from FASTQ-trained errors" << std::endl;
        }
        
        sampler.reset(new Sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns, path_names,
                                  path_ploidies, transcript_expressions, haplotype_transcripts));
    } else {
        // Use the FASTQ-trained sampler
        sampler.reset(new NGSSimulator(*xgidx,
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
                                       // suppresses warnings about fragment length
                                       fragment_length ? fragment_length : std::numeric_limits<double>::max(), 
                                       // eliminates errors from having 0 as stddev without substantial difference
                                       fragment_std_dev ? fragment_std_dev : 0.000001,
                                       error_scale_factor,
                                       !reads_may_contain_Ns,
                                       unsheared_fragments,
                                       seed_val));
    }
    
    // Do common configuration
    sampler->multi_position_annotations = multi_position_annotations;
    sampler->max_tries = max_tries;
    if (!inserted_path_handles.empty()) {
        // Skip paths that we have added ourselves when annotating, so we search
        // further for base-graph path.
        std::function<bool(const path_handle_t&)> annotation_path_filter = [&inserted_path_handles](const path_handle_t& path) {
            return !inserted_path_handles.count(path);
        };
        sampler->annotation_path_filter = std::make_unique<std::function<bool(const path_handle_t&)>>(std::move(annotation_path_filter));
    }
    
    // Generate an Aligner for rescoring
    Aligner aligner(default_score_matrix, default_gap_open, default_gap_extension,
                    default_full_length_bonus, vg::default_gc_content);
                    
    // We define a function to score a using the aligner
    auto rescore = [&] (Alignment& aln) {
        // Score using exact distance.
        aln.set_score(aligner.score_contiguous_alignment(aln, !strip_bonuses, !strip_bonuses));
    };
    
    // And a function to emit either single or paired reads, while recomputing scores.
    auto emit = [&] (Alignment* r1, Alignment* r2) {
        // write the alignment or its string
        if (align_out) {
            // write it out as requested
            
            // We will need scores
            rescore(*r1);
            if (r2) {
                // And we have a paired read
                rescore(*r2);
                alignment_emitter->emit_pair(std::move(*r1), std::move(*r2));
            } else {
                // We have just one read.
                alignment_emitter->emit_single(std::move(*r1));
            }
        } else {
            // Print the sequences of the reads we have.
            #pragma omp critical
            {
                if (fastq_out) {
                    cout << "@" << r1->name() << endl;
                }
                cout << r1->sequence();
                if (fastq_out) {
                    cout << endl << "+" << endl;
                    for (size_t i = 0; i < r1->sequence().size(); ++i) {
                        cout << "I"; // assume all bases are perfect
                    }
                }
                if (r2) {
                    if (fastq_out) {
                        cout << endl << "@" << r2->name() << endl
                             << r2->sequence() << endl
                             << "+" << endl;
                        for (size_t i = 0; i < r2->sequence().size(); ++i) {
                            cout << "I"; // assume all bases are perfect
                        }
                    } else {
                        cout << "\t" << r2->sequence();
                    }
                }
                cout << endl;
            }
        }
    };
    
    // The rest of the process has to split up by the type of sampler in use.
    // TODO: Actually refactor to a common sampling interface.

    if (dynamic_cast<Sampler*>(sampler.get())) {
        // Not everything is bound to the new interface yet, so we need to do
        // actual sampling through this typed pointer.
        Sampler* basic_sampler = dynamic_cast<Sampler*>(sampler.get());
        
        size_t max_iter = 1000;
        int nonce = 1;
        for (int i = 0; i < num_reads; ++i) {
            // For each read we are going to generate
            
            if (fragment_length) {
                // fragment_length is nonzero so make it two paired reads
                auto alns = basic_sampler->alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (alns.front().sequence().size() < read_length
                        || alns.back().sequence().size() < read_length) {
                        // If our read was too short, try again
                        alns = basic_sampler->alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                    }
                }
                
                // write the alignment or its string
                emit(&alns.front(), &alns.back()); 
            } else {
                // Do single-end reads
                auto aln = basic_sampler->alignment_with_error(read_length, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (aln.sequence().size() < read_length) {
                        // If our read is too short, try again
                        auto aln_prime = basic_sampler->alignment_with_error(read_length, base_error, indel_error);
                        if (aln_prime.sequence().size() > aln.sequence().size()) {
                            // But only keep the new try if it is longer
                            aln = aln_prime;
                        }
                    }
                }
                
                // Emit the unpaired alignment
                emit(&aln, nullptr);
            }
        }
        
    } else if (dynamic_cast<NGSSimulator*>(sampler.get())) {
        // Use the trained error rate
        
        // Not everything is bound to the new interface yet, so we need to do
        // actual sampling through this typed pointer.
        NGSSimulator* ngs_sampler = dynamic_cast<NGSSimulator*>(sampler.get());
        
        if (!path_pos_filename.empty()) {
            ngs_sampler->connect_to_position_file(path_pos_filename);
        }
        
        // static scheduling could produce some degradation in speed, but I think it should make
        // the output deterministic (except for ordering)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < num_reads; i++) {
            if (fragment_length) {
                pair<Alignment, Alignment> read_pair = ngs_sampler->sample_read_pair();
                emit(&read_pair.first, &read_pair.second);
            }
            else {
                Alignment read = ngs_sampler->sample_read();
                emit(&read, nullptr);
            }
        }
    } else {
        // We don't know about this sampler type.
        // TODO: Define a real sampler interface that lets you sample.
        throw std::logic_error("Attempted to use sampler type for which sampling is not implemented!");
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_sim("sim", "simulate reads from a graph", TOOLKIT, main_sim);

