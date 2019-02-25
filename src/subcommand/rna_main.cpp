/** \file rna_main.cpp
 *
 * Defines the "vg rna" subcommand, which projects rna transcripts to graph paths and/or GBWT haplotypes.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>

#include "subcommand.hpp"

#include "../transcriptome.hpp"
#include "../stream/vpkg.hpp"
#include "../stream/stream.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_rna(char** argv) {
    cerr << "\nusage: " << argv[0] << " rna [options] <graph.vg> > splice_graph.vg" << endl
         << "options:" << endl
         << "    -n, --transcripts FILE     transcript file(s) in gtf/gff format; may repeat (required)" << endl
         << "    -s, --transcript-tag NAME  use this attribute tag in the gtf/gff file(s) as id [transcript_id]" << endl
         << "    -l, --haplotypes FILE      project transcripts onto haplotypes in GBWT index file" << endl
         << "    -e, --use-embedded-paths   project transcripts onto embedded graph paths" << endl
         << "    -r, --filter-reference     filter reference transcripts" << endl
         << "    -c, --do-not-collapse      do not collapse identical transcripts across haplotypes" << endl
         << "    -a, --add-paths            add transcripts as embedded paths in the graph" << endl
         << "    -b, --write-gbwt FILE      write transcripts as path to GBWT index file" << endl
         << "    -g, --write-gam FILE       write transcripts as alignments to GAM file" << endl
         << "    -f, --write-fasta FILE     write transcripts as sequences to fasta file" << endl
         << "    -t, --threads INT          number of compute threads to use [1]" << endl
         << "    -p, --progress             show progress" << endl
         << "    -h, --help                 print help message" << endl
         << endl;
}

int32_t main_rna(int32_t argc, char** argv) {

    if (argc == 2) {
        help_rna(argv);
        return 1;
    }
    
    vector<string> transcript_filenames;
    string transcript_tag = "transcript_id";
    string haplotypes_filename;
    bool use_embedded_paths = false;
    bool filter_reference_transcript_paths = false;
    bool collapse_transcript_paths = true;
    bool add_transcript_paths = false;
    string gbwt_out_filename = "";
    string gam_out_filename = "";
    string fasta_out_filename = "";
    int32_t num_threads = 1;
    bool show_progress = false;

    int32_t c;
    optind = 2; // force optind past command positional argument

    while (true) {
        static struct option long_options[] =
            {
                {"transcripts",  no_argument, 0, 'n'},
                {"transcript-tag",  no_argument, 0, 's'},
                {"haplotypes",  no_argument, 0, 'l'},
                {"use-embeded-paths",  no_argument, 0, 'e'},
                {"filter-reference",  no_argument, 0, 'r'},
                {"do-not-collapse",  no_argument, 0, 'c'},
                {"add-paths",  no_argument, 0, 'a'},
                {"write-gbwt",  no_argument, 0, 'b'},
                {"write-gam",  no_argument, 0, 'g'},
                {"write-fasta",  no_argument, 0, 'f'},
                {"threads",  no_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "n:s:l:ercab:g:f:t:ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'n':
            transcript_filenames.push_back(optarg);
            break;

        case 's':
            transcript_tag = optarg;
            break;

        case 'l':
            haplotypes_filename = optarg;
            break;

        case 'e':
            use_embedded_paths = true;
            break;

        case 'r':
            filter_reference_transcript_paths = true;
            break;
            
        case 'c':
            collapse_transcript_paths = false;
            break;

        case 'a':
            add_transcript_paths = true;
            break;

        case 'b':
            gbwt_out_filename = optarg;
            break;

        case 'g':
            gam_out_filename = optarg;
            break;

        case 'f':
            fasta_out_filename = optarg;
            break;

        case 't':
            num_threads = stoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_rna(argv);
            exit(1);
            break;

        default:
            abort();
        }
    }

    if (argc < optind + 1) {
        help_rna(argv);
        return 1;
    }

    if (transcript_filenames.empty()) {

        cerr << "[vg rna] ERROR: No transcripts were given. Use --transcripts FILE." << endl;
        return 1;       
    }

    if (haplotypes_filename.empty() && !use_embedded_paths) {

        cerr << "[vg rna] ERROR: No haplotypes or paths were given for transcript projection. Use --haplotypes FILE and/or --use-embeded-paths." << endl;
        return 1;       
    }

    if (show_progress) { cerr << "[vg rna] Parsing graph file ..." << endl; }
    
    // Load the graph
    VG* graph = nullptr;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in, show_progress);
    });

    if (!graph) {
        cerr << "[vg rna] ERROR: Could not load graph." << endl;
        return 1;
    }

    gbwt::GBWT haplotype_index;

    if (!haplotypes_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Parsing haplotype GBWT index file ..." << endl; }

        sdsl::load_from_file(haplotype_index, haplotypes_filename);
    }

    Transcriptome transcriptome;

    transcriptome.num_threads = num_threads;
    transcriptome.transcript_tag = transcript_tag;
    transcriptome.use_embedded_paths = use_embedded_paths;
    transcriptome.collapse_transcript_paths = collapse_transcript_paths;
    transcriptome.filter_reference_transcript_paths = filter_reference_transcript_paths;

    auto time1 = std::chrono::system_clock::now();

    if (show_progress) { cerr << "[vg rna] Parsing and projecting transcripts ..." << endl; }

    for (auto & filename: transcript_filenames) {

        get_input_file(filename, [&](istream& transcript_stream) {
                   transcriptome.add_transcripts(transcript_stream, *graph, haplotype_index);
        });
    }

    auto time2 = std::chrono::system_clock::now();
    cerr << "Time (s): " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << endl;

    if (show_progress) { cerr << "[vg rna] Adding splice-junctions " << (add_transcript_paths ? "and transcript paths " : "") << "to graph ..." << endl; }

    transcriptome.edit_graph(graph, add_transcript_paths);

    auto time3 = std::chrono::system_clock::now();
    cerr << "Time (s): " << std::chrono::duration_cast<std::chrono::seconds>(time3 - time2).count() << endl;

    if (!gbwt_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing " << transcriptome.size() << " transcripts as paths to GBWT index file ..." << endl; }

        gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
        gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(graph->max_node_id(), true)));

        transcriptome.construct_gbwt(&gbwt_builder);
        gbwt_builder.finish();

        gbwt_builder.index.addMetadata();
        gbwt_builder.index.metadata.setHaplotypes(transcriptome.size());
        // gbwt_builder.index.metadata.setSamples();
        // gbwt_builder.index.metadata.setContigs();

        sdsl::store_to_file(gbwt_builder.index, gbwt_out_filename);
    }    

    auto time4 = std::chrono::system_clock::now();
    cerr << "Time (s): " << std::chrono::duration_cast<std::chrono::seconds>(time4 - time3).count() << endl;

    if (!gam_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing " << transcriptome.size() << " transcripts as alignments to GAM file ..." << endl; }

        ofstream gam_ostream;
        gam_ostream.open(gam_out_filename);
        transcriptome.write_gam_alignments(&gam_ostream);
        gam_ostream.close();
    }

    auto time5 = std::chrono::system_clock::now();
    cerr << "Time (s): " << std::chrono::duration_cast<std::chrono::seconds>(time5 - time4).count() << endl;

    if (!fasta_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing " << transcriptome.size() << " transcripts as sequences to fasta file ..." << endl; }

        ofstream fasta_ostream;
        fasta_ostream.open(fasta_out_filename);
        transcriptome.write_fasta_sequences(&fasta_ostream, *graph);
        fasta_ostream.close();
    }    

    auto time6 = std::chrono::system_clock::now();
    cerr << "Time (s): " << std::chrono::duration_cast<std::chrono::seconds>(time6 - time5).count() << endl;

    if (show_progress) { cerr << "[vg rna] Writing graph to stdout ..." << endl; }

    graph->serialize_to_ostream(std::cout);
    delete graph;    

    return 0;
}

// Register subcommand
static Subcommand vg_rna("rna", "project rna transcripts to graph paths and/or GBWT haplotypes", main_rna);

