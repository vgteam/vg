/** \file rna_main.cpp
 *
 * Defines the "vg rna" subcommand.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>

#include "subcommand.hpp"

#include "../transcriptome.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../gbwt_helper.hpp"

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
         << "    -c, --do-not-collapse      do not collapse identical transcripts across haplotypes" << endl
         << "    -d, --remove-non-gene      remove intergenic and intronic regions (removes reference paths if -a or -r)" << endl
         << "    -r, --add-ref-paths        add reference transcripts as embedded paths in the graph" << endl
         << "    -a, --add-non-ref-paths    add non-reference transcripts as embedded paths in the graph" << endl
         << "    -u, --out-ref-paths        output reference transcripts in GBWT, fasta and info" << endl
         << "    -b, --write-gbwt FILE      write transcripts as threads to GBWT index file" << endl
         << "    -f, --write-fasta FILE     write transcripts as sequences to fasta file" << endl
         << "    -i, --write-info FILE      write transcript origin info to tsv file" << endl
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
    bool collapse_transcript_paths = true;
    bool remove_non_transcribed = false;
    bool add_reference_transcript_paths = false;
    bool add_non_reference_transcript_paths = false;
    bool output_reference_transcript_paths = false;
    string gbwt_out_filename = "";
    string fasta_out_filename = "";
    string info_out_filename = "";
    int32_t num_threads = 1;
    bool show_progress = false;

    int32_t c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
            {
                {"transcripts",  no_argument, 0, 'n'},
                {"transcript-tag",  no_argument, 0, 's'},
                {"haplotypes",  no_argument, 0, 'l'},
                {"use-embeded-paths",  no_argument, 0, 'e'},
                {"do-not-collapse",  no_argument, 0, 'c'},
                {"remove-non-gene",  no_argument, 0, 'd'},
                {"add-ref-paths",  no_argument, 0, 'r'},
                {"add-non-ref-paths",  no_argument, 0, 'a'},
                {"out-ref-paths",  no_argument, 0, 'u'},           
                {"write-gbwt",  no_argument, 0, 'b'},
                {"write-fasta",  no_argument, 0, 'f'},
                {"write-info",  no_argument, 0, 'i'},
                {"threads",  no_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "n:s:l:ercdraub:f:i:t:ph?", long_options, &option_index);

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

        case 'c':
            collapse_transcript_paths = false;
            break;

        case 'd':
            remove_non_transcribed = true;
            break;

        case 'r':
            add_reference_transcript_paths = true;
            break;

        case 'a':
            add_non_reference_transcript_paths = true;
            break;

        case 'u':
            output_reference_transcript_paths = true;
            break;

        case 'b':
            gbwt_out_filename = optarg;
            break;

        case 'f':
            fasta_out_filename = optarg;
            break;

        case 'i':
            info_out_filename = optarg;
            break;

        case 't':
            num_threads = stoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'h':
        case '?':
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


    double time_parsing_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Parsing graph file ..." << endl; }

    // Construct transcriptome and parse variation graph.
    Transcriptome transcriptome(get_input_file_name(optind, argc, argv), show_progress);

    unique_ptr<gbwt::GBWT> haplotype_index;

    if (!haplotypes_filename.empty()) {

        // Load haplotype GBWT index.
        if (show_progress) { cerr << "[vg rna] Parsing haplotype GBWT index file ..." << endl; }
        haplotype_index = vg::io::VPKG::load_one<gbwt::GBWT>(haplotypes_filename);
        assert(haplotype_index->bidirectional());

    } else {

        // Construct empty GBWT index if no is given. 
        haplotype_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
    }

    if (show_progress) { cerr << "[vg rna] Graph (and index) parsed in " << gcsa::readTimer() - time_parsing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    transcriptome.num_threads = num_threads;
    transcriptome.transcript_tag = transcript_tag;
    transcriptome.use_embedded_paths = use_embedded_paths;
    transcriptome.use_reference_paths = (add_reference_transcript_paths || output_reference_transcript_paths);
    transcriptome.collapse_transcript_paths = collapse_transcript_paths;


    double time_project_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Parsing and projecting transcripts ..." << endl; }

    // Add transcripts to transcriptome by projecting them onto embedded paths 
    // in a graph and/or haplotypes in a GBWT index. Edit graph with 
    // transcriptome splice-junctions.
    for (auto & filename: transcript_filenames) {

        get_input_file(filename, [&](istream& transcript_stream) {
                   transcriptome.add_transcripts(transcript_stream, *haplotype_index);
        });
    }

    if (show_progress) { cerr << "[vg rna] Transcripts parsed and projected in " << gcsa::readTimer() - time_project_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    // Release and delete GBWT index pointer.
    haplotype_index.reset(nullptr);


    if (remove_non_transcribed) {

        double time_remove_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Removing non-transcribed regions ..." << endl; }

        transcriptome.remove_non_transcribed(!(add_reference_transcript_paths || add_non_reference_transcript_paths));

        if (show_progress) { cerr << "[vg rna] Regions removed in " << gcsa::readTimer() - time_remove_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    double time_sort_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Topological sorting and compacting graph ..." << endl; }
    
    transcriptome.compact_ordered();
    
    if (show_progress) { cerr << "[vg rna] Graph sorted and compacted in " << gcsa::readTimer() - time_sort_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };


    if (add_reference_transcript_paths || add_non_reference_transcript_paths) {

        double time_add_start = gcsa::readTimer();

        if (add_reference_transcript_paths && add_non_reference_transcript_paths) {

            if (show_progress) { cerr << "[vg rna] Adding all transcript paths to graph ..." << endl; }

        } else {

            if (show_progress) { cerr << "[vg rna] Adding " << ((add_reference_transcript_paths) ? "reference" : "non-reference") << " transcript paths to graph ..." << endl; }
        }

        transcriptome.add_paths_to_graph(add_reference_transcript_paths, add_non_reference_transcript_paths, false);

        if (show_progress) { cerr << "[vg rna] Paths added in " << gcsa::readTimer() - time_add_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    double time_writing_start = gcsa::readTimer();

    // Construct and write GBWT index of transcript paths in transcriptome.
    if (!gbwt_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing " << transcriptome.size() << " transcripts as threads to GBWT index file ..." << endl; }

        // Silence GBWT index construction. 
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
        gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(transcriptome.splice_graph().max_node_id(), true)));

        transcriptome.construct_gbwt(&gbwt_builder, output_reference_transcript_paths);

        // Finish contruction and recode index.
        gbwt_builder.finish();

        vg::io::VPKG::save(gbwt_builder.index, gbwt_out_filename);
    }

    // Write transcript path sequences in transcriptome to fasta file.
    if (!fasta_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing " << transcriptome.size() << " transcripts as sequences to fasta file ..." << endl; }

        ofstream fasta_ostream;
        fasta_ostream.open(fasta_out_filename);
        transcriptome.write_sequences(&fasta_ostream, output_reference_transcript_paths);
        fasta_ostream.close();
    }    

    // Write origin info on transcripts in transcriptome to tsv file.
    if (!info_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing origin info on " << transcriptome.size() << " transcripts to tsv file ..." << endl; }

        ofstream info_ostream;
        info_ostream.open(info_out_filename);
        transcriptome.write_info(&info_ostream, output_reference_transcript_paths);
        info_ostream.close();
    }    

    if (show_progress) { cerr << "[vg rna] Writing graph to stdout ..." << endl; }

    // Write spliced variation graph to stdout 
    transcriptome.write_graph(&cout);

    if (show_progress) { cerr << "[vg rna] Graph (and transcripts) written in " << gcsa::readTimer() - time_writing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    return 0;
}

// Register subcommand
static Subcommand vg_rna("rna", "construct spliced variation graphs and transcript paths", main_rna);

