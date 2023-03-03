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
#include "bdsg/packed_graph.hpp"
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/utils.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_rna(char** argv) {
    cerr << "\nusage: " << argv[0] << " rna [options] graph.[vg|pg|hg|gbz] > splicing_graph.[vg|pg|hg]" << endl

         << "\nGeneral options:" << endl

         << "    -t, --threads INT          number of compute threads to use [1]" << endl
         << "    -p, --progress             show progress" << endl
         << "    -h, --help                 print help message" << endl

         << "\nInput options:" << endl

         << "    -n, --transcripts FILE     transcript file(s) in gtf/gff format; may repeat" << endl
         << "    -m, --introns FILE         intron file(s) in bed format; may repeat" << endl
         << "    -y, --feature-type NAME    parse only this feature type in the gtf/gff (parses all if empty) [exon]" << endl
         << "    -s, --transcript-tag NAME  use this attribute tag in the gtf/gff file(s) as id [transcript_id]" << endl
         << "    -l, --haplotypes FILE      project transcripts onto haplotypes in GBWT index file" << endl
         << "    -z, --gbz-format           input graph is in GBZ format (contains both a graph and haplotypes (GBWT index))" << endl

         << "\nConstruction options:" << endl

         << "    -j, --use-hap-ref          use haplotype paths in GBWT index as reference sequences (disables projection)" << endl
         << "    -e, --proj-embed-paths     project transcripts onto embedded haplotype paths" << endl
         << "    -c, --path-collapse TYPE   collapse identical transcript paths across no|haplotype|all paths [haplotype]" << endl
         << "    -k, --max-node-length INT  chop nodes longer than maximum node length (0 disables chopping) [0]" << endl
         << "    -d, --remove-non-gene      remove intergenic and intronic regions (deletes all paths in the graph)" << endl
         << "    -o, --do-not-sort          do not topological sort and compact the graph" << endl
         << "    -r, --add-ref-paths        add reference transcripts as embedded paths in the graph" << endl
         << "    -a, --add-hap-paths        add projected transcripts as embedded paths in the graph" << endl

         << "\nOutput options:" << endl

         << "    -b, --write-gbwt FILE      write pantranscriptome transcript paths as GBWT index file" << endl
         << "    -f, --write-fasta FILE     write pantranscriptome transcript sequences as fasta file" << endl
         << "    -i, --write-info FILE      write pantranscriptome transcript info table as tsv file" << endl
         << "    -q, --out-exclude-ref      exclude reference transcripts from pantranscriptome output" << endl
         << "    -g, --gbwt-bidirectional   use bidirectional paths in GBWT index construction" << endl

         << endl;
}

int32_t main_rna(int32_t argc, char** argv) {

    if (argc == 2) {
        help_rna(argv);
        return 1;
    }
    
    vector<string> transcript_filenames;
    vector<string> intron_filenames;
    string feature_type = "exon";
    string transcript_tag = "transcript_id";
    string haplotypes_filename;
    bool gbz_format = false;
    bool use_hap_ref = false;
    bool proj_emded_paths = false;
    string path_collapse_type = "haplotype";
    uint32_t max_node_length = 0;
    bool remove_non_transcribed_nodes = false;
    bool sort_collapse_graph = true;
    bool add_reference_transcript_paths = false;
    bool add_projected_transcript_paths = false;
    bool exclude_reference_transcripts = false;
    string gbwt_out_filename = "";
    bool gbwt_add_bidirectional = false;
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
                {"introns",  no_argument, 0, 'm'},
                {"feature-type",  no_argument, 0, 'y'},
                {"transcript-tag",  no_argument, 0, 's'},
                {"haplotypes",  no_argument, 0, 'l'},
                {"gbz-format",  no_argument, 0, 'z'},
                {"use-hap-ref",  no_argument, 0, 'j'},
                {"proj-embed-paths",  no_argument, 0, 'e'},
                {"path-collapse",  no_argument, 0, 'c'},
                {"max-node-length",  no_argument, 0, 'k'},
                {"remove-non-gene",  no_argument, 0, 'd'},
                {"do-not-sort",  no_argument, 0, 'o'},
                {"add-ref-paths",  no_argument, 0, 'r'},
                {"add-hap-paths",  no_argument, 0, 'a'},      
                {"write-gbwt",  no_argument, 0, 'b'},
                {"write-fasta",  no_argument, 0, 'f'},
                {"write-info",  no_argument, 0, 'i'},
                {"out-ref-paths",  no_argument, 0, 'u'},
                {"out-exclude-ref",  no_argument, 0, 'q'},
                {"gbwt-bidirectional",  no_argument, 0, 'g'},   
                {"threads",  no_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "n:m:y:s:l:zjec:k:dorab:f:i:uqgt:ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'n':
            transcript_filenames.push_back(optarg);
            break;

        case 'm':
            intron_filenames.push_back(optarg);
            break;

        case 'y':
            feature_type = optarg;
            break;

        case 's':
            transcript_tag = optarg;
            break;

        case 'l':
            haplotypes_filename = optarg;
            break;

        case 'z':
            gbz_format = true;
            break;

        case 'j':
            use_hap_ref = true;
            break;

        case 'e':
            proj_emded_paths = true;
            break;

        case 'c':
            path_collapse_type = optarg;
            break;

        case 'k':
            max_node_length = stoi(optarg);
            break;

        case 'd':
            remove_non_transcribed_nodes = true;
            break;

        case 'o':
            sort_collapse_graph = false;
            break;

        case 'r':
            add_reference_transcript_paths = true;
            break;

        case 'a':
            add_projected_transcript_paths = true;
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

        case 'u':
            exclude_reference_transcripts = false;
            break;

        case 'q':
            exclude_reference_transcripts = true;
            break;

        case 'g':
            gbwt_add_bidirectional = true;
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

    if (transcript_filenames.empty() && intron_filenames.empty()) {

        cerr << "[vg rna] ERROR: No transcripts or introns were given. Use --transcripts FILE and/or --introns FILE." << endl;
        return 1;       
    }

    if (!haplotypes_filename.empty() && gbz_format) {

        cerr << "[vg rna] ERROR: Only one set of haplotypes can be provided (GBZ file contains both a graph and haplotypes). Use either --haplotypes or --gbz-format." << endl;
        return 1;       
    }

    if (remove_non_transcribed_nodes && !add_reference_transcript_paths && !add_projected_transcript_paths) {

        cerr << "[vg rna] WARNING: Reference paths are deleted when removing intergenic and intronic regions. Consider adding transcripts as embedded paths using --add-ref-paths and/or --add-hap-paths." << endl;
    }

    if (path_collapse_type != "no" && path_collapse_type != "haplotype" && path_collapse_type != "all") {

        cerr << "[vg rna] ERROR: Path collapse type (--path-collapse) provided not supported. Options: no, haplotype or all." << endl;
        return 1;
    }

    double time_parsing_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Parsing graph file ..." << endl; }

    string graph_filename = get_input_file_name(optind, argc, argv);

    unique_ptr<MutablePathDeletableHandleGraph> graph(nullptr);
    unique_ptr<gbwt::GBWT> haplotype_index;

    if (!gbz_format) {

        // Load pangenome graph.
        graph = move(vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(graph_filename));
    
        if (!haplotypes_filename.empty()) {

            // Load haplotype GBWT index.
            if (show_progress) { cerr << "[vg rna] Parsing haplotype GBWT index file ..." << endl; }
            haplotype_index = vg::io::VPKG::load_one<gbwt::GBWT>(haplotypes_filename);
            assert(haplotype_index->bidirectional());

        } else {

            // Construct empty GBWT index if non is given. 
            haplotype_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
        }

    } else {

        graph = unique_ptr<MutablePathDeletableHandleGraph>(new bdsg::PackedGraph());

        // Load GBZ file 
        unique_ptr<gbwtgraph::GBZ> gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(graph_filename);

        if (show_progress) { cerr << "[vg rna] Converting graph format ..." << endl; }

        // Convert GBWTGraph to mutable graph type (PackedGraph).
        graph->set_id_increment(gbz->graph.min_node_id());
        handlealgs::copy_handle_graph(&(gbz->graph), graph.get());

        // Copy reference and generic paths to new graph.
        gbz->graph.for_each_path_matching({PathSense::GENERIC, PathSense::REFERENCE}, {}, {}, [&](const path_handle_t& path) {
            
            handlegraph::algorithms::copy_path(&(gbz->graph), path, graph.get());
        });

        haplotype_index = make_unique<gbwt::GBWT>(gbz->index);
    }

    if (graph == nullptr) {
        cerr << "[transcriptome] ERROR: Could not load graph." << endl;
        exit(1);
    }

    // Construct transcriptome and parse graph.
    Transcriptome transcriptome(move(graph));
    assert(graph == nullptr);

    transcriptome.show_progress = show_progress;
    transcriptome.num_threads = num_threads;
    transcriptome.feature_type = feature_type;
    transcriptome.transcript_tag = transcript_tag;
    transcriptome.path_collapse_type = path_collapse_type;

    if (show_progress) { cerr << "[vg rna] Graph " << ((!haplotype_index->empty()) ? "and GBWT index " : "") << "parsed in " << gcsa::readTimer() - time_parsing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };


    if (!intron_filenames.empty()) {

        double time_intron_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Adding intron splice-junctions to graph ..." << endl; }

        vector<istream *> intron_streams;
        intron_streams.reserve(intron_filenames.size());

        for (auto & filename: intron_filenames) {

            auto intron_stream = new ifstream(filename);
            intron_streams.emplace_back(intron_stream);
        }

        // Add introns as novel splice-junctions to graph.
        transcriptome.add_intron_splice_junctions(intron_streams, haplotype_index, true);

        for (auto & intron_stream: intron_streams) {

            delete intron_stream;
        }

        if (show_progress) { cerr << "[vg rna] Introns parsed and graph updated in " << gcsa::readTimer() - time_intron_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    vector<istream *> transcript_streams;

    if (!transcript_filenames.empty()) {

        double time_transcript_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Adding transcript splice-junctions and exon boundaries to graph ..." << endl; }

        transcript_streams.reserve(transcript_filenames.size());

        for (auto & filename: transcript_filenames) {

            auto transcript_stream = new ifstream(filename);
            transcript_streams.emplace_back(transcript_stream);
        }

        // Add transcripts as novel exon boundaries and splice-junctions to graph.
        transcriptome.add_reference_transcripts(transcript_streams, haplotype_index, use_hap_ref, !use_hap_ref);

        if (show_progress) { cerr << "[vg rna] Transcripts parsed and graph updated in " << gcsa::readTimer() - time_transcript_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    if (!transcript_streams.empty() && (!haplotype_index->empty() || proj_emded_paths) && !use_hap_ref) {

        double time_project_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Projecting transcripts to haplotypes ..." << endl; }

        for (auto & transcript_stream: transcript_streams) {

            // Reset transcript file streams.
            transcript_stream->clear();
            transcript_stream->seekg(0);
        }

        // Add transcripts to transcriptome by projecting them onto embedded paths 
        // in a graph and/or haplotypes in a GBWT index.
        transcriptome.add_haplotype_transcripts(transcript_streams, *haplotype_index, proj_emded_paths);

        if (show_progress) { cerr << "[vg rna] Haplotype-specific transcripts constructed in " << gcsa::readTimer() - time_project_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }

    for (auto & transcript_stream: transcript_streams) {

        delete transcript_stream;
    }


    if (remove_non_transcribed_nodes) {

        double time_remove_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Removing non-transcribed regions ..." << endl; }

        transcriptome.remove_non_transcribed_nodes();

        if (show_progress) { cerr << "[vg rna] Regions removed in " << gcsa::readTimer() - time_remove_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    if (max_node_length > 0) {

        double time_chop_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Chopping long nodes ..." << endl; }

        transcriptome.chop_nodes(max_node_length);

        if (show_progress) { cerr << "[vg rna] Nodes chopped in " << gcsa::readTimer() - time_chop_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    if (sort_collapse_graph) {
    
        double time_sort_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Topological sorting graph and compacting node ids ..." << endl; }
        
        if (transcriptome.sort_compact_nodes()) {

            if (show_progress) { cerr << "[vg rna] Graph sorted and compacted in " << gcsa::readTimer() - time_sort_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

        } else {

            if (show_progress) { cerr << "[vg rna] WARNING: Can only sort and compact node ids for a graph in the PackedGraph format" << endl; };            
        }        
    }


    if (add_reference_transcript_paths || add_projected_transcript_paths) {

        double time_add_start = gcsa::readTimer();

        if (add_reference_transcript_paths && add_projected_transcript_paths) {

            if (show_progress) { cerr << "[vg rna] Adding reference and projected transcripts as embedded paths in the graph ..." << endl; }

        } else {

            if (show_progress) { cerr << "[vg rna] Adding " << ((add_reference_transcript_paths) ? "reference" : "projected") << " transcripts as embedded paths in the graph ..." << endl; }
        }

        transcriptome.embed_transcript_paths(add_reference_transcript_paths, add_projected_transcript_paths);

        if (show_progress) { cerr << "[vg rna] Transcript paths added in " << gcsa::readTimer() - time_add_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    double time_writing_start = gcsa::readTimer();

    bool write_pantranscriptome = (!gbwt_out_filename.empty() || !fasta_out_filename.empty() || !info_out_filename.empty());

    if (write_pantranscriptome) {

        if (show_progress) { cerr << "[vg rna] Writing pantranscriptome transcripts to file(s) ..." << endl; }
    }

    // Write transcript paths in transcriptome as GBWT index.
    if (!gbwt_out_filename.empty()) {

        // Silence GBWT index construction. 
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
        gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(transcriptome.graph().max_node_id(), true)), gbwt::DynamicGBWT::INSERT_BATCH_SIZE, gbwt::DynamicGBWT::SAMPLE_INTERVAL);

        transcriptome.add_transcripts_to_gbwt(&gbwt_builder, gbwt_add_bidirectional, exclude_reference_transcripts);

        assert(gbwt_builder.index.hasMetadata());

        // Finish contruction and recode index.
        gbwt_builder.finish();
        save_gbwt(gbwt_builder.index, gbwt_out_filename);
    }

    // Write transcript sequences in transcriptome as fasta file.
    if (!fasta_out_filename.empty()) {

        ofstream fasta_ostream;
        fasta_ostream.open(fasta_out_filename);

        transcriptome.write_transcript_sequences(&fasta_ostream, exclude_reference_transcripts);
     
        fasta_ostream.close();
    }    

    // Write transcript info in transcriptome as tsv file.
    if (!info_out_filename.empty()) {

        ofstream info_ostream;
        info_ostream.open(info_out_filename);

        transcriptome.write_transcript_info(&info_ostream, *haplotype_index, exclude_reference_transcripts);

        info_ostream.close();
    }    

    if (show_progress) { cerr << "[vg rna] Writing splicing graph to stdout ..." << endl; }

    // Write splicing graph to stdout 
    transcriptome.write_graph(&cout);

    if (show_progress) { cerr << "[vg rna] Graph " << (write_pantranscriptome ? "and pantranscriptome " : "") << "written in " << gcsa::readTimer() - time_writing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    return 0;
}

// Register subcommand
static Subcommand vg_rna("rna", "construct splicing graphs and pantranscriptomes", PIPELINE, 3, main_rna);

