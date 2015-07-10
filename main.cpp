#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include "gcsa.h"
#include "json2pb.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "vg_set.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "stream.hpp"
#include "alignment.hpp"
#include "convert.hpp"
#include "google/protobuf/stubs/common.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
         << "Transforms alignments to be relative to particular paths." << endl
         << endl
         << "options:" << endl
         << "    -d, --db-name DIR       use the graph in this database" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -p, --into-path NAME    surject into just this path" << endl
         << "    -i, --into-paths FILE   surject into path names listed in FILE (one per line)" << endl
         << "    -P, --into-prefix NAME  surject into all paths with NAME as their prefix" << endl
        //<< "    -H, --header-from FILE  use the header in the SAM/CRAM/BAM file for the output" << endl
        // todo, reenable
        // << "    -c, --cram-output       write CRAM to stdout (default is vg::Aligment/GAM format)" << endl
        // << "    -f, --reference FILE    use this file when writing CRAM to rebuild sequence header" << endl
         << "    -b, --bam-output        write BAM to stdout" << endl
         << "    -s, --sam-output        write SAM to stdout" << endl
         << "    -C, --compression N     level for compression [0-9]" << endl
         << "    -w, --window N          use N nodes on either side of the alignment to surject (default 5)" << endl;
}

int main_surject(int argc, char** argv) {

    if (argc == 2) {
        help_surject(argv);
        return 1;
    }

    string db_name;
    string path_name;
    string path_prefix;
    string path_file;
    string output_type = "gam";
    string input_type = "gam";
    string header_file;
    int compress_level = 9;
    int window = 5;
    string fasta_filename;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"db-name", required_argument, 0, 'd'},
                {"threads", required_argument, 0, 't'},
                {"into-path", required_argument, 0, 'p'},
                {"into-paths", required_argument, 0, 'i'},
                {"into-prefix", required_argument, 0, 'P'},
                {"cram-output", no_argument, 0, 'c'},
                {"reference", required_argument, 0, 'f'},
                {"bam-output", no_argument, 0, 'b'},
                {"sam-output", no_argument, 0, 's'},
                {"header-from", required_argument, 0, 'H'},
                {"compress", required_argument, 0, 'C'},
                {"window", required_argument, 0, 'w'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hd:p:i:P:cbsH:C:t:w:f:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'd':
            db_name = optarg;
            break;

        case 'p':
            path_name = optarg;
            break;

        case 'i':
            path_file = optarg;
            break;

        case 'P':
            path_prefix = optarg;
            break;

        case 'H':
            header_file = optarg;
            break;

        case 'c':
            output_type = "cram";
            break;

        case 'f':
            fasta_filename = optarg;
            break;

        case 'b':
            output_type = "bam";
            break;

        case 's':
            compress_level = -1;
            output_type = "sam";
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'C':
            compress_level = atoi(optarg);
            break;

        case 'w':
            window = atoi(optarg);
            break;

        case 'h':
        case '?':
            help_surject(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    string file_name = argv[optind];

    Index index;
    // open index
    index.open_read_only(db_name);

    set<string> path_names;
    if (!path_file.empty()){
        // open the file
        ifstream in(path_file);
        string line;
        while (std::getline(in,line)) {
            path_names.insert(line);
        }
    } else {
        path_names.insert(path_name);
    }

    if (input_type == "gam") {
        if (output_type == "gam") {
            int thread_count = get_thread_count();
            vector<vector<Alignment> > buffer;
            buffer.resize(thread_count);
            function<void(Alignment&)> lambda = [&index, &path_names, &buffer, &window](Alignment& src) {
                int tid = omp_get_thread_num();
                Alignment surj;
                string path_name;
                int64_t path_pos;
                index.surject_alignment(src, path_names, surj, path_name, path_pos, window);
                buffer[tid].push_back(surj);
                stream::write_buffered(cout, buffer[tid], 1000);
            };
            if (file_name == "-") {
                stream::for_each_parallel(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each_parallel(in, lambda);
            }
            for (int i = 0; i < thread_count; ++i) {
                stream::write_buffered(cout, buffer[i], 0); // flush
            }
        } else {
            char out_mode[5];
            string out_format = "";
            strcpy(out_mode, "w");
            if (output_type == "bam") { out_format = "b"; }
            else if (output_type == "cram") { out_format = "c"; }
            else { out_format = ""; }
            strcat(out_mode, out_format.c_str());
            if (compress_level >= 0) {
                char tmp[2];
                tmp[0] = compress_level + '0'; tmp[1] = '\0';
                strcat(out_mode, tmp);
            }
            // get the header
            /*
            if (header_file.empty()) {
                cerr << "[vg surject] error: --header-from must be specified for SAM/BAM/CRAM output" << endl;
                return 1;
            }
            */
            string header;
            map<string, int64_t> path_by_id = index.paths_by_id();
            map<string, pair<int64_t, int64_t> > path_layout;
            map<string, int64_t> path_length;
            index.path_layout(path_layout, path_length);
            int thread_count = get_thread_count();
            vector<vector<tuple<string, int64_t, Alignment> > > buffer;
            buffer.resize(thread_count);
            map<string, string> rg_sample;

            // bam/sam/cram output
            samFile* out = 0;
            int buffer_limit = 100;

            bam_hdr_t* hdr = NULL;
            int64_t count = 0;
            omp_lock_t output_lock;
            omp_init_lock(&output_lock);

            // handles buffers, possibly opening the output file if we're on the first record
            auto handle_buffer =
                [&hdr, &header, &path_length, &rg_sample, &buffer_limit,
                 &out_mode, &out, &output_lock, &fasta_filename](vector<tuple<string, int64_t, Alignment> >& buf) {
                if (buf.size() >= buffer_limit) {
                    // do we have enough data to open the file?
#pragma omp critical (hts_header)
                    {
                        if (!hdr) {
                            hdr = hts_string_header(header, path_length, rg_sample);
                            if ((out = sam_open("-", out_mode)) == 0) {
                                /*
                                if (!fasta_filename.empty()) {
                                    string fai_filename = fasta_filename + ".fai";
                                    hts_set_fai_filename(out, fai_filename.c_str());
                                }
                                */
                                cerr << "[vg surject] failed to open stdout for writing HTS output" << endl;
                                exit(1);
                            } else {
                                // write the header
                                if (sam_hdr_write(out, hdr) != 0) {
                                    cerr << "[vg surject] error: failed to write the SAM header" << endl;
                                }
                            }
                        }
                    }
                    // try to get a lock, and force things if we've built up a huge buffer waiting
                    if (omp_test_lock(&output_lock) || buf.size() > 10*buffer_limit) {
                        for (auto& s : buf) {
                            auto& path_nom = get<0>(s);
                            auto& path_pos = get<1>(s);
                            auto& surj = get<2>(s);
                            string cigar = cigar_against_path(surj);
                            bam1_t* b = alignment_to_bam(header,
                                                         surj,
                                                         path_nom,
                                                         path_pos,
                                                         cigar,
                                                         "=",
                                                         path_pos,
                                                         0);
                            int r = 0;
#pragma omp critical (cout)
                            r = sam_write1(out, hdr, b);
                            if (r == 0) { cerr << "[vg surject] error: writing to stdout failed" << endl; exit(1); }
                            bam_destroy1(b);
                        }
                        omp_unset_lock(&output_lock);
                        buf.clear();
                    }
                }
            };

            function<void(Alignment&)> lambda = [&index,
                                                 &path_names,
                                                 &path_length,
                                                 &window,
                                                 &rg_sample,
                                                 &header,
                                                 &out,
                                                 &buffer,
                                                 &count,
                                                 &hdr,
                                                 &out_mode,
                                                 &handle_buffer](Alignment& src) {
                int tid = omp_get_thread_num();
                Alignment surj;
                string path_name;
                int64_t path_pos;
                index.surject_alignment(src, path_names, surj, path_name, path_pos, window);
                if (!surj.path().mapping_size()) {
                    surj = src;
                }
                // record 
                if (!hdr && !surj.read_group().empty() && !surj.sample_name().empty()) {
#pragma omp critical (hts_header)
                    rg_sample[surj.read_group()] = surj.sample_name();
                }

                buffer[tid].push_back(make_tuple(path_name, path_pos, surj));
                handle_buffer(buffer[tid]);

            };

            // now apply the alignment processor to the stream
            if (file_name == "-") {
                stream::for_each_parallel(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each_parallel(in, lambda);
            }
            buffer_limit = 0;
            for (auto& buf : buffer) {
                handle_buffer(buf);
            }
            bam_hdr_destroy(hdr);
            sam_close(out);
            omp_destroy_lock(&output_lock);
        }
    }
    cout.flush();

    return 0;
}

void help_mod(char** argv) {
    cerr << "usage: " << argv[0] << " mod [options] <graph.vg> >[mod.vg]" << endl
         << "Modifies graph, outputs modified on stdout." << endl
         << endl
         << "options:" << endl
         << "    -i, --include-aln FILE  include the paths implied by alignments in the graph" << endl
         << "    -c, --compact-ids       should we sort and compact the id space? (default false)" << endl
         << "    -k, --keep-path NAME    keep only nodes and edges in the path" << endl
         << "    -o, --remove-orphans    remove orphan edges from graph (edge specified but node missing)" << endl
         << "    -p, --prune-complex     remove nodes that are reached by paths of --path-length which" << endl
         << "                            cross more than --edge-max edges" << endl
         << "    -p, --prune-complex     remove nodes that are reached by paths of --length which" << endl
         << "                            cross more than --edge-max edges" << endl
         << "    -S, --prune-subgraphs   remove subgraphs which are shorter than --length" << endl
         << "    -l, --length N          for pruning complex regions and short subgraphs" << endl
         << "    -e, --edge-max N        when pruning complex regions only consider paths which cross" << endl
         << "                            this many potential alternate edges (e.g. if node out-degree is" << endl
         << "                            2, we would count 1 toward --edge-max; for 3 we would count 2)" << endl
         << "    -m, --markers           join all head and tails nodes to marker nodes" << endl
         << "                            ('###' starts and '$$$' ends) of --path-length, for debugging" << endl
         << "    -t, --threads N         for tasks that can be done in parallel, use this many threads" << endl;
}

int main_mod(int argc, char** argv) {

    if (argc == 2) {
        help_mod(argv);
        return 1;
    }

    string path_name;
    bool remove_orphans = false;
    string aln_file;
    bool compact_ids = false;
    bool prune_complex = false;
    int path_length = 0;
    int edge_max = 0;
    bool add_start_and_end_markers = false;
    bool prune_subgraphs = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"include-aln", required_argument, 0, 'i'},
                {"compact-ids", no_argument, 0, 'c'},
                {"keep-path", required_argument, 0, 'k'},
                {"remove-orphans", no_argument, 0, 'o'},
                {"prune-complex", no_argument, 0, 'p'},
                {"prune-subgraphs", no_argument, 0, 'S'},
                {"length", required_argument, 0, 'l'},
                {"edge-max", required_argument, 0, 'e'},
                {"markers", no_argument, 0, 'm'},
                {"threads", no_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:oi:cpl:e:mt:S",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'i':
            aln_file = optarg;
            break;

        case 'c':
            compact_ids = true;
            break;

        case 'k':
            path_name = optarg;
            break;

        case 'o':
            remove_orphans = true;
            break;

        case 'p':
            prune_complex = true;
            break;

        case 'S':
            prune_subgraphs = true;
            break;

        case 'l':
            path_length = atoi(optarg);
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;

        case 'm':
            add_start_and_end_markers = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
            help_mod(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (!path_name.empty()) {
        graph->keep_path(path_name);
    }

    if (remove_orphans) {
        graph->remove_orphan_edges();
    }

    if (!aln_file.empty()) {
        //void edit(const map<int64_t, vector<tuple<Mapping, int64_t, int64_t> > >& mappings);
        
        vector<Path> paths;
        function<void(Alignment&)> lambda = [&graph, &paths](Alignment& aln) {
            Path path = simplify_deletions(aln.path());
            path.set_name(aln.name());
            paths.push_back(path);
        };
        if (aln_file == "-") {
            stream::for_each(std::cin, lambda);
        } else {
            ifstream in;
            in.open(aln_file.c_str());
            stream::for_each(in, lambda);
        }
        // now that everything is sorted, execute the edits
        graph->edit(paths);
        // and optionally compact ids
        if (compact_ids) {
            graph->sort();
            graph->compact_ids();
        }
    }

    if (prune_complex) {
        if (!(path_length > 0 && edge_max > 0)) {
            cerr << "[vg mod]: when pruning complex regions you must specify a --path-length and --edge-max" << endl;
            return 1;
        }
        Node* head_node = NULL;
        Node* tail_node = NULL;
        graph->add_start_and_end_markers(path_length, '#', '$', head_node, tail_node);
        graph->prune_complex(path_length, edge_max, head_node, tail_node);
        // These nodes were created in the graph, so we can destroy them by pointer.
        graph->destroy_node(head_node);
        graph->destroy_node(tail_node);
    }

    if (prune_subgraphs) {
        graph->prune_short_subgraphs(path_length);
    }

    if (add_start_and_end_markers) {
        if (!(path_length > 0)) {
            cerr << "[vg mod]: when adding start and end markers you must provide a --path-length" << endl;
            return 1;
        }
        Node* head_node = NULL;
        Node* tail_node = NULL;
        graph->add_start_and_end_markers(path_length, '#', '$', head_node, tail_node);
    }

    graph->serialize_to_ostream(std::cout);

    delete graph;

    return 0;
}

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options] <graph.vg>" << endl
         << "Simulates reads from the graph(s). Output is a list of reads." << endl
         << endl
         << "options:" << endl
         << "    -l, --read-length N   write reads of length N" << endl
         << "    -n, --num-reads N     simulate N reads" << endl
         << "    -s, --random-seed N   use this specific seed for the PRNG" << endl
         << "    -e, --base-error N    base substitution error rate (default 0.0)" << endl
         << "    -i, --indel-error N   indel error rate (default 0.0)" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    int read_length = 100;
    int num_reads = 1;
    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"read-length", required_argument, 0, 'l'},
                {"num-reads", required_argument, 0, 'n'},
                {"random-seed", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'l':
            read_length = atoi(optarg);
            break;

        case 'n':
            num_reads = atoi(optarg);
            break;

        case 's':
            seed_val = atoi(optarg);
            break;

        case 'e':
            base_error = atof(optarg);
            break;

        case 'i':
            indel_error = atof(optarg);
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

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    int64_t max_id = graph->max_node_id();
    int64_t min_id = graph->min_node_id();

    mt19937 rng;
    rng.seed(seed_val);

    string bases = "ATGC";
    uniform_real_distribution<double> rprob(0, 1);
    uniform_int_distribution<int> rbase(0, 3);

    function<string(const string&)> introduce_read_errors
        = [&rng, &rprob, &rbase, &bases, base_error, indel_error](const string& perfect_read) {

        if (base_error == 0 && indel_error == 0) return perfect_read;
        string read;
        for (auto c : perfect_read) {
            if (rprob(rng) <= base_error) {
                // pick another base than what c is
                char e;
                do {
                    e = bases[rbase(rng)];
                } while (e == c);
                c = e;
            }
            if (rprob(rng) <= indel_error) {
                if (rprob(rng) < 0.5) {
                    read.push_back(bases[rbase(rng)]);
                } // else do nothing, == deletion of base
            } else {
                read.push_back(c);
            }
        }
        return read;
    };

    for (int i = 0; i < num_reads; ++i) {
        string perfect_read = graph->random_read(read_length, rng, min_id, max_id, true);
        // avoid short reads at the end of the graph by retrying
        int iter = 0;
        while (perfect_read.size() < read_length && ++iter < 1000) {
            perfect_read = graph->random_read(read_length, rng, min_id, max_id, true);
            // if we can't make a suitable read in 1000 tries, then maybe the graph is too small?
        }
        if (iter == 1000) {
            cerr << "couldn't simulate read, perhaps the chosen length is too long for this graph?" << endl;
        } else {
            // apply errors
            string readseq = introduce_read_errors(perfect_read);
            cout << readseq << endl;
        }
    }
    delete graph;

    return 0;
}

void help_kmers(char** argv) {
    cerr << "usage: " << argv[0] << " kmers [options] <graph1.vg> [graph2.vg ...] >kmers.tsv" << endl
         << "Generates kmers of the graph(s). Output is: kmer id pos" << endl
         << endl
         << "options:" << endl
         << "    -k, --kmer-size N     print kmers of size N in the graph" << endl
         << "    -e, --edge-max N      only consider paths which cross this many potential alternate edges" << endl
         << "                          (e.g. if node out-degree is 2, we would count 1 toward --edge-max," << endl
         << "                          for 3 we would count 2)" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers in paths (default 1)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -d, --ignore-dups     filter out duplicated kmers" << endl
         << "    -n, --allow-negs      don't filter out relative negative positions of kmers" << endl
         << "    -g, --gcsa-out        output a table suitable for input to GCSA2" << endl
         << "                          kmer, starting position, previous characters," << endl
         << "                          successive characters, successive positions" << endl
         << "    -p, --progress        show progress" << endl;
}

int main_kmers(int argc, char** argv) {

    if (argc == 2) {
        help_kmers(argv);
        return 1;
    }

    int kmer_size = 0;
    int edge_max = 0;
    int kmer_stride = 1;
    bool show_progress = false;
    bool gcsa_out = false;
    bool allow_dups = true;
    bool allow_negs = false;
    // for distributed GCSA2 kmer generation
    int64_t head_id = 0;
    int64_t tail_id = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"kmer-size", required_argument, 0, 'k'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"edge-max", required_argument, 0, 'e'},
                {"threads", required_argument, 0, 't'},
                {"gcsa-out", no_argument, 0, 'g'},
                {"ignore-dups", no_argument, 0, 'd'},
                {"allow-negs", no_argument, 0, 'n'},
                {"progress",  no_argument, 0, 'p'},
                {"head-id", required_argument, 0, 'H'},
                {"tail-id", required_argument, 0, 'T'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:j:pt:e:gdnH:T:",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {

        case 'k':
            kmer_size = atoi(optarg);
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'g':
            gcsa_out = true;
            break;

        case 'd':
            allow_dups = false;
            break;

        case 'n':
            allow_negs = true;
            break;

        case 'p':
            show_progress = true;
            break;

        case 'H':
            head_id = atoi(optarg);
            break;

        case 'T':
            tail_id = atoi(optarg);
            break;

        case 'h':
        case '?':
            help_kmers(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if (edge_max == 0) edge_max = kmer_size + 1;

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);

    graphs.show_progress = show_progress;

    if (gcsa_out) {
        graphs.write_gcsa_out(cout, kmer_size, edge_max, kmer_stride, allow_dups, head_id, tail_id);
    } else {
        function<void(string&, NodeTraversal, int, list<NodeTraversal>&, VG& graph)>
            lambda = [](string& kmer, NodeTraversal n, int p, list<NodeTraversal>& path, VG& graph) {
            // We encode orientation by negating the IDs for backward nodes.
            // Their offsets are from the end of the node in its local forward
            // orientation, and are negated in the output.
            int sign = n.backward ? -1 : 1;            
#pragma omp critical (cout)
            cout << kmer << '\t' << n.node->id() * sign << '\t' << p * sign << '\n';
        };
        graphs.for_each_kmer_parallel(lambda, kmer_size, edge_max, kmer_stride, allow_dups, allow_negs);
    }
    cout.flush();

    return 0;
}

void help_concat(char** argv) {
    cerr << "usage: " << argv[0] << " concat [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Concatenates graphs in order by adding edges from the tail nodes of the" << endl
         << "predecessor to the head nodes of the following graph. Node IDs are" << endl
         << "compacted, so care should be taken if consistent IDs are required." << endl;
}

int main_concat(int argc, char** argv) {

    if (argc == 2) {
        help_concat(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'h':
        case '?':
            help_concat(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    list<VG*> graphs;

    while (optind < argc) {
        VG* graph;
        string file_name = argv[optind++];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        graphs.push_back(graph);
    }

    VG merged;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        merged.append(**g);
    }

    // output
    merged.serialize_to_ostream(std::cout);

    return 0;
}

void help_ids(char** argv) {
    cerr << "usage: " << argv[0] << " ids [options] <graph1.vg> [graph2.vg ...] >new.vg" << endl
         << "options:" << endl
         << "    -c, --compact        minimize the space of integers used by the ids" << endl
         << "    -i, --increment N    increase ids by N" << endl
         << "    -d, --decrement N    decrease ids by N" << endl
         << "    -j, --join           make a joint id space for all the graphs that are supplied" << endl
         << "                         by iterating through the supplied graphs and incrementing" << endl
         << "                         their ids to be non-conflicting" << endl;
}

int main_ids(int argc, char** argv) {

    if (argc == 2) {
        help_ids(argv);
        return 1;
    }

    bool join = false;
    bool compact = false;
    int64_t increment = 0;
    int64_t decrement = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"compact", no_argument, 0, 'c'},
                {"increment", required_argument, 0, 'i'},
                {"decrement", required_argument, 0, 'd'},
                {"join", no_argument, 0, 'j'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hci:d:j",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'c':
            compact = true;
            break;

        case 'i':
            increment = atoi(optarg);
            break;

        case 'd':
            decrement = atoi(optarg);
            break;

        case 'j':
            join = true;
            break;

        case 'h':
        case '?':
            help_ids(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (!join) {
        VG* graph;
        string file_name = argv[optind];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }

        if (compact) {
            graph->compact_ids();
        }

        if (increment != 0) {
            graph->increment_node_ids(increment);
        }

        if (decrement != 0) {
            graph->decrement_node_ids(decrement);
        }

        graph->serialize_to_ostream(std::cout);
        delete graph;
    } else {

        vector<string> graph_file_names;
        while (optind < argc) {
            VG* graph;
            string file_name = argv[optind++];
            graph_file_names.push_back(file_name);
        }

        VGset graphs(graph_file_names);
        graphs.merge_id_space();

    }

    return 0;

}

void help_join(char** argv) {
    cerr << "usage: " << argv[0] << " join [options] <graph1.vg> [graph2.vg ...] >joined.vg" << endl
         << "Joins graphs and sub-graphs into a single variant graph by connecting their" << endl
         << "heads to a single root node with sequence 'N'." << endl
         << "Assumes a single id namespace for all graphs to join." << endl;
}

int main_join(int argc, char** argv) {

    if (argc == 2) {
        help_join(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'h':
        case '?':
            help_join(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    list<VG*> graphs;

    while (optind < argc) {
        VG* graph;
        string file_name = argv[optind++];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        graphs.push_back(graph);
    }

    VG joined;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        joined.extend(**g);
    }

    // combine all subgraphs
    joined.join_heads();

    // output
    joined.serialize_to_ostream(std::cout);

    return 0;
}

void help_stats(char** argv) {
    cerr << "usage: " << argv[0] << " stats [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -z, --size            size of graph" << endl
         << "    -l, --length          length of sequences in graph" << endl
         << "    -s, --subgraphs       describe subgraphs of graph" << endl
         << "    -H, --heads           list the head nodes of the graph" << endl
         << "    -T, --tails           list the tail nodes of the graph" << endl;
}

int main_stats(int argc, char** argv) {

    if (argc == 2) {
        help_stats(argv);
        return 1;
    }

    bool stats_size = false;
    bool stats_length = false;
    bool stats_subgraphs = false;
    bool stats_heads = false;
    bool stats_tails = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"size", no_argument, 0, 'z'},
                {"length", no_argument, 0, 'l'},
                {"subgraphs", no_argument, 0, 's'},
                {"heads", no_argument, 0, 'H'},
                {"tails", no_argument, 0, 'T'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlsHT",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'z':
            stats_size = true;
            break;

        case 'l':
            stats_length = true;
            break;

        case 's':
            stats_subgraphs = true;
            break;

        case 'H':
            stats_heads = true;
            break;

        case 'T':
            stats_tails = true;
            break;

        case 'h':
        case '?':
            help_stats(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (stats_size) {
        cout << "nodes" << "\t" << graph->node_count() << endl
             << "edges" << "\t" << graph->edge_count() << endl;
    }

    if (stats_length) {
        cout << "length" << "\t" << graph->total_length_of_nodes() << endl;
    }

    if (stats_heads) {
        vector<Node*> heads;
        graph->head_nodes(heads);
        cout << "heads" << "\t";
        for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
            cout << (*h)->id() << " ";
        }
        cout << endl;
    }

    if (stats_tails) {
        vector<Node*> tails;
        graph->tail_nodes(tails);
        cout << "tails" << "\t";
        for (vector<Node*>::iterator t = tails.begin(); t != tails.end(); ++t) {
            cout << (*t)->id() << " ";
        }
        cout << endl;
    }

    if (stats_subgraphs) {
        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);
        // these are topologically-sorted
        for (list<VG>::iterator s = subgraphs.begin(); s != subgraphs.end(); ++s) {
            VG& subgraph = *s;
            vector<Node*> heads;
            subgraph.head_nodes(heads);
            int64_t length = subgraph.total_length_of_nodes();
            for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
                cout << (h==heads.begin()?"":",") << (*h)->id();
            }
            cout << "\t" << length << endl;
        }
    }

    delete graph;

    return 0;

}

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -n, --node ID         starting at node with ID" << endl
         << "    -l, --max-length N    generate paths of at most length N" << endl
         << "    -e, --edge-max N      only consider paths which cross this many potential alternate edges" << endl
         << "                          (e.g. if node out-degree is 2, we would count 1 toward --edge-max," << endl
         << "                          for 3 we would count 2)" << endl
         << "    -s, --as-seqs         write each path as a sequence" << endl;
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    int max_length = 0;
    int edge_max = 0;
    int64_t node_id = 0;
    bool as_seqs = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"node", required_argument, 0, 'n'},
                {"max-length", required_argument, 0, 'l'},
                {"edge-max", required_argument, 0, 'e'},
                {"as-seqs", no_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:l:hse:",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'n':
            node_id = atoll(optarg);
            break;

        case 'l':
            max_length = atoi(optarg);
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;

        case 's':
            as_seqs = true;
            break;

        case 'h':
        case '?':
            help_paths(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if (edge_max == 0) edge_max = max_length + 1;

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (max_length == 0) {
        cerr << "error:[vg paths] a --max-length is required when generating paths" << endl;
    }

    function<void(Node*,Path&)> paths_to_seqs = [graph](Node* n, Path& p) {
        string seq = graph->path_sequence(p);
#pragma omp critical(cout)
        cout << seq << endl;
    };

    function<void(Node*,Path&)> paths_to_json = [](Node* n, Path& p) {
        string json2 = pb2json(p);
#pragma omp critical(cout)
        cout<<json2<<endl;
    };

    function<void(Node*, Path&)>* callback = &paths_to_seqs;
    if (!as_seqs) {
        callback = &paths_to_json;
    }

    auto noop = [](NodeTraversal) { }; // don't handle the failed regions of the graph yet

    if (node_id) {
        graph->for_each_kpath_of_node(graph->get_node(node_id), max_length, edge_max,
                                      noop, noop,
                                      *callback);
    } else {
        graph->for_each_kpath_parallel(max_length, edge_max,
                                       noop, noop,
                                       *callback);
    }

    delete graph;

    return 0;

}

void help_find(char** argv) {
    cerr << "usage: " << argv[0] << " find [options] <graph.vg> >sub.vg" << endl
         << "options:" << endl
         << "    -n, --node ID          find node, return 1-hop context as graph" << endl
         << "    -f, --edges-from ID    return edges from node with ID" << endl
         << "    -t, --edges-to ID      return edges from node with ID" << endl
         << "    -k, --kmer STR         return a graph of edges and nodes matching this kmer" << endl
         << "    -T, --table            instead of a graph, return a table of kmers" << endl
         << "    -c, --context STEPS    expand the context of the kmer hit subgraphs" << endl
         << "    -s, --sequence STR     search for sequence STR using --kmer-size kmers" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in sequence (default 1)" << endl
         << "    -z, --kmer-size N      split up --sequence into kmers of size N" << endl
         << "    -C, --kmer-count       report approximate count of kmer (-k) in db" << endl
         << "    -p, --path TARGET      find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]" << endl
         << "    -P, --position-in PATH find the position of the node (specified by -n) in the given path" << endl
         << "    -r, --node-range N:M   get nodes from N to M" << endl
        //<< "    -a, --alignments       write all stored alignments in sorted order (in GAM)" << endl
        //<< "    -m, --mappings         write stored mappings in sorted order (in json)" << endl
         << "    -d, --db-name DIR      use this db (defaults to <graph>.index/)" << endl;
}

int main_find(int argc, char** argv) {

    if (argc == 2) {
        help_find(argv);
        return 1;
    }

    string db_name;
    string sequence;
    int kmer_size=0;
    int kmer_stride = 1;
    vector<string> kmers;
    string output_format;
    int64_t from_id=0, to_id=0;
    vector<int64_t> node_ids;
    int context_size=0;
    bool count_kmers = false;
    bool kmer_table = false;
    string target;
    string path_name;
    string range;
    bool get_alignments = false;
    bool get_mappings = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"node", required_argument, 0, 'n'},
                {"edges-from", required_argument, 0, 'f'},
                {"edges-to", required_argument, 0, 't'},
                {"kmer", required_argument, 0, 'k'},
                {"table", no_argument, 0, 'T'},
                {"sequence", required_argument, 0, 's'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'z'},
                {"output", required_argument, 0, 'o'},
                {"context", required_argument, 0, 'c'},
                {"kmer-count", no_argument, 0, 'C'},
                {"path", required_argument, 0, 'p'},
                {"position-in", required_argument, 0, 'P'},
                {"node-range", required_argument, 0, 'r'},
                {"alignments", no_argument, 0, 'a'},
                {"mappings", no_argument, 0, 'm'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:n:f:t:o:k:hc:s:z:j:CTp:P:r:am",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            db_name = optarg;
            break;

        case 'k':
            kmers.push_back(optarg);
            break;

        case 's':
            sequence = optarg;
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'z':
            kmer_size = atoi(optarg);
            break;

        case 'C':
            count_kmers = true;
            break;

        case 'p':
            target = optarg;
            break;

        case 'P':
            path_name = optarg;
            break;

        case 'c':
            context_size = atoi(optarg);
            break;

        case 'n':
            node_ids.push_back(atoi(optarg));
            break;

        case 'f':
            from_id = atoi(optarg);
            break;

        case 't':
            to_id = atoi(optarg);
            break;

        case 'T':
            kmer_table = true;
            break;

        case 'r':
            range = optarg;
            break;

        case 'a':
            get_alignments = true;
            break;

        case 'm':
            get_mappings = true;
            break;

        case 'o':
            output_format = optarg;
            break;

        case 'h':
        case '?':
            help_find(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if (optind < argc) {
        string file_name = argv[optind];
        if (file_name == "-") {
            if (db_name.empty()) {
                cerr << "error:[vg find] reading variant graph from stdin and no db name (-d) given, exiting" << endl;
                return 1;
            }
        }
        ifstream in;
        if (db_name.empty()) {
            db_name = file_name + ".index";
        }
        in.open(file_name.c_str());
    }

    Index index;
    // open index
    index.open_read_only(db_name);

    if (get_alignments) {
    }

    if (!node_ids.empty() && path_name.empty()) {
        // get the context of the node
        vector<VG> graphs;
        for (auto node_id : node_ids) {
            VG g;
            index.get_context(node_id, g);
            if (context_size > 0) {
                index.expand_context(g, context_size);
            }
            graphs.push_back(g);
        }
        VG result_graph;
        for (auto& graph : graphs) {
            result_graph.extend(graph);
        }
        result_graph.remove_orphan_edges();
        // return it
        result_graph.serialize_to_ostream(cout);
    } else if (from_id != 0) {
        vector<Edge> edges;
        index.get_edges_from(from_id, edges);
        for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
            cout << e->from() << "\t" << e->to() << endl;
        }
    } else if (to_id != 0) {
        vector<Edge> edges;
        index.get_edges_to(to_id, edges);
        for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
            cout << e->from() << "\t" << e->to() << endl;
        }
    }

    if (!node_ids.empty() && !path_name.empty()) {
        int64_t path_id = index.get_path_id(path_name);
        for (auto node_id : node_ids) {
            list<int64_t> path_prev, path_next;
            int64_t prev_pos=0, next_pos=0;
            if (index.get_node_path_relative_position(node_id, path_id,
                                                      path_prev, prev_pos,
                                                      path_next, next_pos)) {

                cout << node_id << "\t" << path_prev.front() << "\t" << prev_pos
                     << "\t" << path_next.back() << "\t" << next_pos << "\t";

                Mapping m = index.path_relative_mapping(node_id, path_id,
                                                        path_prev, prev_pos,
                                                        path_next, next_pos);
                cout << pb2json(m) << endl;
            }
        }
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        VG graph;
        parse_region(target, name, start, end);
        index.get_path(graph, name, start, end);
        if (context_size > 0) {
            index.expand_context(graph, context_size);
        }
        graph.serialize_to_ostream(cout);
    }

    if (!range.empty()) {
        VG graph;
        int64_t id_start=0, id_end=0;
        vector<string> parts = split_delims(range, ":");
        if (parts.size() == 1) {
            cerr << "[vg find] error, format of range must be \"N:M\" where start id is N and end id is M, got " << range << endl;
            exit(1);
        }
        convert(parts.front(), id_start);
        convert(parts.back(), id_end);
        index.get_range(id_start, id_end, graph);
        if (context_size > 0) {
            index.expand_context(graph, context_size);
        }
        graph.remove_orphan_edges();
        graph.serialize_to_ostream(cout);
    }

    // todo cleanup if/else logic to allow only one function
    
    if (!sequence.empty()) {
        set<int> kmer_sizes = index.stored_kmer_sizes();
        if (kmer_sizes.empty()) {
            cerr << "error:[vg find] index does not include kmers, add with vg index -k" << endl;
            return 1;
        }
        if (kmer_size == 0) {
            kmer_size = *kmer_sizes.begin();
        }
        for (int i = 0; i <= sequence.size()-kmer_size; i+=kmer_stride) {
            kmers.push_back(sequence.substr(i,kmer_size));
        }
    }

    if (!kmers.empty()) {
        if (count_kmers) {
            for (auto& kmer : kmers) {
                cout << kmer << "\t" << index.approx_size_of_kmer_matches(kmer) << endl;
            }
        } else if (kmer_table) {
            for (auto& kmer : kmers) {
                map<string, vector<pair<int64_t, int32_t> > > positions;
                index.get_kmer_positions(kmer, positions);
                for (auto& k : positions) {
                    for (auto& p : k.second) {
                        cout << k.first << "\t" << p.first << "\t" << p.second << endl;
                    }
                }
            }
        } else {
            vector<VG> graphs;
            for (auto& kmer : kmers) {
                VG g;
                index.get_kmer_subgraph(kmer, g);
                if (context_size > 0) {
                    index.expand_context(g, context_size);
                }
                graphs.push_back(g);
            }

            VG result_graph;
            for (auto& graph : graphs) {
                result_graph.extend(graph);
            }
            result_graph.remove_orphan_edges();
            result_graph.serialize_to_ostream(cout);
        }
    }

    return 0;

}

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space, and the graph containing the highest-ID node " << endl 
         << "must come first." << endl
         << "general options:" << endl
         << "    -g, --gcsa-out         output a GCSA2 index instead of a rocksdb index" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph" << endl
         << "    -e, --edge-max N       only consider paths which cross this many potential alternate edges" << endl
         << "                           (e.g. if node out-degree is 2, we would count 1 toward --edge-max," << endl
         << "                           for 3 we would count 2)" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in paths (default 1)" << endl
         << "    -d, --db-name PATH     create rocksdb in PATH directory (default: <graph>.index/)" << endl
         << "                           or GCSA2 index in PATH file (default: <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "                           (this is required if you are using multiple graphs files)" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "rocksdb options (ignored with -g):" << endl
         << "    -s, --store-graph      store graph (do this first to build db!)" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -P, --prune KB         remove kmer entries which use more than KB kilobytes" << endl
         << "    -n, --allow-negs       don't filter out relative negative positions of kmers" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -M, --metadata         describe aspects of the db stored in metadata" << endl
         << "    -L, --path-layout      describes the path layout of the graph" << endl
         << "    -S, --set-kmer         assert that the kmer size (-k) is in the db" << endl
        //<< "    -b, --tmp-db-base S    use this base name for temporary indexes" << endl
         << "    -C, --compact          compact the index into a single level (improves performance)" << endl
         << "    -Q, --use-snappy       use snappy compression (faster, larger) rather than zlib" << endl;
         
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    string db_name;
    int kmer_size = 0;
    int edge_max = 0;
    int kmer_stride = 1;
    int prune_kb = -1;
    bool store_graph = false;
    bool dump_index = false;
    bool describe_index = false;
    bool show_progress = false;
    bool set_kmer_size = false;
    bool path_layout = false;
    bool store_alignments = false;
    bool store_mappings = false;
    bool allow_negs = false;
    bool compact = false;
    bool dump_alignments = false;
    bool use_snappy = false;
    bool gcsa_out = false;
    int doubling_steps = gcsa::GCSA::DOUBLING_STEPS; // TODO: add an option for this?

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"kmer-size", required_argument, 0, 'k'},
                {"edge-max", required_argument, 0, 'e'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"store-graph", no_argument, 0, 's'},
                {"store-alignments", no_argument, 0, 'a'},
                {"dump-alignments", no_argument, 0, 'A'},
                {"store-mappings", no_argument, 0, 'm'},
                {"dump", no_argument, 0, 'D'},
                {"metadata", no_argument, 0, 'M'},
                {"set-kmer", no_argument, 0, 'S'},
                {"threads", required_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"prune",  required_argument, 0, 'P'},
                {"path-layout", no_argument, 0, 'L'},
                {"compact", no_argument, 0, 'C'},
                {"allow-negs", no_argument, 0, 'n'},
                {"use-snappy", no_argument, 0, 'Q'},
                {"gcsa-out", no_argument, 0, 'g'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:k:j:pDshMt:b:e:SP:LmaCnAQg",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            db_name = optarg;
            break;

        case 'P':
            prune_kb = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'D':
            dump_index = true;
            break;

        case 'M':
            describe_index = true;
            break;

        case 'L':
            path_layout = true;
            break;

        case 'S':
            set_kmer_size = true;
            break;

        case 's':
            store_graph = true;
            break;

        case 'a':
            store_alignments = true;
            break;

        case 'A':
            dump_alignments = true;
            break;

        case 'm':
            store_mappings = true;
            break;

        case 'n':
            allow_negs = true;
            break;

        case 'C':
            compact = true;
            break;

        case 'Q':
            use_snappy = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;
            
        case 'g':
            gcsa_out = true;
            break;
 
        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if (edge_max == 0) edge_max = kmer_size + 1;

    vector<string> file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        file_names.push_back(file_name);
    }

    if (db_name.empty()) {
        if (file_names.size() > 1) {
            cerr << "error:[vg index] working on multiple graphs and no db name (-d) given, exiting" << endl;
            return 1;
        } else if (file_names.size() == 1) {
            if(gcsa_out) {
                // Name the database for gcsa
                db_name = *file_names.begin() + gcsa::GCSA::EXTENSION;
            } else {
                // Name the database for rocksdb
                db_name = *file_names.begin() + ".index";
            }
        } else {
            cerr << "error:[vg index] no graph or db given, exiting" << endl;
            return 1;
        }
    }
    
    if(kmer_size == 0 && gcsa_out) {
        // gcsa doesn't do anything if we tell it a kmer size of 0.
        cerr << "error:[vg index] kmer size for GCSA2 index must be >0" << endl;
        return 1;
    }
    
    if(kmer_size < 0) {
        cerr << "error:[vg index] kmer size cannot be negative" << endl;
        return 1;
    }
    
    if(kmer_stride <= 0) {
        // kmer strides of 0 (or negative) are silly.
        cerr << "error:[vg index] kmer stride must be positive and nonzero" << endl;
        return 1;
    }

    if(gcsa_out) {
        // We need to make a gcsa index.
    
        // Load up the graphs
        VGset graphs(file_names);

        graphs.show_progress = show_progress;

        // Go get the kmers of the correct size
        vector<gcsa::KMer> kmers;
        graphs.get_gcsa_kmers(kmer_size, edge_max, kmer_stride, kmers, true, 0, 0);
        
        // Handle finding the sink node
        size_t sink_node_id = 0;
        for(auto kmer : kmers) {
            if(gcsa::Key::label(kmer.key) == 0) {
                // This kmer lets us know the sink node.
                sink_node_id = gcsa::Node::id(kmer.from);
                break;
            }
        }
        
        for(auto& kmer : kmers) {
            // Mark kmers that go to the sink node as "sorted", since they have stop
            // characters in them and can't be extended.
            // TODO: Can we just check for the presence of "$" during conversion and not do this serial loop?
            if(gcsa::Node::id(kmer.to) == sink_node_id && gcsa::Node::offset(kmer.to) > 0) {
                kmer.makeSorted();
            }
            
            //cout << kmer << std::endl;
        }
        
        if(show_progress) {
            cerr << "Found " << kmers.size() << " kmer instances" << endl;
        }
        
        // Make the index with the kmers
        gcsa::GCSA gcsa_index(kmers, kmer_size, doubling_steps);
        
        // Save it to the index filename
        sdsl::store_to_file(gcsa_index, db_name);
        
        // Skip all the Snappy stuff we can't do (yet).
        return 0;
    }
    

    Index index;
    index.use_snappy = use_snappy;

    if (compact) {
        index.open_for_write(db_name);
        index.compact();
        index.flush();
        index.close();
    }

    // todo, switch to xg for graph storage
    // index should write and load index/xg or such
    // then a handful of functions used in main.cpp and mapper.cpp need to be rewritten to use the xg index
    if (store_graph && file_names.size() > 0) {
        index.open_for_write(db_name);
        VGset graphs(file_names);
        graphs.show_progress = show_progress;
        graphs.store_in_index(index);
        //index.flush();
        //index.close();
        // reopen to index paths
        // this requires the index to be queryable
        //index.open_for_write(db_name);
        graphs.store_paths_in_index(index);
        index.compact();
        index.flush();
        index.close();
    }

    if (store_alignments && file_names.size() > 0) {
        index.open_for_write(db_name);
        function<void(Alignment&)> lambda = [&index](Alignment& aln) {
            index.put_alignment(aln);
        };
        for (auto& file_name : file_names) {
            if (file_name == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each(in, lambda);
            }
        }
        index.flush();
        index.close();
    }

    if (dump_alignments) {
        vector<Alignment> output_buf;
        index.open_read_only(db_name);
        auto lambda = [&output_buf](const Alignment& aln) {
            output_buf.push_back(aln);
            stream::write_buffered(cout, output_buf, 1000);
        };
        index.for_each_alignment(lambda);
        stream::write_buffered(cout, output_buf, 0);
        index.close();
    }

    if (store_mappings && file_names.size() > 0) {
        index.open_for_write(db_name);
        function<void(Alignment&)> lambda = [&index](Alignment& aln) {
            const Path& path = aln.path();
            for (int i = 0; i < path.mapping_size(); ++i) {
                index.put_mapping(path.mapping(i));
            }
        };
        for (auto& file_name : file_names) {
            if (file_name == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each(in, lambda);
            }
        }
        index.flush();
        index.close();
    }

    if (kmer_size != 0 && file_names.size() > 0) {
        index.open_for_bulk_load(db_name);
        VGset graphs(file_names);
        graphs.show_progress = show_progress;
        graphs.index_kmers(index, kmer_size, edge_max, kmer_stride, allow_negs);
        index.flush();
        index.close();
        // forces compaction
        index.open_for_write(db_name);
        index.flush();
        index.compact();
        index.close();
    }

    if (prune_kb >= 0) {
        if (show_progress) {
            cerr << "pruning kmers > " << prune_kb << " on disk from " << db_name << endl;
        }
        index.open_for_write(db_name);
        index.prune_kmers(prune_kb);
        index.compact();
        index.close();
    }

    if (set_kmer_size) {
        assert(kmer_size != 0);
        index.open_for_write(db_name);
        index.remember_kmer_size(kmer_size);
        index.close();
    }

    if (dump_index) {
        index.open_read_only(db_name);
        index.dump(cout);
        index.close();
    }

    if (describe_index) {
        index.open_read_only(db_name);
        set<int> kmer_sizes = index.stored_kmer_sizes();
        cout << "kmer sizes: ";
        for (auto kmer_size : kmer_sizes) {
            cout << kmer_size << " ";
        }
        cout << endl;
        index.close();
    }

    if (path_layout) {
        index.open_read_only(db_name);
        //index.path_layout();
        map<string, int64_t> path_by_id = index.paths_by_id();
        map<string, pair<int64_t, int64_t> > layout;
        map<string, int64_t> length;
        index.path_layout(layout, length);
        for (auto& p : layout) {
            cout << p.first << " " << p.second.first << " " << p.second.second << " " << length[p.first] << endl;
        }
        index.close();
    }

    return 0;

}

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg> >alignments.gam" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value" << endl
         << "    -j, --json            output alignments in JSON format (default GAM)" << endl;
}

int main_align(int argc, char** argv) {

    string seq;
    string seq_name;

    if (argc == 2) {
        help_align(argv);
        return 1;
    }

    bool print_cigar = false;
    bool output_json = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"seq-name", no_argument, 0, 'Q'},
                {"json", no_argument, 0, 'j'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhQ:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'Q':
            seq_name = optarg;
            break;

        case 'j':
            output_json = true;
            break;
 
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_align(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    Alignment alignment = graph->align(seq);

    if (output_json) {
        cout << pb2json(alignment) << endl;
    } else {
        if (!seq_name.empty()) {
            alignment.set_name(seq_name);
        }
        function<Alignment(uint64_t)> lambda =
            [&alignment] (uint64_t n) {
            return alignment;
        };
        stream::write(cout, 1, lambda);
    }

    delete graph;

    return 0;

}

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] <graph.vg> >alignments.vga" << endl
         << "options:" << endl
         << "    -d, --db-name DIR     use this db (defaults to <graph>.index/)" << endl
         << "                          a graph is not required" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -r, --reads FILE      take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE  align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -f, --fastq FILE      input fastq (possibly compressed), two are allowed, one for each mate" << endl
         << "    -i, --interleaved     fastq is interleaved paired-ended" << endl
         << "    -p, --pair-window N   align to a graph up to N ids away from the mapping location of one mate for the other" << endl
        //<< "    -B, --try-both-mates  attempt to align both reads individually, then used paired end resolution to fix" << endl
         << "    -N, --sample NAME     for --reads input, add this sample" << endl
         << "    -R, --read-group NAME for --reads input, add this read group" << endl
         << "    -k, --kmer-size N     use this kmer size, it must be < kmer size in db (default: from index)" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers to use for seeding (default: kmer size)" << endl
         << "    -E, --min-kmer-entropy N  require shannon entropy of this in order to use kmer (default: no limit)" << endl
         << "    -S, --sens-step N     decrease kmer size by N bp until alignment succeeds (default: 5)" << endl
         << "    -A, --max-attempts N  try to improve sensitivity and align this many times (default: 7)" << endl
         << "    -x, --thread-ex N     grab this many neighboring nodes around each thread for alignment (default: 2)" << endl
         << "    -c, --clusters N      use at most the largest N ordered clusters of the kmer graph for alignment (default: all)" << endl
         << "    -C, --cluster-min N   require at least this many kmer hits in a cluster to attempt alignment (default: 2)" << endl
         << "    -m, --hit-max N       ignore kmers who have >N hits in our index (default: 100)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -F, --prefer-forward  if the forward alignment of the read works, accept it" << endl
         << "    -G, --greedy-accept   if a tested alignment achieves -X score/bp don't try worse seeds" << endl
         << "    -X, --score-per-bp N  accept early alignment if the alignment score per base is > N and -F or -G is set" << endl
         << "    -J, --output-json     output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -B, --band-width N    for very long sequences, align in chunks then merge paths (default 1000bp)" << endl
         << "    -D, --debug           print debugging information about alignment to stderr" << endl;
}

int main_map(int argc, char** argv) {

    if (argc == 2) {
        help_map(argv);
        return 1;
    }

    string seq;
    string seq_name;
    string db_name;
    int kmer_size = 0;
    int kmer_stride = 0;
    int sens_step = 0;
    int best_clusters = 0;
    int cluster_min = 2;
    int max_attempts = 7;
    string read_file;
    string hts_file;
    int hit_max = 100;
    int thread_count = 1;
    int thread_ex = 2;
    bool output_json = false;
    bool debug = false;
    bool prefer_forward = false;
    bool greedy_accept = false;
    float score_per_bp = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_fastq = false;
    int pair_window = 64; // ~11bp/node
    int band_width = 1000; // anything > 1000bp sequences is difficult to align efficiently
    bool try_both_mates_first = false;
    float min_kmer_entropy = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"seq-name", required_argument, 0, 'Q'},
                {"db-name", required_argument, 0, 'd'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'k'},
                {"min-kmer-entropy", required_argument, 0, 'E'},
                {"clusters", required_argument, 0, 'c'},
                {"cluster-min", required_argument, 0, 'C'},
                {"max-attempts", required_argument, 0, 'A'},
                {"reads", required_argument, 0, 'r'},
                {"sample", required_argument, 0, 'N'},
                {"read-group", required_argument, 0, 'R'},
                {"hit-max", required_argument, 0, 'm'},
                {"threads", required_argument, 0, 't'},
                {"prefer-forward", no_argument, 0, 'F'},
                {"greedy-accept", no_argument, 0, 'G'},
                {"score-per-bp", required_argument, 0, 'X'},
                {"sens-step", required_argument, 0, 'S'},
                {"thread-ex", required_argument, 0, 'x'},
                {"output-json", no_argument, 0, 'J'},
                {"hts-input", no_argument, 0, 'b'},
                {"fastq", no_argument, 0, 'f'},
                {"interleaved", no_argument, 0, 'i'},
                {"pair-window", required_argument, 0, 'p'},
                {"band-width", required_argument, 0, 'B'},
                {"debug", no_argument, 0, 'D'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:j:hd:c:r:m:k:t:DX:FS:Jb:R:N:if:p:B:x:GC:A:E:Q:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'Q':
            seq_name = optarg;
            break;

        case 'd':
            db_name = optarg;
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;

        case 'S':
            sens_step = atoi(optarg);
            break;

        case 'c':
            best_clusters = atoi(optarg);
            break;

        case 'C':
            cluster_min = atoi(optarg);
            break;

        case 'E':
            min_kmer_entropy = atof(optarg);
            break;

        case 'A':
            max_attempts = atoi(optarg);
            break;

        case 'm':
            hit_max = atoi(optarg);
            break;

        case 'x':
            thread_ex = atoi(optarg);
            break;

        case 'r':
            read_file = optarg;
            break;

        case 'R':
            read_group = optarg;
            break;

        case 'N':
            sample_name = optarg;
            break;

        case 'b':
            hts_file = optarg;
            break;

        case 'f':
            if (fastq1.empty()) fastq1 = optarg;
            else if (fastq2.empty()) fastq2 = optarg;
            else { cerr << "[vg map] error: more than two fastqs specified" << endl; exit(1); }
            break;

        case 'i':
            interleaved_fastq = true;
            break;

        case 'p':
            pair_window = atoi(optarg);
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'D':
            debug = true;
            break;

        case 'F':
            prefer_forward = true;
            break;

        case 'G':
            greedy_accept = true;
            break;

        case 'X':
            score_per_bp = atof(optarg);
            break;

        case 'J':
            output_json = true;
            break;

        case 'B':
            band_width = atoi(optarg);
            break;
 
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_map(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if (seq.empty() && read_file.empty() && hts_file.empty() && fastq1.empty()) {
        cerr << "error:[vg map] a sequence or read file is required when mapping" << endl;
        return 1;
    }

    string file_name;
    if (optind < argc) {
        file_name = argv[optind];
    }

    if (db_name.empty()) {
        if (file_name.empty()) {
            cerr << "error:[vg map] no graph or db given, exiting" << endl;
            return 1;
        } else {
            db_name = file_name + ".index";
        }
    }

    thread_count = get_thread_count();

    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    vector<vector<Alignment> > output_buffer;
    output_buffer.resize(thread_count);

    Index idx;
    idx.open_read_only(db_name);

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = new Mapper(&idx);
        m->best_clusters = best_clusters;
        m->hit_max = hit_max;
        m->debug = debug;
        if (score_per_bp) m->target_score_per_bp = score_per_bp;
        if (sens_step) m->kmer_sensitivity_step = sens_step;
        m->prefer_forward = prefer_forward;
        m->greedy_accept = greedy_accept;
        m->thread_extension = thread_ex;
        m->cluster_min = cluster_min;
        m->max_attempts = max_attempts;
        m->min_kmer_entropy = min_kmer_entropy;
        mapper[i] = m;
    }

    if (!seq.empty()) {
        int tid = omp_get_thread_num();
        Alignment alignment = mapper[tid]->align(seq, kmer_size, kmer_stride, band_width);
        if (!sample_name.empty()) alignment.set_sample_name(sample_name);
        if (!read_group.empty()) alignment.set_read_group(read_group);
        if (!seq_name.empty()) alignment.set_name(seq_name);
        if (output_json) {
            cout << pb2json(alignment) << endl;
        } else {
            function<Alignment(uint64_t)> lambda =
                [&alignment] (uint64_t n) {
                return alignment;
            };
            stream::write(cout, 1, lambda);
        }
    }

    if (!read_file.empty()) {
        ifstream in(read_file);
        bool more_data = true;
#pragma omp parallel shared(in, more_data)
        {
            string line;
            int tid = omp_get_thread_num();
            while (more_data) {
                line.clear();
#pragma omp critical (readq)
                {
                    more_data = std::getline(in,line);
                }
                if (!line.empty()) {
                    Alignment alignment = mapper[tid]->align(line, kmer_size, kmer_stride, band_width);
                    if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                    if (!read_group.empty()) alignment.set_read_group(read_group);
                    if (output_json) {
                        string json2 = pb2json(alignment);
#pragma omp critical (cout)
                        cout << json2 << "\n";
                    } else {
                        auto& output_buf = output_buffer[tid];
                        output_buf.push_back(alignment);
                        stream::write_buffered(cout, output_buf, 1000);
                    }
                }
            }
        }
    }

    if (!hts_file.empty()) {
        function<void(Alignment&)> lambda =
            [&mapper,
             &output_buffer,
             &output_json,
             &kmer_size,
             &kmer_stride,
             &band_width]
            (Alignment& alignment) {
            int tid = omp_get_thread_num();
            alignment = mapper[tid]->align(alignment, kmer_size, kmer_stride, band_width);
            if (output_json) {
                string json2 = pb2json(alignment);
#pragma omp critical (cout)
                cout << json2 << "\n";
            } else {
                auto& output_buf = output_buffer[tid];
                output_buf.push_back(alignment);
                stream::write_buffered(cout, output_buf, 1000);
            }
        };
        // run
        hts_for_each_parallel(hts_file, lambda);
    }

    if (!fastq1.empty()) {
        if (interleaved_fastq) {
            // paired interleaved
            function<void(Alignment&, Alignment&)> lambda =
                [&mapper,
                 &output_buffer,
                 &output_json,
                 &kmer_size,
                 &kmer_stride,
                 &band_width,
                 &pair_window]
                (Alignment& aln1, Alignment& aln2) {
                int tid = omp_get_thread_num();
                auto alnp = mapper[tid]->align_paired(aln1, aln2, kmer_size, kmer_stride, band_width, pair_window);
                if (output_json) {
                    string json1 = pb2json(alnp.first);
                    string json2 = pb2json(alnp.second);
#pragma omp critical (cout)
                    cout << json1 << "\n" << json2 << "\n";
                } else {
                    auto& output_buf = output_buffer[tid];
                    output_buf.push_back(alnp.first);
                    output_buf.push_back(alnp.second);
                    stream::write_buffered(cout, output_buf, 1000);
                }
            };
            fastq_paired_interleaved_for_each_parallel(fastq1, lambda);
        } else if (fastq2.empty()) {
            // single
            function<void(Alignment&)> lambda =
                [&mapper,
                 &output_buffer,
                 &output_json,
                 &kmer_size,
                 &kmer_stride,
                 &band_width]
                (Alignment& alignment) {
                int tid = omp_get_thread_num();
                alignment = mapper[tid]->align(alignment, kmer_size, kmer_stride, band_width);
                if (output_json) {
                    string json2 = pb2json(alignment);
#pragma omp critical (cout)
                    cout << json2 << "\n";
                } else {
                    auto& output_buf = output_buffer[tid];
                    output_buf.push_back(alignment);
                    stream::write_buffered(cout, output_buf, 1000);
                }
            };
            fastq_unpaired_for_each_parallel(fastq1, lambda);
        } else {
            // paired two-file
            function<void(Alignment&, Alignment&)> lambda =
                [&mapper,
                 &output_buffer,
                 &output_json,
                 &kmer_size,
                 &kmer_stride,
                 &band_width,
                 &pair_window]
                (Alignment& aln1, Alignment& aln2) {
                int tid = omp_get_thread_num();
                auto alnp = mapper[tid]->align_paired(aln1, aln2, kmer_size, kmer_stride, band_width, pair_window);
                if (output_json) {
                    string json1 = pb2json(alnp.first);
                    string json2 = pb2json(alnp.second);
#pragma omp critical (cout)
                    cout << json1 << "\n" << json2 << "\n";
                } else {
                    auto& output_buf = output_buffer[tid];
                    output_buf.push_back(alnp.first);
                    output_buf.push_back(alnp.second);
                    stream::write_buffered(cout, output_buf, 1000);
                }
            };
            fastq_paired_two_files_for_each_parallel(fastq1, fastq2, lambda);
        }
    }

    // clean up
    for (int i = 0; i < thread_count; ++i) {
        delete mapper[i];
        auto& output_buf = output_buffer[i];
        if (!output_json) {
            stream::write_buffered(cout, output_buf, 0);
        }
    }

    cout.flush();

    return 0;

}

void help_view(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] [ <graph.vg> | <graph.json> | <aln.gam> | <read1.fq> [<read2.fq>] ]" << endl
         << "options:" << endl
         << "    -g, --gfa            output GFA format (default)" << endl
         << "    -F, --gfa-in         input GFA format" << endl
         
         << "    -v, --vg             output VG format (input defaults to GFA)" << endl
         << "    -V, --vg-in          input VG format (default)" << endl
         
         << "    -j, --json           output JSON format" << endl
         << "    -J, --json-in        input JSON format" << endl
         
         << "    -G, --gam            output GAM format (vg alignment format: Graph " << endl
         << "                         Alignment/Map)" << endl
         << "    -a, --align-in       input GAM format" << endl
         << "    -A, --aln-graph GAM  add alignments from GAM to the graph" << endl
         
         << "    -d, --dot            output dot format" << endl
         
         << "    -b, --bam            input BAM or other htslib-parseable alignments" << endl
         
         << "    -f, --fastq          input fastq (output defaults to GAM). Takes two " << endl
         << "                         positional file arguments if paired" << endl
         << "    -i, --interleaved    fastq is interleaved paired-ended" << endl;
    //<< "    -p, --paths           extract paths from graph in VG format" << endl;
    // TODO: Can we regularize the option names for input and output types?
}

int main_view(int argc, char** argv) {

    if (argc == 2) {
        help_view(argv);
        return 1;
    }

    // Supported conversions:
    //      TO  vg  json    gfa gam bam fastq   dot
    // FROM
    // vg       Y   Y       Y   N   N   N       Y       
    // json     Y   Y       Y   N   N   N       Y
    // gfa      Y   Y       Y   N   N   N       Y  
    // gam      N   Y       N   N   N   N       N
    // bam      N   N       N   Y   N   N       N
    // fastq    N   N       N   Y   N   N       N
    // dot      N   N       N   N   N   N       N

    string output_type;
    string input_type;
    string alignments;
    string fastq1, fastq2;
    bool interleaved_fastq = false;

    int c;
    optind = 2; // force optind past "view" argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"dot", no_argument, 0, 'd'},
                {"gfa", no_argument, 0, 'g'},
                {"gfa-in", no_argument, 0, 'F'},
                {"json",  no_argument, 0, 'j'},
                {"json-in",  no_argument, 0, 'J'},
                {"vg", no_argument, 0, 'v'},
                {"vg-in", no_argument, 0, 'V'},
                {"paths", no_argument, 0, 'p'},
                {"align-in", no_argument, 0, 'a'},
                {"gam", no_argument, 0, 'G'},
                {"bam", no_argument, 0, 'b'},
                {"fastq", no_argument, 0, 'f'},
                {"interleaved", no_argument, 0, 'i'},
                {"aln-graph", required_argument, 0, 'A'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgFjJhvVpaGbifA:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            output_type = "dot";
            break;
 
        case 'g':
            output_type = "gfa";
            break;
            
        case 'F':
            input_type = "gfa";
            break;
 
        case 'j':
            output_type = "json";
            break;
            
        case 'J':
            input_type = "json";
            break;

        case 'v':
            output_type = "vg";
            if(input_type.empty()) {
                // Default to GFA -> VG
                input_type = "gfa";
            }
            break;
            
        case 'V':
            input_type = "vg";
            break;

        case 'G':
            output_type = "gam";
            break;
            
        case 'a':
            input_type = "gam";
            if(output_type.empty()) {
                // Default to GAM -> JSON
                output_type = "json";
            }
            break;

        case 'b':
            input_type = "bam";
            if(output_type.empty()) {
                // Default to BAM -> GAM, since BAM isn't convertable to our normal default.
                output_type = "gam";
            }
            break;

        case 'p':
            input_type = "paths";
            output_type = "paths";
            break;

        case 'f':
            input_type = "fastq";
            if(output_type.empty()) {
                // Default to FASTQ -> GAM
                output_type = "gam";
            }
            break;

        case 'i':
            interleaved_fastq = true;
            break;

        case 'A':
            alignments = optarg;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_view(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    // If the user specified nothing else, we default to VG in and GFA out.
    if (input_type.empty()) {
        input_type = "vg";
    }
    if (output_type.empty()) {
        output_type = "gfa";
    }
    
    vector<Alignment> alns;
    if (!alignments.empty()) {
        function<void(Alignment&)> lambda = [&alns](Alignment& aln) { alns.push_back(aln); };
        ifstream in;
        in.open(alignments.c_str());
        stream::for_each(in, lambda);
    }

    VG* graph = nullptr;
    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    string file_name = argv[optind];
    if (input_type == "vg") {
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        // VG can convert to any of the graph formats, so keep going
    } else if (input_type == "gfa") {
        if (file_name == "-") {
            graph = new VG;
            graph->from_gfa(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG;
            graph->from_gfa(in);
        }
        // GFA can convert to any of the graph formats, so keep going
    } else if(input_type == "json") {
        // We need to load a JSON graph, which means we need to mess about with FILE pointers, since jansson is a C API.
        
        FILE* json_file;
        
        if (file_name == "-") {
            // Read standard input
            json_file = stdin;
        } else {
            // Open the file for reading
            json_file = fopen(file_name.c_str(), "r");
        }
        
        // Make a new VG that calls this function over and over to read Graphs.
        function<bool(Graph&)> get_next_graph = [&](Graph& subgraph) -> bool {
            // Check if the file ends now, and skip whitespace between records.
            char peeked;
            do {
                peeked = fgetc(json_file);
                if(peeked == EOF) {
                    // File ended or otherwise errored. TODO: check for other
                    // errors and complain.
                    return false;
                }
            } while(isspace(peeked));
            // Put it back
            ungetc(peeked, json_file);
            
            // Now we know we have non-whitespace between here and EOF.
            // If it's not JSON, we want to die. So read it as JSON.
            json2pb(subgraph, json_file);
            
            // We read it successfully!
            return true;
        };
        graph = new VG(get_next_graph, false);
        
    } else if (input_type == "paths") {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG;
        graph->paths.load(in);
        // Paths in is always used with paths out.
    } else if (input_type == "gam") {
        if (output_type == "json") {
            // convert values to printable ones
            function<void(Alignment&)> lambda = [](Alignment& a) {
                //alignment_quality_short_to_char(a);
                cout << pb2json(a) << "\n";
            };
            if (file_name == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each(in, lambda);
            }
            
            cout.flush();
            return 0;
        } else {
            // todo
            cerr << "[vg view] error: GAM can only be converted to JSON" << endl;
            return 1;
        }
    } else if (input_type == "bam") {
        if (output_type == "gam") {
//function<void(const Alignment&)>& lambda) {
            // todo write buffering procedure in alignment.cpp
            vector<Alignment> buf;
            function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                buf.push_back(aln);
                if (buf.size() > 1000) {
                    write_alignments(std::cout, buf);
                    buf.clear();
                }
            };
            hts_for_each(file_name, lambda);
            write_alignments(std::cout, buf);
            buf.clear();
            cout.flush();
            return 0;
        } else if (output_type == "json") {
            // todo
            cerr << "[vg view] error: BAM to JSON conversion not yet implemented" << endl;
            return 0;
        } else {
            cerr << "[vg view] error: BAM can only be converted to GAM" << endl;
            return 1;
        }
    } else if (input_type == "fastq") {
        fastq1 = argv[optind++];
        if (optind < argc) {
            fastq2 = argv[optind];            
        }
        if (output_type == "gam") {
            vector<Alignment> buf;
            if (!interleaved_fastq && fastq2.empty()) {
                function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                    buf.push_back(aln);
                    if (buf.size() > 1000) {
                        write_alignments(std::cout, buf);
                        buf.clear();
                    }
                };
                fastq_unpaired_for_each(fastq1, lambda);
            } else if (interleaved_fastq && fastq2.empty()) {
                function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                    buf.push_back(aln1);
                    buf.push_back(aln2);
                    if (buf.size() > 1000) {
                        write_alignments(std::cout, buf);
                        buf.clear();
                    }
                };
                fastq_paired_interleaved_for_each(fastq1, lambda);
            } else if (!fastq2.empty()) {
                function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                    buf.push_back(aln1);
                    buf.push_back(aln2);
                    if (buf.size() > 1000) {
                        write_alignments(std::cout, buf);
                        buf.clear();
                    }
                };
                fastq_paired_two_files_for_each(fastq1, fastq2, lambda);
            }
            write_alignments(std::cout, buf);
            buf.clear();
        } else {
            // We can't convert fastq to the other graph formats
            cerr << "[vg view] error: FASTQ can only be converted to GAM" << endl;
            return 1;
        }
        cout.flush();
        return 0;
    }

    if(graph == nullptr) {
        // Make sure we didn't forget to implement an input format.
        cerr << "[vg view] error: cannot load graph in " << input_type << " format" << endl;
        return 1;
    }

    // Now we know graph was filled in from the input format. Spit it out in the
    // requested output format.

    if (output_type == "dot") {
        graph->to_dot(std::cout, alns);
    } else if (output_type == "json") {
        cout << pb2json(graph->graph) << endl;
    } else if (output_type == "gfa") {
        graph->to_gfa(std::cout);
    } else if (output_type == "vg") {
        graph->serialize_to_ostream(cout);
    } else if (output_type == "paths") {
        function<void(Path&)> dump_path = [](Path& p) {
            for (int i = 0; i < p.mapping_size(); ++i) {
                cout << p.name() << "\t" << p.mapping(i).position().node_id() << endl;
            }
        };
        graph->paths.for_each(dump_path);
    } else {
        // We somehow got here with a bad output format.
        cerr << "[vg view] error: cannot save a graph in " << output_type << " format" << endl;
        return 1;
    }
    
    cout.flush();
    delete graph;

    return 0;
}

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -P, --ref-paths FILE  write reference paths in protobuf/gzip format to FILE" << endl
         << "    -R, --region REGION   specify a particular chromosome" << endl
         << "    -z, --region-size N   variants per region to parallelize" << endl
         << "    -m, --node-max N      limit the maximum allowable node sequence size" << endl
         << "                          nodes greater than this threshold will be divided" << endl
         << "    -p, --progress        show progress" << endl
         << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl;
}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    string fasta_file_name, vcf_file_name, json_filename;
    string region;
    string output_type = "VG";
    bool progress = false;
    int vars_per_region = 25000;
    int max_node_size = 0;
    string ref_paths_file;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"ref-paths", required_argument, 0, 'P'},
                {"progress",  no_argument, 0, 'p'},
                {"region-size", required_argument, 0, 'z'},
                {"threads", required_argument, 0, 't'},
                {"region", required_argument, 0, 'R'},
                {"node-max", required_argument, 0, 'm'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:phz:t:R:m:P:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'v':
            vcf_file_name = optarg;
            break;

        case 'r':
            fasta_file_name = optarg;
            break;

        case 'P':
            ref_paths_file = optarg;
            break;

        case 'p':
            progress = true;
            break;
 
        case 'z':
            vars_per_region = atoi(optarg);
            break;

        case 'R':
            region = optarg;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'm':
            max_node_size = atoi(optarg);
            break;
            
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_construct(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    vcflib::VariantCallFile variant_file;
    if (!vcf_file_name.empty()) {
        variant_file.open(vcf_file_name);
        if (!variant_file.is_open()) {
            cerr << "error:[vg construct] could not open" << vcf_file_name << endl;
            return 1;
        }
    }

    FastaReference reference;
    if (fasta_file_name.empty()) {
        cerr << "error:[vg construct] a reference is required for graph construction" << endl;
        return 1;
    }
    reference.open(fasta_file_name);

    // store our reference sequence paths
    Paths ref_paths;

    VG graph(variant_file, reference, region, vars_per_region, max_node_size, progress);

    if (!ref_paths_file.empty()) {
        ofstream paths_out(ref_paths_file);
        graph.paths.write(paths_out);
    }

    graph.serialize_to_ostream(std::cout);

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

void vg_help(char** argv) {
    cerr << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl 
         << "  -- construct     graph construction" << endl
         << "  -- view          format conversions for graphs and alignments" << endl
         << "  -- index         index features of the graph in a disk-backed key/value store" << endl
         << "  -- find          use an index to find nodes, edges, kmers, or positions" << endl
         << "  -- paths         traverse paths in the graph" << endl
         << "  -- align         local alignment" << endl
         << "  -- map           global alignment" << endl
         << "  -- stats         metrics describing graph properties" << endl
         << "  -- join          combine graphs via a new head" << endl
         << "  -- ids           manipulate node ids" << endl
         << "  -- concat        concatenate graphs tail-to-head" << endl
         << "  -- kmers         enumerate kmers of the graph" << endl
         << "  -- sim           simulate reads from the graph" << endl
         << "  -- mod           filter, transform, and edit the graph" << endl
         << "  -- surject       map alignments onto specific paths" << endl;
}

int main(int argc, char *argv[])
{

    if (argc == 1) {
        vg_help(argv);
        return 1;
    }

    //omp_set_dynamic(1); // use dynamic scheduling

    string command = argv[1];
    if (command == "construct") {
        return main_construct(argc, argv);
    } else if (command == "view") {
        return main_view(argc, argv);
    } else if (command == "align") {
        return main_align(argc, argv);
    } else if (command == "map") {
        return main_map(argc, argv);
    } else if (command == "index") {
        return main_index(argc, argv);
    } else if (command == "find") {
        return main_find(argc, argv);
    } else if (command == "paths") {
        return main_paths(argc, argv);
    } else if (command == "stats") {
        return main_stats(argc, argv);
    } else if (command == "join") {
        return main_join(argc, argv);
    } else if (command == "ids") {
        return main_ids(argc, argv);
    } else if (command == "concat") {
        return main_concat(argc, argv);
    } else if (command == "kmers") {
        return main_kmers(argc, argv);
    } else if (command == "sim") {
        return main_sim(argc, argv);
    } else if (command == "mod") {
        return main_mod(argc, argv);
    } else if (command == "surject") {
        return main_surject(argc, argv);
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
