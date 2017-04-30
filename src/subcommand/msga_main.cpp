#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include "../stream.hpp"

using namespace vg;
using namespace vg::subcommand;


void help_msga(char** argv) {
    cerr << "usage: " << argv[0] << " msga [options] >graph.vg" << endl
         << "Multiple sequence / graph aligner." << endl
         << endl
         << "options:" << endl
         << "inputs:" << endl
         << "    -f, --from FILE         use sequences in (fasta) FILE" << endl
         << "    -n, --name NAME         include this sequence" << endl
         << "                             (If any --name is specified, use only" << endl
         << "                              specified sequences from FASTA files.)" << endl
         << "    -b, --base NAME         use this sequence as the graph basis if graph is empty" << endl
         << "    -s, --seq SEQUENCE      literally include this sequence" << endl
         << "    -g, --graph FILE        include this graph" << endl
         << "alignment:" << endl
         << "    -k, --min-seed INT      minimum seed (MEM) length [estimated given -e]" << endl
         << "    -c, --hit-max N         ignore kmers or MEMs who have >N hits in our index [512]" << endl
         << "    -e, --seed-chance FLOAT set {-k} such that this fraction of {-k} length hits will by by chance [0.05]" << endl
         << "    -Y, --max-seed INT      ignore seeds longer than this length [0]" << endl
         << "    -r, --reseed-x FLOAT    look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]" << endl
         << "    -u, --try-up-to INT     attempt to align up to the INT best candidate chains of seeds [64]" << endl
         << "    -l, --try-at-least INT  attempt to align up to the INT best candidate chains of seeds [4]" << endl
         << "    -W, --min-chain INT     discard a chain if seeded bases shorter than INT [0]" << endl
         << "    -C, --drop-chain FLOAT  drop chains shorter than FLOAT fraction of the longest overlapping chain [0.4]" << endl
         << "    -P, --min-ident FLOAT   accept alignment only if the alignment identity is >= FLOAT [0]" << endl
         << "    -F, --min-band-mq INT   require mapping quality for each band to be at least this [0]" << endl
         << "    -H, --max-target-x N    skip cluster subgraphs with length > N*read_length [100]" << endl
         << "    -w, --band-width INT    band width for long read alignment [256]" << endl
         << "    -M, --max-multimaps INT when set > 1, thread an optimal alignment through the multimappings of each band [1]" << endl
         << "local alignment parameters:" << endl
         << "    -q, --match INT         use this match score [1]" << endl
         << "    -z, --mismatch INT      use this mismatch penalty [4]" << endl
         << "    -o, --gap-open INT      use this gap open penalty [6]" << endl
         << "    -y, --gap-extend INT    use this gap extension penalty [1]" << endl
         << "    -L, --full-l-bonus INT  the full-length alignment bonus [5]" << endl
         << "index generation:" << endl
         << "    -K, --idx-kmer-size N   use kmers of this size for building the GCSA indexes (default: 16)" << endl
         << "    -O, --idx-no-recomb     index only embedded paths, not recombinations of them" << endl
         << "    -E, --idx-edge-max N    reduce complexity of graph indexed by GCSA using this edge max (default: off)" << endl
         << "    -Q, --idx-prune-subs N  prune subgraphs shorter than this length from input graph to GCSA (default: off)" << endl
         << "    -m, --node-max N        chop nodes to be shorter than this length (default: 2* --idx-kmer-size)" << endl
         << "    -X, --idx-doublings N   use this many doublings when building the GCSA indexes (default: 2)" << endl
         << "graph normalization:" << endl
         << "    -N, --normalize         normalize the graph after assembly" << endl
         << "    -Z, --circularize       the input sequences are from circular genomes, circularize them after inclusion" << endl
         << "generic parameters:" << endl
         << "    -D, --debug             print debugging information about construction to stderr" << endl
         << "    -A, --debug-align       print debugging information about alignment to stderr" << endl
         << "    -S, --align-progress    show a progress bar for each banded alignment" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << endl
         << "Construct a multiple sequence alignment from all sequences in the" << endl
         << "input fasta-format files, graphs, and sequences. Uses the MEM mapping algorithm." << endl
         << endl
         << "Emits the resulting MSA as a (vg-format) graph." << endl;
}

int main_msga(int argc, char** argv) {

    if (argc == 2) {
        help_msga(argv);
        return 1;
    }

    vector<string> fasta_files;
    set<string> seq_names;
    vector<string> sequences;
    vector<string> graph_files;
    string base_seq_name;
    int idx_kmer_size = 16;
    int idx_doublings = 2;
    int hit_max = 100;
    int max_attempts = 10;
    // if we set this above 1, we use a dynamic programming process to determine the
    // optimal alignment through a series of bands based on a proximity metric
    int max_multimaps = 1;
    float min_identity = 0.0;
    int band_width = 256;
    size_t doubling_steps = 3;
    bool debug = false;
    bool debug_align = false;
    size_t node_max = 0;
    int alignment_threads = get_thread_count();
    int edge_max = 0;
    int subgraph_prune = 0;
    bool normalize = false;
    int iter_max = 1;
    bool mem_chaining = true;
    int max_mem_length = 0;
    int min_mem_length = 0;
    int max_target_factor = 100;
    bool idx_path_only = false;
    int match = 1;
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 1;
    int full_length_bonus = 5;
    bool circularize = false;
    float chance_match = 0.05;
    int mem_reseed_length = -1;
    int min_cluster_length = 0;
    float mem_reseed_factor = 1.5;
    float random_match_chance = 0.05;
    int extra_multimaps = 64;
    int min_multimaps = 4;
    float drop_chain = 0.4;
    int max_mapping_quality = 60;
    int method_code = 1;
    int maybe_mq_threshold = 60;
    int min_banded_mq = 0;
    bool use_fast_reseed = true;
    bool show_align_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                {"help", no_argument, 0, 'h'},
                {"from", required_argument, 0, 'f'},
                {"name", required_argument, 0, 'n'},
                {"seq", required_argument, 0, 's'},
                {"graph", required_argument, 0, 'g'},
                {"base", required_argument, 0, 'b'},
                {"idx-kmer-size", required_argument, 0, 'K'},
                {"idx-no-recomb", no_argument, 0, 'O'},
                {"idx-doublings", required_argument, 0, 'X'},
                {"band-width", required_argument, 0, 'w'},
                {"debug", no_argument, 0, 'D'},
                {"debug-align", no_argument, 0, 'A'},
                {"context-depth", required_argument, 0, 'c'},
                {"min-ident", required_argument, 0, 'P'},
                {"min-banded-mq", required_argument, 0, 'F'},
                {"idx-edge-max", required_argument, 0, 'E'},
                {"idx-prune-subs", required_argument, 0, 'Q'},
                {"normalize", no_argument, 0, 'N'},
                {"min-seed", required_argument, 0, 'k'},
                {"max-seed", required_argument, 0, 'Y'},
                {"hit-max", required_argument, 0, 'c'},
                {"seed-chance", required_argument, 0, 'e'},
                {"reseed-x", required_argument, 0, 'r'},
                {"threads", required_argument, 0, 't'},
                {"node-max", required_argument, 0, 'm'},
                {"max-target-x", required_argument, 0, 'H'},
                {"max-multimaps", required_argument, 0, 'M'},
                {"match", required_argument, 0, 'q'},
                {"mismatch", required_argument, 0, 'z'},
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'y'},
                {"full-l-bonus", required_argument, 0, 'L'},
                {"circularize", no_argument, 0, 'Z'},
                {"min-chain", required_argument, 0, 'W'},
                {"try-up-to", required_argument, 0, 'u'},
                {"try-at-least", required_argument, 0, 'l'},
                {"drop-chain", required_argument, 0, 'C'},
                {"align-progress", no_argument, 0, 'S'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hf:n:s:g:b:K:X:w:DAc:P:E:Q:NY:H:t:m:M:q:OI:i:o:y:ZW:z:k:L:e:r:u:l:C:F:S",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'k':
            min_mem_length = atoi(optarg);
            break;

        case 'r':
            mem_reseed_factor = atof(optarg);
            break;

        case 'W':
            min_cluster_length = atoi(optarg);
            break;

        case 'e':
            chance_match = atof(optarg);
            break;

        case 'Y':
            max_mem_length = atoi(optarg);
            break;

        case 'c':
            hit_max = atoi(optarg);
            break;

        case 'M':
            max_multimaps = atoi(optarg);
            break;

        case 'l':
            min_multimaps = atoi(optarg);
            break;

        case 'u':
            extra_multimaps = atoi(optarg);
            break;
            
        case 'H':
            max_target_factor = atoi(optarg);
            break;

        case 'f':
            fasta_files.push_back(optarg);
            break;

        case 'n':
            seq_names.insert(optarg);
            break;

        case 's':
            sequences.push_back(optarg);
            break;

        case 'b':
            base_seq_name = optarg;
            break;

        case 'g':
            if (graph_files.size() != 0) {
                cerr << "[vg msga] Error: graph-graph alignment is not yet implemented." << endl
                     << "We can only use one input graph." << endl;
                return 1;
            }
            graph_files.push_back(optarg);
            break;

        case 'w':
            band_width = atoi(optarg);
            break;

        case 'D':
            debug = true;
            break;

        case 'A':
            debug_align = true;
            break;

        case 'S':
            show_align_progress = true;
            break;

        case 'X':
            doubling_steps = atoi(optarg);
            break;

        case 'K':
            idx_kmer_size = atoi(optarg);
            break;


        case 'O':
            idx_path_only = true;
            break;

        case 'm':
            node_max = atoi(optarg);
            break;

        case 'N':
            normalize = true;
            break;

        case 'Z':
            circularize = true;
            break;

        case 'C':
            drop_chain = atof(optarg);
            break;

        case 'P':
            min_identity = atof(optarg);
            break;

        case 'F':
            min_banded_mq = atoi(optarg);
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            alignment_threads = atoi(optarg);
            break;

        case 'Q':
            subgraph_prune = atoi(optarg);
            break;

        case 'E':
            edge_max = atoi(optarg);
            break;

        case 'q':
            match = atoi(optarg);
            break;

        case 'z':
            mismatch = atoi(optarg);
            break;

        case 'o':
            gap_open = atoi(optarg);
            break;

        case 'y':
            gap_extend = atoi(optarg);
            break;

        case 'L':
            full_length_bonus = atoi(optarg);
            break;

        case 'h':
        case '?':
            help_msga(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    // mapping quality method boilerplate
    MappingQualityMethod mapping_quality_method;
    if (method_code == 0) {
        mapping_quality_method = None;
    }
    else if (method_code == 1) {
        mapping_quality_method = Approx;
    }
    else if (method_code == 2) {
        mapping_quality_method = Exact;
    }
    else {
        cerr << "error:[vg map] unrecognized mapping quality method command line arg '" << method_code << "'" << endl;
        return 1;
    }


    // build the graph or read it in from input
    VG* graph;
    if (graph_files.size() == 1) {
        string file_name = graph_files.front();
        get_input_file(file_name, [&](istream& in) {
            graph = new VG(in);
        });
    } else {
        graph = new VG;
    }

    // we should chop up the inputs into bits
    // then run a kind of alignment/overlap assembly on them
    // to generate the new graph/msa
    // TODO refactor into class

    // map from name to sequence, just a transformation of FASTA records
    map<string, string> strings;

    // open the fasta files, read in the sequences
    vector<string> names_in_order;

    for (auto& fasta_file_name : fasta_files) {
        FastaReference ref;
        ref.open(fasta_file_name);
        if (debug) cerr << "loading " << fasta_file_name << endl;
        for (auto& name : ref.index->sequenceNames) {
            if (!seq_names.empty() && seq_names.count(name) == 0) continue;
            // only use the sequence if we have whitelisted it
            // and also sanitize the input so we have only ATGCN
            strings[name] = nonATGCNtoN(ref.getSequence(name));
            names_in_order.push_back(name);
        }
    }

    // give a label to sequences passed on the command line
    // use the sha1sum, take the head
    // collision avoidance with nonce ensures we get the same names for the same sequences across multiple runs
    for (auto& s : sequences) {
        auto name = sha1head(s, 8);
        int nonce = 0;
        while (strings.find(name) != strings.end()) {
            stringstream ss;
            ss << s << ++nonce;
            name = sha1head(ss.str(), 8);
        }
        strings[name] = nonATGCNtoN(s);
        names_in_order.push_back(name);
    }

    // always go from biggest to smallest in progression
    sort(names_in_order.begin(), names_in_order.end(),
         [&strings](const string& s1, const string& s2) {
             return strings[s1].size() > strings[s2].size(); });

    // align, include, repeat

    if (debug) cerr << "preparing initial graph" << endl;

    size_t max_query_size = pow(2, doubling_steps) * idx_kmer_size;
    // limit max node size
    if (!node_max) node_max = 2*idx_kmer_size;

    // if our graph is empty, we need to take the first sequence and build a graph from it
    if (graph->empty()) {
        auto build_graph = [&graph,&node_max](const string& seq, const string& name) {
            graph->create_node(seq);
            graph->dice_nodes(node_max);
            graph->sort();
            graph->compact_ids();
            // the graph will have a single embedded path in it
            Path& path = *graph->graph.add_path();
            path.set_name(name);
            graph->for_each_node([&path](Node* node) {
                    auto mapping = path.add_mapping();
                    mapping->mutable_position()->set_node_id(node->id());
                    auto edit = mapping->add_edit();
                    edit->set_from_length(node->sequence().size());
                    edit->set_to_length(node->sequence().size());
                });
            graph->paths.extend(path);
            graph->sync_paths(); // sync the paths
        };
        // what's the first sequence?
        if (base_seq_name.empty()) {
            base_seq_name = names_in_order.front();
        }
        // we specified one we wanted to use as the first
        build_graph(strings[base_seq_name], base_seq_name);
    }

#ifdef debug
    cerr << "path going in " << pb2json(graph->graph.path(0)) << endl;
#endif

    // questions:
    // should we preferentially use sequences from fasta files in the order they were given?
    // (considering this a todo)
    // reverse complement?
    Mapper* mapper = nullptr;
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    xg::XG* xgidx = nullptr;
    size_t iter = 0;
    
    // Configure GCSA temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());

    auto rebuild = [&](VG* graph) {
        if (mapper) delete mapper;
        if (xgidx) delete xgidx;
        if (gcsaidx) delete gcsaidx;
        if (lcpidx) delete lcpidx;
    
        //stringstream s; s << iter++ << ".vg";
        graph->sort();
        graph->sync_paths();
        graph->graph.clear_path();
        graph->paths.to_graph(graph->graph);
        graph->rebuild_indexes();

        if (debug) cerr << "building xg index" << endl;
        xgidx = new xg::XG(graph->graph);

        if (debug) cerr << "building GCSA2 index" << endl;
        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if(!debug) gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        
        // Configure its temp directory to the system temp directory
        gcsa::TempFile::setDirectory(find_temp_dir());

        if (edge_max) {
            VG gcsa_graph = *graph; // copy the graph
            // remove complex components
            gcsa_graph.prune_complex_with_head_tail(idx_kmer_size, edge_max);
            if (subgraph_prune) gcsa_graph.prune_short_subgraphs(subgraph_prune);
            // then index
            gcsa_graph.build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
        } else {
            // if no complexity reduction is requested, just build the index
            graph->build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
        }
        mapper = new Mapper(xgidx, gcsaidx, lcpidx);
        { // set mapper variables
            mapper->hit_max = hit_max;
            mapper->max_multimaps = max_multimaps;
            mapper->min_multimaps = min_multimaps;
            mapper->maybe_mq_threshold = maybe_mq_threshold;
            mapper->debug = debug_align;
            mapper->min_identity = min_identity;
            mapper->min_banded_mq = min_banded_mq;
            mapper->drop_chain = drop_chain;
            mapper->min_mem_length = (min_mem_length > 0 ? min_mem_length
                                 : mapper->random_match_length(chance_match));
            mapper->mem_reseed_length = round(mem_reseed_factor * mapper->min_mem_length);
            mapper->min_cluster_length = min_cluster_length;
            if (debug) {
                cerr << "[vg msga] : min_mem_length = " << mapper->min_mem_length
                     << ", mem_reseed_length = " << mapper->mem_reseed_length
                     << ", min_cluster_length = " << mapper->min_cluster_length << endl;
            }
            mapper->fast_reseed = use_fast_reseed;
            mapper->mem_chaining = mem_chaining;
            mapper->max_target_factor = max_target_factor;
            mapper->set_alignment_scores(match, mismatch, gap_open, gap_extend);
            //mapper->adjust_alignments_for_base_quality = qual_adjust_alignments;
            mapper->extra_multimaps = extra_multimaps;
            mapper->mapping_quality_method = mapping_quality_method;
            mapper->full_length_alignment_bonus = full_length_bonus;
            mapper->max_mapping_quality = max_mapping_quality;
            // set up the multi-threaded alignment interface
            mapper->set_alignment_threads(alignment_threads);
            // Set scores after setting threads, since thread count changing resets scores
            mapper->set_alignment_scores(match, mismatch, gap_open, gap_extend);
            mapper->show_progress = show_align_progress;
        }
    };

    // set up the graph for mapping
    rebuild(graph);

    // todo restructure so that we are trying to map everything
    // add alignment score/bp bounds to catch when we get a good alignment
    for (auto& name : names_in_order) {
        //cerr << "do... " << name << " ?" << endl;
        if (!base_seq_name.empty() && name == base_seq_name) continue; // already embedded
        bool incomplete = true; // complete when we've fully included the sequence set
        int iter = 0;
        auto& seq = strings[name];
        //cerr << "doing... " << name << endl;
#ifdef debug
        {
            graph->serialize_to_file("msga-pre-" + name + ".vg");
            ofstream db_out("msga-pre-" + name + ".xg");
            xgidx->serialize(db_out);
            db_out.close();
        }
#endif
        while (incomplete && iter++ < iter_max) {
            stringstream s; s << iter; string iterstr = s.str();
            if (debug) cerr << name << ": adding to graph" << iter << endl;
            vector<Path> paths;
            vector<Alignment> alns;
            int j = 0;
            // align to the graph
            if (debug) cerr << name << ": aligning sequence of " << seq.size() << "bp against " <<
                graph->node_count() << " nodes" << endl;
            Alignment aln = simplify(mapper->align(seq, 0, 0, 0, band_width));
            if (aln.path().mapping_size()) {
                auto aln_seq = graph->path_string(aln.path());
                if (aln_seq != seq) {
                    cerr << "[vg msga] alignment corrupted, failed to obtain correct banded alignment (alignment seq != input seq)" << endl;
                    cerr << "expected " << seq << endl;
                    cerr << "got      " << aln_seq << endl;
                    ofstream f(name + "-failed-alignment-" + convert(j) + ".gam");
                    stream::write(f, 1, (std::function<Alignment(uint64_t)>)([&aln](uint64_t n) { return aln; }));
                    f.close();
                    graph->serialize_to_file(name + "-corrupted-alignment.vg");
                    exit(1);
                }
            } else {
                Edit* edit = aln.mutable_path()->add_mapping()->add_edit();
                edit->set_sequence(aln.sequence());
                edit->set_to_length(aln.sequence().size());
            }
            alns.push_back(aln);
            //if (debug) cerr << pb2json(aln) << endl; // huge in some cases
            paths.push_back(aln.path());
            paths.back().set_name(name); // cache name to trigger inclusion of path elements in graph by edit

            /*
               ofstream f(name + "-pre-edit-" + convert(j) + ".gam");
               stream::write(f, 1, (std::function<Alignment(uint64_t)>)([&aln](uint64_t n) { return aln; }));
               f.close();
               */

            ++j;

            // now take the alignment and modify the graph with it
            if (debug) cerr << name << ": editing graph" << endl;
            //graph->serialize_to_file(name + "-pre-edit.vg");
            graph->edit(paths);
            //if (!graph->is_valid()) cerr << "invalid after edit" << endl;
            //graph->serialize_to_file(name + "-immed-post-edit.vg");
            graph->dice_nodes(node_max);
            //if (!graph->is_valid()) cerr << "invalid after dice" << endl;
            //graph->serialize_to_file(name + "-post-dice.vg");
            if (debug) cerr << name << ": sorting and compacting ids" << endl;
            graph->sort();
            //if (!graph->is_valid()) cerr << "invalid after sort" << endl;
            graph->compact_ids(); // xg can't work unless IDs are compacted.
            //if (!graph->is_valid()) cerr << "invalid after compact" << endl;
            if (circularize) {
                if (debug) cerr << name << ": circularizing" << endl;
                graph->circularize({name});
                //graph->serialize_to_file(name + "-post-circularize.vg");
            }

            // the edit needs to cut nodes at mapping starts and ends
            // thus allowing paths to be included that map directly to entire nodes
            // XXX

            //graph->serialize_to_file(name + "-pre-index.vg");
            // update the paths
            graph->graph.clear_path();
            graph->paths.to_graph(graph->graph);
            // and rebuild the indexes
            rebuild(graph);
            //graph->serialize_to_file(name + "-post-index.vg");

            // verfy validity of path
            bool is_valid = graph->is_valid();
            auto path_seq = graph->path_string(graph->paths.path(name));
            incomplete = !(path_seq == seq) || !is_valid;
            if (incomplete) {
                cerr << "[vg msga] failed to include alignment, retrying " << endl
                    << "expected " << seq << endl
                    << "got " << path_seq << endl
                    << pb2json(aln.path()) << endl
                    << pb2json(graph->paths.path(name)) << endl;
                graph->serialize_to_file(name + "-post-edit.vg");
            }
        }
        // if (debug && !graph->is_valid()) cerr << "graph is invalid" << endl;
        if (incomplete && iter >= iter_max) {
            cerr << "[vg msga] Error: failed to include path " << name << endl;
            exit(1);
        }
    }

    // auto include_paths = [&mapper,
    //      kmer_size,
    //      kmer_stride,
    //      band_width,
    //      debug,
    //      &strings](VG* graph) {
    //          // include the paths in the graph
    //          if (debug) cerr << "including paths" << endl;
    //          for (auto& group : strings) {
    //              auto& name = group.first;
    //              if (debug) cerr << name << ": tracing path through graph" << endl;
    //              auto& seq = group.second;
    //              if (debug) cerr << name << ": aligning sequence of " << seq.size() << "bp" << endl;
    //              Alignment aln = mapper->align(seq, kmer_size, kmer_stride, band_width);
    //              //if (debug) cerr << "alignment score: " << aln.score() << endl;
    //              aln.mutable_path()->set_name(name);
    //              //if (debug) cerr << "alignment: " << pb2json(aln) << endl;
    //              // todo simplify in the mapper itself when merging the banded bits
    //              if (debug) cerr << name << ": labeling" << endl;
    //              graph->include(aln.path());
    //              // now repeat back the path
    //          }
    //      };

    if (normalize) {
        if (debug) cerr << "normalizing graph" << endl;
        graph->remove_non_path();
        graph->normalize();
        graph->dice_nodes(node_max);
        graph->sort();
        graph->compact_ids();
        if (!graph->is_valid()) {
            cerr << "[vg msga] warning! graph is not valid after normalization" << endl;
        }
    }

    // finally, validate the included paths
    set<string> failures;
    for (auto& sp : strings) {
        auto& name = sp.first;
        auto& seq = sp.second;
        if (seq != graph->path_string(graph->paths.path(name))) {
            /*
               cerr << "failed inclusion" << endl
               << "expected " << graph->path_string(graph->paths.path(name)) << endl
               << "got      " << seq << endl;
               */
            failures.insert(name);
        }
    }

    if (!failures.empty()) {
        stringstream ss;
        ss << "vg-msga-failed-include_";
        for (auto& s : failures) {
            cerr << "[vg msga] Error: failed to include path " << s << endl;
            ss << s << "_";
        }
        ss << ".vg";
        graph->serialize_to_file(ss.str());
        exit(1);
    }

    // return the graph
    graph->serialize_to_ostream(std::cout);
    delete graph;

    // todo....
    //
    // strategy for graph/graph alignment
    // ---------------
    // multiple graphs can be aligned by converting them into collections of named sequences
    // e.g. using a strided sampling of a long kmer space
    // of sufficient length to align to a second graph
    //
    // the multi-graph alignment is a graph which contains both of the
    // with homologous sequences merged and the paths from the input graphs retained
    //
    // a progressive approach can be used, where we first attempt to construct a graph using a particular
    // sequence size
    // then, by labeling the first graph where it is shown to map to the new graph, we can retain only
    // the portion which was not included, then attempt to include the remaining fragments using
    // more compute-intensive parameters
    //
    // to limit path complexity, random walks should be used to sample the path space of the first graph
    // we can erode the graph we are aligning as regions of it become completely aligned,
    // so as to avoid over-sampling the graph
    // we already have functionality for this in `vg sim`
    //
    // this is an elaborate but easily-written and flexible approach to aligning large graphs
    // efficiently

    return 0;
}

static Subcommand vg_msga("msga", "multiple sequence graph alignment", main_msga);
