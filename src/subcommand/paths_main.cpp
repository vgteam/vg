/** \file paths_main.cpp
 *
 * Defines the "vg paths" subcommand, which reads paths in the graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include <gbwt/dynamic_gbwt.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options]" << endl
         << "options:" << endl
         << "  input:" << endl
         << "    -v, --vg FILE         use the graph in this vg FILE" << endl
         << "    -x, --xg FILE         use the graph in the XG index FILE" << endl
         << "    -g, --gbwt FILE       use the GBWT index in FILE" << endl
         << "  inspection:" << endl
         << "    -X, --extract-gam     return (as GAM alignments) the stored paths in the graph" << endl
         << "    -V, --extract-vg      return (as path-only .vg) the queried paths (requires -x -g and -q or -Q)" << endl
        //<< "    -X, --extract         return (as a chunked graph) the stored paths in the graph" << endl
         << "    -L, --list            return (as a list of names, one per line) the path names" << endl
         << "    -T, --threads         return the threads (requires GBWT)" << endl
         << "    -q, --threads-by STR  return the threads with the given prefix (requires GBWT)" << endl
         << "    -Q, --paths-by STR    return the paths with the given prefix" << endl;
    //<< "    -s, --as-seqs         write each path as a sequence" << endl;
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    bool as_seqs = false;
    bool extract_as_gam = false;
    bool extract_as_vg = false;
    bool list_paths = false;
    string xg_file;
    string vg_file;
    string gbwt_file;
    string thread_prefix;
    string path_prefix;
    bool extract_threads = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"vg", required_argument, 0, 'v'},
            {"xg", required_argument, 0, 'x'},
            {"gbwt", required_argument, 0, 'g'},
            {"extract-gam", no_argument, 0, 'X'},
            {"extract-vg", no_argument, 0, 'V'},
            {"list", no_argument, 0, 'L'},
            {"max-length", required_argument, 0, 'l'},
            {"as-seqs", no_argument, 0, 's'},
            {"threads-by", required_argument, 0, 'q'},
            {"paths-by", required_argument, 0, 'Q'},
            {"threads", no_argument, 0, 'T'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hs:LXv:x:g:q:Q:VT",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vg_file = optarg;
            break;

        case 'x':
            xg_file = optarg;
            break;

        case 'g':
            gbwt_file = optarg;
            break;
            
        case 'X':
            extract_as_gam = true;
            break;

        case 'V':
            extract_as_vg = true;
            break;

        case 'L':
            list_paths = true;
            break;

        case 's':
            as_seqs = true;
            break;

        case 'q':
            thread_prefix = optarg;
            break;

        case 'Q':
            path_prefix = optarg;
            break;

        case 'T':
            extract_threads = true;
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

    if (!vg_file.empty() && !xg_file.empty()) {
        cerr << "[vg paths] Error: both vg and xg index given" << endl;
        assert(false);
    }

    if (!vg_file.empty()) {
        VG graph;
        if (vg_file == "-") {
            graph.from_istream(std::cin);
        } else {
            ifstream in(vg_file.c_str());
            graph.from_istream(in);
        }
        if (list_paths) {
            graph.paths.for_each_name([&](const string& name) {
                    cout << name << endl;
                });
        } else if (extract_as_gam) {
            vector<Alignment> alns = graph.paths_as_alignments();
            write_alignments(cout, alns);
        } else if (extract_as_vg) {
            cerr << "[vg paths] Error: vg extraction is only defined for prefix queries against a XG/GBWT index pair" << endl;
            assert(false);
        } else if (!thread_prefix.empty()) {
            cerr << "[vg paths] Error: a thread prefix query requires a XG/GBWT index pair" << endl;
            assert(false);
        }
        
    } else if (!xg_file.empty()) {
        xg::XG xgidx;
        ifstream in(xg_file.c_str());
        xgidx.load(in);
        if (list_paths) {
            size_t max_path = xgidx.max_path_rank();
            for (size_t i = 1; i <= max_path; ++i) {
                cout << xgidx.path_name(i) << endl;
            }
        } else if (!thread_prefix.empty() || extract_threads) {
            if (gbwt_file.empty()) {
                cerr << "[vg paths] Error: thread extraction requires a GBWT" << endl;
                assert(false);
            } else if (extract_as_gam == extract_as_vg) {
                cerr << "[vg paths] Error: thread extraction requires -V or -X to specifiy output format" << endl;
                assert(false);
            }
            gbwt::GBWT index;
            sdsl::load_from_file(index, gbwt_file);
            vector<int64_t> thread_ids;
            if (extract_threads) {
                for (gbwt::size_type id = 1; id <= index.sequences()/2; id += 1) {
                    thread_ids.push_back(id);
                }
            } else if (!thread_prefix.empty()) {
                thread_ids = xgidx.threads_named_starting(thread_prefix);
            }
            for (auto& id : thread_ids) {
                //cerr << "thread_id " << id << endl;
                gbwt::vector_type sequence = index.extract(gbwt::Path::encode(id-1, false));
                Path path;
                path.set_name(xgidx.thread_name(id));
                for (auto node : sequence) {
                    Mapping* m = path.add_mapping();
                    Position* p = m->mutable_position();
                    p->set_node_id(gbwt::Node::id(node));
                    p->set_is_reverse(gbwt::Node::is_reverse(node));
                    Edit* e = m->add_edit();
                    size_t len = xgidx.node_length(p->node_id());
                    e->set_to_length(len);
                    e->set_from_length(len);
                }
                if (extract_as_gam) {
                    vector<Alignment> alns;
                    alns.emplace_back(xgidx.path_as_alignment(path));
                    write_alignments(cout, alns);
                } else if (extract_as_vg) {
                    Graph g;
                    *(g.add_path()) = path;
                    vector<Graph> gb = { g };
                    stream::write_buffered(cout, gb, 0);
                }
            }
        } else {
            if (extract_as_gam) {
                auto alns = xgidx.paths_as_alignments();
                write_alignments(cout, alns);
            } else if (!path_prefix.empty()) {
                vector<Path> got = xgidx.paths_by_prefix(path_prefix);
                if (extract_as_gam) {
                    vector<Alignment> alns;
                    for (auto& path : got) {
                        alns.emplace_back(xgidx.path_as_alignment(path));
                    }
                    write_alignments(cout, alns);
                } else if (extract_as_vg) {
                    for(auto& path : got) {
                        Graph g;
                        *(g.add_path()) = xgidx.path(path.name());
                        vector<Graph> gb = { g };
                        stream::write_buffered(cout, gb, 0);
                    }
                }
            }
        }
    } else {
        cerr << "[vg paths] Error: a xg or vg file is required" << endl;
        assert(false);
    }
    
    return 0;

}

// Register subcommand
static Subcommand vg_paths("paths", "traverse paths in the graph", main_paths);

