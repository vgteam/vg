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
#include "../stream.hpp"
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
         << "    -L, --list            return (as a list of names, one per line) the path (or thread) names" << endl
         << "    -E, --lengths         return list of path names (as with -L) but paired with their lengths" << endl
         << "    -T, --threads         operate on threads instead of paths (requires GBWT)" << endl
         << "    -q, --threads-by STR  operate on threads with the given prefix instead of paths (requires GBWT)" << endl
         << "    -Q, --paths-by STR    return the paths with the given prefix" << endl;
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    bool extract_as_gam = false;
    bool extract_as_vg = false;
    bool list_names = false;
    string xg_file;
    string vg_file;
    string gbwt_file;
    string thread_prefix;
    string path_prefix;
    bool extract_threads = false;
    bool list_lengths = false;

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
            {"threads-by", required_argument, 0, 'q'},
            {"paths-by", required_argument, 0, 'Q'},
            {"threads", no_argument, 0, 'T'},
            {"lengths", no_argument, 0, 'E'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hLXv:x:g:q:Q:VTE",
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
            list_names = true;
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

        case 'E':
            list_names = true;
            list_lengths = true;
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
        exit(1);
    }
    
    if (!thread_prefix.empty() && extract_threads) {
        cerr << "[vg paths] Error: cannot extract all threads (-T) and also prefixed threads (-q)" << endl;
        exit(1);
    }
    
    // Load whatever indexes we were given
    unique_ptr<VG> graph;
    if (!vg_file.empty()) {
        // We want a vg
        graph = unique_ptr<VG>(new VG());
        // Load the vg
        get_input_file(vg_file, [&](istream& in) {
            graph->from_istream(in);
        });
    }
    unique_ptr<xg::XG> xg_index;
    if (!xg_file.empty()) {
        // We want an xg
        xg_index = unique_ptr<xg::XG>(new xg::XG());
        // Load the xg
        get_input_file(xg_file, [&](istream& in) {
            xg_index->load(in);
        });
    }
    unique_ptr<gbwt::GBWT> gbwt_index;
    if (!gbwt_file.empty()) {
        // We want a gbwt
        gbwt_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
        // Load the gbwt (TODO: support streams)
        sdsl::load_from_file(*gbwt_index, gbwt_file);
    }
    
    
    if (!thread_prefix.empty() || extract_threads) {
        // We are looking for threads, so we need the GBWT and the xg (which holds the thread name metadata)
        
        if (xg_index.get() == nullptr) {
            cerr << "[vg paths] Error: thread extraction requires an XG for thread metadata" << endl;
            exit(1);
        }
        if (gbwt_index.get() == nullptr) {
            cerr << "[vg paths] Error: thread extraction requires a GBWT" << endl;
            exit(1);
        }
        if (extract_as_gam == extract_as_vg && extract_as_vg == list_names) {
            cerr << "[vg paths] Error: thread extraction requires -V, -X, or -L to specifiy output format" << endl;
            exit(1);
        }
        vector<int64_t> thread_ids;
        if (extract_threads) {
            for (gbwt::size_type id = 1; id <= gbwt_index->sequences()/2; id += 1) {
                thread_ids.push_back(id);
            }
        } else if (!thread_prefix.empty()) {
            thread_ids = xg_index->threads_named_starting(thread_prefix);
        }
        
        // We may need to emit a stream of Alignemnts
        unique_ptr<stream::ProtobufEmitter<Alignment>> gam_emitter;
        // Or we might need to emit a stream of VG Graph objects
        unique_ptr<stream::ProtobufEmitter<Graph>> graph_emitter;
        if (extract_as_gam) {
            // Open up a GAM output stream
            gam_emitter = unique_ptr<stream::ProtobufEmitter<Alignment>>(new stream::ProtobufEmitter<Alignment>(cout));
        } else if (extract_as_vg) {
            // Open up a VG Graph chunk output stream
            graph_emitter = unique_ptr<stream::ProtobufEmitter<Graph>>(new stream::ProtobufEmitter<Graph>(cout));
        }
        
        for (auto& id : thread_ids) {
            // For each matching thread
            
            // Get its name
            auto thread_name = xg_index->thread_name(id); 
        
            if (list_names && !list_lengths) {
                // We are only interested in the name
                cout << thread_name << endl;
                continue;
            }
            
            // Otherwise we need the actual thread data
            gbwt::vector_type sequence = gbwt_index->extract(gbwt::Path::encode(id-1, false));
            Path path;
            path.set_name(thread_name);
            size_t rank = 1;
            for (auto node : sequence) {
                // Put each node in the constructed path
                Mapping* m = path.add_mapping();
                Position* p = m->mutable_position();
                p->set_node_id(gbwt::Node::id(node));
                p->set_is_reverse(gbwt::Node::is_reverse(node));
                Edit* e = m->add_edit();
                size_t len = xg_index->node_length(p->node_id());
                e->set_to_length(len);
                e->set_from_length(len);
                m->set_rank(rank++);
            }
            if (extract_as_gam) {
                // Write as an Alignment
                gam_emitter->write(xg_index->path_as_alignment(path));
            } else if (extract_as_vg) {
                // Write as a Path in a VG
                Graph g;
                *(g.add_path()) = path;
                graph_emitter->write(std::move(g));
            }
            if (list_lengths) {
                cout << path.name() << "\t" << path_to_length(path) << endl;
            }
        }
    } else if (graph.get() != nullptr) {
        // Handle non-thread queries from vg
        
        if (!path_prefix.empty()) {
            cerr << "[vg paths] Error: path prefix not supported for extracting from vg, only for extracting from xg" << endl;
            exit(1);
        }
        
        if (list_names) {
            graph->paths.for_each([&list_lengths](const Path& path) {
                cout << path.name();
                if (list_lengths) {
                    cout << "\t" << path_to_length(path);
                }
                cout << endl;
            });
        } else if (extract_as_gam) {
            vector<Alignment> alns = graph->paths_as_alignments();
            stream::ProtobufEmitter<Alignment> emitter(cout);
            for (auto& aln : alns) {
                emitter.write(std::move(aln));
            }
        } else if (extract_as_vg) {
            cerr << "[vg paths] Error: vg extraction is only defined for prefix queries against a XG/GBWT index pair" << endl;
            exit(1);
        } else {
            cerr << "[vg paths] Error: specify an operation to perform" << endl;
        }
    } else if (xg_index.get() != nullptr) {
        // Handle non-thread queries from xg
        if (list_names) {
            // We aren't looking for threads, but we are looking for names.
            size_t max_path = xg_index->max_path_rank();
            for (size_t i = 1; i <= max_path; ++i) {
                cout << xg_index->path_name(i);
                if (list_lengths) {
                    cout << "\t" << xg_index->path_length(i);
                }
                cout << endl;
            }
        } else if (extract_as_gam) {
            auto alns = xg_index->paths_as_alignments();
            stream::ProtobufEmitter<Alignment> emitter(cout);
            for (auto& aln : alns) {
                emitter.write(std::move(aln));
            }
        } else if (!path_prefix.empty()) {
            vector<Path> got = xg_index->paths_by_prefix(path_prefix);
            if (extract_as_gam) {
                stream::ProtobufEmitter<Alignment> emitter(cout);
                for (auto& path : got) {
                    emitter.write(xg_index->path_as_alignment(path));
                }
            } else if (extract_as_vg) {
                for(auto& path : got) {
                    Graph g;
                    *(g.add_path()) = xg_index->path(path.name());
                    vector<Graph> gb = { g };
                    stream::write_buffered(cout, gb, 0);
                }
            }
        } else {
            cerr << "[vg paths] Error: specify an operation to perform" << endl;
        }
    } else {
        cerr << "[vg paths] Error: an xg (-x) or vg (-v) file is required" << endl;
        exit(1);
    }
    
    return 0;

}

// Register subcommand
static Subcommand vg_paths("paths", "traverse paths in the graph", main_paths);

