// chunk_main.cpp: define the "vg chunk" subcommand, which chunks up vg graphs
// some similar logic to vg find, but focussed on generating multiple chunks
// at once. 

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <regex>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../stream.hpp"
#include "../utility.hpp"
#include "../chunker.hpp"
#include "../region.hpp"
#include "../haplotype_extracter.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;
using namespace xg;

void help_chunk(char** argv) {
    cerr << "usage: " << argv[0] << " chunk [options] > [chunk.vg]" << endl
         << "Splits a graph and/or alignment into smaller chunks" << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE       use this xg index to chunk subgraphs" << endl
         << "    -a, --gam-index FILE     chunk this gam index (made with vg index -N) instead of the graph" << endl
         << "    -g, --gam-and-graph      when used in combination with -a, both gam and graph will be chunked" << endl 
         << "path chunking:" << endl
         << "    -p, --path TARGET        write the chunk in the specified (1-based inclusive)\n"
         << "                             path range TARGET=path[:pos1[-pos2]]" << endl
         << "    -P, --path-list FILE     write chunks for all path regions in (line - separated file). format\n"
         << "                             for each as in -p (all paths chunked unless otherwise specified)" << endl
         << "    -e, --input-bed FILE     write chunks for all (0-based) bed regions" << endl
         << "id range chunking:" << endl
         << "    -r, --node-range N:M     write the chunk for the specified node range\n"
         << "    -R, --node-ranges FILE   write the chunk for each node range in (newline or whitespace separated) file\n"
         << "trace chunking:" << endl
         << "    -n, --start-node INT       start at this node" << endl
         << "    -d, --extend-distance INT  extend search this many nodes [default=50]" << endl
         << "    -T, --trace-file           whitespace-seprated file of node ids to trace chunks from" << endl
         << "general:" << endl
         << "    -s, --chunk-size N       create chunks spanning N bases (or nodes with -r/-R) for all input regions.\n"
         << "    -o, --overlap N          overlap between chunks when using -s [default = 0]" << endl        
         << "    -E, --output-bed FILE    write all created chunks to a bed file" << endl
         << "    -b, --prefix PATH        prefix output chunk paths with PATH. each chunk will have the following name:\n"
         << "                             <PATH>-<i>-<name>-<start>-<length>. <i> is the line# of the chunk in the\n"
         << "                             output bed. [default=./chunk]" << endl
         << "    -c, --context STEPS      expand the context of the chunk this many steps [0]" << endl       
         << "    -t, --threads N          for tasks that can be done in parallel, use this many threads [1]" << endl
         << "    -h, --help" << endl;
}

int main_chunk(int argc, char** argv) {

    if (argc == 2) {
        help_chunk(argv);
        return 1;
    }

    string xg_file;
    string gam_file;
    bool gam_and_graph = false;
    string region_string;
    string path_list_file;
    int chunk_size = 0;
    int overlap = 0;
    string in_bed_file;
    string out_bed_file;
    string out_chunk_prefix = "./chunk";
    int context_steps = 0;
    bool id_range = false;
    string node_range_string;
    string node_ranges_file;
    int threads = 1;
    vg::id_t trace_start = 0;
    int trace_extend = 50;
    string trace_file;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gam-name", required_argument, 0, 'a'},
            {"gam-and-graph", no_argument, 0, 'g'},
            {"path", required_argument, 0, 'p'},
            {"path-names", required_argument, 0, 'P'},
            {"chunk-size", required_argument, 0, 's'},
            {"overlap", required_argument, 0, 'o'},
            {"input-bed", required_argument, 0, 'e'},
            {"output-bed", required_argument, 0, 'E'},
            {"prefix", required_argument, 0, 'b'},
            {"context", required_argument, 0, 'c'},
            {"id-ranges", no_argument, 0, 'r'},
            {"id-range", no_argument, 0, 'R'},
            {"start-node", required_argument, 0, 'n'},
            {"extend-distance", required_argument, 0, 'd'},
            {"trace-file", required_argument, 0, 'T'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:a:gp:P:s:o:e:E:b:c:r:R:n:d:T:t:",
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_file = optarg;
            break;

        case 'a':
            gam_file = optarg;
            break;
            
        case 'g':
            gam_and_graph = true;
            break;            

        case 'p':
            region_string = optarg;
            break;

        case 'P':
            path_list_file = optarg;
            break;

        case 's':
            chunk_size = atoi(optarg);
            break;

        case 'o':
            overlap = atoi(optarg);
            break;

        case 'e':
            in_bed_file = optarg;
            break;

        case 'E':
            out_bed_file = optarg;
            break;

        case 'b':
            out_chunk_prefix = optarg;
            break;

        case'c':
            context_steps = atoi(optarg);
            break;

        case 'r':
            node_range_string = optarg;
            id_range = true;
            break;

        case 'R':
            node_ranges_file = optarg;
            id_range = true;
            break;

        case 'n':
            trace_start = atoi(optarg);
            break;

        case 'd':
            trace_extend = atoi(optarg);
            break;

        case 'T':
            trace_file = optarg;
            break;
            
        case 't':
            threads = atoi(optarg);
            break;

        case 'h':
        case '?':
            help_chunk(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    omp_set_num_threads(threads);            

    // check for invalid options combinations, of which there are many
    if (xg_file.empty()) {
        cerr << "error:[vg chunk] xg index (-x) required" << endl;
        return 1;
    }
    // need at most one of -p, -P, -r, -E, -e,  as an input
    if ((region_string.empty() ? 0 : 1) + (path_list_file.empty() ? 0 : 1) + (in_bed_file.empty() ? 0 : 1) +
        (node_ranges_file.empty() ? 0 : 1) + (node_range_string.empty() ? 0 : 1) +
        (trace_file.empty() ? 0 : 1)> 1 + (trace_start <= 0 ? 0 : 1)) {
        cerr << "error:[vg chunk] at most one of {-p, -P, -r, -E, -e, -n, -T} required to specify input regions" << endl;
        return 1;
    }

    // figure out which outputs we want.  the graph always
    // needs to be chunked, even if only gam output is requested,
    // because we use the graph to get the nodes we're looking for.
    // but we only write the subgraphs to disk if chunk_graph is true. 
    bool chunk_gam = !gam_file.empty();
    bool chunk_graph = gam_and_graph || !chunk_gam;

    bool trace_mode = trace_start > 0 || !trace_file.empty();

    // Load our index
    xg::XG xindex;
    ifstream in(xg_file.c_str());
    if (!in) {
        cerr << "error:[vg chunk] unable to load xg index file" << endl;
        return 1;
    }
    xindex.load(in);
    in.close();

    // This holds the RocksDB index that has all our reads, indexed by the nodes they visit.
    Index gam_index;
    if (chunk_gam) {
        gam_index.open_read_only(gam_file);
    }

    // Hold trace starting information
    vector<int64_t> trace_starts;
    
    // parse the regions into a list
    vector<Region> regions;
    if (!region_string.empty()) {
        Region region;
        parse_region(region_string, region);
        regions.push_back(region);
    }
    else if (!path_list_file.empty()) {
        ifstream pr_stream(path_list_file.c_str());
        if (!pr_stream) {
            cerr << "error:[vg chunk] unable to open path regions: " << path_list_file << endl;
            return 1;
        }
        while (pr_stream) {
            string buf;
            std::getline(pr_stream, buf);
            if (!buf.empty()) {
                Region region;
                parse_region(buf, region);
                regions.push_back(region);
            }
        }
    }
    else if (!in_bed_file.empty()) {
        parse_bed_regions(in_bed_file, regions);
    }
    else if (id_range) {
        istream* range_stream;
        if (!node_range_string.empty()) {
            range_stream = new stringstream(node_range_string);
        } else {
            range_stream = new ifstream(node_ranges_file);
            if (!(*range_stream)) {
                cerr << "error:[vg chunk] unable to open id ranges file: " << node_ranges_file << endl;
                return 1;
            }
        }
        do {
            string range;
            *range_stream >> range;
            if (!range.empty() && isdigit(range[0])) {
                Region region;
                vector<string> parts = split_delims(range, ":");
                if (parts.size() == 1) {
                    convert(parts.front(), region.start);
                    region.end = region.start;
                } else {
                    convert(parts.front(), region.start);
                    convert(parts.back(), region.end);
                }
                regions.push_back(region);
            }
        } while (*range_stream);
        delete range_stream;
    }
    else if (trace_mode) {
        if (trace_start > 0) {
            // one trace
            trace_starts.push_back(trace_start);
        } else {
            // file full of traces
            ifstream trace_stream(trace_file);
            if (!trace_stream) {
                cerr << "[vg chunk] unable to open traces file: " << trace_file << endl;
                return 1;
            }
            while (trace_stream) {
                int64_t ts;
                if (trace_stream >> ts) {
                    trace_starts.push_back(ts);
                }
            }
        }
    }
    else {
        // every path
        size_t max_rank = xindex.max_path_rank();
        for (size_t rank = 1; rank <= max_rank; ++rank) {
            Region region;
            region.seq = xindex.path_name(rank);
            regions.push_back(region);
        }
    }

    // validate and fill in sizes for regions that span entire path
    if (!id_range && !trace_mode) {
        for (auto& region : regions) {
            size_t rank = xindex.path_rank(region.seq);
            if (rank == 0) {
                cerr << "error[vg chunk]: input path " << region.seq << " not found in xg index" << endl;
                return 1;
            }
            if (region.start == 0 || region.end == -1) {
                region.start = 1;
                region.end = xindex.path_length(rank);
            } else if (!id_range) {
                if (region.start < 0 || region.end >= xindex.path_length(rank)) {
                    cerr << "error[vg chunk]: input region " << region.seq << ":" << region.start << "-" << region.end
                         << " is out of bounds of path " << region.seq << " which has length "<< xindex.path_length(rank)
                         << endl;
                    return -1;
                }
            }
        }
    }

    // finally, apply chunk_size and overlap to all regions if they are specified
    if (chunk_size > 0 && !trace_mode) {
        vector<Region> chunked_regions;
        for (auto& region : regions) {
            if (region.end - region.start <= chunk_size) {
                chunked_regions.push_back(region);
            } else {
                for (size_t pos = 1; pos < region.end; pos += chunk_size - overlap) {
                    Region cr = region;
                    cr.start = pos;
                    cr.end = min((int64_t)region.end, (int64_t)(pos + chunk_size - 1));
                    chunked_regions.push_back(cr);
                }
            }
        }
        swap(regions, chunked_regions);
    }

    // now ready to get our chunk on

    // what's the name of chunk i? 
    function<string(int, const Region&, string)> chunk_name =
        [&out_chunk_prefix](int i, const Region& region, string ext) -> string {
        stringstream chunk_name;
        string seq = region.seq.empty() ? "ids" : region.seq;
        chunk_name << out_chunk_prefix << "_" << i << "_" << seq << "_"
        << (region.start - 1) << "_" << region.end << ext;
        return chunk_name.str();
    };

    int num_regions = trace_mode ? trace_starts.size() : regions.size();

    // because we are expanding context, and not cutting nodes, our output
    // chunks are going to cover larger regions that what was asked for.
    // we return this in a bed file. 
    vector<Region> output_regions(num_regions);

    // initialize chunkers
    vector<PathChunker> chunkers(threads);
    for (auto& chunker : chunkers) {
        chunker.xg = &xindex;
    }

    
    // extract chunks in parallel
#pragma omp parallel for
    for (int i = 0; i < num_regions; ++i) {
        int tid = omp_get_thread_num();
        Region& region = regions[i];
        PathChunker& chunker = chunkers[tid];
        VG* subgraph = NULL;
        map<string, int> trace_thread_frequencies;
        if (trace_mode == true) {
            Graph g;
            trace_haplotypes_and_paths(xindex, trace_starts[i], trace_extend,
                                       g, trace_thread_frequencies);
            subgraph = new VG();
            subgraph->extend(g);
            
            // shoehorn trace info into region, just for sake of naming output chunk files
            output_regions[i].seq = "trace";
            output_regions[i].start = trace_starts[i] + 1; // because we -1 in chunk_name()
            output_regions[i].end = trace_extend;
        }
        else if (id_range == false) {
            subgraph = new VG();
            chunker.extract_subgraph(region, context_steps, *subgraph, output_regions[i]);
        } else {
            if (chunk_graph || context_steps > 0) {
                subgraph = new VG();
                output_regions[i].seq = region.seq;                
                chunker.extract_id_range(region.start, region.end, context_steps, *subgraph, output_regions[i]);
            } else {
                // in this case, there's no need to actually build the subgraph, so we don't
                // in order to save time.
                output_regions[i] = region;
            }
        }

        ofstream out_file;
        ostream* out_stream = NULL;
        if (chunk_graph) {
            if ((!region_string.empty() || !node_range_string.empty() || trace_start > 0) &&
                (regions.size() + trace_starts.size() == 1) && chunk_size == 0) {
                // if only one chunk passed in using -p, we output chunk to stdout
                out_stream = &cout;
            } else {
                // otherwise, we write files using prefix-i-seq-start-end convention.
                string name = chunk_name(i, output_regions[i], ".vg");
                out_file.open(name);
                if (!out_file) {
                    cerr << "error[vg chunk]: can't open output chunk file " << name << endl;
                    exit(1);
                }
                out_stream = &out_file;
            }
            
            subgraph->serialize_to_ostream(*out_stream);
        }
        
        // optional gam chunking
        if (chunk_gam) {
            string gam_name = chunk_name(i, output_regions[i], ".gam");
            ofstream out_gam_file(gam_name);
            if (!out_gam_file) {
                cerr << "error[vg chunk]: can't open output gam file " << gam_name << endl;
                exit(1);
            }
            if (subgraph != NULL) {
                chunker.extract_gam_for_subgraph(*subgraph, gam_index, &out_gam_file);
            } else {
                assert(id_range == true);
                chunker.extract_gam_for_id_range(region.start, region.end, gam_index, &out_gam_file);
            }
        }

        // trace annotations
        if (trace_mode) {
            string annot_name = chunk_name(i, output_regions[i], "_trace_annotate.txt");
            ofstream out_annot_file(annot_name);
            if (!out_annot_file) {
                cerr << "error[vg chunk]: can't open output trace annotation file " << annot_name << endl;
                exit(1);
            }
            for (auto tf : trace_thread_frequencies) {
                out_annot_file << tf.first << "\t" << tf.second << endl;
            }
        }

        delete subgraph;
    }
        
    // write a bed file if asked giving a more explicit linking of chunks to files
    if (!out_bed_file.empty()) {
        ofstream obed(out_bed_file);
        if (!obed) {
            cerr << "error[vg chunk]: can't open output bed file: " << out_bed_file << endl;
        }
        for (int i = 0; i < num_regions; ++i) {
            const Region& oregion = output_regions[i];
            string seq = id_range ? "ids" : oregion.seq;
            obed << seq << "\t" << (oregion.start - 1) << "\t" << oregion.end
                 << "\t" << chunk_name(i, oregion, chunk_gam ? ".gam" : ".vg");
            if (trace_mode) {
                obed << "\t" << chunk_name(i, oregion, "_trace_annotate.txt");
            }
            obed << "\n";
        }
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_chunk("chunk", "split graph or alignment into chunks", main_chunk);


