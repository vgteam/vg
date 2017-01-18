// chunk.cpp: define the "vg chunk" subcommand, which chunks up vg graphs
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
         << "    -P, --path-list FILE     write chunks for all path regions in (newline or whitespace\n"
         << "                             separated file). format for each as in -p\n" << endl
         << "    -s, --chunk-size N       create chunks spanning N bases. applies to all regions specified with other flags.\n"
         << "    -o, --overlap N          overlap between chunks when using -s [defaul = 0]" << endl
         << "                             (all paths chunked unless otherwise specified)" << endl
         << "    -r, --input-bed FILE     write chunks for all (0-based) bed regions" << endl
         << "    -R, --output-bed FILE    write all created chunks to a bed file" << endl
         << "    -b, --prefix PATH        prefix output chunk paths with PATH. each chunk will have the following name:\n"
         << "                             <PATH>-<i>-<name>-<start>-<length>. <i> is the line# of the chunk in the\n"
         << "                             output bed. [default=./chunk]" << endl
         << "    -c, --context STEPS      expand the context of the chunk this many steps [10]" << endl
         << "general:" << endl
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
    int context_steps = 10;
    int threads = 1;
    
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
            {"input-bed", required_argument, 0, 'r'},
            {"output-bed", required_argument, 0, 'R'},
            {"prefix", required_argument, 0, 'b'},
            {"context", required_argument, 0, 'c'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:a:gp:P:s:o:r:R:b:c:t:",
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

        case 'r':
            in_bed_file = optarg;
            break;

        case 'R':
            out_bed_file = optarg;
            break;

        case 'b':
            out_chunk_prefix = optarg;
            break;

        case'c':
            context_steps = atoi(optarg);
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
    // need at most one of -p, -P, -r as an input
    if ((region_string.empty() ? 0 : 1) + (path_list_file.empty() ? 0 : 1) + (in_bed_file.empty() ? 0 : 1) > 1) {
        cerr << "error:[vg chunk] at most one of {-p, -P, -r} required to specify input regions" << endl;
        return 1;
    }

    // figure out which outputs we want.  the graph always
    // needs to be chunked, even if only gam output is requested,
    // because we use the graph to get the nodes we're looking for.
    // but we only write the subgraphs to disk if chunk_graph is true. 
    bool chunk_gam = !gam_file.empty();
    bool chunk_graph = gam_and_graph || !chunk_gam;

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
                cerr << "push " << buf << endl;
            }
        }
    }
    else if (!in_bed_file.empty()) {
        parse_bed_regions(in_bed_file, regions);
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
    for (auto& region : regions) {
        size_t rank = xindex.path_rank(region.seq);
        if (rank == 0) {
            cerr << "error[vg chunk]: input path " << region.seq << " not found in xg index" << endl;
            return 1;
        }
        if (region.start == 0 || region.end == -1) {
            region.start = 1;
            region.end = xindex.path_length(rank);
        } else {
            if (region.start < 0 || region.end >= xindex.path_length(rank)) {
                cerr << "error[vg chunk]: input region " << region.seq << ":" << region.start << "-" << region.end
                     << " is out of bounds of path " << region.seq << " which has length "<< xindex.path_length(rank)
                     << endl;
                return -1;
            }
        }
    }

    // finally, apply chunk_size and overlap to all regions if they are specified
    if (chunk_size > 0) {
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
    function<string(int, const Region&, bool)> chunk_name =
        [&out_chunk_prefix](int i, const Region& region, bool gam) -> string {
        stringstream chunk_name;
        chunk_name << out_chunk_prefix << "_" << i << "_" << region.seq << "_"
        << (region.start - 1) << "_" << region.end << (gam ? ".gam" : ".vg");
        return chunk_name.str();
    };

    // because we are expanding context, and not cutting nodes, our output
    // chunks are going to cover larger regions that what was asked for.
    // we return this in a bed file. 
    vector<Region> output_regions(regions.size());

    // initialize chunkers
    vector<PathChunker> chunkers(threads);
    for (auto& chunker : chunkers) {
        chunker.xg = &xindex;
    }

    // extract chunks in parallel
#pragma omp parallel for
    for (int i = 0; i < regions.size(); ++i) {
        int tid = omp_get_thread_num();
        Region& region = regions[i];
        PathChunker& chunker = chunkers[tid];
        VG subgraph;
        output_regions[i].seq = regions[i].seq;
        output_regions[i].start = 1 + chunker.extract_subgraph(region, context_steps, subgraph);
        output_regions[i].end = output_regions[i].start;
        // Is there a better way to get path length? 
        Path output_path = subgraph.paths.path(region.seq);
        for (int j = 0; j < output_path.mapping_size(); ++j) {
            int64_t op_node = output_path.mapping(j).position().node_id();
            output_regions[i].end += xindex.node_length(op_node);
        }

        ofstream out_file;
        ostream* out_stream = NULL;
        if (chunk_graph) {
            if (!region_string.empty() && regions.size() == 1 && chunk_size == 0) {
                // if only one chunk passed in using -p, we output chunk to stdout
                out_stream = &cout;
            } else {
                // otherwise, we write files using prefix-i-seq-start-end convention.
                string name = chunk_name(i, output_regions[i], false);
                out_file.open(name);
                if (!out_file) {
                    cerr << "error[vg chunk]: can't open output chunk file " << name << endl;
                    exit(1);
                }
                out_stream = &out_file;
            }
            
            subgraph.serialize_to_ostream(*out_stream);
        }
        
        // optional gam chunking
        if (chunk_gam) {
            string gam_name = chunk_name(i, output_regions[i], true);
            ofstream out_gam_file(gam_name);
            if (!out_gam_file) {
                cerr << "error[vg chunk]: can't open output gam file " << gam_name << endl;
                exit(1);
            }
            chunker.extract_gam_for_subgraph(subgraph, gam_index, &out_gam_file);
        }
    }
        
    // write a bed file if asked giving a more explicit linking of chunks to files
    if (!out_bed_file.empty()) {
        ofstream obed(out_bed_file);
        if (!obed) {
            cerr << "error[vg chunk]: can't open output bed file: " << out_bed_file << endl;
        }
        for (int i = 0; i < regions.size(); ++i) {
            const Region& oregion = output_regions[i];
            obed << oregion.seq << "\t" << (oregion.start - 1) << "\t" << oregion.end
                 << "\t" << chunk_name(i, oregion, chunk_gam) << "\n";
        }
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_chunk("chunk", "split graph or alignment into chunks", main_chunk);


