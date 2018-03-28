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
         << "Graph chunks are saved to .vg files, read chunks are saved to .gam files, and haplotype annotations are " << endl
         << "saved to .annotate.txt files, of the form <BASENAME>-<i>-<region name or \"ids\">-<start>-<length>.<ext>." << endl
         << "The BASENAME is specified with -b and defaults to \"./chunks\"." << endl
         << "For a single-range chunk (-p or -r), the graph data is sent to standard output instead of a file." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE       use this xg index to chunk subgraphs" << endl
         << "    -G, --gbwt-name FILE     use this GBWT haplotype index for haplotype extraction" << endl
         << "    -a, --gam-index FILE     chunk this gam index (made with vg index -a) instead of the graph" << endl
         << "    -g, --gam-and-graph      when used in combination with -a, both gam and graph will be chunked" << endl 
         << "path chunking:" << endl
         << "    -p, --path TARGET        write the chunk in the specified (0-based inclusive)\n"
         << "                             path range TARGET=path[:pos1[-pos2]] to standard output" << endl
         << "    -P, --path-list FILE     write chunks for all path regions in (line - separated file). format" << endl
         << "                             for each as in -p (all paths chunked unless otherwise specified)" << endl
         << "    -e, --input-bed FILE     write chunks for all (0-based end-exclusive) bed regions" << endl
         << "id range chunking:" << endl
         << "    -r, --node-range N:M     write the chunk for the specified node range to standard output\n"
         << "    -R, --node-ranges FILE   write the chunk for each node range in (newline or whitespace separated) file" << endl
         << "    -n, --n-chunks N         generate this many id-range chunks, which are determined using the xg index" << endl
         << "general:" << endl
         << "    -s, --chunk-size N       create chunks spanning N bases (or nodes with -r/-R) for all input regions." << endl
         << "    -o, --overlap N          overlap between chunks when using -s [0]" << endl        
         << "    -E, --output-bed FILE    write all created chunks to a bed file" << endl
         << "    -b, --prefix BASENAME    write output chunk files with the given base name. Files for chunk i will" << endl
         << "                             be named: <BASENAME>-<i>-<name>-<start>-<length>.<ext> [./chunk]" << endl
         << "    -c, --context-steps N    expand the context of the chunk this many node steps [1]" << endl
         << "    -l, --context-length N   expand the context of the chunk by this many bp [0]" << endl
         << "    -T, --trace              trace haplotype threads in chunks (and only expand forward from input coordinates)." << endl
         << "                             Produces a .annotate.txt file with haplotype frequencies for each chunk." << endl 
         << "    -f, --fully-contained    only return GAM alignments that are fully contained within chunk" << endl
         << "    -A, --search-all         search all nodes of alignment as opposed just the minimum (gam index must be made with -N)" << endl
         << "    -t, --threads N          for tasks that can be done in parallel, use this many threads [1]" << endl
         << "    -h, --help" << endl;
}

int main_chunk(int argc, char** argv) {

    if (argc == 2) {
        help_chunk(argv);
        return 1;
    }

    string xg_file;
    string gbwt_file;
    string gam_file;
    bool gam_and_graph = false;
    string region_string;
    string path_list_file;
    int chunk_size = 0;
    int overlap = 0;
    string in_bed_file;
    string out_bed_file;
    string out_chunk_prefix = "./chunk";
    int context_steps = -1;
    int context_length = 0;
    bool id_range = false;
    string node_range_string;
    string node_ranges_file;
    int threads = 1;
    bool trace = false;
    bool fully_contained = false;
    bool search_all_positions = false;
    int n_chunks = 0;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gbwt-name", required_argument, 0, 'G'},
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
            {"trace", required_argument, 0, 'T'},
            {"fully-contained", no_argument, 0, 'f'},
            {"search-all-positions", no_argument, 0, 'A'},
            {"threads", required_argument, 0, 't'},
            {"n-chunks", required_argument, 0, 'n'},
            {"context-length", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:G:a:gp:P:s:o:e:E:b:c:r:R:TfAt:n:l:",
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_file = optarg;
            break;
        
        case 'G':
            gbwt_file = optarg;
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

        case 'l':
            context_length = atoi(optarg);
            context_steps = (context_steps > 0 ? context_steps : -1);
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
            n_chunks = atoi(optarg);
            id_range = true;
            break;

        case 'T':
            trace = true;
            break;

        case 'f':
            fully_contained = true;
            break;

        case 'A':
            search_all_positions = true;
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

    // need at most one of -n, -p, -P, -e, -r, -R,  as an input
    if ((n_chunks == 0 ? 0 : 1) + (region_string.empty() ? 0 : 1) + (path_list_file.empty() ? 0 : 1) + (in_bed_file.empty() ? 0 : 1) +
        (node_ranges_file.empty() ? 0 : 1) + (node_range_string.empty() ? 0 : 1) > 1) {
        cerr << "error:[vg chunk] at most one of {-n, -p, -P, -e, -r, -R} required to specify input regions" << endl;
        return 1;
    }
    // need -a if using -f
    if (fully_contained && gam_file.empty()) {
        cerr << "error:[vg chunk] gam file must be specified with -a when using -f" << endl;
        return 1;
    }
    // context steps default to 1 if using id_ranges.  otherwise, force user to specify to avoid
    // misunderstandings
    if (context_steps < 0) {
        if (id_range) {
            if (!context_length) {
                context_steps = 1;
            }
        } else {
            cerr << "error:[vg chunk] context expansion steps must be specified with -c/--context when chunking on paths" << endl;
            return 1;
        }
    }

    // figure out which outputs we want.  the graph always
    // needs to be chunked, even if only gam output is requested,
    // because we use the graph to get the nodes we're looking for.
    // but we only write the subgraphs to disk if chunk_graph is true. 
    bool chunk_gam = !gam_file.empty();
    bool chunk_graph = gam_and_graph || !chunk_gam;

    // Load our index
    xg::XG xindex;
    if (chunk_graph || trace || context_steps > 0 || context_length > 0 || !id_range) {
        if (xg_file.empty()) {
            cerr << "error:[vg chunk] xg index (-x) required" << endl;
            return 1;
        }

        ifstream in(xg_file.c_str());
        if (!in) {
            cerr << "error:[vg chunk] unable to load xg index file " << xg_file << endl;
            return 1;
        }
        xindex.load(in);
        in.close();
    }

    // Now load the haplotype data
    unique_ptr<gbwt::GBWT> gbwt_index;
    if (trace && !gbwt_file.empty()) {
        // We are tracing haplotypes, and we want to use the GBWT instead of the old gPBWT.
        gbwt_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
        
        // Open up the index
        ifstream in(gbwt_file.c_str());
        if (!in) {
            cerr << "error:[vg chunk] unable to load gbwt index file " << gbwt_file << endl;
            return 1;
        }

        // And load it
        gbwt_index->load(in);
    }

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
            }
        }
    }
    else if (!in_bed_file.empty()) {
        parse_bed_regions(in_bed_file, regions);
    }
    else if (id_range) {
        if (n_chunks) {
            // determine the ranges from the xg index itself
            // how many nodes per range?
            int nodes_per_chunk = xindex.node_count / n_chunks;
            size_t i = 1;
            // iterate through the node ranks to build the regions
            while (i < xindex.node_count) {
                // make a range from i to i+nodeS_per_range
                vg::id_t a = xindex.rank_to_id(i);
                size_t j = i + nodes_per_chunk;
                if (j > xindex.node_count) j = xindex.node_count;
                vg::id_t b = xindex.rank_to_id(j);
                Region region;
                region.start = a;
                region.end = b;
                regions.push_back(region);
                i = j + 1;
            }
        } else {
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
    if (!id_range) {
        for (auto& region : regions) {
            size_t rank = xindex.path_rank(region.seq);
            if (rank == 0) {
                cerr << "error[vg chunk]: input path " << region.seq << " not found in xg index" << endl;
                return 1;
            }
            region.start = max((int64_t)0, region.start);
            if (region.end == -1) {    
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
    if (chunk_size > 0) {
        vector<Region> chunked_regions;
        for (auto& region : regions) {
            if (region.end - region.start <= chunk_size) {
                chunked_regions.push_back(region);
            } else {
                for (size_t pos = 0; pos < region.end; pos += chunk_size - overlap) {
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
        << region.start << "_" << region.end << ext;
        return chunk_name.str();
    };

    int num_regions = regions.size();

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
        if (id_range == false) {
            subgraph = new VG();
            chunker.extract_subgraph(region, context_steps, context_length,
                                     trace, *subgraph, output_regions[i]);
        } else {
            if (chunk_graph || context_steps > 0) {
                subgraph = new VG();
                output_regions[i].seq = region.seq;
                chunker.extract_id_range(region.start, region.end,
                                         context_steps, context_length, trace,
                                         *subgraph, output_regions[i]);
            } else {
                // in this case, there's no need to actually build the subgraph, so we don't
                // in order to save time.
                output_regions[i] = region;
            }
        }

        // optionally trace our haplotypes
        if (trace && subgraph) {
            int64_t trace_start;
            int64_t trace_end;
            if (id_range) {
                trace_start = output_regions[i].start;
                trace_end = output_regions[i].end;
            } else {
                trace_start = xindex.node_at_path_position(output_regions[i].seq,
                                                           output_regions[i].start);
                trace_end = xindex.node_at_path_position(output_regions[i].seq,
                                                         output_regions[i].end);
            }
            int64_t trace_steps = trace_end - trace_start;
            Graph g;
            trace_haplotypes_and_paths(xindex, gbwt_index.get(), trace_start, trace_steps,
                                       g, trace_thread_frequencies, false);
            subgraph->paths.for_each([&trace_thread_frequencies](const Path& path) {
                    trace_thread_frequencies[path.name()] = 1;});            
            subgraph->extend(g);
        }

        ofstream out_file;
        ostream* out_stream = NULL;
        if (chunk_graph) {
            if ((!region_string.empty() || !node_range_string.empty()) &&
                (regions.size()  == 1) && chunk_size == 0) {
                // If we are going to output only one chunk, it should go to
                // stdout instead of to a file on disk
                out_stream = &cout;
            } else {
                // Otherwise, we write files under the specified prefix, using
                // a prefix-i-seq-start-end convention.
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
                chunker.extract_gam_for_subgraph(*subgraph, gam_index, &out_gam_file,
                                                 fully_contained, search_all_positions);
            } else {
                assert(id_range == true);
                vector<vg::id_t> region_id_range = {region.start, region.end};
                chunker.extract_gam_for_ids(region_id_range, gam_index, &out_gam_file,
                                            true, fully_contained, search_all_positions);
            }
        }

        // trace annotations
        if (trace) {
            // Even if we have only one chunk, the trace annotation data always
            // ends up in a file.
            string annot_name = chunk_name(i, output_regions[i], ".annotate.txt");
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
            obed << seq << "\t" << oregion.start << "\t" << (oregion.end + 1)
                 << "\t" << chunk_name(i, oregion, chunk_gam ? ".gam" : ".vg");
            if (trace) {
                obed << "\t" << chunk_name(i, oregion, ".annotate.txt");
            }
            obed << "\n";
        }
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_chunk("chunk", "split graph or alignment into chunks", main_chunk);


