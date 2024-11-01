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
#include <vg/io/stream.hpp>
#include <vg/io/stream_multiplexer.hpp>
#include <vg/io/vpkg.hpp>
#include "../utility.hpp"
#include "../chunker.hpp"
#include "../stream_index.hpp"
#include "../region.hpp"
#include "../haplotype_extracter.hpp"
#include "../algorithms/sorted_id_ranges.hpp"
#include "../algorithms/find_gbwt.hpp"
#include <bdsg/overlays/overlay_helper.hpp>
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

static string chunk_name(const string& out_chunk_prefix, int i, const Region& region, string ext, int gi = 0, bool components = false);
static int split_gam(istream& gam_stream, size_t chunk_size, const string& out_prefix,
                     size_t gam_buffer_size = 100);
static void check_read(const Alignment& aln, const HandleGraph* graph);
                     

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
         << "    -x, --xg-name FILE       use this graph or xg index to chunk subgraphs" << endl
         << "    -G, --gbwt-name FILE     use this GBWT haplotype index for haplotype extraction (for -T)" << endl
         << "    -a, --gam-name FILE      chunk this gam file instead of the graph (multiple allowed)" << endl
         << "    -g, --gam-and-graph      when used in combination with -a, both gam and graph will be chunked" << endl 
         << "    -F, --in-gaf             input alignment is a sorted bgzipped GAF" << endl 
         << "path chunking:" << endl
         << "    -p, --path TARGET        write the chunk in the specified (0-based inclusive, multiple allowed)\n"
         << "                             path range TARGET=path[:pos1[-pos2]] to standard output" << endl
         << "    -P, --path-list FILE     write chunks for all path regions in (line - separated file). format" << endl
         << "                             for each as in -p (all paths chunked unless otherwise specified)" << endl
         << "    -e, --input-bed FILE     write chunks for all (0-based end-exclusive) bed regions" << endl
         << "    -S, --snarls FILE        write given path-range(s) and all snarls fully contained in them, as alternative to -c" << endl
         << "id range chunking:" << endl
         << "    -r, --node-range N:M     write the chunk for the specified node range to standard output\n"
         << "    -R, --node-ranges FILE   write the chunk for each node range in (newline or whitespace separated) file" << endl
         << "    -n, --n-chunks N         generate this many id-range chunks, which are determined using the xg index" << endl
         << "simple gam chunking:" << endl
         << "    -m, --gam-split-size N   split gam (specified with -a, sort/index not required) up into chunks with at most N reads each" << endl
         << "component chunking:" << endl
         << "    -C, --components         create a chunk for each connected component.  if a targets given with (-p, -P, -r, -R), limit to components containing them" << endl
         << "    -M, --path-components    create a chunk for each path in the graph's connected component" << endl
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
         << "    --no-embedded-haplotypes Don't load haplotypes from the graph. It is possible to -T without any haplotypes available." << endl
         << "    -f, --fully-contained    only return GAM alignments that are fully contained within chunk" << endl
         << "    -O, --output-fmt         Specify output format (vg, pg, hg, gfa).  [pg (vg with -T)]" << endl
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
    vector<string> gam_files;
    bool gam_and_graph = false;
    bool gam_is_gaf = false;
    vector<string> region_strings;
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
    bool trace = false;
    bool no_embedded_haplotypes = false;
    bool fully_contained = false;
    int n_chunks = 0;
    size_t gam_split_size = 0;
    string output_format = "pg";
    bool output_format_set = false;
    bool components = false;
    bool path_components = false;
    string snarl_filename;
    
    #define OPT_NO_EMBEDDED_HAPLOTYPES 1000
    
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
            {"in-gaf", no_argument, 0, 'F'},
            {"path", required_argument, 0, 'p'},
            {"path-names", required_argument, 0, 'P'},
            {"chunk-size", required_argument, 0, 's'},
            {"overlap", required_argument, 0, 'o'},
            {"input-bed", required_argument, 0, 'e'},
            {"snarls", required_argument, 0, 'S'},
            {"output-bed", required_argument, 0, 'E'},
            {"prefix", required_argument, 0, 'b'},
            {"context", required_argument, 0, 'c'},
            {"id-ranges", no_argument, 0, 'r'},
            {"id-range", no_argument, 0, 'R'},
            {"trace", no_argument, 0, 'T'},
            {"no-embedded-haplotypes", no_argument, 0, OPT_NO_EMBEDDED_HAPLOTYPES},
            {"fully-contained", no_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},
            {"n-chunks", required_argument, 0, 'n'},
            {"context-length", required_argument, 0, 'l'},
            {"gam-split-size", required_argument, 0, 'm'},
            {"components", no_argument, 0, 'C'},
            {"path-components", no_argument, 0, 'M'},
            {"output-fmt", required_argument, 0, 'O'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:G:a:gFp:P:s:o:e:S:E:b:c:r:R:Tft:n:l:m:CMO:",
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
            gam_files.push_back(optarg);
            break;
            
        case 'g':
            gam_and_graph = true;
            break;            

        case 'F':
            gam_is_gaf = true;
            break;            

        case 'p':
            region_strings.push_back(optarg);
            break;

        case 'P':
            path_list_file = optarg;
            break;

        case 's':
            chunk_size = parse<int>(optarg);
            break;

        case 'o':
            overlap = parse<int>(optarg);
            break;

        case 'e':
            in_bed_file = optarg;
            break;
            
        case 'S':
            snarl_filename = optarg;
            break;
            
        case 'E':
            out_bed_file = optarg;
            break;

        case 'b':
            out_chunk_prefix = optarg;
            break;

        case'c':
            context_steps = parse<int>(optarg);
            break;

        case 'l':
            context_length = parse<int>(optarg);
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
            n_chunks = parse<int>(optarg);
            id_range = true;
            break;

        case 'm':
            gam_split_size = parse<int>(optarg);
            break;

        case 'C':
            components = true;
            break;

        case 'M':
            components = true;
            path_components = true;
            break;

        case 'T':
            trace = true;
            break;
            
        case OPT_NO_EMBEDDED_HAPLOTYPES:
            no_embedded_haplotypes = true;
            break;

        case 'f':
            fully_contained = true;
            break;

        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg chunk] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
            
        case 'O':
            output_format = optarg;
            output_format_set = true;
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

    // need at most one of -n, -p, -P, -e, -r, -R, -m  as an input
    if ((n_chunks == 0 ? 0 : 1) + (region_strings.empty() ? 0 : 1) + (path_list_file.empty() ? 0 : 1) + (in_bed_file.empty() ? 0 : 1) +
        (node_ranges_file.empty() ? 0 : 1) + (node_range_string.empty() ? 0 : 1) + (gam_split_size == 0 ? 0 : 1)  +
        (path_components ? 1 : 0) > 1) {
        cerr << "error:[vg chunk] at most one of {-n, -p, -P, -e, -r, -R, -m, '-M'} required to specify input regions" << endl;
        return 1;
    }
    // need -a if using -f
    if ((gam_split_size != 0 || fully_contained) && gam_files.empty()) {
        cerr << "error:[vg chunk] read alignment file (gam or gaf) must be specified with -a when using -f or -m" << endl;
        return 1;
    }
    if (components == true && context_steps >= 0) {
        cerr << "error:[vg chunk] context cannot be specified (-c) when splitting into components (-C)" << endl;
        return 1;
    }

    if (!snarl_filename.empty() && context_steps >= 0) {
        cerr << "error:[vg chunk] context cannot be specified (-c) when using snarls (-S)" << endl;
        return 1;
    }
    if (!snarl_filename.empty() && region_strings.empty() && path_list_file.empty() && in_bed_file.empty()) {
        cerr << "error:[vg chunk] snarl chunking can only be used with path regions (-p -P  -e)" << endl;
        return 1;        
    }

    // check the output format
    std::transform(output_format.begin(), output_format.end(), output_format.begin(), ::tolower);
    if (!vg::io::valid_output_format(output_format)) {
        cerr << "error[vg chunk]: invalid output format" << endl;
        return 1;
    }
    if (trace && output_format != "vg") {
        // todo: trace code goes through vg conversion anyway and according to unit tests
        //       fails when not outputting vg
        output_format = "vg";
        if (output_format_set) {
            cerr << "warning[vg chunk]: ignoring -O and setting output format to vg, as required by -T" << endl;
        }
        
    }
    else if (output_format == "vg") {
        cerr << "warning[vg chunk]: the vg-protobuf format is DEPRECATED. you probably want to use PackedGraph (pg) instead" << endl;
    }    
    string output_ext = output_format == "gfa" ? ".gfa"  : ".vg";

    // figure out which outputs we want.  the graph always
    // needs to be chunked, even if only gam output is requested,
    // because we use the graph to get the nodes we're looking for.
    // but we only write the subgraphs to disk if chunk_graph is true. 
    bool chunk_gam = !gam_files.empty() && gam_split_size == 0;
    bool chunk_graph = gam_and_graph || (!chunk_gam && gam_split_size == 0);

    // parse the regions into a list before loading the graph, if we're
    // specifying regions by path name.
    vector<Region> regions;
    if (!region_strings.empty()) {
        for (auto& region_string : region_strings) {
            Region region;
            parse_region(region_string, region);
            regions.push_back(region);
        }
    }
    if (!path_list_file.empty()) {
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
    if (!in_bed_file.empty()) {
        parse_bed_regions(in_bed_file, regions);
    }

    // Load the snarls
    unique_ptr<SnarlManager> snarl_manager;
    if (!snarl_filename.empty()) {
        ifstream snarl_file(snarl_filename.c_str());
        if (!snarl_file) {
            cerr << "error:[vg chunk] Unable to load snarls file: " << snarl_filename << endl;
            return 1;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    }

    // Load our index
    PathPositionHandleGraph* graph = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::ReferencePathOverlayHelper overlay_helper;

    if (chunk_graph || trace || context_steps > 0 || context_length > 0 || (!id_range && gam_split_size == 0) || (id_range && chunk_gam) || components) {
        if (xg_file.empty()) {
            cerr << "error:[vg chunk] graph or xg index (-x) required" << endl;
            return 1;
        }

        ifstream in(xg_file.c_str());
        if (!in) {
            cerr << "error:[vg chunk] unable to load graph / xg index file " << xg_file << endl;
            return 1;
        }
        in.close();

        // To support the regions we were asked for, we might need to ensure
        // the paths they are on are actually indexed for reference style
        // offset lookups.
        std::unordered_set<std::string> ensure_indexed;
        for (auto& region : regions) {
            ensure_indexed.insert(region.seq);
        }
        
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_file);
        graph = overlay_helper.apply(path_handle_graph.get(), ensure_indexed);
        in.close();
    }

    // Now load the haplotype data
    
    unique_ptr<gbwt::GBWT> gbwt_index_holder;
    const gbwt::GBWT* gbwt_index = nullptr; 
    if (trace) {
        // We are tracing haplotypes.
        // We might want a GBWT.
        
        // TODO: Make the tube map stop calling us in trace mode *without* a GBWT?
        
        if (!gbwt_file.empty()) {
            // A GBWT file is specified, so load that
            gbwt_index_holder = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_file);
            if (gbwt_index_holder.get() == nullptr) {
                // Complain if we couldn't get it but were supposed to.
                cerr << "error:[vg::chunk] unable to load gbwt index file " << gbwt_file << endl;
                exit(1);
            }
            gbwt_index = gbwt_index_holder.get();
        }
        
        if (!gbwt_index && !no_embedded_haplotypes) {
            // We didn't get a GBWT from a file, and we are allowed to use the one in the graph, if any.
            gbwt_index = vg::algorithms::find_gbwt(path_handle_graph.get());
        }
        
        // It's OK if gbwt_index is still null here!
    }

    
    // We need an index on the GAM to chunk it (if we're not doing components)
    vector<unique_ptr<GAMIndex>> gam_indexes;
    vector<unique_ptr<tbx_t>> gaf_tbxs;
    vector<unique_ptr<htsFile>> gaf_fps;
    if (chunk_gam && !components) {
        if (gam_is_gaf){
            for (auto gaf_file : gam_files) {
                try {
                    tbx_t *gaf_tbx = NULL;
                    htsFile *gaf_fp = NULL;
                    if (!gaf_file.empty()){
                        gaf_tbx = tbx_index_load3(gaf_file.c_str(), NULL, 0);
                        if ( !gaf_tbx ){
                            cerr << "Could not load .tbi/.csi index of " << gaf_file << endl;
                            exit(1);
                        }
                        int nseq;
                        gaf_fp = hts_open(gaf_file.c_str(),"r");
                        if ( !gaf_fp ) {
                            cerr << "Could not open " << gaf_file << endl;
                            exit(1);
                        }
                        gaf_fps.push_back(unique_ptr<htsFile>(gaf_fp));
                        gaf_tbxs.push_back(unique_ptr<tbx_t>(gaf_tbx));
                    }
                } catch (...) {
                    cerr << "error:[vg chunk] unable to load GAF index file: " << gaf_file << "" << endl;
                    exit(1);
                }
            }
        } else {
            for (auto gam_file : gam_files) {
                try {
                    get_input_file(gam_file + ".gai", [&](istream& index_stream) {
                        gam_indexes.push_back(unique_ptr<GAMIndex>(new GAMIndex()));
                        gam_indexes.back()->load(index_stream);
                    });
                } catch (...) {
                    cerr << "error:[vg chunk] unable to load GAM index file: " << gam_file << ".gai" << endl
                         << "                 note: a GAM index is required when *not* chunking by components with -C or -M" << endl;
                    exit(1);
                }
            }
        }
    }
    
    // If we're chunking on components with a GAM, this map will be used
    // (instead of an index)
    unordered_map<nid_t, int32_t> node_to_component;
    
    if (id_range) {
        if (n_chunks) {
            // Determine the ranges from the source graph itself.
            // how many nodes per range?
            size_t node_count = graph->get_node_count();
            size_t nodes_per_chunk = node_count / n_chunks;
            
            // We need to articulate our chunks in terms of ID ranges, but we
            // have no guarantee that the graph we pull from will be in ID
            // order. An XG probably ought to be in topological order anyway.
            // So we pull all the IDs and sort them in a big vector in order to
            // get the chunk ID breakpoints.
            vector<vg::id_t> rank_to_id(node_count + 1);
            size_t i = 1;
            graph->for_each_handle([&](handle_t handle) {
                rank_to_id[i++] = graph->get_id(handle);
            });
            
            // Sort so we can find the nth ID easily
            std::sort(rank_to_id.begin(), rank_to_id.end());
            
            i = 1;
            // iterate through the node ranks to build the regions
            while (i < node_count) {
                // make a range from i to i+nodeS_per_range
                vg::id_t a = rank_to_id[i];
                size_t j = i + nodes_per_chunk;
                if (j > node_count) j = node_count;
                vg::id_t b = rank_to_id[j];
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
    if (graph != nullptr && path_components) {
        // every reference or generic path (guaranteed to be reference indexed)
        graph->for_each_path_matching({PathSense::REFERENCE, PathSense::GENERIC}, {}, {}, [&](path_handle_t path_handle) {
                Region region;
                region.seq = graph->get_path_name(path_handle);
                if (!Paths::is_alt(region.seq)) {
                    regions.push_back(region);
                }
            });
    }
    
    if (context_steps >= 0 && regions.empty()) {
        cerr << "error:[vg chunk] extracting context (-c) requires a region to take context around" << endl;
        return 1;
    }
    
    // context steps default to 1 if using id_ranges.  otherwise, force user to specify to avoid
    // misunderstandings
    if (context_steps < 0 && gam_split_size == 0) {
        if (id_range) {
            if (!context_length) {
                context_steps = 1;
            }
        } else if (!components && snarl_filename.empty()){
            cerr << "error:[vg chunk] context (-c) or snarls (-S)  must be specified when chunking on paths" << endl;
            return 1;
        }
    }

    // validate and fill in sizes for regions that span entire path
    function<size_t(const string&)> get_path_length = [&](const string& path_name) {
        size_t path_length = 0;
        for (handle_t handle : graph->scan_path(graph->get_path_handle(path_name))) {
            path_length += graph->get_length(handle);
        }
        return path_length;
    };
    if (!id_range) {
        for (auto& region : regions) {
            if (!graph->has_path(region.seq)) {
                cerr << "error[vg chunk]: input path " << region.seq << " not found in graph" << endl;
                return 1;
            }
            region.start = max((int64_t)0, region.start);
            if (region.end == -1) {
                region.end = get_path_length(region.seq) - 1;
            } else if (!id_range && !components) {
                if (region.start < 0 || region.end >= get_path_length(region.seq)) {
                    cerr << "error[vg chunk]: input region " << region.seq << ":" << region.start << "-" << region.end
                         << " is out of bounds of path " << region.seq << " which has length "<< get_path_length(region.seq)
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

    // when using -C for components, regions will be derived from the connected components
    vector<unordered_set<nid_t>> component_ids; 
    if (components == true && regions.empty()) {
        // no regions given, we find our components from scratch and make some dummy regions
        component_ids = handlealgs::weakly_connected_components(graph);
        for (int i = 0; i < component_ids.size(); ++i) {
            Region region;
            region.seq = "";
            region.start = 0;
            region.end = 0;
            regions.push_back(region);
        }
    }

    // now ready to get our chunk on
    if (gam_split_size != 0) {
        if(gam_is_gaf){
            cerr << "error[vg chunk]: GAF file input toggled with -F but, currently, only GAM files can by split. A workaround would be to split the GAF file using split -l/-n which can split text files into chunks." << endl;
            return 1;
        }
        for (size_t gi = 0; gi < gam_files.size(); ++gi) {
            ifstream gam_stream;
            string& gam_file = gam_files[gi];
            // Open the GAM file, whether splitting directly or seeking with an index
            gam_stream.open(gam_file);
            if (!gam_stream) {
                cerr << "error[vg chunk]: unable to open input gam: " << gam_file << endl;
                return 1;
            }
            // just chunk up every N reads in the gam without any path or id logic. Don't do anything else.
            string prefix = gi == 0 ? out_chunk_prefix : out_chunk_prefix + std::to_string(gi);
            split_gam(gam_stream, gam_split_size, prefix);
        }
        return 0;
    }

    int num_regions = regions.size();

    // because we are expanding context, and not cutting nodes, our output
    // chunks are going to cover larger regions that what was asked for.
    // we return this in a bed file. 
    vector<Region> output_regions(num_regions);

    // initialize chunkers
    size_t threads = get_thread_count();
    vector<PathChunker> chunkers(threads);
    for (auto& chunker : chunkers) {
        chunker.graph = graph;
    }
    
    // When chunking GAMs, every thread gets its own cursor to seek into the input GAM.
    // Todo: when operating on multiple gams, we make |threads| X |gams| cursors, even though
    // we only ever use |threads| threads.
    vector<list<ifstream>> gam_streams_vec(gam_files.size());
    vector<vector<GAMIndex::cursor_t>> cursors_vec(gam_files.size());
    
    if (chunk_gam & !gam_is_gaf) {
        for (size_t gam_i = 0; gam_i < gam_streams_vec.size(); ++gam_i) {
            auto& gam_file = gam_files[gam_i];
            auto& gam_streams = gam_streams_vec[gam_i];
            auto& cursors = cursors_vec[gam_i];
            cursors.reserve(threads);
            for (size_t i = 0; i < threads; i++) {
                // Open a stream for every thread
                gam_streams.emplace_back(gam_file);
                if (!gam_streams.back()) {
                    cerr << "error[vg chunk]: unable to open GAM file " << gam_file << endl;
                    return 1;
                }
                // And wrap it in a cursor
                cursors.emplace_back(gam_streams.back());
            }
        }
    }

    // extract chunks in parallel
#pragma omp parallel for
    for (int i = 0; i < num_regions; ++i) {
        int tid = omp_get_thread_num();
        Region& region = regions[i];
        PathChunker& chunker = chunkers[tid];
        unique_ptr<MutablePathMutableHandleGraph> subgraph;
        map<string, int> trace_thread_frequencies;
        if (!component_ids.empty()) {
            subgraph = vg::io::new_output_graph<MutablePathMutableHandleGraph>(output_format);
            chunker.extract_component(component_ids[i], *subgraph, false);
            output_regions[i] = region;
        }
        else if (id_range == false) {
            subgraph = vg::io::new_output_graph<MutablePathMutableHandleGraph>(output_format);
            if (components == true) {
                chunker.extract_path_component(region.seq, *subgraph, output_regions[i]);
            } else if (snarl_manager.get() != nullptr) {
                chunker.extract_snarls(region, *snarl_manager, *subgraph);
                output_regions[i] = region;                
            } else {
                chunker.extract_subgraph(region, context_steps, context_length,
                                         trace, *subgraph, output_regions[i]);
            }
        } else {
            if (chunk_graph || context_steps > 0) {
                subgraph = vg::io::new_output_graph<MutablePathMutableHandleGraph>(output_format);
                output_regions[i].seq = region.seq;
                chunker.extract_id_range(region.start, region.end,
                                         components ? numeric_limits<int64_t>::max() : context_steps,
                                         context_length, trace && !components,
                                         *subgraph, output_regions[i]);
            } else {
                // in this case, there's no need to actually build the subgraph, so we don't
                // in order to save time.
                output_regions[i] = region;
            }
        }

        // optionally trace our haplotypes
        if (trace && subgraph && gbwt_index) {
            int64_t trace_start;
            int64_t trace_steps = 0;
            if (id_range) {
                trace_start = output_regions[i].start;
                trace_steps = output_regions[i].end - trace_start;
            } else {
                path_handle_t path_handle = graph->get_path_handle(output_regions[i].seq);
                step_handle_t trace_start_step = graph->get_step_at_position(path_handle, output_regions[i].start);
                step_handle_t trace_end_step = graph->get_step_at_position(path_handle, output_regions[i].end);
                // make sure we don't loop forever in next loop
                if (output_regions[i].start > output_regions[i].end) {
                    swap(trace_start_step, trace_end_step);
                }
                trace_start = graph->get_id(graph->get_handle_of_step(trace_start_step));
                for (; trace_start_step != trace_end_step; trace_start_step = graph->get_next_step(trace_start_step)) {
                    ++trace_steps;
                }
                // haplotype_extender is forward only.  until it's made bidirectional, try to
                // detect backward paths and trace them backwards.  this will not cover all possible cases though.
                if (graph->get_is_reverse(graph->get_handle_of_step(trace_start_step)) &&
                    graph->get_is_reverse(graph->get_handle_of_step(trace_end_step))) {
                    trace_start = graph->get_id(graph->get_handle_of_step(trace_end_step));
                }
            }
            Graph g;
            trace_haplotypes_and_paths(*graph, *gbwt_index, trace_start, trace_steps,
                                       g, trace_thread_frequencies, false);
            subgraph->for_each_path_handle([&trace_thread_frequencies, &subgraph](path_handle_t path_handle) {
                    trace_thread_frequencies[subgraph->get_path_name(path_handle)] = 1;});
            VG* vg_subgraph = dynamic_cast<VG*>(subgraph.get());
            if (vg_subgraph != nullptr) {
                // our graph is in vg format, just extend it
                vg_subgraph->extend(g);
            } else {
                // our graph is not in vg format.  covert it, extend it, convert it back
                // this can eventually be avoided by handlifying the haplotype tracer
                VG vg;
                handlealgs::copy_path_handle_graph(subgraph.get(), &vg);
                subgraph.reset();
                vg.extend(g);
                subgraph = vg::io::new_output_graph<MutablePathMutableHandleGraph>(output_format);
                handlealgs::copy_path_handle_graph(&vg, subgraph.get());
            }
        }

        ofstream out_file;
        ostream* out_stream = NULL;
        if (chunk_graph) {
            if ((!region_strings.empty() || !node_range_string.empty()) &&
                (regions.size()  == 1) && chunk_size == 0) {
                // If we are going to output only one chunk, it should go to
                // stdout instead of to a file on disk
                out_stream = &cout;
            } else {
                // Otherwise, we write files under the specified prefix, using
                // a prefix-i-seq-start-end convention.
                string name = chunk_name(out_chunk_prefix, i, output_regions[i], output_ext, 0, components);
                out_file.open(name);
                if (!out_file) {
                    cerr << "error[vg chunk]: can't open output chunk file " << name << endl;
                    exit(1);
                }
                out_stream = &out_file;
            }

            assert(subgraph);
            vg::io::save_handle_graph(subgraph.get(), *out_stream);
        }
        
        // optional gam chunking
        if (chunk_gam) {
            if (!components) {
                // Work out the ID ranges to look up
                vector<pair<vg::id_t, vg::id_t>> region_id_ranges;
                if (subgraph) {
                    // Use the regions from the graph
                    region_id_ranges = vg::algorithms::sorted_id_ranges(subgraph.get());
                } else {
                    // Use the region we were asked for
                    region_id_ranges = {{region.start, region.end}};
                }

                if(gam_is_gaf){
                    // use the indexed bgzipped GAFs
                    for (size_t gi = 0; gi < gaf_fps.size(); ++gi) {
                        auto& gaf_fp = gaf_fps[gi];
                        assert(gaf_fp.get() != nullptr);
                        auto& gaf_tbx = gaf_tbxs[gi];
                        assert(gaf_tbx.get() != nullptr);

                        string gaf_name = chunk_name(out_chunk_prefix, i, output_regions[i], ".gaf", gi, components);
                        ofstream out_gaf_file(gaf_name);
                        if (!out_gaf_file) {
                            cerr << "error[vg chunk]: can't open output gaf file " << gaf_name << endl;
                            exit(1);
                        }

                        // loop over ranges and print GAF records
                        for (auto range : region_id_ranges) {
                            string reg = "{node}:" + convert(range.first) + "-" + convert(range.second);
                            hts_itr_t *itr = tbx_itr_querys(gaf_tbx.get(), reg.c_str());
                            kstring_t str = {0,0,0};
                            if ( itr ) {
                                while (tbx_itr_next(gaf_fp.get(), gaf_tbx.get(), itr, &str) >= 0) {
                                    out_gaf_file << str.s << endl;
                                }
                                tbx_itr_destroy(itr);
                            }
                        }
                        out_gaf_file.close();
                    }
                } else {
                    // old way: use the gam index
                    for (size_t gi = 0; gi < gam_indexes.size(); ++gi) {
                        auto& gam_index = gam_indexes[gi];
                        assert(gam_index.get() != nullptr);
                        GAMIndex::cursor_t& cursor = cursors_vec[gi][tid];
            
                        string gam_name = chunk_name(out_chunk_prefix, i, output_regions[i], ".gam", gi, components);
                        ofstream out_gam_file(gam_name);
                        if (!out_gam_file) {
                            cerr << "error[vg chunk]: can't open output gam file " << gam_name << endl;
                            exit(1);
                        }
                        
                        auto emit = vg::io::emit_to<Alignment>(out_gam_file);
                    
                        auto handle_read = [&](const Alignment& aln) {
                            check_read(aln, graph);
                            emit(aln);
                        };
                        
                        gam_index->find(cursor, region_id_ranges, handle_read, fully_contained);
                    }
                }
            } else {
#pragma omp critical (node_to_component)
                {
                    // we're doing components, just use stl map, which we update here
                    subgraph->for_each_handle([&](handle_t sg_handle) {
                            // note, if components overlap, this is arbitrary.  up to user to only use
                            // path components if they are disjoint
                            node_to_component[subgraph->get_id(sg_handle)] = i;
                        });
                }
            }
        }
        // trace annotations
        if (trace) {
            // Even if we have only one chunk, the trace annotation data always
            // ends up in a file.
            string annot_name = chunk_name(out_chunk_prefix, i, output_regions[i], ".annotate.txt", 0, components);
            ofstream out_annot_file(annot_name);
            if (!out_annot_file) {
                cerr << "error[vg chunk]: can't open output trace annotation file " << annot_name << endl;
                exit(1);
            }
            for (auto tf : trace_thread_frequencies) {
                out_annot_file << tf.first << "\t" << tf.second << endl;
            }
        }
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
                 << "\t" << chunk_name(out_chunk_prefix, i, oregion, chunk_gam ? ".gam" : output_ext, 0, components);
            if (trace) {
                obed << "\t" << chunk_name(out_chunk_prefix, i, oregion, ".annotate.txt", 0, components);
            }
            obed << "\n";
        }
    }

    // write out component gams
    if (chunk_gam && components) {
        if(gam_is_gaf){
            cerr << "error[vg chunk]: GAF file input toggled with -F but, currently, only GAM files can by chunked by component. A workaround is to query one chromosome-component as the reference path and all contained snarls using '-p PATHNAME -S SNARLFILE'." << endl;
            return 1;
        }

        // buffer size of each component, total across threads
        static const size_t output_buffer_total_size = 100000;
        // split the buffer into threads
        size_t output_buffer_size = max((size_t)1, output_buffer_total_size / threads);
        // number of thread_buffers
        size_t num_buffers = num_regions * threads;
        
        vector<vector<Alignment>> output_buffers(num_buffers);
        vector<bool> append_buffer(num_regions, false);

        // protect our output buffers
        std::mutex* output_buffer_locks = new std::mutex[num_regions];

        // We may have too many components to keep a buffer open for each one.  So we open them as-needed only when flushing.
        function<void(size_t)> flush_gam_buffer = [&](size_t buffer_idx) {
            size_t comp_number = buffer_idx / threads;
            string gam_name = chunk_name(out_chunk_prefix, comp_number, output_regions[comp_number], ".gam", 0, components);
            {
                std::lock_guard<std::mutex> guard(output_buffer_locks[comp_number]);
                ofstream out_gam_file(gam_name, append_buffer[comp_number] ? std::ios_base::app : std::ios_base::out);
                if (!out_gam_file) {
                    cerr << "error[vg chunk]: can't open output gam file " << gam_name << endl;
                    exit(1);
                }
                vg::io::write_buffered(out_gam_file, output_buffers[buffer_idx], output_buffers[buffer_idx].size());
                append_buffer[comp_number] = true;
            }
            output_buffers[buffer_idx].clear();
        };
        
        function<void(Alignment&)> chunk_gam_callback = [&](Alignment& aln) {
            check_read(aln, graph);
             
            // we're going to lose unmapped reads right here
            if (aln.path().mapping_size() > 0) {
                nid_t aln_node_id = aln.path().mapping(0).position().node_id();
                unordered_map<nid_t, int32_t>::iterator comp_it = node_to_component.find(aln_node_id);                
                if (comp_it != node_to_component.end()) {
                    int32_t aln_component = comp_it->second;
                    size_t buffer_idx = aln_component * threads + omp_get_thread_num();
                    output_buffers[buffer_idx].push_back(aln);
                    if (output_buffers[buffer_idx].size() >= output_buffer_size) {
                        flush_gam_buffer(buffer_idx);
                    }
                }
            }
        };

        for (auto gam_file : gam_files) {
            get_input_file(gam_file, [&](istream& gam_stream) {
                    vg::io::for_each_parallel(gam_stream, chunk_gam_callback);
                });
        }
#pragma omp parallel for
        for (size_t buffer_idx = 0; buffer_idx < num_buffers; ++buffer_idx) {
            if (!output_buffers[buffer_idx].empty()) {
                flush_gam_buffer(buffer_idx);

            }
        }

        delete [] output_buffer_locks;
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_chunk("chunk", "split graph or alignment into chunks", main_chunk);

// Output name of a chunk
string chunk_name(const string& out_chunk_prefix, int i, const Region& region, string ext, int gi, bool components) {
    stringstream chunk_name;
    string seq = region.seq.empty() ? "ids" : region.seq;
    chunk_name << out_chunk_prefix;
    if (gi > 0) {
        chunk_name << "-" << gi;
    }
    if (!components) {
        chunk_name << "_" << i << "_" << seq << "_" << region.start << "_" << region.end;
    } else if (region.seq.empty()) {
        chunk_name << "_" << i;
    } else {
        chunk_name << "_" << region.seq;
    }
    chunk_name << ext;
    return chunk_name.str();
}

// Split out every chunk_size reads into a different file
int split_gam(istream& gam_stream, size_t chunk_size, const string& out_prefix, size_t gam_buffer_size) {
    ofstream out_file;
    size_t count = 0;
    
    // We're going to skip parsing the GAM reads and just treat these as type-tagged opaque messages.
    vg::io::MessageIterator gam_iterator(gam_stream);
    
    // We're going to do multithreaded output of batches of this size.
    // If our reads are paired, these had better be even.
    // We want to use a reasonably substantial size here instead of
    // gam_buffer_size which is pretty small, so we don't spin too much on the
    // OMP task machinery.
    size_t batch_size = std::min((size_t) 1000, chunk_size);
    
    // We use this to merge together compressed data from multiple threads.
    unique_ptr<vg::io::StreamMultiplexer> gam_multiplexer;
    
    // We need to know how many threads there will be.
    size_t thread_count = get_thread_count();
    
    // Each thread needs a place to keep a MessageEmitter
    vector<unique_ptr<vg::io::MessageEmitter>> emitters(thread_count);
    
    // We fill in batch buffers in the main thread and give them away to tasks to write.
    vector<vg::io::MessageIterator::TaggedMessage>* batch_in_progress = nullptr;
    
    #pragma omp parallel shared(gam_multiplexer, emitters)
    {
        #pragma omp single
        {
            while (gam_iterator.has_current()) {
                // There's a read message to process.
                if (count++ % chunk_size == 0) {
                    // We're at read 0, or the first read past the end of a chunk.
                    
                    // Wait for all tasks
                    #pragma omp taskwait
                    
                    for (size_t i = 0; i < thread_count; i++) {
                        // Flush everything tasks have written into the multiplexer
                        if (emitters[i]) {
                            emitters[i]->flush();
                            emitters[i].reset();
                            gam_multiplexer->register_breakpoint(i);
                        }
                    }
                            
                    if (out_file.is_open()) {
                        // Destroy the old multiplexer to flush
                        gam_multiplexer.reset();
                        // And close the file out.
                        out_file.close();
                    }
                    stringstream out_name;
                    out_name << out_prefix << setfill('0') <<setw(6) << (count / chunk_size + 1) << ".gam";
                    out_file.open(out_name.str());
                    if (!out_file) {
                        cerr << "error[vg chunk]: unable to open output gam: " << out_name.str() << endl;
                        exit(1);
                    }
                    // Open a new multiplexer on the new file
                    gam_multiplexer.reset(new vg::io::StreamMultiplexer(out_file, thread_count));
                }
                
                if (!batch_in_progress) {
                    // We need a new batch to fill in. Either we started a new file, or we shipped off our previous batch.
                    batch_in_progress = new vector<vg::io::MessageIterator::TaggedMessage>();
                    batch_in_progress->reserve(batch_size);
                }
                
                
                // Grab the message, paired with its tag, and advance.
                // TODO: Stop copying tag?
                batch_in_progress->emplace_back(std::move(gam_iterator.take()));
                if (batch_in_progress->back().first.empty()) {
                    // This is untagged data; assume it's GAM.
                    batch_in_progress->back().first = "GAM";
                }
                if (!batch_in_progress->back().second) {
                    // This is just a tag alone; throw this away.
                    batch_in_progress->pop_back();
                    count--;
                }
                
                if (batch_in_progress->size() == batch_size || count % chunk_size == 0 || !gam_iterator.has_current()) {
                    // We've hit the batch size, or we've hit the last read that fits in this chunk, or we've hit the end of the input.
                    // Launch a task to deal with this batch
                    #pragma omp task firstprivate(batch_in_progress)
                    {
                        // Get our thread
                        size_t thread = omp_get_thread_num();
                        
                        // Find our GAM emitter
                        auto& emitter_ptr = emitters[thread];
                        if (!emitter_ptr) {
                            // Make it exist if it doesn't yet. Make sure to compress and to pass along the buffer size.
                            emitter_ptr.reset(new vg::io::MessageEmitter(gam_multiplexer->get_thread_stream(thread), true, gam_buffer_size));
                        }
                    
                        for (auto& message : *batch_in_progress) {
                            // Send each message over to the emitter, with the tag, moving the message body
                            emitter_ptr->write(message.first, std::move(*message.second));
                        }
                        
                        // Throw out our assigned batch copy.
                        delete batch_in_progress;
                        
                        if (gam_multiplexer->want_breakpoint(thread)) {
                            // The multiplexer wants our data.
                            // Flush and create a breakpoint.
                            emitter_ptr->flush();
                            gam_multiplexer->register_breakpoint(thread);
                            // TODO: No EOF marker we could remove???
                        }
                    }
                    
                    
                    // Task frees the batch, so null out our pointer to it.
                    batch_in_progress = nullptr;
                }
            }
            
            // Wait for the final tasks.
            #pragma omp taskwait
            
            for (size_t i = 0; i < thread_count; i++) {
                // Flush everything tasks have written into the multiplexer
                if (emitters[i]) {
                    emitters[i]->flush();
                    emitters[i].reset();
                    gam_multiplexer->register_breakpoint(i);
                }
            }
            
            // Get rid of the multiplexer to flush
            gam_multiplexer.reset();
            
            // There will be no final batch, because when we hit EOF we launch a task
            assert(batch_in_progress == nullptr);
        }
    }
    
    if (out_file.is_open()) {
        // Close out the file.
        out_file.close();
    }
    return 0;
}

/// Stop and print an error if the graph exists and the read does not appear to
/// actually be aligned against the graph.
static void check_read(const Alignment& aln, const HandleGraph* graph) { 
    if (!graph) {
        return;
    }
    // Make sure the nodes it visits could be the nodes in the graph.
    AlignmentValidity validity = alignment_is_valid(aln, graph);
    if (!validity) {
        #pragma omp critical (cerr)
        {
            std::cerr << "error:[vg chunk] Alignment " << aln.name() << " cannot be interpreted against this graph: " << validity.message << std::endl;
            std::cerr << "Make sure that you are using the same graph that the reads were mapped to!" << std::endl;
        }
        exit(1);
    }
}

