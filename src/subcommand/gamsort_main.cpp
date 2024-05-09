#include "../stream_sorter.hpp"
#include <vg/io/stream.hpp>
#include "../stream_index.hpp"
#include <getopt.h>
#include "subcommand.hpp"
#include "vg/io/gafkluge.hpp"
#include "alignment.hpp"


/**
* GAM sort main
*/

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_gamsort(char **argv)
{
    cerr << "gamsort: sort a GAM/GAF file, or index a sorted GAM file" << endl
         << "Usage: " << argv[1] << " [Options] gamfile" << endl
         << "Options:" << endl
         << "  -i / --index FILE       produce an index of the sorted GAM file" << endl
         << "  -d / --dumb-sort        use naive sorting algorithm (no tmp files, faster for small GAMs)" << endl
         << "  -s / --shuffle          Shuffle reads by hash (GAM only)" << endl
         << "  -p / --progress         Show progress." << endl
         << "  -G / --gaf-input        Input is a GAF file." << endl
         << "  -c / --chunk-size       Number of reads per chunk when sorting GAFs." << endl
         << "  -t / --threads          Use the specified number of threads." << endl
         << endl;
}

// defines how to compare two GAF records
// first using 'rk1' tag (here, minimum node ID). If tied, use 'rk2' tag (here, maximum node ID)
struct compare_gaf {
    bool operator()(const gafkluge::GafRecord& gaf1, const gafkluge::GafRecord& gaf2) {
        // TODO find a way to not have to convert the node ids to string before and then back to int here?
        long long rk11 = std::stoll(gaf1.opt_fields.find("rk1")->second.second);
        long long rk12 = std::stoll(gaf2.opt_fields.find("rk1")->second.second);
        long long rk21 = std::stoll(gaf1.opt_fields.find("rk2")->second.second);
        long long rk22 = std::stoll(gaf2.opt_fields.find("rk2")->second.second);
        return rk11 < rk12 || (rk11 == rk12 && rk21 < rk22);
    }
};
// defines a pair of a GAF record and the ID of the file it came from (used when merging sorted GAF files)
struct GafFile {
    gafkluge::GafRecord gaf;
    int file_i;
};
// comparator used by the min-heap when merging sorted GAF files
struct greater_gaffile {
    bool operator()(const GafFile& gf1, const GafFile& gf2) {
        // TODO find a way to not have to convert the node ids to string before and then back to int here?
        long long rk11 = std::stoll(gf1.gaf.opt_fields.find("rk1")->second.second);
        long long rk12 = std::stoll(gf2.gaf.opt_fields.find("rk1")->second.second);
        long long rk21 = std::stoll(gf1.gaf.opt_fields.find("rk2")->second.second);
        long long rk22 = std::stoll(gf2.gaf.opt_fields.find("rk2")->second.second);
        return rk11 > rk12 || (rk11 == rk12 && rk21 > rk22);
    }
};

int main_gamsort(int argc, char **argv)
{
    string index_filename;
    bool easy_sort = false;
    bool shuffle = false;
    bool show_progress = false;
    string input_format = "GAM";
    int chunk_size = 1000000; // maximum number reads held in memory
    // We limit the max threads, and only allow thread count to be lowered, to
    // prevent tcmalloc from giving each thread a very large heap for many
    // threads.
    // On my machine we can keep about 4 threads busy.
    size_t num_threads = 4;
    int c;
    optind = 2; // force optind past command positional argument
    while (true)
    {
        static struct option long_options[] =
            {
                {"index", required_argument, 0, 'i'},
                {"dumb-sort", no_argument, 0, 'd'},
                {"shuffle", no_argument, 0, 's'},
                {"progress", no_argument, 0, 'p'},
                {"gaf-input", no_argument, 0, 'g'},
                {"chunk-size", required_argument, 0, 'c'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "i:dshpGt:c:",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            index_filename = optarg;
            break;
        case 'd':
            easy_sort = true;
            break;
        case 's':
            shuffle = true;
            break;
        case 'p':
            show_progress = true;
            break;
        case 'G':
            input_format = "GAF";
            break;
        case 'c':
            chunk_size = parse<int>(optarg);
            break;
        case 't':
            num_threads = min(parse<size_t>(optarg), num_threads);
            break;
        case 'h':
        case '?':
        default:
            help_gamsort(argv);
            exit(1);
        }
    }

    if (argc < 3){
        help_gamsort(argv);
        exit(1);
    }
    
    omp_set_num_threads(num_threads);

    if (input_format == "GAM") {
        if (shuffle && !index_filename.empty()) {
            cerr << "[vg gamsort] Indexing is not allowed when shuffling GAM files." << endl;
            exit(1);
        }
        get_input_file(optind, argc, argv, [&](istream& gam_in) {

            GAMSorter gs(shuffle ? GAMSorter::Order::RANDOM : GAMSorter::Order::BY_GRAPH_POSITION, show_progress);

            // Do a normal GAMSorter sort
            unique_ptr<GAMIndex> index;
        
            if (!index_filename.empty()) {
                // Make an index
                index = unique_ptr<GAMIndex>(new GAMIndex());
            }
        
            if (easy_sort) {
                // Sort in a single pass in memory
                gs.easy_sort(gam_in, cout, index.get());
            } else {
                // Sort using fan-in-limited temp file merging
                gs.stream_sort(gam_in, cout, index.get());
            }
        
            if (index.get() != nullptr) {
                // Save the index
                ofstream index_out(index_filename);
                index->save(index_out);
            }
        });
    } else if (input_format == "GAF") {
        if (shuffle) {
            // TODO: Implement shuffling for GAF files by making the
            // comparators switch modes and hashing the record strings.
            // TODO: Is there a way to be less duplicative with the
            // StreamSorter?
            cerr << "[vg gamsort] Shuffling is not implemented for GAF files." << endl;
            exit(1);
        }

        std::string input_gaf_filename = get_input_file_name(optind, argc, argv);

        // where to store the chunk of GAF records that will be sorted, then written to disk,
        // (then later merged with the other sorted chunks)
        std::vector<gafkluge::GafRecord> current_gaf_chunk;
        int count = 0; // read count
        int chunk_id = 0; // ID of the current chunk
        std::vector<std::string> chunk_files; // names of the chunk files
        
        // read input GAF file
        htsFile* in = hts_open(input_gaf_filename.c_str(), "r");
        if (in == NULL) {
            cerr << "[vg gamsort] couldn't open " << input_gaf_filename << endl; exit(1);
        }
        kstring_t s_buffer = KS_INITIALIZE;
        gafkluge::GafRecord gaf;

        string chunk_outf = temp_file::create();
        if(show_progress){
            cerr << "Preparing temporary chunk " << chunk_outf << "..." << endl;
        }
        
        while (vg::io::get_next_record_from_gaf(nullptr, nullptr, in, s_buffer, gaf) == true) {
            // find the minimum and maximum node IDs
            nid_t min_node = std::numeric_limits<nid_t>::max();
            nid_t max_node = 0;
            for (size_t i = 0; i < gaf.path.size(); ++i) {
                const auto& gaf_step = gaf.path[i];
                assert(gaf_step.is_stable == false);
                assert(gaf_step.is_interval == false);
                nid_t nodeid = std::stol(gaf_step.name);
                if (min_node > nodeid){
                    min_node = nodeid;
                }
                if (max_node < nodeid){
                    max_node = nodeid;
                }
            }
            // write them as new GAF tags 'rk1' and 'rk2'
            // they'll get written in the temporary chunks to avoid having
            // to find them again when merging them
            gaf.opt_fields["rk1"] = make_pair('i', std::to_string(min_node));
            gaf.opt_fields["rk2"] = make_pair('i', std::to_string(max_node));
            current_gaf_chunk.push_back(gaf);
            count++;

            // if we've read enough reads, sort them and write to disk
            if(count == chunk_size){
                // sort by minimum node id
                if(show_progress){
                    cerr << "   Sorting chunk..." << endl;
                }
                std::stable_sort(current_gaf_chunk.begin(), current_gaf_chunk.end(), compare_gaf());
                // write to temporary file
                if(show_progress){
                    cerr << "   Writing chunk..." << endl;
                }
                std::ofstream out_file(chunk_outf);
                for (int ii=0; ii<current_gaf_chunk.size(); ii++){
                    out_file << current_gaf_chunk[ii] << endl;
                } 
                out_file.close();
                chunk_files.push_back(chunk_outf);
                // init next chunk
                current_gaf_chunk.clear();
                count = 0;
                chunk_id++;
                chunk_outf = temp_file::create();
                if(show_progress){
                    cerr << "Preparing temporary chunk " << chunk_outf << "..." << endl;
                }
            } 
        }
        hts_close(in);

        // write the current last chunk too, if it has any reads
        if(count > 0){
            // sort by minimum node id
            if(show_progress){
                cerr << "   Sorting chunk..." << endl;
            }
            std::stable_sort(current_gaf_chunk.begin(), current_gaf_chunk.end(), compare_gaf());
            // write to temporary file
            if(show_progress){
                cerr << "   Writing chunk..." << endl;
            }
            std::ofstream out_file(chunk_outf);
            for (int ii=0; ii<current_gaf_chunk.size(); ii++){
                out_file << current_gaf_chunk[ii] << endl;
            } 
            out_file.close();
            chunk_files.push_back(chunk_outf);
        } 

        // merge the chunks of sorted reads
        // open all the files TODO don't do that if too many files
        if(show_progress){
            cerr << "Merging " << chunk_files.size() << " files..." << endl;
        }
        
        std::vector<htsFile*> opened_files;
        std::vector<bool> more_in_file;
        std::vector<kstring_t> opened_file_buffers;
        // heap with the current GAF record of each file
        std::priority_queue<GafFile, std::vector<GafFile>, greater_gaffile > opened_records;

        std::string line;

        // open the temp GAF files and read the first record
        GafFile gf;
        for(int ii=0; ii < chunk_files.size(); ii++){
            htsFile* in = hts_open(chunk_files[ii].c_str(), "r");
            if (in == NULL) {
                cerr << "[vg::alignment.cpp] couldn't open " << input_gaf_filename << endl; exit(1);
            }
            opened_file_buffers.push_back(KS_INITIALIZE);
            opened_files.push_back(in);
            if(vg::io::get_next_record_from_gaf(nullptr, nullptr, opened_files.back(), opened_file_buffers.back(), gaf)){
                gf.gaf = gaf;
                gf.file_i = ii;
                opened_records.push(gf);
            }
        }

        while(opened_records.size() > 0){
            // which file will have the smallest record (i.e. to output first)
            gf = opened_records.top();
            // output smallest record
            cout << gf.gaf << endl;
            opened_records.pop();
            if(vg::io::get_next_record_from_gaf(nullptr, nullptr, opened_files[gf.file_i], opened_file_buffers[gf.file_i], gf.gaf)){
                opened_records.push(gf);
            }
        }

    }
    return 0;
}

static Subcommand vg_gamsort("gamsort", "Sort a GAM/GAF file or index a sorted GAM file.", main_gamsort);
