#include "../gamsorter.hpp"
#include "../gam_index.hpp"
#include "../stream.hpp"
#include <getopt.h>
#include "subcommand.hpp"
#include "../index.hpp"

/**
* GAM sort main
*/

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_gamsort(char **argv)
{
    cerr << "gamsort: sort a GAM file, or index a sorted GAM file" << endl
         << "Usage: " << argv[1] << " [Options] gamfile" << endl
         << "Options:" << endl
         << "  -s / --sorted           Input GAM is already sorted." << endl
         << "  -i / --index FILE       produce an index of the sorted GAM file" << endl
         << "  -d / --dumb-sort        use naive sorting algorithm (no tmp files, faster for small GAMs)" << endl
         << "  -r / --rocks            Just use the old RocksDB-style indexing scheme for sorting." << endl
         << "  -a / --aln-index        Create the old RocksDB-style node-to-alignment index." << endl
         << "  -p / --progress         Show progress." << endl
         << "  -t / --threads          Use the specified number of threads." << endl
         << endl;
}

int main_gamsort(int argc, char **argv)
{
    string index_filename;
    bool dumb_sort = false;
    bool is_sorted = false;
    bool just_use_rocks = false;
    bool do_aln_index = false;
    bool show_progress = false;
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
                {"rocks", no_argument, 0, 'r'},
                {"aln-index", no_argument, 0, 'a'},
                {"is-sorted", no_argument, 0, 's'},
                {"progress", no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "i:dhraspt:",
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
            dumb_sort = true;
            break;
        case 's':
            is_sorted = true;
            break;
        case 'r':
            just_use_rocks = true;
            break;
        case 'a':
            do_aln_index = true;
            break;
        case 'p':
            show_progress = true;
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

    get_input_file(optind, argc, argv, [&](istream& gam_in) {

        GAMSorter gs(show_progress);

        if (just_use_rocks) {
            // Do the sort the old way - write a big ol'
            // RocksDB index of alignments, then dump them
            // from that DB. Loses unmapped reads.
            
            if (index_filename.empty()) {
                cerr << "error:[vg gamsort]: Cannot index using RocksDB without an index filename (-i)" << endl;
                exit(1);
            }
            
            string dbname = index_filename;
            Index index;

            // Index the alignments in RocksDB
            index.open_for_bulk_load(dbname);
            int64_t aln_idx = 0;
            function<void(Alignment&)> lambda_reader = [&index](Alignment& aln) {
                    index.put_alignment(aln);
            };
            stream::for_each_parallel(gam_in, lambda_reader);
            
            // Print them out again in order
            vector<Alignment> output_buf;
            auto lambda_writer = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 1000);
            };
            index.for_each_alignment(lambda_writer);
            stream::write_buffered(cout, output_buf, 0);
            
            index.flush();
            index.close();
        } else {
            // Do a normal GAMSorter sort
            unique_ptr<GAMIndex> index;
            
            if (!index_filename.empty()) {
                // Make an index
                index = unique_ptr<GAMIndex>(new GAMIndex());
            }
            
            if (dumb_sort) {
                // Sort in a single pass in memory
                gs.dumb_sort(gam_in, cout, index.get());
            } else {
                // Sort using fan-in-limited temp file merging 
                gs.stream_sort(gam_in, cout, index.get());
            }
            
            if (index.get() != nullptr) {
                // Save the index
                ofstream index_out(index_filename);
                index->save(index_out);
            }
        }
    });

    return 0;
}

static Subcommand vg_gamsort("gamsort", "Sort a GAM file or index a sorted GAM file.", main_gamsort);
