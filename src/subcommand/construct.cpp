// construct.cpp: define the "vg construct" subcommand.

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../stream.hpp"
#include "../constructor.hpp"
#include "../region.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -a, --alt-paths       save paths for alts of variants by variant ID" << endl
         << "    -R, --region REGION   specify a particular chromosome or 1-based inclusive region" << endl
         << "    -C, --region-is-chrom don't attempt to parse the region (use when the reference" << endl
         << "                          sequence name could be inadvertently parsed as a region)" << endl
         << "    -z, --region-size N   variants per region to parallelize" << endl
         << "    -m, --node-max N      limit the maximum allowable node sequence size (defaults to 1000)" << endl
         << "                          nodes greater than this threshold will be divided" << endl
         << "                          Note: nodes larger than ~1024 bp can't be GCSA2-indexed" << endl
         << "    -p, --progress        show progress" << endl
         << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl
         << "    -f, --flat-alts N     don't chop up alternate alleles from input vcf" << endl;
}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    // Make a constructor to fill in
    Constructor constructor;

    // We also parse some arguments separately.
    string fasta_file_name, vcf_file_name, json_filename;
    string region;
    bool region_is_chrom = false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"alt-paths", no_argument, 0, 'a'},
                {"progress",  no_argument, 0, 'p'},
                {"region-size", required_argument, 0, 'z'},
                {"threads", required_argument, 0, 't'},
                {"region", required_argument, 0, 'R'},
                {"region-is-chrom", no_argument, 0, 'C'},
                {"node-max", required_argument, 0, 'm'},\
                {"flat-alts", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:phz:t:R:m:as:Cf",
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

        case 'a':
            constructor.alt_paths = true;
            break;

        case 'p':
            constructor.show_progress = true;
            break;

        case 'z':
            constructor.vars_per_chunk = atoi(optarg);
            break;

        case 'R':
            region = optarg;
            break;

        case 'C':
            region_is_chrom = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'm':
            constructor.max_node_size = atoi(optarg);
            break;

        case 'f':
            constructor.flat = true;
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
    
    if (constructor.max_node_size == 0) {
        // Make sure we can actually make nodes
        cerr << "error:[vg construct] max node size cannot be 0" << endl;
        exit(1);
    }
    
    if (!region.empty()) {
        // We want to limit to a certain region
        if (!region_is_chrom) {
            // We are allowed to parse the region.
            // Break out sequence name and region bounds
            string seq_name;
            int start_pos = 0, stop_pos = 0;
            parse_region(region,
                         seq_name,
                         start_pos,
                         stop_pos);
                         
            if (start_pos != 0 && stop_pos != 0) {
                // These are 0-based, so if both are nonzero we got a real set of coordinates
                if (constructor.show_progress) {
                    cerr << "Restricting to " << seq_name << " from " << start_pos << " to " << stop_pos << endl;
                }
                constructor.allowed_vcf_names.insert(seq_name);
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                constructor.allowed_vcf_regions[seq_name] = make_pair(start_pos - 1, stop_pos);
            } else if (start_pos == 0 && stop_pos == 0) {
                // We just got a name
                cerr << "Restricting to " << seq_name << " from 1 to end" << endl;
                constructor.allowed_vcf_names.insert(seq_name);
            } else {
                // This doesn't make sense. Does it have like one coordinate?
                cerr << "error:[vg construct] could not parse " << region << endl;
                exit(1);
            }
         } else {
            // We have been told not to parse the region
            cerr << "Restricting to " << region << " from 1 to end" << endl;
            constructor.allowed_vcf_names.insert(region);
         }
    }

    vcflib::VariantCallFile variant_file;
    if (!vcf_file_name.empty()) {
        // Make sure the file exists. Otherwise Tabix++ may exit with a non-
        // helpful message.

        // We can't invoke stat woithout a place for it to write. But all we
        // really want is its return value.
        struct stat temp;
        if(stat(vcf_file_name.c_str(), &temp)) {
            cerr << "error:[vg construct] file \"" << vcf_file_name << "\" not found" << endl;
            return 1;
        }
        variant_file.open(vcf_file_name);
        if (!variant_file.is_open()) {
            cerr << "error:[vg construct] could not open" << vcf_file_name << endl;
            return 1;
        }
    }

    if (fasta_file_name.empty()) {
        cerr << "error:[vg construct] a reference is required for graph construction" << endl;
        return 1;
    }
    FastaReference reference;
    reference.open(fasta_file_name);

    // We need a callback to handle pieces of graph as they are produced.
    auto callback = [&](Graph& big_chunk) {
        // TODO: these chunks may be too big to (de)serialize directly. For now,
        // just serialize them directly anyway.
        #pragma omp critical (cout)
        stream::write(cout, 1, std::function<Graph(uint64_t)>([&](uint64_t chunk_number) -> Graph {
            assert(chunk_number == 0);
            // Just spit out our one chunk
            return big_chunk;
        }));
    };

    // Construct the graph. Make sure to make our FASTA and VCF into vectors.
    constructor.construct_graph({&reference}, {&variant_file}, callback);

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("construct", "graph construction", main_construct);

