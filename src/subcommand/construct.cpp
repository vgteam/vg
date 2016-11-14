// construct.cpp: define the "vg construct" subcommand.

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -P, --ref-paths FILE  write reference paths in protobuf/gzip format to FILE" << endl
         << "    -B, --phase-blocks    save paths for phased blocks with the ref paths" << endl
         << "    -a, --alt-paths       save paths for alts of variants by variant ID" << endl
         << "    -R, --region REGION   specify a particular chromosome" << endl
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

    string fasta_file_name, vcf_file_name, json_filename;
    string region;
    bool region_is_chrom = false;
    string output_type = "VG";
    bool progress = false;
    int vars_per_region = 25000;
    int max_node_size = 1000;
    string ref_paths_file;
    bool flat_alts = false;
    // Should we make paths out of phasing blocks in the called samples?
    bool load_phasing_paths = false;
    // Should we make alt paths for variants?
    bool load_alt_paths = false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                // TODO: change the long option here?
                {"ref-paths", required_argument, 0, 'P'},
                {"phase-blocks", no_argument, 0, 'B'},
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
        c = getopt_long (argc, argv, "v:r:phz:t:R:m:P:Bas:Cf",
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

        case 'P':
            ref_paths_file = optarg;
            break;

        case 'B':
            load_phasing_paths = true;
            break;

        case 'a':
            load_alt_paths = true;
            break;

        case 'p':
            progress = true;
            break;

        case 'z':
            vars_per_region = atoi(optarg);
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
            max_node_size = atoi(optarg);
            break;

        case 'f':
            flat_alts = true;
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

    if(load_phasing_paths && ref_paths_file.empty()) {
        cerr << "error:[vg construct] cannot save phasing paths without a paths file name" << endl;
        return 1;
    }

    if (fasta_file_name.empty()) {
        cerr << "error:[vg construct] a reference is required for graph construction" << endl;
        return 1;
    }
    FastaReference reference;
    reference.open(fasta_file_name);

    // store our reference sequence paths
    // TODO: use this. Maybe dump paths here instead of in the graph?
    Paths ref_paths;

    VG graph(variant_file, reference, region, region_is_chrom, vars_per_region,
             max_node_size, flat_alts, load_phasing_paths, load_alt_paths, progress);

    if (!ref_paths_file.empty()) {
        ofstream paths_out(ref_paths_file);
        graph.paths.write(paths_out);
        if(load_phasing_paths) {
            // Keep only the non-phasing paths in the graph. If you keep too
            // many paths in a graph, you'll make chunks that are too large.
            // TODO: dynamically deliniate the chunks in the serializer so you
            // won't write vg files you can't read.

            set<string> non_phase_paths;
            string phase_prefix = "_phase";
            graph.paths.for_each_name([&](string path_name) {
                    if(!equal(phase_prefix.begin(), phase_prefix.end(), path_name.begin())) {
                        // Path is not a phase path
                        non_phase_paths.insert(path_name);
                    }
                });

            // Keep only the non-phase paths
            graph.paths.keep_paths(non_phase_paths);
        }
    }

    graph.serialize_to_ostream(std::cout);

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("construct", "graph construction", main_construct);

