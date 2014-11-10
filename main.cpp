#include <iostream>
#include <fstream>
#include <getopt.h>
#include "pb2json.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"
#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

void main_help(char** argv) {
    cerr << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl 
         << "  -- construct     graph construction" << endl
         << "  -- view          conversion (protobuf/json/GFA)" << endl;
}

void view_help(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] <file>" << endl;
}

void construct_help(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options]" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -P, --protobuf        output VG protobuf format (default)" << endl
         << "    -G, --gfa             output GFA format" << endl
         << "    -J, --json            output VG JSON format" << endl;
}

int view_main(int argc, char** argv) {
    view_help(argv);
    return 0; // not implemented
}

int construct_main(int argc, char** argv) {

    if (argc == 2) {
        construct_help(argv);
        return 1;
    }

    string fasta_file_name, vcf_file_name;
    bool output_VG=true, output_GFA=false, output_JSON=false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"protobuf",  no_argument, 0, 'P'},
                {"gfa",  no_argument, 0, 'G'},
                {"json",  no_argument, 0, 'J'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:PGJ",
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
            output_VG = true;
            break;
 
        case 'G':
            output_GFA = true;
            break;
 
        case 'J':
            output_JSON = true;
            break;
 
        case '?':
            /* getopt_long already printed an error message. */
            construct_help(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    // check that we are only outputting a single stream type
    bool single_output_type = true;
    if (output_VG) {
        if (output_GFA || output_JSON) {
            single_output_type = false;
        }
    } else if (output_GFA) {
        if (output_VG || output_JSON) {
            single_output_type = false;
        }
    } else if (output_JSON) {
        if (output_VG || output_GFA) {
            single_output_type = false;
        }
    }

    if (!single_output_type) {
        cerr << "error:[vg construct] can only output a single graph format (one of VG, GFA, JSON)" << endl;
        return 1;
    }

    // set up our inputs

    vcf::VariantCallFile variant_file;
    if (vcf_file_name.empty()) {
        cerr << "error:[vg construct] a VCF file is required for graph construction" << endl;
        return 1;
    }
    variant_file.open(vcf_file_name);
    if (!variant_file.is_open()) {
        cerr << "error:[vg construct] could not open" << vcf_file_name << endl;
        return 1;
    }

    FastaReference reference;
    if (fasta_file_name.empty()) {
        cerr << "error:[vg construct] a reference is required for graph construction" << endl;
        return 1;
    }
    reference.open(fasta_file_name);

    VariantGraph graph(variant_file, reference);

    char *json2 = pb2json(graph);
	cout<<json2<<endl;
	free(json2);

    return 0;
}

int main(int argc, char *argv[])
{

    if (argc == 1) {
        main_help(argv);
        return 1;
    }

    string command = argv[1];
    if (command == "construct") {
        return construct_main(argc, argv);
    } else if (command == "view") {
        return view_main(argc, argv);
    }

    return 0;

}
