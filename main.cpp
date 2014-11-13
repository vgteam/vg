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
         << "  -- view          conversion (protobuf/json/GFA)" << endl
         << "  -- align         alignment" << endl;
}

void align_help(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg" << endl
        //<< "    -p, --print-cigar     output graph cigar for alignments" << endl
         << "    -j, --json            output alignments in JSON format (default)" << endl;
}

void view_help(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -d, --dot             output dot format (default)" << endl
         << "    -g, --gfa             output GFA format" << endl
         << "    -j, --json            output VG JSON format" << endl;
}

void construct_help(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options]" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -p, --protobuf        output VG protobuf format (default)" << endl
         << "    -g, --gfa             output GFA format" << endl
         << "    -j, --json            output VG JSON format" << endl;
}


int align_main(int argc, char** argv) {

    string seq;

    if (argc == 2) {
        align_help(argv);
        return 1;
    }

    bool print_cigar = false;
    bool output_json = true;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:j",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 's':
            seq = optarg;
            break;

            /*
        case 'p':
            print_cigar = true;
            break;
            */

        case 'j':
            output_json = true;
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

    VariantGraph* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VariantGraph(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VariantGraph(in);
    }

    graph->create_alignable_graph();
    Alignment alignment = graph->align(seq);

    char *json2 = pb2json(alignment);
    cout<<json2<<endl;
    free(json2);

    delete graph;

    return 0; // not implemented

}

int view_main(int argc, char** argv) {

    if (argc == 2) {
        view_help(argv);
        return 1;
    }

    string output_type = "dot";

    int c;
    optind = 2; // force optind past "view" argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"dot", no_argument, 0, 'd'},
                {"gfa", no_argument, 0, 'g'},
                {"json",  no_argument, 0, 'j'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgj",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            output_type = "dot";
            break;
 
        case 'g':
            output_type = "GFA";
            break;
 
        case 'j':
            output_type = "JSON";
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

    VariantGraph* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VariantGraph(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VariantGraph(in);
    }

    if (output_type == "dot") {
        graph->to_dot(std::cout);
    } else if (output_type == "JSON") {
        char *json2 = pb2json(graph->graph);
        cout<<json2<<endl;
        free(json2);
    }

    delete graph;

    return 0; // not implemented
}

int construct_main(int argc, char** argv) {

    if (argc == 2) {
        construct_help(argv);
        return 1;
    }

    string fasta_file_name, vcf_file_name;
    string output_type = "VG";

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"protobuf",  no_argument, 0, 'p'},
                {"gfa",  no_argument, 0, 'g'},
                {"json",  no_argument, 0, 'j'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:pgj",
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

        case 'p':
            output_type = "VG";
            break;
 
        case 'g':
            output_type = "GFA";
            break;
 
        case 'j':
            output_type = "JSON";
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

    if (output_type == "VG") {
        //ofstream of("test.vg");
        graph.graph.SerializeToOstream(&std::cout);
    } else if (output_type == "JSON") {
        char *json2 = pb2json(graph.graph);
        cout<<json2<<endl;
        free(json2);
    }

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
    } else if (command == "align") {
        return align_main(argc, argv);
    }

    return 0;

}
