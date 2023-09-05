/** \file translate_main.cpp
 *
 * Defines the "vg translate" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../translator.hpp"
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_translate(char** argv) {
    cerr << "usage: " << argv[0] << " translate [options] translation" << endl
         << "Translate alignments or paths using the translation map." << endl
         << endl
         << "options:" << endl
         << "    -p, --paths FILE      project the input paths into the from-graph" << endl
         << "    -a, --alns FILE       project the input alignments into the from-graph" << endl
         << "    -l, --loci FILE       project the input locus descriptions into the from-graph" << endl
         << "    -m, --mapping JSON    print the from-mapping corresponding to the given JSON mapping" << endl
         << "    -P, --position JSON   print the from-position corresponding to the given JSON position" << endl
         << "    -o, --overlay FILE    overlay this translation on top of the one we are given" << endl;
}

int main_translate(int argc, char** argv) {

    if (argc <= 2) {
        help_translate(argv);
        return 1;
    }

    string position_string;
    string mapping_string;
    string path_file;
    string aln_file;
    string loci_file;
    string overlay_file;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"position", required_argument, 0, 'P'},
            {"mapping", required_argument, 0, 'm'},
            {"paths", required_argument, 0, 'p'},
            {"alns", required_argument, 0, 'a'},
            {"loci", required_argument, 0, 'l'},
            {"overlay", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:m:P:a:o:l:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'P':
            position_string = optarg;
            break;

        case 'm':
            mapping_string = optarg;
            break;

        case 'p':
            path_file = optarg;
            break;

        case 'a':
            aln_file = optarg;
            break;

        case 'l':
            loci_file = optarg;
            break;

        case 'o':
            overlay_file = optarg;
            break;

        case 'h':
        case '?':
            help_translate(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    Translator* translator;
    get_input_file(optind, argc, argv, [&](istream& in) {
        translator = new Translator(in);
    });

    // test the position translation
    if (!position_string.empty()) {
        Position position;
        json2pb(position, position_string.c_str(), position_string.size());
        cout << pb2json(translator->translate(position)) << endl;
    }

    // test the mapping translation
    if (!mapping_string.empty()) {
        Mapping mapping;
        json2pb(mapping, mapping_string.c_str(), mapping_string.size());
        cout << pb2json(translator->translate(mapping)) << endl;
    }

    if (!path_file.empty()) {
        vector<Path> buffer;
        function<void(Path&)> lambda = [&](Path& path) {
            buffer.push_back(translator->translate(path));
            vg::io::write_buffered(cout, buffer, 100);
        };
        ifstream path_in(path_file);
        vg::io::for_each(path_in, lambda);
        vg::io::write_buffered(cout, buffer, 0);
    } else if (!aln_file.empty()) {
        vector<Alignment> buffer;
        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            buffer.push_back(translator->translate(aln));
            vg::io::write_buffered(cout, buffer, 100);
        };
        ifstream aln_in(aln_file);
        vg::io::for_each(aln_in, lambda);
        vg::io::write_buffered(cout, buffer, 0);
    } else if (!loci_file.empty()) {
        vector<Locus> buffer;
        function<void(Locus&)> lambda = [&](Locus& locus) {
            buffer.push_back(translator->translate(locus));
            vg::io::write_buffered(cout, buffer, 100);
        };
        ifstream loci_in(loci_file);
        vg::io::for_each(loci_in, lambda);
        vg::io::write_buffered(cout, buffer, 0);
    }

    if (!overlay_file.empty()) {
        vector<Translation> buffer;
        function<void(Translation&)> lambda = [&](Translation& trans) {
            buffer.push_back(translator->overlay(trans));
            vg::io::write_buffered(cout, buffer, 100);
        };
        ifstream overlay_in(overlay_file);
        vg::io::for_each(overlay_in, lambda);
        vg::io::write_buffered(cout, buffer, 0);
    }

    return 0;
}


// Register subcommand
static Subcommand vg_version("translate", "project alignments and paths through a graph translation", DEPRECATED, main_translate);

