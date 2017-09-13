#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <sys/stat.h>
#include "gcsa/gcsa.h"
#include "gcsa/algorithms.h"
#include "json2pb.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "vg_set.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "stream.hpp"
#include "alignment.hpp"
#include "convert.hpp"
#include "caller.hpp"
#include "deconstructor.hpp"
#include "filter.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"
#include "version.hpp"
#include "genotyper.hpp"
#include "bubbles.hpp"
#include "readfilter.hpp"
#include "distributions.hpp"
#include "unittest/driver.hpp"
// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"
#include "flow_sort.hpp"


using namespace std;
using namespace google::protobuf;
using namespace vg;

    void help_sv(char** argv){
        cerr << "usage: " << argv[0] << " sv [options] <aln.gam>" << endl
            << "options: " << endl
            << " -g --graph <graph>.vg " << endl
            << " -m --mask <vcf>.vcf" << endl
            << endl;
    }

void help_locify(char** argv){
    cerr << "usage: " << argv[0] << " locify [options] " << endl
         << "    -l, --loci FILE      input loci over which to locify the alignments" << endl
         << "    -a, --aln-idx DIR    use this rocksdb alignment index (from vg index -N)" << endl
         << "    -x, --xg-idx FILE    use this xg index" << endl
         << "    -n, --name-alleles   generate names for each allele rather than using full Paths" << endl
         << "    -f, --forwardize     flip alignments on the reverse strand to the forward" << endl
         << "    -s, --sorted-loci FILE  write the non-nested loci out in their sorted order" << endl
         << "    -b, --n-best N       keep only the N-best alleles by alignment support" << endl
         << "    -o, --out-loci FILE  rewrite the loci with only N-best alleles kept" << endl;
        // TODO -- add some basic filters that are useful downstream in whatshap
}

int main_locify(int argc, char** argv){
    string gam_idx_name;
    string loci_file;
    Index gam_idx;
    string xg_idx_name;
    bool name_alleles = false;
    bool forwardize = false;
    string loci_out, sorted_loci;
    int n_best = 0;

    if (argc <= 2){
        help_locify(argv);
        exit(1);
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"gam-idx", required_argument, 0, 'g'},
            {"loci", required_argument, 0, 'l'},
            {"xg-idx", required_argument, 0, 'x'},
            {"name-alleles", no_argument, 0, 'n'},
            {"forwardize", no_argument, 0, 'f'},
            {"sorted-loci", required_argument, 0, 's'},
            {"loci-out", required_argument, 0, 'o'},
            {"n-best", required_argument, 0, 'b'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:x:g:nfo:b:s:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gam_idx_name = optarg;
            break;

        case 'l':
            loci_file = optarg;
            break;

        case 'x':
            xg_idx_name = optarg;
            break;

        case 'n':
            name_alleles = true;
            break;

        case 'f':
            forwardize = true;
            break;

        case 'o':
            loci_out = optarg;
            break;

        case 's':
            sorted_loci = optarg;
            break;

        case 'b':
            n_best = atoi(optarg);
            name_alleles = true;
            break;

        case 'h':
        case '?':
            help_locify(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (!gam_idx_name.empty()) {
        gam_idx.open_read_only(gam_idx_name);
    }

    if (xg_idx_name.empty()) {
        cerr << "[vg locify] Error: no xg index provided" << endl;
        return 1;
    }
    ifstream xgstream(xg_idx_name);
    xg::XG xgidx(xgstream);

    std::function<vector<string>(string, char)> strsplit = [&](string x, char delim){

        vector<string> ret;
        stringstream ss;
        std::string tok;
        while (getline(ss, tok, delim)){
            ret.push_back(tok);
        }
        return ret;

    };

    vector<string> locus_names;
    map<string, map<string, int > > locus_allele_names;
    map<string, Alignment> alignments_with_loci;
    map<pos_t, set<string> > pos_to_loci;
    map<string, set<pos_t> > locus_to_pos;
    map<string, map<int, int> > locus_allele_support;
    map<string, vector<int> > locus_to_best_n_alleles;
    map<string, set<int> > locus_to_keep;
    int count = 0;

    std::function<void(Locus&)> lambda = [&](Locus& l){
        locus_names.push_back(l.name());
        set<vg::id_t> nodes_in_locus;
        for (int i = 0; i < l.allele_size(); ++i) {
            auto& allele = l.allele(i);
            for (int j = 0; j < allele.mapping_size(); ++j) {
                auto& position = allele.mapping(j).position();
                nodes_in_locus.insert(position.node_id());
            }
            // for position in mapping
            map<pos_t, int> ref_positions;
            map<pos_t, Edit> edits;
            decompose(allele, ref_positions, edits);
            // warning: uses only reference positions!!!
            for (auto& pos : ref_positions) {
                pos_to_loci[pos.first].insert(l.name());
                locus_to_pos[l.name()].insert(pos.first);
            }
        }
        // void for_alignment_in_range(int64_t id1, int64_t id2, std::function<void(const Alignment&)> lambda);
        std::function<void(const Alignment&)> fill_alns = [&](const Alignment& a){
            // TODO reverse complementing alleles ?
            // overlap is stranded
            //matching
            // find the most-matching allele
            map<double, vector<int> > matches;
            for (int i = 0; i < l.allele_size(); ++i) {
                auto& allele = l.allele(i);
                matches[overlap(a.path(), allele)].push_back(i);
            }
            assert(l.allele_size());
            int best = matches.rbegin()->second.front();
            Locus matching;
            matching.set_name(l.name());
            if (name_alleles) {
                //map<string, map<string, int > > locus_allele_names;
                auto& allele = l.allele(best);
                string s;
                allele.SerializeToString(&s);
                auto& l_names = locus_allele_names[l.name()];
                auto f = l_names.find(s);
                int name_int = 0;
                if (f == l_names.end()) {
                    int next_id = l_names.size() + 1;
                    l_names[s] = next_id;
                    name_int = next_id;
                } else {
                    name_int = f->second;
                }
                string allele_name = vg::convert(name_int);
                Path p;
                p.set_name(allele_name);
                *matching.add_allele() = p;
                if (n_best) {
                    // record support for this allele
                    // we'll use to filter the locus records later
                    locus_allele_support[l.name()][name_int]++;
                }
            } else {
                *matching.add_allele() = l.allele(best);
                // TODO get quality score relative to this specific allele / alignment
                // record in the alignment we'll save
            }
            if (alignments_with_loci.find(a.name()) == alignments_with_loci.end()) {
                alignments_with_loci[a.name()] = a;
            }
            Alignment& aln = alignments_with_loci[a.name()];
            *aln.add_locus() = matching;
        };
        vector<vg::id_t> nodes_vec;
        for (auto& id : nodes_in_locus) nodes_vec.push_back(id);
        gam_idx.for_alignment_to_nodes(nodes_vec, fill_alns);
    };

    if (!loci_file.empty()){
        ifstream ifi(loci_file);
        stream::for_each(ifi, lambda);
    } else {
        cerr << "[vg locify] Warning: empty locus file given, could not annotate alignments with loci." << endl;
    }

    // find the non-nested loci
    vector<string> non_nested_loci;
    for (auto& name : locus_names) {
        // is it nested?
        auto& positions = locus_to_pos[name];
        int min_loci = 0;
        for (auto& pos : positions) {
            auto& loci = pos_to_loci[pos];
            min_loci = (min_loci == 0 ? (int)loci.size() : min(min_loci, (int)loci.size()));
        }
        if (min_loci == 1) {
            // not fully contained in any other locus
            non_nested_loci.push_back(name);
        }
    }

    // filter out the non-best alleles
    if (n_best) {
        // find the n-best
        for (auto& supp : locus_allele_support) {
            auto& name = supp.first;
            auto& alleles = supp.second;
            map<int, int> ranked;
            for (auto& allele : alleles) {
                ranked[allele.second] = allele.first;
            }
            auto& to_keep = locus_to_keep[name];
            for (auto r = ranked.rbegin(); r != ranked.rend(); ++r) {
                to_keep.insert(r->second);
                if (to_keep.size() == n_best) {
                    break;
                }
            }
        }
        // filter out non-n-best from the alignments
        for (auto& a : alignments_with_loci) {
            auto& aln = a.second;
            vector<Locus> kept;
            for (int i = 0; i < aln.locus_size(); ++i) {
                auto& allele = aln.locus(i).allele(0);
                if (locus_to_keep[aln.locus(i).name()].count(atoi(allele.name().c_str()))) {
                    kept.push_back(aln.locus(i));
                }
            }
            aln.clear_locus();
            for (auto& l : kept) {
                *aln.add_locus() = l;
            }
        }
    }

    if (n_best && !loci_out.empty()) {
        // filter out non-n-best from the loci
        if (!loci_file.empty()){
            ofstream outloci(loci_out);
            vector<Locus> buffer;
            std::function<void(Locus&)> lambda = [&](Locus& l){
                // remove the alleles which are to filter
                //map<string, map<string, int > > locus_allele_names;
                auto& allele_names = locus_allele_names[l.name()];
                auto& to_keep = locus_to_keep[l.name()];
                vector<Path> alleles_to_keep;
                for (int i = 0; i < l.allele_size(); ++i) {
                    auto allele = l.allele(i);
                    string s; allele.SerializeToString(&s);
                    auto& name = allele_names[s];
                    if (to_keep.count(name)) {
                        allele.set_name(vg::convert(name));
                        alleles_to_keep.push_back(allele);
                    }
                }
                l.clear_allele();
                for (auto& allele : alleles_to_keep) {
                    *l.add_allele() = allele;
                }
                buffer.push_back(l);
                stream::write_buffered(outloci, buffer, 100);
            };
            ifstream ifi(loci_file);
            stream::for_each(ifi, lambda);
            stream::write_buffered(outloci, buffer, 0);
            outloci.close();
        } else {
            cerr << "[vg locify] Warning: empty locus file given, could not update loci." << endl;
        }
    }

    // sort them using... ? ids?
    sort(non_nested_loci.begin(), non_nested_loci.end(),
         [&locus_to_pos](const string& s1, const string& s2) {
             return *locus_to_pos[s1].begin() < *locus_to_pos[s2].begin();
         });

    if (!sorted_loci.empty()) {
        ofstream outsorted(sorted_loci);
        for (auto& name : non_nested_loci) {
            outsorted << name << endl;
        }
        outsorted.close();
    }

    vector<Alignment> output_buf;
    for (auto& aln : alignments_with_loci) {
        // TODO order the loci by their order in the alignments
        if (forwardize) {
            if (aln.second.path().mapping_size() && aln.second.path().mapping(0).position().is_reverse()) {
                output_buf.push_back(reverse_complement_alignment(aln.second,
                                                                  [&xgidx](int64_t id) { return xgidx.node_length(id); }));
            } else {
                output_buf.push_back(aln.second);
            }
        } else {
            output_buf.push_back(aln.second);
        }
        stream::write_buffered(cout, output_buf, 100);
    }
    stream::write_buffered(cout, output_buf, 0);        
    
    return 0;
}

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] -p <PATH> <my_graph>.vg" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "--path / -p     REQUIRED: A reference path to deconstruct against." << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    //cerr << "WARNING: EXPERIMENTAL" << endl;
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    string graphname;
    string outfile = "";
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {0, 0, 0, 0}

            };

            int option_index = 0;
            c = getopt_long (argc, argv, "hp:",
                    long_options, &option_index);

            // Detect the end of the options.
            if (c == -1)
                break;

            switch (c)
            {
                case 'p':
                    refpaths = split(optarg, ",");
                    break;
                case '?':
                case 'h':
                    help_deconstruct(argv);
                    return 1;
                default:
                    help_deconstruct(argv);
                    abort();
            }

        }
        graphname = argv[optind];
        vg::VG* graph;
        if (!graphname.empty()){
            ifstream gstream(graphname);
            graph = new vg::VG(gstream);
        }

        // load graph

        // Deconstruct
        Deconstructor dd;
        dd.deconstruct(refpaths, graph);
    return 0;
}


void help_sort(char** argv){
    cerr << "usage: " << argv[0] << " sort [options] -i <input_file> -r <reference_name> > sorted.vg " << endl
         << "options: " << endl
         << "           -g, --gfa              input in GFA format" << endl
         << "           -i, --in               input file" << endl
         << "           -r, --ref              reference name" << endl
         << "           -w, --without-grooming no grooming mode" << endl
         << "           -f, --fast             sort using Eades algorithm, otherwise max-flow sorting is used" << endl   
         << endl;
}

int main_sort(int argc, char *argv[]) {

    //default input format is vg
    bool gfa_input = false;
    string file_name = "";
    string reference_name = "";
    bool without_grooming = false;
    bool use_fast_algorithm = false;
    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"gfa", no_argument, 0, 'g'},
                {"in", required_argument, 0, 'i'},
                {"ref", required_argument, 0, 'r'},
                {"without-grooming", no_argument, 0, 'w'},
                {"fast", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "i:r:gwf",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gfa_input = true;
            break;
        case 'r':
            reference_name = optarg;
            break;
        case 'i':
            file_name = optarg;
            break;
        case 'w':
            without_grooming = true;
            break;
        case 'f':
            use_fast_algorithm = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_sort(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }
  
    if (reference_name.empty() || file_name.empty()) {
        help_sort(argv);
        exit(1);
    }
    
    ifstream in;
    std::unique_ptr<VG> graph;
    {
        in.open(file_name.c_str());        
        if (gfa_input) {
            graph.reset(new VG());
            graph->from_gfa(in);
        } else {
            graph.reset(new VG(in));
        }
    }
    FlowSort flow_sort(*graph.get());
    if (use_fast_algorithm) {
        flow_sort.fast_linear_sort(reference_name, !without_grooming);
    } else {
        flow_sort.max_flow_sort(reference_name);
    }
    
    graph->serialize_to_ostream(std::cout);
    in.close();
    return 0;
}

// No help_test is necessary because the unit testing library takes care of
// complaining about missing options.

int main_test(int argc, char** argv){
    // Forward arguments along to the main unit test driver
    return vg::unittest::run_unit_tests(argc, argv);
}

void vg_help(char** argv) {
    cerr << "vg: variation graph tool, version " << VG_VERSION_STRING << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl;
         
     vg::subcommand::Subcommand::for_each([](const vg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
         
     // Also announce all the old-style hardcoded commands
     cerr << "  -- deconstruct   convert a graph into VCF relative to a reference." << endl
         << "  -- genotype      compute genotypes from aligned reads" << endl
         << "  -- sort          sort variant graph using max flow algorithm or Eades fast heuristic algorithm" << endl
         << "  -- test          run unit tests" << endl;
}

int main(int argc, char *argv[])
{

    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        vg_help(argv);
        return 1;
    }
    
    auto* subcommand = vg::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    }
    
    // Otherwise, fall abck on the old chain of if statements.

    //omp_set_dynamic(1); // use dynamic scheduling

    string command = argv[1];
    if (command == "deconstruct"){
        return main_deconstruct(argc, argv);
    } else if (command == "test") {
        return main_test(argc, argv);
    } else if (command == "locify"){
        return main_locify(argc, argv);
    } else if (command == "sort") {
        return main_sort(argc, argv);
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
