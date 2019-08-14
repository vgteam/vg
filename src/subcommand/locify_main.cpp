/** \file locify_main.cpp
 *
 * Defines the "vg locify" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../index.hpp"
#include "../convert.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

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
            n_best = parse<int>(optarg);
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
    unique_ptr<PathPositionHandleGraph> xgidx = vg::io::VPKG::load_one<PathPositionHandleGraph>(xgstream);

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
        vg::io::for_each(ifi, lambda);
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
                vg::io::write_buffered(outloci, buffer, 100);
            };
            ifstream ifi(loci_file);
            vg::io::for_each(ifi, lambda);
            vg::io::write_buffered(outloci, buffer, 0);
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
                                                                  [&xgidx](int64_t id) { return xgidx->get_length(xgidx->get_handle(id)); }));
            } else {
                output_buf.push_back(aln.second);
            }
        } else {
            output_buf.push_back(aln.second);
        }
        vg::io::write_buffered(cout, output_buf, 100);
    }
    vg::io::write_buffered(cout, output_buf, 0);        
    
    return 0;
}

// Register subcommand
static Subcommand vg_locify("locify", "find loci", main_locify);

