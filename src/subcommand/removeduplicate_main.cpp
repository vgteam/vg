/*
 * Define the "vg remove_duplicate" subcommand, which remove the duplicate PCRs
 *
 */
#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>
#include <cstdint>

#include <subcommand.hpp>

#include "../algorithms/alignment_path_offsets.hpp"
#include "../multipath_alignment.hpp"
#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../stream_sorter.hpp"
#include <BooPHF.h>
#include <bitset>
#include "../hash_map.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <iostream>
#include <vector>
#include <algorithm>


class Custom_string_hasher {
public:


    uint64_t operator()(const string &s, uint64_t seed = 0) const {
        size_t hash = seed;
        hash_combine(hash, s);
        return hash;


    }
};

vector <pair<long long, long long>> make_coalesced_sorted_intervals(const Alignment &aln) {
    /// return a vector of pairs of id-ts of the node ids that are themselves sorted and properly coalesced
    vector<long long> sorted_node_ids;
    vector <pair<long long, long long>> intervals;
    for (const Mapping mapping: aln.path().mapping()) {
        sorted_node_ids.push_back(mapping.position().node_id());
    }

    if (sorted_node_ids.empty()) {
        return intervals;
    }
    // We have a sorted vector of node_ids
    sort(sorted_node_ids.begin(), sorted_node_ids.end());


    long long start = sorted_node_ids[0];
    long long end = sorted_node_ids[0];
    for (long long i = 1; i < sorted_node_ids.size(); i++) {
        if (sorted_node_ids[i] == end + 1) {
            end = sorted_node_ids[i];
        } else {
            intervals.push_back(make_pair(start, end));
            start = end = sorted_node_ids[i];
        }
    }

    intervals.push_back(make_pair(start, end));
    return intervals;
}


vector <pos_t> alignment_position(const Alignment &aln) {
    /// return a vector of pos_t variables which the index of the vector points to the read offset and the value points to the node pos_t
    vector <pos_t> positions;
    for (const Mapping &mapping: aln.path().mapping()) {
        long long node_offset = mapping.position().offset(); // This is the position of the first Edit
        long long node_id = mapping.position().node_id(); // The id of the node we are working with
        bool is_reverse = mapping.position().is_reverse();

        for (const Edit &edit: mapping.edit()) {

            // Match or Mismatch
            if (edit.from_length() == edit.to_length()) {

                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset + i));
                }
                node_offset += edit.from_length();

            } else if (edit.to_length() == 0 && edit.from_length() > edit.to_length()) { // Deletion
                node_offset += edit.from_length();

            } else if (edit.from_length() < edit.to_length()) { // insertion
                node_offset += edit.from_length();
                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset));
                }
            }


        }
    }
    return positions;

}

bool check_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have one equal position on one base
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);
    for (size_t i = 1; i < first_alignment_pos.size(); ++i) {
        if (id(first_alignment_pos[i]) == id(second_alignment_pos[i]) &&
            offset(first_alignment_pos[i]) == offset(second_alignment_pos[i]) &&
            (is_rev(first_alignment_pos[i]) == is_rev(second_alignment_pos[i]))) {
            size_t t = first_alignment_pos.size() - 1;
            // We check if the alignments end on the same position
            if (id(first_alignment_pos[t]) == id(second_alignment_pos[t]) &&
                offset(first_alignment_pos[t]) == offset(second_alignment_pos[t]) &&
                (is_rev(first_alignment_pos[t]) == is_rev(second_alignment_pos[t])))
                return true;
        }

    }
    return false;


}

bool check_pair_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have same begining and end
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);
    if (id(first_alignment_pos.front()) == id(second_alignment_pos.front()) &&
        offset(first_alignment_pos.front()) == offset(second_alignment_pos.front()) &&
        (is_rev(first_alignment_pos.front()) == is_rev(second_alignment_pos.front()))) {
//        size_t t = first_alignment_pos.size() - 1;
        // We check if the alignments end on the same position
        if (id(first_alignment_pos.back()) == id(second_alignment_pos.back()) &&
            offset(first_alignment_pos.back()) == offset(second_alignment_pos.back()) &&
            (is_rev(first_alignment_pos.back()) == is_rev(second_alignment_pos.back())))
            return true;
    }

    return false;


}

string name_id(const Alignment &aln) {
    /// Make a unique name for each alignment base on the name and if they have next/prev fragment
    string id;
    if (aln.has_fragment_next()) {
        id = aln.name() + "/1";
    } else if (aln.has_fragment_prev()) {
        id = aln.name() + "/2";
    } else {
        id = aln.name();
    }
    return id;

}


void help_rmvdup(char **argv) {
// TODO: add whatever option is needed to this list. Change long_option and getopt_long if want to add an option.
// TODO: see what input and output file formats is possible
    cerr << "usage: " << argv[0] << " rmvdup [options] inputfile.gam > output.gam " << endl
         << "Remove duplicate PCRs from the input file. A gam index file (.gam.gai) must exists." << endl
         << "  -o, --output_type               prints the pairs of (sequence, name) as strings in the output" << endl
         << "    -t, --threads N            number of threads to use" << endl
         << "   -g --get_duplicates            prints the duplicate candidates in the output(can use with the -o)"
         << endl;

}


typedef boomphf::mphf <string, Custom_string_hasher> boophf_t;

int main_rmvdup(int argc, char *argv[]) {
    string filename;
    bool output_t = false;
    bool print_duplicates = false;
    int threads = 1;

    int c;
    optind = 2;  // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
                {"help",           no_argument,       0, 'h'},
                {"output_type",    no_argument,       0, 'o'},
                {"get_duplicates", no_argument,       0, 'g'},
                {"threads",        required_argument, 0, 't'},
                {0,                0,                 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hogt:",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            case 'o':
                output_t = true;
                break;

            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
                break;

            case 'g':
                print_duplicates = true;
                break;

            case 'h':
            case '?':
                help_rmvdup(argv);
                exit(1);
                break;

            default:
                help_rmvdup(argv);
                abort();
        }
    }

    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }


    string sorted_gam_name = get_input_file_name(optind, argc, argv);
    unique_ptr <GAMIndex> gam_index;
    if (!sorted_gam_name.empty()) {
        // Load the GAM index
        gam_index = unique_ptr<GAMIndex>(new GAMIndex());
        get_input_file(sorted_gam_name + ".gai", [&](istream &in) {
            // We get it form the appropriate .gai, which must exist
            gam_index->load(in);
        });
    }

    vector <string> keys;
    function<void(Alignment & )> fill_hash = [&](const Alignment &aln) {
        keys.push_back(name_id(aln));
    };

    get_input_file(sorted_gam_name, [&](istream &in) {
        vg::io::for_each(in, fill_hash);
    });

    vector<bool> checked(keys.size(), false);

    boophf_t *bphf = new boomphf::mphf<string, Custom_string_hasher>(keys.size(), keys, threads, 2.0, false, false);

    std::unique_ptr<vg::io::ProtobufEmitter<Alignment>> emitter;
    emitter = std::unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
    function<void(Alignment & )> pcr_removal = [&](Alignment &aln) {
        if (gam_index.get() != nullptr) {
            // I mark all the reads that have to get deleted as duplicates. This means one read from each duplicate set remains unmarked.
            // This way we can remove all reads with duplicate flag and not worry about deleting them all
            if (!checked[bphf->lookup(name_id(aln))]) {
                // make the alignment nodes list that can be use as input if .find function of the gam_index
                vector <pair<long long, long long>> intervals = make_coalesced_sorted_intervals(aln);
                // Find all alignments that share at least one node with the current working alignment
                get_input_file(sorted_gam_name, [&](istream &input_gam) {
                    vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);
                    // find all sharing nodes alignments and call the function to handle the result
                    gam_index->find(gam_cursor, intervals, [&](const Alignment &share_aln) {

                        if (!checked[bphf->lookup(name_id(share_aln))]) {
                            if (name_id(aln) != name_id(share_aln)) {
                                if (check_duplicate(aln, share_aln)) {
                                    if (print_duplicates) {
#pragma omp critical (cerr)
                                        if (output_t)
                                            cout << share_aln.sequence() << endl;
                                        else
                                            emitter->write(std::move(const_cast<Alignment &>(share_aln)));
                                    }
                                    checked[bphf->lookup(name_id(share_aln))] = true;
                                }

                            }

                        }
                    });

                });
                // Writing the non-duplicate PCRs to the file
                if (!checked[bphf->lookup(name_id(aln))]) {
#pragma omp critical (cerr)
                    checked[bphf->lookup(name_id(aln))] = true;
                    if (!print_duplicates) {
                        if (output_t)
                            cout << aln.sequence() << "\t" << aln.name() << endl;
                        else
                            emitter->write(std::move(aln));
                    }
                }


            }

        }
    };

//    vector<unordered_map <string, Alignment>> memory(threads);
    vector<unordered_set<std::string>> pair1_duplicates_set(threads);


    function<void(Alignment &, Alignment & )> pcr_removal_pair_end_double = [&](Alignment &aln, Alignment &aln_pair) {

        int thread_number = omp_get_thread_num();

        if (gam_index.get() != nullptr) {
            if (!checked[bphf->lookup(name_id(aln))]) {
                checked[bphf->lookup(name_id(aln))] = true;
                checked[bphf->lookup(name_id(aln_pair))] = true;

                vector <pair<long long, long long>> intervals_pair1 = make_coalesced_sorted_intervals(aln);


                // Set of duplicate alns with the first pair
//                unordered_set <string> pair1_duplicates_set;

                get_input_file(sorted_gam_name, [&](istream &input_gam) {

                    vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);

                    // This finds all other alns that share node/nodes with the "aln"
                    gam_index->find(gam_cursor, intervals_pair1, [&](const Alignment &share_aln) {

                        if (check_duplicate(aln, share_aln)) {
                            pair1_duplicates_set[thread_number].insert(
                                    share_aln.has_fragment_prev() ? share_aln.fragment_prev().name()
                                                                  : share_aln.fragment_next().name());
                        }


                    });


                });


                vector <pair<long long, long long>> intervals_pair2 = make_coalesced_sorted_intervals(aln_pair);
                get_input_file(sorted_gam_name, [&](istream &input_gam) {
                    vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);

                    gam_index->find(gam_cursor, intervals_pair2, [&](const Alignment &share_aln) {
                        if (check_duplicate(aln_pair, share_aln)) {

                            // The pair is in the set
                            if (pair1_duplicates_set[thread_number].find(share_aln.name()) != pair1_duplicates_set[thread_number].end()) {


                                // This means we find pairs that are duplicates with the main pairs so we check them for being duplicate

                                // we add /1 to the prev one because we are naming the other pair
                                string pair_name_id = share_aln.has_fragment_prev() ?
                                                      share_aln.fragment_prev().name() + "/1" :
                                                      share_aln.fragment_next().name() + "/2";
//                                        cout << "share aln " << name_id(share_aln) << endl;
//                                        cout << "pair share aln " << pair_name_id << endl;
                                checked[bphf->lookup(name_id(share_aln))] = true;
                                checked[bphf->lookup(pair_name_id)] = true;

                                if (print_duplicates) {
#pragma omp critical (cerr)
                                    if (output_t){
                                        cout << "Pair 1" << share_aln.sequence() << endl;
                                        cout << "Pair 2" << aln.sequence() << endl;
                                    }


                                    else{
                                        emitter->write(std::move(const_cast<Alignment &>(share_aln)));
                                        emitter->write(std::move(const_cast<Alignment &>(aln)));
                                    }

                                }


                            }


                        }

                    });

                });

#pragma omp critical (cerr)
//                        cout << name_id(aln) << endl;
//                checked[bphf->lookup(name_id(aln))] = true;
//                checked[bphf->lookup(name_id(aln_pair))] = true;
                if (!print_duplicates) {
                    if (output_t)
                        cout << "Pair1 " << aln.sequence() << "\t" << aln.name() << "Pair2 "
                             << aln_pair.sequence() << "\t" << aln_pair.name() << endl;
                    else {
                        emitter->write(std::move(aln));
                        emitter->write(std::move(aln_pair));
                    }

                }




            }

        }


    };


//
//    // This is a memory for finding pairs. The name of the aln as the key and the aln as the value
//    vector<unordered_map <string, Alignment>> memory(threads);
////    unordered_map <string, Alignment> temp_memory;
//    // This is the function that works on pair_end data
//    function<void(Alignment & )> pcr_removal_pair_end = [&](Alignment &aln) {
//        int thread_number = omp_get_thread_num();
//        // TODO: Check this. I have to clear all the used alignments here and then use the memory again. right?
////        unordered_map<string, Alignment>().swap(memory[thread_number]);
//
//
//        if (gam_index.get() != nullptr) {
//            // If this alignment is not already checked for being duplicate
//            if (!checked[bphf->lookup(name_id(aln))]) {
//                // This is when the read has a pair
//                if (aln.has_fragment_prev() || aln.has_fragment_next()) {
//                    string pair_name = aln.has_fragment_prev() ? aln.fragment_prev().name()
//                                                               : aln.fragment_next().name();
//                    if (memory[thread_number].find(aln.name()) == memory[thread_number].end()) {
//                        // This means the we don't have the pair yet
//                        memory[thread_number][pair_name] = aln;
//                    } else {
//                        // We have the pair of our alignment
//                        Alignment aln_pair = memory[thread_number][aln.name()];
//                        memory[thread_number].erase(aln.name());
//
//                        // We have to find all the alignments that share nodes with both pairs and find if they are both end of a pair
//                        // TODO: For now I check the equality of pairs using check_duplicate function which check
//                        //  if they end the same and if they have at least one more same base. This could gets better.
//
//                        vector <pair<long long, long long>> intervals_pair1 = make_coalesced_sorted_intervals(aln);
//
//
//                        // Set of duplicate alns with the first pair
//                        unordered_set <string> pair1_duplicates_set;
//                        get_input_file(sorted_gam_name, [&](istream &input_gam) {
//
//                            vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);
//
//                            // This finds all other alns that share node/nodes with the "aln"
//                            gam_index->find(gam_cursor, intervals_pair1, [&](const Alignment &share_aln) {
//                                if (check_duplicate(aln, share_aln)) {
//                                    pair1_duplicates_set.insert(
//                                            share_aln.has_fragment_prev() ? share_aln.fragment_prev().name()
//                                                                          : share_aln.fragment_next().name());
//                                }
//
//                            });
//
//
//                        });
//
//                        vector <pair<long long, long long>> intervals_pair2 = make_coalesced_sorted_intervals(aln_pair);
//                        get_input_file(sorted_gam_name, [&](istream &input_gam) {
//                            vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);
//
//                            gam_index->find(gam_cursor, intervals_pair2, [&](const Alignment &share_aln) {
//                                if (check_duplicate(aln_pair, share_aln)) {
//
//                                    // The pair is in the set
//                                    if (pair1_duplicates_set.find(share_aln.name()) != pair1_duplicates_set.end()) {
//
//
//                                        // This means we find pairs that are duplicates with the main pairs so we check them for being duplicate
//
//                                        // we add /1 to the prev one because we are naming the other pair
//                                        string pair_name_id = share_aln.has_fragment_prev() ?
//                                                              share_aln.fragment_prev().name() + "/1" :
//                                                              share_aln.fragment_next().name() + "/2";
////                                        cout << "share aln " << name_id(share_aln) << endl;
////                                        cout << "pair share aln " << pair_name_id << endl;
//                                        // TODO: WE dont get here!
//                                        checked[bphf->lookup(name_id(share_aln))] = true;
//                                        checked[bphf->lookup(pair_name_id)] = true;
//
//                                        // TODO: for now I am just printing the second pair because I am not saving all
//                                        //  the alignments for the first pair. I am just saving the name for efficient mem usage.
//                                        if (print_duplicates) {
//#pragma omp critical (cerr)
//                                            if (output_t)
//                                                cout << "End Pair" << share_aln.sequence() << endl;
//                                            else
//                                                emitter->write(std::move(const_cast<Alignment &>(share_aln)));
//                                        }
//
//
//                                    }
//
//
//                                }
//
//                            });
//
//                        });
//
//                        // Writing the non-duplicate PCRs to the file
//#pragma omp critical (cerr)
////                        cout << name_id(aln) << endl;
//                        checked[bphf->lookup(name_id(aln))] = true;
//                        checked[bphf->lookup(name_id(aln_pair))] = true;
//                        if (!print_duplicates) {
//                            if (output_t)
//                                cout << "Pair1 " << aln.sequence() << "\t" << aln.name() << "Pair2 "
//                                     << aln_pair.sequence() << "\t" << aln_pair.name() << endl;
//                            else {
//                                emitter->write(std::move(aln));
//                                emitter->write(std::move(aln_pair));
//                            }
//
//                        }
//
//
//                    }
//
//                } else { // When the read is not a pair_end read
//                    pcr_removal(aln);
//
//                }
//
//
//            }
//
//        }
//
//
//    };

    if (gam_index.get() != nullptr) {
        get_input_file(sorted_gam_name, [&](istream &in) {
//            vg::io::for_each(in, pcr_removal);
//            vg::io::for_each_parallel_shuffled(in, pcr_removal_pair_end);
            vg::io::for_each_parallel_shuffled_double(in, pcr_removal_pair_end_double);



        });
    }


    return 0;

}


// Register subcommand
static Subcommand vg_removeduplicate("rmvdup", "Remove duplicate PCRs from the input file", WIDGET, main_rmvdup);