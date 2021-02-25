// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#include "deconstructor.hpp"
#include "traversal_finder.hpp"

//#define debug

using namespace std;


namespace vg {
Deconstructor::Deconstructor(){

}
Deconstructor::~Deconstructor(){
}

/**
 * Takes in a vector of snarltraversals
 * returns their sequences as a vector<string>
 * returns a boolean hasRef
 * if a reference path is present, hasRef is set to true and the first
 * string in the vector is the reference allele
 * otherwise, hasRef is set to false and all strings are alt alleles.
 */
vector<int> Deconstructor::get_alleles(vcflib::Variant& v, const vector<SnarlTraversal>& travs, int ref_path_idx,
                                       char prev_char, bool use_start) {

    assert(ref_path_idx >=0 && ref_path_idx < travs.size());

    // map strings to allele numbers
    // (we are using the traversal finder in such a way that duplicate alleles can get returned
    // in order to be able to preserve the path names)
    map<string, int> allele_idx;
    size_t cur_alt = 1;

    // go from traversals number (offset in travs) to allele number
    vector<int> trav_to_allele(travs.size());

    // compute the allele as a string
    auto trav_to_string = [&](const SnarlTraversal& trav) {
        string allele;
        // we skip the snarl endpoints
        for (int j = 1; j < trav.visit_size() - 1; ++j) {
            const string& node_sequence = graph->get_sequence(graph->get_handle(trav.visit(j).node_id()));
            allele += trav.visit(j).backward() ? reverse_complement(node_sequence) : node_sequence;
        }
        return toUppercase(allele);
    };
       

    // set the reference allele
    string ref_allele = trav_to_string(travs.at(ref_path_idx));
    allele_idx[ref_allele] = 0;
    trav_to_allele[ref_path_idx] = 0;
    bool substitution = true;
        
    // set the other alleles (they can end up as 0 alleles too if their strings match the reference)
    for (int i = 0; i < travs.size(); ++i) {
        if (i != ref_path_idx) {
            string allele = trav_to_string(travs[i]);
            auto ai_it = allele_idx.find(allele);
            if (ai_it == allele_idx.end()) {
                // make a new allele for this string
                allele_idx[allele] = cur_alt;
                trav_to_allele.at(i) = cur_alt;
                ++cur_alt;
                substitution = substitution && allele.size() == ref_allele.size();
            } else {
                // allele string has been seen, map this traversal to it
                trav_to_allele.at(i) = ai_it->second;
            }
        }
    }

    // fill in the variant
    v.alleles.resize(allele_idx.size());
    assert(allele_idx.size() > 0);
    v.alt.resize(allele_idx.size() - 1);

    for (auto ai_pair : allele_idx) {
        string allele_string = ai_pair.first;
        if (!use_start) {
            reverse_complement_in_place(allele_string);
        }
        if (!substitution) {
            allele_string = string(1, prev_char) + allele_string;
        }
        v.alleles[ai_pair.second] = allele_string;
        if (ai_pair.second > 0) {
            v.alt[ai_pair.second - 1] = allele_string;
        } else {
            v.ref = allele_string;
        }
    }

    // shift our variant back if it's an indel
    if (!substitution) {
        assert(v.position >= 2);
        --v.position;
    }

    v.updateAlleleIndexes();

    return trav_to_allele;
}

void Deconstructor::get_genotypes(vcflib::Variant& v, const vector<string>& names,
                                  const vector<int>& trav_to_allele) {
    assert(names.size() == trav_to_allele.size());
    // set up our variant fields
    v.format.push_back("GT");
    if (path_to_sample) {
        v.format.push_back("PI");
    }

    // get a list of traversals for every vcf sample
    // (this will be 1:1 unless we're using the path_to_sample name map)
    map<string, vector<int> > sample_to_traversals;
    for (int i = 0; i < names.size(); ++i) {
        string sample_name;
        if (path_to_sample && path_to_sample->count(names[i])) {
            sample_name = path_to_sample->find(names[i])->second;
        } else {
            sample_name = names[i];
        }
        if (sample_names.count(sample_name)) {
            if (sample_to_traversals.count(sample_name)) {
                sample_to_traversals[sample_name].push_back(i);
            } else {
                sample_to_traversals[sample_name] = {i};
            }
        }
    }

    // write out the genotype for each sample
    // if we're mapping a vg path name to its prefix for the sample name, we stick some information about the full
    // path name in the PI part of format
    set<string> conflicts;
    for (auto& sample_name : sample_names) {
        if (sample_to_traversals.count(sample_name)) {
            const vector<int>& travs = sample_to_traversals[sample_name];
            assert(!travs.empty());
            vector<int> chosen_travs;
            bool conflict;
            std::tie(chosen_travs, conflict) = choose_traversals(travs, trav_to_allele, names);
            if (conflict) {
                conflicts.insert(sample_name);
            }
            string genotype = std::to_string(trav_to_allele[chosen_travs[0]]);
            for (int i = 1; i < chosen_travs.size(); ++i) {
                genotype += "/" + std::to_string(trav_to_allele[chosen_travs[i]]);
            }
            v.samples[sample_name]["GT"] = {genotype};
            if (path_to_sample) {
                for (auto trav : travs) {
                    v.samples[sample_name]["PI"].push_back(names[trav] + "=" + std::to_string(trav_to_allele[trav]));
                }
            }
        } else {
            v.samples[sample_name]["GT"] = {"."};
            if (path_to_sample) {
                v.samples[sample_name]["PI"] = {"."};
            }
        }
    }
    for (auto& conflict_sample : conflicts) {
        v.info["CONFLICT"].push_back(conflict_sample);
    }
}

pair<vector<int>, bool> Deconstructor::choose_traversals(const vector<int>& travs, const vector<int>& trav_to_allele,
                                                         const vector<string>& trav_to_name) {
    assert(!travs.empty());
    // count the number of times each allele comes up in a traversal
    vector<int> allele_frequencies(trav_to_allele.size(), 0);
    for (auto trav : travs) {
        ++allele_frequencies[trav_to_allele[trav]];
    }
    // sort on frquency
    function<bool(int, int)> comp = [&] (int trav1, int trav2) {
        if (allele_frequencies[trav_to_allele[trav1]] < allele_frequencies[trav_to_allele[trav2]]) {
            return true;
        } else if (allele_frequencies[trav_to_allele[trav1]] == allele_frequencies[trav_to_allele[trav2]]) {
            // prefer non-ref when possible
            if (trav_to_allele[trav1] == 0 && trav_to_allele[trav2] != 0) {
                return true;
            }
            // or break tie using lex order on path name
            else {
                return trav_to_name[trav1] < trav_to_name[trav2];
            }
        } else {
            return false;
        }
    };
    vector<int> sorted_travs = travs;
    std::sort(sorted_travs.begin(), sorted_travs.end());
    
    // find the <ploidy> most frequent traversals
    vector<int> most_frequent_travs;
    for (int i = sorted_travs.size() - 1; i >= 0 && most_frequent_travs.size() < ploidy; --i) {
        most_frequent_travs.push_back(sorted_travs[i]);
    }

    // check if there's a conflict
    size_t zero_count = std::count(allele_frequencies.begin(), allele_frequencies.end(), 0);
    bool conflict = allele_frequencies.size() - zero_count > ploidy;

    return make_pair(most_frequent_travs, conflict);
}
    
bool Deconstructor::deconstruct_site(const Snarl* snarl) {

    auto contents = snarl_manager->shallow_contents(snarl, *graph, false);
    if (contents.first.empty()) {
        // Nothing but the boundary nodes in this snarl
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping empty site " << pb2json(*snarl) << endl;
#endif
        return false;
    }
#ifdef debug
#pragma omp crtiical (cerr)
    cerr << "Computing traversals of site " << pb2json(*snarl) << endl;
#endif

    // find every traversal that runs through a path in the graph
    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs;
    path_travs = path_trav_finder->find_path_traversals(*snarl);
    vector<string> path_trav_names;

    // pick out the traversal corresponding to a reference path, breaking ties consistently
    string ref_trav_name;
    for (int i = 0; i < path_travs.first.size(); ++i) {
        string path_trav_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
#ifdef debug
#pragma omp critical (cerr)
        {
            cerr << "Traversal " << i << ": name=" << path_trav_name << ", size=" << path_travs.first[i].visit_size()
                 << ", start=" << graph->get_position_of_step(path_travs.second[i].first)
                 << ", end=" << graph->get_position_of_step(path_travs.second[i].second) << endl
                 << " trav=" << pb2json(path_travs.first[i]) << endl;
        }
#endif
        if (ref_paths.count(path_trav_name) &&
            (ref_trav_name.empty() || path_trav_name < ref_trav_name)) {
            ref_trav_name = path_trav_name;
        }
        path_trav_names.push_back(path_trav_name);
    }

    // remember all the reference traversals (there can be more than one only in the case of a
    // cycle in the reference path
    vector<int> ref_travs;
    if (!ref_trav_name.empty()) {
        for (int i = 0; i < path_travs.first.size(); ++i) {
            string path_trav_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
            if (path_trav_name == ref_trav_name) {
                ref_travs.push_back(i);
            }
        }
    }

    // there's no reference path through the snarl, so we can't make a variant
    // (todo: should we try to detect this before computing traversals?)
    if (ref_travs.empty()) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site becuase no reference traversal was found " << pb2json(*snarl) << endl;
#endif
        return false;
    }

    // add in the exhaustive traversals
    if (!path_restricted) {
        // exhaustive traversal can't do all snarls
        if (snarl->type() != ULTRABUBBLE) {
            return false;
        }
        if (!check_max_nodes(snarl)) {
#pragma omp critical (cerr)
            cerr << "Warning: Skipping site because it is too complex for exhaustive traversal enumeration: " << pb2json(*snarl) << endl << "         Consider using -e to traverse embedded paths" << endl;
            return false;
        }
        vector<SnarlTraversal> exhaustive_travs = explicit_exhaustive_traversals(snarl);
        // happens when there was a nested non-ultrabubble snarl
        if (exhaustive_travs.empty()) {
            return false;
        }
        path_travs.first.insert(path_travs.first.end(), exhaustive_travs.begin(), exhaustive_travs.end());
        for (int i = 0; i < exhaustive_travs.size(); ++i) {
            // dummy names so we can use the same code as the named path traversals above
            path_trav_names.push_back(" >>" + std::to_string(i));
            // dummy handles so we can use the same code as the named path traversals above
            path_travs.second.push_back(make_pair(step_handle_t(), step_handle_t()));
            
        }
    }
    
    // there's not alt path through the snarl, so we can't make an interesting variant
    if (path_travs.first.size() < 2) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site because to alt traversal was found " << pb2json(*snarl) << endl;
#endif
        return false;
    }

    // we write a variant for every reference traversal
    for (auto ref_trav_idx : ref_travs) {

        const SnarlTraversal& ref_trav = path_travs.first[ref_trav_idx];
        
        vcflib::Variant v;
        v.setVariantCallFile(outvcf);
        v.quality = 23;

        // write variant's sequenceName (VCF contig)
        v.sequenceName = ref_trav_name;

        // Map our snarl endpoints to oriented positions in the embedded path in the graph
        step_handle_t start_step = path_travs.second[ref_trav_idx].first;
        step_handle_t end_step = path_travs.second[ref_trav_idx].second;
        handle_t start_handle = graph->get_handle_of_step(start_step);
        handle_t end_handle = graph->get_handle_of_step(end_step);
        size_t start_pos = graph->get_position_of_step(start_step);
        size_t end_pos = graph->get_position_of_step(end_step);
        bool use_start = start_pos < end_pos;
        handle_t first_path_handle = use_start ? start_handle : end_handle;
        size_t first_path_pos = use_start ? start_pos : end_pos;
        // Get the first visit of our snarl traversal
        const Visit& first_trav_visit = use_start ? ref_trav.visit(0) : ref_trav.visit(ref_trav.visit_size() - 1);

        char prev_char;
        if ((use_start && first_trav_visit.backward() == graph->get_is_reverse(first_path_handle)) ||
            (!use_start && first_trav_visit.backward() != graph->get_is_reverse(first_path_handle))) {
            // Our path and traversal have consistent orientation.  leave off the end of the start node going forward
            first_path_pos += graph->get_length(first_path_handle);
            prev_char = ::toupper(graph->get_sequence(first_path_handle)[graph->get_length(first_path_handle) - 1]);
        } else {
            // They are flipped: leave off the beginning of the start node going backward
            prev_char = reverse_complement(::toupper(graph->get_sequence(first_path_handle)[0]));
        }
        
        // shift from 0-based to 1-based for VCF
        first_path_pos += 1;

        v.position = first_path_pos;

        v.id = std::to_string(snarl->start().node_id()) + "_" + std::to_string(snarl->end().node_id());
        
        // Convert the snarl traversals to strings and add them to the variant
        vector<int> trav_to_allele = get_alleles(v, path_travs.first, ref_trav_idx, prev_char, use_start);

        // Fill in the genotypes
        if (path_restricted) {
            get_genotypes(v, path_trav_names, trav_to_allele);
        }

        // we only bother printing out sites with at least 1 non-reference allele
        if (!std::all_of(trav_to_allele.begin(), trav_to_allele.end(), [](int i) { return i == 0; })) {
#pragma omp critical (cout)
            {
                cout << v << endl;
            }
        }
    }
    return true;
}

/**
 * Convenience wrapper function for deconstruction of multiple paths.
 */
void Deconstructor::deconstruct(vector<string> ref_paths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                                bool path_restricted_traversals, int ploidy, bool include_nested,
                                const unordered_map<string, string>* path_to_sample) {

    this->graph = graph;
    this->snarl_manager = snarl_manager;
    this->path_restricted = path_restricted_traversals;
    this->ploidy =ploidy;
    this->path_to_sample = path_to_sample;
    this->ref_paths = set<string>(ref_paths.begin(), ref_paths.end());
    assert(path_to_sample == nullptr || path_restricted);
    
    // Keep track of the non-reference paths in the graph.  They'll be our sample names
    sample_names.clear();
    graph->for_each_path_handle([&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (!this->ref_paths.count(path_name)) {
                // rely on the given map.  if a path isn't in it, it'll be ignored
                if (path_to_sample) {
                    if (path_to_sample->count(path_name)) {
                        sample_names.insert(path_to_sample->find(path_name)->second);
                    }
                    // if we have the map, we only consider paths there-in
                }
                else {
                    // no name mapping, just use every path as is
                    sample_names.insert(path_name);
                }
            }
        });
    
    // print the VCF header
    stringstream stream;
    stream << "##fileformat=VCFv4.2" << endl;
    if (path_restricted) {
        stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    }
    if (path_to_sample) {
        stream << "##FORMAT=<ID=PI,Number=.,Type=String,Description=\"Path information. Original vg path name for sample as well as its allele (can be many paths per sample)\">" << endl;
        stream << "##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in the graph with conflicting alleles (details in PI field)\">" << endl;
    }
    for(auto& refpath : ref_paths) {
        size_t path_len = 0;
        path_handle_t path_handle = graph->get_path_handle(refpath);
        for (handle_t handle : graph->scan_path(path_handle)) {
            path_len += graph->get_length(handle);
        }
        stream << "##contig=<ID=" << refpath << ",length=" << path_len << ">" << endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (path_restricted) {
        for (auto& sample_name : sample_names) {
            stream << "\t" << sample_name;
        }
    }
    stream << endl;
    
    string hstr = stream.str();
    assert(outvcf.openForOutput(hstr));
    cout << outvcf.header << endl;

    // create the traversal finder
    map<string, const Alignment*> reads_by_name;
    path_trav_finder = unique_ptr<PathTraversalFinder>(new PathTraversalFinder(*graph,
                                                                               *snarl_manager));
    
    if (!path_restricted) {
        trav_finder = unique_ptr<TraversalFinder>(new ExhaustiveTraversalFinder(*graph,
                                                                                *snarl_manager,
                                                                                true));

    }
    
    // Do the top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
            vector<const Snarl*> todo(1, snarl);
            vector<const Snarl*> next;
            while (!todo.empty()) {
                for (auto next_snarl : todo) {
                    // if we can't make a variant from the snarl due to not finding
                    // paths through it, we try again on the children
                    // note: we may want to push the parallelism down a bit 
                    if (!deconstruct_site(next_snarl) || include_nested) {
                        const vector<const Snarl*>& children = snarl_manager->children_of(next_snarl);
                        next.insert(next.end(), children.begin(), children.end());
                    }
                }
                swap(todo, next);
                next.clear();
            }
        });
}

bool Deconstructor::check_max_nodes(const Snarl* snarl)  {
    unordered_set<id_t> nodeset = snarl_manager->deep_contents(snarl, *graph, false).first;
    int node_count = 0;
    for (auto node_id : nodeset) {
        handle_t node = graph->get_handle(node_id);
        if (graph->get_degree(node, true) > 1 || graph->get_degree(node, false) > 1) {
            ++node_count;
            if (node_count > max_nodes_for_exhaustive) {
                return false;
            }
        }
    }
    return true;
};

vector<SnarlTraversal> Deconstructor::explicit_exhaustive_traversals(const Snarl* snarl){
    vector<SnarlTraversal> out_travs;
    bool ultra_all_the_way_down = true;
    function<void(const SnarlTraversal&, const Snarl&)> extend_trav =
        [&](const SnarlTraversal& trav, const Snarl& nested_snarl) {
        // exhaustive traversal finder is limited.  if we find something
        // that's not an ultrabubble, not much we can do
        if (nested_snarl.type() != ULTRABUBBLE) {
            ultra_all_the_way_down = false;
            return;
        }
        vector<SnarlTraversal> nested_travs = trav_finder->find_traversals(nested_snarl);
        for (auto& nested_trav : nested_travs) {
            SnarlTraversal extended_trav = trav;
            bool is_explicit = true;
            for (int i = 0; i < nested_trav.visit_size(); ++i) {
                if (nested_trav.visit(i).node_id() != 0) {
                    Visit* visit = extended_trav.add_visit();
                    *visit = nested_trav.visit(i);
                } else {
                    extend_trav(extended_trav, nested_trav.visit(i).snarl());
                    is_explicit = false;
                }
            }
            if (is_explicit) {
                out_travs.push_back(extended_trav);
            }
        }
    };
    SnarlTraversal trav;
    extend_trav(trav, *snarl);
    if (!ultra_all_the_way_down) {
        out_travs.clear();
    }        
    return out_travs;
}

}

