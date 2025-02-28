#include "deconstructor.hpp"
#include "traversal_finder.hpp"
#include <gbwtgraph/gbwtgraph.h>
#include "traversal_clusters.hpp"

//#define debug

using namespace std;


namespace vg {
Deconstructor::Deconstructor() : VCFOutputCaller("") {
}
Deconstructor::~Deconstructor(){
}

/**
 * Takes in a vector of snarltraversals, an index of the ref path among them, a
 * vector of flags for traversals to actually use, the character before all the
 * traversals, and a flag for whether the start should be used (???).
 *
 * Returns a vector where entires are which allele number a traversal in travs
 * ought to become in the VCF. If a traversal is flagged off, it gets a -1.
 */
vector<int> Deconstructor::get_alleles(vcflib::Variant& v,
                                       const vector<Traversal>& travs,
                                       const vector<pair<step_handle_t, step_handle_t>>& trav_steps,
                                       int ref_path_idx,
                                       const vector<vector<int>>& trav_clusters,
                                       char prev_char, bool use_start) const {

    assert(ref_path_idx >=0 && ref_path_idx < travs.size());

    // map strings to allele numbers (and their traversal)
    // (we are using the traversal finder in such a way that duplicate alleles can get returned
    // in order to be able to preserve the path names)
    map<string, pair<int, int>> allele_idx;
    size_t cur_alt = 1;

    // go from traversals number (offset in travs) to allele number
    vector<int> trav_to_allele(travs.size());

    // compute the allele as a string
    auto trav_to_string = [&](const Traversal& trav) {
        string allele;
        // hack to support star alleles
        if (trav.size() == 0) {
            allele = "*";
        } else {
            // we skip the snarl endpoints
            for (int j = 1; j < trav.size() - 1; ++j) {
                allele += toUppercase(graph->get_sequence(trav[j]));
            }
        }
        return allele;
    };

    // set the reference allele
    string ref_allele = trav_to_string(travs.at(ref_path_idx));
    allele_idx[ref_allele] = make_pair(0, ref_path_idx);
    trav_to_allele[ref_path_idx] = 0;
    bool substitution = true;
        
    // set the other alleles (they can end up as 0 alleles too if their strings match the reference)
    // note that we have one (unique) allele per cluster, so we take advantage of that here
    for (const vector<int>& cluster : trav_clusters) {
        string allele = trav_to_string(travs[cluster.front()]);
        for (const int& i : cluster) {
            if (i != ref_path_idx) {
                auto ai_it = allele_idx.find(allele);
                if (ai_it == allele_idx.end()) {
                    // make a new allele for this string
                    allele_idx[allele] = make_pair(cur_alt, i);
                    trav_to_allele.at(i) = cur_alt;
                    ++cur_alt;
                    substitution = substitution && allele.size() == ref_allele.size();
                } else {
                    // allele string has been seen, map this traversal to it
                    trav_to_allele.at(i) = ai_it->second.first;
                }
            } else {
                trav_to_allele.at(i) = -1; // HACK! negative allele indexes are ignored
            }
        }
    }

    // fill in the variant
    v.alleles.resize(allele_idx.size());
    assert(allele_idx.size() > 0);
    v.alt.resize(allele_idx.size() - 1);

    // if we should flip the traversals
    bool reversed = !use_start;
    //if (reversed) {
    //    cerr << "it's reversed!!!" << endl;
    //}

    // we're going to go through the allele traversals
    // for each, we should find the set of supporting paths
    // we will use them to untangle the reference positions
    // using the same context mapping typically used for record positions
    // AP field (allele positions)
    // same pattern as AT except the steps are listed [>|<][id]_[start]_[end]+
    // id is the node id and start/end give the 1-based, half-open reference coordinates of the node
    // we establish these by making a reference context per node in the reference allele through the site
    // for each step where the reference touches the node more than once
    // (in cases where it's only once, this is trivial)
    // then we compare each reference matching step inside
    // (todo: this should be done only for top level bubbles)

    vector<int> allele_idx_unfolded(allele_idx.size());
    for (auto& ai_pair : allele_idx) {
        int allele_no = ai_pair.second.first;
        int allele_trav_no = ai_pair.second.second;
        allele_idx_unfolded[allele_no] = allele_trav_no;
    }

    // record the alleles in the VCF record
    for (auto ai_pair : allele_idx) {
        string allele_string = ai_pair.first;
        int allele_no = ai_pair.second.first;
        int allele_trav_no = ai_pair.second.second;
        if (reversed) {
            reverse_complement_in_place(allele_string);
        }
        if (!substitution) {
            allele_string = string(1, prev_char) + allele_string;
        }
        v.alleles[allele_no] = allele_string;
        if (allele_no > 0) {
            v.alt[allele_no - 1] = allele_string;
        } else {
            v.ref = allele_string;
        }
    }

    if (untangle_allele_traversals) {

        // set up for reference position context mapping across allele traversals
        path_handle_t ref_path = graph->get_path_handle_of_step(trav_steps.at(ref_path_idx).first);
        unordered_map<nid_t, vector<pair<uint64_t, dac_vector<>>>> ref_dup_nodes;
        unordered_map<nid_t, nid_t> ref_simple_pos;
        {
            auto& trav = travs.at(ref_path_idx);
            for (size_t i = 0; i < trav.size(); ++i) {
                size_t j = !reversed ? i : trav.size() - 1 - i;
                const handle_t& handle = trav[j];
                nid_t node_id = graph->get_id(handle);
                if (ref_simple_pos.find(node_id) != ref_simple_pos.end()) continue;
                if (ref_dup_nodes.find(node_id) != ref_dup_nodes.end()) continue;
                handle_t h = graph->get_handle(node_id);
                // count reference occurrences on node
                step_handle_t ref_step;
                uint64_t ref_step_count = 0;
                graph->for_each_step_on_handle(
                    h, [&](const step_handle_t& step) {
                        auto p = graph->get_path_handle_of_step(step);
                        if (p == ref_path) {
                            ++ref_step_count;
                            ref_step = step;
                        }
                    });
                if (ref_step_count > 1) {
                    //ref_dup_nodes[node_id] = make_pair(ref_context);
                    ref_dup_nodes[node_id] = {};
                    auto& contexts = ref_dup_nodes[node_id];
                    //vector<pair<uint64_t, vector<nid_t>>> contexts;
                    graph->for_each_step_on_handle(
                        h, [&](const step_handle_t& step) {
                            auto p = graph->get_path_handle_of_step(step);
                            if (p == ref_path) {
                                contexts.emplace_back();
                                auto& c = contexts.back();
                                c.first = graph->get_position_of_step(step);
                                c.second = get_context(step, step);
                            }
                        });
                    //
                } else if (ref_step_count == 1) {
                    auto pos = graph->get_position_of_step(ref_step) + 1;
                    auto len = graph->get_length(graph->get_handle_of_step(ref_step));
                    ref_simple_pos[node_id] = pos;
                }
            }
        }

        // set up the UT field for our reference-relative traversal position untangling
        auto& ut_field = v.info["UT"];
        ut_field.resize(allele_idx.size());

#pragma omp parallel for schedule(dynamic,1)
        for (size_t i = 0; i < allele_idx_unfolded.size(); i++) {
            int allele_no = i;
            int allele_trav_no = allele_idx_unfolded[i];
            auto start_step = trav_steps.at(allele_trav_no).first;
            auto end_step = trav_steps.at(allele_trav_no).second;
            auto start_pos = graph->get_position_of_step(start_step);
            auto end_pos = graph->get_position_of_step(end_step);
            bool flip_path = start_pos > end_pos;
            if (flip_path) {
                std::swap(start_step, end_step);
            }
            path_handle_t path = graph->get_path_handle_of_step(start_step);
            std::vector<step_handle_t> steps;
            for (auto s = start_step; ;
                 s = graph->get_next_step(s)) {
                steps.push_back(s);
                if (s == end_step) break;
                if (!graph->has_next_step(s)) break;
            }
            if (steps.front() != start_step || steps.back() != end_step) {
                //cerr << "warning!" << endl;
                // something went wrong
                ut_field[allele_no] = ".";
                continue;
            }
            if (flip_path) {
                std::reverse(steps.begin(), steps.end());
            }
            // update the traversal info
            stringstream trav_pos_info;
            //string& trav_pos_info = ut_field[allele_no];
            for (auto& step : steps) {
                handle_t h = graph->get_handle_of_step(step);
                nid_t node_id = graph->get_id(h);
                bool step_rev = graph->get_is_reverse(h) != reversed;
                trav_pos_info << (reversed ? "<" : ">") << node_id << "_";
                if (allele_no == 0) { // special case the reference allele
                    auto pos = graph->get_position_of_step(step) + 1;
                    auto len = graph->get_length(graph->get_handle_of_step(step));
                    trav_pos_info << pos << "_" << pos+len;
                } else { // for non-reference alleles
                    auto f = ref_simple_pos.find(node_id);
                    if (f != ref_simple_pos.end()) {
                        // we have a single reference position at this node
                        auto pstart = f->second + 1;
                        auto pend = pstart + graph->get_length(h);
                        trav_pos_info << pstart << "_" << pend;
                    } else {
                        auto d = ref_dup_nodes.find(node_id);
                        if (d == ref_dup_nodes.end()) {
                            // no reference position at this node
                            trav_pos_info << "._.";
                        } else {
                            // multiple reference positions at this node
                            // compare the reference contexts of each step to our
                            // path context to determine reference position assignment
                            // ... first we get our path's context
                            auto path_context = get_context(step, step);
                            auto& ref_contexts = d->second;
                            //cerr << "path context size " << path_context.size() << endl;
                            // check vs. the contexts
                            double best_jaccard = -1;
                            uint64_t best_pos = 0;
                            for (auto& c : ref_contexts) {
                                auto& ref_context = c.second;
                                auto& ref_pos = c.first;
                                double j = jaccard_coefficient(ref_context, path_context);
                                if (j > best_jaccard) {
                                    best_jaccard = j;
                                    best_pos = ref_pos;
                                }
                            }
                            trav_pos_info << best_pos+1 << "_" << best_pos+1+graph->get_length(h);
                        }
                    }
                }
            }
            // save the untangled traversal field
            ut_field[allele_no] = trav_pos_info.str();
        }
    } else {
        // init the traversal info
        v.info["AT"].resize(allele_idx.size());
        for (auto ai_pair : allele_idx) {
            string allele_string = ai_pair.first;
            int allele_no = ai_pair.second.first;
            int allele_trav_no = ai_pair.second.second;
            // update the traversal info
            add_allele_path_to_info(graph, v, allele_no, travs[allele_trav_no], reversed, !substitution);
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
                                  const vector<int>& trav_to_allele,
                                  const vector<pair<double, int64_t>>& trav_to_cluster_info) const {  
    assert(names.size() == trav_to_allele.size());
    // set up our variant fields
    v.format.push_back("GT");
    if (show_path_info && path_to_sample_phase) {
        v.format.push_back("PI");
    }
    if (this->cluster_threshold < 1.0) {
        v.format.push_back("TS");
        v.format.push_back("TL");
    }
    
    // get a list of traversals for every vcf sample
    // (this will be 1:1 unless we're using the path_to_sample name map)
    map<string, vector<int> > sample_to_traversals;
    // phasing information from the gbwt where applicable
    vector<int> gbwt_phases(trav_to_allele.size(), -1);
    for (int i = 0; i < names.size(); ++i) {
        string sample_name = PathMetadata::parse_sample_name(names[i]);
        // for backward compatibility
        if (sample_name.empty()) {
            sample_name = names[i];
        }
        auto phase = PathMetadata::parse_haplotype(names[i]);
        if (!sample_name.empty() && phase == PathMetadata::NO_HAPLOTYPE) {
            // THis probably won't fit in an int. Use 0 instead.
            phase = 0;
        }
        gbwt_phases[i] = (int)phase;
        if (sample_names.count(sample_name)) {
            sample_to_traversals[sample_name].push_back(i);
        }
    }

    // write out the genotype for each sample
    // if we're mapping a vg path name to its prefix for the sample name, we stick some information about the full
    // path name in the PI part of format
    set<string> conflicts;
    for (const auto& sample_name : sample_names) {
        if (sample_to_traversals.count(sample_name)) {
            const vector<int>& travs = sample_to_traversals[sample_name];
            assert(!travs.empty());
            vector<int> chosen_travs;
            bool conflict;
            std::tie(chosen_travs, conflict) = choose_traversals(sample_name, travs, trav_to_allele, names, gbwt_phases);
            if (conflict) {
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Conflict for sample " << sample_name << endl;
#endif
                conflicts.insert(sample_name);
            }
            string genotype;
            for (int i = 0; i < chosen_travs.size(); ++i) {
                if (i > 0) {
                    // TODO check flag for phasing
                    genotype += (path_to_sample_phase || gbwt_trav_finder.get())  ? "|" : "/";
                }
                genotype += (chosen_travs[i] != -1 && (!conflict || keep_conflicted_genotypes))
                    ? std::to_string(trav_to_allele[chosen_travs[i]]) : ".";
                if (this->cluster_threshold < 1.0) {
                    if (*genotype.rbegin() == '.') {
                        v.samples[sample_name]["TS"].push_back(".");
                        v.samples[sample_name]["TL"].push_back(".");
                    } else {
                        ostringstream ss;
                        ss.precision(3);
                        ss << trav_to_cluster_info[chosen_travs[i]].first;
                        v.samples[sample_name]["TS"].push_back(ss.str());
                        v.samples[sample_name]["TL"].push_back(std::to_string(trav_to_cluster_info[chosen_travs[i]].second));
                    }
                }
            }
            v.samples[sample_name]["GT"] = {genotype};
            if (show_path_info && path_to_sample_phase) {
                for (auto trav : travs) {
                    auto allele = trav_to_allele[trav];
                    if (allele != -1) {
                        v.samples[sample_name]["PI"].push_back(names[trav] + "=" + std::to_string(allele));
                    }
                }
            }
        } else {
            string blank_gt = ".";
            if (gbwt_sample_to_phase_range.count(sample_name)) {
                // note: this is the same logic used when filling in actual genotypes
                int min_phase, max_phase;
                std::tie(min_phase, max_phase) = gbwt_sample_to_phase_range.at(sample_name);
                // shift left by 1 unless min phase is 0
                int sample_ploidy = min_phase == 0 ? max_phase + 1 : max_phase;
                assert(sample_ploidy > 0);
                for (int phase = 1; phase < sample_ploidy; ++phase) {
                    blank_gt += "|.";
                }
            }
            v.samples[sample_name]["GT"] = {blank_gt};
            if (show_path_info && path_to_sample_phase) {
                v.samples[sample_name]["PI"] = {blank_gt};
            }
        }
    }
    for (auto& conflict_sample : conflicts) {
        v.info["CONFLICT"].push_back(conflict_sample);
    }
}

pair<vector<int>, bool> Deconstructor::choose_traversals(const string& sample_name,
                                                         const vector<int>& travs, const vector<int>& trav_to_allele,
                                                         const vector<string>& trav_to_name,
                                                         const vector<int>& gbwt_phases) const {

    assert(trav_to_name.size() == trav_to_allele.size());
    assert(gbwt_phases.size() == trav_to_name.size());
    assert(!travs.empty());
    // count the number of times each allele comes up in a traversal
    vector<int> allele_frequencies(std::max(0, *max_element(trav_to_allele.begin(), trav_to_allele.end()) + 1), 0);
    for (auto trav : travs) {
        // we always want to choose alt over ref when possible in sorting logic below, so
        // cap ref frequency at 1
        int allele = trav_to_allele.at(trav);
        // avoid counting filtered alleles (== -1)
        if (allele >= 0 && (allele > 0 || allele_frequencies[allele] == 0)) {
            ++allele_frequencies[allele];
        }
    }
    // sort on frquency
    function<bool(int, int)> comp = [&] (int trav1, int trav2) {
        auto& trav1_allele = trav_to_allele[trav1];
        auto& trav2_allele = trav_to_allele[trav2];
        // avoid the filtered allele state (== -1)
        auto trav1_af = trav1_allele >= 0 ? allele_frequencies[trav1_allele] : 0;
        auto trav2_af = trav2_allele >= 0 ? allele_frequencies[trav2_allele] : 0;
        if (trav1_af < trav2_af) {
            return true;
        } else if (trav1_af == trav2_af) {
            // prefer non-ref when possible
            if (trav1_allele < trav2_allele) {
                return true;
            } else if (trav1_allele == trav2_allele) {
                return trav_to_name[trav1] < trav_to_name[trav2];
            }
        }
        return false;
    };
    vector<int> sorted_travs = travs;
    std::sort(sorted_travs.begin(), sorted_travs.end(), comp);
    // find the <ploidy> most frequent traversals
    vector<int> most_frequent_travs;

    // try to pull out unique phases if available
    bool has_phasing = gbwt_sample_to_phase_range.count(sample_name) &&
        std::any_of(gbwt_phases.begin(), gbwt_phases.end(), [](int i) { return i >= 0; });
    //|| path_to_sample_phase;
    bool phasing_conflict = false;
    int sample_ploidy = 1;
    int min_phase = 1;
    int max_phase = 1;
    if (has_phasing || path_to_sample_phase) {
        if (has_phasing) {
            // override ploidy with information about all phases found in input
            std::tie(min_phase, max_phase) = gbwt_sample_to_phase_range.at(sample_name);
            // shift left by 1 unless min phase is 0
            sample_ploidy = min_phase == 0 ? max_phase + 1 : max_phase;
            assert(sample_ploidy > 0);
        } else {
            sample_ploidy = sample_ploidys.at(sample_name);
        }        
        set<int> used_phases;
        for (int i = sorted_travs.size() - 1; i >= 0 && most_frequent_travs.size() < sample_ploidy; --i) {
            int phase = gbwt_phases.at(sorted_travs[i]);
            if (!used_phases.count(phase)) {
                if (trav_to_allele[sorted_travs[i]] >= 0) {
                    most_frequent_travs.push_back(sorted_travs[i]);
                    used_phases.insert(phase);
                }
            } else if (strict_conflict_checking) {
                phasing_conflict = true;
            }
        }
    } else {
        for (int i = sorted_travs.size() - 1; i >= 0 && most_frequent_travs.size() < sample_ploidy; --i) {
            if (trav_to_allele[sorted_travs[i]] >= 0) {
                most_frequent_travs.push_back(sorted_travs[i]);
            }
        }
    }

    // sort by phase
    if (has_phasing) {
        std::sort(most_frequent_travs.begin(), most_frequent_travs.end(),
                  [&](int t1, int t2) {return gbwt_phases.at(t1) < gbwt_phases.at(t2);});
        if (max_phase > 0) {
            // pad out by phase
            assert(most_frequent_travs.empty() || gbwt_phases.at(most_frequent_travs.back()) <= max_phase);
            assert(max_phase < 1000);
            // we normally expect to have phases 1,2,3, ...
            // in this case, we shift them all back, otherwise leave 0-based
            int offset = min_phase != 0 ? -1 : 0;
            vector<int> padded_travs(max_phase + 1 + offset, -1);
            for (auto ft : most_frequent_travs) {
                int phase = gbwt_phases.at(ft) + offset;
                padded_travs.at(phase) = ft;
            }
            swap(padded_travs, most_frequent_travs);
        }
    } else if (path_to_sample_phase) {
        std::sort(most_frequent_travs.begin(), most_frequent_travs.end(),
                  [&](int t1, int t2) {return gbwt_phases.at(t1) < gbwt_phases.at(t2);});
        vector<int> padded_travs(sample_ploidy, -1);
        for (auto ft : most_frequent_travs) {
            int phase = gbwt_phases.at(ft);
            //std::cerr << "on phase " << phase << std::endl;
            padded_travs.at(phase) = ft;
        }
        swap(padded_travs, most_frequent_travs);
    }
    // check if there's a conflict
    size_t zero_count = std::count(allele_frequencies.begin(), allele_frequencies.end(), 0);
    bool conflict = phasing_conflict || allele_frequencies.size() - zero_count > sample_ploidy;
    return make_pair(most_frequent_travs, conflict);
}

vector<nid_t> Deconstructor::get_context(
    step_handle_t start_step,
    step_handle_t end_step) const {
    if (graph->get_position_of_step(start_step)
        > graph->get_position_of_step(end_step)) {
        std::swap(start_step, end_step);
    }
    // by definition, our start and end are shared among all traversals
    // we establish a graph context sketch including the nodes traversed in the bubble
    // and flanks upstream and downstream of path_jaccard_window bp
    vector<nid_t> context;
    step_handle_t curr = start_step;
    const int check_distance = this->path_jaccard_window; // how far we look in each direction
    int distance_checked = 0;
    while (distance_checked < check_distance && graph->has_previous_step(curr)) {
        curr = graph->get_previous_step(curr);
        auto h = graph->get_handle_of_step(curr);
        context.push_back(graph->get_id(h));
        distance_checked += graph->get_length(h);
    }
    // add the nodes in the bubble
    if (start_step != end_step) {
        curr = start_step;
        context.push_back(graph->get_id(graph->get_handle_of_step(curr)));
        while (graph->has_next_step(curr) &&
               curr != end_step) {
            curr = graph->get_next_step(curr);
            context.push_back(graph->get_id(graph->get_handle_of_step(curr)));
        }
    }
    distance_checked = 0;
    curr = end_step;
    while (distance_checked < check_distance && graph->has_next_step(curr)) {
        curr = graph->get_next_step(curr);
        auto h = graph->get_handle_of_step(curr);
        context.push_back(graph->get_id(h));
        distance_checked += graph->get_length(h);
    }
    std::sort(context.begin(), context.end());
    return context;
}


void Deconstructor::get_traversals(const handle_t& snarl_start, const handle_t& snarl_end,
                                   vector<Traversal>& out_travs,
                                   vector<string>& out_trav_path_names,
                                   vector<pair<step_handle_t, step_handle_t>>& out_trav_steps) const {
    // empty snarl check
    vector<handle_t> next_handles;
    graph->follow_edges(snarl_start, false, [&](handle_t handle) {
            next_handles.push_back(handle);
        });
    if (next_handles.size() == 1 && next_handles.back() == snarl_end) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping empty site " << graph_interval_to_string(graph, snarl_start, snarl_end) << endl;
#endif        
        return;
    }
    
#ifdef debug
#pragma omp crtiical (cerr)
    cerr << "Computing traversals of site " << graph_interval_to_string(graph, snarl_start, snarl_end) << endl;
#endif

    // find every traversal that runs through a path in the graph
    std::tie(out_travs, out_trav_steps) = path_trav_finder->find_path_traversals(snarl_start, snarl_end);
    for (const pair<step_handle_t, step_handle_t>& trav_steps : out_trav_steps) {
        out_trav_path_names.push_back(graph->get_path_name(graph->get_path_handle_of_step(trav_steps.first)));
    }
    
    // add in the gbwt traversals
    // after this, all traversals are treated the same, with metadata embedded in their names
    if (gbwt_trav_finder.get() != nullptr) {
        const gbwt::GBWT& gbwt_index = gbwt_trav_finder->get_gbwt();
        pair<vector<Traversal>, vector<gbwt::size_type>> thread_travs = gbwt_trav_finder->find_path_traversals(snarl_start, snarl_end);
        for (int i = 0; i < thread_travs.first.size(); ++i) {
            // We need to get a bunch of metadata about the path, but the GBWT
            // we have might not even have structured path names stored.
            gbwt::size_type path_id = gbwt::Path::id(thread_travs.second[i]);
            if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || path_id >= gbwt_index.metadata.paths()) {
                continue;
            }
            
            PathSense sense = gbwtgraph::get_path_sense(gbwt_index, path_id, gbwt_reference_samples);
            
            if (sense == PathSense::HAPLOTYPE) {
                // we count on convention of reference as embedded path above, so only use haplotype paths here.
                // todo: would be nice to be more flexible...
                string path_name = PathMetadata::create_path_name(
                    sense,
                    gbwtgraph::get_path_sample_name(gbwt_index, path_id, sense),
                    gbwtgraph::get_path_locus_name(gbwt_index, path_id, sense),
                    gbwtgraph::get_path_haplotype(gbwt_index, path_id, sense),
                    gbwtgraph::get_path_phase_block(gbwt_index, path_id, sense),
                    gbwtgraph::get_path_subrange(gbwt_index, path_id, sense));
                out_trav_path_names.push_back(path_name);
                out_travs.push_back(std::move(thread_travs.first[i]));
            }
        }
    }
}

unordered_map<string, vector<int>> Deconstructor::add_star_traversals(vector<Traversal>& travs,
                                                                      vector<string>& names,
                                                                      vector<vector<int>>& trav_clusters,
                                                                      vector<pair<double, int64_t>>& trav_cluster_info,
                                                                      const unordered_map<string, vector<int>>& parent_haplotypes) const {    
    // todo: refactor this into general genotyping code
    unordered_map<string, vector<int>> sample_to_haps;

    // find out what's in the traversals
    assert(names.size() == travs.size());
    for (int64_t i = 0; i < names.size(); ++i) {
        string sample_name = PathMetadata::parse_sample_name(names[i]);
        // for backward compatibility
        if (sample_name.empty()) {
            sample_name = names[i];
        }
        auto phase = PathMetadata::parse_haplotype(names[i]);
        if (!sample_name.empty() && phase == PathMetadata::NO_HAPLOTYPE) {
            // THis probably won't fit in an int. Use 0 instead.
            phase = 0;
        }
        sample_to_haps[sample_name].push_back(phase);
    }

    // find everything that's in parent_haplotyes but not the travefsals,
    // and add in dummy start-alleles for them
    for (const auto& parent_sample_haps : parent_haplotypes) {
        string parent_sample_name = PathMetadata::parse_sample_name(parent_sample_haps.first);
        if (parent_sample_name.empty()) {
            parent_sample_name = parent_sample_haps.first;
        }
        if (!this->sample_names.count(parent_sample_name)) {
            // dont' bother for purely reference samples -- we don't need to force and allele for them.
            continue;
        }
        for (int parent_hap : parent_sample_haps.second) {
            bool found = false;
            if (sample_to_haps.count(parent_sample_haps.first)) {
                // note: this is brute-force search, but number of haplotypes usually tiny.
                for (int hap : sample_to_haps[parent_sample_haps.first]) {
                    if (parent_hap == hap) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                travs.push_back(Traversal());
                names.push_back(PathMetadata::create_path_name(PathSense::REFERENCE,
                                                               parent_sample_haps.first,
                                                               "star",
                                                               parent_hap,
                                                               PathMetadata::NO_PHASE_BLOCK,
                                                               PathMetadata::NO_SUBRANGE));
                sample_to_haps[parent_sample_haps.first].push_back(parent_hap);
                trav_clusters.push_back({(int)travs.size() - 1});
                trav_cluster_info.push_back(make_pair(0, 0));
            }
        }
    }
    
    return sample_to_haps;
}


bool Deconstructor::deconstruct_site(const handle_t& snarl_start, const handle_t& snarl_end,
                                     const NestingInfo* in_nesting_info,
                                     vector<NestingInfo>* out_nesting_infos) const {


    vector<Traversal> travs;
    vector<string> trav_path_names;
    // note that this vector (unlike the above two) is for embedded paths only (not GBWT)
    vector<pair<step_handle_t, step_handle_t>> trav_steps;

    // compute all the traversals from embedded paths and gbwt
    this->get_traversals(snarl_start, snarl_end, travs, trav_path_names, trav_steps);
    int64_t trav_count = travs.size();
    int64_t trav_step_count = trav_steps.size();

    if (travs.empty()) {
        return false;        
    }
    
    // pick out the traversal corresponding to an embedded reference path, breaking ties consistently
    string ref_trav_name;
    string parent_ref_trav_name;
    if (in_nesting_info != nullptr && in_nesting_info->has_ref) {
        parent_ref_trav_name = graph->get_path_name(graph->get_path_handle_of_step(in_nesting_info->parent_path_interval.first));
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Using nesting information to set reference to " << parent_ref_trav_name << endl;
#endif
        // remember it for the vcf header
        this->off_ref_paths[omp_get_thread_num()].insert(graph->get_path_handle_of_step(in_nesting_info->parent_path_interval.first));
    }
    for (int i = 0; i < travs.size(); ++i) {
        const string& path_trav_name = trav_path_names[i];
#ifdef debug
#pragma omp critical (cerr)
        {
            cerr << "Traversal " << i << ": name=" << path_trav_name << ", size=" << travs[i].size();
            if (i < trav_steps.size()) {
                cerr << ", start=" << graph->get_position_of_step(trav_steps[i].first)
                     << ", end=" << graph->get_position_of_step(trav_steps[i].second) << endl;
            }
            cerr << " trav=" << traversal_to_string(graph, travs[i]) << endl;
        }
#endif
        bool ref_path_check;
        if (!parent_ref_trav_name.empty()) {
            // the reference was specified by the parent
            ref_path_check = path_trav_name == parent_ref_trav_name;
        } else {
            // the reference comes from the global options
            ref_path_check = ref_paths.count(path_trav_name);
        }
        if (ref_path_check &&
            (ref_trav_name.empty() || path_trav_name < ref_trav_name)) {
            ref_trav_name = path_trav_name;
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Setting ref_trav_name " << ref_trav_name << (in_nesting_info ? " using nesting info" : "") << endl;
#endif
        }
    }
    
    // remember all the reference traversals (there can be more than one only in the case of a
    // cycle in the reference path

    // in case of cycles, we need our allele traversals to be associated to the correct reference position
    // this is done with the path jaccard metric over all overlapping reference paths the given path_jaccard_window size
    vector<int> ref_travs;
    // hacky subpath support -- gets added to variant on output
    vector<int64_t> ref_offsets;
    if (!ref_trav_name.empty()) {
        for (int i = 0; i < trav_steps.size(); ++i) {
            const string& path_trav_name = trav_path_names.at(i);
            subrange_t subrange ;
            Paths::strip_subrange(path_trav_name, &subrange);
            int64_t sub_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
            if (path_trav_name == ref_trav_name) {
                ref_travs.push_back(i);
                ref_offsets.push_back(sub_offset);
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Adding ref_trav idx=" << i << " offset=" << sub_offset << " because " << path_trav_name << " == " << ref_trav_name << endl;
#endif                
            }
        }
    }

    // there's no reference path through the snarl, so we can't make a variant
    // (todo: should we try to detect this before computing traversals?)
    if (ref_travs.empty()) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site because no reference traversal was found " << graph_interval_to_string(graph, snarl_start, snarl_end) << endl;
#endif
        return false;
    }    
    
    // there's not alt path through the snarl, so we can't make an interesting variant
    if (travs.size() < 2) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site because to alt traversal was found " << graph_interval_to_string(graph, snarl_start, snarl_end) << endl;
#endif
        return false;
    }

    if (ref_travs.size() > 1 && this->nested_decomposition) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Multiple ref traversals not yet supported with nested decomposition: removing all but first" << endl;
#endif
        size_t min_start_pos = numeric_limits<size_t>::max();
        int64_t first_ref_trav;
        for (int64_t i = 0; i < ref_travs.size(); ++i) {            
            auto& ref_trav_idx = ref_travs[i];
            step_handle_t start_step = trav_steps[ref_trav_idx].first;
            step_handle_t end_step = trav_steps[ref_trav_idx].second;
            size_t ref_trav_pos = min(graph->get_position_of_step(start_step), graph->get_position_of_step(end_step));
            if (ref_trav_pos < min_start_pos) {
                min_start_pos = ref_trav_pos;
                first_ref_trav = i;
            }
        }
        ref_travs = {ref_travs[first_ref_trav]};        
    }

    // XXX CHECKME this assumes there is only one reference path here, and that multiple traversals are due to cycles
    
    // we collect windows around the reference traversals
    // to compare with equivalent windows from the alternate allele paths
    // we will associate these 1:1 with reference traversals

    // remember that path_travs := pair<vector<Traversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs;

    // map from each path_trav index to the ref_trav index it best maps to
    vector<int> path_trav_to_ref_trav;

    if (ref_travs.size() > 1 && this->path_jaccard_window) {
        path_trav_to_ref_trav.resize(trav_steps.size());
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Multiple ref traversals!" << endl;
#endif
        vector<vector<nid_t>> ref_contexts(ref_travs.size());
#pragma omp parallel for schedule(dynamic,1)
        for (size_t i = 0; i < ref_travs.size(); ++i) {
            ref_contexts[i] = get_context(trav_steps[ref_travs[i]].first, trav_steps[ref_travs[i]].second);
        }

        // now for each traversal, we compute and equivalent context and match it to a ref context
        // using a jaccard metric over node ids
#pragma omp parallel for schedule(dynamic,1)
        for (size_t i = 0; i < trav_steps.size(); ++i) {
            vector<nid_t> context = get_context(trav_steps[i].first, trav_steps[i].second);
            // map jaccard metric to the index of the ref_trav
            vector<pair<double, int>> ref_mappings;
            for (uint64_t j = 0; j < ref_travs.size(); ++j) {
                ref_mappings.push_back(make_pair(
                                           jaccard_coefficient(
                                               ref_contexts[j],
                                               context),
                                           ref_travs[j]));
            }
            std::stable_sort(ref_mappings.begin(), ref_mappings.end());
            // the best is the last, which has the highest jaccard
            path_trav_to_ref_trav[i] = ref_mappings.back().second;
        }
    }

    // we write a variant for every reference traversal
    // (optionally) selecting the subset of path traversals that are 1:1
    for (size_t i = 0; i < ref_travs.size(); ++i) {
        // we zap these to their original size, as the nesting logic can add
        // dummy traversals and these are reference-specific (and so need to be cleaned each iteration here)
        travs.resize(trav_count);
        trav_path_names.resize(trav_count);
        trav_steps.resize(trav_step_count);
        auto& ref_trav_idx = ref_travs[i];
        auto& ref_trav_offset = ref_offsets[i];

        const Traversal& ref_trav = travs[ref_trav_idx];

        vcflib::Variant v;
        v.quality = 60;

        // in VCF we usually just want the contig
        string contig_name = PathMetadata::parse_locus_name(ref_trav_name);
        if (contig_name == PathMetadata::NO_LOCUS_NAME) {
            contig_name = ref_trav_name;
        } else if (long_ref_contig) {
            // the sample name isn't unique enough, so put a full ugly name in the vcf
            if (PathMetadata::parse_sense(ref_trav_name) == PathSense::GENERIC) {
                contig_name = ref_trav_name;
            } else {
                contig_name = PathMetadata::create_path_name(PathSense::REFERENCE,
                                                             PathMetadata::parse_sample_name(ref_trav_name),
                                                             contig_name,
                                                             PathMetadata::parse_haplotype(ref_trav_name),
                                                             PathMetadata::NO_PHASE_BLOCK,
                                                             PathMetadata::NO_SUBRANGE);
            }
        }
            
        // write variant's sequenceName (VCF contig)
        v.sequenceName = contig_name;

        // Map our snarl endpoints to oriented positions in the embedded path in the graph
        handle_t first_path_handle;
        size_t first_path_pos;
        bool use_start;
        step_handle_t start_step = trav_steps[ref_trav_idx].first;
        step_handle_t end_step = trav_steps[ref_trav_idx].second;
        handle_t start_handle = graph->get_handle_of_step(start_step);
        handle_t end_handle = graph->get_handle_of_step(end_step);
        size_t start_pos = graph->get_position_of_step(start_step);
        size_t end_pos = graph->get_position_of_step(end_step);
        use_start = start_pos < end_pos;
        first_path_handle = use_start ? start_handle : end_handle;
        first_path_pos = use_start ? start_pos : end_pos;
            
        // Get the first visit of our snarl traversal
        const handle_t& first_trav_handle = use_start ? ref_trav.front() : ref_trav.back();

        char prev_char;
        if ((use_start && graph->get_is_reverse(first_trav_handle) == graph->get_is_reverse(first_path_handle)) ||
            (!use_start && graph->get_is_reverse(first_trav_handle) != graph->get_is_reverse(first_path_handle))) {
            // Our path and traversal have consistent orientation.  leave off the end of the start node going forward
            first_path_pos += graph->get_length(first_path_handle);
            prev_char = ::toupper(graph->get_sequence(first_path_handle)[graph->get_length(first_path_handle) - 1]);
        } else {
            // They are flipped: leave off the beginning of the start node going backward
            prev_char = reverse_complement(::toupper(graph->get_sequence(first_path_handle)[0]));
        }
            
        // shift from 0-based to 1-based for VCF
        first_path_pos += 1;

        v.position = first_path_pos + ref_trav_offset;

        v.id = print_snarl(graph, snarl_start, snarl_end);
            
        // Convert the snarl traversals to strings and add them to the variant
        vector<bool> use_trav(travs.size());
        if (path_trav_to_ref_trav.size()) {
            for (uint64_t i = 0; i < use_trav.size(); ++i) {
                use_trav[i] = (ref_trav_idx == path_trav_to_ref_trav[i]);
            }
        } else {
            for (uint64_t i = 0; i < use_trav.size(); ++i) {
                use_trav[i] = true;
            }
        }

        if (std::none_of(use_trav.begin(), use_trav.end(), [](bool b) {return b;})) {
            // no alts were jaccard-assigned to this reference, so abort before an assertion gets tripped
            continue;
        }

        // Sort the traversals for clustering
        vector<int> sorted_travs = get_traversal_order(graph, travs, trav_path_names, ref_travs, ref_trav_idx, use_trav);

        // jaccard clustering (using handles for now) on traversals
        vector<pair<double, int64_t>> trav_cluster_info;
        vector<int> child_snarl_to_trav;
        vector<vector<int>> trav_clusters = cluster_traversals(graph, travs, sorted_travs,
                                                               (in_nesting_info ? in_nesting_info->child_snarls :
                                                                vector<pair<handle_t, handle_t>>()),
                                                               cluster_threshold,
                                                               trav_cluster_info,
                                                               child_snarl_to_trav);

#ifdef debug
        cerr << "cluster priority";
        for (const auto& t: sorted_travs) {
            cerr << " " << t;
        }
        cerr << endl;
        for (const auto& tc : trav_clusters) {
            cerr << "traversal cluster: { ";
            for (const auto& t: tc) {
                cerr << t << "(" << trav_cluster_info[t].first << "," << trav_cluster_info[t].second << ") ";
            }
            cerr << " }" << endl;
        }
#endif

        unordered_map<string, vector<int>> sample_to_haps;                     
        if (in_nesting_info != nullptr) {
            // if the reference traversal is also an alt traversal, we pop out an extra copy
            // todo: this is a hack add add off-reference support while keeping the current
            // logic where the reference traversal is always distinct from the alts. this step
            // could be avoided, but it would come at the cost of some detailed refactoring of the
            // allele getting code...
            string ref_sample_name = PathMetadata::parse_sample_name(trav_path_names[ref_trav_idx]);
            if (this->sample_names.count(ref_sample_name)) {
                int alt_trav_copy = travs.size();
                travs.push_back(travs[ref_trav_idx]);
                trav_path_names.push_back(trav_path_names[ref_trav_idx]);
                trav_cluster_info.push_back(make_pair(0, 0));
                if (trav_steps.size() == travs.size()) {
                    trav_steps.push_back(trav_steps[ref_trav_idx]);
                }
                bool found_cluster = false;
                for (vector<int>& cluster : trav_clusters) {
                    if (cluster[0] == ref_trav_idx) {
                        found_cluster =true;
                        cluster.push_back(alt_trav_copy);
                        break;
                    }
                }
                assert(found_cluster == true);
            }
            
            // add in the star alleles -- these are alleles that were genotyped in the parent but not
            // the current allele, and are treated as *'s in VCF.
            if (this->star_allele) {
                sample_to_haps = add_star_traversals(travs, trav_path_names, trav_clusters, trav_cluster_info,
                                                     in_nesting_info->sample_to_haplotypes);
            }
            
        }

        vector<int> trav_to_allele = get_alleles(v, travs, trav_steps,
                                                 ref_trav_idx,
                                                 trav_clusters,
                                                 prev_char, use_start);

      
#ifdef debug
        assert(trav_to_allele.size() == travs.size());
        cerr << "trav_to_allele =";
        for (const auto& tta : trav_to_allele) {
            cerr << " " << tta;
        }
        cerr << endl;
#endif           

        // Fill in the genotypes
        get_genotypes(v, trav_path_names, trav_to_allele, trav_cluster_info);

        // Fill in some nesting-specific (site-level) tags
        NestingInfo ref_info; // since in_nesting_info is const, we put top-level stuff here       
        if (this->nested_decomposition) {
            if (in_nesting_info != nullptr && in_nesting_info->has_ref == true) {
                // if we're a child, just take what's passed in
                ref_info.parent_allele = in_nesting_info->parent_allele;
                ref_info.parent_len = in_nesting_info->parent_len;
                ref_info.parent_ref_len = in_nesting_info->parent_ref_len;
                ref_info.lv0_ref_name = in_nesting_info->lv0_ref_name;
                ref_info.lv0_ref_start = in_nesting_info->lv0_ref_start;
                ref_info.lv0_ref_len = in_nesting_info->lv0_ref_len;
                ref_info.lv0_alt_len = in_nesting_info->lv0_alt_len;
            } else {
                // if we're a root, compute values from the prsent site
                // todo: should they just be left undefined?
                ref_info.parent_allele = 0;
                ref_info.parent_len = v.alleles[0].length();
                ref_info.parent_ref_len = v.alleles[0].length();
                ref_info.lv0_ref_name = v.sequenceName;
                ref_info.lv0_ref_start = v.position;
                ref_info.lv0_ref_len = v.alleles[0].length();
                ref_info.lv0_alt_len = v.alleles[ref_info.parent_allele].length();
            }            
            v.info["PA"].push_back(std::to_string(ref_info.parent_allele));
            v.info["PL"].push_back(std::to_string(ref_info.parent_len));
            v.info["PR"].push_back(std::to_string(ref_info.parent_ref_len));            
            v.info["RC"].push_back(ref_info.lv0_ref_name);
            v.info["RS"].push_back(std::to_string(ref_info.lv0_ref_start));
            v.info["RD"].push_back(std::to_string(ref_info.lv0_ref_len));
            v.info["RL"].push_back(std::to_string(ref_info.lv0_alt_len));
        }

        if (i == 0 && out_nesting_infos != nullptr) {
            // we pass some information down to the children
            // todo: do/can we consider all the diferent reference intervals?
            //       right now, the info passed is hopefully coarse-grained enough not to matter?
            assert(in_nesting_info != nullptr &&
                   in_nesting_info->child_snarls.size() == out_nesting_infos->size());

            for (int64_t j = 0; j < out_nesting_infos->size(); ++j) {
                out_nesting_infos->at(j).child_snarls.clear();
                out_nesting_infos->at(j).has_ref = false;
                if (child_snarl_to_trav[j] >= 0) {
                    if (child_snarl_to_trav[j] < trav_steps.size()) {
                        NestingInfo& child_info = out_nesting_infos->at(j);
                        child_info.has_ref = true;
                        child_info.parent_path_interval = trav_steps[child_snarl_to_trav[j]];
                        child_info.sample_to_haplotypes = sample_to_haps;
                        child_info.parent_allele = trav_to_allele[child_snarl_to_trav[j]] >= 0 ?
                            trav_to_allele[child_snarl_to_trav[j]] : 0;
                        child_info.parent_len = v.alleles[child_info.parent_allele].length();
                        child_info.parent_ref_len = v.alleles[0].length();
                        child_info.lv0_ref_name = ref_info.lv0_ref_name;
                        child_info.lv0_ref_start = ref_info.lv0_ref_start;
                        child_info.lv0_ref_len = ref_info.lv0_ref_len;
                        if (in_nesting_info == nullptr || in_nesting_info->has_ref == false) {
                            // we're the parent of root, so we want to set this here
                            child_info.lv0_alt_len = child_info.parent_len;
                        } else {
                            child_info.lv0_alt_len = ref_info.lv0_alt_len;
                        }
                    }
                }
            }
        }

        // we only bother printing out sites with at least 1 non-reference allele
        if (!std::all_of(trav_to_allele.begin(), trav_to_allele.end(), [](int i) { return (i == 0 || i == -1); })) {
            // run vcffixup to add some basic INFO like AC
            vcf_fixup(v);
            bool added = add_variant(v);
            if (!added) {
                stringstream ss;
                ss << v;
                cerr << "Warning [vg deconstruct]: Skipping variant at " << v.sequenceName << ":" << v.position
                     << " with ID=" << v.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
                     << VCFOutputCaller::max_vcf_line_length << endl;
                return false;            
            }
        }
    }
    return true;
}

string Deconstructor::get_vcf_header() {
    // Keep track of the non-reference paths in the graph.  They'll be our sample names
    ref_samples.clear();
    set<size_t> ref_haplotypes;
    for (const string& ref_path_name : ref_paths) {
        ref_samples.insert(PathMetadata::parse_sample_name(ref_path_name));
        ref_haplotypes.insert(PathMetadata::parse_haplotype(ref_path_name));
    }
    if (!long_ref_contig) {
        long_ref_contig = ref_samples.size() > 1 || ref_haplotypes.size() > 1 || nested_decomposition;
    }
    this->long_ref_contig = long_ref_contig;
    sample_names.clear();
    unordered_map<string, set<int>> sample_to_haps;

    // find sample names from non-reference paths
    graph->for_each_path_handle([&](const path_handle_t& path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (!this->ref_paths.count(path_name)) {
            string sample_name = graph->get_sample_name(path_handle);
            // for backward compatibility
            if (sample_name == PathMetadata::NO_SAMPLE_NAME) {
                sample_name = path_name;
            }
            if (!ref_samples.count(sample_name)) {
                size_t haplotype = graph->get_haplotype(path_handle);
                if (haplotype == PathMetadata::NO_HAPLOTYPE) {
                    haplotype = 0;
                }
                if (haplotype > 10) {
                    cerr << "Warning [vg deconstruct]: Suspiciously large haplotype, " << haplotype
                         << ", parsed from path, " << path_name << ": This will leed to giant GT entries.";
                }
                sample_to_haps[sample_name].insert((int)haplotype);
                sample_names.insert(sample_name);
            }
        }
    });
    
    // add in the GBWT sample names
    if (gbwt) {
        // add in sample names from the gbwt
        for (size_t i = 0; i < gbwt->metadata.paths(); i++) {
            PathSense sense = gbwtgraph::get_path_sense(*gbwt, i, gbwt_reference_samples);
            if (sense == PathSense::HAPLOTYPE) {
                string path_name = PathMetadata::create_path_name(
                    sense,
                    gbwtgraph::get_path_sample_name(*gbwt, i, sense),
                    gbwtgraph::get_path_locus_name(*gbwt, i, sense),
                    gbwtgraph::get_path_haplotype(*gbwt, i, sense),
                    gbwtgraph::get_path_phase_block(*gbwt, i, sense),
                    gbwtgraph::get_path_subrange(*gbwt, i, sense));
                if (!this->ref_paths.count(path_name)) {
                    string sample_name = gbwtgraph::get_path_sample_name(*gbwt, i, sense);
                    if (!ref_samples.count(sample_name)) {
                        auto phase = gbwtgraph::get_path_haplotype(*gbwt, i, sense);
                        if (phase == PathMetadata::NO_HAPLOTYPE) {
                            // Default to 0.
                            phase = 0;
                        }
                        if (phase > 10) {
                            cerr << "Warning [vg deconstruct]: Suspiciously large haplotype, " << phase
                                 << ", parsed from GBWT thread, " << path_name << ": This will leed to giant GT entries.";
                        }
                        sample_to_haps[sample_name].insert((int)phase);
                        sample_names.insert(sample_name);
                    }
                }
            }
        }
    }

    if (sample_to_haps.empty()) {
        cerr << "Error [vg deconstruct]: No paths other than selected reference(s) found in the graph, "
             << "so no alt alleles can be generated. Note that exhaustive path-free traversal finding "
             << "is no longer supported, and vg deconstruct now only works on embedded paths and GBWT "
             << "threads." << endl;
        exit(1);
    }

    // find some stats about the haplotypes for each sample    
    gbwt_sample_to_phase_range.clear();
    sample_ploidys.clear();
    for (auto& sample_haps : sample_to_haps) {
        sample_ploidys[sample_haps.first] = sample_haps.second.size();
        gbwt_sample_to_phase_range[sample_haps.first] = make_pair(*sample_haps.second.begin(), *sample_haps.second.rbegin());
    }
        
    // print the VCF header
    stringstream stream;
    stream << "##fileformat=VCFv4.2" << endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    if (show_path_info && path_to_sample_phase) {
        stream << "##FORMAT=<ID=PI,Number=.,Type=String,Description=\"Path information. Original vg path name for sample as well as its allele (can be many paths per sample)\">" << endl;
    }
    if (path_to_sample_phase || gbwt) {
        stream << "##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in the graph with conflicting alleles";
        if (!gbwt && show_path_info) {
            stream << " (details in PI field)";
        }
        stream << "\">" << endl;
    }
    if (path_to_sample_phase && cluster_threshold < 1) {
        stream << "##FORMAT=<ID=TS,Number=1,Type=Float,Description=\"Similarity between the sample's actual path and its allele\">"
               << endl;
        stream << "##FORMAT=<ID=TL,Number=1,Type=Integer,Description=\"Length difference between the sample's actual path and its allele\">"
               << endl;

    }
    stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl;
    stream << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl;
    stream << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl;
    stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl;
    if (include_nested) {
        stream << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << endl;
        stream << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    }
    if (this->nested_decomposition) {
        stream << "##INFO=<ID=PA,Number=1,Type=Integer,Description=\"Allele number of the reference allele in Parent Snarl\">" << endl;
        stream << "##INFO=<ID=PL,Number=1,Type=Integer,Description=\"Length of the reference allele in Parent Snarl\">" << endl;
        stream << "##INFO=<ID=PR,Number=1,Type=Integer,Description=\"Length of 0th allele in the Parent Snarl\">" << endl;
        stream << "##INFO=<ID=RC,Number=1,Type=String,Description=\"Reference chromosome name of top-level containing site\">" << endl;
        stream << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Reference start position of top-level containing site\">" << endl;
        stream << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"Reference end position name of top-level containing site\">" << endl;
        stream << "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Length of the top-level allele in which this site nests\">" << endl;
    }
    if (untangle_allele_traversals) {
        stream << "##INFO=<ID=UT,Number=R,Type=String,Description=\"Untangled allele Traversal with reference node start and end positions, format: [>|<][id]_[start|.]_[end|.], with '.' indicating non-reference nodes.\">" << endl;
    } else {
        stream << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    }
    
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto& sample_name : sample_names) {
        stream << "\t" << sample_name;
    }
    stream << endl;
    return stream.str();
}

string Deconstructor::add_contigs_to_vcf_header(const string& vcf_header) const {

    vector<string> header_lines = split_delims(vcf_header, "\n");

    stringstream patched_header;
    for (int64_t i = 0; i < header_lines.size() - 1; ++i) {
        patched_header << header_lines[i] << "\n";
    }

    set<string> all_ref_paths = this->ref_paths;
    
    // add in the off-ref paths that nested deconstruction may have found
    for (const unordered_set<path_handle_t>& off_ref_path_set : this->off_ref_paths) {
        for (const path_handle_t& off_ref_path : off_ref_path_set) {
            all_ref_paths.insert(graph->get_path_name(off_ref_path));
        }
    }
    
    map<string, int64_t> ref_path_to_length;
    for(auto& refpath : all_ref_paths) {
        assert(graph->has_path(refpath));
        int64_t path_len = 0;
        path_handle_t path_handle = graph->get_path_handle(refpath);
        for (handle_t handle : graph->scan_path(path_handle)) {
            path_len += graph->get_length(handle);
        }
        string locus_name = graph->get_locus_name(path_handle);
        if (locus_name == PathMetadata::NO_LOCUS_NAME) {
            locus_name = refpath;
        } else if (long_ref_contig) {
            // the sample name isn't unique enough, so put a full ugly name in the vcf
            if (graph->get_sense(path_handle) == PathSense::GENERIC) {
                locus_name = graph->get_path_name(path_handle);
            } else {
                locus_name = PathMetadata::create_path_name(PathSense::REFERENCE,
                                                            graph->get_sample_name(path_handle),
                                                            locus_name,
                                                            graph->get_haplotype(path_handle),
                                                            PathMetadata::NO_PHASE_BLOCK,
                                                            PathMetadata::NO_SUBRANGE);
            }
        }            

        subrange_t subrange = graph->get_subrange(path_handle);
        int64_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
        ref_path_to_length[locus_name] = std::max(ref_path_to_length[locus_name], path_len + offset);
    }
    for (auto& ref_path_len : ref_path_to_length) {
        patched_header << "##contig=<ID=" << ref_path_len.first << ",length=" << ref_path_len.second << ">" << endl;
    }

    assert(header_lines.back().substr(0, 6) == "#CHROM");
    patched_header << header_lines.back();
    return patched_header.str();
}

void Deconstructor::deconstruct_graph(SnarlManager* snarl_manager) {

    vector<const Snarl*> snarls;
    vector<const Snarl*> queue;

    // read all our snarls into a list
    snarl_manager->for_each_top_level_snarl([&](const Snarl* snarl) {
        queue.push_back(snarl);
    });
    if (include_nested) {
        while (!queue.empty()) {
            const Snarl* snarl = queue.back();
            queue.pop_back();
            snarls.push_back(snarl);
            const vector<const Snarl*>& children = snarl_manager->children_of(snarl);
            queue.insert(queue.end(), children.begin(), children.end());
        }
    } else {
        swap(snarls, queue);
    }

    // process the whole shebang in parallel
#pragma omp parallel for schedule(dynamic,1)
    for (size_t i = 0; i < snarls.size(); i++) {
        deconstruct_site(graph->get_handle(snarls[i]->start().node_id(), snarls[i]->start().backward()),
                         graph->get_handle(snarls[i]->end().node_id(), snarls[i]->end().backward()));
    }
}

void Deconstructor::deconstruct_graph_top_down(SnarlManager* snarl_manager) {
    // logic copied from vg call (graph_caller.cpp)
    
    size_t thread_count = get_thread_count();
    this->off_ref_paths.clear();
    this->off_ref_paths.resize(get_thread_count());
    // Used to recurse on children of parents that can't be called
    vector<vector<pair<const Snarl*, NestingInfo>>> snarl_queue(thread_count);

    // Run the deconstructor on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl, NestingInfo nesting_info) {
        if (!snarl_manager->is_trivial(snarl, *graph)) {
            const vector<const Snarl*>& children = snarl_manager->children_of(snarl);
            assert(nesting_info.child_snarls.empty());
            for (const Snarl* child : children) {
                nesting_info.child_snarls.push_back(make_pair(graph->get_handle(child->start().node_id(), child->start().backward()),
                                                              graph->get_handle(child->end().node_id(), child->end().backward())));
                                                    
            }
            vector<NestingInfo> out_nesting_infos(children.size());
            bool was_deconstructed = deconstruct_site(graph->get_handle(snarl->start().node_id(), snarl->start().backward()),
                                                      graph->get_handle(snarl->end().node_id(), snarl->end().backward()),
                                                      include_nested ? &nesting_info : nullptr,
                                                      include_nested ? &out_nesting_infos : nullptr);
            if (include_nested || !was_deconstructed) {
                vector<pair<const Snarl*, NestingInfo>>& thread_queue = snarl_queue[omp_get_thread_num()];
                for (int64_t i = 0; i < children.size(); ++i) {
                    thread_queue.push_back(make_pair(children[i], out_nesting_infos[i]));
                }

            }            
        }
    };

    // Start with the top level snarls
    // (note, can't do for_each_top_level_snarl_parallel() because interface wont take nesting info)
    vector<pair<const Snarl*, NestingInfo>> top_level_snarls;
    snarl_manager->for_each_top_level_snarl([&](const Snarl* snarl) {
        NestingInfo nesting_info;
        nesting_info.has_ref = false;
        top_level_snarls.push_back(make_pair(snarl, nesting_info));
    });
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < top_level_snarls.size(); ++i) {
        process_snarl(top_level_snarls[i].first, top_level_snarls[i].second);
    }

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(snarl_queue.begin(), snarl_queue.end(),
                        [](const vector<pair<const Snarl*, NestingInfo>>& snarl_vec) {return snarl_vec.empty();})) {
        vector<pair<const Snarl*, NestingInfo>> cur_queue;
        for (vector<pair<const Snarl*, NestingInfo>>& thread_queue : snarl_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i].first, cur_queue[i].second);
        }
    }
}

/**
 * Convenience wrapper function for deconstruction of multiple paths.
 */
void Deconstructor::deconstruct(vector<string> ref_paths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                                bool include_nested,
                                int context_jaccard_window,
                                bool untangle_traversals,
                                bool keep_conflicted,
                                bool strict_conflicts,
                                bool long_ref_contig,
                                double cluster_threshold,
                                gbwt::GBWT* gbwt,
                                bool nested_decomposition,
                                bool star_allele) {

    this->graph = graph;
    this->ref_paths = set<string>(ref_paths.begin(), ref_paths.end());
    this->include_nested = include_nested || nested_decomposition;
    this->path_jaccard_window = context_jaccard_window;
    this->untangle_allele_traversals = untangle_traversals;
    this->keep_conflicted_genotypes = keep_conflicted;
    this->strict_conflict_checking = strict_conflicts;
    this->long_ref_contig = long_ref_contig;
    if (gbwt) {
        this->gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(*gbwt);
    }
    this->cluster_threshold = cluster_threshold;
    this->gbwt = gbwt;
    this->nested_decomposition = nested_decomposition;
    this->star_allele = star_allele;

    // the need to use nesting is due to a problem with omp tasks and shared state
    // which results in extremely high memory costs (ex. ~10x RAM for 2 threads vs. 1)
    omp_set_nested(1);
    omp_set_max_active_levels(3);
    
    // create the traversal finder
    map<string, const Alignment*> reads_by_name;
    path_trav_finder = unique_ptr<PathTraversalFinder>(new PathTraversalFinder(*graph));
        
    if (gbwt != nullptr) {
        gbwt_trav_finder = unique_ptr<GBWTTraversalFinder>(new GBWTTraversalFinder(*graph, *gbwt));
    }

    string hstr = this->get_vcf_header();
    assert(output_vcf.openForOutput(hstr));

    if (nested_decomposition) {
        deconstruct_graph_top_down(snarl_manager);
    } else {
        deconstruct_graph(snarl_manager);
    }

    string patched_header = this->add_contigs_to_vcf_header(output_vcf.header);
    cout << patched_header << endl;

    // write variants in sorted order
    write_variants(cout, snarl_manager);
}


}

