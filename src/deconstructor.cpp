#include "deconstructor.hpp"
#include "traversal_finder.hpp"
#include <gbwtgraph/gbwtgraph.h>
#include "rgfa.hpp"

//#define debug

using namespace std;


namespace vg {
Deconstructor::Deconstructor() : VCFOutputCaller(""),
                                 exhaustive_jaccard_warning(false){
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
                                       const pair<vector<SnarlTraversal>,
                                       vector<pair<step_handle_t, step_handle_t>>>& path_travs,
                                       int ref_path_idx,
                                       const vector<bool>& use_trav,
                                       char prev_char, bool use_start) const {

    auto& travs = path_travs.first;
    assert(ref_path_idx >=0 && ref_path_idx < travs.size());

    // map strings to allele numbers (and their traversal)
    // (we are using the traversal finder in such a way that duplicate alleles can get returned
    // in order to be able to preserve the path names)
    map<string, pair<int, int>> allele_idx;
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
    allele_idx[ref_allele] = make_pair(0, ref_path_idx);
    trav_to_allele[ref_path_idx] = 0;
    bool substitution = true;
        
    // set the other alleles (they can end up as 0 alleles too if their strings match the reference)
    for (int i = 0; i < travs.size(); ++i) {
        if (i != ref_path_idx) {
            if (use_trav[i]) {
                string allele = trav_to_string(travs[i]);
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
        path_handle_t ref_path = graph->get_path_handle_of_step(path_travs.second.at(ref_path_idx).first);
        unordered_map<nid_t, vector<pair<uint64_t, dac_vector<>>>> ref_dup_nodes;
        unordered_map<nid_t, nid_t> ref_simple_pos;
        {
            auto& trav = travs.at(ref_path_idx);
            for (size_t i = 0; i < trav.visit_size(); ++i) {
                size_t j = !reversed ? i : trav.visit_size() - 1 - i;
                const Visit& visit = trav.visit(j);
                nid_t node_id = visit.node_id();
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
            auto start_step = path_travs.second.at(allele_trav_no).first;
            auto end_step = path_travs.second.at(allele_trav_no).second;
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
                                double j = context_jaccard(ref_context, path_context);
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
            add_allele_path_to_info(v, allele_no, travs[allele_trav_no], reversed, !substitution);
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
                                  const vector<int>& trav_to_allele) const {
    assert(names.size() == trav_to_allele.size());
    // set up our variant fields
    v.format.push_back("GT");
    if (show_path_info && path_to_sample_phase && path_restricted) {
        v.format.push_back("PI");
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
                auto& phase_range = gbwt_sample_to_phase_range.at(sample_name);
                for (int phase = phase_range.first + 1; phase <= phase_range.second; ++phase) {
                    blank_gt += "|.";
                }
            }
            v.samples[sample_name]["GT"] = {blank_gt};
            if (show_path_info && path_to_sample_phase && path_restricted) {
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
    int sample_ploidy = ploidy;
    int min_phase = 1;
    int max_phase = ploidy;
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


// todo refactor if we need to reuse elsewhere in vg
// implemented inline for development
// assumes sorted input
double Deconstructor::context_jaccard(
    const vector<nid_t>& target,
    const vector<nid_t>& query) const {
    size_t node_isec = 0;
    std::set_intersection(target.begin(), target.end(),
                          query.begin(), query.end(),
                          count_back_inserter<nid_t>(node_isec));
    size_t node_union = 0;
    std::set_union(target.begin(), target.end(),
                   query.begin(), query.end(),
                   count_back_inserter<nid_t>(node_union));
    return (double)node_isec / (double)node_union;
}

double Deconstructor::context_jaccard(
    const dac_vector<>& target,
    const vector<nid_t>& query) const {
    size_t node_isec = 0;
    std::set_intersection(target.begin(), target.end(),
                          query.begin(), query.end(),
                          count_back_inserter<nid_t>(node_isec));
    size_t node_union = 0;
    std::set_union(target.begin(), target.end(),
                   query.begin(), query.end(),
                   count_back_inserter<nid_t>(node_union));
    return (double)node_isec / (double)node_union;
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

vector<nid_t> Deconstructor::get_context(
    const pair<vector<SnarlTraversal>,
               vector<pair<step_handle_t, step_handle_t>>>& path_travs,
    const int& trav_idx) const {
    step_handle_t start_step = path_travs.second[trav_idx].first;
    step_handle_t end_step = path_travs.second[trav_idx].second;
    return get_context(start_step, end_step);
}

bool Deconstructor::deconstruct_site(const Snarl* snarl) const {

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
    for (const pair<step_handle_t, step_handle_t>& trav_ends : path_travs.second) {
        path_trav_names.push_back(graph->get_path_name(graph->get_path_handle_of_step(trav_ends.first)));
    }

    // pick out the traversal corresponding to a reference path, breaking ties consistently
    string ref_trav_name;
    for (int i = 0; i < path_travs.first.size(); ++i) {
        const string& path_trav_name = path_trav_names.at(i);
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
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Setting ref_trav_name " << ref_trav_name << endl;
#endif
        }
    }
    
    // add in the gbwt traversals
    // after this, all traversals are treated the same, with metadata embedded in their names
    int64_t first_gbwt_trav_idx = path_trav_names.size();
    vector<gbwt::size_type> gbwt_path_ids;
    if (gbwt_trav_finder.get() != nullptr) {
        const gbwt::GBWT& gbwt_index = gbwt_trav_finder->get_gbwt();
        pair<vector<SnarlTraversal>, vector<gbwt::size_type>> thread_travs = gbwt_trav_finder->find_path_traversals(*snarl);
        for (int i = 0; i < thread_travs.first.size(); ++i) {
            // We need to get a bunch of metadata about the path, but the GBWT
            // we have might not even have structured path names stored.
            gbwt::size_type path_id = gbwt::Path::id(thread_travs.second[i]);
            if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || path_id >= gbwt_index.metadata.paths()) {
                continue;
            }
            
            gbwt_path_ids.push_back(path_id);
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
                path_trav_names.push_back(path_name);
                path_travs.first.push_back(thread_travs.first[i]);
                // dummy handles so we can use the same code as the named path traversals above
                path_travs.second.push_back(make_pair(step_handle_t(), step_handle_t()));
            }
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
        for (int i = 0; i < path_travs.first.size(); ++i) {
            const string& path_trav_name = path_trav_names.at(i);
            subrange_t subrange ;
            Paths::strip_subrange(path_trav_name, &subrange);
            int64_t sub_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
            if (path_trav_name == ref_trav_name) {
                ref_travs.push_back(i);
                ref_offsets.push_back(sub_offset);
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Adding ref_tav idx=" << i << " offset=" << sub_offset << " because " << path_trav_name << " == " << ref_trav_name << endl;
#endif                
            }
        }
    }

    // there's no reference path through the snarl, so we can't make a variant
    // (todo: should we try to detect this before computing traversals?)
    if (ref_travs.empty()) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site because no reference traversal was found " << pb2json(*snarl) << endl;
#endif
        return false;
    }    

    bool exhaustive = !path_restricted && gbwt_trav_finder.get() == nullptr;
    if (exhaustive) {        
        // add in the exhaustive traversals
        vector<SnarlTraversal> additional_travs;
                        
        // exhaustive traversal can't do all snarls
        if (snarl->type() != ULTRABUBBLE) {
            return false;
        }
        if (!check_max_nodes(snarl)) {
#pragma omp critical (cerr)
            cerr << "Warning: Skipping site because it is too complex for exhaustive traversal enumeration: " << pb2json(*snarl) << endl << "         Consider using -e to traverse embedded paths" << endl;
            return false;
        }
        additional_travs = explicit_exhaustive_traversals(snarl);
         
        // happens when there was a nested non-ultrabubble snarl
        if (additional_travs.empty()) {
            return false;
        }
        path_travs.first.insert(path_travs.first.end(), additional_travs.begin(), additional_travs.end());
        for (int i = 0; i < additional_travs.size(); ++i) {
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

    // XXX CHECKME this assumes there is only one reference path here, and that multiple traversals are due to cycles
    
    // we collect windows around the reference traversals
    // to compare with equivalent windows from the alternate allele paths
    // we will associate these 1:1 with reference traversals

    // remember that path_travs := pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs;

    // map from each path_trav index to the ref_trav index it best maps to
    vector<int> path_trav_to_ref_trav;

    if (ref_travs.size() > 1 && this->path_jaccard_window && exhaustive && !exhaustive_jaccard_warning) {
#pragma omp critical (cerr)
        cerr << "warning [vg deconstruct]: Conext Jaccard logic for multiple references disabled with exhaustive traversals. Use -e, -g or GBZ input to switch to path-based traversals only (recommended)." << endl;
        exhaustive_jaccard_warning = true;
    }
    if (ref_travs.size() > 1 && this->path_jaccard_window && !exhaustive) {
        path_trav_to_ref_trav.resize(path_travs.first.size());
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Multiple ref traversals!" << endl;
#endif
        {
            vector<vector<nid_t>> ref_contexts(ref_travs.size());
#pragma omp parallel for schedule(dynamic,1)
            for (size_t i = 0; i < ref_travs.size(); ++i) {
                auto& trav_id = ref_travs[i];
                ref_contexts[i] = get_context(path_travs, trav_id);
            }
            // now for each traversal, we compute and equivalent context and match it to a ref context
            // using a jaccard metric over node ids
#pragma omp parallel for schedule(dynamic,1)
            for (size_t i = 0; i < path_travs.first.size(); ++i) {
                vector<nid_t> context = get_context(path_travs, i);
                // map jaccard metric to the index of the ref_trav
                vector<pair<double, int>> ref_mappings;
                for (uint64_t j = 0; j < ref_travs.size(); ++j) {
                    ref_mappings.push_back(make_pair(
                                               context_jaccard(
                                                   ref_contexts[j],
                                                   context),
                                               ref_travs[j]));
                }
                std::sort(ref_mappings.begin(), ref_mappings.end());
                // the best is the last, which has the highest jaccard
                path_trav_to_ref_trav[i] = ref_mappings.back().second;
            }
        }
    }

    // we write a variant for every reference traversal
    // (optionally) selecting the subset of path traversals that are 1:1
//#pragma omp parallel for
    for (size_t i = 0; i < ref_travs.size(); ++i) {
//#pragma omp task firstprivate(i)
        {
            auto& ref_trav_idx = ref_travs[i];
            auto& ref_trav_offset = ref_offsets[i];

            const SnarlTraversal& ref_trav = path_travs.first[ref_trav_idx];

            string fixed_trav_name = RGFACover::revert_rgfa_path_name(ref_trav_name, true);

            vcflib::Variant v;
            v.quality = 60;

            // in VCF we usually just want the contig
            string contig_name = PathMetadata::parse_locus_name(fixed_trav_name);
            if (contig_name == PathMetadata::NO_LOCUS_NAME) {
                contig_name = fixed_trav_name;
            } else if (long_ref_contig) {
                // the sample name isn't unique enough, so put a full ugly name in the vcf
                if (PathMetadata::parse_sense(fixed_trav_name) == PathSense::GENERIC) {
                    contig_name = fixed_trav_name;
                } else {
                    contig_name = PathMetadata::create_path_name(PathSense::REFERENCE,
                                                                 PathMetadata::parse_sample_name(fixed_trav_name),
                                                                 contig_name,
                                                                 PathMetadata::parse_haplotype(fixed_trav_name),
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
            assert(ref_trav_idx < first_gbwt_trav_idx);
            step_handle_t start_step = path_travs.second[ref_trav_idx].first;
            step_handle_t end_step = path_travs.second[ref_trav_idx].second;
            handle_t start_handle = graph->get_handle_of_step(start_step);
            handle_t end_handle = graph->get_handle_of_step(end_step);
            size_t start_pos = graph->get_position_of_step(start_step);
            size_t end_pos = graph->get_position_of_step(end_step);
            use_start = start_pos < end_pos;
            first_path_handle = use_start ? start_handle : end_handle;
            first_path_pos = use_start ? start_pos : end_pos;
            
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

            v.position = first_path_pos + ref_trav_offset;

            v.id = print_snarl(*snarl);
            
            // Convert the snarl traversals to strings and add them to the variant
            vector<bool> use_trav(path_travs.first.size());
            if (path_trav_to_ref_trav.size()) {
                for (uint64_t i = 0; i < use_trav.size(); ++i) {
                    use_trav[i] = (ref_trav_idx == path_trav_to_ref_trav[i]);
                }
            } else {
                for (uint64_t i = 0; i < use_trav.size(); ++i) {
                    use_trav[i] = true;
                }
            }

            vector<int> trav_to_allele = get_alleles(v, path_travs,
                                                     ref_trav_idx,
                                                     use_trav,
                                                     prev_char, use_start);

            // Fill in the genotypes
            if (path_restricted || gbwt_trav_finder.get()) {
                get_genotypes(v, path_trav_names, trav_to_allele);
            }

            // we only bother printing out sites with at least 1 non-reference allele
            if (!std::all_of(trav_to_allele.begin(), trav_to_allele.end(), [](int i) { return (i == 0 || i == -1); })) {
                if (path_restricted || gbwt_trav_finder.get()) {
                    // run vcffixup to add some basic INFO like AC
                    vcf_fixup(v);
                }
                add_variant(v);
            }
        }
    }
//#pragma omp taskwait
    return true;
}

/**
 * Convenience wrapper function for deconstruction of multiple paths.
 */
void Deconstructor::deconstruct(vector<string> ref_paths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                                bool path_restricted_traversals,
                                int ploidy,
                                bool include_nested,
                                int context_jaccard_window,
                                bool untangle_traversals,
                                bool keep_conflicted,
                                bool strict_conflicts,
                                bool long_ref_contig,
                                gbwt::GBWT* gbwt) {

    this->graph = graph;
    this->snarl_manager = snarl_manager;
    this->path_restricted = path_restricted_traversals;
    this->ploidy = ploidy;
    this->ref_paths = set<string>(ref_paths.begin(), ref_paths.end());
    this->include_nested = include_nested;
    this->path_jaccard_window = context_jaccard_window;
    this->untangle_allele_traversals = untangle_traversals;
    this->keep_conflicted_genotypes = keep_conflicted;
    this->strict_conflict_checking = strict_conflicts;
    if (gbwt) {
        this->gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(*gbwt);
    }

    // the need to use nesting is due to a problem with omp tasks and shared state
    // which results in extremely high memory costs (ex. ~10x RAM for 2 threads vs. 1)
    omp_set_nested(1);
    omp_set_max_active_levels(3);

    // Keep track of the non-reference paths in the graph.  They'll be our sample names
    ref_samples.clear();
    set<size_t> ref_haplotypes;
    for (const string& ref_path_name : ref_paths) {
        ref_samples.insert(PathMetadata::parse_sample_name(ref_path_name));
        ref_haplotypes.insert(PathMetadata::parse_haplotype(ref_path_name));
    }
    if (!long_ref_contig) {
        long_ref_contig = ref_samples.size() > 1 || ref_haplotypes.size() > 1;
    }
    this->long_ref_contig = long_ref_contig;
    sample_names.clear();
    unordered_map<string, set<int>> sample_to_haps;

    // find sample names from non-reference paths
    graph->for_each_path_handle([&](const path_handle_t& path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (!this->ref_paths.count(path_name)) {
            path_name = RGFACover::revert_rgfa_path_name(path_name, true);
            string sample_name = PathMetadata::parse_sample_name(path_name);
            // for backward compatibility
            if (sample_name == PathMetadata::NO_SAMPLE_NAME) {
                sample_name = path_name;
            }
            if (!ref_samples.count(sample_name)) {
                size_t haplotype = PathMetadata::parse_haplotype(path_name);
                if (haplotype == PathMetadata::NO_HAPLOTYPE) {
                    haplotype = 0;
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
                        sample_to_haps[sample_name].insert((int)phase);
                        sample_names.insert(sample_name);
                    }
                }
            }
        }
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
    if (path_restricted || gbwt) {
        stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    }
    if (show_path_info && path_to_sample_phase && path_restricted) {
        stream << "##FORMAT=<ID=PI,Number=.,Type=String,Description=\"Path information. Original vg path name for sample as well as its allele (can be many paths per sample)\">" << endl;
    }
    if (path_to_sample_phase || gbwt) {
        stream << "##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in the graph with conflicting alleles";
        if (!gbwt && show_path_info) {
            stream << " (details in PI field)";
        }
        stream << "\">" << endl;
    }
    if (path_restricted || gbwt) {
        stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl;
        stream << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl;
        stream << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl;
        stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl;
    }
    if (include_nested) {
        stream << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << endl;
        stream << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    }
    if (untangle_allele_traversals) {
        stream << "##INFO=<ID=UT,Number=R,Type=String,Description=\"Untangled allele Traversal with reference node start and end positions, format: [>|<][id]_[start|.]_[end|.], with '.' indicating non-reference nodes.\">" << endl;
    } else {
        stream << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    }
    set<string> gbwt_ref_paths;
    map<string, int64_t> ref_path_to_length;
    for(auto& refpath : ref_paths) {
        if (graph->has_path(refpath)) {
            int64_t path_len = 0;
            path_handle_t path_handle = graph->get_path_handle(refpath);
            for (handle_t handle : graph->scan_path(path_handle)) {
                path_len += graph->get_length(handle);
            }
            string fixed_refpath = RGFACover::revert_rgfa_path_name(refpath, true);
            string locus_name = PathMetadata::parse_locus_name(fixed_refpath);
            if (locus_name == PathMetadata::NO_LOCUS_NAME) {
                locus_name = fixed_refpath;
            } else if (long_ref_contig) {
                // the sample name isn't unique enough, so put a full ugly name in the vcf
                if (PathMetadata::parse_sense(fixed_refpath) == PathSense::GENERIC) {
                    locus_name = fixed_refpath;
                } else {
                    locus_name = PathMetadata::create_path_name(PathSense::REFERENCE,
                                                                PathMetadata::parse_sample_name(fixed_refpath),
                                                                locus_name,
                                                                PathMetadata::parse_haplotype(fixed_refpath),
                                                                PathMetadata::NO_PHASE_BLOCK,
                                                                PathMetadata::NO_SUBRANGE);
                }
            }            

            subrange_t subrange = graph->get_subrange(path_handle);
            int64_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
            ref_path_to_length[locus_name] = std::max(ref_path_to_length[locus_name], path_len + offset);
        } else {
            gbwt_ref_paths.insert(refpath);
        }       
    }
    for (auto& ref_path_len : ref_path_to_length) {
        stream << "##contig=<ID=" << ref_path_len.first << ",length=" << ref_path_len.second << ">" << endl;
    }
    if (!gbwt_ref_paths.empty()) {
        unordered_map<string, vector<gbwt::size_type>> gbwt_name_to_ids;
        for (size_t i = 0; i < gbwt->metadata.paths(); i++) {
            // Collect all the GBWT path IDs for each sample and contig.
            gbwt_name_to_ids[compose_short_path_name(*gbwt, i)].push_back(i);
        }
        for (const string& refpath : gbwt_ref_paths) {
            // For each sample and contig name that is a GBWT ref path
            vector<gbwt::size_type>& thread_ids = gbwt_name_to_ids.at(refpath);
            size_t path_len = 0;
            for (gbwt::size_type thread_id : thread_ids) {
                // For each actual path in the GBWT for that sample-and-contig,
                // we need to see how long it extends the space of the sample
                // and contig.
                
                // TODO: These are probably all guaranteed to be haplotype sense?
                PathSense sense = gbwtgraph::get_path_sense(*gbwt, thread_id, gbwt_reference_samples);
                subrange_t subrange = gbwtgraph::get_path_subrange(*gbwt, thread_id, sense);
                
                // TODO: when importing GFAs we might cram the start of a walk
                // into the GBWT count field. But we don't ever guarantee that
                // we've done that so it might not be visible as a subrange
                // here. Fix that somehow???
                size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                size_t len = path_to_length(extract_gbwt_path(*graph, *gbwt, thread_id));
                path_len = std::max(path_len, offset + len);
            }
            stream << "##contig=<ID=" << refpath << ",length=" << path_len << ">" << endl;
        }
    }
    
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (path_restricted || gbwt) {
        for (auto& sample_name : sample_names) {
            stream << "\t" << sample_name;
        }
    }
    stream << endl;
    
    string hstr = stream.str();
    assert(output_vcf.openForOutput(hstr));
    cout << output_vcf.header << endl;

    // create the traversal finder
    map<string, const Alignment*> reads_by_name;
    path_trav_finder = unique_ptr<PathTraversalFinder>(new PathTraversalFinder(*graph,
                                                                               *snarl_manager));
    
    if (!path_restricted && !gbwt) {
        trav_finder = unique_ptr<TraversalFinder>(new ExhaustiveTraversalFinder(*graph,
                                                                                *snarl_manager,
                                                                                true));

    }
    
    if (gbwt != nullptr) {
        gbwt_trav_finder = unique_ptr<GBWTTraversalFinder>(new GBWTTraversalFinder(*graph, *gbwt));
    }

    vector<const Snarl*> snarls_todo;
    // Do the top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl([&](const Snarl* snarl) {
            vector<const Snarl*> todo(1, snarl);
            vector<const Snarl*> next;
            while (!todo.empty()) {
                for (auto next_snarl : todo) {
                    // if we can't make a variant from the snarl due to not finding
                    // paths through it, we try again on the children
                    // note: we may want to push the parallelism down a bit
#pragma omp critical (snarls_todo)
                    snarls_todo.push_back(next_snarl);
                    if (include_nested) {
                        // n.b. we no longer attempt to deconstruct the site to determine if we nest
                        const vector<const Snarl*>& children = snarl_manager->children_of(next_snarl);
                        next.insert(next.end(), children.begin(), children.end());
                    }
                }
                swap(todo, next);
                next.clear();
            }
        });

//#pragma omp parallel
//#pragma omp single
    {
#pragma omp parallel for schedule(dynamic,1)
        for (size_t i = 0; i < snarls_todo.size(); i++) {
//#pragma omp task firstprivate(i)
            {
                auto& snarl = snarls_todo[i];
                deconstruct_site(snarl);
            }
        }
    }
//#pragma omp taskwait

    // write variants in sorted order
    write_variants(cout, snarl_manager);
}

bool Deconstructor::check_max_nodes(const Snarl* snarl) const  {
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

vector<SnarlTraversal> Deconstructor::explicit_exhaustive_traversals(const Snarl* snarl) const {
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

