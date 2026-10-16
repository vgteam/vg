#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "vcf_record.hpp"
#include "path.hpp"
#include "traversal_clusters.hpp"

//#define debug

namespace vg {

bool visits_to_walk(const HandleGraph& graph, const vector<NodeVisit>& visits, Traversal& walk) {
    // cluster_traversals asserts size() >= 2, and a visit to a child site has no single handle.
    // Both are real inputs here (the "*" placeholder, and alleles that step over child sites), so
    // refuse rather than fabricate something.
    if (visits.size() < 2) {
        return false;
    }
    walk.clear();
    walk.reserve(visits.size());
    for (const NodeVisit& visit : visits) {
        if (visit.first <= 0) {
            return false;
        }
        walk.push_back(graph.get_handle(visit.first, visit.second));
    }
    return true;
}

namespace {
/// Parse a VCF FORMAT value without throwing and without exiting.  vg::parse<double> exits on
/// failure and the 2-argument vg::parse can throw; merge_similar_alleles runs inside an OpenMP
/// region, where an escaping exception is std::terminate rather than something a caller can handle,
/// and exit() would abandon whatever the other threads had already buffered.  Missing values
/// ("." and "") are ordinary input here.
bool parse_vcf_double(const string& field, double& value) {
    try {
        size_t after;
        value = std::stod(field, &after);
        return after == field.size();
    } catch (const std::exception&) {
        return false;
    }
}
}

int64_t allele_core_length(const vector<string>& alleles) {
    vector<const string*> seqs;
    for (const string& a : alleles) {
        if (a != "*") {
            seqs.push_back(&a);
        }
    }
    if (seqs.empty()) {
        return 0;
    }
    size_t min_len = seqs[0]->length();
    size_t max_len = 0;
    for (const string* s : seqs) {
        min_len = std::min(min_len, s->length());
        max_len = std::max(max_len, s->length());
    }
    // The prefix and the suffix may not overlap, exactly as in flatten_common_allele_ends: a shared
    // region can only be counted once, or {"AAAA","AAAAA"} would come out at -3 instead of 1.
    // Case-insensitive to match flatten's own toupper and deconstruct's toUppercase.
    auto shared = [&](size_t skip, bool from_back) {
        auto at = [&](const string* s, size_t i) {
            return std::toupper((*s)[from_back ? s->length() - 1 - i : i]);
        };
        size_t n = 0;
        while (skip + n < min_len) {
            int ch = at(seqs[0], n);
            bool match = true;
            for (size_t j = 1; j < seqs.size() && match; ++j) {
                match = at(seqs[j], n) == ch;
            }
            if (!match) {
                break;
            }
            ++n;
        }
        return n;
    };
    size_t prefix = shared(0, false);
    size_t suffix = shared(prefix, true);
    // non-negative structurally, not by clamping: the loop caps give prefix + suffix <= min_len <=
    // max_len
    return (int64_t)(max_len - prefix - suffix);
}

bool merge_similar_alleles(const PathPositionHandleGraph& graph, size_t allele_count,
                           const function<bool(size_t allele, Traversal& walk)>& allele_walk,
                           vector<int>& site_genotype, const string& sample_name,
                           vcflib::Variant& out_variant, GLLayout gl_layout, double threshold,
                           int64_t min_len) {
    if (!(threshold < 1.0)) {
        return false;
    }
    // we only collapse a genotype that actually calls two distinct ALTs.  This also keeps the -a
    // padding block above (which adds uncalled alleles) out of scope: those alleles are advertised,
    // not called, and rewriting them without a genotype change would be a silent surprise.
    set<int> called_alts;
    for (int g : site_genotype) {
        if (g > 0) {
            called_alts.insert(g);
        }
    }
    if (called_alts.size() < 2) {
        return false;
    }
    // Per-site gate, decided over the alleles this record actually emits.  NOT over the traversal
    // finder's candidate list: that is up to max_yens_traversals (50) speculative paths, most of
    // which never become an allele, so gating on them lets an invisible branch with no reads and no
    // AT entry decide whether merging happens.  deconstruct's equivalent gate is decided over the
    // set that becomes ITS alleles -- the reference plus everything get_traversal_order kept.  Both
    // tools gate on what they emit, but those sets differ (we see only the called genotype's
    // alleles, deconstruct sees every haplotype), so the two can disagree at a site whose uncalled
    // haplotypes are much larger than its called ones.
    // The quantity is CORE LENGTH (see allele_core_length): the longest allele once the prefix and
    // suffix shared by every allele are stripped.  Raw string length would answer differently from
    // vg deconstruct on the same variant, because this record has been flattened down to an anchor
    // base and deconstruct's has not.
    if (min_len > 0 &&
        allele_core_length(out_variant.alleles) < min_len) {
        return false;
    }

    // ALT-vs-ALT only.  Absorbing an ALT into allele 0 would empty out_variant.alt and the record
    // would then be dropped entirely by the caller, turning a het call into no call at all.
    // (vg deconstruct does fold near-reference alleles into the reference cluster and drop the
    //  record; that is deliberate there and deliberately not copied here.)
    vector<Traversal> alt_travs;
    vector<int> alt_to_allele;
    for (size_t i = 1; i < allele_count; ++i) {
        if (!called_alts.count((int)i)) {
            continue;
        }
        Traversal trav;
        if (!allele_walk(i, trav)) {
            // star placeholder or a child-snarl visit: leave this allele alone
            continue;
        }
        alt_travs.push_back(std::move(trav));
        alt_to_allele.push_back((int)i);
    }
    if (alt_travs.size() < 2) {
        return false;
    }

    // same clustering call deconstruct makes, so the metric, the endpoint pruning and the
    // >= comparison are inherited rather than reimplemented.
    //
    // Cluster in descending allele-depth order, so each cluster's head -- the allele that survives,
    // and the one MAT's similarity is measured against -- is its best-supported member.  Identity
    // order would instead inherit the traversal finder's ranking, which
    // FlowCaller::call_snarl_internal (and NestedFlowCaller's copy of it) switches to
    // length-weighted average flow once a snarl's interior passes the average-support threshold.
    // That ranking can put a short, lightly-supported allele ahead of a long, heavily-supported
    // one, and merging into it emits the minority sequence as a homozygous call carrying the pooled
    // depth.
    vector<int> order(alt_travs.size());
    std::iota(order.begin(), order.end(), 0);
    {
        auto& sample_fields = out_variant.samples[sample_name];
        auto ad_it = sample_fields.find("AD");
        if (ad_it != sample_fields.end() && ad_it->second.size() == out_variant.alleles.size()) {
            vector<double> ad(alt_travs.size(), 0);
            bool usable = true;
            for (size_t k = 0; k < alt_travs.size() && usable; ++k) {
                usable = parse_vcf_double(ad_it->second.at(alt_to_allele[k]), ad[k]);
            }
            if (usable) {
                // stable, so equal depths keep the finder's own ranking
                std::stable_sort(order.begin(), order.end(),
                                 [&](int a, int b) { return ad[a] > ad[b]; });
            }
        }
    }
    // The VCF reference allele sets the scale a pure deletion is measured against.  It is
    // site_traversals[0] and is deliberately not clustered (the loop above starts at 1).
    // The nullptr fallback is currently unreachable: the only producer of a Visit without a node is
    // NestedFlowCaller, and --bottom-up is rejected with -L.  It is kept because the consequence of
    // being wrong about that is a crash, and because falling back to pairwise scoring can only
    // merge less, never more.
    Traversal ref_trav;
    const Traversal* site_ref_trav = nullptr;
    if (allele_count > 0 && allele_walk(0, ref_trav)) {
        site_ref_trav = &ref_trav;
    }
    vector<pair<double, int64_t>> cluster_info;
    vector<int> unused_child_snarl_mapping;
    vector<vector<int>> clusters = cluster_traversals(&graph, alt_travs, order,
                                                      vector<pair<handle_t, handle_t>>(),
                                                      threshold,
                                                      cluster_info, unused_child_snarl_mapping,
                                                      site_ref_trav);

    // merge_to[a] == a for a surviving allele, else the allele it collapses into
    vector<int> merge_to(out_variant.alleles.size());
    std::iota(merge_to.begin(), merge_to.end(), 0);
    vector<string> mat_entries;
    bool merged_any = false;
    for (const vector<int>& cluster : clusters) {
        if (cluster.size() < 2) {
            continue;
        }
        int survivor = alt_to_allele[cluster.front()];
        for (size_t j = 1; j < cluster.size(); ++j) {
            int absorbed = alt_to_allele[cluster[j]];
            merge_to[absorbed] = survivor;
            merged_any = true;
            stringstream ss;
            ss.precision(3);
            ss << absorbed << ">" << survivor << ":" << cluster_info[cluster[j]].first;
            mat_entries.push_back(ss.str());
        }
    }
    if (!merged_any) {
        return false;
    }

    // dense renumbering of the survivors, preserving order
    vector<int> new_index(merge_to.size(), -1);
    int next = 0;
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_index[a] = next++;
        }
    }
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] != (int)a) {
            new_index[a] = new_index[merge_to[a]];
        }
    }
    int n_new = next;

    // alleles / alt
    vector<string> new_alleles(n_new);
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_alleles[new_index[a]] = out_variant.alleles[a];
        }
    }
    out_variant.alleles = new_alleles;
    out_variant.alt.assign(new_alleles.begin() + 1, new_alleles.end());

    // AT is Number=R, so it is indexed by allele just like alleles
    auto at_it = out_variant.info.find("AT");
    if (at_it != out_variant.info.end() && at_it->second.size() == merge_to.size()) {
        vector<string> new_at(n_new);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            if (merge_to[a] == (int)a) {
                new_at[new_index[a]] = at_it->second[a];
            }
        }
        at_it->second = new_at;
    }

    auto& sample = out_variant.samples[sample_name];

    // AD is a per-allele count, so the absorbed allele's reads move onto the survivor.  sum(AD) is
    // therefore unchanged by the merge; it can slightly over-count when the merged alleles share
    // interior nodes, whose depth was proportionally split between them.
    auto ad_it = sample.find("AD");
    if (ad_it != sample.end() && ad_it->second.size() == merge_to.size()) {
        vector<double> summed(n_new, 0);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            double v = 0;
            // treat an unparseable entry as 0 rather than bailing: the merge is already committed,
            // and dropping AD would leave the record with no per-allele depth at all
            parse_vcf_double(ad_it->second[a], v);
            summed[new_index[a]] += v;
        }
        vector<string> new_ad(n_new);
        for (int a = 0; a < n_new; ++a) {
            new_ad[a] = std::to_string((int64_t)std::llround(summed[a]));
        }
        ad_it->second = new_ad;
        // MAD is the min allele depth over the called alleles; recompute so it agrees with the AD
        // and GT printed beside it.  FILTER is left as computed pre-merge, and since the new MAD is
        // >= the old one, a depth filter can only over-filter, never under-filter.
        auto mad_it = sample.find("MAD");
        if (mad_it != sample.end() && mad_it->second.size() == 1) {
            double min_ad = -1;
            for (int g : site_genotype) {
                if (g >= 0 && g < (int)merge_to.size()) {
                    double v = summed[new_index[g]];
                    if (min_ad < 0 || v < min_ad) {
                        min_ad = v;
                    }
                }
            }
            if (min_ad >= 0) {
                mad_it->second[0] = std::to_string((int64_t)std::llround(min_ad));
            }
        }
    }

    // GL is Number=G, and its order depends on the caller that wrote it (see GLLayout), so the
    // layout is passed in rather than assumed.
    auto gl_it = sample.find("GL");
    if (gl_it != sample.end()) {
        size_t n_old = merge_to.size();
        // Only the diploid layout can occur: a merge needs two different called ALTs, so the
        // ploidy is 2.
        assert(gl_it->second.size() == n_old * (n_old + 1) / 2);
        bool gl_usable = true;
        vector<double> old_gl(gl_it->second.size(), 0.0);
        for (size_t g = 0; g < gl_it->second.size() && gl_usable; ++g) {
            gl_usable = parse_vcf_double(gl_it->second[g], old_gl[g]);
        }
        vector<double> folded;
        if (gl_usable) {
            folded = fold_genotype_likelihoods(old_gl, new_index, (size_t)n_new, gl_layout);
        }
        if (gl_usable) {
            vector<string> new_gl(folded.size());
            for (size_t i = 0; i < folded.size(); ++i) {
                new_gl[i] = std::to_string(folded[i]);
            }
            gl_it->second = new_gl;
        } else {
            // A value we cannot parse.  Leaving GL alone would emit a Number=G field whose length
            // disagrees with the new allele count, so drop it rather than lie.
            sample.erase(gl_it);
            auto& fmt = out_variant.format;
            fmt.erase(std::remove(fmt.begin(), fmt.end(), string("GL")), fmt.end());
        }
    }
    // GQ and GP are deliberately untouched: they come from the caller's own CallInfo, computed over
    // its candidate set rather than from the emitted GL, so recomputing them here would silently
    // swap one statistic for another.

    // GT, from the renumbered genotype
    for (int& g : site_genotype) {
        if (g >= 0 && g < (int)merge_to.size()) {
            g = new_index[g];
        }
    }
    stringstream vcf_gt;
    for (size_t i = 0; i < site_genotype.size(); ++i) {
        if (site_genotype[i] == MISSING_ALLELE_MARKER) {
            vcf_gt << ".";
        } else {
            vcf_gt << site_genotype[i];
        }
        if (i != site_genotype.size() - 1) {
            vcf_gt << "/";
        }
    }
    sample["GT"] = {vcf_gt.str()};

    // record what happened: without this a merged 1/1 is indistinguishable from a real hom-alt
    out_variant.info["MAT"] = mat_entries;

    out_variant.updateAlleleIndexes();
    return true;
}

void add_allele_path_to_info(vcflib::Variant& v, int allele, const vector<NodeVisit>& visits,
                             bool reversed, const NodeTranslation* translation) {
    auto& trav_info = v.info["AT"];
    assert(allele < trav_info.size());

    const NodeVisit* prev_visit = nullptr;
    NodeTranslation::const_iterator prev_trans;
    
    for (size_t i = 0; i < visits.size(); ++i) {
        size_t j = !reversed ? i : visits.size() - 1 - i;
        const NodeVisit& visit = visits[j];
        nid_t node_id = visit.first;
        string node_name = std::to_string(node_id);
        bool skip = false;
        if (translation) {
            auto i = translation->find(node_id);
            if (i == translation->end()) {
                throw runtime_error("Error [vg deconstruct]: Unable to find node " + node_name + " in translation file");
            }
            if (prev_visit) {
                nid_t prev_node_id = prev_visit->first;
                if (prev_trans->second.first == i->second.first && node_id != prev_node_id) {
                    // here is a case where we have two consecutive nodes that map back to
                    // the same source node.
                    // todo: check if translation node properly covered
                    skip = true;
                }
            }
            node_name = i->second.first;
            prev_trans = i;
        }

        if (!skip) {
            bool vrev = visit.second != reversed;
            trav_info[allele] += (vrev ? "<" : ">");
            trav_info[allele] += node_name;
        }
        prev_visit = &visit;
    }
    if (trav_info[allele].empty()) {
        // note: * alleles get empty traversals
        trav_info[allele] = ".";
    }
}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(
    const PathPositionHandleGraph& graph, const handle_t& start, const handle_t& end,
    const string& ref_path_name) {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    handle_t start_handle = start;
    map<size_t, step_handle_t> start_steps;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                start_steps[graph.get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = end;
    map<size_t, step_handle_t> end_steps;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                end_steps[graph.get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    // just because we found a pair of steps on our path that correspond to the snarl ends, doesn't
    // mean the path threads the snarl.  verify that we can actaully walk, either forwards or backwards
    // along the path from the start node and hit then end node in the right orientation. 
    bool start_rev = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != graph.get_is_reverse(start);
    bool end_rev = graph.get_is_reverse(graph.get_handle_of_step(end_step)) != graph.get_is_reverse(end);
    bool found_end = start_rev == end_rev && start_rev == start_steps.begin()->first > end_steps.begin()->first;
        
    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        found_end = false;
        // try each start step
        for (auto i = start_steps.begin(); i != start_steps.end() && !found_end; ++i) {
            start_step = i->second;
            bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != graph.get_is_reverse(start);
            if (scan_backward) {
                // if we're going backward, we expect to reach the end backward
                end_handle = graph.flip(end);
            }            
            if (scan_backward) {
                for (step_handle_t cur_step = start_step; graph.has_previous_step(cur_step) && !found_end;
                     cur_step = graph.get_previous_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            } else {
                for (step_handle_t cur_step = start_step; graph.has_next_step(cur_step) && !found_end;
                     cur_step = graph.get_next_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            }
        }
    }
    int64_t start_position = start_steps.begin()->first;
    step_handle_t out_start_step = start_step;
    int64_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : graph.get_position_of_step(end_step);
    step_handle_t out_end_step = end_step == end_steps.begin()->second ? end_steps.begin()->second : end_step;
    bool backward = end_position < start_position;
    

    if (!found_end) {
        // oops, once of the above checks failed.  we tell caller we coudlnt find by hacking in a -1
        // coordinate.
        start_position = -1;
        end_position = -1;
    }

    if (backward) {
        return make_tuple(end_position, start_position, backward, out_end_step, out_start_step);
    } else {
        return make_tuple(start_position, end_position, backward, out_start_step, out_end_step);
    }
}

pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const handle_t& start,
                                       const handle_t& end, const string& ref_path_name,
                                       int64_t ref_path_offset) {
    const string basepath_name = Paths::strip_subrange(ref_path_name);
    const int64_t position = base_path_position(
        ref_path_name, get<0>(get_ref_interval(graph, start, end, ref_path_name)) + ref_path_offset);
    return make_pair(basepath_name, position);
}

void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) {
    if (variant.alt.size() == 0) {
        return;
    }

    // find the minimum allele length to make sure we don't delete an entire allele
    size_t min_allele_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_allele_len = std::min(min_allele_len, variant.alleles[i].length());
    }

    // With min_allele_len 0 the decrement below would wrap max_flatten_len, a size_t, and the
    // backward pass would read past the end of an empty allele, such as a pure deletion's ALT.
    if (min_allele_len == 0) {
        return;
    }

    // the maximum number of bases we want ot zip up, applying override if provided
    size_t max_flatten_len = len_override > 0 ? len_override : min_allele_len;
    
    // want to leave at least one in the reference position
    if (max_flatten_len == min_allele_len) {
        --max_flatten_len;
    }
    
    bool match = true;
    int shared_prefix_len = 0;
    for (int i = 0; i < max_flatten_len && match; ++i) {
        char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
        for (int j = 1; j < variant.alleles.size() && match; ++j) {
            char c2 = std::toupper(variant.alleles[j][!backward ? i : variant.alleles[j].length() - 1 - i]);
            match = c1 == c2;
        }
        if (match) {
            ++shared_prefix_len;
        }
    }

    if (!backward) {
        variant.position += shared_prefix_len;
    }
    for (int i = 0; i < variant.alleles.size(); ++i) {
        if (!backward) {
            variant.alleles[i] = variant.alleles[i].substr(shared_prefix_len);
        } else {
            variant.alleles[i] = variant.alleles[i].substr(0, variant.alleles[i].length() - shared_prefix_len);
        }
        if (i == 0) {
            variant.ref = variant.alleles[i];
        } else {
            variant.alt[i - 1] = variant.alleles[i];
        }
    }
}

int64_t base_path_position(const string& ref_path_name, int64_t along_path) {
    subrange_t subrange;
    Paths::strip_subrange(ref_path_name, &subrange);
    const int64_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : (int64_t)subrange.first;
    return along_path + 1 + basepath_offset;
}

string site_name(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                 const NodeTranslation* translation, bool in_brackets) {
    return site_name(graph.get_id(start), graph.get_is_reverse(start), graph.get_id(end),
                     graph.get_is_reverse(end), translation, in_brackets);
}

size_t record_key_of(const string& site_id) {
    return std::hash<string>{}(site_id);
}

size_t record_key_of(const HandleGraph& graph, const handle_t& start, const handle_t& end,
                     const NodeTranslation* translation) {
    return record_key_of(site_name(graph, start, end, translation, false));
}

string site_name(nid_t start_node_id, bool start_backward, nid_t end_node_id, bool end_backward,
                 const NodeTranslation* translation, bool in_brackets) {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    string start_node = std::to_string(start_node_id);
    string end_node = std::to_string(end_node_id);
    if (translation) {
        auto i = translation->find(start_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + start_node + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + end_node + " in translation file");
        }
        end_node = i->second.first;
    }
    // Built in place rather than through a stringstream: names are printed for every site of the
    // graph when the VCF is written, and for each ancestor of every record.
    string name;
    name.reserve(start_node.size() + end_node.size() + 4);
    if (in_brackets) {
        name += '(';
    }
    name += start_backward ? '<' : '>';
    name += start_node;
    name += end_backward ? '<' : '>';
    name += end_node;
    if (in_brackets) {
        name += ')';
    }
    return name;
}

SiteRecord build_site_record(const PathPositionHandleGraph& graph, const SiteToWrite& site,
                             const SiteAlleles& alleles, const SiteHooks& hooks,
                             const RecordOptions& options) {
    const vector<int>& genotype = site.genotype;
    const int ref_trav_idx = site.ref_trav_idx;

    SiteRecord record;
    vcflib::Variant& out_variant = record.variant;

    // The traversal of each VCF allele, or -1 for the "*" placeholder.
    vector<int> site_trav = {ref_trav_idx};
    vector<int>& site_genotype = record.genotype;
    auto ref_gt_it = std::find(genotype.begin(), genotype.end(), ref_trav_idx);
    out_variant.ref = alleles.spell(ref_trav_idx,
                                    ref_gt_it != genotype.end() ? ref_gt_it - genotype.begin() : 0);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;
    // Which VCF allele each called traversal became. Alleles with the same sequence are merged, and
    // only called traversals are written.
    map<int, int>& trav_to_allele = record.trav_to_allele;
    allele_to_gt[out_variant.ref] = 0;
    trav_to_allele[ref_trav_idx] = 0;
    int star_allele_idx = -1;  // index for star allele in allele_to_gt, if needed
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == STAR_ALLELE_MARKER) {
            // Star allele: haplotype spans this site but has no defined traversal here
            if (star_allele_idx < 0) {
                // Add star allele to allele list
                star_allele_idx = allele_to_gt.size();
                allele_to_gt["*"] = star_allele_idx;
                // Add a placeholder traversal (written as "." in AT)
                site_trav.push_back(-1);
            }
            site_genotype.push_back(star_allele_idx);
        } else if (genotype[i] == MISSING_ALLELE_MARKER) {
            // Missing allele: parent doesn't traverse this child, output as '.' in VCF
            site_genotype.push_back(MISSING_ALLELE_MARKER);
        } else if (genotype[i] == ref_trav_idx
                   || (hooks.same_as_reference && hooks.same_as_reference(genotype[i]))) {
            // The reference traversal, or one that takes the same route through this site and
            // differs only inside child sites, whose own records report those differences.
            site_genotype.push_back(0);
            trav_to_allele[genotype[i]] = 0;
        } else {
            string allele_string = alleles.spell(genotype[i], i);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_trav.push_back(genotype[i]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
            trav_to_allele[genotype[i]] = site_genotype.back();
        }
    }

    // add on fixed number of uncalled traversals if we're making a ref-call
    // with genotype_snarls set to true
    if (site.genotype_snarls && site_trav.size() <= 1) {
        // note: we're adding all the strings here and sorting to make this deterministic
        // at the cost of speed
        map<string, int> allele_map;
        for (int i = 0; i < site.traversal_count; ++i) {
            // todo: verify index below.  it's for uncalled traversals so not important tho
            string allele_string = alleles.spell(i, max(0, (int)genotype.size() - 1));
            if (!allele_map.count(allele_string)) {
                allele_map[allele_string] = i;
            }
        }
        // pick out the first "max_uncalled_alleles" traversals to add
        int i = 0;
        for (auto ai = allele_map.begin(); i < options.max_uncalled_alleles && ai != allele_map.end(); ++i, ++ai) {
            if (!allele_to_gt.count(ai->first)) {
                allele_to_gt[ai->first] = allele_to_gt.size();
                site_trav.push_back(ai->second);
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    
    // init the traversal info
    out_variant.info["AT"].resize(allele_to_gt.size());

    for (auto& allele_gt : allele_to_gt) {
#ifdef debug
        cerr << "allele " << allele_gt.first << " -> gt " << allele_gt.second << endl;
#endif
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;

        // update the traversal info
        const int trav = site_trav.at(allele_gt.second);
        add_allele_path_to_info(out_variant, allele_gt.second,
                                trav >= 0 ? alleles.visits(trav) : vector<NodeVisit>(), false,
                                options.translation);
    }

    // resolve subpath naming
    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(site.ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    // in VCF we usually just want a contig
    string contig_name = PathMetadata::parse_locus_name(basepath_name);
    if (contig_name != PathMetadata::NO_LOCUS_NAME) {
        basepath_name = contig_name;
    }
    // fill out the rest of the variant    
    out_variant.sequenceName = basepath_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get<0>(get_ref_interval(graph, site.start, site.end, site.ref_path_name))
                           + site.ref_offset + 1 + basepath_offset;
    // Kept before flattening moves it: the position of the reference traversal's first base, from
    // which block offsets are measured.
    record.unflattened_position = out_variant.position;
    out_variant.id = site.id != nullptr
                         ? *site.id
                         : site_name(graph, site.start, site.end, options.translation, false);
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[options.sample_name]["GT"];
    
    stringstream vcf_gt;
    if (!genotype.empty()) {
        for (int i = 0; i < site_genotype.size(); ++i) {
            if (site_genotype[i] == MISSING_ALLELE_MARKER) {
                vcf_gt << ".";
            } else {
                vcf_gt << site_genotype[i];
            }
            if (i != site_genotype.size() - 1) {
                vcf_gt << "/";
            }
        }
    } else {
        for (int i = 0; i < site.ploidy; ++i) {
            vcf_gt << ".";
            if (i != site.ploidy - 1) {
                vcf_gt << "/";
            }
        }
    }
                    
    genotype_vector.push_back(vcf_gt.str());

    // Phase the genotype here, where `trav_to_allele` maps traversals to this record's alleles.
    int64_t phase_set_to_write = -1;
    if (!genotype.empty() && hooks.phase) {
        phase_set_to_write = hooks.phase(site_genotype, trav_to_allele, genotype_vector[0]);
    }

    // add some support info
    if (hooks.fill_info) {
        hooks.fill_info(site_trav, site_genotype, out_variant);
    }

    // PS last in FORMAT.
    if (phase_set_to_write >= 0) {
        out_variant.format.push_back("PS");
        out_variant.samples[options.sample_name]["PS"].push_back(std::to_string(phase_set_to_write));
    }

    // if genotype_snarls, then we only flatten up to the snarl endpoints
    // (this is when we are in genotyping mode and want consistent calls regardless of the sample)
    int64_t flatten_len_s = 0;
    int64_t flatten_len_e = 0;
    if (site.genotype_snarls) {
        flatten_len_s = graph.get_length(site.start);
        assert(flatten_len_s >= 0);
        flatten_len_e = graph.get_length(site.end);
    }
    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true, flatten_len_e);
    flatten_common_allele_ends(out_variant, false, flatten_len_s);

    // Merge near-identical called ALT alleles, turning 1/2 into 1/1. After the caller's fields, so
    // the genotyper saw every candidate, and after flattening, so the surviving allele's string,
    // POS and REF are the same as without merging. The missing-allele fixup below runs after it and
    // is unaffected, since merging never empties alt.
    record.alleles_merged = merge_similar_alleles(
        graph, site_trav.size(),
        [&](size_t allele, Traversal& walk) {
            const int trav = site_trav[allele];
            return visits_to_walk(graph, trav >= 0 ? alleles.visits(trav) : vector<NodeVisit>(),
                                  walk);
        },
        site_genotype, options.sample_name, out_variant, hooks.gl_layout,
        options.allele_merge_threshold, options.allele_merge_min_len);

    // If genotype contains missing allele but no ALT, add * as ALT to emit valid VCF
    // This happens when one parent haplotype doesn't traverse a nested child snarl
    bool has_missing = std::find(site_genotype.begin(), site_genotype.end(), MISSING_ALLELE_MARKER) != site_genotype.end();
    if (has_missing && out_variant.alt.empty()) {
        out_variant.alt.push_back("*");
        out_variant.alleles.push_back("*");
        out_variant.info["AT"].push_back(".");
        // AD is Number=R, written by the caller's fields for the alleles before this one was added,
        // so it gets an entry for the new allele. GL needs none: it is left out whenever the
        // genotype has a marker, as MISSING_ALLELE_MARKER is.
        auto ad_it = out_variant.samples[options.sample_name].find("AD");
        if (ad_it != out_variant.samples[options.sample_name].end()) {
            ad_it->second.push_back("0");
        }
    }

    return record;
}

string block_site_name(const string& id) {
    size_t underscore = id.rfind('_');
    return underscore == string::npos ? id : id.substr(0, underscore);
}

string nested_strand_genotype(int allele, int strand) {
    const string a = std::to_string(allele);
    return strand == 0 ? a + "|." : "." + ("|" + a);
}

vector<NodeVisit> visits_of(const SnarlTraversal& trav) {
    vector<NodeVisit> visits;
    visits.reserve(trav.visit_size());
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        visits.emplace_back(visit.node_id(), visit.backward());
    }
    return visits;
}

vector<NodeVisit> visits_of(const HandleGraph& graph, const Traversal& walk) {
    vector<NodeVisit> visits;
    visits.reserve(walk.size());
    for (const handle_t& handle : walk) {
        visits.emplace_back(graph.get_id(handle), graph.get_is_reverse(handle));
    }
    return visits;
}

}
