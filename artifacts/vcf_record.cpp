#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "vcf_record.hpp"
#include "../src/path.hpp"
#include "../src/site_emitter.hpp"
#include "../src/traversal_clusters.hpp"

namespace vg {

bool visits_to_walk(const HandleGraph& graph, const vector<NodeVisit>& visits, Traversal& walk) {
    // cluster_traversals needs at least two handles, and a visit to a child site has no single
    // handle.
    if (visits.size() < 2) {
        return false;
    }
    walk.clear();
    walk.reserve(visits.size());
    for (const auto& [node_id, backward] : visits) {
        if (node_id <= 0) {
            return false;
        }
        walk.push_back(graph.get_handle(node_id, backward));
    }
    return true;
}

/// Return the GT text of an unphased genotype in VCF allele numbers, with "." for
/// MISSING_ALLELE_MARKER.
string unphased_gt_text(const vector<int>& site_genotype) {
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
    return vcf_gt.str();
}

namespace {

/// Parse a VCF FORMAT value, or return false if it is not a number.
///
/// Missing values ("." and "") are ordinary input here. This never throws or exits, since
/// merge_similar_alleles runs inside an OpenMP region, where an escaping exception terminates the
/// program.
bool parse_vcf_double(const string& field, double& value) {
    try {
        size_t after;
        value = std::stod(field, &after);
        return after == field.size();
    } catch (const std::logic_error&) {
        // std::stod throws std::invalid_argument or std::out_of_range.
        return false;
    }
}

/// Return the order in which to cluster the called ALTs, deepest first.
///
/// `alt_to_allele` gives the VCF allele of each ALT, and the order is of indices into it. ALTs of
/// equal depth, or all of them if AD cannot be read, keep the order of `alt_to_allele`.
vector<int> order_by_depth(vcflib::Variant& variant, const string& sample_name,
                           const vector<int>& alt_to_allele) {
    vector<int> order(alt_to_allele.size());
    std::iota(order.begin(), order.end(), 0);
    auto& sample_fields = variant.samples[sample_name];
    auto ad_it = sample_fields.find("AD");
    if (ad_it == sample_fields.end() || ad_it->second.size() != variant.alleles.size()) {
        return order;
    }
    vector<double> ad(alt_to_allele.size(), 0);
    for (size_t k = 0; k < alt_to_allele.size(); ++k) {
        if (!parse_vcf_double(ad_it->second.at(alt_to_allele[k]), ad[k])) {
            return order;
        }
    }
    std::ranges::stable_sort(order, [&](int a, int b) { return ad[a] > ad[b]; });
    return order;
}

/// Return, for each of `allele_count` alleles, the allele it merges into, or itself if it survives.
///
/// Each of `clusters` that has two or more ALTs merges into its first ALT. The clusters hold
/// indices into `alt_to_allele`, which gives the VCF allele of each, and `cluster_info` gives the
/// similarity of each ALT to the first of its cluster. Each merge is added to `mat_entries` as
/// "absorbed>survivor:similarity", for INFO/MAT.
vector<int> merge_targets(const vector<vector<int>>& clusters, const vector<int>& alt_to_allele,
                          const vector<pair<double, int64_t>>& cluster_info, size_t allele_count,
                          vector<string>& mat_entries) {
    vector<int> merge_to(allele_count);
    std::iota(merge_to.begin(), merge_to.end(), 0);
    for (const vector<int>& cluster : clusters) {
        if (cluster.size() < 2) {
            continue;
        }
        int survivor = alt_to_allele[cluster.front()];
        for (size_t j = 1; j < cluster.size(); ++j) {
            int absorbed = alt_to_allele[cluster[j]];
            merge_to[absorbed] = survivor;
            stringstream ss;
            ss.precision(3);
            ss << absorbed << ">" << survivor << ":" << cluster_info[cluster[j]].first;
            mat_entries.push_back(ss.str());
        }
    }
    return merge_to;
}

/**
 * A renumbering of the alleles of a record that merges some alleles into others.
 *
 * The surviving alleles keep their order, and are numbered from 0 with no gaps. Each merged
 * allele takes the new number of the allele it merged into.
 */
class AlleleRenumbering {
public:
    /// Make the renumbering in which each old allele `a` merges into `merge_to[a]`.
    ///
    /// A surviving allele merges into itself, and every other allele merges into a survivor.
    explicit AlleleRenumbering(vector<int> merge_to)
        : merge_to(std::move(merge_to)), new_index(this->merge_to.size(), -1) {
        for (size_t a = 0; a < this->merge_to.size(); ++a) {
            if (survives(a)) {
                new_index[a] = new_count++;
            }
        }
        for (size_t a = 0; a < this->merge_to.size(); ++a) {
            if (!survives(a)) {
                new_index[a] = new_index[this->merge_to[a]];
            }
        }
    }

    /// Get how many alleles the record had before merging.
    size_t old_size() const { return merge_to.size(); }

    /// Get how many alleles the record has after merging.
    int new_size() const { return new_count; }

    /// Get the new number of each old allele.
    const vector<int>& new_numbers() const { return new_index; }

    /// Return whether old allele `a` survives the merge.
    bool survives(size_t a) const { return merge_to[a] == (int)a; }

    /// Return the entries of the surviving alleles, in their new order, from a column that has
    /// an entry for each old allele.
    vector<string> survivors_of(const vector<string>& column) const {
        vector<string> kept(new_count);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            if (survives(a)) {
                kept[new_index[a]] = column[a];
            }
        }
        return kept;
    }

private:
    /// For each old allele, the old allele it merges into.
    vector<int> merge_to;
    /// For each old allele, its new number.
    vector<int> new_index;
    /// How many alleles survive.
    int new_count = 0;
};

/// Move the reads of each merged allele onto its survivor in AD, and recompute MAD to match.
///
/// The total of AD stays the same. It can over-count a little where the merged alleles share
/// nodes, whose depth was split between them. `site_genotype` is in the old allele numbers.
void merge_allele_depths(map<string, vector<string>>& sample, const AlleleRenumbering& renumbering,
                         const vector<int>& site_genotype) {
    auto ad_it = sample.find("AD");
    if (ad_it == sample.end() || ad_it->second.size() != renumbering.old_size()) {
        return;
    }
    const vector<int>& new_index = renumbering.new_numbers();
    vector<double> summed(renumbering.new_size(), 0);
    for (size_t a = 0; a < renumbering.old_size(); ++a) {
        double v = 0;
        // An entry that cannot be read counts as 0. The merge is already decided, and dropping
        // AD would leave the record with no depth for each allele.
        parse_vcf_double(ad_it->second[a], v);
        summed[new_index[a]] += v;
    }
    vector<string> new_ad(renumbering.new_size());
    for (int a = 0; a < renumbering.new_size(); ++a) {
        new_ad[a] = std::to_string((int64_t)std::llround(summed[a]));
    }
    ad_it->second = new_ad;
    // MAD is the least depth of a called allele. FILTER keeps the value it had before the merge,
    // and MAD can only grow, so a depth filter can only over-filter.
    auto mad_it = sample.find("MAD");
    if (mad_it == sample.end() || mad_it->second.size() != 1) {
        return;
    }
    double min_ad = -1;
    for (int g : site_genotype) {
        if (g >= 0 && g < (int)renumbering.old_size()) {
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

/// Fold GL onto the merged alleles, or drop GL if a value in it cannot be read.
///
/// `gl_layout` is the order the caller wrote GL in. GL is dropped rather than left alone, since
/// its length would disagree with the new number of alleles.
void fold_allele_likelihoods(vcflib::Variant& variant, map<string, vector<string>>& sample,
                             const AlleleRenumbering& renumbering, GLLayout gl_layout) {
    auto gl_it = sample.find("GL");
    if (gl_it == sample.end()) {
        return;
    }
    size_t n_old = renumbering.old_size();
    // A merge needs two different called ALTs, so the genotype is diploid.
    assert(gl_it->second.size() == n_old * (n_old + 1) / 2);
    vector<double> old_gl(gl_it->second.size(), 0.0);
    for (size_t g = 0; g < gl_it->second.size(); ++g) {
        if (!parse_vcf_double(gl_it->second[g], old_gl[g])) {
            sample.erase(gl_it);
            std::erase(variant.format, string("GL"));
            return;
        }
    }
    vector<double> folded = fold_genotype_likelihoods(old_gl, renumbering.new_numbers(),
                                                      (size_t)renumbering.new_size(), gl_layout);
    vector<string> new_gl(folded.size());
    for (size_t i = 0; i < folded.size(); ++i) {
        new_gl[i] = std::to_string(folded[i]);
    }
    gl_it->second = new_gl;
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
    // Count the shared bases from one end, case-insensitively, as flatten_common_allele_ends
    // does. The prefix and the suffix may not overlap, or {"AAAA","AAAAA"} would measure -3
    // instead of 1.
    auto shared = [&](size_t skip, bool from_back) {
        auto at = [&](const string* s, size_t i) {
            return std::toupper((*s)[from_back ? s->length() - 1 - i : i]);
        };
        size_t n = 0;
        while (skip + n < min_len) {
            if (int ch = at(seqs[0], n); !std::all_of(seqs.begin() + 1, seqs.end(),
                             [&](const string* s) { return at(s, n) == ch; })) {
                break;
            }
            ++n;
        }
        return n;
    };
    size_t prefix = shared(0, false);
    size_t suffix = shared(prefix, true);
    // prefix + suffix <= min_len <= max_len, so this is not negative.
    return (int64_t)(max_len - prefix - suffix);
}

bool merge_similar_alleles(const PathPositionHandleGraph& graph,
                           const vector<vector<NodeVisit>>& allele_visits,
                           vector<int>& site_genotype, const string& sample_name,
                           vcflib::Variant& out_variant, GLLayout gl_layout, double threshold,
                           int64_t min_len) {
    if (!(threshold < 1.0)) {
        return false;
    }
    // Merge only in a genotype that calls two different ALTs. The alleles that were not called,
    // added to a reference call written anyway, are left as they are.
    set<int> called_alts;
    for (int g : site_genotype) {
        if (g > 0) {
            called_alts.insert(g);
        }
    }
    if (called_alts.size() < 2) {
        return false;
    }
    // Measure the alleles this record writes, by core length, as vg deconstruct measures the
    // alleles it writes, so that the two agree on a variant written by both. The traversal
    // finder's other candidates have no reads and no AT entry, so they do not count. Raw string
    // length would disagree, since this record is trimmed to an anchor base and deconstruct's is
    // not.
    if (min_len > 0 &&
        allele_core_length(out_variant.alleles) < min_len) {
        return false;
    }

    // Cluster only the called ALTs. Merging an ALT into the reference allele could leave the
    // record with no ALT, which the caller would drop, turning a het call into no call.
    vector<Traversal> alt_travs;
    vector<int> alt_to_allele;
    for (size_t i = 1; i < allele_visits.size(); ++i) {
        if (!called_alts.contains((int)i)) {
            continue;
        }
        Traversal trav;
        if (!visits_to_walk(graph, allele_visits[i], trav)) {
            continue;
        }
        alt_travs.push_back(std::move(trav));
        alt_to_allele.push_back((int)i);
    }
    if (alt_travs.size() < 2) {
        return false;
    }

    // Cluster the ALTs deepest first, as vg deconstruct does, so that each cluster merges into its
    // best-supported member, and MAT measures against it. The order of the traversal finder can
    // put a short, lightly supported allele first, and merging into it would write the minority
    // sequence as a homozygous call with the pooled depth.
    const vector<int> order = order_by_depth(out_variant, sample_name, alt_to_allele);
    // The reference allele sets the scale a pure deletion is measured against. Without a walk
    // for it, cluster_traversals scores pairs alone, which can only merge less.
    Traversal ref_trav;
    const Traversal* site_ref_trav = nullptr;
    if (!allele_visits.empty() && visits_to_walk(graph, allele_visits[0], ref_trav)) {
        site_ref_trav = &ref_trav;
    }
    vector<pair<double, int64_t>> cluster_info;
    vector<int> unused_child_snarl_mapping;
    vector<vector<int>> clusters = cluster_traversals(&graph, alt_travs, order,
                                                      vector<pair<handle_t, handle_t>>(),
                                                      threshold,
                                                      cluster_info, unused_child_snarl_mapping,
                                                      site_ref_trav);

    vector<string> mat_entries;
    vector<int> merge_to = merge_targets(clusters, alt_to_allele, cluster_info,
                                         out_variant.alleles.size(), mat_entries);
    if (mat_entries.empty()) {
        return false;
    }
    const AlleleRenumbering renumbering(std::move(merge_to));

    // Rewrite each field that has an entry per allele.
    out_variant.alleles = renumbering.survivors_of(out_variant.alleles);
    out_variant.alt.assign(out_variant.alleles.begin() + 1, out_variant.alleles.end());
    if (auto at_it = out_variant.info.find("AT");
        at_it != out_variant.info.end() && at_it->second.size() == renumbering.old_size()) {
        at_it->second = renumbering.survivors_of(at_it->second);
    }
    auto& sample = out_variant.samples[sample_name];
    merge_allele_depths(sample, renumbering, site_genotype);
    fold_allele_likelihoods(out_variant, sample, renumbering, gl_layout);
    // GQ and GP stay as they are. The caller computed them over its own candidates, not from GL.

    for (int& g : site_genotype) {
        if (g >= 0 && g < (int)renumbering.old_size()) {
            g = renumbering.new_numbers()[g];
        }
    }
    sample["GT"] = {unphased_gt_text(site_genotype)};

    // Record the merge, since a merged 1/1 otherwise looks like a real homozygous ALT call.
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
            auto found = translation->find(node_id);
            if (found == translation->end()) {
                throw runtime_error("Error [vg deconstruct]: Unable to find node " + node_name + " in translation file");
            }
            // Two consecutive nodes made from the same source node are written once.
            // TODO: Check that the pieces cover the source node.
            if (prev_visit && prev_trans->second.first == found->second.first
                && node_id != prev_visit->first) {
                skip = true;
            }
            node_name = found->second.first;
            prev_trans = found;
        }

        if (!skip) {
            bool vrev = visit.second != reversed;
            trav_info[allele] += (vrev ? "<" : ">");
            trav_info[allele] += node_name;
        }
        prev_visit = &visit;
    }
    if (trav_info[allele].empty()) {
        // A "*" allele visits nothing.
        trav_info[allele] = ".";
    }
}

namespace {

/// Return the steps of path `path_handle` at `handle`, by their positions on the path.
map<size_t, step_handle_t> path_steps_at(const PathPositionHandleGraph& graph,
                                         const handle_t& handle, path_handle_t path_handle) {
    map<size_t, step_handle_t> steps;
    graph.for_each_step_on_handle(handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                steps[graph.get_position_of_step(step)] = step;
            }
        });
    return steps;
}

/// Walk the path from step `from`, backward along it if `backward` is set, to the first step that
/// visits `target`, and put that step in `found`. Return false if no step does.
///
/// Both `from` and the last step of the path in the direction of the walk are compared.
bool find_step_along_path(const PathPositionHandleGraph& graph, step_handle_t from,
                          const handle_t& target, bool backward, step_handle_t& found) {
    step_handle_t cur = from;
    while (true) {
        if (graph.get_handle_of_step(cur) == target) {
            found = cur;
            return true;
        }
        if (backward ? !graph.has_previous_step(cur) : !graph.has_next_step(cur)) {
            return false;
        }
        cur = backward ? graph.get_previous_step(cur) : graph.get_next_step(cur);
    }
}

}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(
    const PathPositionHandleGraph& graph, const handle_t& start, const handle_t& end,
    const string& ref_path_name) {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    map<size_t, step_handle_t> start_steps = path_steps_at(graph, start, path_handle);
    map<size_t, step_handle_t> end_steps = path_steps_at(graph, end, path_handle);

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    // The path visits both bounds, but it threads the site only if it goes from the start to the
    // end in the orientations of the bounds, either forward or backward along the path.
    bool start_rev = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != graph.get_is_reverse(start);
    bool end_rev = graph.get_is_reverse(graph.get_handle_of_step(end_step)) != graph.get_is_reverse(end);
    bool found_end = start_rev == end_rev
                     && start_rev == (start_steps.begin()->first > end_steps.begin()->first);

    if (start_steps.size() > 1 || end_steps.size() > 1) {
        // The path visits a bound more than once, as on a cycle. Scan the path from each step at
        // the start for a step at the end.
        found_end = false;
        for (const auto& [position, step] : start_steps) {
            start_step = step;
            bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != graph.get_is_reverse(start);
            // Going backward, we reach the end backward.
            const handle_t end_handle = scan_backward ? graph.flip(end) : end;
            if (find_step_along_path(graph, start_step, end_handle, scan_backward, end_step)) {
                found_end = true;
                break;
            }
        }
    }
    int64_t start_position = graph.get_position_of_step(start_step);
    step_handle_t out_start_step = start_step;
    int64_t end_position = graph.get_position_of_step(end_step);
    step_handle_t out_end_step = end_step;
    bool backward = end_position < start_position;

    if (!found_end) {
        // The path does not thread the site, which the positions say by being -1.
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
    if (variant.alt.empty()) {
        return;
    }

    // Find the shortest allele, so that we never trim an allele away entirely.
    size_t min_allele_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_allele_len = std::min(min_allele_len, variant.alleles[i].length());
    }

    // An empty allele, such as the ALT of a pure deletion, has nothing to trim, and the
    // decrement below would wrap.
    if (min_allele_len == 0) {
        return;
    }

    // Trim at most `len_override` bases, and never more than all but one base of the shortest
    // allele.
    size_t max_flatten_len = min_allele_len - 1;
    if (len_override > 0) {
        max_flatten_len = std::min(max_flatten_len, len_override);
    }

    // Count the bases every allele shares at the end being trimmed, case-insensitively.
    int shared_prefix_len = 0;
    for (int i = 0; i < max_flatten_len; ++i) {
        char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
        if (!std::all_of(variant.alleles.begin() + 1, variant.alleles.end(),
                         [&](const string& allele) {
                char c2 = std::toupper(allele[!backward ? i : allele.length() - 1 - i]);
                return c1 == c2;
            })) {
            break;
        }
        ++shared_prefix_len;
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
            throw runtime_error("Error [vg]: Unable to find node " + start_node + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [vg]: Unable to find node " + end_node + " in translation file");
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

/// Set the contig, position and ID of a record for `site`.
void place_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                  const NodeTranslation* translation, SiteRecord& record) {
    vcflib::Variant& out_variant = record.variant;
    // Name the contig by the locus of the base path, where the path name has one.
    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(site.ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    if (string contig_name = PathMetadata::parse_locus_name(basepath_name);
        contig_name != PathMetadata::NO_LOCUS_NAME) {
        basepath_name = contig_name;
    }
    out_variant.sequenceName = basepath_name;
    // Add 1 to make the position 1-based, as VCF is.
    out_variant.position = get<0>(get_ref_interval(graph, site.start, site.end, site.ref_path_name))
                           + site.ref_offset + 1 + basepath_offset;
    record.unflattened_position = out_variant.position;
    out_variant.id = site.id != nullptr
                         ? *site.id
                         : site_name(graph, site.start, site.end, translation, false);
}

/// Return the GT text of a genotype with no call: `ploidy` "." separated by "/".
string no_call_gt_text(int ploidy) {
    stringstream vcf_gt;
    for (int i = 0; i < ploidy; ++i) {
        vcf_gt << ".";
        if (i != ploidy - 1) {
            vcf_gt << "/";
        }
    }
    return vcf_gt.str();
}

/// Add "*" as the ALT of a record whose genotype has "." and no ALT, so that it is valid VCF.
///
/// This happens where a strand of the parent does not visit a nested site.
void add_star_for_missing(vcflib::Variant& out_variant, const vector<int>& site_genotype,
                          const string& sample_name) {
    if (std::ranges::find(site_genotype, MISSING_ALLELE_MARKER) == site_genotype.end()
        || !out_variant.alt.empty()) {
        return;
    }
    out_variant.alt.emplace_back("*");
    out_variant.alleles.emplace_back("*");
    out_variant.info["AT"].emplace_back(".");
    // AD has an entry per allele, so it gets one for "*". GL needs none, since the caller leaves
    // it out of a genotype with an allele marker.
    auto ad_it = out_variant.samples[sample_name].find("AD");
    if (ad_it != out_variant.samples[sample_name].end()) {
        ad_it->second.emplace_back("0");
    }
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
