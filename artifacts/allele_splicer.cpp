#include <algorithm>

#include "allele_splicer.hpp"

namespace vg::multipass {

AlleleSplicer::AlleleSplicer(const HandleGraph& graph, const EditSnarlDecompositionOverlay& sites,
                             const PhaseTable& phase_table) :
    graph(graph), sites(sites), phase_table(phase_table) {
    // Nothing to do
}

void AlleleSplicer::index(StagedSiteTable& staged_sites, bool phasing) {
    children.clear();
    if (!phasing || !phase_table.has_rendered()) {
        return;
    }
    auto file = [&](const StagedSite& site) {
        if (site.placement.parent_record_key != 0 && !site.dropped) {
            children[site.placement.parent_record_key].push_back(&site);
        }
    };
    for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
        for (const StagedSite& site : staged_sites.queue(q)) {
            file(site);
        }
    }
    for (const StagedSite& site : staged_sites.unrendered()) {
        file(site);
    }
}

void AlleleSplicer::clear() {
    children.clear();
}

bool AlleleSplicer::has_children(size_t record_key) const {
    return children.count(record_key) != 0;
}

bool AlleleSplicer::splice(const StagedSite& site, SplicedAlleles& out) const {
    if (children.empty() || !children.count(site.record_key)) {
        return false;
    }
    const LinkageCollector::PhaseCall* phase = phase_table.rendered(site.record_key);
    if (phase == nullptr) {
        return false;
    }
    bool changed = false;
    const int strands = phase->ploidy == 1 ? 1 : 2;
    int first_trav = -1;
    for (int strand = 0; strand < strands; ++strand) {
        const int trav = strand == 0 ? phase->trav_first : phase->trav_second;
        if (trav < 0 || trav >= (int)site.alleles.walks.size() || trav == site.alleles.ref_allele
            || trav == first_trav) {
            // Leave the reference allele as it is. Where both strands carry the same allele,
            // splice it from strand 0 alone.
            continue;
        }
        first_trav = trav;
        Traversal walk;
        StepSequences steps;
        if (!spliced_allele(site, trav, strand, strands == 1, walk, steps)) {
            continue;
        }
        if (!changed) {
            out.travs = site.alleles.walks;
            out.sequences.assign(site.alleles.walks.size(), StepSequences());
            for (size_t t = 0; t < site.alleles.walks.size(); ++t) {
                out.sequences[t] = site.alleles.step_sequences(graph, t);
            }
            changed = true;
        }
        out.travs[trav] = std::move(walk);
        out.sequences[trav] = std::move(steps);
    }
    return changed;
}

bool AlleleSplicer::chosen_replacement(const StagedSite& child, int strand, bool haploid,
                                       const Traversal& parent_walk,
                                       const unordered_map<nid_t, int>& first_visit,
                                       Replacement& rep) const {
    const LinkageCollector::PhaseCall* phase = phase_table.rendered(child.record_key);
    if (phase == nullptr) {
        return false;
    }
    int child_trav;
    if (phase->ploidy == 1) {
        // The child has one copy, on the parent's strand that carries it.
        if (!haploid && phase->nested_strand != strand) {
            return false;
        }
        child_trav = phase->trav_first;
    } else {
        child_trav = strand == 0 ? phase->trav_first : phase->trav_second;
    }
    if (child_trav < 0 || child_trav >= (int)child.alleles.walks.size()) {
        return false;
    }
    rep.inside_node = sites.is_edit(child.net);
    rep.key = child.record_key;
    if (!spliced_allele(child, child_trav, strand, phase->ploidy == 1, rep.walk, rep.steps)) {
        rep.walk = child.alleles.walks[child_trav];
        if (!child.alleles.spelled_by_walk(child_trav)) {
            rep.steps = child.alleles.step_sequences(graph, child_trav);
        } else {
            rep.steps.clear();
            for (const handle_t& h : rep.walk) {
                rep.steps.push_back(graph.get_sequence(h));
            }
        }
    }
    if (rep.walk.empty()) {
        return false;
    }
    // Find where the parent's walk crosses the child, and which way.
    auto at_start = first_visit.find(graph.get_id(rep.walk.front()));
    auto at_end = first_visit.find(graph.get_id(rep.walk.back()));
    if (at_start == first_visit.end() || at_end == first_visit.end()) {
        return false;
    }
    bool turned = false;
    if (rep.inside_node) {
        // The edit is placed on the spliced walk later. Here find only which way the parent
        // reads it.
        rep.begin = at_start->second;
        rep.end = rep.begin + 1;
        turned = graph.get_is_reverse(parent_walk[rep.begin])
                 != graph.get_is_reverse(rep.walk.front());
    } else if (at_start->second <= at_end->second
               && parent_walk[at_start->second] == rep.walk.front()) {
        rep.begin = at_start->second;
        rep.end = at_end->second + 1;
    } else if (at_end->second <= at_start->second
               && parent_walk[at_end->second] == graph.flip(rep.walk.back())) {
        rep.begin = at_end->second;
        rep.end = at_start->second + 1;
        turned = true;
    } else {
        return false;
    }
    if (turned) {
        // Read the child's allele the way the parent's walk reads it.
        std::ranges::reverse(rep.walk);
        std::ranges::reverse(rep.steps);
        for (size_t i = 0; i < rep.walk.size(); ++i) {
            rep.walk[i] = graph.flip(rep.walk[i]);
            rep.steps[i] = reverse_complement(rep.steps[i]);
        }
    }
    return true;
}

bool AlleleSplicer::splice_routes(const Traversal& parent_walk, const StepSequences* own,
                                  const vector<Replacement>& routes, Traversal& walk,
                                  StepSequences& steps) const {
    bool changed = false;
    size_t r = 0;
    for (int i = 0; i < (int)parent_walk.size();) {
        if (r < routes.size() && routes[r].begin < i) {
            // Skip a route that overlaps one already spliced.
            ++r;
            continue;
        }
        if (r < routes.size() && routes[r].begin == i) {
            // Replace the parent's route through the child with the child's allele.
            const Replacement& rep = routes[r];
            bool same = rep.end - rep.begin == (int)rep.walk.size();
            for (int k = 0; same && k < (int)rep.walk.size(); ++k) {
                same = parent_walk[rep.begin + k] == rep.walk[k]
                       && rep.steps[k] == (own != nullptr ? own->at(rep.begin + k)
                                                          : graph.get_sequence(rep.walk[k]));
            }
            changed = changed || !same;
            walk.insert(walk.end(), rep.walk.begin(), rep.walk.end());
            steps.insert(steps.end(), rep.steps.begin(), rep.steps.end());
            i = rep.end;
            ++r;
            continue;
        }
        walk.push_back(parent_walk[i]);
        steps.push_back(own != nullptr ? own->at(i) : graph.get_sequence(parent_walk[i]));
        ++i;
    }
    return changed;
}

/// Find where `child_walk` lies on `walk`, starting at the first visit to its first node.
///
/// @return The index in `walk` where it starts, or `walk.size()` if it does not lie there.
static size_t find_on_walk(const HandleGraph& graph, const Traversal& walk,
                           const Traversal& child_walk) {
    const nid_t first = graph.get_id(child_walk.front());
    size_t i = 0;
    while (i < walk.size() && graph.get_id(walk[i]) != first) {
        ++i;
    }
    if (i + child_walk.size() <= walk.size()
        && std::equal(child_walk.begin(), child_walk.end(), walk.begin() + i)) {
        return i;
    }
    return walk.size();
}

vector<AlleleSplicer::StepEdit> AlleleSplicer::step_edits_of(const vector<Replacement>& edits,
                                                             const Traversal& walk) const {
    vector<StepEdit> step_edits;
    for (const Replacement& edit : edits) {
        // Find the child's walk on the spliced walk, read either way.
        Traversal turned_walk;
        StepSequences turned_steps;
        for (size_t k = edit.walk.size(); k-- > 0;) {
            turned_walk.push_back(graph.flip(edit.walk[k]));
            turned_steps.push_back(reverse_complement(edit.steps[k]));
        }
        const StepSequences* child_steps = &edit.steps;
        size_t i = find_on_walk(graph, walk, edit.walk);
        if (i == walk.size()) {
            child_steps = &turned_steps;
            i = find_on_walk(graph, walk, turned_walk);
        }
        if (i == walk.size()) {
            continue;
        }
        for (size_t k = 0; k < edit.walk.size(); ++k) {
            const string node_seq = graph.get_sequence(walk[i + k]);
            const string& child_seq = (*child_steps)[k];
            if (child_seq == node_seq) {
                continue;
            }
            // The child's change on this node is what is left without the longest shared prefix
            // and suffix.
            size_t prefix = 0;
            while (prefix < child_seq.size() && prefix < node_seq.size()
                   && child_seq[prefix] == node_seq[prefix]) {
                ++prefix;
            }
            size_t suffix = 0;
            while (suffix < child_seq.size() - prefix && suffix < node_seq.size() - prefix
                   && child_seq[child_seq.size() - 1 - suffix]
                          == node_seq[node_seq.size() - 1 - suffix]) {
                ++suffix;
            }
            step_edits.push_back(StepEdit{i + k, prefix, node_seq.size() - prefix - suffix,
                                          child_seq.substr(prefix,
                                                           child_seq.size() - prefix - suffix),
                                          edit.key});
        }
    }
    return step_edits;
}

bool AlleleSplicer::apply_step_edits(vector<StepEdit> step_edits, StepSequences& steps) {
    std::sort(step_edits.begin(), step_edits.end(),
              [](const StepEdit& a, const StepEdit& b) {
                  // Right to left. At one offset, put a replacement before an insertion in front
                  // of it.
                  return std::tie(a.step, b.at, b.removed, a.key)
                         < std::tie(b.step, a.at, a.removed, b.key);
              });
    bool changed = false;
    for (size_t e = 0; e < step_edits.size(); ++e) {
        const StepEdit& edit = step_edits[e];
        if (e > 0 && step_edits[e - 1].step == edit.step
            && edit.at + edit.removed > step_edits[e - 1].at
            && !(edit.removed == 0 && step_edits[e - 1].removed == 0
                 && edit.at == step_edits[e - 1].at)) {
            continue;
        }
        string& spelled = steps[edit.step];
        if (edit.at + edit.removed > spelled.size()) {
            continue;
        }
        if (spelled.compare(edit.at, edit.removed, edit.added) != 0) {
            spelled.replace(edit.at, edit.removed, edit.added);
            changed = true;
        }
    }
    return changed;
}

bool AlleleSplicer::spliced_allele(const StagedSite& site, int trav, int strand, bool haploid,
                                   Traversal& walk, StepSequences& steps) const {
    auto kids = children.find(site.record_key);
    if (kids == children.end()) {
        return false;
    }
    const Traversal& parent_walk = site.alleles.walks[trav];
    // Index where each node is first visited, which is where a child is entered.
    unordered_map<nid_t, int> first_visit;
    for (int i = (int)parent_walk.size() - 1; i >= 0; --i) {
        first_visit[graph.get_id(parent_walk[i])] = i;
    }
    // Split the replacements into routes and edits inside nodes. Sort each, so that the result
    // does not depend on the order the children were staged in.
    vector<Replacement> routes;
    vector<Replacement> edits;
    for (const StagedSite* child : kids->second) {
        Replacement rep;
        if (chosen_replacement(*child, strand, haploid, parent_walk, first_visit, rep)) {
            (rep.inside_node ? edits : routes).push_back(std::move(rep));
        }
    }
    if (routes.empty() && edits.empty()) {
        return false;
    }
    std::sort(routes.begin(), routes.end(), [](const Replacement& a, const Replacement& b) {
        return std::tie(a.begin, a.end, a.key) < std::tie(b.begin, b.end, b.key);
    });
    std::sort(edits.begin(), edits.end(),
              [](const Replacement& a, const Replacement& b) { return a.key < b.key; });

    walk.clear();
    steps.clear();
    // Splice the routes in walk order, then make the edits of the children inside nodes on the
    // spliced walk.
    const StepSequences own = site.alleles.step_sequences(graph, trav);
    const bool routes_changed =
        splice_routes(parent_walk, own.empty() ? nullptr : &own, routes, walk, steps);
    const bool edits_changed = apply_step_edits(step_edits_of(edits, walk), steps);
    if (!routes_changed && !edits_changed) {
        walk.clear();
        steps.clear();
        return false;
    }
    return true;
}

}
