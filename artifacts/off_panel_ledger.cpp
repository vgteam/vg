#include <fstream>

#include "off_panel_ledger.hpp"
#include "../src/utility.hpp"

namespace vg::multipass {

OffPanelLedger::OffPanelLedger(const PathPositionHandleGraph& graph,
                               const ReferencePathTable& ref_paths) :
    graph(graph), ref_paths(ref_paths) {
    // Nothing to do
}

void OffPanelLedger::configure(const SiteReadSource* reads, string dump,
                               const EditCallBar& bar) {
    candidate_source = reads;
    dump_path = std::move(dump);
    call_bar = bar;
}

void OffPanelLedger::open(size_t threads) {
    outcomes_filed = std::make_unique<EditOutcomes>(threads);
}

pair<string, int64_t> OffPanelLedger::reference_position(nid_t node, uint32_t offset) const {
    pair<string, int64_t> found;
    const handle_t handle = graph.get_handle(node);
    graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
        const string path_name = graph.get_path_name(graph.get_path_handle_of_step(step));
        if (!ref_paths.contains(path_name)) {
            return true;
        }
        const size_t len = graph.get_length(handle);
        const bool reverse = graph.get_is_reverse(graph.get_handle_of_step(step));
        const int64_t along = (int64_t)graph.get_position_of_step(step)
                              + (int64_t)(reverse ? len - 1 - offset : offset)
                              + (int64_t)ref_paths.offset_of(path_name);
        string contig = PathMetadata::parse_locus_name(Paths::strip_subrange(path_name));
        found = make_pair(contig == PathMetadata::NO_LOCUS_NAME ? Paths::strip_subrange(path_name)
                                                                : contig,
                          base_path_position(path_name, along));
        return false;
    });
    return found;
}

void OffPanelLedger::report(size_t nodes_after_pass) const {
    size_t top = 0;
    size_t nested = 0;
    size_t staged = 0;
    size_t unplaced = 0;
    size_t indels_staged = 0;
    for (const EditOutcome& o : outcomes_filed->sorted()) {
        if (o.placed == "top") {
            ++top;
        } else if (o.placed == "nested") {
            ++nested;
        } else {
            ++unplaced;
        }
        if (o.staged) {
            ++staged;
            if (o.kind == 'D' || o.kind == 'I') {
                ++indels_staged;
            }
        }
    }
    const vector<EditCandidate>& candidates = candidate_source->counted_edit_candidates();
    cerr << "[vg call] off-panel: " << candidates.size()
         << " candidates counted; the call bar's genotyped as " << top << " top-level and "
         << nested << " nested sites, " << unplaced << " left unplaced; " << staged
         << " staged with their ALT; " << nodes_after_pass
         << " nodes' top-level sites left to after the direct pass" << endl;
    const auto indels = std::ranges::count_if(candidates,
                                              [](const EditCandidate& c) { return c.is_indel(); });
    cerr << "[vg call] off-panel: " << indels << " of the candidates are indels; "
         << indels_staged << " indels staged with their ALT" << endl;
}

/// Write the columns of the dump that describe the candidate itself, from REF to the mean ALT
/// mapping quality.
static void write_candidate_columns(ostream& out, const EditCandidate& c) {
    if (c.is_indel()) {
        // Start an indel's REF and ALT with its anchor base, the base before it.
        out << c.ref << (c.kind == 'D' ? c.sequence : "") << "\t" << c.ref
            << (c.kind == 'I' ? c.sequence : "");
    } else if (c.kind == 'M') {
        out << c.sequence << "\t";
        for (size_t a = 0; a < c.alts.size(); ++a) {
            out << (a > 0 ? "," : "") << c.alts[a];
        }
    } else {
        out << c.ref << "\t" << c.alt;
    }
    out << "\t" << c.alt_fragments << "\t" << c.ref_fragments << "\t" << c.other_fragments << "\t"
        << c.discordant_fragments << "\t" << c.fraction() << "\t" << c.mean_alt_mapq;
}

/// Write the columns of the dump that say what became of a candidate, from `placed` to `dr`.
///
/// @param outcome The outcome of the candidate, or null if the call bar did not admit it.
static void write_outcome_columns(ostream& out, const EditCandidate& c, const EditOutcome* outcome,
                                  const EditCallBar& bar) {
    if (outcome != nullptr) {
        out << "\t" << outcome->placed << "\t" << (outcome->site.empty() ? "." : outcome->site)
            << "\t";
        for (size_t i = 0; i < outcome->genotype.size(); ++i) {
            out << (i ? "/" : "") << outcome->genotype[i];
        }
        if (outcome->genotype.empty()) {
            out << ".";
        }
        out << "\t" << outcome->gq << "\t" << (outcome->staged ? 1 : 0);
    } else {
        out << "\t" << (bar.admits(c) ? "unplaced" : "below_bar") << "\t.\t.\t.\t0";
    }
    out << "\t" << c.homopolymer << "\t";
    if (outcome != nullptr && !outcome->genotype.empty()) {
        out << outcome->depth_ratio;
    } else {
        out << ".";
    }
}

void OffPanelLedger::write(size_t nodes_after_pass) const {
    if (!is_enabled()) {
        return;
    }
    report(nodes_after_pass);
    if (dump_path.empty()) {
        return;
    }
    ofstream out(dump_path);
    if (!out) {
        cerr << "error [vg call]: could not open " << dump_path << " for the off-panel dump" << endl;
        return;
    }
    out << "#contig\tpos\tnode\toffset\tref\talt\talt_frags\tref_frags\tother_frags"
           "\tdiscordant_frags\taf\tmean_alt_mapq\tplaced\tsite\tgt\tgq\tstaged\thomopolymer\tdr\n";
    const vector<EditOutcome> outcomes =
        outcomes_filed != nullptr ? outcomes_filed->sorted() : vector<EditOutcome>();
    // Walk the outcomes alongside the candidates, which are in the same order.
    auto outcome = outcomes.begin();
    for (const EditCandidate& c : candidate_source->counted_edit_candidates()) {
        const auto [contig, position] = reference_position(c.node, c.offset);
        out << (contig.empty() ? "." : contig) << "\t";
        if (contig.empty()) {
            out << ".";
        } else {
            out << position;
        }
        out << "\t" << c.node << "\t" << c.offset << "\t";
        write_candidate_columns(out, c);
        while (outcome != outcomes.end()
               && tie(outcome->node, outcome->offset, outcome->kind, outcome->alt,
                      outcome->sequence)
                      < tie(c.node, c.offset, c.kind, c.alt, c.sequence)) {
            ++outcome;
        }
        const bool has_outcome = outcome != outcomes.end() && outcome->node == c.node
                                 && outcome->offset == c.offset && outcome->kind == c.kind
                                 && outcome->alt == c.alt && outcome->sequence == c.sequence;
        write_outcome_columns(out, c, has_outcome ? std::to_address(outcome) : nullptr, call_bar);
        out << "\n";
    }
}

}
