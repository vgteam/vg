/**
 * \file surjecting_alignment_emitter.cpp
 * Implementation for SurjectingAlignmentEmitter
 */


#include "surjecting_alignment_emitter.hpp"
#include "hts_alignment_emitter.hpp"

#include <map>

namespace vg {

using namespace std;

SurjectingAlignmentEmitter::SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph, unordered_set<path_handle_t> paths,
    unique_ptr<AlignmentEmitter>&& backing, bool prune_suspicious_anchors, bool add_graph_alignment_tag, bool report_supplementary) : surjector(graph), paths(paths), backing(std::move(backing)) {
    
    // Configure the surjector
    surjector.prune_suspicious_anchors = prune_suspicious_anchors;
    surjector.annotate_with_graph_alignment = add_graph_alignment_tag;
    surjector.report_supplementary = report_supplementary;
}

void SurjectingAlignmentEmitter::surject_alignments_in_place(vector<Alignment>& alns) const {
    for (size_t i = 0, n = alns.size(); i < n; ++i) {
        // Surject each alignment and annotate with surjected path position
        auto& aln = alns[i];
        auto surjected = surjector.surject(aln, paths, surject_subpath_global);
        aln = std::move(surjected.front());
        for (size_t j = 1; j < surjected.size(); ++j) {
            alns.emplace_back(std::move(surjected[j]));
        }
    }
}

void SurjectingAlignmentEmitter::surject_paired_alignments_in_place(vector<Alignment>& alns1, vector<Alignment>& alns2,
                                                                    vector<Alignment>& supplementary_alns) const {
    for (size_t i = 0; i < alns1.size(); ++i) {
        auto surjected1 = surjector.surject(alns1[i], paths, surject_subpath_global);
        auto surjected2 = surjector.surject(alns2[i], paths, surject_subpath_global);
        if (surjected1.size() > 1 || surjected2.size() > 1) {

            // find the primaries
            size_t primary_idx1 = -1, primary_idx2 = -1;
            for (size_t j = 0; j < surjected1.size(); ++j) {
                if (!has_annotation(surjected1[j], "supplementary") || !get_annotation<bool>(surjected1[j], "supplementary")) {
                    primary_idx1 = j;
                    break;
                }
            }
            for (size_t j = 0; j < surjected2.size(); ++j) {
                if (!has_annotation(surjected2[j], "supplementary") || !get_annotation<bool>(surjected2[j], "supplementary")) {
                    primary_idx2 = j;
                    break;
                }
            }

            // annotate supplementaries with primary mate info
            const auto& primary_pos1 = surjected1[primary_idx1].refpos(0);
            const auto& primary_pos2 = surjected2[primary_idx2].refpos(0);
            for (size_t j = 0; j < surjected1.size(); ++j) {
                if (j == primary_idx1) {
                    continue;
                }
                supplementary_alns.emplace_back(std::move(surjected1[j]));
                set_annotation(supplementary_alns.back(), "mate_info", 
                               mate_info(primary_pos2.name(), primary_pos2.offset(), primary_pos2.is_reverse(), false));
            }
            for (size_t j = 0; j < surjected2.size(); ++j) {
                if (j == primary_idx1) {
                    continue;
                }
                supplementary_alns.emplace_back(std::move(surjected2[j]));
                set_annotation(supplementary_alns.back(), "mate_info", 
                               mate_info(primary_pos1.name(), primary_pos1.offset(), primary_pos1.is_reverse(), true));
            }

            alns1[i] = std::move(surjected1[primary_idx1]);
            alns2[i] = std::move(surjected2[primary_idx2]);
        }
        else {
            alns1[i] = std::move(surjected1.front());
            alns2[i] = std::move(surjected2.front());
        }
    }
}

void SurjectingAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln_batch_caught(aln_batch);
    // Surject it in place
    surject_alignments_in_place(aln_batch_caught);
    // Forward it along
    backing->emit_singles(std::move(aln_batch_caught));
}

void SurjectingAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns_batch_caught(alns_batch);
    for (auto& mappings : alns_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_singles(std::move(alns_batch_caught));
}

void SurjectingAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln1_batch_caught(aln1_batch);
    vector<Alignment> aln2_batch_caught(aln2_batch);
    // Surject non-supplementary in place and gather supplementary
    vector<Alignment> supplementary_alns;
    surject_paired_alignments_in_place(aln1_batch_caught, aln2_batch_caught, supplementary_alns);
    // Forward them along
    backing->emit_pairs(std::move(aln1_batch_caught), std::move(aln2_batch_caught), std::move(tlen_limit_batch));
    backing->emit_singles(std::move(supplementary_alns));
}

void SurjectingAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns1_batch_caught(alns1_batch);
    vector<vector<Alignment>> alns2_batch_caught(alns2_batch);
    vector<vector<Alignment>> supplementary_batch;
    for (size_t i = 0; i < alns1_batch_caught.size(); ++i) {
        supplementary_batch.emplace_back();
        surject_paired_alignments_in_place(alns1_batch_caught[i], alns2_batch_caught[i], supplementary_batch.back());
    }
    // Forward it along
    backing->emit_mapped_pairs(std::move(alns1_batch_caught), std::move(alns2_batch_caught), std::move(tlen_limit_batch));
    backing->emit_mapped_singles(std::move(supplementary_batch));
}

void SurjectingAlignmentEmitter::emit_extra_message(const std::string& tag, std::string&& data) {
    backing->emit_extra_message(tag, std::move(data));
}

}
