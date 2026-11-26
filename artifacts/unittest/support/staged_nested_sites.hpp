#ifndef VG_UNITTEST_STAGED_NESTED_SITES_HPP_INCLUDED
#define VG_UNITTEST_STAGED_NESTED_SITES_HPP_INCLUDED

/** \file
 * The sites of the nested fixture graph, staged as the render has them, with a phase table and a
 * splicer over them, for tests of code that reads the called strands of a site.
 */

#include <vector>

#include "../../../src/integrated_snarl_finder.hpp"
#include "../../allele_splicer.hpp"
#include "../../edit_snarl_decomposition_overlay.hpp"
#include "../../phase_table.hpp"
#include "../../staged_site.hpp"
#include "decomposition_pair.hpp"
#include "nested_graph.hpp"

namespace vg {
namespace unittest {

using namespace std;

/**
 * The three nested sites of `NestedGraph`, and a site inside node 4, staged under one another.
 *
 * The top site, 1 to 8, has the reference walk through every node but 6, the walk from 1 straight
 * to 8, and the walk through node 6. The middle site, 2 to 7, has the reference walk and the walk
 * through node 6. The inner site, 3 to 5, has the reference walk and the walk that leaves out
 * node 4. The site inside node 4 has the node and the node with its second base changed to A.
 *
 * Give the sites their phases with `phase`, then call `index`.
 */
struct StagedNestedSites {
    NestedGraph nested;
    IntegratedSnarlFinder finder;
    DecompositionPair pair;
    EditSnarlDecompositionOverlay overlay;
    multipass::PhaseTable phases;
    multipass::StagedSiteTable table;
    multipass::AlleleSplicer splicer;

    const size_t top_key = 1;
    const size_t middle_key = 2;
    const size_t inner_key = 3;
    const size_t edit_key = 4;

    StagedNestedSites() :
        finder(nested.graph), pair(nested.graph, finder), overlay(pair.adapter, nested.graph),
        splicer(nested.graph, overlay, phases) {
        table.start(1);
        table.add_top_level(snarl(top_key, 0, 1, 8,
                                  {{1, 2, 3, 4, 5, 7, 8}, {1, 8}, {1, 2, 6, 7, 8}}));
        table.add_nested(snarl(middle_key, top_key, 2, 7, {{2, 3, 4, 5, 7}, {2, 6, 7}}));
        table.add_nested(snarl(inner_key, middle_key, 3, 5, {{3, 4, 5}, {3, 5}}));
        multipass::StagedSite edit;
        edit.id = "edit";
        edit.record_key = edit_key;
        edit.placement.parent_record_key = inner_key;
        edit.net = overlay.add_edit(nested.handle(4), 1, nested.handle(4), 2);
        edit.alleles.walks = {{nested.handle(4)}, {nested.handle(4)}};
        edit.alleles.ref_allele = 0;
        edit.alleles.edits = {{}, {AlleleEdit{0, 1, 1, "A"}}};
        table.add_nested(std::move(edit));
    }

    /// Give the site with record key `key` a phase with traversal `first` on strand 0 and
    /// `second` on strand 1.
    ///
    /// At `ploidy` 1 the site has `first` alone, on strand `nested_strand` of its parent.
    void phase(size_t key, int first, int second, size_t ploidy = 2, int nested_strand = -1) {
        multipass::PhaseTable::PhaseCall call;
        call.record_key = key;
        call.contig = "ref";
        call.trav_first = first;
        call.trav_second = second;
        call.ploidy = ploidy;
        call.nested_strand = (int8_t)nested_strand;
        phases.calls().push_back(call);
    }

    /// Hand the nested sites to the render queue and index them, as the render does before it
    /// writes any record.
    void index() {
        phases.freeze_for_render(true);
        table.gather_nested();
        table.hand_off();
        splicer.index(table, true);
    }

    /// Get the staged site with record key `key`, once `index` has run.
    multipass::StagedSite& site(size_t key) {
        for (multipass::StagedSite& staged : table.queue(0)) {
            if (staged.record_key == key) {
                return staged;
            }
        }
        throw std::out_of_range("no staged site has the record key");
    }

private:
    /// Make the staged site of the snarl from node `start` to node `end`, with the walks that
    /// visit `walks` forward, the first being the reference.
    multipass::StagedSite snarl(size_t key, size_t parent_key, nid_t start, nid_t end,
                                const vector<vector<nid_t>>& walks) {
        multipass::StagedSite site;
        site.id = ">" + std::to_string(start) + ">" + std::to_string(end);
        site.record_key = key;
        site.placement.parent_record_key = parent_key;
        site.net = find_snarl(overlay, nested.graph, start, end);
        for (const vector<nid_t>& walk : walks) {
            Traversal visits;
            for (nid_t id : walk) {
                visits.push_back(nested.handle(id));
            }
            site.alleles.walks.push_back(std::move(visits));
        }
        site.alleles.ref_allele = 0;
        return site;
    }
};

}
}

#endif
