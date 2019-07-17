#include "snarl_caller.hpp"

namespace vg {

SnarlCaller::SnarlCaller(TraversalGenotyper& traversal_genotyper,
                         SnarlManager& snarl_manager,
                         ostream& out_stream) :
    traversal_genotyper(traversal_genotyper), snarl_manager(snarl_manager), out_stream(out_stream) {
}

SnarlCaller::~SnarlCaller() {
}

void SnarlCaller::call_top_level_snarls(bool recurse_on_fail) {

    header();
    
    // Used to recurse on children of parents that can't be called
    vector<const Snarl*> snarl_queue;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {
        bool was_called = call_snarl(*snarl);
        if (!was_called && recurse_on_fail) {
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
            snarl_queue.insert(snarl_queue.end(), children.begin(), children.end());
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_snarl_parallel(process_snarl);

    // Then recurse on any children the snarl caller failed to handle
    while (!snarl_queue.empty()) {
        vector<const Snarl*> cur_queue;
        std::swap(snarl_queue, cur_queue);
#pragma omp parallel for
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }

}

}

