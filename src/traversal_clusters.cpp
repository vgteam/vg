#include "traversal_clusters.hpp"

namespace vg {

vector<vector<int>> cluster_traversals(const PathHandleGraph* graph,
                                       const vector<Traversal>& traversals,
                                       const vector<int64_t>& traversal_order,
                                       double min_jaccard) {

    assert(traversal_order.size() == traversals.size());
    
    // the values are indexes in the input traversals vector
    // the "reference" traversal of each cluster (to which distance is computed)
    // is always its first element
    vector<vector<int>> clusters;

    for (const int64_t& i : traversal_order) {
        const Traversal& trav = traversals[i];
        double min_jaccard = numeric_limits<double>::max();
        int64_t min_cluster_idx = -1;
        for (int64_t j = 0; j < clusters.size(); ++j) {
            const Traversal& cluster_trav = traversals[clusters[j][0]];
            double jac = jaccard_coefficient(trav, cluster_trav);
            if (jac < min_jaccard) {
                min_jaccard = jac;
                min_cluster_idx = j;
                if (jac == 0) {
                    break;
                }                    
            }        
        }
        if (min_cluster_idx >= 0) {
            // we've found a suitably similar cluster, add it do that
            clusters[min_cluster_idx].push_back(i);
        } else {
            // there's no cluster close enough, need to start a new one
            clusters.push_back({i});
        }        
    }

    return clusters;
}


}
