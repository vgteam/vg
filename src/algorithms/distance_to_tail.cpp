// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#include "distance_to_tail.hpp"
#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

bool is_tail_node(handle_t h, const HandleGraph* g) {
	bool no_right_edges = true;
	h = g->forward(h);
    g->follow_edges(h, false, [&](const handle_t& ignored) {
        // We found a right edge!
        no_right_edges = false;
        // We only need one
        return false;
    });
    return no_right_edges;
}

int32_t distance_to_tail(handle_t h, int32_t limit, const HandleGraph* graph) {
	unordered_set<handle_t> seen;
	return distance_to_tail(h, limit, 0, seen, graph);
}

int32_t distance_to_tail(handle_t h, int32_t limit, int32_t dist, unordered_set<handle_t>& seen, const HandleGraph* graph) {
	if (seen.count(h)) return -1;
	seen.insert(h);
	if (limit <= 0) {
		return -1;
	}

	if(is_tail_node(h, graph)) {
		return dist;
	}

	int32_t t = -1; 

	graph->follow_edges(h, false, [&](const handle_t& current) {
		int32_t l = graph->get_length(current);
		t = distance_to_tail(current, limit-l, dist+l, seen, graph);
		if (t != -1) {
			return false;
		}
		else {
			return true;
		}
	});
	return t;
}
}
}
