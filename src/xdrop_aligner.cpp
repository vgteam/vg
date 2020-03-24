
/**
 * @file xdrop_aligner.hpp
 * @author Hajime Suzuki
 * @date 2018/03/23
 */
#include <cstdio>
#include <assert.h>
#include <utility>
#include "mem.hpp"
#include "xdrop_aligner.hpp"
#include "proto_handle_graph.hpp"
#include "reverse_graph.hpp"
#include "extra_node_graph.hpp"
#include "algorithms/topological_sort.hpp"
#include "convert_handle.hpp"

#define DZ_FULL_LENGTH_BONUS

// We require these particular values for this enum because we index arrays with it.
enum { MISMATCH = 1, MATCH = 2, INS = 3, DEL = 4 };

// Tell Dozeu to use those values
#define DZ_CIGAR_OP 0x04030201

#include <dozeu/dozeu.h>

// To turn on debugging:
// #define DEBUG

#ifdef DEBUG
#include <dozeu/log.h>
#endif

// #include "json2pb.h"			// I don't know whether it's needed

using namespace vg;

XdropAligner::XdropAligner(XdropAligner const &rhs)
{
	dz = dz_init(
		rhs.dz->matrix,
		(uint16_t)rhs.dz->giv[0],
		(uint16_t)rhs.dz->gev[0],
		(uint16_t)rhs.dz->max_gap_len,
		(uint16_t)rhs.dz->bonus
	);
    
    aa_match = rhs.aa_match;
}

XdropAligner& XdropAligner::operator=(XdropAligner const &rhs)
{
	if(this == &rhs) { return(*this); }
	dz = dz_init(
		rhs.dz->matrix,
		(uint16_t)rhs.dz->giv[0],
		(uint16_t)rhs.dz->gev[0],
		(uint16_t)rhs.dz->max_gap_len,
		(uint16_t)rhs.dz->bonus
	);
    aa_match = rhs.aa_match;
	return(*this);
}

XdropAligner::XdropAligner(XdropAligner&& rhs)
{
	dz_destroy(dz);
	dz = rhs.dz;		// move
	rhs.dz = nullptr;
    aa_match = rhs.aa_match;
}

XdropAligner& XdropAligner::operator=(XdropAligner&& rhs)
{
	if(this == &rhs) { return(*this); }
	dz_destroy(dz);
	dz = rhs.dz;
	rhs.dz = nullptr;
    aa_match = rhs.aa_match;
	return(*this);
}

XdropAligner::XdropAligner()
{
	dz = nullptr;
    aa_match = 0;
}

XdropAligner::XdropAligner(
	int8_t const *_score_matrix,
	int8_t _gap_open,
	int8_t _gap_extension,
	int32_t _full_length_bonus,
	uint32_t _max_gap_length)
{
	assert(_gap_open - _gap_extension >= 0);
	assert(_gap_extension > 0);
	assert(_full_length_bonus >= 0);
	assert(_max_gap_length > 0);

	uint64_t go = _gap_open - _gap_extension, ge = _gap_extension;
	dz = dz_init((int8_t const *)_score_matrix, go, ge, _max_gap_length, _full_length_bonus);
    aa_match = _score_matrix[0];
	// bench_init(bench);
}

XdropAligner::XdropAligner(
	int8_t _match,
	int8_t _mismatch,
	int8_t _gap_open,
	int8_t _gap_extension,
	int32_t _full_length_bonus,
	uint32_t _max_gap_length)
{
	assert(_match > 0);
	assert(_mismatch > 0);
	assert(_gap_open - _gap_extension >= 0);
	assert(_gap_extension > 0);
	assert(_full_length_bonus >= 0);
	assert(_max_gap_length > 0);

	int8_t const M = _match, X = -((int8_t)_mismatch);
	int8_t _score_matrix[16] = {
		M, X, X, X,
		X, M, X, X,
		X, X, M, X,
		X, X, X, M
	};
	uint64_t go = _gap_open - _gap_extension, ge = _gap_extension;
	dz = dz_init((int8_t const *)_score_matrix, go, ge, _max_gap_length, _full_length_bonus);
    aa_match = M;
	// bench_init(bench);
}

XdropAligner::~XdropAligner(void)
{
	dz_destroy(dz);
	dz = nullptr;
	// fprintf(stderr, "xdrop: time(%lu), count(%lu)\n", bench_get(bench) / 1000, bench_get_count(bench));
}

static inline char comp(char x)
{
	switch(x) {
		case 'a': case 'A': return('T');
		case 'c': case 'C': return('G');
		case 'g': case 'G': return('C');
		case 't': case 'T': return('A');
		default: return('N');
	}
}

/* edge packing / extraction macros */
#define _edge(_s, _d)		( ((int64_t)(_d)<<32) | ((int64_t)(_s)) )
#define _src_index(_x)		( (_x) & 0xffffffff )
#define _dst_index(_x)		( (_x)>>32 )

void XdropAligner::build_id_index_table(OrderedGraph const &graph)
{
	// construct node_id -> index map
	id_to_index.clear();							// unordered_map< id_t, uint64_t >
	for(size_t i = 0; i < graph.order.size(); i++) {
		id_to_index[graph.graph.get_id(graph.order[i])] = (uint32_t)i;
		debug("i(%lu), id(%ld), length(%lu)", i, graph.graph.get_id(graph.order[i]), graph.graph.get_sequence(graph.order[i]).length());
	}
	return;
}

void XdropAligner::build_index_edge_table(OrderedGraph const &graph, uint32_t const seed_node_index, bool left_to_right)
{
	// build (src_index, dst_index) array
	index_edges.clear();
	index_edges_head.clear();
    
    // TODO: reserve space for all the edges.
    // We no longer have an easy way to count edges in the handle graph interface...

    graph.graph.for_each_edge([&](handlegraph::edge_t handle_pair) {
        // For each edge
        
        
        if (graph.graph.get_is_reverse(handle_pair.first)) {
            // Flip reverse-reverse to forward-forward orientation
            std::swap<handle_t>(handle_pair.first, handle_pair.second);
            handle_pair.first = graph.graph.flip(handle_pair.first);
            handle_pair.second = graph.graph.flip(handle_pair.second);
        }
        
        // Make sure it isn't an edge we can't deal with (must be forward to forward).
        assert(!graph.graph.get_is_reverse(handle_pair.first));
        assert(!graph.graph.get_is_reverse(handle_pair.second));
        
        auto from_id = graph.graph.get_id(handle_pair.first);
        auto to_id = graph.graph.get_id(handle_pair.second);
        
        // Record the edge in the normal index.
        // If left_to_right is true, use function 0, which puts to in the high bits
        // If it is false, put from in the high bits
		
        index_edges.push_back(edge[!left_to_right](id_to_index[from_id],
            id_to_index[to_id]));
            
        debug("append edge, %u (id %lu) -> %u (id %lu)", id_to_index[from_id], from_id, id_to_index[to_id], to_id);
        
        // if left_to_right is true, use function 0, so compare with <
        // Then invert that, so we check (index of from_id) > seed_node_index
        // If left_to_right is false, this is (index of to_id) < seed_node_index
		if(!compare[!left_to_right](id_to_index[left_to_right ? from_id : to_id], seed_node_index)) {
            // We do not need to record this edge in the flipped index.
            return true;
        }
        
        // If left_to_right is true, use function 1, so put from in the high bits
        // If left_to_right is false, put to in the high bits.
		index_edges_head.push_back(edge[left_to_right](id_to_index[from_id], id_to_index[to_id]));	// reversed
		// fprintf(stderr, "append head edge, %u -> %u\n", id_to_index[from_id], id_to_index[to_id]);
        
        return true;
    });

    // Sort index_edges in ascending order if left_to_right is true, and
    // descending order otherwise.
    //
    // If left_to_right is true, we put the edge's to index in the high bits,
    // so this will order edges ascending by to node and then from node.
    //
    // If left_to_right is false, we put the edge's from index in the high
    // bits, so this will order edges descending by from node and then to node.
    sort(index_edges.begin(), index_edges.end(), compare[!left_to_right]);
    
    // Sort index_edges_head in descending order if left_to_right is true, and
    // ascending order otherwise.
    //
    // If left_to_right is true, we put the from index in the high bits, so
    // this will order edges in descending order by from and then to.
    //
    // If left_to_right is false, we put the to index in the high bits, so this
    // will order edges in asecnding order by to and then from.
	sort(index_edges_head.begin(), index_edges_head.end(), compare[left_to_right]);
	return;
}

XdropAligner::graph_pos_s XdropAligner::calculate_seed_position(OrderedGraph const &graph, vector<MaximalExactMatch> const &mems, size_t query_length, bool direction)
{
	/*
	 * seed selection:
	 * the most upstream for forward one, the most downstream one for reverse one.
	 * FIXME: better selection strategy (or multiple trial) would be possible.
	 *
	 * node_id extraction:
	 * use the most downstrem one for forward mapping to avoid any *artifact*
	 * caused by seed position. (redundant multiple mapping may reported due to redundant seeds;
	 * having almost the same position, length, and sequence but pass slightly different paths)
	 * GSSW aligner avoids such artifacts because it does not depends on any seed positions on
	 * computing the DP matrix (it is true local alignment).
	 */
    graph_pos_s pos;
	MaximalExactMatch const &seed = mems[direction ? (mems.size() - 1) : 0];		// here mems is never empty

	auto seed_pos = direction ? seed.nodes.front() : seed.nodes.back();
	size_t node_id = gcsa::Node::id(seed_pos);
	size_t node_offset = gcsa::Node::offset(seed_pos);
	pos.node_index = id_to_index[node_id];

	// calc ref_offset
	handle_t n = graph.order[pos.node_index];
	pos.ref_offset = direction ? (graph.graph.get_length(n) - node_offset) : node_offset;

	// calc query_offset (FIXME: is there O(1) solution?)
	pos.query_offset = query_length - (seed.end - seed.begin);
	for(std::string::const_iterator p = seed.end; *p != '\0'; pos.query_offset--, p++) {}
	pos.query_offset = direction ? query_length - pos.query_offset : pos.query_offset;
	// fprintf(stderr, "calc_seed_pos, direction(%d), rpos(%lu, %u), qpos(%u), len(%lu)\n", direction, pos.node_index, pos.ref_offset, pos.query_offset, seed.end - seed.begin);
	return(pos);
}

XdropAligner::graph_pos_s XdropAligner::calculate_max_position(OrderedGraph const &graph, graph_pos_s const &seed_pos, size_t max_node_index, bool direction)
{
	// save node id
	graph_pos_s pos;
	pos.node_index = max_node_index;
    
    // Find the node
    handle_t n = graph.order[pos.node_index];

    assert(forefronts[max_node_index]->mcap != nullptr);
    if (forefronts[max_node_index]->mcap == nullptr) {
        // No alignment. Not safe to call dz_calc_max_pos.
        // Bail out early.
        pos.ref_offset = direction ? 0 : graph.graph.get_length(n);
        pos.query_offset = seed_pos.query_offset;
        return(pos);
    }

	// calc max position on the node
	uint64_t max_pos = (uint64_t)dz_calc_max_pos(dz, forefronts[max_node_index]);

	// ref-side offset fixup
	int32_t rpos = (int32_t)(max_pos>>32);

	
	pos.ref_offset = direction ? -rpos : (graph.graph.get_length(n) - rpos);

	// query-side offset fixup
	int32_t qpos = max_pos & 0xffffffff;
	pos.query_offset = seed_pos.query_offset + (direction ? -qpos : qpos);
	// fprintf(stderr, "calc_max_pos, rpos(%lu, %u), qpos(%u, %u)\n", pos.node_index, pos.ref_offset, seed_pos.query_offset, pos.query_offset);
	return(pos);
}

XdropAligner::graph_pos_s XdropAligner::scan_seed_position(OrderedGraph const &graph, std::string const &query_seq, bool direction)
{
	uint64_t const qlen = query_seq.length(), scan_len = qlen < 15 ? qlen : 15;		// FIXME: scan_len should be variable

	debug("scan seeds, direction(%u), qlen(%lu), scan_len(%lu)", direction, qlen, scan_len);
	for(vector<uint64_t>::const_reverse_iterator p = index_edges.rbegin(); p != index_edges.rend(); p++) {
		// debug("i(%lu), %lu -> %lu", p - index_edges.rbegin(), _dst_index(*p), _src_index(*p));
	}
	dz_query_s const *packed_query = (direction
		? dz_pack_query_reverse(dz, &query_seq.c_str()[0],               scan_len)
		: dz_pack_query_forward(dz, &query_seq.c_str()[qlen - scan_len], scan_len)
	);

	int64_t inc = direction ? -1 : 1;
	size_t max_node_index = direction ? graph.order.size() - 1 : 0, prev_node_index = max_node_index - inc;
	for(vector<uint64_t>::const_reverse_iterator p = index_edges.rbegin(); p != index_edges.rend();) {
		size_t node_index = _src_index(*p), n_incoming_edges = 0;
		debug("i(%lu), node_index(%lu), max_node_index(%lu), prev_node_index(%lu)", p - index_edges.rbegin(), node_index, max_node_index, prev_node_index);

		// fill ends
		for(size_t empty_node_index = prev_node_index + inc; compare[direction](empty_node_index, node_index); empty_node_index += inc) {
			debug("fill root for empty_node_index(%lu)", empty_node_index);
			auto seq = graph.graph.get_sequence(graph.order[empty_node_index]);
			forefronts[empty_node_index] = dz_scan(dz,
				packed_query,
				dz_root(dz), 1,
				&seq.c_str()[direction ? seq.length() : 0],
				direction ? -seq.length() : seq.length(),
				empty_node_index
			);
		}
		prev_node_index = node_index;


		// count #incoming edges
		while(p + n_incoming_edges != index_edges.rend() && _src_index(p[n_incoming_edges]) == node_index) { n_incoming_edges++; }
		assert(n_incoming_edges > 0);
		debug("n_incoming_edges(%lu)", n_incoming_edges);

		// pack incoming edges
		dz_forefront_s const *incoming_forefronts[n_incoming_edges + 1];
		for(size_t i = 0; i < n_incoming_edges; i++) {
			dz_forefront_s const *t = forefronts[_dst_index(p[i])];
			assert(t != nullptr);
			incoming_forefronts[i] = t;
		}

		// forward edge_index_base
		p += n_incoming_edges;

		handle_t n = graph.order[node_index];
		auto ref_seq = graph.graph.get_sequence(n);
		forefronts[node_index] = dz_scan(dz,
			packed_query,
			incoming_forefronts, n_incoming_edges,
			&ref_seq.c_str()[direction ? ref_seq.length() : 0],
			direction ? -ref_seq.length() : ref_seq.length(),
			node_index
		);
		debug("i(%lu), forefront(%p), max(%d)", p - index_edges.rbegin(), forefronts[node_index], forefronts[node_index]->max);

		if(forefronts[node_index]->max + (direction & dz_geq(forefronts[node_index])) > forefronts[max_node_index]->max) {
			max_node_index = node_index;
		}

	}

	graph_pos_s pos;
	pos.node_index = 0;
	pos.ref_offset = 0;
	pos.query_offset = direction ? scan_len : qlen - scan_len;
	graph_pos_s p = calculate_max_position(graph, pos, max_node_index, direction);
	debug("node_index(%lu), ref_offset(%d), query_offset(%d), max(%d)", p.node_index, p.ref_offset, p.query_offset, forefronts[max_node_index]->max);
	return(p);
}

size_t XdropAligner::extend(
	OrderedGraph const &graph,
	vector<uint64_t>::const_iterator begin,			// must be sorted ascending order for forward extension, descending order in reverse
	vector<uint64_t>::const_iterator end,
	dz_query_s const *packed_query,
	size_t seed_node_index,
	uint64_t seed_offset,							// offset: 0-------->L for both forward and reverse
	bool right_to_left)									// true for a right-to-left pass with left-to-right traceback, false otherwise
{
	// get root node
	auto root_seq = graph.graph.get_sequence(graph.order[seed_node_index]);

	// load position and length
	uint64_t rpos = seed_offset;
	int64_t rlen = (right_to_left ? 0 : root_seq.length()) - seed_offset;
    
    debug("extend pass: %s over %lu forefronts from %lu", right_to_left ? "right-to-left" : "left-to-right", forefronts.size(), seed_node_index);
    
	debug("rpos(%lu), rlen(%ld)", rpos, rlen);
	forefronts[seed_node_index] = dz_extend(dz,
		packed_query,
		dz_root(dz), 1,
		&root_seq.c_str()[rpos], rlen, seed_node_index
	);

	int64_t inc = right_to_left ? -1 : 1;
	size_t max_node_index = seed_node_index;
    size_t prev_node_index = seed_node_index;
	debug("root: node_index(%lu, %ld), ptr(%p), score(%d)", seed_node_index, graph.graph.get_id(graph.order[seed_node_index]), forefronts[seed_node_index], forefronts[seed_node_index]->max);

    debug("edge count: %ld", end - begin);

	for(vector<uint64_t>::const_iterator p = begin; p != end;) {
		size_t node_index = _dst_index(*p);
        size_t n_incoming_edges = 0;
        size_t n_incoming_forefronts = 0;

		for(size_t empty_node_index = prev_node_index + inc; compare[right_to_left](empty_node_index, node_index); empty_node_index += inc) {
			forefronts[empty_node_index] = nullptr;
		}
		prev_node_index = node_index;

		// pack incoming edges
		while(p + n_incoming_edges != end && _dst_index(p[n_incoming_edges]) == node_index) { n_incoming_edges++; }
		assert(n_incoming_edges > 0);

		// pack incoming edges
		dz_forefront_s const *incoming_forefronts[n_incoming_edges + 1];
		for(size_t k = 0; k < n_incoming_edges; k++) {
			int64_t s = _src_index(p[k]), d = _dst_index(p[k]), r = seed_node_index;
            // fiddly way to check that s and d are on the same side of r
			if((s - r) * (d - r) < 0) { continue; }

			dz_forefront_s const *t = forefronts[_src_index(p[k])];
			if(t == nullptr || dz_is_terminated(t)) { continue; }	// skip terminated
			incoming_forefronts[n_incoming_forefronts++] = t;		// copy onto stack
		}
		debug("node_index(%lu, %ld), n_incoming_edges(%lu, %lu)", node_index, graph.graph.get_id(graph.order[node_index]), n_incoming_edges, n_incoming_forefronts);

		// forward edge_index_base
		p += n_incoming_edges;
		if(n_incoming_forefronts == 0) { forefronts[node_index] = nullptr; continue; }

		handle_t n = graph.order[node_index];
		auto ref_seq = graph.graph.get_sequence(n);
		forefronts[node_index] = dz_extend(dz,
			packed_query,
			incoming_forefronts, n_incoming_forefronts,
			&ref_seq.c_str()[right_to_left ? ref_seq.length() : 0],
			right_to_left ? -ref_seq.length() : ref_seq.length(),
			node_index
		);
		if(forefronts[node_index]->max + (right_to_left & dz_geq(forefronts[node_index])) > forefronts[max_node_index]->max) {
			max_node_index = node_index;
		}
		debug("node_index(%lu, %ld), n_incoming_edges(%lu, %lu), forefront(%p), range[%u, %u), term(%u), max(%d), curr(%d, %d, %u)",
			node_index, graph.graph.get_id(n), n_incoming_edges, n_incoming_forefronts,
			forefronts[node_index], forefronts[node_index]->r.spos, forefronts[node_index]->r.epos, dz_is_terminated(forefronts[node_index]),
			forefronts[max_node_index]->max, forefronts[node_index]->max, forefronts[node_index]->inc,
            right_to_left & dz_geq(forefronts[node_index]));
	}
    
    // Get max query pos
    assert(max_node_index <= forefronts.size());
    assert(forefronts[max_node_index] != nullptr);
    
    if (forefronts[max_node_index]->mcap == nullptr) {
        // We got no optimal alignment. This must be non-null to dz_calc_max_qpos.
        debug("null mcap in winning forefront; no alignment");
        return(max_node_index);
    } else {
    
        uint64_t query_max_pos = dz_calc_max_qpos(dz, forefronts[max_node_index]);
        uint64_t ref_node_max_pos = dz_calc_max_rpos(dz, forefronts[max_node_index]);
        
        debug("max(%p), score(%d), qpos(%lu), rpos(%lu)", forefronts[max_node_index], forefronts[max_node_index]->max, query_max_pos, ref_node_max_pos);
        return(max_node_index);
    }
}

// append an edit at the end of the current mapping array, returns forwarded length on the query
size_t XdropAligner::push_edit(
	Mapping *mapping,
	uint8_t op,
	char const *alt,
	size_t len)
{
	/* see aligner.cpp:gssw_mapping_to_alignment */
	#define _add_edit(_from_len, _to_len, _subseq) { \
		Edit *e = mapping->add_edit(); \
		e->set_from_length((_from_len)); \
		e->set_to_length((_to_len)); \
		/* expect a branch dependent on a compile-time NULL will be eliminated */ \
		if((_subseq) != nullptr) { e->set_sequence((char const *)(_subseq), (size_t)(_to_len)); } \
	}

	debug("push_edit: %lu%c in %s", len, "-XMID"[op], pb2json(mapping->position()).c_str());
	if(op == MISMATCH) {
		// break down into multiple SNVs
		for(size_t i = 0; i < len; i++) { _add_edit(1, 1, &alt[i]); }
	} else if (len > 0) {
        // Only add an edit if the operation has nonzero length
		alt = (op == INS) ? alt : nullptr;
		size_t rlen = (op & 0x01) ? 0 : len;
		size_t qlen = (op & 0x02) ? len : 0;
		_add_edit(rlen, qlen, alt); len = qlen;
	}
	return(len);

	#undef _add_edit
}

void XdropAligner::calculate_and_save_alignment(
	Alignment &alignment,
	OrderedGraph const &graph,
	graph_pos_s const &head_pos,
	size_t tail_node_index,
	bool left_to_right)
{
	// clear existing alignment (no matter if any significant path is not obtained)
	alignment.clear_path();
	alignment.set_score(forefronts[tail_node_index]->max);
	if(forefronts[tail_node_index]->max == 0) { 
        // No alignment scoring anything other than 0, or not safe to try and trace back.
        // Emit a full-length insertion
        debug("no alignment; emit full length insertion");
        Mapping* m = alignment.mutable_path()->add_mapping();
        handle_t start = graph.order.front();
        m->mutable_position()->set_node_id(graph.graph.get_id(start));
        m->mutable_position()->set_is_reverse(graph.graph.get_is_reverse(start));
        m->mutable_position()->set_offset(0);
        m->set_rank(1);
        Edit* e = m->add_edit();
        e->set_from_length(0);
        e->set_to_length(alignment.sequence().size());
        e->set_sequence(alignment.sequence());
        return;
    }
    
    // If we have a nonzero score we should have a nonempty mcap on the winning forefront
    assert(forefronts[tail_node_index]->mcap != nullptr);
    
	// traceback. This produces an alignment in the same order as we traversed the nodes when filling them in.
	dz_alignment_s const *aln = dz_trace(dz, forefronts[tail_node_index]);
	if(aln == nullptr || aln->path_length == 0) {
        // No alignment actually computed.
        // Emit a full-length insertion
        debug("no traceback; emit full length insertion");
        Mapping* m = alignment.mutable_path()->add_mapping();
        handle_t start = graph.order.front();
        m->mutable_position()->set_node_id(graph.graph.get_id(start));
        m->mutable_position()->set_is_reverse(graph.graph.get_is_reverse(start));
        m->mutable_position()->set_offset(0);
        m->set_rank(1);
        Edit* e = m->add_edit();
        e->set_from_length(0);
        e->set_to_length(alignment.sequence().size());
        e->set_sequence(alignment.sequence());
        return;
    }
    
    #ifdef DEBUG
    // Dump the Dozeu alignment
    // Make sure to translate CIGAR numbers to characters
    string translated_path; 
    for (auto* op = aln->path; *op != 0; ++op) {
        translated_path.push_back("-XMID"[*op]);
    }
    debug("ref_length(%u), query_length(%u), span_length(%d), path_length(%d), score(%d), path(%s)", aln->ref_length, aln->query_length, aln->span_length, aln->path_length, aln->score, translated_path.c_str());
    debug("matches(%u), mismatches(%u), inserts(%u), deletes(%u)", aln->match_count, aln->mismatch_count, aln->ins_count, aln->del_count);
    for(size_t i = 0; i < aln->span_length; i++) {
        dz_path_span_s const *s = &aln->span[i];
        
        translated_path.clear();
        for (auto* op = &aln->path[s->offset]; op != &aln->path[s->offset] + (s[1].offset - s[0].offset); ++op) {
            translated_path.push_back("-XMID"[*op]);
        }
        
        debug("node_id(%u), subpath_length(%u:%u-%u), subpath(%s)",
            s->id,
            s[1].offset - s[0].offset,
            s[0].offset, s[1].offset,
            translated_path.c_str()
        );
    }
    #endif
    
    // Make sure it ends where we started traceback from.
	assert(aln->span[aln->span_length - 1].id == tail_node_index);
    
    // Check the length and make sure it is right.
    if (aln->query_length > alignment.sequence().size()) {
        cerr << "[vg xdrop_aligner.cpp] Error: dozeu alignment query_length longer than sequence" << endl;
        exit(1);
    }
    // aln->query_length can be shorter than alignment.sequence().size() if we
    // didn't traceback from the very last base of the query, or if we didn't
    // pack the whole query because of an offset.

	#define _push_mapping(_id) ({ \
		handle_t n = graph.order[(_id)]; \
		Mapping *mapping = path->add_mapping(); \
		mapping->set_rank(path->mapping_size()); \
		Position *position = mapping->mutable_position(); \
		position->set_node_id(graph.graph.get_id(n)); \
        position->set_is_reverse(graph.graph.get_is_reverse(n)); \
		position->set_offset(ref_offset); ref_offset = 0; \
		mapping; \
	})
	// #define _push_mapping(_id) ({ fprintf(stderr, "%u", (_id)); nullptr; })
	#define _push_op(_m, _op, _len) { \
		query_offset += push_edit(_m, (_op), &query[query_offset], (_len)); \
	}
	#define _append_op(_m, _op, _init) { \
		if((state & 0xff) == (_op)) { state += 0x100; } \
		else { _push_op(_m, state & 0xff, state>>8); state = (_op) | ((_init)<<8); } \
	}
	#define _flush_op(_m, _next_op) { \
		_push_op(_m, state & 0xff, state>>8); state = (_next_op); \
	}
    
    
    
    // Work out what region of the unpacked query sequence has a real alignment.
    // query_min_pos is first, query_max_pos is past end
    // It will be surrounded by inserts of the rest of the query.
    //
    // Account for the head_pos.query_offset.
    // If left_to_right, it is number of leading bases in alignment.sequence() not packed and not visible to dozeu.
    // Else, it is number of leading bases packd and visible to dozeu.
    // Also, if not left_to_right, alignment.sequence() and dozeu's sequence run in opposite directions.
    uint64_t query_max_pos;
    uint64_t query_min_pos;
    if (left_to_right) {
        // We end where we end in the Dozeu sequence, plus the offset from the real sequence.
        query_max_pos = aln->query_length + head_pos.query_offset;
        // We begin the distance before there accounted for by the alignment
        query_min_pos = query_max_pos - aln->match_count - aln->mismatch_count - aln->ins_count;
    } else {
        // Total packed length minus where we end in the Dozeu sequence counts from the left edge of aln.sequence
        query_min_pos = head_pos.query_offset - aln->query_length;
        // Then we advance by the number of query bases used
        query_max_pos = query_min_pos + aln->match_count + aln->mismatch_count + aln->ins_count;
    }
    debug("aligned query region: %lu-%lu", query_min_pos, query_max_pos);
    
    
	// extract query (again).
    // We're going to go through it in alignment.sequence() order, no matter the order Dozeu ran in.
	std::string const &query_seq = alignment.sequence();
	char const *query = query_seq.c_str();
    // Start a cursor at 0. If there is a leading insert, we will emit it and add the length into the cursor.
	size_t query_offset = 0;

	// set score and pos
	alignment.set_score(aln->score);
	alignment.set_identity((double)aln->match_count / (double)query_seq.length());
	alignment.set_query_position(0);		// always zero?

	// convert the result to protobuf object
	Path *path = alignment.mutable_path();
	Mapping *m = nullptr;
	if(left_to_right) {
        // The order that Dozeu gave us the alignment in (and in which we
        // filled the nodes) is left to right (i.e. the order we want to emit
        // the alignment in)
		uint64_t ref_offset = head_pos.ref_offset;
        if (query_min_pos != 0) {
            debug("leading insert of %ld bp will be required", query_min_pos);
        }
        uint64_t state = query_min_pos<<8;

		handle_t n = graph.order[aln->span[aln->span_length - 1].id];
        
		debug("rid(%u, %ld), ref_length(%lu), ref_offset(%lu), query_length(%u), query_init_length(%lu)", aln->span[aln->span_length - 1].id, graph.graph.get_id(n), graph.graph.get_length(n), ref_offset, aln->query_length, state>>8);

		state |= state == 0 ? MATCH : INS;
		for(size_t i = 0, path_offset = aln->span[0].offset; i < aln->span_length; i++) {
            debug("accounted for query up to %lu/%lu", query_offset, query_max_pos);
			dz_path_span_s const *span = &aln->span[i];
			debug("i(%lu), rid(%u, %ld), ref_length(%lu), path_offset(%lu), span->offset(%lu)", i, span->id, graph.graph.get_id(graph.order[aln->span[i].id]), graph.graph.get_length(graph.order[aln->span[i].id]), (uint64_t)path_offset, (uint64_t)span->offset);

			for(m = _push_mapping(span->id); path_offset < span[1].offset; path_offset++) {
                _append_op(m, aln->path[path_offset], 1);
			}
			_flush_op(m, aln->path[path_offset]);
		}
        // We should have made it all the way through our region of the query.
        debug("accounted for query up to %lu/%lu", query_offset, query_max_pos);
		if(m != nullptr && query_seq.length() != query_offset) {
            // We have extra query sequence that dozeu didn't actually align.
            // Treat it as a trailing insert.
            // TODO: how do we know it should be trailing?
            debug("trailing insert of %ld bp to make up length difference", query_seq.length() - query_offset);
			_push_op(m, INS, query_seq.length() - query_offset);
		}
		debug("rv: (%ld, %u) -> (%ld, %u), score(%d), %s\n", graph.graph.get_id(graph.order[aln->span[aln->span_length - 1].id]), head_pos.ref_offset, graph.graph.get_id(graph.order[aln->span[0].id]), aln->span[1].offset, aln->score, alignment.sequence().c_str());
	} else {
        // The order that Dozeu gave us the alignment in (and in which we
        // filled the nodes) is right to left (i.e. backwards). We have to flip
        // it. Also we packed the query sequence in backward, so to follow
        // along it forward, we go through the Dozeu alignment backward.
		
		uint64_t ref_offset = -((int32_t)aln->rrem);
        if (query_min_pos != 0) {
            debug("leading insert of %ld bp will be required", query_min_pos);
        }
        uint64_t state = query_min_pos<<8;
        
        handle_t n = graph.order[aln->span[aln->span_length - 1].id];
        
		debug("rid(%u, %ld), ref_length(%lu), ref_offset(%lu), query_length(%lu), query_aln_length(%u), query_init_length(%lu)", aln->span[aln->span_length - 1].id, graph.graph.get_id(n), graph.graph.get_length(n), ref_offset, query_seq.length(), aln->query_length, state>>8);

		state |= state == 0 ? MATCH : INS;
		for(size_t i = aln->span_length, path_offset = aln->path_length; i > 0; i--) {
            debug("accounted for query up to %lu/%lu", query_offset, query_max_pos);
			dz_path_span_s const *span = &aln->span[i - 1];
			debug("i(%lu), rid(%u, %ld), ref_length(%lu), path_offset(%lu), span->offset(%lu)", i, span->id, graph.graph.get_id(graph.order[aln->span[i - 1].id]), graph.graph.get_length(graph.order[aln->span[i - 1].id]), (uint64_t)path_offset, (uint64_t)span->offset);

			for(m = _push_mapping(span->id); path_offset > span->offset; path_offset--) {
				_append_op(m, aln->path[path_offset - 1], 1);
			}
			_flush_op(m, aln->path[path_offset - 1]);
		}
        // We should have made it all the way through our region of the query.
        debug("accounted for query up to %lu/%lu", query_offset, query_max_pos);
		if(m != nullptr && query_seq.length() != query_offset) {
            // We have extra query sequence that dozeu didn't actually align.
            // Treat it as a trailing insert.
            // TODO: how do we know it should be trailing?
            debug("trailing insert of %ld bp to make up length difference", query_seq.length() - query_offset);
			_push_op(m, INS, query_seq.length() - query_offset);
		}
		debug("fw: (%ld, %u) -> (%ld, %u), score(%d), %s", graph.graph.get_id(graph.order[aln->span[aln->span_length - 1].id]), -((int32_t)aln->rrem), graph.graph.get_id(graph.order[aln->span[0].id]), aln->span[1].offset, aln->score, alignment.sequence().c_str());
	}
	return;

	#undef _push_mapping
	#undef _push_op
	#undef _append_op
	#undef _flush_op
}

#if 0
void XdropAligner::debug_print(Alignment const &alignment, OrderedGraph const &graph, MaximalExactMatch const &seed, bool reverse_complemented)
{
	uint64_t seed_pos = gcsa::Node::offset(seed.nodes.front());
	uint64_t rlen = graph.graph.get_length(graph.order[id_to_index[gcsa::Node::id(seed_pos)]]);
	auto sequence = graph.graph.get_sequence(graph.order[id_to_index[gcsa::Node::id(seed_pos)]]);
    char const *rseq = sequence.c_str();
	uint64_t qlen = alignment.sequence().length(), qpos = calculate_query_seed_pos(alignment, seed);
	char const *qseq = alignment.sequence().c_str();
	fprintf(stderr, "xdrop_aligner::align, rev(%d), ptr(%p, %p), (%u, %u, %lu), (%d, %d), %s\n",
		reverse_complemented,
		&(*seed.begin), &(*alignment.sequence().begin()),
		qpos, seed.end - seed.begin, rlen,
		gcsa::Node::id(seed_pos), gcsa::Node::offset(seed_pos),
		alignment.sequence().c_str());

	if(reverse_complemented) {
		for(uint64_t i = qpos; i > 0; i--) { fprintf(stderr, "%c", comp(qseq[i - 1])); } fprintf(stderr, "\n");
		for(uint64_t i = rlen - gcsa::Node::offset(seed_pos); i > 0; i--) { fprintf(stderr, "%c", comp(rseq[i - 1])); } fprintf(stderr, "\n");
	} else {
		for(uint64_t i = qpos; i < qlen; i++) { fprintf(stderr, "%c", qseq[i]); } fprintf(stderr, "\n");
		for(uint64_t i = gcsa::Node::offset(seed_pos); i < rlen; i++) { fprintf(stderr, "%c", rseq[i]); } fprintf(stderr, "\n");
	}
	return;
}
#endif

/*
 * align query: forward-backward banded alignment
 *
 * First we find a "head" position, on the upstream side of the backing graph. If we have MEMs we do it by extending an alignment back from the most backing-downstream MEM; if we don't have MEMs then we hunt for a good match see dourdelves.
 *
 * Then we extend the head seed backing-downstream, and trace that back to find the optimal alignment.
 */
void
XdropAligner::align(
	Alignment &alignment,
	OrderedGraph const &graph,
	vector<MaximalExactMatch> const &mems,
	bool reverse_complemented)
{
	// fprintf(stderr, "called, direction(%u)\n", reverse_complemented);
	// bench_start(bench);
	// the indices are iterated by size_t (64-bit unsigned), though
	assert(graph.graph.get_node_count() < UINT32_MAX);
    // TODO: also check the edge count without iterating over all of them.

	// debug_print(alignment, graph, mems[0], reverse_complemented);

	// compute direction (currently just copied), FIXME: direction (and position) may contradict the MEMs when the function is called via the unfold -> dagify path
	bool direction = reverse_complemented;

	// extract query
	std::string const &query_seq = alignment.sequence();
	uint64_t const qlen = query_seq.length();

	// construct node_id -> index mapping table
	build_id_index_table(graph);
	forefronts.resize(graph.order.size());		// vector< void * >

	// extract seed node
	graph_pos_s head_pos;
	if(mems.empty()) {
		// seeds are not available here; probably called from mate_rescue
		build_index_edge_table(graph, direction ? graph.order.size() : 0, direction);			// we need edge information before we traverse the graph for scan seeds
		head_pos = scan_seed_position(graph, query_seq, direction);								// scan seed position mems is empty
	} else {
		// ordinary extension DP
		graph_pos_s seed_pos = calculate_seed_position(graph, mems, qlen, direction);	// we need seed to build edge table (for semi-global extension)
		build_index_edge_table(graph, seed_pos.node_index, direction);

		// pack query (upward)
		dz_query_s const *packed_query_seq_up = (direction
			? dz_pack_query_reverse(dz, &query_seq.c_str()[0],                            seed_pos.query_offset)
			: dz_pack_query_forward(dz, &query_seq.c_str()[seed_pos.query_offset], qlen - seed_pos.query_offset)
		);

		// upward extension
		head_pos = calculate_max_position(graph, seed_pos,
			extend(graph, index_edges_head.begin(), index_edges_head.end(),
				packed_query_seq_up, seed_pos.node_index, seed_pos.ref_offset,
				direction
			),
			direction
		);
	}
	// fprintf(stderr, "head_node_index(%lu), rpos(%lu, %u), qpos(%u), direction(%d)\n", head_pos.node_index, head_pos.node_index, head_pos.ref_offset, head_pos.query_offset, direction);
    
    // Now that we have determined head_pos, do the downward alignment from there, and the traceback.
    align_downward(alignment, graph, head_pos, reverse_complemented);
    
    #ifdef DEBUG
		if(mems.empty()) { fprintf(stderr, "rescue: score(%d)\n", alignment.score()); }
	#endif
    // bench_end(bench);
}
    
void
XdropAligner::align_downward(
	Alignment &alignment,
	OrderedGraph const &graph,
    graph_pos_s const &head_pos,
	bool left_to_right)
{ 

    // extract query
	std::string const &query_seq = alignment.sequence();
	uint64_t const qlen = query_seq.length();

	// pack query (downward)
	dz_query_s const *packed_query_seq_dn = (left_to_right
		? dz_pack_query_forward(dz, &query_seq.c_str()[head_pos.query_offset], qlen - head_pos.query_offset)
		: dz_pack_query_reverse(dz, &query_seq.c_str()[0],                     head_pos.query_offset)
	);

	// get head node index
	vector<uint64_t>::const_iterator begin = index_edges.begin(), end = index_edges.end();
	while(begin != end && !compare[left_to_right](_dst_index(*begin), head_pos.node_index)) { begin++; }

	// downward extension
	calculate_and_save_alignment(alignment, graph, head_pos,
		extend(graph, begin, end,
			packed_query_seq_dn, head_pos.node_index, head_pos.ref_offset,
			!left_to_right
		),
		left_to_right
	);
	dz_flush(dz);
}

void
XdropAligner::align(
	Alignment &alignment,
	Graph const &graph,
	vector<MaximalExactMatch> const &mems,
	bool reverse_complemented)
{

    // Wrap the Protobuf graph up in a ProtoHandleGraph
    ProtoHandleGraph wrapper(&graph);
    
    // Make the topological order, which we get by assuming the input Graph is
    // already in topological order.
    vector<handle_t> order;
    
    // Fill in the topological order as just the Protobuf graph's order.
    order.reserve(graph.node_size());
    for (size_t i = 0; i < graph.node_size(); i++) {
        order.push_back(wrapper.get_handle_by_index(i));
    }
    
    // Attach order to graph
    OrderedGraph ordered = {wrapper, order};
    
    // Call the backing implementation
    align(alignment, ordered, mems, reverse_complemented);
}

void
XdropAligner::align_pinned(
    Alignment& alignment, 
    const HandleGraph& g, 
    bool pin_left)
{
    // Compute our own topological order
    vector<handle_t> order = algorithms::topological_order(&g);
    
    // Align with it
    align_pinned(alignment, g, order, pin_left);
}

void
XdropAligner::align_pinned(
    Alignment& alignment, 
    const HandleGraph& g, 
    const vector<handle_t>& order,
    bool pin_left)
{
    
    if (!pin_left) {
        // If not pin_left, flip around so we can use our left-pinning-only implementation.
        // TODO: can we get around this somehow with proper use of reverse_complemented?
        
        // Flip the graph
        ReverseGraph rc(&g, true);
        
        // Flip the sequence
        string original_sequence(std::move(*alignment.mutable_sequence()));
        alignment.set_sequence(reverse_complement(original_sequence));
    
        // Fill the topological order with promoted versions of the backing graph's handles, in reverse order
        vector<handle_t> reverse_order;
        reverse_order.reserve(order.size());
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            reverse_order.push_back(rc.from_backing(*it));
        }
    
        align_pinned(alignment, rc, reverse_order, true);
        
        // Put back in the right order
        reverse_complement_path_in_place(alignment.mutable_path(), [&](id_t id) {
            // Get the length of this node
            return rc.get_length(rc.get_handle(id));
        });
        
        // Because we used a reverse complementing overlay, we have the wrong
        // idea of which orientation of each node is "forward".
        //
        // TODO: if we supported reversing edges in our internal edge indices
        // we could just flip the handles and avoid the reverse complement
        // graph and this fixup...
        for (size_t i = 0; i < alignment.path().mapping_size(); i++) {
            // Flip the orientation of each mapping
            auto mapping = alignment.mutable_path()->mutable_mapping(i);
            mapping->mutable_position()->set_is_reverse(!mapping->position().is_reverse());
        }
        
        *alignment.mutable_sequence() = std::move(original_sequence);
    } else {
        // Otherwise we have the actual implementation, for left pinning.
        
        if (order.empty()) {
            // Can't do anything with no nodes in the graph.
            return;
        }
        
        // Dozeu needs a seed position to start at, but that position doesn't necessarily actually become a match.
        // So just seed position 0 on the first node and read base 0 and start the DP from that matrix position.
       
        // Create a fixed seed match to start from.
        // Match up base 0 on node 0 in the topological sort with query base 0.
        graph_pos_s head_pos = {0, 0, 0}; 
    
        // Attach order to graph
        OrderedGraph ordered = {g, order};
        
        // construct node_id -> index mapping table
        build_id_index_table(ordered);
        forefronts.resize(ordered.order.size());		// vector< void * >
        
        // Index and order the edges for a left-to-right pass
        build_index_edge_table(ordered, head_pos.node_index, true);
    
        // Do the left-to-right alignment from the fixed head_pos seed, and then do the traceback.
        align_downward(alignment, ordered, head_pos, true);
        
        // We just pass the alignment we got through, because we didn't touch the graph.
    }
}

/**
 * end of xdrop_aligner.cpp
 */
