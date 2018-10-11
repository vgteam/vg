
/**
 * @file xdrop_aligner.hpp
 * @author Hajime Suzuki
 * @date 2018/03/23
 */
#include <cstdio>
#include <assert.h>
#include "mem.hpp"
#include "xdrop_aligner.hpp"

#define DZ_FULL_LENGTH_BONUS
#define DZ_CIGAR_OP				0x04030201
enum { MISMATCH = 1, MATCH = 2, INS = 3, DEL = 4 };
#include "../deps/dozeu/dozeu.h"

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
	return(*this);
}

XdropAligner::XdropAligner(XdropAligner&& rhs)
{
	dz_destroy(dz);
	dz = rhs.dz;		// move
	rhs.dz = nullptr;
}

XdropAligner& XdropAligner::operator=(XdropAligner&& rhs)
{
	if(this == &rhs) { return(*this); }
	dz_destroy(dz);
	dz = rhs.dz;
	rhs.dz = nullptr;
	return(*this);
}

XdropAligner::XdropAligner()
{
	dz = nullptr;
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

void XdropAligner::build_id_index_table(Graph const &graph)
{
	// construct node_id -> index map
	id_to_index.clear();							// unordered_map< id_t, uint64_t >
	for(size_t i = 0; i < graph.node_size(); i++) {
		Node const &n = graph.node(i);
		id_to_index[n.id()] = (uint32_t)i;
		// debug("i(%lu), id(%ld), length(%lu)", i, n.id(), n.sequence().length());
	}
	return;
}

void XdropAligner::build_index_edge_table(Graph const &graph, uint32_t const seed_node_index, bool direction)
{
	// build (src_index, dst_index) array
	index_edges.clear();
	index_edges_head.clear();

	index_edges.reserve(graph.edge_size());			// vector< uint64_t >
	index_edges_head.reserve(graph.edge_size());	// vector< uint64_t >
	for(size_t i = 0; i < graph.edge_size(); i++) {
		Edge const &e = graph.edge(i);
		index_edges.push_back(edge[!direction](id_to_index[e.from()], id_to_index[e.to()]));
		// debug("append edge, %u -> %u", id_to_index[e.from()], id_to_index[e.to()]);

		if(!compare[!direction](id_to_index[direction ? e.from() : e.to()], seed_node_index)) { continue; }
		index_edges_head.push_back(edge[direction](id_to_index[e.from()], id_to_index[e.to()]));	// reversed
		// fprintf(stderr, "append head edge, %u -> %u\n", id_to_index[e.from()], id_to_index[e.to()]);
	}

	sort(index_edges.begin(), index_edges.end(), compare[!direction]);
	sort(index_edges_head.begin(), index_edges_head.end(), compare[direction]);
	return;
}

struct graph_pos_s XdropAligner::calculate_seed_position(Graph const &graph, vector<MaximalExactMatch> const &mems, size_t query_length, bool direction)
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
    struct graph_pos_s pos;
	MaximalExactMatch const &seed = mems[direction ? (mems.size() - 1) : 0];		// here mems is never empty

	auto seed_pos = direction ? seed.nodes.front() : seed.nodes.back();
	size_t node_id = gcsa::Node::id(seed_pos);
	size_t node_offset = gcsa::Node::offset(seed_pos);
	pos.node_index = id_to_index[node_id];

	// calc ref_offset
	Node const &n = graph.node(pos.node_index);
	pos.ref_offset = direction ? n.sequence().length() - node_offset : node_offset;

	// calc query_offset (FIXME: is there O(1) solution?)
	pos.query_offset = query_length - (seed.end - seed.begin);
	for(std::string::const_iterator p = seed.end; *p != '\0'; pos.query_offset--, p++) {}
	pos.query_offset = direction ? query_length - pos.query_offset : pos.query_offset;
	// fprintf(stderr, "calc_seed_pos, direction(%d), rpos(%lu, %u), qpos(%u), len(%lu)\n", direction, pos.node_index, pos.ref_offset, pos.query_offset, seed.end - seed.begin);
	return(pos);
}

struct graph_pos_s XdropAligner::calculate_max_position(Graph const &graph, struct graph_pos_s const &seed_pos, size_t max_node_index, bool direction)
{
	// save node id
	struct graph_pos_s pos;
	pos.node_index = max_node_index;

	// calc max position on the node
	uint64_t max_pos = (uint64_t)dz_calc_max_pos(dz, forefronts[max_node_index]);

	// ref-side offset fixup
	int32_t rpos = (int32_t)(max_pos>>32);

	Node const &n = graph.node(pos.node_index);
	pos.ref_offset = direction ? -rpos : (n.sequence().length() - rpos);

	// query-side offset fixup
	int32_t qpos = max_pos & 0xffffffff;
	pos.query_offset = seed_pos.query_offset + (direction ? -qpos : qpos);
	// fprintf(stderr, "calc_max_pos, rpos(%lu, %u), qpos(%u, %u)\n", pos.node_index, pos.ref_offset, seed_pos.query_offset, pos.query_offset);
	return(pos);
}

struct graph_pos_s XdropAligner::scan_seed_position(Graph const &graph, std::string const &query_seq, bool direction)
{
	uint64_t const qlen = query_seq.length(), scan_len = qlen < 15 ? qlen : 15;		// FIXME: scan_len should be variable

	debug("scan seeds, direction(%u), qlen(%lu), scan_len(%lu)", direction, qlen, scan_len);
	for(vector<uint64_t>::const_reverse_iterator p = index_edges.rbegin(); p != index_edges.rend(); p++) {
		// debug("i(%lu), %lu -> %lu", p - index_edges.rbegin(), _dst_index(*p), _src_index(*p));
	}
	struct dz_query_s const *packed_query = (direction
		? dz_pack_query_reverse(dz, &query_seq.c_str()[0],               scan_len)
		: dz_pack_query_forward(dz, &query_seq.c_str()[qlen - scan_len], scan_len)
	);

	int64_t inc = direction ? -1 : 1;
	size_t max_node_index = direction ? graph.node_size() - 1 : 0, prev_node_index = max_node_index - inc;
	for(vector<uint64_t>::const_reverse_iterator p = index_edges.rbegin(); p != index_edges.rend();) {
		size_t node_index = _src_index(*p), n_incoming_edges = 0;
		debug("i(%lu), node_index(%lu), max_node_index(%lu), prev_node_index(%lu)", p - index_edges.rbegin(), node_index, max_node_index, prev_node_index);

		// fill ends
		for(size_t empty_node_index = prev_node_index + inc; compare[direction](empty_node_index, node_index); empty_node_index += inc) {
			debug("fill root for empty_node_index(%lu)", empty_node_index);
			std::string const &seq = graph.node(empty_node_index).sequence();
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
		struct dz_forefront_s const *incoming_forefronts[n_incoming_edges + 1];
		for(size_t i = 0; i < n_incoming_edges; i++) {
			struct dz_forefront_s const *t = forefronts[_dst_index(p[i])];
			assert(t != nullptr);
			incoming_forefronts[i] = t;
		}

		// forward edge_index_base
		p += n_incoming_edges;

		Node const &n = graph.node(node_index);
		std::string const &ref_seq = n.sequence();
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

	struct graph_pos_s pos;
	pos.node_index = 0;
	pos.ref_offset = 0;
	pos.query_offset = direction ? scan_len : qlen - scan_len;
	struct graph_pos_s p = calculate_max_position(graph, pos, max_node_index, direction);
	debug("node_index(%lu), ref_offset(%d), query_offset(%d), max(%d)", p.node_index, p.ref_offset, p.query_offset, forefronts[max_node_index]->max);
	return(p);
}

size_t XdropAligner::extend(
	Graph const &graph,
	vector<uint64_t>::const_iterator begin,			// must be sorted ascending order for forward extension, descending order in reverse
	vector<uint64_t>::const_iterator end,
	struct dz_query_s const *packed_query,
	size_t seed_node_index,
	uint64_t seed_offset,							// offset: 0-------->L for both forward and reverse
	bool direction)									// true for forward, false for reverse
{
	// get root node
	std::string const &root_seq = graph.node(seed_node_index).sequence();

	// load position and length
	uint64_t rpos = seed_offset;
	int64_t rlen = (direction ? 0 : root_seq.length()) - seed_offset;

	debug("rpos(%lu), rlen(%ld)", rpos, rlen);
	forefronts[seed_node_index] = dz_extend(dz,
		packed_query,
		dz_root(dz), 1,
		&root_seq.c_str()[rpos], rlen, seed_node_index
	);

	int64_t inc = direction ? -1 : 1;
	size_t max_node_index = seed_node_index, prev_node_index = seed_node_index;
	debug("root: node_index(%lu, %ld), ptr(%p), score(%d)", seed_node_index, graph.node(seed_node_index).id(), forefronts[seed_node_index], forefronts[seed_node_index]->max);

	for(vector<uint64_t>::const_iterator p = begin; p != end;) {
		size_t node_index = _dst_index(*p), n_incoming_edges = 0, n_incoming_forefronts = 0;

		for(size_t empty_node_index = prev_node_index + inc; compare[direction](empty_node_index, node_index); empty_node_index += inc) {
			forefronts[empty_node_index] = nullptr;
		}
		prev_node_index = node_index;

		// pack incoming edges
		while(p + n_incoming_edges != end && _dst_index(p[n_incoming_edges]) == node_index) { n_incoming_edges++; }
		assert(n_incoming_edges > 0);

		// pack incoming edges
		struct dz_forefront_s const *incoming_forefronts[n_incoming_edges + 1];
		for(size_t k = 0; k < n_incoming_edges; k++) {
			int64_t s = _src_index(p[k]), d = _dst_index(p[k]), r = seed_node_index;
			if((s - r) * (d - r) < 0) { continue; }

			struct dz_forefront_s const *t = forefronts[_src_index(p[k])];
			if(t == nullptr || dz_is_terminated(t)) { continue; }	// skip terminated
			incoming_forefronts[n_incoming_forefronts++] = t;		// copy onto stack
		}
		debug("node_index(%lu, %ld), n_incoming_edges(%lu, %lu)", node_index, graph.node(node_index).id(), n_incoming_edges, n_incoming_forefronts);

		// forward edge_index_base
		p += n_incoming_edges;
		if(n_incoming_forefronts == 0) { forefronts[node_index] = nullptr; continue; }

		Node const &n = graph.node(node_index);
		std::string const &ref_seq = n.sequence();
		forefronts[node_index] = dz_extend(dz,
			packed_query,
			incoming_forefronts, n_incoming_forefronts,
			&ref_seq.c_str()[direction ? ref_seq.length() : 0],
			direction ? -ref_seq.length() : ref_seq.length(),
			node_index
		);
		if(forefronts[node_index]->max + (direction & dz_geq(forefronts[node_index])) > forefronts[max_node_index]->max) {
			max_node_index = node_index;
		}
		debug("node_index(%lu, %ld), n_incoming_edges(%lu, %lu), forefront(%p), range[%u, %u), term(%u), max(%d), curr(%d, %d, %u)",
			node_index, graph.node(node_index).id(), n_incoming_edges, n_incoming_forefronts,
			forefronts[node_index], forefronts[node_index]->r.spos, forefronts[node_index]->r.epos, dz_is_terminated(forefronts[node_index]),
			forefronts[max_node_index]->max, forefronts[node_index]->max, forefronts[node_index]->inc, direction & dz_geq(forefronts[node_index]));
	}
	debug("max(%p), score(%d)", forefronts[max_node_index], forefronts[max_node_index]->max);
	return(max_node_index);
}

// append an edit at the end of the current mapping array, returns forwarded length on the query
size_t XdropAligner::push_edit(
	Mapping *mapping,
	uint8_t op,
	char const *alt,
	size_t len)
{
	/* see gssw_aligner.cpp:gssw_mapping_to_alignment */
	#define _add_edit(_from_len, _to_len, _subseq) { \
		Edit *e = mapping->add_edit(); \
		e->set_from_length((_from_len)); \
		e->set_to_length((_to_len)); \
		/* expect a branch dependent on a compile-time NULL will be eliminated */ \
		if((_subseq) != nullptr) { e->set_sequence((char const *)(_subseq), (size_t)(_to_len)); } \
	}

	// fprintf(stderr, ":%lu%c", len, "-XMID"[op]);
	if(op == MISMATCH) {
		// break down into multiple SNVs
		for(size_t i = 0; i < len; i++) { _add_edit(1, 1, &alt[i]); }
	} else {
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
	Graph const &graph,
	struct graph_pos_s const &head_pos,
	size_t tail_node_index,
	bool direction)
{
	// clear existing alignment (no matter if any significant path is not obtained)
	alignment.clear_path();
	alignment.set_score(forefronts[tail_node_index]->max);
	if(forefronts[tail_node_index]->max == 0) { return; }

	// traceback
	struct dz_alignment_s const *aln = dz_trace(dz, forefronts[tail_node_index]);
	if(aln == nullptr || aln->path_length == 0) { return; }
	assert(aln->span[aln->span_length - 1].id == tail_node_index);
    if (aln->query_length > alignment.sequence().size()) {
        cerr << "[vg xdrop_aligner.cpp] Error: dozeu alignment query_length longer than sequence" << endl;
        alignment.set_score(0);
        return;
    }

	#define _push_mapping(_id) ({ \
		Node const &n = graph.node((_id)); \
		Mapping *mapping = path->add_mapping(); \
		mapping->set_rank(path->mapping_size()); \
		Position *position = mapping->mutable_position(); \
		position->set_node_id(n.id()); \
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

	// extract query (again)
	std::string const &query_seq = alignment.sequence();
	char const *query = query_seq.c_str();
	size_t query_offset = 0;

	// set score and pos
	alignment.set_score(aln->score);
	alignment.set_identity((double)aln->match_count / (double)query_seq.length());
	alignment.set_query_position(0);		// always zero?

	// convert the result to protobuf object
	Path *path = alignment.mutable_path();
	Mapping *m = nullptr;
	if(direction) {
		uint64_t ref_offset = head_pos.ref_offset, state = head_pos.query_offset<<8;

		Node const &n = graph.node(aln->span[aln->span_length - 1].id);
		// fprintf(stderr, "rid(%u, %ld), ref_length(%lu), ref_offset(%lu), query_length(%u), query_init_length(%lu)\n", aln->span[aln->span_length - 1].id, n.id(), n.sequence().length(), ref_offset, aln->query_length, state>>8);

		state |= state == 0 ? MATCH : INS;
		for(size_t i = 0, path_offset = aln->span[0].offset; i < aln->span_length; i++) {
			struct dz_path_span_s const *span = &aln->span[i];
			// fprintf(stderr, "i(%lu), rid(%u, %ld), ref_length(%lu), path_offset(%lu), span->offset(%lu)\n", i, span->id, graph.node(aln->span[i].id).id(), graph.node(aln->span[i].id).sequence().length(), (uint64_t)path_offset, (uint64_t)span->offset);

			for(m = _push_mapping(span->id); path_offset < span[1].offset; path_offset++) {
				_append_op(m, aln->path[path_offset], 1);
			}
			_flush_op(m, aln->path[path_offset]);
		}
		if(m != nullptr && query_seq.length() != query_offset) {
			_push_op(m, INS, query_seq.length() - query_offset);
		}
		// fprintf(stderr, "rv: (%ld, %u) -> (%ld, %u), score(%d), %s\n", graph.node(aln->span[aln->span_length - 1].id).id(), head_pos.ref_offset, graph.node(aln->span[0].id).id(), aln->span[1].offset, aln->score, alignment.sequence().c_str());
	} else {
		Node const &n = graph.node(aln->span[aln->span_length - 1].id);
		uint64_t ref_offset = -((int32_t)aln->rrem), state = (head_pos.query_offset - aln->query_length)<<8;
		// fprintf(stderr, "rid(%u, %ld), ref_length(%lu), ref_offset(%lu), query_length(%lu), query_aln_length(%u), query_init_length(%lu)\n", aln->span[aln->span_length - 1].id, n.id(), n.sequence().length(), ref_offset, query_seq.length(), aln->query_length, state>>8);

		state |= state == 0 ? MATCH : INS;
		for(size_t i = aln->span_length, path_offset = aln->path_length; i > 0; i--) {
			struct dz_path_span_s const *span = &aln->span[i - 1];
			// fprintf(stderr, "i(%lu), rid(%u, %ld), ref_length(%lu), path_offset(%lu), span->offset(%lu)\n", i, span->id, graph.node(aln->span[i - 1].id).id(), graph.node(aln->span[i - 1].id).sequence().length(), (uint64_t)path_offset, (uint64_t)span->offset);

			for(m = _push_mapping(span->id); path_offset > span->offset; path_offset--) {
				_append_op(m, aln->path[path_offset - 1], 1);
			}
			_flush_op(m, aln->path[path_offset - 1]);
		}
		if(m != nullptr && query_seq.length() != query_offset) {
			_push_op(m, INS, query_seq.length() - query_offset);
		}
		// fprintf(stderr, "fw: (%ld, %u) -> (%ld, %u), score(%d), %s\n", graph.node(aln->span[aln->span_length - 1].id).id(), -((int32_t)aln->rrem), graph.node(aln->span[0].id).id(), aln->span[1].offset, aln->score, alignment.sequence().c_str());
	}
	return;

	#undef _push_mapping
	#undef _push_op
	#undef _append_op
	#undef _flush_op
}

#if 0
void XdropAligner::debug_print(Alignment const &alignment, Graph const &graph, MaximalExactMatch const &seed, bool reverse_complemented)
{
	uint64_t seed_pos = gcsa::Node::offset(seed.nodes.front());
	uint64_t rlen = graph.node(id_to_index[gcsa::Node::id(seed_pos)]).sequence().length();
	char const *rseq = graph.node(id_to_index[gcsa::Node::id(seed_pos)]).sequence().c_str();
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
 * First the head (the most upstream) seed in MEMs is selected and extended downward to detect the downstream breakpoint.
 * Next the alignment path is generated by second upward extension from the downstream breakpoint.
 */
void
XdropAligner::align(
	Alignment &alignment,
	Graph const &graph,
	vector<MaximalExactMatch> const &mems,
	bool reverse_complemented)
{
	// fprintf(stderr, "called, direction(%u)\n", reverse_complemented);
	// bench_start(bench);
	// the indices are iterated by size_t (64-bit unsigned), though
	assert(graph.node_size() < UINT32_MAX);
	assert(graph.edge_size() < UINT32_MAX);

	// debug_print(alignment, graph, mems[0], reverse_complemented);

	// compute direction (currently just copied), FIXME: direction (and position) may contradict to the MEMs when the function is called via the unfold -> dagify path
	bool direction = reverse_complemented;

	// extract query
	std::string const &query_seq = alignment.sequence();
	uint64_t const qlen = query_seq.length();

	// construct node_id -> index mapping table
	build_id_index_table(graph);
	forefronts.reserve(graph.node_size());		// vector< void * >

	// extract seed node
	struct graph_pos_s head_pos;
	if(mems.empty()) {
		// seeds are not available here; probably called from mate_rescue
		build_index_edge_table(graph, direction ? graph.node_size() : 0, direction);			// we need edge information before we traverse the graph for scan seeds
		head_pos = scan_seed_position(graph, query_seq, direction);								// scan seed position mems is empty
	} else {
		// ordinary extension DP
		struct graph_pos_s seed_pos = calculate_seed_position(graph, mems, qlen, direction);	// we need seed to build edge table (for semi-global extension)
		build_index_edge_table(graph, seed_pos.node_index, direction);

		// pack query (upward)
		struct dz_query_s const *packed_query_seq_up = (direction
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

	// pack query (downward)
	struct dz_query_s const *packed_query_seq_dn = (direction
		? dz_pack_query_forward(dz, &query_seq.c_str()[head_pos.query_offset], qlen - head_pos.query_offset)
		: dz_pack_query_reverse(dz, &query_seq.c_str()[0],                     head_pos.query_offset)
	);

	// get head node index
	vector<uint64_t>::const_iterator begin = index_edges.begin(), end = index_edges.end();
	while(begin != end && !compare[direction](_dst_index(*begin), head_pos.node_index)) { begin++; }

	// downward extension
	calculate_and_save_alignment(alignment, graph, head_pos,
		extend(graph, begin, end,
			packed_query_seq_dn, head_pos.node_index, head_pos.ref_offset,
			!direction
		),
		direction
	);
	dz_flush(dz);
	// bench_end(bench);

	#ifdef DEBUG
		if(mems.empty()) { fprintf(stderr, "rescue: score(%d)\n", alignment.score()); }
	#endif
	return;
}

/**
 * end of xdrop_aligner.cpp
 */
