
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
    score_matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    for (size_t i = 0; i < 16; ++i) {
        score_matrix[i] = rhs.score_matrix[i];
    }
    gap_open = rhs.gap_open;
    gap_extend = rhs.gap_extend;
    full_length_bonus = rhs.full_length_bonus;
    max_gap_length = rhs.max_gap_length;
}

XdropAligner& XdropAligner::operator=(XdropAligner const &rhs)
{
	if (this == &rhs) {
        return *this;
    }
    free(score_matrix);
    score_matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    for (size_t i = 0; i < 16; ++i) {
        score_matrix[i] = rhs.score_matrix[i];
    }
    gap_open = rhs.gap_open;
    gap_extend = rhs.gap_extend;
    full_length_bonus = rhs.full_length_bonus;
    max_gap_length = rhs.max_gap_length;
	return *this;
}

XdropAligner::XdropAligner(XdropAligner&& rhs)
{
    score_matrix = rhs.score_matrix;
    rhs.score_matrix = nullptr;
    gap_open = rhs.gap_open;
    gap_extend = rhs.gap_extend;
    full_length_bonus = rhs.full_length_bonus;
    max_gap_length = rhs.max_gap_length;
}

XdropAligner& XdropAligner::operator=(XdropAligner&& rhs)
{
	if (this == &rhs) {
        return *this;
        
    }
    score_matrix = rhs.score_matrix;
    rhs.score_matrix = nullptr;
    gap_open = rhs.gap_open;
    gap_extend = rhs.gap_extend;
    full_length_bonus = rhs.full_length_bonus;
    max_gap_length = rhs.max_gap_length;
	return *this;
}

XdropAligner::XdropAligner()
{
	
}

XdropAligner::XdropAligner(int8_t const *_score_matrix, int8_t _gap_open, int8_t _gap_extension,
                           int32_t _full_length_bonus, uint32_t _max_gap_length) :
    gap_open(_gap_open - _gap_extension), gap_extend(_gap_extension), max_gap_length(_max_gap_length), full_length_bonus(_full_length_bonus)
{
	assert(_gap_open - _gap_extension >= 0);
	assert(_gap_extension > 0);
	assert(_full_length_bonus >= 0);
	assert(_max_gap_length > 0);

    score_matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    for (size_t i = 0; i < 16; ++i) {
        score_matrix[i] = _score_matrix[i];
    }
	// bench_init(bench);
}

XdropAligner::XdropAligner(int8_t _match, int8_t _mismatch, int8_t _gap_open, int8_t _gap_extension,
                           int32_t _full_length_bonus, uint32_t _max_gap_length) :
    gap_open(_gap_open - _gap_extension), gap_extend(_gap_extension), max_gap_length(_max_gap_length), full_length_bonus(_full_length_bonus)
{
	assert(_match > 0);
	assert(_mismatch > 0);
	assert(_gap_open - _gap_extension >= 0);
	assert(_gap_extension > 0);
	assert(_full_length_bonus >= 0);
	assert(_max_gap_length > 0);
    
    score_matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    for (size_t i = 0; i < 16; ++i) {
        if (i % 5 == 0) {
            score_matrix[i] = _match;
        }
        else {
            score_matrix[i] = -_mismatch;
        }
    }
    // bench_init(bench);
}

XdropAligner::~XdropAligner(void)
{
    free(score_matrix);
	// fprintf(stderr, "xdrop: time(%lu), count(%lu)\n", bench_get(bench) / 1000, bench_get_count(bench));
}

XdropAligner::OrderedGraph::OrderedGraph(const HandleGraph& graph, const vector<handle_t>& order) : graph(graph), order(order) {
    for (size_t i = 0; i < order.size(); ++i) {
        index_of[order[i]] = i;
    }
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


XdropAligner::graph_pos_s XdropAligner::calculate_seed_position(const OrderedGraph& graph, const vector<MaximalExactMatch>& mems,
                                                                size_t query_length, bool direction) const
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
	const MaximalExactMatch& seed = mems[direction ? (mems.size() - 1) : 0];		// here mems is never empty

	auto seed_pos = direction ? seed.nodes.front() : seed.nodes.back();
	size_t node_id = gcsa::Node::id(seed_pos);
	size_t node_offset = gcsa::Node::offset(seed_pos);
	pos.node_index = graph.index_of.at(graph.graph.get_handle(node_id));

	// calc ref_offset
	handle_t n = graph.order[pos.node_index];
	pos.ref_offset = direction ? (graph.graph.get_length(n) - node_offset) : node_offset;

	// calc query_offset (FIXME: is there O(1) solution?)
    // TODO: i think we need access to the original Alignment to make this O(1)
	pos.query_offset = query_length - (seed.end - seed.begin);
	for (auto p = seed.end; *p != '\0'; pos.query_offset--, p++) {
        // nothing to do in loop, just measuring
    }
	pos.query_offset = direction ? query_length - pos.query_offset : pos.query_offset;
	// fprintf(stderr, "calc_seed_pos, direction(%d), rpos(%lu, %u), qpos(%u), len(%lu)\n", direction, pos.node_index, pos.ref_offset, pos.query_offset, seed.end - seed.begin);
	return pos;
}

XdropAligner::graph_pos_s XdropAligner::calculate_max_position(const OrderedGraph& graph, const graph_pos_s& seed_pos, size_t max_node_index,
                                                               bool direction, dz_s* dz, const vector<const dz_forefront_s*>& forefronts) const
{
	// save node id
	graph_pos_s pos;
	pos.node_index = max_node_index;
    
    // Find the node
    handle_t n = graph.order[pos.node_index];

    assert(forefronts.at(max_node_index)->mcap != nullptr);

	// calc max position on the node
	uint64_t max_pos = (uint64_t) dz_calc_max_pos(dz, forefronts[max_node_index]);

	// ref-side offset fixup
	int32_t rpos = (int32_t)(max_pos>>32);

	pos.ref_offset = direction ? -rpos : (graph.graph.get_length(n) - rpos);

	// query-side offset fixup
	int32_t qpos = max_pos & 0xffffffff;
	pos.query_offset = seed_pos.query_offset + (direction ? -qpos : qpos);
	// fprintf(stderr, "calc_max_pos, rpos(%lu, %u), qpos(%u, %u)\n", pos.node_index, pos.ref_offset, seed_pos.query_offset, pos.query_offset);
	return pos;
}

XdropAligner::graph_pos_s XdropAligner::scan_seed_position(const OrderedGraph& graph, const string& query_seq, bool direction,
                                                           dz_s* dz, vector<const dz_forefront_s*>& forefronts) const
{
	const uint64_t qlen = query_seq.length(), scan_len = qlen < 15 ? qlen : 15;		// FIXME: scan_len should be variable

	const dz_query_s* packed_query = (direction
		? dz_pack_query_reverse(dz, &query_seq.c_str()[0],               scan_len)
		: dz_pack_query_forward(dz, &query_seq.c_str()[qlen - scan_len], scan_len)
	);

	int64_t inc = direction ? -1 : 1;
    int64_t max_idx  = direction ? graph.order.size() - 1 : 0;
    for (int64_t i = max_idx; i >= 0 && i < graph.order.size(); i += inc) {
        
        handle_t handle = graph.order.at(i);
        
        vector<const dz_forefront_s*> incoming_forefronts;
        graph.graph.follow_edges(handle, !direction, [&](const handle_t& incoming) {
            const dz_forefront_s* inc_ff = forefronts[graph.index_of.at(incoming)];
            assert(inc_ff != nullptr);
            incoming_forefronts.push_back(inc_ff);
        });
        
        auto seq = graph.graph.get_sequence(handle);
        if (incoming_forefronts.empty()) {
            forefronts[i] = dz_scan(dz, packed_query, dz_root(dz), 1,
                                    &seq.c_str()[direction ? seq.size() : 0],
                                    direction ? -seq.size() : seq.size(), i);
        }
        else {
            forefronts[i] = dz_scan(dz,packed_query, incoming_forefronts.data(), incoming_forefronts.size(),
                                    &seq.c_str()[direction ? seq.size() : 0],
                                    direction ? -seq.size() : seq.size(), i);
        }
        
        if(forefronts[i]->max + (direction & dz_geq(forefronts[i])) > forefronts[max_idx]->max) {
            max_idx = i;
        }
    }

	graph_pos_s pos;
	pos.node_index = 0;
	pos.ref_offset = 0;
	pos.query_offset = direction ? scan_len : qlen - scan_len;
	graph_pos_s p = calculate_max_position(graph, pos, max_idx, direction, dz, forefronts);
	debug("node_index(%lu), ref_offset(%d), query_offset(%d), max(%d)", p.node_index, p.ref_offset, p.query_offset, forefronts[max_node_index]->max);
	return p;
}

size_t XdropAligner::extend(const OrderedGraph& graph, const dz_query_s* packed_query,
                            const vector<graph_pos_s>& seed_positions, bool right_to_left,
                            dz_s* dz, vector<const dz_forefront_s*>& forefronts) const
{
    // seed_offset: 0-------->L for both forward and reverse
    // right_to_left: true for a right-to-left pass with left-to-right traceback, false otherwise
    
    // ensure that the forefronts are reset
    for (size_t i = 0; i < forefronts.size(); ++i) {
        forefronts[i] = nullptr;
    }
    
    // how far into the topological order we can start
    size_t start_idx = right_to_left ? 0 : graph.order.size();
    
    // seed an alignment at each of the seed positions
    for (const graph_pos_s& seed_pos : seed_positions) {
        
        // get root node
        auto root_seq = graph.graph.get_sequence(graph.order[seed_pos.node_index]);
        
        // load position and length
        int64_t rlen = (right_to_left ? 0 : root_seq.size()) - seed_pos.ref_offset;
        
        debug("extend pass: %s over %lu forefronts from %lu", right_to_left ? "right-to-left" : "left-to-right", forefronts.size(), seed_pos.node_index);
        
        debug("rpos(%lu), rlen(%ld)", seed_pos.ref_offset, rlen);
        forefronts[seed_pos.node_index] = dz_extend(dz, packed_query, dz_root(dz), 1,
                                                    root_seq.c_str() + seed_pos.ref_offset,
                                                    rlen, seed_pos.node_index);
        
        // push the start index out as far as we can
        if (right_to_left) {
            start_idx = max(start_idx, seed_pos.node_index);
        }
        else {
            start_idx = min(start_idx, seed_pos.node_index);
        }
    }

	size_t max_idx = start_idx;
	debug("root: node_index(%lu, %ld), ptr(%p), score(%d)", start_idx, graph.graph.get_id(graph.order[start_idx]), forefronts[start_idx], forefronts[start_idx]->max);
    
    int64_t inc = right_to_left ? -1 : 1;
    for (int64_t i = start_idx + inc; i < graph.order.size() && i >= 0; i += inc) {

        handle_t handle = graph.order.at(i);
        
        vector<const dz_forefront_s*> incoming_forefronts;
        graph.graph.follow_edges(handle, !right_to_left, [&](const handle_t& incoming) {
            const dz_forefront_s* inc_ff = forefronts[graph.index_of.at(incoming)];
            if (inc_ff) {
                incoming_forefronts.push_back(inc_ff);
            }
        });
        
        if (!incoming_forefronts.empty()) {
            
            // TODO: if there were multiple seed positions and we didn't choose head nodes, we
            // can end up clobbering them here, seems like it might be fragile if anyone develops this again...
            
            auto ref_seq = graph.graph.get_sequence(handle);
            forefronts[i] = dz_extend(dz, packed_query, incoming_forefronts.data(), incoming_forefronts.size(),
                                      &ref_seq.c_str()[right_to_left ? ref_seq.length() : 0],
                                      right_to_left ? -ref_seq.length() : ref_seq.length(), i);
        }
        
        if (forefronts[i] != nullptr) {
            if (forefronts[i]->max + (right_to_left & dz_geq(forefronts[i])) > forefronts[max_idx]->max) {
                max_idx = i;
            }
        }
    }
    
    // Get max query pos
    assert(max_idx <= forefronts.size());
    assert(forefronts[max_idx] != nullptr);
    
    if (forefronts[max_idx]->mcap != nullptr) {
        
        uint64_t query_max_pos = dz_calc_max_qpos(dz, forefronts[max_idx]);
        uint64_t ref_node_max_pos = dz_calc_max_rpos(dz, forefronts[max_idx]);
        
        debug("max(%p), score(%d), qpos(%lu), rpos(%lu)", forefronts[max_idx], forefronts[max_idx]->max, query_max_pos, ref_node_max_pos);
    }
    return max_idx;
}

// append an edit at the end of the current mapping array, returns forwarded length on the query
size_t XdropAligner::push_edit(Mapping *mapping, uint8_t op, char const *alt, size_t len) const
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
        // TODO: why not make this one substitution?
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

void XdropAligner::calculate_and_save_alignment(Alignment &alignment, const OrderedGraph& graph, const vector<graph_pos_s>& head_positions,
                                                size_t tail_node_index, bool left_to_right, dz_s* dz, const vector<const dz_forefront_s*>& forefronts) const
{
	
    
    
    
    // clear existing alignment (no matter if any significant path is not obtained)
	alignment.clear_path();
	alignment.set_score(forefronts.at(tail_node_index)->max);
	if(forefronts.at(tail_node_index)->max == 0) {
        // No alignment scoring anything other than 0, or not safe to try and trace back.
        // Emit a full-length insertion
        debug("no alignment; emit full length insertion");
        Mapping* m = alignment.mutable_path()->add_mapping();
        handle_t start = graph.order[head_positions.front().node_index];
        m->mutable_position()->set_node_id(graph.graph.get_id(start));
        m->mutable_position()->set_is_reverse(graph.graph.get_is_reverse(start));
        m->mutable_position()->set_offset(head_positions.front().ref_offset);
        m->set_rank(1);
        Edit* e = m->add_edit();
        e->set_from_length(0);
        e->set_to_length(alignment.sequence().size());
        e->set_sequence(alignment.sequence());
        return;
    }
        
    // If we have a nonzero score we should have a nonempty mcap on the winning forefront
    assert(forefronts.at(tail_node_index)->mcap != nullptr);
    
	// traceback. This produces an alignment in the same order as we traversed the nodes when filling them in.
	const dz_alignment_s* aln = dz_trace(dz, forefronts.at(tail_node_index));
	if(aln == nullptr || aln->path_length == 0) {
        // No alignment actually computed.
        // Emit a full-length insertion
        debug("no traceback; emit full length insertion");
        Mapping* m = alignment.mutable_path()->add_mapping();
        handle_t start = graph.order[head_positions.front().node_index];
        m->mutable_position()->set_node_id(graph.graph.get_id(start));
        m->mutable_position()->set_is_reverse(graph.graph.get_is_reverse(start));
        m->mutable_position()->set_offset(head_positions.front().ref_offset);
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
    
    // figure out which of the head positions we ended up using
    graph_pos_s head_pos;
    for (const graph_pos_s& pos : head_positions) {
        if (pos.node_index == aln->span[0].id) {
            head_pos = pos;
            break;
        }
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
	const string& query_seq = alignment.sequence();
	const char* query = query_seq.c_str();
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
void XdropAligner::debug_print(const Alignment& alignment, const OrderedGraph& graph, const MaximalExactMatch& seed, bool reverse_complemented) const
{
	uint64_t seed_pos = gcsa::Node::offset(seed.nodes.front());
	uint64_t rlen = graph.graph.get_length(graph.order[graph.index_of.at(gcsa::Node::id(seed_pos))];
	auto sequence = graph.graph.get_sequence(graph.order[graph.index_of(gcsa::Node::id(seed_pos))]);
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
void XdropAligner::align(Alignment &alignment, const HandleGraph& graph, const vector<handle_t>& topological_order,
                         const vector<MaximalExactMatch>& mems, bool reverse_complemented) const
{

    const OrderedGraph ordered_graph(graph, topological_order);
    
    dz_s* dz = dz_init(score_matrix, gap_open, gap_extend, max_gap_length, full_length_bonus);
    
	// debug_print(alignment, graph, mems[0], reverse_complemented);

	// compute direction (currently just copied), FIXME: direction (and position) may contradict the MEMs when the function is called via the unfold -> dagify path
	bool direction = reverse_complemented;

	// extract query
	const string& query_seq = alignment.sequence();
	const uint64_t qlen = query_seq.length();

	// construct node_id -> index mapping table
    vector<const dz_forefront_s*> forefronts(ordered_graph.order.size(), nullptr);
    
	// extract seed node
	graph_pos_s head_pos;
	if(mems.empty()) {
		// seeds are not available here; probably called from mate_rescue
        
        // scan seed position mems is empty
		head_pos = scan_seed_position(ordered_graph, query_seq, direction, dz, forefronts);
	}
    else {
		// ordinary extension DP
        
        // we need seed to build edge table (for semi-global extension)
		graph_pos_s seed_pos = calculate_seed_position(ordered_graph, mems, qlen, direction);

		// pack query (upward)
		dz_query_s const *packed_query_seq_up = (direction
			? dz_pack_query_reverse(dz, &query_seq.c_str()[0],                     seed_pos.query_offset)
			: dz_pack_query_forward(dz, &query_seq.c_str()[seed_pos.query_offset], qlen - seed_pos.query_offset)
		);

		// upward extension
		head_pos = calculate_max_position(ordered_graph, seed_pos,
                                          extend(ordered_graph, packed_query_seq_up,
                                                 {seed_pos},
                                                 direction, dz, forefronts),
                                          direction, dz, forefronts);
	}
	// fprintf(stderr, "head_node_index(%lu), rpos(%lu, %u), qpos(%u), direction(%d)\n", head_pos.node_index, head_pos.node_index, head_pos.ref_offset, head_pos.query_offset, direction);
    
    // Now that we have determined head_pos, do the downward alignment from there, and the traceback.
    align_downward(alignment, ordered_graph, {head_pos}, reverse_complemented, dz, forefronts);
    
    #ifdef DEBUG
		if (mems.empty()) {
            fprintf(stderr, "rescue: score(%d)\n", alignment.score());
        }
	#endif
    
    dz_destroy(dz);
    // bench_end(bench);
}
    
void XdropAligner::align_downward(Alignment &alignment, const OrderedGraph& graph, const vector<graph_pos_s>& head_positions,
                                  bool left_to_right, dz_s* dz, vector<const dz_forefront_s*>& forefronts) const
{ 

    // we're now allowing multiple graph start positions, but not multiple read start positions
    for (size_t i = 1; i < head_positions.size(); ++i) {
        assert(head_positions.at(i).query_offset == head_positions.front().query_offset);
    }
    
    // extract query
	const string& query_seq = alignment.sequence();
	const uint64_t qlen = query_seq.length();

	// pack query (downward)
	const dz_query_s* packed_query_seq_dn = (left_to_right
		? dz_pack_query_forward(dz, &query_seq.c_str()[head_positions.front().query_offset], qlen - head_positions.front().query_offset)
		: dz_pack_query_reverse(dz, &query_seq.c_str()[0], head_positions.front().query_offset)
	);

	// downward extension
	calculate_and_save_alignment(alignment, graph, head_positions,
                                 extend(graph, packed_query_seq_dn,
                                        head_positions,
                                        !left_to_right, dz, forefronts),
                                 left_to_right, dz, forefronts);
    
	dz_flush(dz);
}

void XdropAligner::align(Alignment &alignment, const HandleGraph& graph, const vector<MaximalExactMatch>& mems,
                         bool reverse_complemented) const
{
    align(alignment, graph, algorithms::topological_order(&graph), mems, reverse_complemented);
}

void XdropAligner::align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left) const
{
    // Compute our own topological order
    vector<handle_t> order = algorithms::topological_order(&g);
    
    // Align with it
    align_pinned(alignment, g, order, pin_left);
}

void XdropAligner::align_pinned(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order, bool pin_left) const
{
    if (order.empty()) {
        // Can't do anything with no nodes in the graph.
        return;
    }
    
    // Make a new dozeu instance
    dz_s* dz = dz_init(score_matrix, gap_open, gap_extend, max_gap_length, full_length_bonus);
    
    // Dozeu needs a seed position to start at, but that position doesn't necessarily actually become a match.
    
    // Find all of the tips that we'd want to pin at
    vector<graph_pos_s> head_positions;
    for (size_t i = 0; i < order.size(); ++i) {
        handle_t handle = order[i];
        // check if this is a tip in the correct direction so that we'd want to pin on it
        bool do_pinning = g.follow_edges(handle, pin_left, [](const handle_t& neighbor) { return false; });
        
        if (do_pinning) {
            head_positions.emplace_back();
            head_positions.back().node_index = i;
            if (pin_left) {
                head_positions.back().ref_offset = 0;
                head_positions.back().query_offset = 0;
            }
            else {
                head_positions.back().ref_offset = g.get_length(handle);
                head_positions.back().query_offset = alignment.sequence().size();
            }
        }
    }
    
    
    // Attach order to graph
    OrderedGraph ordered(g, order);
    
    // construct node_id -> index mapping table
    vector<const dz_forefront_s*> forefronts(ordered.order.size(), nullptr);
    
    // Do the left-to-right alignment from the fixed head_pos seed, and then do the traceback.
    align_downward(alignment, ordered, head_positions, pin_left, dz, forefronts);

}

/**
 * end of xdrop_aligner.cpp
 */
