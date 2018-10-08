
/**
 * @file xdrop_aligner.hpp
 * @author Hajime Suzuki
 * @date 2018/03/23
 */
#ifndef VG_XDROP_ALIGNER_HPP_INCLUDED
#define VG_XDROP_ALIGNER_HPP_INCLUDED

#include <algorithm>
#include <cstdint>			/* int8_t, ... */
#include <functional>
#include <unordered_map>
#include <vector>

#include "vg.pb.h"
#include "types.hpp"
#include "mem.hpp"

// #define BENCH
// #include "bench.h"

struct dz_s;
struct dz_forefront_s;
struct dz_query_s;

namespace vg {

	struct graph_pos_s {
		size_t node_index;
		uint32_t ref_offset, query_offset;
	};

	class XdropAligner {
	private:
		// context (contains memory arena and constants) and working buffers
		struct dz_s *dz;
		std::unordered_map< id_t, uint64_t > id_to_index;			// can be lighter? index in lower 32bit and graph_id -> mem_id mapping (inverse of trans mapping)
		std::vector< uint64_t > index_edges, index_edges_head;		// (int32_t, int32_t) tuple; FIXME: index_edges and index_edges_head are partly duplicated
		std::vector< struct dz_forefront_s const * > forefronts;

		// forward and reverse comparators; [0] for forward and [1] for reverse (FIXME: can we embed them in the vtable?)
		std::function<bool (uint64_t const &, uint64_t const &)> const compare[2] = {
			[](uint64_t const &x, uint64_t const &y) -> bool { return((int64_t)x < (int64_t)y); },
			[](uint64_t const &x, uint64_t const &y) -> bool { return((int64_t)x > (int64_t)y); }
		};
		std::function<uint64_t (uint64_t const &, uint64_t const &)> const edge[2] = {
			[](uint64_t const &from, uint64_t const &to) -> uint64_t { return((to<<32) | from); },
			[](uint64_t const &from, uint64_t const &to) -> uint64_t { return((from<<32) | to); }
		};

		// working buffer init functions
		void build_id_index_table(Graph const &graph);
		void build_index_edge_table(Graph const &graph, uint32_t const seed_node_index, bool direction);

		// position handling -> (node_index, ref_offset, query_offset): struct graph_pos_s
		// MaximalExactMatch const &select_root_seed(vector<MaximalExactMatch> const &mems);
		struct graph_pos_s calculate_seed_position(Graph const &graph, vector<MaximalExactMatch> const &mems, size_t query_length, bool direction);
		struct graph_pos_s calculate_max_position(Graph const &graph, struct graph_pos_s const &seed_pos, size_t max_node_index, bool direction);
		struct graph_pos_s scan_seed_position(Graph const &graph, std::string const &query_seq, bool direction);

		size_t push_edit(Mapping *mapping, uint8_t op, char const *alt, size_t len);

		// extension -> max_node_index: size_t
		size_t extend(Graph const &graph, vector<uint64_t>::const_iterator begin, vector<uint64_t>::const_iterator end, struct dz_query_s const *packed_query, size_t seed_node_index, uint64_t seed_offset, bool direction);
		void calculate_and_save_alignment(Alignment &alignment, Graph const &graph, struct graph_pos_s const &head_pos, size_t tail_node_index, bool direction);

		// void debug_print(Alignment const &alignment, Graph const &graph, MaximalExactMatch const &seed, bool reverse_complemented);
		// bench_t bench;

	public:
		// default_* defined in vg::, see gssw_aligner.hpp
		XdropAligner();
		XdropAligner(XdropAligner const &);
		XdropAligner& operator=(XdropAligner const &);
		XdropAligner(XdropAligner&&);
		XdropAligner& operator=(XdropAligner&&);
		XdropAligner(int8_t _match,
			int8_t _mismatch,
			int8_t _gap_open,
			int8_t _gap_extension,
			int32_t _full_length_bonus,
			uint32_t _max_gap_length);
		XdropAligner(int8_t const *_score_matrix,
			int8_t _gap_open,
			int8_t _gap_extension,
			int32_t _full_length_bonus,
			uint32_t _max_gap_length);
		~XdropAligner(void);

		// copied from gssw_aligner.hpp
		void align(Alignment &alignment, Graph const &graph, const vector<MaximalExactMatch> &mems, bool reverse_complemented);
	};
} // end of namespace vg

#endif
/**
 * end of xdrop_aligner.hpp
 */
